% Load data
load('problem4data.mat');
load('problem4truth.mat');
dt  = 0.1;
Qk  = diag([0.1, (5*pi/180)])^2;
Rk  = diag([1 1 1])^2;

% set number of particles to N = 1000
Ns  = 1000;

% threshold for resampling (as defined by problem)
resamplingTheshold = Ns/2;

% get initial particle distirbution
% sample from uniform distribution
mintheta = 0;
maxtheta = 2*pi;

% initialize x0,y0,theta0
x0       = unifrnd(minx, maxx, [Ns,1]);
y0       = unifrnd(miny, maxy, [Ns,1]);
theta0   = unifrnd(mintheta, maxtheta, [Ns,1]);

% define initial particles
Xi0 = [x0'; 
       y0'; 
       theta0'];
   
% Execute Filtering
[xhat,P] = executeParicleFiltering(...
            Xi0, Rk, Qk, Ns, encoder, beacons, sonar, robot, resamplingTheshold);
        
% plot results
plotComparison(xhat, P, robot)

%% subfunctions 

function [xhat,P] = executeParicleFiltering(...
                Xi0, Rk, Qk, Ns, encoder, beacons, sonar, robot, resamplingTheshold)
% number of states
nx = 3;
nv = 2;
nz = 3;
% get the initial weights
% set to 1/N for each
wik0     = repmat(1/Ns, Ns, 1);
% get number of steps
numSteps = numel(encoder);
% history
xhat     = zeros(nx, 1,numSteps);
P        = zeros(nx,nx,numSteps);
% initialize first value for weights
wik      = wik0;
% inialize particles
Xi       = Xi0;
hFig     = figure("Name","Robot Tracking");
xtrue12k    = [];  
for k=1:numSteps
    % first extract encoder inputs for timestep k
    uk = encoder(k).u; 
    
    % sample process noise from gaussian distribution
    vk = zeros(nv,Ns);
    for ii = 1:Ns
        % sample process noise from a normal distribtuon
        vk(:,ii) = mvnrnd(zeros(nv,1), Qk, 1)';
        % propagate particles forward from dynamics model
        % using current process noise and current 
        Xi(:,ii) = robotDynamics(Xi(:,ii), vk(:,ii), uk);
    end
    % iterate through each particle and evaluate max 3 ranges for each
    hki = zeros(3,Ns);
    for ii=1:Ns
        hki(:,ii) = findSmallestRanges(Xi(:,ii), beacons);
    end
    % extract sonar measurements for time step k
    zk          = sonar(k).z;
    
    % fine wprime from notes (note this is the natural log of wkprime)
    lnWkiPrime  = zeros(Ns,1);
    for ii=1:Ns
        lnWkiPrime(ii) = log( wik(ii) ) ...
            - 0.5*(zk -  hki(:,ii))'*inv(Rk)*(zk -  hki(:,ii));
%         lnWkiPrime(ii) = log( wik(ii) ) ...
%             + 0.5*(zk -  hki(:,ii))'*inv(Rk)*(zk -  hki(:,ii));

    end
    
    % log of lnWkiprime is max(lnWkiprime)
    lnWkiStar     = max( lnWkiPrime );
    % wtilde is the exp{ log(wmiprime(i)) - log(wkiStar) }
    wkiTilde      = exp( lnWkiPrime - lnWkiStar );
    
    % find new weights by normalizing samples
    wikp1         = wkiTilde/sum(wkiTilde);
    % find the effective number of particles
    Nhat          = 1/sum(wikp1.^2);
    % perform resampling if needed
    if Nhat <= resamplingTheshold
        Xi_new    = zeros(nx,Ns);
        wik_new   = zeros(Ns,1);
        % find m such that wikp1(1:m-1) <= eta < wikp1(1:m-1)
        % iterate through each paticle
        for ii=1:Ns
            % sample from unifrom distribution to get eta
            eta     = unifrnd(0,1);
            % for each patricle, increment m until we meet the criteria
            for m=1:Ns
                if (sum(wikp1(1:m-1)) <= eta) && (eta < sum(wikp1(1:m)))
                      % if the criteria is met, X_new(:,i) = X_new(:,m) 
                      Xi_new(:,ii)     = Xi(:,m);
                      % wik_new gets 1/Ns
                      wik_new(ii)      = 1/Ns;
                    break;
                end
            end
        end
        % delete old weights and Xi's
        Xi      = Xi_new;
        wikp1   = wik_new;
    end
    % update weights for next iteration
    wik     = wikp1;
    % evaluate xhat(k)
    xhatk   = zeros(nx,1);
    for ii=1:Ns
        xhatk = xhatk + wik(ii)*Xi(:,ii);
    end
    % evaluate P(k)
    Pk = zeros(nx,nx);
    for ii=1:Ns
        Pk = Pk + wik(ii)*(Xi(:,ii) - xhatk)*(Xi(:,ii) - xhatk)';
    end
    % update time history for best estimate and continue
    xhat(:,:,k) = xhatk;
    P(:,:,k)    = Pk; 
    % get true state for plotting
    xtrue12k = [xtrue12k robot(k).x];
    plotEstimate(k , xhat, xtrue12k, Xi, hFig)
end
end
function  [hX] = findSmallestRanges(Xi, beacons)
    % @brief: findSmallestRanges
    % helper function that calculates a vector which returns the 3
    % largest points 
    % for the given particle find the minimum ranges to the beacon
    % returns the 3 smallest values from Xi from the beacons
    
    Nb = length(beacons);
    
    ranges = [];
    for ii=1:Nb
        % x, y elements are first two elements of particle filter
        rangei  = norm(Xi(1:2) - beacons(ii,:)');
        ranges  = [ranges; rangei];
    end
    % sort ranges in ascending order
    ranges  = sort(ranges, 1, 'ascend');
    % hi is the 3 smallest values to the particle
    hX      = ranges(1:3);
end
function [Xi_new] = robotDynamics(Xi, vk, uk)
    % propagates forward robot dynamics
    % adds process noise to input 
    % Provides non-linear map to propagate 
    % Xi(k) -> Xi(k+1)
    % apply noise to measurement
    uk      = uk + vk;
    % extract ds/dtheta
    ds      = uk(1);
    dtheta  = uk(2);
    % get current value of theta
    theta   = Xi(3);
    
    % provide nonlinear map to propagate model forward  
    % use current value of theta to evaluate sin/cos functions
    dX      = [ds*cos(theta + dtheta);
               ds*sin(theta + dtheta);
               dtheta];
    
    % get new particle
    Xi_new = Xi + dX;
end
function plotEstimate(k, xhat, xtrue, Xi, h)
    % use same figure handle or create a new one
    figure(h);
    plot(Xi(1,:),Xi(2,:),'bo','LineWidth', 1');
    hold on;
    plot(xhat(1,1:k), xhat(2,1:k),'r--o','LineWidth', 2);
    hold on;
    plot(xtrue(1,1:k), xtrue(2,1:k),'g--o','LineWidth', 3);
    hold off;
    xlim([0, 10])
    ylim([0, 10])
    grid on;
    title('Propgation of Robot Dynamics')
    ylabel('$y$','Interpreter','latex')
    xlabel('$x$','Interpreter','latex')
    legend({'$\chi_i(k)$', ...
            '$\hat{x}(k)$',...
            '$x_{truth}(k)$'},...
            'Interpreter','latex', 'Location','southeast')
end
function plotComparison(xhat, P, robot)
    % @brief: plotComparison 
    % plot comparison for robot position and truth values
    figure("Name","Robot Tracking Comparison");
    steps = 1:numel(robot);
    
    xTrueHist   = [];
    xhatHist    = [];
    Pkhist      = [];
    for ii=1:numel(robot)
        xTrueHist   = [xTrueHist robot(ii).x];
        xhatHist    = [xhatHist xhat(:,:,ii)];
        Pkhist      = [Pkhist  diag(P(:,:,ii))];
    end
    subplot(1,3,1)
    plot(steps, xhatHist(1,:),'r--','LineWidth', 2);
    hold on;
    plot(steps, xTrueHist(1,:),'b-','LineWidth', 2);
    hold on;
    plot(steps, xhatHist(1,:) + sqrt(Pkhist(1,:)),'k--')
    hold on;
    plot(steps, xhatHist(1,:) - sqrt(Pkhist(1,:)),'k--')
    ylabel('$x$','Interpreter','latex')
    xlabel('$k$','Interpreter','latex')
    legend({'$\hat{x}$', ...
            '$x_{truth}$',...
            '$\hat{x} \pm \sigma_{\hat{x}}$'},...
            'Interpreter','latex', 'Location','southeast')
    grid on;
    subplot(1,3,2)
    title('Comparison of Estimates vs. Truth')
    plot(steps, xhatHist(2,:),'r--','LineWidth', 2);
    hold on;
    plot(steps, xTrueHist(2,:),'-b','LineWidth', 2);
    hold on;
    plot(steps, xhatHist(2,:) + sqrt(Pkhist(2,:)),'k--')
    hold on;
    plot(steps, xhatHist(2,:) - sqrt(Pkhist(2,:)),'k--')
    grid on;
    ylabel('$y$','Interpreter','latex')
    xlabel('$k$','Interpreter','latex')
    legend({'$\hat{y}$', ...
            '$y_{truth}$',...
            '$\hat{y} \pm \sigma_{\hat{y}}$'},...
            'Interpreter','latex', 'Location','southeast')
    
    subplot(1,3,3)
    plot(steps, xhatHist(3,:),'r--','LineWidth', 2');
    hold on;
    plot(steps, xTrueHist(3,:),'b-','LineWidth', 2');
    hold on;
    plot(steps, xhatHist(3,:) + sqrt(Pkhist(3,:)),'k--')
    hold on;
    plot(steps, xhatHist(3,:) - sqrt(Pkhist(3,:)),'k--')
    grid on;
    ylabel('$\theta$','Interpreter','latex')
    xlabel('$k$','Interpreter','latex')
    legend({'$\hat{\theta}$', ...
            '$\theta_{truth}$',...
            '$\hat{\theta} \pm \sigma_{\hat{\theta}}$'},...
            'Interpreter','latex', 'Location','southeast')
end