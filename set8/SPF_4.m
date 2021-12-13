%% SPF 4: PARTICLE FILTER 

clear 
% clc 

load problem4data.mat 
load problem4truth.mat 

% rng(0)

% noise 
Q = diag( [ 0.1, 5*pi/180 ] )^2;    % robot wheel encoders  
R = diag( [ 1 1 1 ] )^2;            % robot sonar 

% grab truth states 
x_truth = [];  t = []; 
for i = 1:length(robot)
    x_truth = [ x_truth; robot(i).x' ]; 
    t       = [ t; robot(i).t ]; 
end 

%% initialize 

% # particles 
Ns  = 1000;

% state size 
nx = 3; 

% draw initial particles from initial uniform probability density,
% initialize weights equally 
r0     = unifrnd(minx, maxx, [Ns, 2]); 
theta0 = rand(Ns, 1) * 2*pi; 
w_k0   = ones(Ns, 1) / Ns; 

% start figure 
fname = 'Robot Particle Filtering'; 
hf_map = figure('name', fname); 
    sgtitle(fname) 
    hold on; 
    xlim([minx - 1, maxx + 2]) 
    ylim([miny - 1, maxy + 2])
    xlabel('X (m)'); ylabel('Y (m)') 

%% PARTICLE FILTER 

XX_k  = [ r0, theta0 ]; 
w_k   = w_k0; 

x_hat = []; 
P     = []; 

% measurement index 
for k = 1 : length(encoder) 
    
    % particle filter 
    [x_khat, P_k, XX_k, w_k] = particle_filter(k, w_k, Q, R, Ns, XX_k, beacons, encoder, sonar, nx); 
    
    % save outputs 
    x_hat    = [x_hat; x_khat]; 
    P(:,:,k) = P_k; 
    
    % update figure 
    gcf = hf_map; 
        cla 
        scatter(XX_k(:,1), XX_k(:,2), 4, 'g'); 
        scatter(x_truth(1:k,1), x_truth(1:k,2), 12, 'b', 'filled');
        scatter(x_hat(:,1), x_hat(:,2) , 12, 'r', 'filled'); 

        legend('particles', 'truth', 'est', 'location', 'southeast') 
        pause(0.02)
    
end

%% final analysis plot 

% extract std devs 
x_sigma = sqrt(squeeze(P(1,1,:))); 
y_sigma = sqrt(squeeze(P(2,2,:))); 
theta_sigma = sqrt(squeeze(P(3,3,:))); 

fname = 'Robot Particle Filtering: Truth vs. Estimate';
n = 3; p = 1; 
pos = [100 100 600 600]; 
figure('name', fname, 'position', pos) 
hold on; grid on; 

    % x compare 
    subplot(n,p,1) 
    hold on; 
        plot_lines(1, t, x_truth, x_hat, x_sigma)
        title('X Compare') 
        legend('truth', 'est', 'est +/- \sigma', 'location', 'best')  
        ylabel('X (m)') 

    % x diff compare 
    subplot(n,p,2) 
    hold on; 
        plot_lines(2, t, x_truth, x_hat, x_sigma)
        title('Y Compare')  
        ylabel('Y (m)') 
    
    % y compare 
    subplot(n,p,3)
    hold on; 
        plot_lines(3, t, x_truth, x_hat, x_sigma)
        title('\Theta Compare') 
        xlabel('Time (s)') 
        ylabel('rad') 
        
    sgtitle(fname)
    
%% PART B     
% b. Compare your state estimate to the true state, stored in problem4truth.mat. In
% particular, consider the following points:
% i. Does your state estimate time history change much if you run the filter more
% than once? Why might that be good or bad? Hint: try running the filter a few
% times with only 100 particles.
% ii. Why do clusters of particles sometimes persist in incorrect locations on the
% map?
% iii. Why would it be difficult to implement this filter as an extended Kalman
% Filter?


    

%% subfunctions 

function plot_lines(i, t, x_truth, x_hat, x_sigma)

    plot(t, x_truth(:,i), 'b', 'linewidth', 2); 
    plot(t, x_hat(:,i), 'r--', 'linewidth', 2); 
    plot(t, x_hat(:,i) + x_sigma, 'r--'); 
    plot(t, x_hat(:,i) - x_sigma, 'r--'); 

end 

function [x_khatp1, P_kp1, XX_kp1, w_kp1] = particle_filter(k, w_k, Q, R, Ns, XX_k, beacons, encoder, sonar, nx) 

    % extract coder command 
    uk = encoder(k).u;      uk = uk'; 
    vk = covdraw(Q, Ns);    vk = vk'; 
    
    % propagate state 
    XX_kp1 = robot_dyn(uk, vk, Q, Ns, XX_k); 
    
    % measurement model 
    Z_mdl = Z_mdl_fn(XX_kp1, beacons, Ns); 

    % Calculate innovation 
    sonar_k = sonar(k).z'; 
    nu_k    = Z_mdl - sonar_k; 
    
    % update weights 
    w_kp1 = update_weights(Ns, nu_k, R, w_k); 

    % evaluate effective # of particles 
    w_sq_sum = sum(w_kp1.^2); 
    Ns_hat = 1 / w_sq_sum; 

    % resample if necessary 
    if Ns_hat < Ns / 2
        [XX_kp1, w_kp1] = resample(XX_kp1, w_kp1, Ns); 
    end 

    % Weighted state and covariance 
    x_khatp1 = sum(w_kp1 .* XX_kp1); 
    xtilde   = XX_kp1 - x_khatp1; 
    P_kp1    = xtilde' * (xtilde .* w_kp1); 
    
    
end 

function w_kp1 = update_weights(Ns, nu_k, R, w_k)

    % Recalculate weights 
    w_kp1 = zeros(Ns, 1); 
    for i = 1:Ns 
        
        nu_ki = nu_k(i,:)'; 
        nR = length(R); 

        % pdf 
        p_ki = exp( -1/2 * nu_ki' * R^-1 * nu_ki ); 
        
        % log of recalculated weight 
        w_kp1_ln(i) = log(p_ki) + log(w_k(i)); 

    end 
    
    % Update according to log likelihood 
    for i = 1:Ns 
        w_kp1(i) = exp( w_kp1_ln(i) - max(w_kp1_ln) ); 
    end 

    % Normalize weights 
    w_kp1 = w_kp1 ./ sum(w_kp1); 

end 

function XX_kp1 = robot_dyn(uk, vk, Q, Ns, XX_k)

    % add noise to distance and angle cmds 
    uk = uk + vk; 

    % determine xa and xb change 
    ds_k     = uk(:,1);     % distance delta 
    dtheta_k = uk(:,2);     % angle delta 
    theta_k  = XX_k(:,3);   % OG angle 

    dxa_k = ds_k .* cos(theta_k + dtheta_k); 
    dxb_k = ds_k .* sin(theta_k + dtheta_k); 

    % propagate dynamics 
    x_k   = XX_k(:,1:2); 
    x_kp1 = x_k + [dxa_k, dxb_k]; 
    theta_kp1 = theta_k + dtheta_k ; 

    % propagated state 
    XX_kp1 = [x_kp1, theta_kp1]; 

end 

function Z_mdl_min3 = Z_mdl_fn(XX_kp1, beacons, Ns)

    % for beacon index 
    for i = 1:5
        dxa = XX_kp1(:,1) - beacons(i,1); 
        dxb = XX_kp1(:,2) - beacons(i,2); 
        Z_mdl_all(:,i) = sqrt(dxa.^2 + dxb.^2); 
    end 

    % Create measurement model with min 3 ranges 
    Z_mdl_min3 = zeros(Ns, 3); 
    for i = 1:Ns 
        min3 = sort(Z_mdl_all(i,:)); 
        min3 = min3(1:3); 
        Z_mdl_min3(i,:) = min3; 
    end 

end 

function [XX_kp1, w_kp1] = resample(XX_kp1, w_kp1, Ns)

    % Cumulative distribution function 
    XX_kp1_new = XX_kp1; 
    w_kp1_cdf = cumsum(w_kp1); 

    % for each particle 
    for pi = 1:Ns 

        % choose random number [0,1]
        n_rand = rand;  

        % loop stops right before n_rand exceeds w_kp1_cdf threshold 
        wi = 1; 
        while n_rand > w_kp1_cdf(wi) && wi < Ns 
            wi = wi + 1; 
        end 
        XX_kp1_new(pi,:) = XX_kp1(wi,:); 
        
%         % for each patricle, increment m until we meet the criteria
%         for m = 1:Ns
%             if (sum(w_kp1(1:m-1)) <= n_rand) && (n_rand < sum(w_kp1(1:m)))
%                   % if the criteria is met, X_new(:,i) = X_new(:,m) 
%                   XX_kp1_new(pi,:)     = XX_kp1(m,:);
%                   % wik_new gets 1/Ns
% %                   wik_new(ii)      = 1/Ns;
%                 break;
%             end
%         end

    end 

    XX_kp1 = XX_kp1_new; 
    w_kp1  = ones(Ns,1) / Ns; 

    % Normalize weights 
    w_kp1 = w_kp1 ./ sum(w_kp1); 

end 





