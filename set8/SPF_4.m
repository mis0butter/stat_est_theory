clear 
% clc 

load problem4data.mat 
load problem4truth.mat 

rng(0)

% noise on robot wheel encoders 
Q = diag( [ 0.1, 5*pi/180 ] )^2; 

% noise on robot sonar 
R = diag( [ 1 1 1 ] )^2; 

% # particles 
Ns  = 1000;

% time step 
dt = encoder(2).t - encoder(1).t; 

% state size 
nx = 3; 

% grab truth states 
truth = []; 
for i = 1:length(robot)
    truth = [ truth; robot(i).x' ]; 
end 

%% initialize 

% draw initial particles from initial uniform probability density,
% initialize weights equally 
r0     = unifrnd(minx, maxx, [Ns, 2]); 
theta0 = rand(Ns, 1) * 2*pi; 
w_k0   = ones(Ns, 1) / Ns; 

XX_k  = [ r0, theta0 ]; 
w_k   = w_k0; 

x_hat = []; 
P     = []; 

% start figure 
fname = 'Robot Particle Filtering'; 
hf = figure('name', fname); 
    title(fname) 
    hold on; 
    xlim([minx - 1, maxx + 2]) 
    ylim([miny - 1, maxy + 2])

%% PARTICLE FILTER 

% measurement index 
for k = 1 : length(encoder) 
    
    % particle filter 
    [x_khat, P_k, XX_k, w_k] = particle_filter(k, w_k, Q, R, Ns, XX_k, beacons, encoder, sonar, nx); 
    
    % save outputs 
    x_hat    = [x_hat; x_khat]; 
    P(:,:,k) = P_k; 
    
    % update figure 
    gcf = hf; 
    cla 
    scatter(XX_k(:,1), XX_k(:,2), 'g'); 
    scatter(truth(1:k,1), truth(1:k,2), 'b');
    scatter(x_hat(:,1), x_hat(:,2) , 'r'); 

    legend('particles', 'truth', 'est', 'location', 'southeast') 
    pause(0.05)
    
end 

%% subfunctions 

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

    % evaluate effective # of particles 
    w_sq_sum = 0; 
    for i = 1:Ns 
        w_sq_sum = w_sq_sum + w_kp1(i)^2; 
    end 
    Ns_hat = 1 / w_sq_sum; 

    % resample if necessary 
    if Ns_hat < Ns / 2
        [XX_kp1, w_kp1] = resample(XX_kp1, w_kp1, Ns); 
    end 

    % Compute weighted estimate 
    x_khatp1 = zeros(1, nx); 
    P_kp1    = zeros(nx); 
    for i = 1:Ns 
        x_khatp1 = x_khatp1 + w_kp1(i) * XX_kp1(i,:);
        xtilde   = (XX_kp1(i,:) - x_khatp1)'; 
        P_kp1    = P_kp1 + w_kp1(i) * xtilde * xtilde';       % outer product 
    end 

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

    end 

    XX_kp1 = XX_kp1_new; 
    w_kp1  = ones(Ns,1) / Ns; 

    % Normalize weights 
    w_kp1 = w_kp1 ./ sum(w_kp1); 

end 




