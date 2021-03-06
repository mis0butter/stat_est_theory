
%% 

rho0 = 0; 
theta0 = 0; 

% polar covariance 
o_rho = 100; 
o_theta = 0.5 * pi/180; 
R = diag([o_rho^2; o_theta^2]); 

% create noise 
N = 100; 
w = mvnrnd([0; 0], R, N); 

% create measurements 
z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

%% 

N = 100; 
xbar = ones(N, 1); 
Pxx  = diag(ones(N,1)); 
R = chol(Pxx); 

w1 = mvnrnd(xbar, Pxx, N); 
w2 = rand(N,1)' * R; 

