clear
% clc

rng(0)

%% part 1: find transformation into cartesian coordinates 

w = sym('w', [2 1]); 
z = sym('z', [2 1]); 
x = sym('x', [2 1]); 
y = sym('y', [2 1]); 
syms rho theta 
syms sigma_rho sigma_theta 

% rho = sqrt( x(1)^2 + x(2)^2 );
% theta = atan(x(2), x(1)); 

y(1) = rho * cos(theta); 
y(2) = rho * sin(theta); 

R = [ sigma_rho^2 0; 0 sigma_theta^2 ]; 

%% part 2: find covariance matrix 

% linearize 
y_lin = jacobian(y, [rho theta]); 

% covariance of linearized y 
Rc = y_lin * R * y_lin.'; 
Rc_fun = matlabFunction(Rc); 

%% part 3: simulate N = 100 realizations 

rho0 = 10^5; 
theta0 = 45*pi/180; % deg --> rad 

z1 = rho0; 
z2 = theta0; 

yc_true = [ z1 * cos(z2); z1 * sin(z2) ]; 

% polar covariance 
rho_sigma = 100; 
theta_sigma = 0.5 * pi/180; 
R = diag([rho_sigma^2; theta_sigma^2]); 

% create noise 
N = 100; 
w = mvnrnd([0; 0], R, N); 

% create measurements 
z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

% transform polar measurements to cartesian 
for i = 1:length(w)
    yc(:,i) = [ z1(i) * cos(z2(i)); z1(i) * sin(z2(i)) ]; 
end 

yc = yc'; 

%% part 4: analyze whether the errors w_c are zero mean 

w_c = yc - yc_true'; 
disp('w_c is not zero-mean') 

N
mean(w_c)

%% part 5: covariance test statistic 

Rc_double = Rc_fun( rho0, rho_sigma, theta_sigma , theta0 ); 

% calculate error 
for i = 1:length(w_c)
    err(i,:) = w_c(i,:) * inv(Rc_double) * w_c(i,:)'; 
end 

% calculate mean of error as test statistic 
err_mean = mean(err); 

N = length(w); 
Nx = length(x); 
Nz = length(z); 

% NEED STATISTICS TOOLBOX 
a = .01; 
r1 = chi2inv( a/2, N * Nz) / N; 
r2 = chi2inv( 1 - a/2, N * Nz ) / N; 

if err_mean > r1 && err_mean < r2 
    disp('Covariance matrix can be accepted as correct!') 
else
    disp('Covariance matrix cannot be accepted as correct') 
end 







