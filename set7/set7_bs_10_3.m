clear; clc

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

y(1) = z(1) * cos(z(2)) + w(1); 
y(2) = z(1) * sin(z(2)) + w(2); 

R = [ sigma_rho^2 0; 0 sigma_theta^2 ]; 

% linearize 
y_lin = jacobian(y, z); 

% covariance of linearized y 
Rc = y_lin * R * y_lin.'; 

%% part 2: find covariance matrix 
rho0 = 10^5; 
theta0 = 45*pi/180; % deg --> rad 

z1 = rho0; 
z2 = theta0; 

y_true = [ z1 * cos(z2); z1 * sin(z2) ]; 

%% part 3: simulate N = 100 realizations 
o_rho = 100; 
o_theta = 0.5 * pi/180; 
R = diag([o_rho^2; o_theta^2]); 

w = mvnrnd([0; 0], R, 100); 

z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

clear y 
for i = 1:length(w)
    y(:,i) = [ z1(i) * cos(z2(i)); z1(i) * sin(z2(i)) ]; 
end 

y = y'; 

%% part 4: analyze whether the errors w_c are zero mean 
w_c = y - y_true'; 

%% part 5: covariance test statistic 

Rc_subs = Rc; 
Rc_subs = subs(Rc_subs, z(1), rho0); 
Rc_subs = subs(Rc_subs, z(2), theta0); 
Rc_subs = subs(Rc_subs, sigma_rho, o_rho); 
Rc_subs = subs(Rc_subs, sigma_theta, o_theta); 

Rc_double = double(Rc_subs); 

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







