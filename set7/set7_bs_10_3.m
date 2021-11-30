clear;

rng(0)

%% part 1: find transformation into cartesian coordinates 
w = sym('w', [2 1]); 
z = sym('z', [2 1]); 
x = sym('x', [2 1]); 
y = sym('y', [2 1]); 
syms rho theta 

% rho = sqrt( x(1)^2 + x(2)^2 );
% theta = atan(x(2), x(1)); 

y(1) = rho * cos(theta) + w(1); 
y(2) = rho * sin(theta) + w(2); 

% linearize 
y_lin = jacobian(y, [rho theta]); 

% covariance of linearized y 
cov(y_lin) 

%% part 2: find covariance matrix 
rho0 = 10^5; 
theta0 = 45; % deg 

z1 = rho0; 
z2 = theta0; 

y_true = [ z1 * cosd(z2); z1 * sind(z2) ]; 

%% part 3: simulate N = 100 realizations 
R = diag([100^2; 0.5^2]); 

w = mvnrnd([0; 0], R, 100); 

z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

for i = 1:length(w)
    y(:,i) = [ z1(i) * cosd(z2(i)); z1(i) * sind(z2(i)) ]; 
end 

y = y'; 

%% part 4: analyze whether the errors w_c are zero mean 
mean(y)' - y_true

%% part 5: covaraince test statistic 








