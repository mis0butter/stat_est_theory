clear
% clc

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
Rc_fn = matlabFunction(Rc); 

%% run sims 

mean_wc_arr = []; 
for N = 100 : 10 : 10000
    
    mean_wc = calc_wc(N, Rc_fn); 
    
    mean_wc_arr = [mean_wc_arr; mean_wc]; 
    
end 

%% 

plot(mean_wc_arr) 
xlabel('N') 
title('w_c mean') 
legend('w_c(1)', 'w_c(2)')

%% functions 

function mean_wc = calc_wc(N, Rc_fn)

% Rc_fn = @(rho,sigma_rho,sigma_theta,theta)reshape([sigma_rho.^2.*cos(theta).^2+rho.^2.*sigma_theta.^2.*sin(theta).^2,sigma_rho.^2.*cos(theta).*sin(theta)-rho.^2.*sigma_theta.^2.*cos(theta).*sin(theta),sigma_rho.^2.*cos(theta).*sin(theta)-rho.^2.*sigma_theta.^2.*cos(theta).*sin(theta),sigma_rho.^2.*sin(theta).^2+rho.^2.*sigma_theta.^2.*cos(theta).^2],[2,2]); 


%% part 3: simulate N = 100 realizations 

rho0 = 10^5; 
theta0 = 45*pi/180; % deg --> rad 

z1 = rho0; 
z2 = theta0; 

yc_true = [ z1 * cos(z2); z1 * sin(z2) ]; 

o_rho = 100; 
o_theta = 0.5 * pi/180; 
R = diag([o_rho^2; o_theta^2]); 

w = mvnrnd([0; 0], R, N); 

z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

for i = 1:length(w)
    yc(:,i) = [ z1(i) * cos(z2(i)); z1(i) * sin(z2(i)) ]; 
end 

yc = yc'; 

%% part 4: analyze whether the errors w_c are zero mean 

w_c = yc - yc_true'; 
disp('w_c is not zero-mean') 

mean_wc = mean(w_c); 

end 