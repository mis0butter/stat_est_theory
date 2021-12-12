clear
clc 

%%

syms r_bar w_r theta_bar w_theta
syms sigma_r sigma_theta 

% polar measurements 
r = r_bar + w_r; 
theta = theta_bar + w_theta; 

% polar covariance 
R_i = diag([ sigma_r^2 sigma_theta^2 ]); 
R_fn = matlabFunction(R_i); 

% cartesian coordinates 
x = r * cos(theta); 
y = r * sin(theta); 

% expected cartesian 
x_bar = r_bar * cos(theta_bar); 
y_bar = r_bar * sin(theta_bar); 
x_bar_fn = matlabFunction(x_bar); 
y_bar_fn = matlabFunction(y_bar); 

% linearized cartesian about expected values 
dx_bar = jacobian(x_bar, [r_bar theta_bar]); 
dy_bar = jacobian(y_bar, [r_bar theta_bar]); 
dx_bar_fn = matlabFunction(dx_bar); 
dy_bar_fn = matlabFunction(dy_bar); 

% linearized cartesian covariance 
Rc = [dx_bar; dy_bar] * R_i * [dx_bar; dy_bar].'; 
Rc_fn = matlabFunction(Rc); 

% inputs part i 
r_bar_val_i = 76; 
theta_bar_val_i = -3*pi/180; 
sigma_r_val_i = 1; 
sigma_theta_val_i = pi/180; 

% compute 
H = [ dx_bar_fn(r_bar_val_i, theta_bar_val_i); 
      dy_bar_fn(r_bar_val_i, theta_bar_val_i) ]; 
dr = [ r_bar_val_i - r_bar_val_i ]; 
dtheta = [ theta_bar_val_i - theta_bar_val_i ]; 
xy_bar = [ x_bar_fn(r_bar_val_i, theta_bar_val_i); y_bar_fn(r_bar_val_i, theta_bar_val_i) ]; 
z_lin_i  = xy_bar + H * [dr; dtheta]; 
Rc_val_i = Rc_fn(r_bar_val_i, sigma_r_val_i, sigma_theta_val_i, theta_bar_val_i); 

% inputs part ii 
r_bar_val_ii = 76; 
theta_bar_val_ii = -3*pi/180; 
sigma_r_val_ii = 1; 
sigma_theta_val_ii = 15*pi/180; 

% compute 
H = [ dx_bar_fn(r_bar_val_ii, theta_bar_val_ii); 
      dy_bar_fn(r_bar_val_ii, theta_bar_val_ii) ]; 
dr = [ r_bar_val_ii - r_bar_val_ii ]; 
dtheta = [ theta_bar_val_ii - theta_bar_val_ii ]; 
xy_bar = [ x_bar_fn(r_bar_val_ii, theta_bar_val_ii); y_bar_fn(r_bar_val_ii, theta_bar_val_ii) ]; 
z_lin_ii  = xy_bar + H * [dr; dtheta]; 
Rc_val_ii = Rc_fn(r_bar_val_ii, sigma_r_val_ii, sigma_theta_val_ii, theta_bar_val_ii); 

%% part b: unscented transform 

% cholesky factorize 
R_i = R_fn(sigma_r_val_i, sigma_theta_val_i); 

[z_bar_i, Pzz_i] = unscented_transform(r_bar_val_i, theta_bar_val_i, R_i); 

% cholesky factorize 
R_ii = R_fn(sigma_r_val_ii, sigma_theta_val_ii); 

[z_bar_ii, Pzz_ii] = unscented_transform(r_bar_val_ii, theta_bar_val_ii, R_ii); 

%% part c: large random vectors 

N = 100000; 
w = mvnrnd([r_bar_val_i; theta_bar_val_i], R_i, N); 

z = [ w(:,1) .* cos(w(:,2)), w(:,1) .* sin(w(:,2)) ];  
disp('Part i inputs:') 
fprintf('Mean of %d samples: x = %g, y = %g \n', N, mean(z(:,1)), mean(z(:,2)));
fprintf('Linearized mean: x = %g, y = %g \n', z_lin_i(1), z_lin_i(2)); 
fprintf('UT mean: x = %g, y = %g \n\n', z_bar_i(1), z_bar_i(2)); 

N = 100000; 
w = mvnrnd([r_bar_val_ii; theta_bar_val_ii], R_ii, N); 

z = [ w(:,1) .* cos(w(:,2)), w(:,1) .* sin(w(:,2)) ];  

disp('Part ii inputs: ') 
fprintf('Mean of %d samples: x = %g, y = %g \n', N, mean(z(:,1)), mean(z(:,2))); 
fprintf('Linearized mean: x = %g, y = %g \n', z_lin_ii(1), z_lin_ii(2)); 
fprintf('UT mean: x = %g, y = %g \n', z_bar_ii(1), z_bar_ii(2)); 


%% subfunctions 

function [z_bar, Pzz] = unscented_transform(r_bar_val, theta_bar_val, R) 
s = chol(R)'; 

% obtain lambda 
a = 10^-3; 
b = 2; 
k = 0; 
nx = 2; 
nv = 0; 
nz = 2; 
% lambda = a^2 * (nx + 1/2) - nx; 
lambda = a^2 * (nx + nv + k) - (nx + nv); 

% build sigma points (still polar coordinates) 
rtheta_bar = [r_bar_val; theta_bar_val]; 
R = rtheta_bar' ; 
for i = 1 : nx
    Ri = rtheta_bar + sqrt( nx + lambda ) * s(:,i); 
    R = [ R; Ri' ]; 
end 
for i = nx + 1 : 2*nx
    Ri = rtheta_bar - sqrt( nx + lambda) * s(:, i-nx); 
    R = [ R; Ri' ]; 
end 

% push sigma points through measurement model --> 
z = @(sigma) [sigma(1) * cos(sigma(2)); sigma(1) * sin(sigma(2))]; 
for i = 1 : 2*nx + 1
    Z(i,:) = z(R(i,:)); 
end 

% determine weighting 
w_0m = lambda / (nx + lambda); 
w_im = 1 / ( 2*(nx + lambda) ); 
w_0c = lambda / (nx + lambda) + 1 - a^2 - b; 
w_ic = 1 / ( 2*(nx + lambda) ); 

% recombine with the weighting 
z_bar = w_0m * Z(1,:); 
for i = 2 : 2*nx+1
    z_bar = z_bar + w_ic * Z(i,:); 
end  

% covariance 
ztilde = Z - z_bar; ztilde = ztilde'; 
Pzz = zeros(size(nz)); 
for i = 1:length(ztilde) 
    if i == 1 
        w = w_0c; 
    else
        w = w_ic; 
    end 
    Pzz = Pzz + w * ztilde(:,i) * ztilde(:,i)'; 
end 

end 













