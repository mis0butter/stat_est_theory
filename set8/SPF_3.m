
clear; 
% clc

% sim dt 
dt = 0.1; 

% process noise covariance 
Q = diag([0.25 0.25 3 40*pi/180 0.1 0.1])^2 / dt; 

% measurement vector:
% z = [b_min b_max r_min]'; 
%   b_min = minimum bearing (rad)
%   b_max = maximum bearing (rad) 
%   r_min = minimum range (m) 

% measurement covariance 
R = diag( [2*pi/180, 2*pi/180, 0.1] )^2; 

% initial state est and covar 
x0 = [90; 4.25; 13; pi; 5; 2]; 
P0 = diag([2 5 1 pi/4 4 2])^2; 

xa = [ x0; zeros(6,1) ]; 
Pa = [ P0, zeros(6); zeros(6), Q ]; 

% for lambda 
a = 10^-3; 
b = 2; 
k = 0; 
nx = 6; 
nz = 3; 
nv = 6; 

% load LIDAR data 
load problem3data.mat 

%% DYNAMICS PROPAGATION 

% cholesky factorize Pa 
Sx = chol(Pa)'; 

% build sigma points 
lambda_x = a^2 * (nx + nv + k) - (nx + nv); 
XX = build_SP(xa, Pa, nx, nv, lambda_x); 

% propagate sigma points 
XX_prop = []; 
for i = 1:length(XX)
    
    [tvec, XX_i_prop] = ode45(@dyn_car, [0 dt], XX(i,:)); 
    XX_i_prop = XX_i_prop(end,:); 
    XX_prop   = [ XX_prop; XX_i_prop ]; 
    
end 

% combine sigma points 
[x_bar, P_bar] = combine_SP(nx, nv, lambda_x, a, b, XX_prop); 

%     % determine weights 
%     w_0m = lambda_x / (nx + nv + lambda_x); 
%     w_im = 1 / (2*(nx + nv + lambda_x)); 
%     w_0c = lambda_x / (nx + nv + lambda_x) + 1 - a^2 + b; 
%     w_ic = w_im; 
%     
%     N_SP = 2*(nx+nv) + 1; 
% 
%     % use predicted sigma points to calculate predicted state x_bar 
%     x_bar = zeros(1, nx); 
%     for i = 1 : N_SP
% 
%         if i == 1;      wx = w_0m; 
%         else;           wx = w_im; 
%         end
% 
%         % build xa_bar 
%         xi_bar = wx * XX_prop(i, 1:nx); 
%         x_bar = x_bar + xi_bar; 
% 
%     end  
% 
%     % predict Pa_bar 
%     P_bar = zeros(nx); 
%     for i = 1 : N_SP
% 
%         if i == 1;      wP = w_0c; 
%         else;           wP = w_ic; 
%         end
% 
%         % build Pa_bar 
%         xtilde = [ XX_prop(i, 1:nx) - x_bar ]'; 
%         Pi_bar = wP * (xtilde) * (xtilde)'; 
%         P_bar = P_bar + Pi_bar; 
% 
%     end 

%% MEASUREMENT UPDATE 

% push sigma points through measurement model. Use min bearing, max
% bearing, and min range 
i    = 1; 
zi   = lidar(i).z; 
z_ri = zi(:,1);     % range 
z_bi = zi(:,2);     % bearing 
z    = [ min(z_bi); max(z_bi); min(z_ri) ]; 

% stack states and covariance 
xa_bar = [x_bar; zeros(nz,1)]; 
Pa_bar = [P_bar, zeros(nx, nz); zeros(nz, nx), R]; 

% build sigma points 
lambda_z = a^2 * (nx + nz + k) - (nx + nz); 
% XX_bar = build_SP(xa_bar, Pa_bar, nz, nx, lambda_z); 

    % cholesky factorize Pa 
    Sx_bar = chol(Pa_bar)'; 

    % build sigma points. REMEMBER: there will be 2*(nx+nv)+1 sigma points 
    XX_bar(1,:) = xa_bar'; 
    for i = 1 : nx + nz 

        XXi = xa_bar' + sqrt( nx + nz + lambda_z ) * Sx_bar(:,i)'; 
        XX_bar  = [ XX_bar; XXi ]; 

    end 
    for i = nx + nz + 1 : 2*(nx + nz)

        XXi = xa_bar' - sqrt( nx + nz + lambda_z) * Sx_bar(:, i - nx - nz)'; 
        XX_bar  = [ XX_bar; XXi ]; 

    end 

% push sigma points through measurement model 
ZZ_bar = []; 
for i = 1:length(XX_bar)
    ZZ_i_bar = h_car(XX_bar(i,:)); 
    ZZ_bar   = [ ZZ_bar; ZZ_i_bar' ]; 
end 

% determine weights 
w_0m = lambda_z / (nx + nz + lambda_z); 
w_im = 1 / (2*(nx + nz + lambda_z)); 
w_0c = lambda_z / (nx + nz + lambda_z) + 1 - a^2 + b; 
w_ic = w_im; 

N_SP = 2*(nx+nz) + 1; 

% use predicted sigma points to calculate predicted state x_bar 
z_bar = zeros(1, nz); 
for i = 1 : N_SP

    if i == 1;      wz = w_0m; 
    else;           wz = w_im; 
    end

    % build xa_bar 
    zi_bar = wz * ZZ_bar(i, 1:nz); 
    z_bar = z_bar + zi_bar; 

end 
z_bar = z_bar'; 

% combine sigma points 
[z_bar, Pzz_bar] = combine_SP(nz, nx, lambda_z, a, b, ZZ_bar); 




%% subfunctions 

function XX = build_SP(xa, Pa, nx, nv, lambda) 

    % cholesky factorize Pa 
    Sx = chol(Pa)'; 

    % build sigma points. REMEMBER: there will be 2*(nx+nv)+1 sigma points 
    XX(1,:) = xa'; 
    for i = 1 : nx + nv 

        XXi = xa' + sqrt( nx + nv + lambda ) * Sx(:,i)'; 
        XX  = [ XX; XXi ]; 

    end 
    for i = nx + nv + 1 : 2*(nx + nv)

        XXi = xa' - sqrt( nx + nv + lambda) * Sx(:, i - nx - nv)'; 
        XX  = [ XX; XXi ]; 

    end 

end 

function [x_bar, P_bar] = combine_SP(nx, nv, lambda, a, b, XX_prop) 

    % determine weights 
    w_0m = lambda / (nx + nv + lambda); 
    w_im = 1 / (2*(nx + nv + lambda)); 
    w_0c = lambda / (nx + nv + lambda) + 1 - a^2 + b; 
    w_ic = w_im; 
    
    N_SP = 2*(nx+nv) + 1; 

    % use predicted sigma points to calculate predicted state x_bar 
    x_bar = zeros(1, nx); 
    for i = 1 : N_SP

        if i == 1;      wx = w_0m; 
        else;           wx = w_im; 
        end

        % build xa_bar 
        xi_bar = wx * XX_prop(i, 1:nx); 
        x_bar = x_bar + xi_bar; 

    end  

    % predict Pa_bar 
    P_bar = zeros(nx); 
    for i = 1 : N_SP

        if i == 1;      wP = w_0c; 
        else;           wP = w_ic; 
        end

        % build Pa_bar 
        xtilde = [ XX_prop(i, 1:nx) - x_bar ]'; 
        Pi_bar = wP * (xtilde) * (xtilde)'; 
        P_bar = P_bar + Pi_bar; 

    end 
    
    x_bar = x_bar';

end 

