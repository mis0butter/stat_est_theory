close all; 
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

x_hat = x0; 
P     = P0; 

% for lambda 
a = 10^-3; 
b = 2; 
k = 0; 
nx = 6; 
nz = 3; 
nv = 6; 

% load LIDAR data 
load problem3data.mat 
N = length(lidar); 

% initialize output arrays 
x_hat_arr = [x_hat']; 
P_cell = {P}; 

% create plot 
fname = 'Car Plot'; 
h = figure('name', fname); 
    xlim([0 100]); ylim([-50 50]); 
    plotcar(x_hat, '-', h)
%     pause(0.05)
%     drawnow 

%% UNSCENTED KALMAN FILTER 

% iterate through measurements 
for j = 1 : N 

% augmented state and covariance 
xa = [ x_hat; zeros(6,1) ]; 
Pa = [ P, zeros(6); zeros(6), Q ]; 

%% DYNAMICS PROPAGATION 

% cholesky factorize Pa 
Sx = chol(Pa)'; 

% build sigma points 
lambda_xv = a^2 * (nx + nv + k) - (nx + nv); 
XX = build_SP(xa, Pa, nx, nv, lambda_xv); 

% propagate sigma points 
XX_prop = []; 
for i = 1:length(XX)
    
    [tvec, XX_i_prop] = ode45(@dyn_car, [0 dt], XX(i,:)); 
    XX_i_prop = XX_i_prop(end,:); 
    XX_prop   = [ XX_prop; XX_i_prop ]; 
    
end 

% combine sigma points 
[x_bar, P_bar] = combine_SP(nx, nv, lambda_xv, a, b, XX_prop); 

%% MEASUREMENT UPDATE 

% push sigma points through measurement model. Use min bearing, max
% bearing, and min range 
zj   = lidar(j).z; 
z_rj = zj(:,1);     % range 
z_bj = zj(:,2);     % bearing 
z    = [ min(z_bj); max(z_bj); min(z_rj) ]; 

% stack states and covariance 
xa_bar = [x_bar; zeros(nz,1)]; 
Pa_bar = [P_bar, zeros(nx, nz); zeros(nz, nx), R]; 

% build sigma points 
lambda_xz = a^2 * (nx + nz + k) - (nx + nz); 
XX_bar = build_SP(xa_bar, Pa_bar, nz, nx, lambda_xz); 

% push sigma points through measurement model 
ZZ_bar = []; 
for i = 1:length(XX_bar)
    ZZ_i_bar = h_car(XX_bar(i,:)); 
    ZZ_bar   = [ ZZ_bar; ZZ_i_bar' ]; 
end 

% combine sigma points 
[z_bar, Pzz] = combine_SP(nz, nx, lambda_xz, a, b, ZZ_bar); 

% calculate Pxz 
Pxz = calc_Pxz(nx, nz, lambda_xz, a, b, XX_bar, x_bar, ZZ_bar, z_bar); 

% LMMSE update!!! 
x_hat = x_bar + Pxz * Pzz^(-1) * [z - z_bar]; 
P     = P_bar - Pxz * Pzz^(-1) * Pxz'; 

% save state and covariance
x_hat_arr = [x_hat_arr; x_hat']; 
P_cell{j+1} = {P}; 

% update plot 
plotcar(x_hat, '-', h)

end 

%% ANALYSIS 

load problem3truth.mat


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

function Pxz = calc_Pxz(nx, nz, lambda_xz, a, b, XX_bar, x_bar, ZZ_bar, z_bar)

    % determine weights 
    w_0m = lambda_xz / (nx + nz + lambda_xz); 
    w_im = 1 / (2*(nx + nz + lambda_xz)); 
    w_0c = lambda_xz / (nx + nz + lambda_xz) + 1 - a^2 + b; 
    w_ic = w_im; 

    Pxz = zeros(size(nx, nz)); 
    N_SP = 2*(nx + nz) + 1; 
    for i = 1 : N_SP

        if i == 1;      wP = w_0c; 
        else;           wP = w_ic; 
        end

        % build Pa_bar 
        xtilde = [ XX_bar(i, 1:nx)' - x_bar ]; 
        ztilde = [ ZZ_bar(i, 1:nz)' - z_bar ]; 

        Pi  = wP * (xtilde) * (ztilde)'; 
        Pxz = Pxz + Pi; 

    end 

end 

function [x_bar, Pxx_bar] = combine_SP(nx, nv, lambda, a, b, XX_prop) 

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
    Pxx_bar = zeros(nx); 
    for i = 1 : N_SP

        if i == 1;      wP = w_0c; 
        else;           wP = w_ic; 
        end

        % build Pa_bar 
        xtilde = [ XX_prop(i, 1:nx) - x_bar ]'; 
        Pi_bar = wP * (xtilde) * (xtilde)'; 
        Pxx_bar = Pxx_bar + Pi_bar; 

    end 
    
    x_bar = x_bar';

end 

