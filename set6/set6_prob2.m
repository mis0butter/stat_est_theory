%% problem set-up 

% There are two Kalman filtering problems defined by the two MATLAB scripts
% kf example03a.m and kf example03b.m. Solve each of these Kalman filtering 
% problems in two ways. First, use a standard Kalman filter. Second, use a 
% square-root information filter (SRIF). The two filters should work about 
% equally well for the problem in kf example03a.m, but the Kalman filter 
% will not work as well as the SRIF for the problem in kf example03b.m. 
% This is true because of the small R(k) value, which causes the computed 
% covariance matrix to be ill-conditioned. Compare the state estimates and 
% the covariances for the two filters at the terminal time. Is there a 
% significant difference for the second filtering problem but not for the 
% first? Note: this is a relatively benign case. The improvement due to use 
% of the SRIF is not extremely significant. There are, however, known 
% practical situations where a Kalman filter completely breaks down while 
% an SRIF functions well.

clear; clc; close all 

%% KALMAN FILTER AND SQUARE-ROOT INFORMATION FILTER: example 03a 

clear; 
disp('EXAMPLE 03A') 
kf_example03a; 

[xhat_arr_kf, Pxx_arr, Pzz_arr, P_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk ); 
% plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)

uk = [0; 0]; 
nx = length(xhat0); 
nv = length(Qk); 
[xhat_arr_srif, Rxx_arr] = srif( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk, uk, nx, nv); 

xhat_kf   = xhat_arr_kf(end,:); 
xhat_srif = xhat_arr_srif(end,:); 
P_kf      = P_cell{end}; 
Rxx       = Rxx_arr{end}; 
P_srif    = inv(Rxx) * inv(Rxx)'; 

disp('KF terminal estimate: ')
xhat_kf 

disp('SRIF terminal estimate: ') 
xhat_srif 

disp('KF terminal covariance: ') 
P_kf 

disp('SRIF terminal covariance: ') 
P_srif 



%% KALMAN FILTER AND SQURE-ROOT INFORMATION FILTER: example 03b 

clear; 
disp('EXAMPLE 03B') 
kf_example03b; 

[xhat_arr_kf, Pxx_arr, Pzz_arr, P_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk ); 
% plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)

uk = [0; 0]; 
nx = length(xhat0); 
nv = length(Qk); 
[xhat_arr_srif, Rxx_arr] = srif( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk, uk, nx, nv); 

xhat_kf   = xhat_arr_kf(end,:); 
xhat_srif = xhat_arr_srif(end,:); 
P_kf      = P_cell{end}; 
Rxx       = Rxx_arr{end}; 
P_srif    = inv(Rxx) * inv(Rxx)'; 

disp('KF terminal estimate: ')
xhat_kf 

disp('SRIF terminal estimate: ') 
xhat_srif 

disp('KF terminal covariance: ') 
P_kf 

disp('SRIF terminal covariance: ') 
P_srif 

%% results 

disp('There is a significant difference for the state estimates for the second person, but not the first') 

%% subfunctions 

function [xhat_arr, Rxx_cell] = srif( ... 
    xhat0, P0, zhist, Fk, Gk, Qk, Hk, Rk, uk, nx, nv)

% Initialize 
xhat_arr = []; 
Rxx_cell = {}; 

% START AT k = 0: 
I = inv(P0); 
Rxx = chol(I); 
Rvv = chol(inv(Qk)); 
zv = zeros(nv, 1); 
zx = Rxx * xhat0; 

for i = 1:length(zhist)

    % PROPAGATION STEP 
    % a) QR factorize 
    A = [Rvv, zeros(2,3); -Rxx * inv(Fk) * Gk, Rxx * inv(Fk)]; 
    [QA, RA] = qr(A); 

    % b) orthonormal transformation 
    B = [ zv; zx + Rxx * inv(Fk) * Gk * uk ];  
    [zv_zx_bar] = QA' * B; 

    % c) extract Rxx_bar(k+1) and zx_bar(k+1) 
    zv_bar = zv_zx_bar(1:nv); 
    zx_bar = zv_zx_bar(nv+1:end); 
    Rxx_bar = RA(nv+1:end, nv+1:end); 

    % MEASUREMENT UPDATE: 
    % a) Cholesky factorize R 
    Ra = chol(Rk); 

    % b) Transform z(k+1) and H(k+1) 
    za = inv(Ra)' * zhist(i); 
    Ha = inv(Ra)' * Hk; 

    % c) perform another QR factorization: 
    [QB, RB] = qr([Rxx_bar; Ha]); 

    % d) transform as 
    [zx_zr] = QB' * [zx_bar; za];    

    % e) extract Rxx_bar(k+1) and zx(k+1) 
    zx = zx_zr(1:nx); 
    Rxx = RB(1:nx, :); 
    
    xhat = inv(Rxx) * zx; 
    xhat_arr = [xhat_arr; xhat']; 
    Rxx_cell{i} = Rxx;  
    
end 

end 

function plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)
% plot xhat and Pxx 

thist0 = [ 0; thist ]; 

% plot 
ftitle = 'States and Covariances'; 
figure('name', ftitle); 
    subplot(2,1,1) 
        plot( thist0, xhat_arr(:,1), '.' ); hold on; grid on; 
        plot( thist0, xhat_arr(:,1) + sqrt( Pxx_arr ), 'r--'); 
        plot( thist0, xhat_arr(:,1) - sqrt( Pxx_arr ), 'r--'); 
        title('$\hat{x}$(1)', 'interpreter', 'latex'); 
        legend('$\hat{x}$', '$ \hat{x} \pm \sigma_{xx}$', 'interpreter', 'latex', 'location', 'best'); 
        ylabel('state units'); 
    subplot(2,1,2) 
        plot( thist0, xhat_arr(:,2), '.' ); hold on; grid on; 
        plot( thist0, xhat_arr(:,2) + sqrt( Pzz_arr ), 'r--'); 
        plot( thist0, xhat_arr(:,2) - sqrt( Pzz_arr ), 'r--'); 
        title('$\hat{x}$(2)', 'interpreter', 'latex'); 
        legend('$\hat{x}$', '$ \hat{x} \pm \sigma_{zz}$', 'interpreter', 'latex', 'location', 'best'); 
        ylabel('state units'); 
    xlabel('time'); 
    sgtitle(ftitle); 

% print final values 
disp('xhat(50) =')
disp(xhat_arr(end,:))

disp('P(50) =')
disp(P_cell{end})

end 

function [xhat_arr, Pxx_arr, Pzz_arr, P_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gk, Qk, Hk, Rk )

% initialize for k = 0 
xhat = xhat0; 
P    = P0; 

% Initialize saved output arrays 
xbar_arr = []; 
Pbar_arr = []; 
xhat_arr = [xhat']; 
P_cell   = {P}; 
Pxx_arr  = [P(1,1)]; 
Pzz_arr  = [P(2,2)]; 

% Propagate and filter through all measurements 
for k = 0 : length(zhist)-1

    % propagate state and covar 
    xbar = Fk * xhat;                       % a priori state est 
    Pbar = Fk * P * Fk' + Gk * Qk * Gk';    % a posteriori covar est 

    % update 
    v = zhist(k+1) - Hk * xbar;             % innovation 
    S = Hk * Pbar * Hk' + Rk;               % innovation covariance 
    W = Pbar * Hk' * inv(S);                % Kalman gain 
    xhat = xbar + W * v;                    % a posteriori state est 
    P = Pbar - W * S * W';                  % a posteriori covar est 
    
    % next step 
    k = k + 1; 
    
    % save states and covariances 
    xbar_arr = [xbar_arr; xbar']; 
    Pbar_arr = [Pbar_arr; Pbar]; 
    xhat_arr = [xhat_arr; xhat']; 
    P_cell   = {P_cell; P}; 
    Pxx_arr  = [Pxx_arr; P(1,1)]; 
    Pzz_arr  = [Pzz_arr; P(2,2)]; 

end 

end 








