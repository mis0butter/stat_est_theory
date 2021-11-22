%% problem set-up 

% There are two Kalman filtering problems defined by the two MATLAB scripts
% kf example03a.m and kf example03b.m. Solve each of these Kalman filtering problems
% in two ways. First, use a standard Kalman filter. Second, use a square-root information filter (SRIF). The two filters should work about equally well for the problem in
% kf example03a.m, but the Kalman filter will not work as well as the SRIF for the problem in kf example03b.m. This is true because of the small R(k) value, which causes the
% computed covariance matrix to be ill-conditioned. Compare the state estimates and the
% covariances for the two filters at the terminal time. Is there a significant difference for
% the second filtering problem but not for the first?
% Note: this is a relatively benign case. The improvement due to use of the SRIF is not
% extremely significant. There are, however, known practical situations where a Kalman
% filter completely breaks down while an SRIF functions well.

clear; clc; close all 

%% KALMAN FILTER 


kf_example03a; 
[xhat_arr, Pxx_arr, Pzz_arr, P_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk ); 
plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)

clear; 
kf_example03b; 
[xhat_arr, Pxx_arr, Pzz_arr, P_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk ); 
plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)

%% results 

%% subfunctions KALMAN FILTER 

function plot_xhat_Pxx(thist, xhat_arr, Pxx_arr, Pzz_arr, P_cell)

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








