% Implement a Kalman filter for a stochastic linear time invariant (SLTI) 
% system in the standard form used in class (with Γ(k) =/= I). 

% The problem matrices and the measurement data, z(k) for k = 1, ..., 50, 
% can be loaded into your Matlab workspace by running the Matlab script 
% kf_example02a.m. 

% Hand in plots of the two elements of xˆ(k) vs. time and of the predicted 
% standard deviations of xˆ(k) vs. time, i.e., of sqrt([P(k)_11]) and 
% sqrt([P(k)_22]).

% Plot each element of xˆ(k) and its corresponding standard deviation 
% together on the same graph. 

% Use symbols on the plot at each of the 51 points and do not connect the 
% symbols by lines (type “help plot” in order to learn how to do this). 
% Also, hand in numerical values for the terminal values of xˆ(50) and P(50).

clear; clc; 
% close all 

kf_example02a; 
Gk = Gammak; 

%% KALMAN FILTER 

% initialize for k = 0 
xhat = xhat0; 
P    = P0; 

% Initialize saved output arrays 
xbar_arr = []; 
Pbar_arr = []; 
xhat_arr = [xhat']; 
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
    Pxx_arr  = [Pxx_arr; P(1,1)]; 
    Pzz_arr  = [Pzz_arr; P(2,2)]; 

end 

%% plot 

ftitle = 'States and Covariances'; 
figure('name', ftitle); 
    subplot(2,1,1) 
        plot( xhat_arr(:,1), '.' ); hold on; grid on; 
        plot( sqrt( Pxx_arr ), '.'); 
        title('$\hat{x}$(1)', 'interpreter', 'latex'); 
        legend('$\hat{x}$', '$\sqrt{(P_{xx})}$', 'interpreter', 'latex'); 
    subplot(2,1,2) 
        plot( xhat_arr(:,2), '.' ); hold on; grid on; 
        plot( sqrt( Pzz_arr ), '.' ); 
        title('$\hat{x}$(2)', 'interpreter', 'latex'); 
        legend('$\hat{x}$', '$\sqrt{(P_{zz})}$', 'interpreter', 'latex'); 
    sgtitle(ftitle); 







