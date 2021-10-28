% Repeat Problem 3, except use the problem matrices and measurement data 
% that are defined by the Matlab script kf example02b.m. Notice that the R 
% and Q values are different for this problem and that there is a different 
% measurement time history. 

% Run your Kalman filter two additional times using the two alternate Q 
% values that are mentioned in the comments in the file kf example02b.m. 
% It is uncertain which is the correct Q value.

% Decide which is the best value in the following way: Calculate err(ν(k)) 
% for k = 1, 2, ..., 50 for each of your runs. Compute the average of these 
% 50 values. This average times 50, i.e., {err_ν(1) + err_ν(2) + ... + 
% err_ν(50)}, will be a sample of chi-square distribution of degree 50 if 
% the filter model is correct. Develop upper and lower limits between which 
% the average {err_ν(1) + err_ν(2) + ... + err_ν(50)}/50 must lie 99% of 
% the time if the Kalman filter model is correct, and test your averages 
% for each of the three candidate Q values. Which is the most reasonable? 
% Look at the state estimate differences between the best filter and the
% other two filters. Compute the RMS value of the difference time history 
% for each state vector element. Do the averaging over the last 40 points. 
% Are these differences significant compared to the computed state 
% estimation error standard deviations for the best filter?

clear; clc; 
close all;  

kf_example02b

% Qk_a = Qk 
% Qk_b = first alternate Qk 
% Qk_c = second alternate Qk 
Gk = Gammak; 

%% KALMAN FILTER 

% Qk_a 
[xhat_a, Pxx_a, Pzz_a, P_a, xbar_a, Pbar_a] = kf( ... 
    xhat0, P0, zhist, Fk, Gk, Qk_a, Hk, Rk ); 
e_a = xbar_a - xhat_a; 
e_a_mean = mean(e_a); 

% Qk_b 
[xhat_b, Pxx_b, Pzz_b, P_b, xbar_b, Pbar_b] = kf( ... 
    xhat0, P0, zhist, Fk, Gk, Qk_b, Hk, Rk ); 


% Qk_c 
[xhat_c, Pxx_c, Pzz_c, P_c, xbar_c, Pbar_c] = kf( ... 
    xhat0, P0, zhist, Fk, Gk, Qk_c, Hk, Rk ); 



%% subfunctions KALMAN FILTER 

function [xhat_arr, Pxx_arr, Pzz_arr, P_cell, xbar_arr, Pbar_cell] = kf( ... 
    xhat0, P0, zhist, Fk, Gk, Qk, Hk, Rk )

% initialize for k = 0 
xhat = xhat0; 
P    = P0; 

% Initialize saved output arrays 
xbar_arr  = [xhat']; 
Pbar_cell = {P}; 
xhat_arr  = [xhat']; 
P_cell    = {P}; 
Pxx_arr   = [P(1,1)]; 
Pzz_arr   = [P(2,2)]; 

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
    Pbar_cell = {Pbar_cell; Pbar}; 
    xhat_arr = [xhat_arr; xhat']; 
    P_cell   = {P_cell; P}; 
    Pxx_arr  = [Pxx_arr; P(1,1)]; 
    Pzz_arr  = [Pzz_arr; P(2,2)]; 

end 

end 



