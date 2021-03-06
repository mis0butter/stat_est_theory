% This problem refers back to Problems 5 & 6 of Problem Set 5. Re-do Problem 6 using the
% Kalman filter problem matrices of Problem 5. Note that in Problem 5 there are 3 possible
% Q values. Run your truth-model simulation using the largest of the three Q values, but
% run your Kalman filter using the smallest of the three Q values. What does this do to
% your consistency evaluation? Is this what you expect?

%% problem set-up 

clear; clc 
close all 

kf_example02b

% Qk_a = Qk 
% Qk_b = first alternate Qk  
% Qk_c = second alternate Qk 
Gk = Gammak; 


n = 10000; 
Nx = 2; 
Nz = 1; 

% NEED STATISTICS TOOLBOX 
a = .01; 
r1 = chi2inv( a/2, n * Nx) / n; 
r2 = chi2inv( 1 - a/2, n * Nx ) / n; 

%% LARGE Q MC, SMALL Q KF 

rng(0) 

% monte carlo and KF 
txt = 'MC and KF using large Q for MC and small Q for KF:'; 
disp(txt); 
[xtilde10_vec, err10_vec, xtilde35_vec, err35_vec] = ... 
    mcltisim_kf(Fk, Gammak, Hk, Qk_a, Qk_c, Rk, xhat0, P0, n); 

xtilde10_cumsum = cumsum(xtilde10_vec); 
xtilde35_cumsum = cumsum(xtilde35_vec); 
for k = 1:length(xtilde10_vec) 
    xtilde10_cummean(k) = xtilde10_cumsum(k) / k; 
    xtilde35_cummean(k) = xtilde35_cumsum(k) / k; 
end 

hf = figure(); 
    subplot(2,1,1) 
        plot(xtilde10_cummean); hold on; grid on; 
        plot(xtilde35_cummean, '--'); 
        legend('xtilde10', 'xtilde35'); 
        title(txt); 
        
mean_err10 = mean(err10_vec); 
mean_err35 = mean(err35_vec); 
if r1 < mean_err10 && mean_err10 < r2 
    disp('nu10 is consistent!') 
else 
    disp('nu10 is not consistent');  
end 
if r1 < mean_err35 && mean_err35 < r2 
    disp('nu35 is consistent!') 
else 
    disp('nu35 is not consistent');  
end 

        
%% LARGE Q MC, LARGE Q KF 

rng(0)

% monte carlo and KF 
txt = 'MC and KF using large Q for MC and large Q for KF:';  
disp(txt); 
[xtilde10_vec, err10_vec, xtilde35_vec, err35_vec] = ... 
    mcltisim_kf(Fk, Gammak, Hk, Qk_a, Qk_a, Rk, xhat0, P0, n); 

xtilde10_cumsum = cumsum(xtilde10_vec); 
xtilde35_cumsum = cumsum(xtilde35_vec); 
for k = 1:length(xtilde10_vec) 
    xtilde10_cummean(k) = xtilde10_cumsum(k) / k; 
    xtilde35_cummean(k) = xtilde35_cumsum(k) / k; 
end 

figure(hf); 
    subplot(2,1,2) 
        plot(xtilde10_cummean); hold on; grid on; 
        plot(xtilde35_cummean, '--'); 
        legend('xtilde10', 'xtilde35'); 
        title(txt); 
        
mean_err10 = mean(err10_vec); 
mean_err35 = mean(err35_vec); 
if r1 < mean_err10 && mean_err10 < r2 
    disp('\nu10 is consistent!') 
else 
    disp('\nu10 is not consistent');  
end 
if r1 < mean_err35 && mean_err35 < r2 
    disp('\nu35 is consistent!') 
else 
    disp('\nu35 is not consistent');  
end 

%% The 




%% subfunctions 

function [xtilde10_vec, err10_vec, xtilde35_vec, err35_vec] = ... 
    mcltisim_kf(Fk, Gammak, Hk, Qk_a, Qk_c, Rk, xhat0, P0, n)

xtilde10_vec = []; 
xtilde35_vec = []; 
xtilde_vec   = []; 
err10_vec    = []; 
err35_vec    = []; 

for i_mc = 1 : n % run 50 or 10000 monte carlos 

    kmax = 50; % zhist --> 10 or 35 
    [xhist, zhist] = mcltisim(Fk, Gammak, Hk, Qk_a, Rk, xhat0, P0, kmax); 

    % Pass into kalman filter 
    [xhat, Pxx, Pzz, P_cell] = kf( ... 
        xhat0, P0, zhist, Fk, Gammak, Qk_c, Hk, Rk );

    % xtilde = truth - estimate (xhat0 in 1st element) 
    xtilde10_vec = [xtilde10_vec; xhist(11,:) - xhat(11,:)]; 
    xtilde35_vec = [xtilde35_vec; xhist(36,:) - xhat(36,:)]; 
    xtilde_vec = [xtilde_vec; xhist - xhat]; 
    
    xtilde10 = xhist(11,:) - xhat(11,:); 
    err10 = xtilde10 * inv(P_cell{11}) * xtilde10'; 
    err10_vec = [err10_vec; err10]; 
    
    xtilde35 = xhist(36,:) - xhat(36,:); 
    err35 = xtilde35 * inv(P_cell{36}) * xtilde35'; 
    err35_vec = [err35_vec; err35]; 
    
    if i_mc == 50
        disp('i = 50') 
        sprintf('xtilde10 x1 mean = %g, x2 mean = %g', mean(xtilde10_vec(:,1)), mean(xtilde10_vec(:,2)) )
        sprintf('xtilde10 x1 cov = %g, x2 cov = %g', cov(xtilde10_vec(:,1)), cov(xtilde10_vec(:,2)) )
        sprintf('xtilde35 x1 mean = %g, x2 mean = %g', mean(xtilde35_vec(:,1)), mean(xtilde35_vec(:,2)) )
        sprintf('xtilde35 x1 cov = %g, x2 cov = %g', cov(xtilde35_vec(:,1)), cov(xtilde35_vec(:,2)) )
    elseif i_mc == n
        sprintf('i = %d', n) 
        sprintf('xtilde10 x1 mean = %g, x2 mean = %g', mean(xtilde10_vec(:,1)), mean(xtilde10_vec(:,2)) )
        sprintf('xtilde10 x1 cov = %g, x2 cov = %g', cov(xtilde10_vec(:,1)), cov(xtilde10_vec(:,2)) )
        sprintf('xtilde35 x1 mean = %g, x2 mean = %g', mean(xtilde35_vec(:,1)), mean(xtilde35_vec(:,2)) )
        sprintf('xtilde35 x1 cov = %g, x2 cov = %g', cov(xtilde35_vec(:,1)), cov(xtilde35_vec(:,2)) )
    end 

end 

end 

function [xhist, zhist] = mcltisim(F, Gamma, H, Q, R, xbar0, P0, kmax)
% ----------------------------------------------------------------------- %
% ltisim : Monte-Carlo simulation of a linear time invariant system.
%
% Performs a truth-model Monte-Carlo simulation for the discrete-time
% stochastic system model:
%   x(k+1) = F*x(k) + Gamma*v(k)
%   z(k) = H*x(k+1) + w(k)
%
% Where v(k) and w(k) are uncorrelated, zero-mean, white-noise Gaussian random
% processes with covariances E[v(k)*v(k)???] = Q and E[w(k)*w(k)???] = R. The
% simulation starts from an initial x(0) that is drawn from a Gaussian
% distribution with mean xbar0 and covariance P0. The simulation starts at
% time k = 0 and runs until time k = kmax.
%
% INPUTS
%   F ----------- nx-by-nx state transition matrix
%   Gamma ------- nx-by-nv process noise gain matrix
%   H ----------- nz-by-nx measurement sensitivity matrix
%   Q ----------- nv-by-nv symmetric positive definite process noise covariance
%       matrix.
%   R ----------- nz-by-nz symmetric positive definite measurement noise
%       covariance matrix.
%   xbar0 ------- nx-by-1 mean of probability distribution for initial state
%   P0 ---------- nx-by-nx symmetric positive definite covariance matrix
%       associated with the probability distribution of the initial
%       state.
%   kmax -------- Maximum discrete-time index of the simulation
%
% OUTPUTS
%   xhist ------- (kmax+1)-by-nx matrix whose kth row is equal to x(k-1)???. Thus,
%       xhist = [x(0), x(1), ..., x(kmax)]???.
%   zhist ------- kmax-by-nz matrix whose kth row is equal to z(k)???. Thus, zhist
%       = [z(1), z(2), ..., z(kmax)]. Note that the state vector
%       xhist(k+1,:)??? and the measurement vector zhist(k,:)???
%       correspond to the same time.
% ----------------------------------------------------------------------- %

Nx = length(xbar0); 
Nz = length(R); 

% P0 --> P 
Ra_p  = chol(P0); 
% build realization off of xhat0 --> xhat0 becomes mean of normally
% distributed x 
x0 = Ra_p' * randn(Nx, 1) + xbar0; 

% Q --> v 
Ra_v = chol(Q);                   % Cholesky of Covariance Matrix
% v    = Ra_v' * randn(Nx, 1); 

% R --> w 
Ra_w = chol(R); 
% w    = Ra_w' * randn(Nz,1); 

% Hints: Use the Cholesky factorization function chol to generate your initial truth state
% and the truth process and measurement noise vectors. Also use the Gaussian random
% number generation function, randn, which samples a Gaussian distribution with a mean
% equal to zero and covariance equal to the identity matrix. This function???s outputs are
% uncorrelated from call to call. You can check whether you have used chol and randn
% to correctly define the various error and noise vectors by re-computing their covariances
% based on your knowledge of how chol and randn work.

% v(k) and w(k) are uncorrelated, zero-mean, white-noise Gaussian random
% processes with covariances E[v(k)*v(k)???] = Q and E[w(k)*w(k)???] = R. The
% simulation starts from an initial x(0) that is drawn from a Gaussian
% distribution with mean xbar0 and covariance P0. The simulation starts at
% time k = 0 and runs until time k = kmax.

x = x0; 

for k = 1 : kmax
    
    v = Ra_v' * randn(1, 1); 
    w = Ra_w' * randn(Nz,1); 

    x(:,k+1) = F*x(:,k) + Gamma * v; 
    z(:,k)   = H*x(:,k+1) + w; 

end 

% Assign outputs 
xhist = x'; 
zhist = z'; 

end 



%% subfunctions KALMAN FILTER 

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
    P_cell{k} = P; 
    Pxx_arr  = [Pxx_arr; P(1,1)]; 
    Pzz_arr  = [Pzz_arr; P(2,2)]; 

end 

end 

