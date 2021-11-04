% Develop a truth model simulation of a stochastic linear system and use it to test the
% consistency of a Kalman filter by doing Monte Carlo runs. Your truth model simulation
% should adhere to the following interface, which you can cut-and-paste from the pdf into
% Matlab:

% Use this truth model simulation and a Kalman filter to do a Monte Carlo study of the
% consistency of the Kalman filter. Test the Kalman filter’s consistency on the Kalman
% filtering problem described in Problem 3. Use Monte Carlo techniques to calculate
% E[x˜(10)], E[x˜(10)x˜T(10)], E[x˜(35)], and E[x˜(35)x˜T(35)] . Show that E[x˜(10)] and
% E[x˜(35)] approach zero as the number of Monte Carlo runs gets to be large. 
% Recall that x˜(k) = x(k) − xˆ(k). Show that E[x˜(10)x˜T(10)] approaches P(10) and that
% E[x˜(35)x˜T(35)] approaches P(35) as the number of Monte Carlo runs gets to be large.
% Do one Monte Carlo run that uses 50 simulations and do another that uses 10,000 simulations and compare your results. As part of the record of your work include your calculated
% values for E[x˜(10)], E[x˜(10)x˜T(10)], E[x˜(35)], and E[x˜(35)x˜T(35)] and any code that
% you used to perform this study.

% Hints: Use the Cholesky factorization function chol to generate your initial truth state
% and the truth process and measurement noise vectors. Also use the Gaussian random
% number generation function, randn, which samples a Gaussian distribution with a mean
% equal to zero and covariance equal to the identity matrix. This function’s outputs are
% uncorrelated from call to call. You can check whether you have used chol and randn
% to correctly define the various error and noise vectors by re-computing their covariances
% based on your knowledge of how chol and randn work.

clear
% clc 

kf_example02a; clear thist; clear zhist; 

xtilde10 = []; 
xtilde35 = []; 
xtilde = []; 

for i = 1 : 10000 % run 50 or 10000 monte carlos 

    kmax = 50; % zhist --> 10 or 35 
    [xhist, zhist] = mcltisim(Fk, Gammak, Hk, Qk, Rk, xhat0, P0, kmax); 

    % Pass into kalman filter 
    [xhat, Pxx, Pzz, P_cell] = kf( ... 
        xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk );

    % xtilde = truth - estimate 
    xtilde10 = [xtilde10; xhist(10,:) - xhat(10,:)]; 
    xtilde35 = [xtilde35; xhist(35,:) - xhat(35,:)]; 
    xtilde = [xtilde; xhist - xhat]; 
    
    if i == 50
        disp('i = 50') 
        mean(xtilde10) 
        mean(xtilde35) 
        cov(xhist(10,:) - xhat(10,:)); 
    elseif i == 10000 
        disp('i = 10000') 
        mean(xtilde10) 
        mean(xtilde35) 
    end 

end 







%% subfunctions 

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
% processes with covariances E[v(k)*v(k)’] = Q and E[w(k)*w(k)’] = R. The
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
%   xhist ------- (kmax+1)-by-nx matrix whose kth row is equal to x(k-1)’. Thus,
%       xhist = [x(0), x(1), ..., x(kmax)]’.
%   zhist ------- kmax-by-nz matrix whose kth row is equal to z(k)’. Thus, zhist
%       = [z(1), z(2), ..., z(kmax)]. Note that the state vector
%       xhist(k+1,:)’ and the measurement vector zhist(k,:)’
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
% equal to zero and covariance equal to the identity matrix. This function’s outputs are
% uncorrelated from call to call. You can check whether you have used chol and randn
% to correctly define the various error and noise vectors by re-computing their covariances
% based on your knowledge of how chol and randn work.

% v(k) and w(k) are uncorrelated, zero-mean, white-noise Gaussian random
% processes with covariances E[v(k)*v(k)’] = Q and E[w(k)*w(k)’] = R. The
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
    P_cell   = {P_cell; P}; 
    Pxx_arr  = [Pxx_arr; P(1,1)]; 
    Pzz_arr  = [Pzz_arr; P(2,2)]; 

end 

end 

