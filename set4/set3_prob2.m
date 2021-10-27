% HW Set 3, Prob 2
% Stat Estimation Theory Fall 2021 
% Junette Hsin 

% z = [ n_z x 1] = [ 5 x 1 ]
z = [   5.505751578229795
        9.754581690707774
        -3.380491358445003
        -1.240740004272318
        1.402850573448269   ]; 

% H = [ n_z x n_x ] = [ 5 x 3 ] 
H = [   -0.809498694424876 -0.809498871800032 -1.618995978525818
        -2.944284161994896 -2.944284358048384 -5.888569324509237
        1.438380292815098 1.438381712125249 2.876762701564763
        0.325190539456198 0.325190831040572 0.650382205584935
        -0.754928319169703 -0.754928121358650 -1.509856684243494 ]; 

% R = [n_z x n_z ] = [ 5 x 5 ] 
R = [   9.599563453667082 1.695803491633456 3.374804722838387 -2.036285115123790 -1.912094793736365
        1.695803491633456 22.490148367908073 2.451654084712032 9.799176417244743 7.078780745346799
        3.374804722838387 2.451654084712032 12.073999822869279 1.904627078020524 -3.804148586359607
        -2.036285115123790 9.799176417244743 1.904627078020524 6.270434397020487 4.070651550638881
        -1.912094793736365 7.078780745346799 -3.804148586359607 4.070651550638881 5.320048568546946     ]; 

[x_hat, Ro_tilde] = sribls(H, z, R);     

%% function scribls 

function [x_hat, Ro_tilde] = sribls(H_prime, z_prime, R)
% ------------------------------------------------------------------------
% sribls : Square-root-information-based least squares routine. Given the
% linear measurement model:
%   zprime = Hprime*x + w, w ~ N(0,R)
% sribls returns xhat, the least squares (and Maximum Likelihood)
% estimate of x, and Rotilde, the associated square root
% information matrix.
%
% INPUTS
%   - H_prime ----- n_z-by-n_x measurement sensitivity matrix.
%   - z_prime ----- n_z-by-1 measurement vector.
%   - R ---------- n_z-by-n_z measurement noise covariance matrix.
%
% OUTPUTS
%   - x_hat ------- n_x-by-1 Maximum Likelihood estimate of x.
%   - Ro_tilde ---- n_x-by-n_x square root information matrix, where P =
%       inv(Rotilde’*Rotilde) is the estimation error covariance matrix.
%
% ------------------------------------------------------------------------
% References: Todd Humphreys 
%
% Author: Junette Hsin 
% ------------------------------------------------------------------------

% Perform change of coordinates to normalize: Cholesky factorize R 
Ra = chol(R); 

% Change coordinates 
z = inv(Ra') * z_prime; 
H = inv(Ra') * H_prime; 

% QR factorize H 
[Q_tilde, R_tilde] = qr(H); 

% Cost function 
Ro_tilde = R_tilde(1:3, 1:3); 
z_tilde  = Q_tilde' * z; 
zo_tilde = z_tilde(1:3); 

x_hat = inv(Ro_tilde) * zo_tilde; 

end 