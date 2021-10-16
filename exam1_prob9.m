% Exam 1, Problem 9 
% Stat Estimation Theory Fall 2021 
% Junette Hsin 

z = [ 
    367.426067023724
    23779.4593915456
    12846.8299987297
    -18662.4933331644
    7984.52275474047
    907.737690145488 ]; 

H = [ 
    122.747974553825 -700.205624002536 -288.72922150566
    -6161.41287306561 -17321.8364059845 -11741.6232165692
    -8745.71054008378 1414.28877217844 -3665.71330902301
    8529.78921867569 6245.78023146386 7387.79003058372
    -7900.66437226995 5784.61586018735 -1058.01203361061
    2945.79594212046 -6989.02671091266 -2021.60425818948 ]; 

R = [ 
    2.36519215943924 0.057845076600491 1.01307108063044 -0.558431761588577 1.81387196612649 -0.394688155222653
    0.057845076600491 1.85476678558121 -0.611482029811377 1.69343840421604 2.62478404586119 2.34080310702718
    1.01307108063044 -0.611482029811377 8.56587413452779 1.37856896419802 2.75188091884609 -2.38121782107832
    -0.558431761588577 1.69343840421604 1.37856896419802 5.17266218453676 3.29441874285897 1.78429853967157
    1.81387196612649 2.62478404586119 2.75188091884609 3.29441874285897 9.17143948110111 0.447330044307135
    -0.394688155222653 2.34080310702718 -2.38121782107832 1.78429853967157 0.447330044307135 4.62657818768543 ]; 

[x_hat_scribls, Ro_tilde] = sribls(H, z, R) 

[x_hat_nls] = nls(H, z, R); 

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
%   - R ----------- n_z-by-n_z measurement noise covariance matrix.
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

%% Normal least squares 

function [x_hat] = nls(H, z, R)
% ------------------------------------------------------------------------
% nls : Normal least squares routine. Given the
% linear measurement model:
%   z = H * x + w, w ~ N(0,R)
% nls returns xhat, the least squares (and Maximum Likelihood)
% estimate of x. 
%
% INPUTS
%   - H ----- n_z-by-n_x measurement sensitivity matrix.
%   - z ----- n_z-by-1 measurement vector.
%   - R ----- n_z-by-n_z measurement noise covariance matrix.
%
% OUTPUTS
%   - x_hat ------- n_x-by-1 Maximum Likelihood estimate of x.
%
% ------------------------------------------------------------------------
% References: Todd Humphreys 
%
% Author: Junette Hsin 
% ------------------------------------------------------------------------

% Normal least squares 
x_hat = inv(H' * inv(R) * H) * H' * inv(R) * z; 

end 