%  Calculate the smoothed estimates for the problem in kf example03a.m. Compare xˆ(10)
% with x∗(10) and compare P(10) with P∗(10). Is P∗(10) ≤ P(10)? Do the smoothed state
% time history estimate plots look “smoother” than the filtered state time history estimate
% plots?

%% SRIF (forward dynamics) 

clear; clc 

disp('EXAMPLE 03A') 
kf_example03a; 

uk = [0; 0]; 
nx = length(xhat0); 
nv = length(Qk); 
[xhat_arr_srif, Rxx_arr, zx_arr, zv_bar_arr, Rvv_bar_cell, Rvx_bar_cell] = ... 
    srif( xhat0, P0, zhist, Fk, Gammak, Qk, Hk, Rk, uk, nx, nv); 


%% smoother (backward dynamics) 

zx_star = zx_arr(end,:)'; 
Rxx_star = Rxx_arr{end}; 
% wx_star = wx; 

% START AT k = N 
x_star = inv(Rxx_star) * zx_star; 
P_star = inv(Rxx_star) * inv(Rxx_star)'; 
Rvv = chol(inv(Qk)); 
zv  = zeros(nv, 1); 

% initialize 
N = length(zhist); 
x_star_arr = zeros(N, nx); 
P_star_cell = cell(N,1); 

for k = N-1 : -1 : 1
    
    zx_star = zx_arr(k+1, :)'; 
    Rxx_star = Rxx_arr{k+1}; 
    
    Rvv_bar = Rvv_bar_cell{k+1}; 
    Rvx_bar = Rvx_bar_cell{k+1}; 
    
    A = [ Rvv_bar + Rvx_bar * Gammak,   Rvx_bar * Fk; 
          Rxx_star * Gammak,            Rxx_star * Fk ]; 
    [QA, RA] = qr(A); 
    
    R_QR = QA' * A; 
    Rxx_star = R_QR(nv+1:end, nv+1:end); 
    
    zv_bar = zv_bar_arr(k+1,:)'; 
    z_star = QA' * [ zv_bar; zx_star ]; 
    zx_star = z_star(nv+1:end); 
    
    % extract state and covariance 
    x_star = inv(Rxx_star) * zx_star; 
    P_star = inv(Rxx_star) * inv(Rxx_star)'; 
    
    % save outputs 
    x_star_arr(k,:) = [x_star_arr; x_star]; 
    P_star_cell{k} = P_star; 
    
end 


%% subfunctions 

function [xhat_arr, Rxx_cell, zx_arr, zv_bar_arr, Rvv_bar_cell, Rvx_bar_cell] = ... 
    srif( xhat0, P0, zhist, Fk, Gk, Qk, Hk, Rk, uk, nx, nv)

% Initialize 
xhat_arr = []; 
Rxx_cell = {}; 
zx_arr = []; 
zv_bar_arr = []; 
Rvv_bar_cell = {}; 
Rvx_bar_cell = {}; 

% START AT k = 0: 
I = inv(P0); 
Rxx = chol(I); 
Rvv = chol(inv(Qk)); 
zv = zeros(nv, 1); 
zx = Rxx * xhat0; 

for k = 1:length(zhist)

    % PROPAGATION STEP 
    % a) QR factorize 
    A = [Rvv, zeros(2,3); -Rxx * inv(Fk) * Gk, Rxx * inv(Fk)]; 
    [QA, RA] = qr(A); 
    
    % [ Rvv_bar(k) [2x2], Rvx_bar(k+1) [2x3]; 
    %   0 [3x2]         , Rxx_bar(k+1) [3x3] ] = RA 

    % b) orthonormal transformation 
    B = [ zv; zx + Rxx * inv(Fk) * Gk * uk ];  
    [zv_zx_bar] = QA' * B; 

    % c) extract Rxx_bar(k+1) and zx_bar(k+1) 
    zv_bar = zv_zx_bar(1:nv); 
    zx_bar = zv_zx_bar(nv+1:end); 
    Rxx_bar = RA(nv+1:end, nv+1:end); 
    
    % extract Rvv and Rvx bars 
    Rvv_bar = RA(1:nv, 1:nv); 
    Rvx_bar = RA(1:nv, nv+1:end); 

    % MEASUREMENT UPDATE: 
    % a) Cholesky factorize R 
    Ra = chol(Rk); 

    % b) Transform z(k+1) and H(k+1) 
    za = inv(Ra)' * zhist(k); 
    Ha = inv(Ra)' * Hk; 

    % c) perform another QR factorization: 
    [QB, RB] = qr([Rxx_bar; Ha]); 

    % d) transform as 
    [zx_zr] = QB' * [zx_bar; za];    

    % e) extract Rxx_bar(k+1) and zx(k+1) 
    zx = zx_zr(1:nx); 
    Rxx = RB(1:nx, :); 
    
    xhat = inv(Rxx) * zx; 
    
    zx_arr = [zx_arr; zx']; 
    zv_bar_arr = [zv_bar_arr; zv_bar']; 
    xhat_arr = [xhat_arr; xhat']; 
    Rxx_cell{k} = Rxx;  
    Rvv_bar_cell{k} = Rvv_bar; 
    Rvx_bar_cell{k} = Rvx_bar; 
    
end 

end 