%% set 4, prob 4 
clear; clc 

% loads rhoahist, rhobhist, and thist 
% load radarmeasdata_missle.mat
load radarmeasdata_missile_new.mat

global la lb 
la = 3.5e5; 
lb = 4.0e5; 

sigma_thetaa = 0.01; 

R_aj = zeros(1); 
R_aj(1,1) = sigma_thetaa^2; 

% Build full R matrix 
R = zeros(length(thist)); 
for j = 1:length(thist)
    R(j,j) = R_aj; 
end 

Ra = chol(R); 

% build full zhist 
zhist = []; 
for j = 1:length(thist) 
    zhist = [ zhist; thetaahist(j) ]; 
end 
        
%% best initial estimate from rhoandtheta

xg0_OG = [2009.37317814612
          899.930080622043
          2250.33006742093
          1598.79191840315 ]; 
% xg0_OG = find_xg0(rhoahist, rhobhist, thist, 3, 25) 
xg0 = xg0_OG; 

%% Jacobian H 

x = sym('x', [4 1]); 
syms la_sym lb_sym tj g 

y1 = x(1) + x(2)*tj; 
dy_1a = la_sym - y1; 
dy_1b = lb_sym - y1; 
dy_2 = x(3) + tj * x(4) - 4.9*tj^2; 

h_thetaa = atan2( dy_2, dy_1a ); 

% inputs: la, lb, tj, x1, x2, x3, x4 
Hhist_j = matlabFunction( [ jacobian(h_thetaa, x) ] ); 

% if no symbolic toolbox - here is Hhist_j copied from comand window 
% Hhist_j = @(la_sym,tj,x1,x2,x3,x4)[-(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3))./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),((imag(tj)./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))-real(tj).*(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).*1.0./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2).*(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2)./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),-(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),-((real(tj)./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))+imag(tj).*(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).*1.0./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2).*(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2)./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2)]; 

%% First cost function 

[Jg, h, H, dx] = cost_fn(xg0, thist, zhist, Ra); 

% first a step 
a = 1; 

% First step-size adjusted cost function 
xg = xg0 + a * dx; 
[Jgnew, h, H, ~] = cost_fn(xg, thist, zhist, Ra); 

% Gauss-Newton dx 
% dx = inv((H' * H)) * H' * (z - h); 

%% The while loop: Jgnew > Jg 

Jg_i = []; 
while norm(dx) > 1e-10
    
    while Jgnew >= Jg 

        % Next a 
        a = a/2; 
        if a < 0.001
            break; end 

        % Step size-adjusted guess and cost fn 
        xg = xg0 + a * dx; 
        [Jgnew, h, H, ~] = cost_fn(xg, thist, zhist, Ra); 

    end 
    
    %% While loop: "New" first guess - saved from last iteration 
        
%     if a < eps
%         break; end 
    
    xg0 = xg; 
    Jg  = Jgnew; 
    
    % Gauss-Newton dx (H, z, and h saved from last iteration) 
    z  = inv(Ra') * zhist; 
    dx = inv((H' * H)) * H' * (z - h); 

    % first a step 
    a = 1; 

    % "new" step-size adjusted guess 
    xg = xg0 + a * dx; 
    [Jgnew, h, H, ~] = cost_fn(xg, thist, zhist, Ra); 
    
    Jg_i = [Jg_i; Jg]; 

end 

xg0_sol = xg0; 

%% output 

% original initial guess 
xg0_OG 

% Gauss-Newton approximated solution 
xg0_sol

% covariance 
Pxx = inv(H' * H)

%% subfunctions 

function h = h_NL(x, t) 
% Nonlinear measurement h 
    
global la lb 

    % Initialize h 
    h = []; 
    
    for i = 1:length(t) 

        % y1 = y10 + v10*t = x1 + x2*t         
        y1 = x(1) + x(2)*t(i); 
        dy_1a = la - y1; 
        dy_1b = lb - y1; 
        
        % y2 = y20 + v20*t - 0.5 * 9.8 * t^2 
        dy_2 = x(3) + x(4)*t(i) - 4.9*t(i)^2; 

        h_thetaa = atan2( dy_2, dy_1a ); 

        % Build nonlinear h from guess 
        h = [h; h_thetaa]; 
        
    end 

end 
  
function H = Hhist(x, thist) 
% Full jacobian of h 

    global la lb 

    % Copied from symbolic toolbox output 
    Hhist_j = @(la_sym,tj,x1,x2,x3,x4)[-(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3))./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),((imag(tj)./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))-real(tj).*(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).*1.0./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2).*(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2)./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),-(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2),-((real(tj)./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1))+imag(tj).*(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).*1.0./(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2).*(imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2)./((imag(tj.^2).*(-4.9e+1./1.0e+1)+imag(tj.*x4)+real(tj.*x2)-real(la_sym)+imag(x3)+real(x1)).^2+(real(tj.^2).*(4.9e+1./1.0e+1)+imag(tj.*x2)-real(tj.*x4)-imag(la_sym)+imag(x1)-real(x3)).^2)]; 
    
    H = []; 
    for j = 1:length(thist)

        H = [ H; Hhist_j(la, thist(j), x(1), x(2), x(3), x(4)) ]; 

    end 
    
end 

function [Jg, h, H, dx] = cost_fn(xg, thist, zhist, Ra)

    % Normalized NL at guess 
    h = inv(Ra') * h_NL(xg, thist); 

    % Normalized jacobian at guess 
    H = inv(Ra') * Hhist(xg, thist); 

    % Normalized measurement 
    z = inv(Ra') * zhist; 

    % Gauss-Newton dx 
    dx = inv((H' * H)) * H' * (z - h); 

    % Cost function 
    Jg = norm(z - h); 

end 

function xg0_OG = find_xg0(rhoahist, rhobhist, thist, i, f)

clear x 

% "initial" measurements 
% i = 3; 
p_ai = rhoahist(i); 
p_bi = rhobhist(i); 

global la lb 

y_1i = 1/( 2*lb - 2*la ) * ( p_ai^2 - la^2 - p_bi^2 + lb^2); 
y_2i = sqrt( p_ai^2 - ( la - y_1i )^2 ); 

% last measurements 
% f = 26; 
p_af = rhoahist(f); 
p_bf = rhobhist(f); 

y_1f = 1/( 2*lb - 2*la ) * ( p_af^2 - la^2 - p_bf^2 + lb^2 ); 
y_2f = sqrt( p_af^2 - ( la - y_1f )^2 ); 

% guessing x1 (y10) and x2 (v10) 
% y_1s = (1)*y10 + (ts)*v10 
% y_1f = (1)*y10 + (tf)*v10 
ti = thist(i); tf = thist(f); 
x = pinv( [ 1 ti; 1 tf ] ) * [y_1i; y_1f]; 
y_10 = x(1); 
v_10 = x(2); 

% guessing x3 (y20) and x4 (v20) 
% y_2s = (1)*y20 + (ts)*v20 - 4.9ts^2
% y_2f = (1)*y20 + (tf)*v20 - 4.9tf^2
x = pinv( [ 1 ti; 1 tf ] ) * ( [ y_2i; y_2f ] + 4.9 * [ ti^2; tf^2 ] ); 
y_20 = x(1); 
v_20 = x(2); 

% SANITY CHECK linear algebra 
t = [ 0; -0.5 * 9.8 * ti^2; 0; -0.5 * 9.8 * tf^2 ]; 
y = [y_1i; y_2i; y_1f; y_2f]; 
A = [1, ti, 0, 0; 
     0, 0, 1, ti; 
     1, tf, 0, 0; 
     0, 0, 1, tf ]; 
x = pinv( A ) * (y - t); 

% First guess 
xg0_OG = [y_10; v_10; y_20; v_20]; 

end 


  
      
