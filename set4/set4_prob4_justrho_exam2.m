%% set 4, prob 4 
clear; clc 

% loads rhoahist, rhobhist, and thist 
% load radarmeasdata_missle.mat
load radarmeasdata_missile_new.mat

global la lb 
la = 3.5e5; 
lb = 4.0e5; 

sigma_rhoa = 10; 
sigma_rhob = 30; 

R_j = [sigma_rhoa^2 0; 0 sigma_rhob^2]; 

% Build full R matrix 
R = zeros(length(thist)); 
for j = 1:length(thist)
    R(2*j-1 : 2*j, 2*j-1 : 2*j) = R_j; 
end 

Ra = chol(R); 

% build full zhist 
zhist = []; 
for j = 1:length(thist) 
    zhist = [ zhist; rhoahist(j); rhobhist(j) ]; 
end 

%% Initial condition guessing 

xg0_arr = []; 
% First guess 
for i = 1:5
    for f = 23:28 
        sprintf('i = %d, f = %d', i, f)
        xg0 = find_xg0(rhoahist, rhobhist, thist, i, f); 
        xg0_arr = [xg0_arr; xg0']; 
    end
end

y = [xg0_arr(:,1), xg0_arr(:,3)]; 
ftitle = 'IC guessing: first 5 and last 5 range meas'; 
figure('name', ftitle); 
    subplot(2,1,1) 
        plot(y(:,1), y(:,2),'.')
        grid on; hold on; 
        yline(0, 'r') 
        xlabel('y1'); ylabel('y2'); 
        bigger_ylim; bigger_xlim 
    subplot(2,1,2) 
        plot(xg0_arr(:,3), xg0_arr(:,4), '.'); 
        grid on; hold on; 
        xlabel('v1'); ylabel('v2'); 
        bigger_ylim; bigger_xlim 
    title(ftitle); 
        
%% this one looks good 

xg0_OG = find_xg0(rhoahist, rhobhist, thist, 3, 25); 
xg0 = xg0_OG; 

%% Jacobian H 

x = sym('x', [4 1]); 
syms la_sym lb_sym tj g 

y1 = x(1) + x(2)*tj; 
dy_1a = la_sym - y1; 
dy_1b = lb_sym - y1; 
dy_2 = x(3) + tj * x(4) - 4.9*tj^2; 

ha = sqrt( dy_1a^2 + dy_2^2 ); 
hb = sqrt( dy_1b^2 + dy_2^2 ); 

% inputs: la, lb, tj, x1, x2, x3, x4 
Hhist_j = matlabFunction( [ jacobian(ha, x); jacobian(hb, x) ] ); 

% if no symbolic toolbox - here is Hhist_j copied from comand window 
% Hhist_j = @(la_sym,lb_sym,tj,x1,x2,x3,x4)reshape([(1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(la_sym.*-2.0+x1.*2.0+tj.*x2.*2.0))./2.0,(1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(lb_sym.*-2.0+x1.*2.0+tj.*x2.*2.0))./2.0,tj.*1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(-la_sym+x1+tj.*x2),tj.*1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(-lb_sym+x1+tj.*x2),(1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3.*2.0+tj.*x4.*2.0-tj.^2.*(4.9e+1./5.0)))./2.0,(1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3.*2.0+tj.*x4.*2.0-tj.^2.*(4.9e+1./5.0)))./2.0,tj.*1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)),tj.*1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1))],[2,4]); 

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
        
        % Range measurements 
        p_a = sqrt( dy_1a^2 + dy_2^2 ); 
        p_b = sqrt( dy_1b^2 + dy_2^2 ); 
        
        % Build nonlinear h from guess 
        h = [h; p_a; p_b]; 
        
    end 

end 
  
function H = Hhist(x, thist) 
% Full jacobian of h 

    global la lb 

    % Copied from symbolic toolbox output 
    Hhist_j = @(la_sym,lb_sym,tj,x1,x2,x3,x4)reshape([(1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(la_sym.*-2.0+x1.*2.0+tj.*x2.*2.0))./2.0,(1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(lb_sym.*-2.0+x1.*2.0+tj.*x2.*2.0))./2.0,tj.*1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(-la_sym+x1+tj.*x2),tj.*1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(-lb_sym+x1+tj.*x2),(1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3.*2.0+tj.*x4.*2.0-tj.^2.*(4.9e+1./5.0)))./2.0,(1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3.*2.0+tj.*x4.*2.0-tj.^2.*(4.9e+1./5.0)))./2.0,tj.*1.0./sqrt((-la_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)),tj.*1.0./sqrt((-lb_sym+x1+tj.*x2).^2+(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1)).^2).*(x3+tj.*x4-tj.^2.*(4.9e+1./1.0e+1))],[2,4]); 

    H = []; 
    for j = 1:length(thist)

        H = [ H; Hhist_j(la, lb, thist(j), x(1), x(2), x(3), x(4)) ]; 

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


  
      
