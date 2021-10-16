%% set 4, prob 4 
clear;clc 

% loads rhoahist, rhobhist, and thist 
load radarmeasdata_missle.mat

o_pa = 10; 
o_pb = 30; 

R = [o_pa^2 0; 0 o_pb^2]; 

Ra = chol(R); 

%% Initial condition guessing 
% p_ai = [ ( l_a - y_1i )^2 + (y_2i)^2 ]^(1/2)
% p_bi = [ ( l_b - y_1i )^2 + (y_2i)^2 ]^(1/2) 

% second measurements 
p_as = rhoahist(2); 
p_bs = rhobhist(2); 

global la lb 
la = 4.1e5; 
lb = 4.4e5; 

y_1s = 1/( 2*la - 2*lb ) * ( p_bs^2 - p_as^2 - lb^2 + la^2 ); 
y_2s = sqrt( p_bs^2 - ( lb - y_1s )^2 ); 

% last measurements 
p_af = rhoahist(end); 
p_bf = rhobhist(end); 

y_1f = 1/( 2*la - 2*lb ) * ( p_bf^2 - p_af^2 - lb^2 + la^2 ); 
y_2f = sqrt( p_bf^2 - ( lb - y_1f )^2 ); 

% guessing x1 (y10) and x2 (v10) 
% y_1s = (1)*y10 + (ts)*v10 
% y_1f = (1)*y10 + (tf)*v10 
ts = 5; tf = 115; 
x = inv( [ 1 ts; 1 tf ] ) * [y_1s; y_1f]; 
y_10 = x(1); 
v_10 = x(2); 

% guessing x3 (y20) and x4 (v20) 
% y_2s = (1)*y20 + (ts)*v20 - 4.9ts^2
% y_2f = (1)*y20 + (tf)*v20 - 4.9tf^2
x = inv( [ 1 ts; 1 tf ] ) * ( [ y_2s; y_2f ] + 4.9 * [ ts^2; tf^2 ] ); 
y_20 = x(1); 
v_20 = x(2); 

%% First cost function 
xg0_OG = [y_10; v_10; y_20; v_20]; 
xg0 = xg0_0G; 

% Normalized NL at guess 
h = inv(Ra') * h_NL(xg0, thist); 

% Normalized jacobian at guess 
H = inv(Ra') * Hhist(xg0, thist); 

% Normalized measurement 
z = inv(Ra') * zhist; 

% Gauss-Newton dx 
dx = inv((H' * H)) * H' * (z - h); 

% Cost function 
Jg = norm(z - h); 

% first a step 
a = 1; 

%% First new cost function 
xg = xg0 + a * dx'; 

% Normalized NL at guess 
h = inv(Ra') * h_NL(xg, thist); 

% Normalized jacobian at guess 
H = inv(Ra') * Hhist(xg, thist); 

% Normalized measurement 
z = inv(Ra') * zhist; 

% Cost function 
Jgnew = norm(z - h); 

% Gauss-Newton dx 
dx = inv((H' * H)) * H' * (z - h); 

%% The while loop: Jgnew > Jg 

while norm(dx) > 0.00001 
    
    while Jgnew >= Jg 

        % Next a 
        a = a/2; 

        % First guess + dx 
        xg = xg0 + a * dx'; 

        % Normalized NL at guess 
        h = inv(Ra') * h_NL(xg, thist); 

        % Normalized jacobian at guess 
        H = inv(Ra') * Hhist(xg, thist); 

        % Normalized measurement 
        z = inv(Ra') * zhist; 

        % Cost function 
        Jgnew = norm(z - h); 

    end 
    
    %% While loop: "New" first guess - saved from last iteration 
    
    xg0 = xg; 
    
    Jg = Jgnew; 
    
    % Gauss-Newton dx (H, z, and h saved from last iteration) 
    dx = inv((H' * H)) * H' * (z - h); 

    % first a step 
    a = 1; 

    %% While loop: "new" first guess + dx 
    xg = xg0 + a * dx'; 

    % Normalized NL at guess 
    h = inv(Ra') * h_NL(xg, thist); 

    % Normalized jacobian at guess 
    H = inv(Ra') * Hhist(xg, thist); 

    % Normalized measurement 
    z = inv(Ra') * zhist; 

    % Cost function 
    Jgnew = norm(z - h); 

end 

%% output 

% original initial guess 
xg0_OG 

% solution to initial guess 
xg0

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

function H = Hhist(x, t) 
% Jacobian of h 

    Hrow = @(x, t) [ cos( x(2)*t +x(3) ) , ... 
                    -x(1)*sin( x(2)*t + x(3) )*t , ... 
                    -x(1) * sin( x(2)*t + x(3) ) ]; 
                
    H = [ Hrow(x, t(1)) ]; 

    for j = 2:11 
        H = [ H ; Hrow(x, t(j)) ]; 
    end 
    
end 


  
      
      
