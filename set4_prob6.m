%% set 4, prob 6 
clear; clc 

z = 0; 

% First guess 
xg0_OG = 1.5; 
xg0 = xg0_OG; 

%% First cost function and "new" cost function 

% First NL, jacobian, and cost fn at guess 
[Jg, h, H, dx] = cost_fn(xg0, z); 

% first a step 
a = 1; 

% First new cost function 
xg = xg0 + a * dx; 

% First new NL, jacobian, and cost fn at guess 
[Jgnew, h, H, ~] = cost_fn(xg, z); 

%% The while loop: Jgnew > Jg 

while norm(dx) > 0.000001 
    
    while Jgnew >= Jg 

        % Next a 
        a = a/2; 
        
        if a < eps 
            break
        end 
        
        [Jgnew, h, H, ~] = cost_fn(xg, z); 

    end 
    
    %% While loop: "New" first guess - saved from last iteration 
        
    if a < eps 
        break
    end 
    
    xg0 = xg; 
    Jg = Jgnew; 
    
    % Gauss-Newton dx (H, z, and h saved from last iteration) 
    dx = inv((H' * H)) * H' * (z - h); 

    % first a step 
    a = 1; 

    %% While loop: "new" first guess + dx 
    xg = xg0 + a * dx; 
 
    [Jgnew, h, H, ~] = cost_fn(xg, z); 

end 

%% output 

% original initial guess 
xg0_OG 

% solution to initial guess 
xg0

% covariance 
Pxx = inv(H' * H)

%% subfunctions 

function h = h_NL(x) 
% Nonlinear measurement h 
    
    h = atan(x); 

end 
  
function H = H_NL(x) 
% Full jacobian of h 

    H = (sec(x)).^2; 
    
end 

function [Jg, h, H, dx] = cost_fn(xg, z)

    % NL at guess 
    h = h_NL(xg); 

    % jacobian at guess 
    H = H_NL(xg); 

    % Cost function 
    Jg = norm(z - h); 
    
    % Gauss-Newton dx 
    dx = inv((H' * H)) * H' * (z - h); 

end 


  
      
