%% set 4, prob 1 

clear; clc

% nonlinear equations f1 and f2  
fn = @(x) ... 
    [ x(1) + x(2) + x(1)*x(2) + 5   ; 
      x(1)^2 + 2*x(2) - x(2)^2 - 2  ]; 

% jacobian of f1 and f2 
J_fn = @(x) ... 
    [ 1+x(2),   1+x(1)      ; 
      2*x(1),   2-2*x(2)    ];

  
%% initial conditions 

% first guess 
xg = [4; -4]; 
fprintf('\n\n FIRST GUESS: xg = %g, %g \n', xg(1), xg(2)); 
fprintf('initial guess f norm = %g \n', norm(fn(xg))) 

iter_6(xg, J_fn, fn)

% second guess 
xg = [6; 0]; 
fprintf('\n\n SECOND GUESS: xg = %g, %g \n', xg(1), xg(2)); 
fprintf('initial guess f norm = %g \n', norm(fn(xg))) 

iter_6(xg, J_fn, fn)

% third guess 
xg = [-5; 5]; 
fprintf('\n\n THIRD GUESS: xg = %g, %g \n', xg(1), xg(2)); 
fprintf('initial guess f norm = %g \n', norm(fn(xg))) 

iter_6(xg, J_fn, fn)

%% subfunctions 

% first 6 iterates 

function iter_6(xg, J_fn, fn)
    for i = 1:6 

        fprintf('\n ITER %d \n\n', i); 
        fprintf('xg = \n'); 
        xg = xg - inv(J_fn(xg)) * fn(xg); 
        disp(xg); 
        fprintf('f norm = %g \n', norm(fn(xg))); 

    end 
end 

