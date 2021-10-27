clear; clc

J_fn = @(x) ... 
    [ 1+x(2),   1+x(1)      ; 
      2*x(1),   2-2*x(2)    ];

fn = @(x) ... 
    [ x(1) + x(2) + x(1)*x(2) + 5   ; 
      x(1)^2 + 2*x(2) - x(2)^2 - 2  ]; 

% first guess 
xg = [-5; 5]; 
fprintf('initial guess f norm') 
norm(fn(xg))

% first 6 iterates 
for i = 1:6 
    
    fprintf('\n iter %d \n', i); 
    fprintf('xg = \n'); 
    xg = xg - inv(J_fn(xg)) * fn(xg); 
    disp(xg); 
    fprintf('f norm = '); 
    disp((norm(fn(xg))))

end 

