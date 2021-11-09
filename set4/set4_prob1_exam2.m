clear; clc

J_fn = @(x) ... 
    [ 1+x(2),   1+x(1)      ; 
      2*x(1),   2-2*x(2)    ];

fn = @(x) ... 
    [ x(1) + x(2) + x(1)*x(2) + 5   ; 
      x(1)^2 + 2*x(2) - x(2)^2 - 2  ]; 

% first guess 
xg = [-10; -10]; 
fprintf('initial guess f norm') 
norm(fn(xg))

xg_arr = []; 

for i = -10 : 10  
    for j = -10 : 10 
        
        xg = [i; j]; 
        
        k = 0; 
        % first 6 iterates 
        while norm(fn(xg)) > 0.00000001 

            k = k + 1; 
            xg = xg - inv(J_fn(xg)) * fn(xg); 
            disp((norm(fn(xg))))

        end 

        fprintf('for i = %d, j = %d, converged xg = \n', i, j); 
        disp(xg); 
        fprintf('f norm = %g \n', norm(fn(xg))); 

        xg_arr = [xg_arr; i, j, xg(1), xg(2)]; 
        
    end
end

%% unique converged values 

tol   = 1e-5; 
xg_u  = rmmissing(xg_arr); 
xg_u(:,3) = round(xg_u(:,3), 5);
xg_u(:,4) = round(xg_u(:,4), 5); 

% unique x1 
row_u = unique(xg_u(:,3)); 

for j = 1:length(row_u) 
    for i = 1:length(xg_u) 
        if xg_u(i,3) == row_u(j)
            xg_u(i,5) = j; 
        end
    end 
end 

% first unique points 
temp = xg_u(:,3) == row_u(1); 
idx  = find(temp, 1, 'first'); 
u1   = xg_u(idx, 3:4); 

% second unique points 
temp = xg_u(:,3) == row_u(2); 
idx  = find(temp, 1, 'first'); 
u2   = xg_u(idx, 3:4); 

%% plot 

figure()
    plot(u1(1), u1(2), 'bp', 'linewidth', 5); grid on; hold on; 
    plot(u2(1), u2(2), 'rd', 'linewidth', 5); 
    xlim([-10 10])
    ylim([-10 10]) 
    
    
    
    
    
    
    

