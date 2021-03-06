set5_prob3 

%% set 5 prob 4 

sys = ss(Fk, [ Gammak ], Hk, 0, -1); 
[km, L, Pbar_ss, W_ss] = kalman(sys, Qk, Rk); 

% update 
S = Hk * Pbar_ss * Hk' + Rk;               % innovation covariance 
P_ss = Pbar_ss - W_ss * S * W_ss';         % a posteriori covar est 

disp('W steady state = ')
disp(W)

disp('A posteriori covariance steady state = ')
disp(P_ss)

%% plot 

ftitle = 'Covariances Comparison'; 
figure('name', 'ftitle'); 
    subplot(2,1,1) 
        plot( thist0, Pxx_arr); hold on; grid on; 
        yline(P_ss(1,1), 'r--'); 
        legend( 'P_{11}', 'P_{ss}' ); 
        title('P_{11} vs P_{ss}') 
        bigger_ylim 
    subplot(2,1,2) 
        plot( thist0, Pzz_arr); hold on; grid on; 
        yline(P_ss(2,2), 'r--'); 
        legend( 'P_{22}', 'P_{ss}' ); 
        title('P_{22} vs P_{ss}')
        bigger_ylim 
    xlabel('Time') 

disp('A posteriori covariance converges to steady-state covariance'); 
    
%% Stability 
        
disp('Eigenvalues of [I - W_ss * H] * F ')
disp( eig( (eye(2) - W_ss * Hk) * Fk ) )

disp('Eigenvalues have complex magnitudes less than 1; error transition matrix is stable'); 
