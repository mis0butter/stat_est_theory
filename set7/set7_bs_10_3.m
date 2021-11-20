clear;

rng(0)

rho0 = 10^5; 
theta0 = 45; % deg 

z1 = rho0; 
z2 = theta0; 

y_true = [ z1 * cosd(z2); z1 * sind(z2) ]; 

R = diag([100^2; 0.5^2]); 

w = mvnrnd([0; 0], R, 100); 

z1 = rho0 + w(:,1); 
z2 = theta0 + w(:,2); 

for i = 1:length(w)
    y(:,i) = [ z1(i) * cosd(z2(i)); z1(i) * sind(z2(i)) ]; 
end 

y = y'; 

mean(y) - y_true' 
