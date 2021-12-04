clear; clc 

% altitude + radius (m) 
r = ( 780 + 6378 ) * 1000; 

% GM 
mu = 3.986005e14; 

% Keplerian orbit  
v_mag = sqrt(mu/r); 

% initial position and velocity 
r0 = [1; 0; 0] * r; 
v0 = [0; 1; 0] * v_mag ; 

% period 
P = 2 * pi * sqrt( r^3 / mu ); 

% sample time 
% dt = 1; 
thist = linspace(0, P, round(2*P)); 
dt = thist(2) - thist(1); 

% input 
uk = zeros(3,1); 
% uk = rand(length(thist), 3); 

vk = zeros(3,1); 

%% propagate orbit 

xk = [r0; v0]; 
xhist = []; 
Fhist = []; 
Ghist = []; 
for i = 1:length(thist) 

    tk = thist(i); 
    [xkp1, ~, ~] = propagateOrbit(tk, dt, xk, uk, vk, mu); 
%     [xkp1, ~, ~] = propagateOrbit(tk, dt, xk, uk(i,:)', vk, mu); 

    xk = xkp1; 
    xhist = [xhist; xk']; 
%     Fhist = [Fhist; Fk']; 
%     Ghist = [Ghist; GAMMAk']; 
    
end 

%% Final A 

A = Afun(xhist(end,:)', mu); 

%% 

close all; 
plot3(xhist(:,1), xhist(:,2), xhist(:,3))

%% 

tk = 0; 

uk = zeros(3,1); 

x0 = [r0; v0]; 
[x1, Fk, GAMMAk] = propagateOrbit(tk, dt, x0, uk, vk, mu); 

r0_bar = [r0(1) + 10; 0; 0]; 
v0_bar = [0; v0(2) + 10; 0]; 
x0bar = [r0_bar; v0_bar]; 
[x1bar, Fk_bar, ~] = propagateOrbit(tk, dt, x0bar, uk, vk, mu); 

RHS = Fk * (x0 - x0bar); 
LHS = x1 - x1bar; 

disp('norm(LHS) - norm(RHS)')
norm(LHS) - norm(RHS) 

disp('norm(LHS - RHS)') 
norm(LHS - RHS)


