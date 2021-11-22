clear; clc 

% altitude + radius (km) 
r = ( 780 + 6378 ) * 1000; 

% GM 
mu = 3.986005e14; 

% Keplerian orbit; initial velocity can be found via online calculator 
r0 = [1; 0; 0] * r; 
v0 = [0; 1; 0] * 7.4623 * 1000; 

% find semi-major axis 
oe = rv2oe([r0; v0], mu); 
a = oe(1); 

% period 
P = a^(3/2); 

% sample time 
T = 1; 
t = 0 : T : 7000; 

% input 
uk = zeros(3,1); 
vk = zeros(3,1); 

%% propagate orbit 

xk = [r0; v0]; 
xhist = []; 
Fhist = []; 
Ghist = []; 
for i = 1:length(t) 

    tk = t(i); 
    [xkp1, ~, ~] = propagateOrbit(tk, T, xk, uk, vk, mu); 

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

x0 = [r0; v0]; 
[x0p1, Fk, GAMMAk] = propagateOrbit(tk, T, x0, uk, vk, mu); 

x0bar = [r0-1; v0]; 
[x0p1bar, ~, ~] = propagateOrbit(tk, T, x0bar, uk, vk, mu); 

RHS = Fk * (x0 - x0bar); 
LHS = x0p1 - x0p1bar; 

norm(LHS) - norm(RHS) 

