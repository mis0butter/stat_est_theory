% altitude + radius (km) 
r = ( 780 + 6378 ) * 1000; 

% Keplerian orbit; initial velocity can be found via online calculator 
r0 = [1; 0; 0] * r; 
v0 = [0; 1; 0] * 7.4623 * 1000; 

T = 1; 
t = 0 : T : 7000; 

uk = zeros(3,1); 
vk = zeros(3,1); 

mu = 3.986005e14; 

xk = [r0; v0]; 
xhist = []; 
Fhist = []; 
Ghist = []; 
for i = 1:length(t) 

    tk = t(i); 
    [xkp1, Fk, GAMMAk] = propagateOrbit(tk, T, xk, uk, vk, mu); 

    xk = xkp1; 
    xhist = [xhist; xk']; 
    Fhist = [Fhist; Fk']; 
    Ghist = [Ghist; GAMMAk']; 
    
end

close all; 
plot3(xhist(:,1), xhist(:,2), xhist(:,3))

