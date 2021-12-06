clear; clc

rng(0)

%% part 4 

syms xa xb pa pb 

h = ( (xa - pa)^2 + (xb - pb)^2 )^(1/2); 

H_fn = matlabFunction(jacobian(h, [pa pb])); 

syms sigma_r 

R = diag([sigma_r^2 sigma_r^2]); 
R_fn = matlabFunction(R); 

% P = ( H.' * R^-1 * H )^-1;

%% part 5 

p0a = 0; 
p0b = 0; 
x1a = -10^4; 
x1b = 10^5; 
x2a = 10^4; 
x2b = 10^5; 

% H = [ H_fn(p0a, p0b, x1a, x2b); 
%       H_fn(p0a, p0b, x2a, x2b) ]; 
%   
% R = sigma_r^2 * eye(length(H)); 
% P = inv( H.' * inv(R) * H); 

d1a = p0a - x1a; 
d2a = p0a - x2a; 
d1b = p0b - x1b; 
d2b = p0b - x2b; 
r1  = sqrt( ( x1a - p0a )^2 + ( x1b - p0b )^2 ); 
r2  = sqrt( ( x2a - p0a )^2 + ( x2b - p0b )^2 ); 

T = [ d1a^2/r1^2 + d2a^2/r2^2,      d1a*d1b/r1^2 + d2a*d2b/r2^2; 
          d1a*d1b/r1^2 + d2a*d2b/r2^2,  d1b^2/r1^2 + d2b^2/r2^2]; 
      
H = [ d1a / r1, d1b / r1; d2a / r2, d2b / r2 ]; 
      
GDOP = sqrt( trace(inv(T)) ); 




