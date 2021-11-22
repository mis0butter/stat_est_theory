function [oe, evec] = rv2oe(rv, mu)
% ------------------------------------------------------------------------
% Inputs 
%   rv = [6x1] position and velocity states vector 
% 
% Outputs 
%   oe = [6x1] orbital elements: a, e, i, w, Omega, nu
%           a       = semimajor axis 
%           e       = eccentricity 
%           i       = inclination 
%           w       = argument of perigee 
%           Omega   = right ascension of ascending node 
%           nu      = true anomaly 
% ------------------------------------------------------------------------

r = rv(1:3); 
v = rv(4:6); 

% angular momentum 
h       = cross(r,v); 

% node vector 
nhat    = cross( [0 0 1], h ); 

% eccentricity 
evec    = ( (norm(v)^2 - mu/norm(r))*r - dot(r,v)*v ) / mu; 
e       = norm(evec); 

% specific mechanical energy 
energy  = norm(v)^2/2 - mu/norm(r); 

% semi-major axis and p
if abs(e-1.0)>eps
   a = -mu/(2*energy); 
   p = a*(1-e^2); 
else
   p = norm(h)^2/mu; 
   a = inf; 
end

% inclination 
i = acos(h(3)/norm(h)); 

% right ascension of ascending node (check for equatorial orbit) 
if i > 0.000001
    Omega = acos( nhat(1)/norm(nhat) ); 
else
    Omega = 0; 
end
if isnan(Omega)
    Omega = 0; 
end
% if n(2)<0
%    Omega = 360-Omega; 
% end

% argument of perigee 
if e > 0.000001
    w = acos(dot(nhat,evec)/(norm(nhat)*e)); 
else
    w = 0; 
end
if isnan(w)
    w = 0; 
end
% if e(3)<0
%    argp = 360-argp
% end

% true anomaly 
nu = acos( dot(evec,r) / (e*norm(r)) );  
% if dot(r,v)<0
%    nu = 360 - nu
% end

oe = [a, e, i, w, Omega, nu]; 

end