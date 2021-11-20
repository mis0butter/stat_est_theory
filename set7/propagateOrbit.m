function [xkp1,Fk,GAMMAk]=propagateOrbit(tk,T,xk,uk,vk,mu)
% propagateOrbit   Propagates the nx-by-1 orbital elements vector xk to xkp1
%                  via numerical integration.  Also solves for the linearized
%                  state transition matrix Fk = F(tkp1,tk) and the linearized
%                  process noise matrix GAMMAk = GAMMA(tkp1,tk).
%
% INPUTS
%
% tk           The time at which xk is defined (sec).
%
% T            Update interval (sec).          
% 
% xk           The nx-by-1 vector of orbital elements at time tk. 
%              xk = [rsk; vsk] where rsk and vsk are the SV's 
%              position (m) and velocity (m/s) at time tk expressed 
%              in ECI coordinates.
%
% uk           The nu-by-1 vector of control accelerations for the
%              interval tk to tkp1 expressed in ECI coordinates.
%
% vk           The nv-by-1 vector of disturbance accelerations for the
%              interval tk to tkp1 expressed in ECI coordinates. 
%
% mu           Earth's gravitational parameter (G*Mearth) (m^3/s^2)
%
% OUTPUTS
% xkp1         The nx-by-1 vector of orbital elements at time tkp1. 
%
% Fk           The linearized state transition matrix: Fk = F(tkp1,tk).
%
% GAMMAk       The linearized process noise matrix: 
%              GAMMAk = GAMMA(tkp1,tk).
%
%+------------------------------------------------------------------+
% References:  
%
%
% Author:  Todd Humphreys
%+==================================================================+

%----- Setup
nx = 6;
nu = 3;
nv = 3;
% Set options for numerical integration
options = odeset('reltol',1e-1);
% Initialize F and GAMMA 
F0 = eye(nx); 
F0 = reshape(F0, [], 1); 
GAMMA0 = zeros(nx,nv); 
GAMMA0 = reshape(GAMMA0, [], 1); 
% The nx*(nx+nv+1)-by-1 augmented state vector Xk is of the form 
% X = [x; phi1; ...; phinx; gamma1; ...; gammanv], 
% where x = [rs; vs] is the nx-by-1 vector of orbital parameters, phij is
% the nx-by-1 jth column of F(t,tk), and gammaj is the nx-by-1 jth column of
% GAMMA(t,tk).
Xk = [xk; F0; GAMMA0]; 

%----- Call the numerical integration routine (select one)
%[TVec,XMat] = ode113(@odePHI,[tk,tk+T/2,tk+T],Xk,options);
[TVec,XMat] = ode45(@odePHI,[tk,tk+T/2,tk+T],Xk,options);

%----- Unpack Xkp1 for output
Xkp1 = XMat(end,:)';
xkp1 = Xkp1(1:nx);
Fvec = Xkp1(nx+1:(nx+1)*nx);
Fk = reshape(Fvec,nx,nx);
GAMMAvec = Xkp1((nx+1)*nx + 1:end);
GAMMAk = reshape(GAMMAvec,nx,nv);

%----- odePHI is a nested function that has access to all input paramters
function [Xdot] = odePHI(t,X)

%-- Calculate the time derivative of the nx-by-1 state vector
x = X(1:nx,1);
rs = x(1:3); vs = x(4:6);
xdot = zeros(nx,1);
xdot(1:3) = X(4:6); 
xdot(4:6) = -mu * rs / norm(rs)^3; 

%-- Calculate A(t) and the time derivatives Fdot and GAMMAdot
A = Afun(x, mu);
D = [zeros(3); eye(3)];
Fvec = X(nx+1:(nx+1)*nx);
F = reshape(Fvec,nx,nx);
GAMMAvec = X((nx+1)*nx + 1:end);
GAMMA = reshape(GAMMAvec, [nx nv]);
Fdot = A * F; 
GAMMAdot = A * GAMMA; 
Xdot = [xdot;Fdot(:);GAMMAdot(:)];

end

end
  
  
  