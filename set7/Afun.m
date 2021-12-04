function [A] = Afun(rv,mu)
%  AFUN    Generate the linearized state matrix A for the orbit 
%          propagation problem.  
%
%  INPUTS
%  xbar         nx-by-1 state about which the nonlinear dynamics function
%               f will be linearized, where xdot = f(x,t).  The state
%               vector elements are :
%               x = [rs;vs]
%               with rs the 3-by-1 position of the SV expressed in ECI
%               coordinates (J2000), and vs is the 3-by-1 velocity of the SV
%               with respect to the ECI coordinate system and expressed in
%               ECI coordinates.
%                            
%  mu           Earth's gravitational parameter (G*Mearth) (m^3/s^2)
%
%  OUTPUTS
%  A            The nx-by-nx linearized state matrix A where 
%               A = dftilde/dxtilde_k(t)
%+------------------------------------------------------------------+
% References:
%
%
% Author:  Todd Humphreys
%+==================================================================+
  
% do I even need this 
nx = size(rv,1);  % total elements in state vector x

% % Symbolically find A 
% r = sym('r', [3 1]); 
% v = sym('v', [3 1]); 
% dr = v; 
% rmag = norm(r); 
% dv = -mu*r / rmag^3; 
% 
% % linearize dynamics 
% A = jacobian(dr, r); 
% B = jacobian(dr, v); 
% C = jacobian(dv, r); 
% D = jacobian(dv, v); 
% 
% F = [A B; C D]; 
% F_fn = matlabFunction(F); 

F_fn = ... 
    @(r1,r2,r3)reshape([0.0,0.0,0.0,1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(3.0./2.0).*-3.986005e+14+r1.*abs(r1).*sign(r1).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,r2.*abs(r1).*sign(r1).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,r3.*abs(r1).*sign(r1).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,0.0,0.0,0.0,r1.*abs(r2).*sign(r2).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(3.0./2.0).*-3.986005e+14+r2.*abs(r2).*sign(r2).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,r3.*abs(r2).*sign(r2).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,0.0,0.0,0.0,r1.*abs(r3).*sign(r3).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,r2.*abs(r3).*sign(r3).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(3.0./2.0).*-3.986005e+14+r3.*abs(r3).*sign(r3).*1.0./(abs(r1).^2+abs(r2).^2+abs(r3).^2).^(5.0./2.0).*1.1958015e+15,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0],[6,6]) ; 

% substitute 
F = F_fn(rv(1), rv(2), rv(3)); 

% output linearized dyamics 
A = F; 
  
end 

  
    
  
  