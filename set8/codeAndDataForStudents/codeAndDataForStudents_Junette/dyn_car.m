function xdot = dyn_car(t, xaug)
%
%Calculates the state derivative for a constant-velocity car modeled as a
%box and driven at the back axle (back of the box).
%
%INPUTS:
%   xaug - current vehicle state, augmented to include process noise.
%   Consists of:
%       x - vehicle x-position (back axle, m)
%       y - vehicle y-position (back axle, m)
%       s - vehicle forward speed (no-slip assumed)
%       t - vehicle heading (rad.)
%       l - vehicle length (m)
%       w - vehicle width (m)
%       ex - white noise driving x-position
%       ey - white noise driving y-position
%       es - white noise driving speed
%       eh - white noise driving heading
%       el - white noise driving vehicle length
%       ew - white noise driving vehicle width
%
%OUTPUTS:
%   xdot - state derivative

%% calculate A 

% x = sym('x', [12 1]); 
% syms t 
% dx = sym('dx', [12 1]); 
% 
% dx(1) = x(3) * cos(t) + x(7); 
% dx(2) = x(3) * sin (t) + x(8); 
% dx(3) = x(9); 
% dx(4) = x(10); 
% dx(5) = x(11); 
% dx(6) = x(12); 
% dx(7:12) = zeros(6,1); 

%% state derivative 

%extract speed from the state vector
s = xaug(3);
%extract heading from the state vector
t = xaug(4);

%extract the process noise from the augmented state vector
ex = xaug(7);
ey = xaug(8);
es = xaug(9);
eh = xaug(10);
el = xaug(11);
ew = xaug(12);

%calculate the state derivative
xdot = zeros(12, 1);
xdot(1) = s*cos(t) + ex;
xdot(2) = s*sin(t) + ey;
%assume speed is "constant," driven only by noise
xdot(3) = es;
%assume heading is "constant," driven only by noise
xdot(4) = eh;
%assume size is "constant," driven only by noise
xdot(5) = el;
xdot(6) = ew;
%assume each of the process noise inputs are zero order hold
xdot(7:12) = 0;

return;
