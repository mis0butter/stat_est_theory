function plotcar(x, pstr, h)
%
%Plots the rectangular car specified by the state x into a figure of handle
%h.
%
%INPUTS:
%   x - state of the rectangular car to be plotted.  State elements are:
%       x - vehicle x-position (back axle, m)
%       y - vehicle y-position (back axle, m)
%       s - vehicle forward speed (no-slip assumed)
%       t - vehicle heading (rad.)
%       l - vehicle length (m)
%       w - vehicle width (m)
%   pstr - (OPTIONAL) plot string.  Describes the line and marker style to
%       be used for the car rectangle.  If not provided, function defaults
%       to 'b-'
%   h - (OPTIONAL) handle of the figure in which the car is to be plotted.
%       If not specified, a new figure is spawned for drawing.
%
%OUTPUTS:
%   Draws a figure with the specified rectangular car.

if (nargin < 2)
    pstr = 'b-';
end
if (nargin < 3)
    h = figure;
end

figure(h);
%calculate the car corners
Rot = [cos(x(4)), -sin(x(4)); ...
       sin(x(4)), cos(x(4))];
%the corners of the car: (front/back) (driver/passenger)
fd = x(1:2) + Rot*[x(5); x(6)/2];
fp = x(1:2) + Rot*[x(5); -x(6)/2];
bd = x(1:2) + Rot*[0; x(6)/2];
bp = x(1:2) + Rot*[0; -x(6)/2];

hold on;
plot([fd(1); fp(1)], [fd(2); fp(2)], pstr);
plot([fd(1); bd(1)], [fd(2); bd(2)], pstr);
plot([bd(1); bp(1)], [bd(2); bp(2)], pstr);
plot([bp(1); fp(1)], [bp(2); fp(2)], pstr);

return;
