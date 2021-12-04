function z = h_car(xaug)
%
%Computes measurements to a car modeled geometrically as a rectangle.
%Measurements are computed via raytracing to the corners and sides of the
%rectangle.  Measurements computed are smallest and largest bearings, and
%range to closest point.
%
%INPUTS:
%   xaug - car state vector augmented with process noise.  Vector consists
%   of:
%       x - vehicle x-position (back axle, m)
%       y - vehicle y-position (back axle, m)
%       s - vehicle forward speed (no-slip assumed)
%       t - vehicle heading (rad.)
%       l - vehicle length (m)
%       w - vehicle width (m)
%       ebs - extra white noise corrupting smallest bearing measurement
%       ebl - extra white noise corrupting largest bearing measurement
%       er - extra white noise corrupting range measurement
%
%OUTPUTS:
%   z - the vector of measurements [bmin, bmax, rmin]: smallest bearing,
%       largest bearing, and smallest range

%extract elements from the augmented state:
xc = xaug(1);
yc = xaug(2);
hdg = xaug(4);
len = xaug(5);
wid = xaug(6);

%extract measurement noise from augmented state
ebs = xaug(7);
ebl = xaug(8);
ebr = xaug(9);

%calculate the box/car corners
Rot = [cos(hdg), -sin(hdg); ...
    sin(hdg), cos(hdg)];
%the corners of the car: (front/back) (driver/passenger)
fd = [xc; yc] + Rot*[len; wid/2];
fp = [xc; yc] + Rot*[len; -wid/2];
bd = [xc; yc] + Rot*[0; wid/2];
bp = [xc; yc] + Rot*[0; -wid/2];
corners = [fd, fp, bd, bp];

%start by finding the angles and ranges to the four corners
rc = sqrt(sum(corners.^2, 1));
ac = atan2(corners(2, :), corners(1, :));
%detect and correct wraparound problems
if (abs(max(ac) - min(ac)) > 2*pi)
    loc = find(ac < 0);
    ac(loc) = ac(loc) + 2*pi;
end

%next find the (one or) two exposed faces of the rectangle
%the segments are between the minimum angle pt, the minimum range pt, and 
%the max angle pt.
segs = zeros(2, 3);
[mina, locmin] = min(ac);
segs(:, 1) = corners(:, locmin);
[minr, locr] = min(rc);
mida = ac(locr);
segs(:, 2) = corners(:, locr);
[maxa, locmax] = max(ac);
segs(:, 3) = corners(:, locmax);

%calculate the vector direction of each segment
if (locmin == locr)
    %only 1 segment visible
    %update minr by solving a perpendicular distance problem
    segs = segs(:, 2:3);
    
    ps1 = segs(:, 1);
    dirs1 = segs(:, 2) - segs(:, 1);
    dirs1 = dirs1 / norm(dirs1);
    pss = [0; 0];
    dirss = [0, -1; 1, 0]*dirs1;
    b = ps1 - pss;
    A = [-dirs1, dirss];
    w = A\b;
    minr = abs(w(2));
elseif (locmax == locr)
    segs = segs(:, 1:2);
    ps1 = segs(:, 1);
    dirs1 = segs(:, 2) - segs(:, 1);
    dirs1 = dirs1 / norm(dirs1);
    pss = [0; 0];
    dirss = [0, -1; 1, 0]*dirs1;
    b = ps1 - pss;
    A = [-dirs1, dirss];
    w = A\b;
    minr = abs(w(2));
end
%if 2 segments are visible, the closest point is a corner

%calculate the meta measurements from the range measurements
z = zeros(3, 1);
%1. bearing of most CW edge
z(1) = min(ac) + ebs;
%2. bearing of most CCW edge
z(2) = max(ac) + ebl;
%3. range to closest point
z(3) = minr + ebr;

return;
