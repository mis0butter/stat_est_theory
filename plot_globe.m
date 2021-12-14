function plot_globe(pos_ecef, scan_pixel,scan_pixel_dots,time_mode_change,E_rg_d_mode_change)

 % Textured 3D Earth example
%
% Ryan Gray
% 8 Sep 2004
% Revised 9 March 2006, 31 Jan 2006, 16 Oct 2013

%% Options

space_color = 'k';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 0.5; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = 0; % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.

image_file = 'earth.jpg';

% Mean spherical earth

erad    = 6371008.7714; % equatorial radius (meters)
prad    = 6371008.7714; % polar radius (meters)
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure


hold on;

% Turn off the normal axes
set(gcf, 'Color', space_color)
set(gca, 'NextPlot','add', 'Visible','off');


axis equal;
axis auto;

% Set initial view

view(0,30);

axis vis3d;

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);
% globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 'none');

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%% Texturemap the globe

% Load Earth image for texture map

cdata = imread(image_file,'jpg');

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Keep the mesh edges.

set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'Clipping','off')
hold on;
% keep clipping off for better zooming, and match colors to DIS plots
plot3(pos_ecef(1,:),pos_ecef(2,:),pos_ecef(3,:),'r','linewidth', 2,'Clipping','off');
if exist('scan_pixel')
    plot3(scan_pixel(1,:), scan_pixel(2,:), scan_pixel(3,:), 'r','Clipping','off')
end
if nargin > 3 && ~isempty(scan_pixel_dots)
    % if available, plot dots for times when scan is in progress
    plot3(scan_pixel_dots(1,:), scan_pixel_dots(2,:), scan_pixel_dots(3,:), 'g.','Clipping','off')
end
if nargin > 5 && ~isempty(time_mode_change)
    % if available, put time labels on mode changes
    text(E_rg_d_mode_change(1,:),...
        E_rg_d_mode_change(2,:),...
        E_rg_d_mode_change(3,:),...
        num2str(time_mode_change,'%-.4f'),...
        'Color','w')
end    
hold off;

end