function [meanCar,covCar] = cartesianTransform(rho,theta,sigma_rho,sigma_theta)

meanCar = [rho*cos(theta) rho*sin(theta)]';

covCar = [sigma_rho^2*cos(theta)^2 + sigma_theta^2*rho^2*sin(theta)^2 ...;
    .5*(sigma_rho^2 - sigma_theta^2*rho^2)*sin(2*theta);
    .5*(sigma_rho^2 - sigma_theta^2*rho^2)*sin(2*theta) ...;
    sigma_rho^2*sin(theta)^2 + sigma_theta^2*rho^2*cos(theta)^2];


end

