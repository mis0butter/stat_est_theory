function [zbar,total] = cartesianTransformUnscented(rho,theta,sigma_rho,sigma_theta,alpha,beta,kappa)

% example call: [zbar,total] = cartesianTransformUnscented(rho,theta,sigma_rho,sigma_theta,alpha,beta,kappa)

nx = 2;
S = chol(blkdiag(sigma_rho^2,sigma_theta^2))';

tempMean = [rho theta]';
sps = 2*nx + 1;
lambda1 = alpha^2*(nx + kappa) - nx; 
scale = sqrt(nx + lambda1); 

x0 = tempMean;

% positive pass and negative pass
for i = 1:nx
    xplus(:,i) = x0 + scale*S(:,i);
    xminus(:,i) = x0 - scale*S(:,i);
end

sigmaPoints = [x0 xplus xminus];
z_i = [sigmaPoints(1,:).*cos(sigmaPoints(2,:)); sigmaPoints(1,:).*sin(sigmaPoints(2,:))];

% weighting
W0_m = (lambda1/(nx+lambda1));
Wi_m = 1/(2*(nx + lambda1));
temp2 = Wi_m*ones(1,nx*2);
W_m = [W0_m temp2];

W0_c = W0_m + 1 - alpha^2 + beta;
% Wi_c = Wi_m;
W_c = [W0_c temp2];

zbar = 0;
for i = 1:length(W_c)
    weighted_zi = W_m(i)*z_i(:,i);
    zbar = zbar + weighted_zi;
end

total = 0; 
for i = 1:length(z_i)
    temp4 = W_c(i)*(z_i(:,i)-zbar)*(z_i(:,i)-zbar)';
    total = total + temp4;
end



