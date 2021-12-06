% mmExample_temp.m
% Multiple-model estimation example: The solar panel deployment problem

clear;clc;

%% Simulation parameters

dt = 0.1;
Nsim = 1000;
tkhist = [0:Nsim-1]'*dt;
% nmod is the number of models considered
nmod = 3;  
nx = 2;
nz = 1;
nu = 1;

%----- Random number seed
n_seed = round(sum(clock*100));
rng(n_seed);
rng(0); % set rng to 0 for now 

%----- Storage matrices
xtrue_hist = zeros(Nsim,nx);
ztrue_hist = zeros(Nsim,nz);
xhat_hist = zeros(Nsim,nx);
P_hist = zeros(Nsim,nx);    
% xhat Mhist is an array of all the mode-conditioned estimates
xhat_Mhist = zeros(Nsim,nx,nmod); 
P_Mhist = zeros(Nsim,nx,nx,nmod);
mu_hist = zeros(Nsim,nmod);
nu_Mhist = zeros(Nsim,nz,nmod);
Lambda_Mhist = zeros(Nsim,nmod);
F_M = zeros(nx,nx,nmod);
G_M = zeros(nx,nu,nmod);

%----- Set up multiple system models
cSet = [0.1;0.5;1];  hzSet = [1;2;3];
avec = cSet./hzSet;
bvec = 1./hzSet;
for j=1:nmod
    a = avec(j);
    b = bvec(j);
    % Discrete-time state transition matrix
    A = [0 1; 0 -a]; 
    B = [0; b]; 
    F_M(:,:,j) = expm(A*dt); 
    % Discrete-time control matrix

    % integral function 
    fcn = @(sigma) expm(A*sigma) * B; 
    Gk = integral(fcn, 0, dt, 'ArrayValued', true); 
    
    G_M(:,:,j) = Gk; 

    % using ode45 
%     options = odeset('reltol',1e-10); 
%     G0 = reshape(eye(2), [], 1); 
%     G0 = [0; 0]; 
%     [tvec, Gk] = ode45(@ode_Gk, [0 dt], G0, options );

end
 
%----- Truth model
% switchHist gives the model switching time history and switchIndex gives
% the indices at and after which each model obtains.  For the static case,
% switchHist is either 1, 2, or 3, and switchIndex = 1;
% switchHist = [3];
% switchIndex = [1];

switchHist = [1, 2, 3, 2, 3 ]; 
switchIndex = [1, 201, 401, 601, 801]; 
mu_lowerbound = 1e-8; 

Qk = 0.001*diag([0.1 1]); Rq = chol(Qk);
vtruekhist = (Rq'*randn(nx,Nsim))';
Rkp1 = 0.15; Rr = chol(Rkp1);
wtruekhist = (Rr'*randn(nz,Nsim))';
utruekhist = 2*randn(Nsim,nu);
Hkp1 = [1 0];
x1 = [0;0.1];

%----- Generate truth-model states and measurements
xk = x1;
ii = 1;
for k=1:Nsim-1
  if k==switchIndex(ii)
    Fk = F_M(:,:,switchHist(ii));
    Gk = G_M(:,:,switchHist(ii));
    if ii<length(switchHist)
      ii = ii+1;
    end
  end
  kp1 = k + 1;
  uk = utruekhist(k,:)';
  vk = vtruekhist(k,:)';
  wkp1 = wtruekhist(kp1,:)';
  xkp1 = Fk*xk + Gk*uk + vk;
  zkp1 = Hkp1*xkp1 + wkp1;
  xtrue_hist(k,:) = xk';
  ztrue_hist(kp1,:) = zkp1';
  xk = xkp1;
end

%----- Run the Multiple-model filter
% Set up initial states and covariances
P1 = 100*eye(nx);

xhat1 = x1 + (chol(P1))'*randn(nx,1);
for j=1:nmod
  P_Mhist(1,:,:,j) = P1;
  % Initialize the mode probabilities as equally probable
  mu_hist(1,j) = 1/nmod; 
  xhat_Mhist(1,:,j) = xhat1';
end
for k=1:Nsim-1
  kp1 = k + 1;
  zkp1 = ztrue_hist(kp1,:)'; 
  uk = utruekhist(k,:)';
  for j=1:nmod
    % Propagation step
    Fkj = F_M(:,:,j);
    Gkj = G_M(:,:,j);
    xhatkj = xhat_Mhist(k,:,j)';
    Pkj = squeeze(P_Mhist(k,:,:,j));
    xbarkp1j = Fkj * xhatkj + Gkj * uk; 
    Pbarkp1j = Fkj * Pkj * Fkj' + Qk; 
    
    % Measurement update
    nukp1j = zkp1 - Hkp1*xbarkp1j;         
    Skp1j = Hkp1 * Pbarkp1j * Hkp1' + Rkp1; 
    invSkp1j = inv(Skp1j); 
    Wkp1j = Pbarkp1j * Hkp1' * invSkp1j; 
    xhatkp1j = xbarkp1j + Wkp1j * nukp1j; 
    Pkp1j = Pbarkp1j - Wkp1j * Skp1j * Wkp1j'; 
    
    % Calculate the likelihood of the current innovation
    Lambdakp1j = normpdf(nukp1j, 0, Skp1j); 
    
    % Store state estimate and covariance, innovation, and likelihood
    xhat_Mhist(kp1,:,j) = xhatkp1j';
    P_Mhist(kp1,:,:,j) = Pkp1j;
    nu_Mhist(kp1,:,j) = nukp1j';
    Lambda_Mhist(kp1,j) = Lambdakp1j;
  end
  
  % Update the mode probabilities
  mukvec = mu_hist(k,:)';
  Lambdakp1vec = Lambda_Mhist(kp1,:)';
  % mukp1vec is the nmod-by-1 vector that contains the mode probabilities
  % for the nmod modes at time index kp1; it must satisfy sum(mukp1vec) = 1.
  mukp1vec = mukvec .* Lambdakp1vec / dot(mukvec, Lambdakp1vec); 
  
  % implement lower bound 
  for ind = 1:length(mukp1vec) 
      if mukp1vec(ind) < mu_lowerbound 
          mukp1vec(ind) = mu_lowerbound; 
      end 
  end 
  
  mukp1vec = mukp1vec / sum(mukp1vec); 
  
  mu_hist(kp1,:) = mukp1vec';
  
  % Calculate the combined state estimate and covariance
  xMdum = reshape(squeeze(xhat_Mhist(kp1,:,:)),nx,nmod);
  mudum = mu_hist(kp1,:)';
  % Calculate the MMSE estimate averaged over all the modes
  xhatkp1 = (xMdum*mudum);
  % Calculate the estimation error covariance of xhatkp1
  Pkp1 = zeros(nx,nx);
  for j=1:nmod
    Pkp1j = squeeze(P_Mhist(kp1,:,:,j));
    xhatkp1j = xhat_Mhist(kp1,:,j)';
    mukp1j = mu_hist(kp1,j);
    Pkp1 = Pkp1 + mukp1j * (Pkp1j + (xhatkp1j - xhatkp1) * (xhatkp1j - xhatkp1)'); 
  end
  xhat_hist(kp1,:) = xhatkp1'; 
  P_hist(kp1,:) = diag(Pkp1)';  
end

%----- Display results
figure(1);clf;
plot(tkhist,mu_hist);shg;
xlabel('Time (s)');
ylabel('Probability');
legend('\mu_1', '\mu_2','\mu_3');
title('Mode probability time histories');
ylim([0 1.1]);

figure(2);clf;
subplot(211)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtrue_hist(iidum,1) - xhat_hist(iidum,1));
hold on;
plot(tkhist(iidum), sqrt(P_hist(iidum,1)), 'k');
plot(tkhist(iidum), -sqrt(P_hist(iidum,1)), 'k');
ylabel('\Delta \theta_z (rad)');
title('Estimation errors and covariances');

subplot(212)
iidum = 1:Nsim-1;
plot(tkhist(iidum), xtrue_hist(iidum,2) - xhat_hist(iidum,2));
hold on;
plot(tkhist(iidum), sqrt(P_hist(iidum,2)), 'k');
plot(tkhist(iidum), -sqrt(P_hist(iidum,2)), 'k');
ylabel('\Delta \omega_z (rad/s)');
xlabel('Time (s)');

%% subfunctions 

function xdot = ode_Gk(t, x) 

    global A B 

%     x = reshape(x, [2 2]); 

    xdot = expm(A*t) * B;
    
%     xdot = reshape(xdot, [], 1); 

end 

  
  


