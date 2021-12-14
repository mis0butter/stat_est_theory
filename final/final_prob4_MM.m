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

end
 
%----- Truth model
% switchHist gives the model switching time history and switchIndex gives
% the indices at and after which each model obtains.  For the static case,
% switchHist is either 1, 2, or 3, and switchIndex = 1;
% switchHist = [3];
% switchIndex = [1];

switchHist = [1 3 2 3 3];
% switchHist = [1 1 1 1 1]; 
% switchHist = [1 3 2 1 3]; 
% switchHist = [1 2 1 2 3]; 
switchIndex = [1, 201, 401, 601, 801]; 

Qk = 0.001*diag([0.1 1]); Rq = chol(Qk);
vtruekhist = (Rq'*randn(nx,Nsim))';
Rkp1 = 0.1; Rr = chol(Rkp1);
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
P1 = 10*eye(nx);

N_mu = 100; 
mu_lb_vec = linspace(1e-10, 1e-2, N_mu)'; 

% N_mu = 1; 
% mu_lb_vec = 1e-5; 
    
% select alpha 
a = 0.01; 

% compute bounds 
r1 = chi2inv(a/2, (Nsim-2) * nx) / (Nsim-2); 
r2 = chi2inv( 1 - a/2, (Nsim-2) * nx) / (Nsim-2); 

err_nu_mean_arr = []; 

for i = 1 : N_mu 
    
    mu_lb = mu_lb_vec(i); 
    
    % MM filter 
    [ xhat_hist, P_hist, mu_hist ] = MM_filter( P1, nx, x1, nmod, Nsim, ... 
        ztrue_hist, utruekhist, F_M, G_M, Qk, Hkp1, Rkp1, mu_lb) ; 

    % chi-squared test 
    xtilde = (xtrue_hist - xhat_hist); 

    % compute err nu mean 
    for i = 1:length(xtilde)
        err_nu(i,:) = xtilde(i,:) * inv(P_hist(:,:,i)) * xtilde(i,:)'; 
    end 
%     err_nu = err_nu(~isnan(err_nu)); 
    err_nu_mean = mean(err_nu(2:end-1)); 
    
    if r1 < err_nu_mean && err_nu_mean < r2 
        disp('found consistent mu') 
        mu_lb 
        err_nu_mean 
        r1
        r2 
        break 
    end 
    
    err_nu_mean_arr = [err_nu_mean_arr; err_nu_mean]; 

end 

%% 

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
    plot(tkhist(iidum), sqrt(squeeze(P_hist(1,1,iidum))), 'k');
    plot(tkhist(iidum), -sqrt(squeeze(P_hist(1,1,iidum))), 'k');
    ylabel('\Delta \theta_z (rad)');
    title('Estimation errors and covariances');

    subplot(212)
    iidum = 1:Nsim-1;
    plot(tkhist(iidum), xtrue_hist(iidum,2) - xhat_hist(iidum,2));
    hold on;
    plot(tkhist(iidum), sqrt(squeeze(P_hist(2,2,iidum))), 'k');
    plot(tkhist(iidum), -sqrt(squeeze(P_hist(2,2,iidum))), 'k');
    ylabel('\Delta \omega_z (rad/s)');
    xlabel('Time (s)');

%% Monte carlo chi-squared tests 

% for monte carlo multiple realizations: 
%   nu_bar_k = sum_epsilon i=1 to n ( xtilde(k) Pxx(k)^-1 & xtilde(k) )
% what do you prize: accuracy or consistency? 
% xtilde(k)' * P^-1(k) * xtilde(k) = nx 

% prize consistency 
% if mu lower bound is too high, then the filter is not confident in its 
% estimate of the model, and thus is not efficient (and also biased). 
% you want the lowest lower bound that gives you consistency 
% filter is consistent if the the filter is unbiased and efficient 

%% subfunctions 


function [ xhat_hist, P_hist, mu_hist ] = MM_filter( P1, nx, x1, nmod, Nsim, ... 
    ztrue_hist, utruekhist, F_M, G_M, Qk, Hkp1, Rkp1, mu_lowerbound) 

%     xhat1 = x1 + (chol(P1))'*randn(nx,1);
    xhat1 = x1 + chol(P1) * [ 1.30562310353041
         0.983969531303835 ]; 
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
    %     Lambdakp1j = normpdf(nukp1j, 0, Skp1j); 
        Lambdakp1j = mvnpdf(nukp1j, 0, Skp1j); 

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

    %   % implement lower bound 
    %   for ind = 1:length(mukp1vec) 
    %       if mukp1vec(ind) < mu_lowerbound 
    %           mukp1vec(ind) = mu_lowerbound; 
    %       end 
    %   end 
    %   mukp1vec = mukp1vec / sum(mukp1vec); 

    %% run lower bound thingy twice 
        % lower bound logical mask
%         lb_mask = mukp1vec <= mu_lowerbound;
%         % above bound mask
%         ab_mask = ~lb_mask;
%         % normalization constraint
%         constraint = 1 - mu_lowerbound*sum(lb_mask);
%         if sum(lb_mask) > 0
%             mukp1vec(lb_mask) = mu_lowerbound;
%             mukp1vec(ab_mask) = mukp1vec(ab_mask)*constraint/sum(mukp1vec(ab_mask));
%             sanityCheck = sum(mukp1vec(ab_mask));
%         end
%         % lower bound logical mask
%         lb_mask = mukp1vec <= mu_lowerbound; 
%         % above bound mask
%         ab_mask = ~lb_mask;
%         % normalization constraint
%         constraint = 1 - mu_lowerbound*sum(lb_mask);
%         if sum(lb_mask) > 0
%             mukp1vec(lb_mask) = mu_lowerbound;
%             mukp1vec(ab_mask) = mukp1vec(ab_mask)*constraint/sum(mukp1vec(ab_mask));
%             sanityCheck = sum(mukp1vec(ab_mask));
%         end
   
    modeSwitchingFlag = true; 
    % enable model switching behavior
    if modeSwitchingFlag == true
        % lower bound logical mask
        [~, pos_bestState] = max(mukp1vec);
        bestState = xhat_Mhist(kp1,:,pos_bestState);
        lb_mask = mukp1vec <= mu_lowerbound;
        % above bound mask
        ab_mask = ~lb_mask;
        % normalization constraint
        constraint = 1 - mu_lowerbound*sum(lb_mask);
        % check bleeding twice
        if sum(lb_mask) > 0
            mukp1vec(lb_mask) = mu_lowerbound;
            mukp1vec(ab_mask) = mukp1vec(ab_mask)*constraint/sum(mukp1vec(ab_mask));
            %            sanityCheck = sum(mukp1vec(ab_mask));
        end
        % check again to prevent bleeding which occurs when the model probabilties after adjustment fall below threshold
        lb_mask2 = mukp1vec <= mu_lowerbound;
        ab_mask2 = ~lb_mask2;
        constraint2 = 1 - mu_lowerbound*sum(lb_mask2);
        if sum(lb_mask2) > 0
            mukp1vec(lb_mask2) = mu_lowerbound;
            mukp1vec(ab_mask2) = mukp1vec(ab_mask2)*constraint/sum(mukp1vec(ab_mask2));
        end
        % overwite state for model probabilites that are below threshold
        % find best state
        indexVec = 1:nmod;
        pos2replace = indexVec(lb_mask2);
        for i = 1:sum(lb_mask2)
            xhatk_hist(kp1,:,pos2replace(i)) = bestState;
        end
    end



%% 

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
      P_hist(:,:,kp1) = Pkp1;  
    end

end 

  
  


