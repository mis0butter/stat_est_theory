function [xhathist,Phist,sigmahist,nhist,enuhist] = ...
                    pll_lkf(delt,zhist,xhat0,P0,qC,sigma)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function performs linear Kalman filtering for the
%  phase-locked loop problem.  The linear Kalman filtering assumes
%  that the phase angle can be "measured" by taking the 2-argument
%  arctangent of the in-phase and quadrature measurements.  It does
%  the filtering using a square-root information filter format.
%  
%
%  Inputs:
%
%    delt        The time in seconds between RF front-end output samples.
%
%    zhist       = [I_1,Q_1;I_2,Q_2;I_3,Q_3;...;...
%                I_k,Q_k], the kx2 time history of the simulated
%                in-phase and quadrature base-band mixed signals.
%
%    xhat0       The 3x1 initial state estimate.  The first element is
%                the initial carrier phase estimate in radians, the 
%                second element is the initial Doppler shift estimate in 
%                rad/sec, and the third element is the initial Doppler 
%                shift rate estimate in rad/sec^2.
%
%    P0          The 3x3 symmetric positive definite initial state
%                estimation error covariance matrix.
%
%    qC          The white-noise intensity of the carrier phase jerk
%                (rate of change of Doppler shift rate) in units
%                of rad^2/sec^5.
%
%    sigma       The receiver measurement noise standard deviation in
%                receiver output amplitude units.
%  
%  Outputs:
%
%    xhathist    The (k+1)x3 array that contains the time history of the
%                state vector estimates.  xhathist(j,:)' is the 
%                PLL state vector estimate at time j*delt.
%
%    Phist       The 3x3x(k+1) array that contains the time history
%                of the estimation error covariance matrices.
%                Phist(:,:,j) is the 3x3 estimation error covariance
%                matrix for the estimate xhathist(j,:)'.
%
%    sigmahist   The (k+1)x3 array that contains the filter standard
%                deviation time histories.  sigmahist(j,:) = ...
%                sqrt(diag(Phist(:,:,j)))'.
%
%    nhist       The (k+1)x1 vector of the number of cycles needed to
%                "unwind" the 2*pi phase ambiguity of the 2-argument
%                arctangent function.  The first value is zero because
%                there is no measurement update at the first sample 
%                time.
%
%    enuhist     The (k+1)x1 vector of the normalized innovation
%                statistic.  The first value is zero because there is no
%                measurement update at the first sample time.
%

%
%  Get the problem dimensions and initialize the output arrays.
%
   k = size(zhist,1);
   kp1 = k + 1;
   xhathist = zeros(kp1,3);
   Phist = zeros(3,3,kp1);
   sigmahist = zeros(kp1,3);
   nhist = zeros(kp1,1);
   enuhist = zeros(kp1,1);
   xhathist(1,:) = xhat0';
   Phist(:,:,1) = P0;
   sigmahist(1,:) = sqrt(diag(P0))';
%
%  Initialize the state information equation for the SRIF
%
   Rscriptj = ????;
   zscriptj = ????;
%
%  Determine the a priori square root information matrix for the
%  process noise.
%
   Qrootdum = chol([ (1/20), (1/8), (1/6);...
                      (1/8), (1/3), (1/2);...
                      (1/6), (1/2),     1])';
   Qroot = diag([(delt^2);delt;1]*sqrt(qC*delt))*Qrootdum;
   Rww = inv(Qroot);
%
%  Set up the inverse state transition matrix and its inverse multiplied
%  by Gamma, which equals the inverse of the state transition matrix
%  in this case because Gamma equals the identity matrix.
%
   Finv = eye(3);
   Finv(1,2) = -delt;
   Finv(1,3) = 0.5*(delt^2);
   Finv(2,3) = -delt;
   FinvGamma = Finv; 
%
%  Determine the carrier phase measurements and their approximate
%  noise standard deviation.
%
   zphihist = atan2(zhist(:,2),zhist(:,1));
   Aplusnoisesqhist = sum(zhist.^2,2);
   Asqmean = mean(Aplusnoisesqhist) - 2*(sigma^2);
   sigmaphi = sigma/sqrt(Asqmean);
%
%  This loop performs one SRIF propagation step and one measurement
%  update step per iteration.
%
   for j = 1:k
      zscriptjm1 = zscriptj;
      Rscriptjm1 = Rscriptj;
%
%  Propagation from stage j-1 to stage j.
%
      ????;
      ????;
      ????;
      ????;
      Rscripttilj = ????;
      zscripttilj = ????;
%
%  Ad hoc determination of nj.
%
      xtilj = Rscripttilj\zscripttilj;
      zphij = zphihist(j,1);
      nj = ????;
      jp1 = j + 1;
      nhist(jp1,1) = nj;
%
%  Measurement update.
%
      ????;
      ????;
      ????;
      ????;
      Rscriptj = ????;
      zscriptj = ????;
      Rscriptjinv = inv(Rscriptj);
      xhatj = Rscriptjinv*zscriptj;
      xhathist(jp1,:) = xhatj';
      Pj = Rscriptjinv*(Rscriptjinv');
      Phist(:,:,jp1) = Pj;
      sigmahist(jp1,:) = sqrt(diag(Pj))';
      enuhist(jp1,1) = ????;
   end