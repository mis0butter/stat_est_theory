function [xhathist,Phist,sigmahist,enuhist] = ...
                    pll_iekf(delt,zhist,xhat0,P0,qC,qA,sigma,Niter)
%
%  Copyright (c) 2002 Mark L. Psiaki.  All rights reserved.  
% 
%  This function performs iterated extended Kalman filtering for the
%  phase-locked loop problem.  It uses both the in-phase and the
%  quadrature base-band signals in its nonlinear measurement update,
%  and it estimates the carrier amplitude in addition to the
%  carrier phase states.
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
%    xhat0       The 4x1 initial state estimate.  The first element is
%                the initial carrier phase estimate in radians, the 
%                second element is the initial Doppler shift estimate in 
%                rad/sec, the third element is the initial Doppler shift
%                rate estimate in rad/sec^2, and the fourth element
%                is the initial carrier amplitude estimate in receiver
%                output amplitude units.
%
%    P0          The 4x4 symmetric positive definite initial state
%                estimation error covariance matrix.
%
%    qC          The white-noise intensity of the carrier phase jerk
%                (rate of change of Doppler shift rate) in units
%                of rad^2/sec^5.
%
%    qA          The white-noise intensity of the amplitude rate
%                in receiver output amplitude units^2/sec.
%
%    sigma       The receiver measurement noise standard deviation in
%                receiver output amplitude units.
%    
%    Niter       The number of Gauss-Newton interations to do during
%                the nonlinear measurement update.  Niter = 1 will
%                yield the same results as the function pll_ekf.m.
%  
%  Outputs:
%
%    xhathist    The (k+1)x4 array that contains the time history of the
%                state vector estimates.  xhathist(j,:)' is the 
%                PLL state vector estimate at time j*delt.
%
%    Phist       The 4x4x(k+1) array that contains the time history
%                of the estimation error covariance matrices.
%                Phist(:,:,j) is the 4x4 estimation error covariance
%                matrix for the estimate xhathist(j,:)'.
%
%    sigmahist   The (k+1)x4 array that contains the filter standard
%                deviation time histories.  sigmahist(j,:) = ...
%                sqrt(diag(Phist(:,:,j)))'.
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
   xhathist = zeros(kp1,4);
   Phist = zeros(4,4,kp1);
   sigmahist = zeros(kp1,4);
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
   Qrootdum = diag([(delt^2);delt;1]*sqrt(qC*delt))*Qrootdum;
   Qroot = zeros(4,4);
   Qroot(1:3,1:3) = Qrootdum;
   Qroot(4,4) = sqrt(qA*delt);
   Rvv = inv(Qroot);
%
%  Set up the inverse state transition matrix and its inverse multiplied
%  by Gamma, which equals the inverse of the state transition matrix
%  in this case because Gamma equals the identity matrix.
%
   Finv = eye(4);
   Finv(1,2) = -delt;
   Finv(1,3) = 0.5*(delt^2);
   Finv(2,3) = -delt;
   FinvGamma = Finv; 
%
%  This loop performs one SRIF propagation step and one extended 
%  measurement update step per iteration.
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
%  Measurement update.
%
      xlinearizej = ????;
      for iiter = 1:Niter
         cosphij = cos(xlinearizej(1,1));
         sinphij = sin(xlinearizej(1,1));
         Aj = xlinearizej(4,1);
         Ha = ????;
         delza = ????;
         ????;
         ????;
         ????;
         ????;
         Rscriptj = ????;
         zscriptj = ????;
         Rscriptjinv = inv(Rscriptj);
         xhatj = Rscriptjinv*zscriptj;
         xlinearizej = ????;
      end
      jp1 = j + 1;
      xhathist(jp1,:) = xhatj';
      Pj = Rscriptjinv*(Rscriptjinv');
      Phist(:,:,jp1) = Pj;
      sigmahist(jp1,:) = sqrt(diag(Pj))';
      zscriptrj = zdum(5:6,1);
      enuhist(jp1,1) = ????;
   end