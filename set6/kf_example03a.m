% kf_example03a.m
%  
%  This Matlab script loads into your work space the problem matrices
%  and data for the linear Kalman filter problem of an assignment.
%
   Fk     = [   1,     2,    2;...
                0,     1,    2;...
                0,     0,    1];                         % for all k
   Gammak = [  0.53333333333333,  1.33333333333333;...
               0.40000000000000,  2.00000000000000;...
               0.20000000000000,  2.00000000000000];     % for all k
   Hk     = [  0.5,  1.0,  0.5];                         % for all k
%
   Qk     = [ 100,  -20;...
              -20,   30];                                % for all k
   Rk     =    1;                                        % for all k
%
   xhat0   = [ 1300;...
                 15;...
                -10];
   P0      = [  20,    4,  -10;...
                 4,   30,   -3;...
               -10,   -3,   10];
%
%  Note that the following arrays define the sample times 
%  [t(1);t(2);t(3);...;t(50)] and the measurements 
%  [z(1);z(2);z(3);...;z(50)]
%
   thist = [1:50]'*2;
   zhist = 1.0e+005*[ 0.00647190560160;...
                      0.00580832801273;...
                      0.00445132359833;...
                      0.00253439244633;...
                     -0.00006278142685;...
                     -0.00332403265430;...
                     -0.00725898293689;...
                     -0.01180616838013;...
                     -0.01711976133027;...
                     -0.02302569007292;...
                     -0.02952043775875;...
                     -0.03636698683843;...
                     -0.04368628867131;...
                     -0.05167905879917;...
                     -0.06046381044327;...
                     -0.07008479240219;...
                     -0.08082039731799;...
                     -0.09241848863219;...
                     -0.10510576215007;...
                     -0.11872649827803;...
                     -0.13325609057895;...
                     -0.14873717916064;...
                     -0.16496798873743;...
                     -0.18204340978200;...
                     -0.20037744227150;...
                     -0.21982932904536;...
                     -0.24058051614791;...
                     -0.26271969768579;...
                     -0.28641763273904;...
                     -0.31198229882915;...
                     -0.33944951777224;...
                     -0.36899542646973;...
                     -0.40056657189951;...
                     -0.43407598481054;...
                     -0.46923175348785;...
                     -0.50590903055982;...
                     -0.54411147877176;...
                     -0.58389778906594;...
                     -0.62541641352699;...
                     -0.66836734281539;...
                     -0.71256546708227;...
                     -0.75796349721355;...
                     -0.80439072948236;...
                     -0.85175259478702;...
                     -0.90028360330410;...
                     -0.94993913510295;...
                     -1.00063533043694;...
                     -1.05276741145956;...
                     -1.10640222282830;...
                     -1.16157538876315];
