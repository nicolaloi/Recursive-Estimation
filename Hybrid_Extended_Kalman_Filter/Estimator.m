function [posEst,linVelEst,oriEst,windEst,driftEst,...
          posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   ut: u_t, thrust command
%                   ur: u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization

const = EstimatorConst();
if (tm == 0)
    % Do the initialization of your estimator here!
    
    % initial state mean
    posEst    = zeros(1,2); % 1x2 matrix
    linVelEst = zeros(1,2); % 1x2 matrix
    oriEst    = 0; % 1x1 matrix
    windEst   = 0; % 1x1 matrix
    driftEst  = 0; % 1x1 matrix
    
    % initial state variance
    posVar    = (const.StartRadiusBound^2)/4*ones(1,2); % 1x2 matrix https://math.stackexchange.com/questions/3992430/variance-and-covariance-of-a-circles-coordinates
    linVelVar = zeros(1,2); % 1x2 matrix
    oriVar    = (1/12)*(-const.RotationStartBound - (const.RotationStartBound))^2; % 1x1 matrix
    windVar   = (1/12)*(-const.WindAngleStartBound - (const.WindAngleStartBound))^2; % 1x1 matrix
    driftVar  = (1/12)*(-const.GyroDriftStartBound - (const.GyroDriftStartBound))^2; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar, linVelVar, oriVar, windVar, driftVar]);
    
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst]';
    
    % time of last update
    estState.tm = tm;
else
    % = = = = = = = = = = = = = = = = %
    % --------- CONSTANTS ----------- %
    % = = = = = = = = = = = = = = = = %
    
    Cdh = estConst.dragCoefficientHydr; 
    Cda = estConst.dragCoefficientAir;
    Cr  = estConst.rudderCoefficient;
    Cw  = estConst.windVel;

    xa = estConst.pos_radioA(1); ya = estConst.pos_radioA(2);
    xb = estConst.pos_radioB(1); yb = estConst.pos_radioB(2);
    xc = estConst.pos_radioC(1); yc = estConst.pos_radioC(2);
    
    ut = actuate(1); ur = actuate(2);

               %Q_d             %Q_r               %Q_rho                %Q_b
    Qc = diag([const.DragNoise, const.RudderNoise, const.WindAngleNoise, const.GyroDriftNoise]);
    
              %sigma_a^2        %sigma_b^2        sigma_c^2         %sigma_g^2       %sigma_n^2
    R = diag([const.DistNoiseA, const.DistNoiseB, const.DistNoiseC, const.GyroNoise, const.CompassNoise]);

    %% Estimator iteration.
    % get time since last estimator update
    dt = tm - estState.tm;
    estState.tm = tm; % update measurement update time

    % = = = = = = = = = = = = = = = = %
    % -------- PRIOR UPDATE --------- %
    % = = = = = = = = = = = = = = = = %

    %State dynamics
    
    sx_dot = @(x) cos(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2))-...
                  Cda*(x(3)-Cw*cos(x(6)))*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);
              
    sy_dot = @(x) sin(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2))-...
                  Cda*(x(4)-Cw*sin(x(6)))*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);     
    
    q = @(x) [     x(3);
                   x(4);
              sx_dot(x);
              sy_dot(x);
                  Cr*ur;
                      0;
                      0];


    %Variance dynamics

    A33 = @(x) -cos(x(5))*Cdh*2*x(3)-Cda*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2)...
                -Cda*(x(3)-Cw*cos(x(6)))*(x(3)-Cw*cos(x(6)))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);

    A43 = @(x) -sin(x(5))*Cdh*2*x(3)-Cda*(x(4)-Cw*sin(x(6)))*(x(3)-Cw*cos(x(6)))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);

    A34 = @(x) -cos(x(5))*Cdh*2*x(4)-Cda*(x(3)-Cw*cos(x(6)))*(x(4)-Cw*sin(x(6)))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);

    A44 = @(x) -sin(x(5))*Cdh*2*x(4)-Cda*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2)...
                -Cda*(x(4)-Cw*sin(x(6)))*(x(4)-Cw*sin(x(6)))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);

    A35 = @(x) -sin(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2));

    A45 = @(x)  cos(x(5))*(tanh(ut)-Cdh*(x(3)^2+x(4)^2));

    A36 = @(x) -Cda*Cw*sin(x(6))*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2)...
                -Cda*(x(3)-Cw*cos(x(6)))*(Cw*sin(x(6))*(x(3)-Cw*cos(x(6)))-Cw*cos(x(6))*(x(4)-Cw*sin(x(6))))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);

    A46 = @(x) Cda*Cw*cos(x(6))*sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2)...
                -Cda*(x(4)-Cw*sin(x(6)))*(Cw*sin(x(6))*(x(3)-Cw*cos(x(6)))-Cw*cos(x(6))*(x(4)-Cw*sin(x(6))))/sqrt((x(3)-Cw*cos(x(6)))^2+(x(4)-Cw*sin(x(6)))^2);
          
            
             %px    %py    %sx      %sy    %phi    %rho    %b
    A = @(x) [0,     0,     1,       0,      0,      0,     0;  % px influenced by sx
              0,     0,     0,       1,      0,      0,     0;  % py influenced by sy
              0,     0, A33(x), A34(x), A35(x), A36(x),     0;  % sx influenced by sx,sy,phi,rho
              0,     0, A43(x), A44(x), A45(x), A46(x),     0;  % sy influenced by sx,sy,phi,rho
              0,     0,      0,      0,      0,      0,     0;  % phi boat orientation not influenced by any state
              0,     0,      0,      0,      0,      0,     0;  % rho wind orientation not influenced by any state
              0,     0,      0,      0,      0,      0,     0]; % b drift not influenced by any state

    % diff sx in vd
    L31 = @(x) -Cdh*cos(x(5))*(x(3)^2+x(4)^2);
    % diff sy in vd
    L41 = @(x) -Cdh*sin(x(5))*(x(3)^2+x(4)^2);
    % diff phi in vr
    L52 = @(x) Cr*ur;

                %vd      %vr    %vrho   %vb
    L = @(x) [    0,       0,      0,     0;  % px not influenced by any disturbance
                  0,       0,      0,     0;  % py not influenced by any disturbance
              L31(x),      0,      0,     0;  % sy influenced by vd
              L41(x),      0,      0,     0;  % sy influenced by vd
                   0, L52(x),      0,     0;  % phi boat orientation influenced by vr
                   0,      0,      1,     0;  % rho wind orientation influenced by vrho
                   0,      0,      0,     1]; % b drift influenced by vb


      
    %solve ODE for x 
    [~, xp] = ode45(@(t,xp) q(xp), [0 dt], estState.xm);
    xp = xp(end,:)';
    
    %solve ODE for P
    len_Pm = size(estState.Pm,1);
    
    Pm_dot = @(x,Pp) A(x)*reshape(Pp,[len_Pm,len_Pm]) +reshape(Pp,[len_Pm,len_Pm])*A(x)' +L(x)*Qc*L(x)';

    [~, Pp_vector] = ode45(@(t,Pp_vector) reshape(Pm_dot(xp, Pp_vector),[len_Pm^2,1]), [0 dt], reshape(estState.Pm,[len_Pm^2,1]));
    Pp = reshape(Pp_vector(end,:), size(estState.Pm));
    
  
    % = = = = = = = = = = = = = = = = %
    % ----- MEASUREMENT UPDATE ------ %
    % = = = = = = = = = = = = = = = = %
    
    h = [sqrt((xp(1) - xa)^2+(xp(2) - ya)^2);  %za
         sqrt((xp(1) - xb)^2+(xp(2) - yb)^2);  %zb
         sqrt((xp(1) - xc)^2+(xp(2) - yc)^2);  %zc
                               xp(5) + xp(7);  %zg
                                      xp(5)];  %zn
     
    % diff za in x(1)
    H11 = (xp(1)-xa)/sqrt((xp(1) - xa)^2+(xp(2) - ya)^2);
    % diff za in x(2)
    H12 = (xp(2)-ya)/sqrt((xp(1) - xa)^2+(xp(2) - ya)^2);
    
    % diff zb in x(1)
    H21 = (xp(1)-xb)/sqrt((xp(1) - xb)^2+(xp(2) - yb)^2);
    % diff zb in x(2)
    H22 = (xp(2)-yb)/sqrt((xp(1) - xa)^2+(xp(2) - ya)^2);
    
    % diff zc in x(1)
    H31 = (xp(1)-xc)/sqrt((xp(1) - xc)^2+(xp(2) - yc)^2);   
    % diff zc in x(2)
    H32 = (xp(2)-yc)/sqrt((xp(1) - xa)^2+(xp(2) - ya)^2);
    
         %px  %py  %sx  %sy  %phi  %rho  %b
    H = [H11, H12,   0,   0,   0,   0,   0;  %za influenced by px,py
         H21, H22,   0,   0,   0,   0,   0;  %zb influenced by px,py
         H31, H32,   0,   0,   0,   0,   0;  %zc influenced by px,py
           0,   0,   0,   0,   1,   0,   1;  %zg influenced by phi,b 
           0,   0,   0,   0,   1,   0,   0]; %zn influenced by phi

        %wa  %wb  %wc  %wg  %wn
    M = [1,   0,   0,   0,   0;  %za influenced by wa
         0,   1,   0,   0,   0;  %zb influenced by wb
         0,   0,   1,   0,   0;  %zc influenced by wc
         0,   0,   0,   1,   0;  %zg influenced by wg
         0,   0,   0,   0,   1]; %zn influenced by wn
    
    if isinf(sense(3))
        sense(3) = [];
        h(3)     = [];
        H(3,:)   = [];
        M(3,:)   = [];
        M(:,3)   = [];
        R(3,:)   = [];
        R(:,3)   = [];
    end
    
    K = (Pp*H')/(H*Pp*H'+M*R*M');
    estState.xm = xp +K*(sense'-h);
    estState.Pm = (eye(size(K,1)) -K*H)*Pp;
    
    % Get resulting estimates and variances
    % Output quantities
    posEst = estState.xm(1:2);
    linVelEst = estState.xm(3:4);
    oriEst = estState.xm(5);
    windEst = estState.xm(6);
    driftEst = estState.xm(7);
    
    posVar = [estState.Pm(1,1),estState.Pm(2,2)];
    linVelVar = [estState.Pm(3,3),estState.Pm(4,4)];
    oriVar = estState.Pm(5,5);
    windVar = estState.Pm(6,6);
    driftVar = estState.Pm(7,7);
end
end