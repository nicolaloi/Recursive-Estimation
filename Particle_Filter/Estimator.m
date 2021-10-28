function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles)
%                       corresponding to:
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles)
%                       corresponding to:
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

persistent kappa_hit

% Set number of particles:
N_particles = 11000; % obviously, you will need more particles than 10.
K_roughening = 0.01;


%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    kappa_hit = false; %flag used to activate the wall resampling
    
    phi_zero = - estConst.phi_0 + (2 * estConst.phi_0) * rand(1,N_particles);   % the initial heading phi_zero is uniformly distributed in [-pi/4, pi/4]
    kappa_zero = - estConst.l + 2 * estConst.l * rand(1,N_particles);
    
    is_in_circleA = round(rand(1,N_particles));    % 1 means that we start in circle A, otherwise we start in circle B
    % is_in_circleB = not(is_in_circleA);
    
    x_A = estConst.pA(1);
    y_A = estConst.pA(2);
    x_B = estConst.pB(1);
    y_B = estConst.pB(2);
    
    r = estConst.d * sqrt(rand(1,N_particles));
    phi = 2 * pi * rand(1,N_particles);
    
    postParticles.x_r = r .* cos(phi) + is_in_circleA * x_A + (1 - is_in_circleA) * x_B;   % 1xN_particles matrix
    postParticles.y_r = r .* sin(phi) + is_in_circleA * y_A + (1 - is_in_circleA) * y_B;   % 1xN_particles matrix
    postParticles.phi = phi_zero; % 1xN_particles matrix
    postParticles.kappa = kappa_zero; % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!


% Prior Update:
v_f = - estConst.sigma_f/2 + estConst.sigma_f * rand(1,N_particles);
v_phi = - estConst.sigma_phi/2 + estConst.sigma_phi * rand(1,N_particles);

priorParticles.x_r = prevPostParticles.x_r + (act(1) + v_f) .* cos(prevPostParticles.phi);
priorParticles.y_r = prevPostParticles.y_r + (act(1) + v_f) .* sin(prevPostParticles.phi);
priorParticles.phi = prevPostParticles.phi + act(2) + v_phi;
priorParticles.kappa = prevPostParticles.kappa;


% Posterior Update:
distance = zeros(1, N_particles);
w = zeros(1, N_particles);

contour = estConst.contour;
direction = [cos(priorParticles.phi'),sin(priorParticles.phi')];

for k = 1:N_particles
    
    contour(8,1) = prevPostParticles.kappa(k);
    contour(9,1) = prevPostParticles.kappa(k);
    [distance(k),flag] = getDistance([priorParticles.x_r(k),priorParticles.y_r(k)],direction(k,:),contour);
    
    if(flag)
        kappa_hit = true;
    end
    
    w(k) = sens - distance(k);
end

p_w = getProbability(w, estConst.epsilon);


% It can happen that all the particles have zero measurement likelihood and
% this leads to having all particle weights equal to zero.

if(sum(p_w)==0)
    
    postParticles.x_r = priorParticles.x_r;
    postParticles.y_r = priorParticles.y_r;
    postParticles.phi = priorParticles.phi;
    postParticles.kappa = priorParticles.kappa;
    
else
    
    beta = p_w./sum(p_w); %alpha;
    [postParticles.x_r, postParticles.y_r, postParticles.phi, postParticles.kappa] = resampleParticles(beta, priorParticles.x_r, priorParticles.y_r, priorParticles.phi, priorParticles.kappa, kappa_hit);
    
end

[postParticles.x_r, postParticles.y_r, postParticles.phi, postParticles.kappa] = roughening(postParticles.x_r, postParticles.y_r, postParticles.phi, postParticles.kappa, K_roughening, kappa_hit);

% END ESTIMATOR


%------------------------------------------------------%
%------------------------------------------------------%
%------------------- FUNCTIONS ------------------------%
%------------------------------------------------------%
%------------------------------------------------------%

function [dist,flag] = getDistance(x,vr,contour)
    flag = false;
    dist = 0;
    x2 = x+5*vr;

    flip = [0 -1;1 0];

    %parametrisation
    a = circshift(contour,-1,1)'- contour';%2x10
    b = flip*(x2'-x'); %2x1
    c = x'-contour'; %2x10

    %Cramers Rule
    divisor = sum(a.*b);
    s = sum(c.*b)./divisor;
    t = -sum(a.*(flip*c))./divisor;
    t_min = 42;
    for i = 1:10
        if(s(i) >= 0 && s(i) <= 1 && t(i) >= 0 && t(i) < t_min)
            dist = t(i)*5;
            t_min = t(i);
            if i == 8
                %Hit the unknown wall
                flag = true;
            end
        end
    end
end


function p_w = getProbability(w, eps)

    p_w = zeros(size(w));
    
    % Since the probability density function of the measurement noise w_k is
    % symmetric around the y-axis, we can take the absolute value

    w = abs(w);

    % from 0 to 2*eps
    p_w( (w >= 0) & (w < 2*eps) ) = -1/(5*(eps^2))*w( (w >= 0) & (w < 2*eps) ) + 2/(5*eps);        % equation of the line passing through the points (0, 2/(5*epsilon)) and (2*epsilon, 0)

    % from 2*eps to 2.5*eps
    p_w( (w >= 2*eps) & (w < 2.5*eps) ) = 2/(5*(eps^2))*w( (w >= 2*eps) & (w < 2.5*eps) ) - 4/(5*eps);   % equation of the line passing through the points (2*epsilon, 0) and (2.5*epsilon, 1/(5*epsilon))

    % from 2.5*eps to 3*eps
    p_w( (w >= 2.5*eps) & (w < 3*eps) ) = -2/(5*(eps^2))*w( (w >= 2.5*eps) & (w < 3*eps) ) + 6/(5*eps);  % equation of the line passing through the points (2.5*epsilon, 1/(5*epsilon)) and (3*epsilon, 0)

end


function [x_r, y_r, phi, kappa] = resampleParticles(beta, x_r_prior, y_r_prior, phi_prior, kappa_prior, flag)

    particles = size(x_r_prior, 2);

    beta = [0 beta];
    cumulative_sum_beta = zeros(1, particles + 1);
    summ = 0;

    for i = 1:(particles + 1)
        summ = summ + beta(i);
        cumulative_sum_beta(i) = summ;
    end

    [~, ~, bin_indices] = histcounts(rand(1,particles + 1), cumulative_sum_beta);


    x_r = zeros(1, particles);
    y_r = zeros(1, particles);
    phi = zeros(1, particles);

    if (flag)
        kappa = zeros(1, particles);
    end

    for i = 1:particles
        index_particle_k = bin_indices(i);

        x_r(i) = x_r_prior(index_particle_k);
        y_r(i) = y_r_prior(index_particle_k);
        phi(i) = phi_prior(index_particle_k);
        if (flag)
            kappa(i) = kappa_prior(index_particle_k);
        end
    end

    if (~flag)
        kappa = kappa_prior;
    end

end

function [x_r, y_r, phi, kappa] = roughening(x_r, y_r, phi, kappa, K, flag)

    d = 4;      % dimension of the state-space
    n_particles = size(x_r, 2);

    % The maximum inter-sample variability (E_i) is given by the difference
    % between the maximum and the minimum values for each state variable
    E_1 = max(x_r) - min(x_r);
    E_2 = max(y_r) - min(y_r);
    E_3 = max(phi) - min(phi);


    sigma_1 = K * E_1 * n_particles^(-1/d);
    sigma_2 = K * E_2 * n_particles^(-1/d);
    sigma_3 = K * E_3 * n_particles^(-1/d);


    x_r = x_r + sigma_1 * randn(1, n_particles);
    y_r = y_r + sigma_2 * randn(1, n_particles);
    phi = phi + sigma_3 * randn(1, n_particles);

    if (flag)
        E_4 = max(kappa) - min(kappa);
        sigma_4 = K * E_4 * n_particles^(-1/d);
        kappa = kappa + sigma_4 * randn(1, n_particles);
    end

end

end % end estimator