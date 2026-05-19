function A = jacobian_ccr4bp_H(t, x_vec, params)
%==========================================================================
%
% Computes Jacobian matrix A = df/dx for CCR4BP Hamiltonian system.
% Used in STM propagation: dot(Phi) = A*Phi.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  t         -   time                                               -
%  x_vec     -   state vector (4x1)                                 -
%  params    -   parameter struct                                   -
%
% OUTPUT:
%
%  A         -   Jacobian matrix (4x4)                              -
%==========================================================================
    
    % Inputs
    switch params.frame
        case 'synodic_m1m2'
            mu = params.mu;
            mu_3 = params.mu_pert;
            r_13 = params.r_pert;
            Omega_3 = params.Omega_pert;
            theta_3_0 = params.theta_pert_0;
            theta_3_t = theta_pert(t, params);
            
            x_1 = -mu; y_1 = 0; r_1_vec = [x_1; y_1];
            x_2 = 1-mu; y_2 = 0; r_2_vec = [x_2; y_2];

            x_3 = -mu + r_13*cos(theta_3_t);
            y_3 = r_13*sin(theta_3_t);
            r_3_vec = [x_3; y_3];

        case 'synodic_m1m3'
        otherwise
            error('Invalid frame')
    end

    % Extract values from x_vec
    x = x_vec(1);
    y = x_vec(2);
    
    % Calculate some elements
    H_xx = - (1-mu)*drm1_d2('xx', x, y, - mu, 0)...
        - mu*drm1_d2('xx', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('xx', x, y, x_3, y_3);
    
    H_yy = - (1-mu)*drm1_d2('yy', x, y, - mu, 0)...
        - mu*drm1_d2('yy', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('yy', x, y, x_3, y_3);

    H_xy = - (1-mu)*drm1_d2('xy', x, y, - mu, 0)...
        - mu*drm1_d2('xy', x, y, 1- mu, 0)...
        -mu_3*drm1_d2('xy', x, y, x_3, y_3);

    % Complete matrix A
    A = [0, 1, 1, 0;
        -1, 0, 0, 1;
        -H_xx, -H_xy, 0, 1;
        -H_xy, -H_yy, -1, 0];


end