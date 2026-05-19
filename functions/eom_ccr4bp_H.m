function dx_dt_vec =  eom_ccr4bp_H(t, x_vec, params)
%==========================================================================
%
% Computes planar CCR4BP Hamiltonian equations of motion in synodic frame.
% State vector is defined as x = [x y p_x p_y]^T.
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
%  dx_dt_vec -   time derivative of state (4x1)                     -
%==========================================================================

    % Initialization
    x_vec = x_vec(:);
    dx_dt_vec = zeros(size(x_vec));

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
    p_x = x_vec(3);
    p_y = x_vec(4);

    % Calculate important quantities
    r_vec = [x; y];
    r_14 = norm(r_vec-r_1_vec);
    r_24 = norm(r_vec-r_2_vec);
    r_34 = norm(r_vec-r_3_vec);


    % Hamilton Eqs
    dx_dt_vec(1) = p_x + y;
    dx_dt_vec(2) = p_y - x;
    dx_dt_vec(3) = p_y - (1-mu)*(x-x_1)/r_14^3 - mu*(x - x_2)/r_24^3 ...
    -mu_3*(x - x_3)/r_34^3 - mu_3*cos(theta_3_t)/r_13^2;
    dx_dt_vec(4) = -p_x - (1-mu)*(y-y_1)/r_14^3 - mu*(y - y_2)/r_24^3 ...
    -mu_3*(y - y_3)/r_34^3 - mu_3*sin(theta_3_t)/r_13^2;
    
end