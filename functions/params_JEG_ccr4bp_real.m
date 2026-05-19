function params = params_JEG_ccr4bp_real()
%==========================================================================
%
% Defines parameters for the Jupiter–Europa–Ganymede CCR4BP model.
% Synodic frame is defined with respect to the primary pair (m1,m2).
% 1 - Jupiter
% 2 - Europa
% 3  - Ganymede
%
% Author: G. Montseny
% Date: April 09, 2026
%

%
% OUTPUT:
%
%  params     -   parameter struct for CCR4BP model                  -
%==========================================================================
    
    %======================================================================
    % Frame definition
    %======================================================================
    
    params.frame = 'synodic_m1m2';     % Synodic frame (Jupiter–Europa)
    
    %======================================================================
    % Mass parameters (normalized)
    %======================================================================
    
    params.m1 = 1.89813e27; % kg
    params.m2 = 4.79984e22; % kg
    params.m3 = 1.4819e23; % kg 

    params.mu    = params.m2/(params.m1 + params.m2);   % m2 / (m1 + m2)
    params.mu_pert  = params.m3/(params.m1 + params.m2);   % m3 / (m1 + m2)
    
    %======================================================================
    % Geometry (normalized distances)
    %======================================================================
    
    params.r_JE = 670900; % in km
    params.r_JG = 1070400; % in km
    params.r_pert  = params.r_JG/params.r_JE;   % distance of perturbing body orbit
    
    %======================================================================
    % Perturbing body motion
    %======================================================================
    params.T2 = 3.551181; % in days
    params.T3 = 7.15455296; % in days
    params.Omega_pert    = params.T2/params.T3;      % angular rate in synodic frame
    params.theta_pert_0  = 0;      % initial phase
    
    %======================================================================
    % Numerical settings
    %======================================================================
    
    params.ode.options = odeset( ...
        'RelTol', 1e-13, ...
        'AbsTol', 1e-13);
    
end