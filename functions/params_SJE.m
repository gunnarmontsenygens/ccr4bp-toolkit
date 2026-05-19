function params = params_SJE()
%==========================================================================
%
% Defines parameters for the Saturn–Janus–Epimetheus CCR4BP model.
% Synodic frame is defined with respect to the primary pair (m1,m2).
% 1 - Saturn
% 2 - Europa
% 3  - Ganymede
%
% Author: G. Montseny
% Date: April 30, 2026
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
    
    params.m1 = 5.683e26; % kg
    params.m2 = 1.98e18; % kg
    params.m3 = 5.6e17; % kg 

    params.mu    = params.m2/(params.m1 + params.m2);   % m2 / (m1 + m2)
    params.mu_pert  = params.m3/(params.m1 + params.m2);   % m3 / (m1 + m2)
    
    %======================================================================
    % Geometry (normalized distances)
    %======================================================================
    
    params.r_SJ = 151472; % in km
    params.r_SE = 151422; % in km
    params.r_pert  = params.r_SE/params.r_SJ;   % distance of perturbing body orbit
    
    %======================================================================
    % Perturbing body motion
    %======================================================================
    %params.T2 = 3.551181; % in days
    %params.T3 = 7.15455296; % in days
    params.Omega_pert    = 1; % angular rate in synodic frame
    params.theta_pert_0  = pi;      % initial phase
    
    %======================================================================
    % Numerical settings
    %======================================================================
    
    params.ode.options = odeset( ...
        'RelTol', 1e-13, ...
        'AbsTol', 1e-13);
    
end