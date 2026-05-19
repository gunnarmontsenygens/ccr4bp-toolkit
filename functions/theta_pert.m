function theta_t = theta_pert(t, params)
%==========================================================================
%
% Computes angular position of perturbing body in synodic frame.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  t         -   time                                               -
%  params    -   parameter struct                                   -
%
% OUTPUT:
%
%  theta_t   -   angular position (scalar)                          -
%==========================================================================
    
    % Inputs
    Omega_pert = params.Omega_pert;
    theta_pert_0 = params.theta_pert_0;

    % Formula
    theta_t = (Omega_pert - 1)*t + theta_pert_0;

end