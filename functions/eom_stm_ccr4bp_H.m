function dY_dt_vec = eom_stm_ccr4bp_H(t, Y_vec, params)
%==========================================================================
%
% Computes CCR4BP Hamiltonian equations of motion together with STM dynamics.
% Augmented state is defined as Y = [x; vec(Phi)], where x is 4x1 and Phi is 4x4.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  t         -   time                                               -
%  Y_vec     -   augmented state (20x1)                             -
%                 Y_vec = [x(4x1); vec(Phi)(16x1)]
%  params    -   parameter struct                                   -
%
% OUTPUT:
%
%  dY_dt_vec -   time derivative of augmented state (20x1)          -
%==========================================================================

    % Initialization
    Y_vec = Y_vec(:);
    
    % Extract state and STM vector
    x_vec = Y_vec(1:4);
    Phi_vec = Y_vec(5:20);

    % EoM
    dx_dt_vec = eom_ccr4bp_H(t, x_vec, params);

    % STM
    Phi_mtx = reshape(Phi_vec, 4, 4);
    A_t = jacobian_ccr4bp_H(t, x_vec, params);
    dPhi_dt_mtx = A_t*Phi_mtx;

    % Put vectors back into Y
    dx_dt_vec = dx_dt_vec(:);
    dPhi_dt_vec = dPhi_dt_mtx(:);
    dY_dt_vec = [dx_dt_vec; dPhi_dt_vec];
end