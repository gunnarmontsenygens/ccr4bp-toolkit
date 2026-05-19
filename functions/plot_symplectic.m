function plot_symplectic(traj)
%==========================================================================
%
% Plots determinant and symplecticity error of STM along trajectory.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  traj      -   trajectory struct containing symplecticity results   -
%
% OUTPUT:
%
%  None
%==========================================================================
    
    % Extract variables
    t_vec = traj.t_vec;
    det_vec = traj.det_vec;
    sympl_err_vec = traj.sympl_err_vec;

    % Plot determinant variation
    % figure();
    % plot(t_vec, det_vec - mean(det_vec));
    % grid on;
    % xlabel('$t$','Interpreter','latex');
    % ylabel('$\Delta \det(\Phi(t))$','Interpreter','latex');

    % Plot symplecticity error
    figure();
    plot(t_vec, sympl_err_vec);
    grid on;
    xlabel('$t$','Interpreter','latex');
    ylabel('$\|\Phi(t)^T J \Phi(t)-J\|_{F}$','Interpreter','latex');

end