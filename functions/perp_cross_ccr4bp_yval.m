function [X_vec, traj, hist] = perp_cross_ccr4bp_yval(X_0_vec, t_star, params, tol, max_iter, t_0, y_val)
%==========================================================================
%
% Computes an initial condition that yields a perpendicular crossing of the
% x-axis at a prescribed time t_star in the planar CCR4BP using a
% single-shooting Newton correction.
%
% The state is expressed in Hamiltonian coordinates:
%   X = [x; y; p_x; p_y]
%
% A reduced 2x2 Newton system is used by fixing:
%
%   y_0   = 0
%   p_x,0 = 0
%
% and correcting only the free variables:
%
%   z = [x_0; p_y,0]
%
% so that
%
%   X_0(z) = [x_0; 0; 0; p_y,0]
%
% The shooting function enforces a perpendicular crossing of the x-axis:
%
%   f(z) = C * phi(t_star, X_0(z)) = 0
%
% where:
%   phi(t, X_0)  - Hamiltonian flow map of the planar CCR4BP
%   C            - selection matrix extracting [y; y + p_x]
%
% Explicitly:
%
%   f(z) = [ y(t_star; X_0(z));
%            y(t_star; X_0(z)) + p_x(t_star; X_0(z)) ] = 0
%
% Since in the standard rotating-frame Hamiltonian formulation
%
%   x_dot = p_x + y
%
% the second condition enforces:
%
%   x_dot(t_star; X_0(z)) = 0
%
% Therefore, f(z) = 0 corresponds to:
%
%   y(t_star; X_0(z))     = 0
%   x_dot(t_star; X_0(z)) = 0
%
% i.e. a perpendicular crossing of the x-axis at time t_star.
%
% The Jacobian is computed using the STM:
%
%   Df(z) = C * Phi(t_star, X_0(z)) * dX_0/dz
%
% where:
%
%   dX_0/dz = [1 0;
%              0 0;
%              0 0;
%              0 1]
%
%--------------------------------------------------------------------------
% INPUT:               Description                                   Units
%
%  X_0_vec   - initial state guess in Hamiltonian coordinates (4x1)    -
%              only x_0 and p_y,0 are used
%  t_star    - prescribed crossing time                                -
%  params    - parameter struct                                        -
%  tol       - convergence tolerance (default: 1e-10)                  -
%  max_iter  - maximum Newton iterations (default: 20)                 -
%
% OUTPUT:
%
%  X_vec            - corrected initial state [x_0; 0; 0; p_y,0]       -
%  traj             - trajectory struct from final integration         -
%  hist.z_vec       - iterate history (2 x N_iter+1)                   -
%  hist.X_vec       - initial state history (4 x N_iter+1)             -
%  hist.f_vec       - residual history (2 x N_iter)                    -
%  hist.f_norm_vec  - residual norm history (N_iter x 1)               -
%  hist.Delta_z_vec - Newton step history (2 x N_iter)                 -
%  hist.iter        - number of Newton iterations performed            -
%  hist.converged   - convergence flag                                 -
%
%--------------------------------------------------------------------------
% METHOD:
%
% At each Newton iteration:
%
%   1) Build X_0(z) = [x_0; 0; 0; p_y,0]
%   2) Integrate the CCR4BP dynamics and STM from t = 0 to t_star
%   3) Evaluate the shooting function:
%        f = C * X(t_star)
%   4) Evaluate the Jacobian:
%        Df = C * Phi(t_star) * dX_0/dz
%   5) Solve the 2x2 Newton correction:
%        Df * Delta_z = -f
%   6) Update the free variables:
%        z <- z + Delta_z
%
%--------------------------------------------------------------------------
% NOTES:
%
% - This reduced formulation enforces y_0 = 0 and p_x,0 = 0 by construction.
%
% - The second row of C assumes the Hamiltonian convention x_dot = p_x + y.
%   If a different Hamiltonian convention is used, C must be modified
%   accordingly.
%
%==========================================================================
    
    %----------------------------------------------------------------------
    % Defaults
    %----------------------------------------------------------------------
    if nargin < 4 || isempty(tol)
        tol = 1e-10;
    end

    if nargin < 5 || isempty(max_iter)
        max_iter = 20;
    end

    if nargin < 6 || isempty(t_0)
        t_0 = 0;
    end

    if nargin < 7 || isempty(y_val)
        y_val = 0;
    end

    %----------------------------------------------------------------------
    % Initialization
    %----------------------------------------------------------------------
    X_0_vec = X_0_vec(:);

    % Free variables: z = [x0; py0]
    z_vec = [X_0_vec(1); X_0_vec(4)];

    hist.z_vec = zeros(2, max_iter + 1);
    hist.X_vec = zeros(4, max_iter + 1);
    hist.f_vec = zeros(2, max_iter);
    hist.f_norm_vec = zeros(max_iter, 1);
    hist.Delta_z_vec = zeros(2, max_iter);
    hist.iter = 0;
    hist.converged = false;
    hist.z_vec(:,1) = z_vec;
    hist.X_vec(:,1) = [z_vec(1); y_val; 0; z_vec(2)];

    C = [0, 1, 0, 0;
         0, 1, 1, 0];

    dX0_dz_mtx = [1, 0;
                  0, 0;
                  0, 0;
                  0, 1];

    %----------------------------------------------------------------------
    % Newton loop
    %----------------------------------------------------------------------
    for k = 1:max_iter

        % (i) Build current initial condition from free variables
        X0_vec = [z_vec(1);
                  y_val;
                  0;
                  z_vec(2)];

        % (ii) Integrate trajectory + STM to prescribed crossing time
        traj = integrate_ccr4bp([t_0, t_0 + t_star], X0_vec, params);

        Xf_vec = traj.x_f_vec;
        Phi_f_mtx = traj.Phi_f_mtx;

        % (iii) Evaluate residual
        f_vec = [Xf_vec(2) - y_val;
                 Xf_vec(2) + Xf_vec(3)];

        % (iv) Store history
        hist.f_vec(:,k) = f_vec;
        hist.f_norm_vec(k) = norm(f_vec, 2);
        hist.iter = k;

        % (v) Convergence check
        if norm(f_vec, 2) < tol
            hist.converged = true;
            X_vec = X0_vec;

            hist.z_vec = hist.z_vec(:,1:k);
            hist.X_vec = hist.X_vec(:,1:k);
            hist.f_vec = hist.f_vec(:,1:k);
            hist.f_norm_vec = hist.f_norm_vec(1:k);
            hist.Delta_z_vec = hist.Delta_z_vec(:,1:max(k-1,1));
            return
        end

        % (vi) Evaluate Jacobian
        Df_mtx = C * Phi_f_mtx * dX0_dz_mtx;

        % (vii) Newton correction
        Delta_z_vec = -(Df_mtx \ f_vec);

        % (viii) Update
        %z_vec = z_vec + Delta_z_vec;

        %Limit step size
        
        max_step = 1e-3;

        if norm(Delta_z_vec) > max_step

            Delta_z_vec = Delta_z_vec * max_step/norm(Delta_z_vec);

        end

        % Damped Newton line search

        alpha = 1.0;

        alpha_min = 1e-5;

        f_norm_old = norm(f_vec, 2);

        while alpha > alpha_min

            z_trial_vec = z_vec + alpha*Delta_z_vec;

            X0_trial_vec = [z_trial_vec(1);

                            y_val;

                            0;

                            z_trial_vec(2)];

            traj_trial = integrate_ccr4bp([t_0, t_0 + t_star], X0_trial_vec, params);

            f_trial_vec = [traj_trial.x_f_vec(2) - y_val;
                           traj_trial.x_f_vec(2) + traj_trial.x_f_vec(3)];

            f_norm_trial = norm(f_trial_vec, 2);

            if f_norm_trial < f_norm_old

                break

            end

            alpha = alpha/2;

        end

        if alpha <= alpha_min

            warning('Line search failed at iteration %d.', k)

            break

        end

        % Update with damped step

        z_vec = z_vec + alpha*Delta_z_vec;

        % Save history
        hist.Delta_z_vec(:,k) = Delta_z_vec;
        hist.z_vec(:,k+1) = z_vec;
        hist.X_vec(:,k+1) = [z_vec(1); y_val; 0; z_vec(2)];
    end

    %----------------------------------------------------------------------
    % No convergence
    %----------------------------------------------------------------------
    X_vec = [z_vec(1); y_val; 0; z_vec(2)];

    hist.z_vec = hist.z_vec(:,1:max_iter+1);
    hist.X_vec = hist.X_vec(:,1:max_iter+1);
    hist.f_vec = hist.f_vec(:,1:max_iter);
    hist.f_norm_vec = hist.f_norm_vec(1:max_iter);
    hist.Delta_z_vec = hist.Delta_z_vec(:,1:max_iter);

    warning('perp_cross_ccr4bp did not converge in %d iterations.', max_iter)

end