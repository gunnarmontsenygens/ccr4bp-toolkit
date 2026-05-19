function [z_vec, traj, hist] = single_shooting_ccr4bp(z_0_vec, params, tol, max_iter)
%==========================================================================
%
% Computes a periodic orbit in the planar CCR4BP using single-shooting
% Newton correction (lecture-notes formulation).
%
% State is expressed in Hamiltonian coordinates:
%   X = [x; y; px; py]
%
% Unknown vector:
%   z = [x0; y0; px0; py0; T]
%
%--------------------------------------------------------------------------
% EQUATIONS USED
%--------------------------------------------------------------------------
%
% Periodicity condition:
%
%   f(z) = phi(T, X0) - X0 = 0
%
% Phase constraint:
%
%   g(z) = y0 = 0
%
% Full residual:
%
%   h(z) = [ phi(T, X0) - X0 ;
%            y0               ] = 0
%
% Jacobian:
%
%   Dh(z) = [ Phi(T,0) - I    |   Xdot(T) ;
%             0  1  0  0       |     0     ]
%
% where:
%
%   Phi(T,0) = State Transition Matrix (STM)
%   Xdot(T)  = f(T, X(T)) from equations of motion
%
% Newton step:
%
%   Dh(z) * Delta_z = -h(z)
%
% Update:
%
%   z_{k+1} = z_k + Delta_z
%
%--------------------------------------------------------------------------
%
% Author: G. Montseny
%==========================================================================
%
%
% INPUT:               Description                                   Units
%
%  z_0_vec    -   initial guess [x0; y0; px0; py0; T]               -
%  params     -   parameter struct                                 -
%  tol        -   convergence tolerance                            -
%  max_iter   -   maximum number of Newton iterations              -
%
% OUTPUT:
%
%  z_vec      -   corrected solution vector [X0; T]                -
%  traj       -   trajectory struct from final integration         -
%  hist       -   iteration history struct                         -
%==========================================================================
    
    % DEFAULTS
    if nargin < 3 || isempty(tol)
        tol = 1e-10;
    end
    if nargin < 4 || isempty(max_iter)
        max_iter = 20;
    end

    % INITIALIZATION
    z_vec = z_0_vec(:);

    hist.z_vec = zeros(5, max_iter + 1);
    hist.h_vec = zeros(5, max_iter);
    hist.h_norm_vec = zeros(max_iter, 1);
    hist.Delta_z_vec = zeros(5, max_iter);
    hist.iter = 0;
    hist.converged = false;
    hist.z_vec(:,1) = z_vec;

    % NEWTON LOOP
    for k = 1:max_iter

        % (i) Extract variables
        X0_vec = z_vec(1:4);
        T = z_vec(5);

        % (ii) Integrate trajectory + STM
        traj = integrate_ccr4bp([0, T], X0_vec, params);

        Xf_vec = traj.x_f_vec;
        Phi_f_mtx = traj.Phi_f_mtx;

        % (iii) Residual
        f_vec = Xf_vec - X0_vec;
        g_val = X0_vec(2);

        h_vec = [f_vec; g_val];

        % (iv) Store history
        hist.h_vec(:,k) = h_vec;
        hist.h_norm_vec(k) = norm(h_vec,2);
        hist.iter = k;

        % (v) Convergence check
        if norm(h_vec,2) < tol
            hist.converged = true;
            hist.z_vec = hist.z_vec(:,1:k);
            hist.h_vec = hist.h_vec(:,1:k);
            hist.h_norm_vec = hist.h_norm_vec(1:k);
            hist.Delta_z_vec = hist.Delta_z_vec(:,1:max(k-1,1));
            return
        end

        % (vi) Jacobian
        Df_DX0 = Phi_f_mtx - eye(4);
        Xdot_f_vec = eom_ccr4bp_H(T, Xf_vec, params);

        Dg_Dz = [0, 1, 0, 0, 0];

        Dh_mtx = [Df_DX0, Xdot_f_vec;
                  Dg_Dz];

        % (vii) Newton step
        if rcond(Dh_mtx) < 1e-12
            Delta_z = -pinv(Dh_mtx) * h_vec;
        else
            Delta_z = -(Dh_mtx \ h_vec);
        end

        % (viii) Update
        z_vec = z_vec + Delta_z;

        hist.Delta_z_vec(:,k) = Delta_z;
        hist.z_vec(:,k+1) = z_vec;
    end

    % NO CONVERGENCE
    hist.z_vec = hist.z_vec(:,1:max_iter+1);
    hist.h_vec = hist.h_vec(:,1:max_iter);
    hist.h_norm_vec = hist.h_norm_vec(1:max_iter);
    hist.Delta_z_vec = hist.Delta_z_vec(:,1:max_iter);

    warning('single_shooting_ccr4bp did not converge in %d iterations.', max_iter)

end