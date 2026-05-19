function traj = check_symplectic(traj)
%==========================================================================
%
% Computes symplectic error of STM along trajectory using Frobenius norm.
% Also computes determinant of STM and symplectic defect matrix at each step.
% All outputs are stored in the traj struct.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  traj.Y_vec_hist - augmented state history (N x 20)               -
%
% OUTPUT:
%
%  traj.sympl_err_vec - symplecticity error at each time step (N x 1) -
%  traj.det_vec       - determinant of STM at each step (N x 1)       -
%  traj.E_mtx_hist    - symplectic defect matrices (N x 4 x 4)        -
%==========================================================================
    
    % Extract augmented state history from traj
    Y_vec_hist = traj.Y_vec_hist;

    % Number of time steps
    N = size(Y_vec_hist,1);

    % Symplectic matrix
    J = [0, 0, 1, 0;
         0, 0, 0, 1;
        -1, 0, 0, 0;
         0,-1, 0, 0];

    % Initialization
    sympl_err_vec = zeros(N,1);
    det_vec = zeros(N,1);
    E_mtx_hist = zeros(N,4,4);

    % Loop through trajectory
    for k = 1:N

        % Extract STM
        Phi_k_mtx = reshape(Y_vec_hist(k,5:20),4,4);

        % Compute defect
        E_k_mtx = Phi_k_mtx.'*J*Phi_k_mtx - J;

        % Store defect matrix
        E_mtx_hist(k,:,:) = E_k_mtx;

        % Frobenius norm
        sympl_err_vec(k) = norm(E_k_mtx,'fro');

        % Determinant
        det_vec(k) = det(Phi_k_mtx);

    end

    % Store outputs in traj struct
    traj.sympl_err_vec = sympl_err_vec;
    traj.det_vec = det_vec;
    traj.E_mtx_hist = E_mtx_hist;

end