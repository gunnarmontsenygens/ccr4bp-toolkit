function traj = integrate_ccr4bp(tspan, x_0_vec, params, event_fun)
%==========================================================================
%
% Integrates planar CCR4BP Hamiltonian equations of motion together with the STM.
% Augmented state is defined as Y = [x; vec(Phi)].
% All outputs are returned together in the struct traj.
%
% Author: G. Montseny
% Date: April 09, 2026
%
%
% INPUT:               Description                                   Units
%
%  tspan      -   integration interval [t0 tf]                        -
%  x_0_vec    -   initial state (4x1)                                 -
%  params     -   parameter struct                                    -
%  event_fun  -   (optional) event function handle                    -
%
% OUTPUT:
%
%  traj.t_vec        -   time vector                                  -
%  traj.Y_vec_hist   -   augmented state history (N x 20)             -
%                        [x(4); vec(Phi)(16)]
%  traj.x_vec_hist   -   state history (N x 4)                        -
%  traj.Phi_mtx_hist -   STM history (N x 4 x 4)                      -
%  traj.x_f_vec      -   final state (4x1)                            -
%  traj.Phi_f_mtx    -   final STM (4x4)                              -
%  traj.t_e          -   event times                                  -
%  traj.Y_e_vec      -   state at events                              -
%  traj.i_e          -   event indices                                -
%==========================================================================
    
    % Check if there's any events
    if nargin < 4
        event_fun = [];
    end
    
    % Extract base ODE options from params 
    ode_options = params.ode.options;
    
    % Implement event in case it's not empy
    if ~isempty(event_fun)
        ode_options = odeset(ode_options, ...
            'Events', @(t,Y) event_fun(t,Y,params));
    end
    
    % Initialization
    x_0_vec = x_0_vec(:);
    Phi_0_mtx = eye(4);
    Y_0_vec = [x_0_vec; Phi_0_mtx(:)];
    
    % Integration
    if isempty(event_fun)

        [t_vec, Y_vec_hist] = ode113( ...
            @(t,Y) eom_stm_ccr4bp_H(t, Y, params), ...
            tspan, Y_0_vec, ode_options);

        % Empty event outputs
        t_e = [];
        Y_e_vec = [];
        i_e = [];

    else

        [t_vec, Y_vec_hist, t_e, Y_e_vec, i_e] = ode113( ...
            @(t,Y) eom_stm_ccr4bp_H(t, Y, params), ...
            tspan, Y_0_vec, ode_options);

    end

    % Recover state history
    x_vec_hist = Y_vec_hist(:,1:4);

    % Recover STM history as (N x 4 x 4)
    N = size(Y_vec_hist,1);
    Phi_mtx_hist = zeros(N,4,4);
    
    for k = 1:N
        Phi_mtx_hist(k,:,:) = reshape(Y_vec_hist(k,5:20),4,4);
    end

    % Recover final states
    x_f_vec = Y_vec_hist(end,1:4).';
    Phi_f_mtx = reshape(Y_vec_hist(end,5:20),4,4);

    % Store outputs in traj struct
    traj.t_vec = t_vec;
    traj.Y_vec_hist = Y_vec_hist;
    traj.x_vec_hist = x_vec_hist;
    traj.Phi_mtx_hist = Phi_mtx_hist;
    traj.x_f_vec = x_f_vec;
    traj.Phi_f_mtx = Phi_f_mtx;
    traj.t_e = t_e;
    traj.Y_e_vec = Y_e_vec;
    traj.i_e = i_e;

end