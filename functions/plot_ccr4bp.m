function fig = plot_ccr4bp(traj, fig, plot_opts)
%==========================================================================
%
% Plots x-y trajectory from traj struct returned by integrate_ccr4bp.
% If a figure handle is provided, the trajectory is plotted in that figure.
% Otherwise, a new figure is created.
%
% Author: G. Montseny
% Date: April 23, 2026
%
%
% INPUT:               Description                                   Units
%
%  traj      -   trajectory struct from integrate_ccr4bp               -
%  fig       -   (optional) figure handle                              -
%  plot_opts -   (optional) plotting options struct with fields:       -
%
%                .LineStyle    line style              (default '-')
%                .LineWidth    line width              (default 1.5)
%                .Color        line color              (default [0 0.4470 0.7410])
%                .DisplayName  legend entry            (default '')
%                .Title        plot title              (default '')
%                .XLabel       x-axis label            (default '$x$')
%                .YLabel       y-axis label            (default '$y$')
%                .Interpreter  text interpreter        (default 'latex')
%                .AxisEqual    true/false              (default true)
%                .Grid         true/false              (default true)
%                .HoldOn       true/false              (default true)
%
% OUTPUT:
%
%  fig       -   figure handle                                         -
%==========================================================================
    
    %------------------------------------------------------------------
    % Handle optional inputs
    %------------------------------------------------------------------
    if nargin < 2 || isempty(fig)
        fig = figure();
    else
        figure(fig);
    end

    if nargin < 3 || isempty(plot_opts)
        plot_opts = struct();
    end

    %------------------------------------------------------------------
    % Default plotting options
    %------------------------------------------------------------------
    if ~isfield(plot_opts, 'LineStyle')
        plot_opts.LineStyle = '-';
    end

    if ~isfield(plot_opts, 'LineWidth')
        plot_opts.LineWidth = 1;
    end

    if ~isfield(plot_opts, 'Color')
        plot_opts.Color = [0 0.4470 0.7410];
    end

    if ~isfield(plot_opts, 'DisplayName')
        plot_opts.DisplayName = '';
    end

    if ~isfield(plot_opts, 'Title')
        plot_opts.Title = '';
    end

    if ~isfield(plot_opts, 'XLabel')
        plot_opts.XLabel = '$x$';
    end

    if ~isfield(plot_opts, 'YLabel')
        plot_opts.YLabel = '$y$';
    end

    if ~isfield(plot_opts, 'Interpreter')
        plot_opts.Interpreter = 'latex';
    end

    if ~isfield(plot_opts, 'AxisEqual')
        plot_opts.AxisEqual = true;
    end

    if ~isfield(plot_opts, 'Grid')
        plot_opts.Grid = true;
    end

    if ~isfield(plot_opts, 'HoldOn')
        plot_opts.HoldOn = true;
    end

    %------------------------------------------------------------------
    % Extract x-y data
    %------------------------------------------------------------------
    x_vec = traj.x_vec_hist(:,1);
    y_vec = traj.x_vec_hist(:,2);

    %------------------------------------------------------------------
    % Plot
    %------------------------------------------------------------------
    if plot_opts.HoldOn
        hold on;
    end

    plot(x_vec, y_vec, ...
        'LineStyle', plot_opts.LineStyle, ...
        'LineWidth', plot_opts.LineWidth, ...
        'Color', plot_opts.Color, ...
        'DisplayName', plot_opts.DisplayName);

    xlabel(plot_opts.XLabel, 'Interpreter', plot_opts.Interpreter);
    ylabel(plot_opts.YLabel, 'Interpreter', plot_opts.Interpreter);

    if ~isempty(plot_opts.Title)
        title(plot_opts.Title, 'Interpreter', plot_opts.Interpreter);
    end

    if plot_opts.Grid
        grid on;
    else
        grid off;
    end

    if plot_opts.AxisEqual
        axis equal;
    end

    if ~isempty(plot_opts.DisplayName)
        legend show;
    end

    box on;

end