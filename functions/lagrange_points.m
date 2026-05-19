function L_pts = lagrange_points(mu)
%==========================================================================
%
% Computes L1, L2, L3, L4, L5 for the planar CR3BP
%
% INPUT:
%   mu - mass parameter
%
% OUTPUT:
%   L_pts - struct with fields L1, L2, L3, L4, L5 (each 1x2 vector)
%
%==========================================================================

    options = optimset('TolX',1e-12,'TolFun',1e-12);

    x_m2 = 1 - mu;

    % Brackets for the collinear points
    x_L1 = fzero(@(x) V_x(x,0), [x_m2 - 0.1, x_m2 - 1e-8], options);
    x_L2 = fzero(@(x) V_x(x,0), [x_m2 + 1e-8, x_m2 + 0.1], options);
    x_L3 = fzero(@(x) V_x(x,0), [-1.5, -mu - 1e-8], options);

    % Triangular points
    L4 = [0.5 - mu,  sqrt(3)/2];
    L5 = [0.5 - mu, -sqrt(3)/2];

    % Store output
    L_pts.L1 = [x_L1, 0];
    L_pts.L2 = [x_L2, 0];
    L_pts.L3 = [x_L3, 0];
    L_pts.L4 = L4;
    L_pts.L5 = L5;

    function dVdx = V_x(x,y)
        dVdx = x - (1-mu)*(x+mu)/((x+mu)^2 + y^2)^(3/2)- mu*(x-1+mu)/((x-1+mu)^2 + y^2)^(3/2);
    end
end