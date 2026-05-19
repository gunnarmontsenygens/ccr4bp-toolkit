function D = drm1_d2(dvar, x, y, a, b)
%==========================================================================
%
% Computes second partial derivatives of 1/r potential term.
% Used for constructing the Hessian of the gravitational potential.
%
% Author: G. Montseny
% Date: April XX, 2026
%
%
% INPUT:               Description                                   Units
%
%  dvar      -   derivative type ('xx','yy','xy')                    -
%  x,y       -   evaluation point                                   -
%  a,b       -   source point                                       -
%
% OUTPUT:
%
%  D         -   second derivative value (scalar)                    -
%==========================================================================

    r = sqrt((x-a)^2 + (y-b)^2);

    switch dvar
        case 'xy'
            D = 3*(x-a)*(y-b)/r^5;
        case 'xx'
            D = -1/r^3+3*(x-a)^2/r^5;
        case 'yy'
            D = -1/r^3+3*(y-b)^2/r^5;
        otherwise
            error('Pick correct derivation variables');
    end
end