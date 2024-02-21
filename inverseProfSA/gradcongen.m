% Copyright: Nestor Cardozo 2009. All rights reserved

function [history] = gradcongen(xprof,yprof,xo,yo,teta,xs,yu,tetau,ftblock,xtf,ytf,ramp,ps,tra,slip,sinc,maxit)

% 2D
% COMPUTE BEST FIT PARAMETERS: CONSTRAINED OPTIMIZATION USING SIMULATED ANNEALING
% MATLAB FUNCTION simulannealbnd
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE AND
%       NEGATIVE FOR NORMAL FAULTS

% INPUT: Topographic profile (xprof,yprof), location of bed intersections
% (xo,yo), dips at intersections (teta), x location of undeformed
% stratigraphy (xs), tops of undeformed stratigraphy (yu),
% dip of undeformed stratigraphy (tetau), side of fault where
% undeformed stratigraphy is given (ftblock): footwall (0), hangingwall (1)
% current tip (xtf, ytf), ramp angle (ramp), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc), max number of
% iterations(maxit).

% OUTPUT: history structure. Type history.x to get evolution of parameter
% values. Type history.fval to get evolution of objective function (chi
% square values)

% Set up shared variables with outfuns
history.x = [];
history.fval = [];

% -------------------------------------------------------------------------------
% Set parameters to search and construct initial estimate vector x0
% WARNING: All searched parameters must be within the same order of 
% magnitude. For the example in model.jpg, multiply slip, xtf, and ytf by
% 1e-2
% -------------------------------------------------------------------------------

% Uncomment line below if searching for Slip, Trishear angle, PS, Ramp
% angle, and final tip
x0= [slip*1e-2 tra ps ramp xtf*1e-2 ytf*1e-2];

% Uncomment line below if searching for Slip, Trishear angle, PS, and Ramp
% angle
%x0= [slip*1e-2 tra ps ramp];

% Uncomment line below if searching for Slip, Trishear angle, and PS
%x0 = [slip*1e-2 tra ps];

% -------------------------------------------------------------------------------
% Set lower lb and upper ub limits of the parameters
% REMEMBER TO USE SCALED VALUES OF THE PARAMETERS
% -------------------------------------------------------------------------------

% Uncomment two lines below if searching for Slip, Trishear angle, PS, Ramp
% angle, and final tip
lb = [1. 40.*pi/180. 1. 30.*pi/180 3.5 1.5];
ub = [3. 80.*pi/180. 3. 50.*pi/180 4.5 2.5];

% Uncomment two lines below if searching for Slip, Trishear angle, PS, and
% Ramp angle
%lb = [1. 40.*pi/180. 1. 30.*pi/180.];
%ub = [3. 80.*pi/180. 3. 50.*pi/180.];

% Uncomment two lines below if searching for Slip, Trishear angle, and PS
%lb = [1. 40.*pi/180. 1.];
%ub = [3. 80.*pi/180. 3.];

options = saoptimset('OutputFcn',@outfuns,'Display','iter','MaxIter',maxit);
simulannealbnd(@objfun,x0,lb,ub,options);

    function [stop,options,optchanged]  = outfuns(options,optimValues,flag)
        stop = false;
        optchanged = false;
        switch flag
            case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; optimValues.x];
            otherwise
        end
    end

    % objective function
    function f = objfun(x)
        % ------------------------------------------------------------------
        % Unknown parameters:
        % REMEMBER TO RETURN THE PARAMETERS TO THEIR ACTUAL VALUE
        % BEFORE COMPUTING THE OBJECTIVE FUNCTION
        % For the example in model.jpg, multiply slip, xtf, and ytf by 1e2
        % ------------------------------------------------------------------
        slip = x(1)*1e2;
        tra = x(2);   
        ps = x(3);
        % Uncomment line below if searching for ramp angle
        ramp = x(4);
        % Uncomment two lines below if searching for final tip
        xtf = x(5)*1e2;
        ytf = x(6)*1e2;
        % Compute objective function
        f = backtrishear(xprof,yprof,xo,yo,teta,xs,yu,tetau,ftblock,xtf,ytf,ramp,ps,tra,slip,sinc);
    end
end