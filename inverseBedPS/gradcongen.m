% Copyright: Nestor Cardozo 2009. All rights reserved

function [history] = gradcongen(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc,maxit)

% 2D
% COMPUTE BEST FIT PARAMETERS: CONSTRAINED OPTIMIZATION USING PATTERN SEARCH
% MATLAB FUNCTION patternsearch
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE 
%          AND NEGATIVE FOR NORMAL FAULTS

% INPUT: bed (xp,yp), final tip (xtf, ytf), ramp angle (ramp), P/S (ps),
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
% magnitude. For Lin et al. (2007) example, multiply slip, xtf, and ytf by 1e-1
% -------------------------------------------------------------------------------

% Uncomment line below if searching for Slip, Trishear angle, PS, Ramp angle, and final tip
x0= [slip*1e-1 tra ps ramp xtf*1e-1 ytf*1e-1];

% Uncomment line below if searching for Slip, Trishear angle, PS, and Ramp angle
%x0= [slip*1e-1 tra ps ramp];

% Uncomment line below if searching for Slip, Trishear angle, and PS
%x0 = [slip*1e-1 tra ps];

% -------------------------------------------------------------------------------
% Set lower lb and upper ub limits of the parameters
% REMEMBER TO USE SCALED VALUES OF THE PARAMETERS
% -------------------------------------------------------------------------------

% Uncomment two lines below if searching for Slip, Trishear angle, PS, Ramp angle, and final tip
lb = [0.5 20.*pi/180. 1. 25.*pi/180 0.75 0.25];
ub = [1.0 100.*pi/180. 4. 45.*pi/180 1.6 1.];

% Uncomment two lines below if searching for Slip, Trishear angle, PS, and Ramp angle
%lb = [0.5 20.*pi/180. 1. 25.*pi/180.];
%ub = [1.0 100.*pi/180. 4. 45.*pi/180.];

% Uncomment two lines below if searching for Slip, Trishear angle, and PS
%lb = [0.5 20.*pi/180. 1.];
%ub = [1.0 100.*pi/180. 4.];

options = psoptimset('OutputFcn',@outfuns,'Display','iter','MaxIter',maxit);
patternsearch(@objfun,x0,[],[],[],[],lb,ub,[],options);

    function [stop,options,optchanged]  = outfuns(optimValues,options,flag)
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
        % For Lin et al. (2007) example, multiply slip, xtf, and ytf by 1e1
        % ------------------------------------------------------------------
        slip = x(1)*1e1;
        tra = x(2);   
        ps = x(3);
        % Uncomment line below if searching for ramp angle
        ramp = x(4);
        % Uncomment two lines below if searching for final tip
        xtf = x(5)*1e1;
        ytf = x(6)*1e1;
        % Compute objective function
        f = backtrishear(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc);
    end
end