% Copyright: Nestor Cardozo 2012: Circular fault model and parallel shear
% in the backlimb

function [history] = gradcongen(xp,yp,ccx,ccy,ccr,maxarc,ps,tra,slip,sinc,maxit)

% 2D
% COMPUTE BEST FIT PARAMETERS: CONSTRAINED OPTIMIZATION USING SIMULATED ANNEALING
% MATLAB FUNCTION simulannealbnd
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE 
%          AND NEGATIVE FOR NORMAL FAULTS

% INPUT: bed (xp,yp), x and y of center of curvature (ccx, ccy), radius of
% curvature (ccr), maximum central angle of fault (maxarc), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc), 
% and max number of iterations(maxit).

% OUTPUT: history structure. Type history.x to get evolution of parameter
% values. Type history.fval to get evolution of objective function (chi
% square values)

% Set up shared variables with outfuns
history.x = [];
history.fval = [];

% -------------------------------------------------------------------------------
% Set parameters to search and construct initial estimate vector x0
% WARNING: All searched parameters must be within the same order of 
% magnitude. 
% -------------------------------------------------------------------------------

% Uncomment line below if searching for Slip, Trishear angle, PS, center 
% and radius of curvature
x0= [slip*1e-2 tra ps ccx*1e-2 ccy*1e-2 ccr*1e-2];

% -------------------------------------------------------------------------------
% Set lower lb and upper ub limits of the parameters
% REMEMBER TO USE SCALED VALUES OF THE PARAMETERS
% -------------------------------------------------------------------------------

% Uncomment two lines below if searching for Slip, Trishear angle, PS
% center and radius of curvature
lb = [1.0 40.*pi/180. 1. 2.9 2.9 2.9];
ub = [3.0 80.*pi/180. 3. 3.1 3.1 3.1];

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
        % ------------------------------------------------------------------
        slip = x(1)*1e2;
        tra = x(2);   
        ps = x(3);
        % Uncomment 3 lines below if searching for center and radius of curvature
        ccx = x(4)*1e2;
        ccy = x(5)*1e2;
        ccr = x(6)*1e2;
        % Compute objective function
        f = backtrishear(xp,yp,ccx,ccy,ccr,maxarc,ps,tra,slip,sinc);
    end
end