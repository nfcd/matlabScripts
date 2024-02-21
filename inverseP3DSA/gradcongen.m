% Copyright: Nestor Cardozo 2009. All rights reserved

function [history] = gradcongen(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake,maxit)

% COMPUTE BEST FIT PARAMETERS: CONSTRAINED OPTIMIZATION USING SIMULATED ANNEALING
% MATLAB FUNCTION simulannealbnd

% Pseudo 3D algorithm of Cristallini and Allmendinger
% Journal of Structural Geology v. 23, pp. 1883-1899 (2001)

% INPUT: beds (xp,yp,zp), current tip: south (xtsf, ytsf, ztxf), north (xtnf, ytnf, ztnf),
% ramp angle (ramp), P/S south (pss), north (psn),trishear angle south (tras), north(tran),
% slip south (sls), north (sln), number of slip increments (ninc), slip rake (slrake), 
% maximum number of iterations (maxit).

% NOTICE: 1. ENTER ANGLES IN RADIANS
%         2. SLIP SHOULD BE POSTIVE FOR REVERSE FAULTS AND NEGATIVE FOR NORMAL FAULTS

% xp, yp, zp are vectors with size (1, numberOfPoints)

% Set up shared variables with outfuns
history.x = [];
history.fval = [];

% -------------------------------------------------------------------------------
% Set parameters to search and construct initial estimate vector x0
% WARNING: All searched parameters must be within the same order of 
% magnitude. Slip and x, y, and z tip location are scaled by 1e-2
% -------------------------------------------------------------------------------
% uncomment line below if searching for slip, TA, PS, rake, ramp angle, and
% fault tips
%x0= [sls*1e-2 sln*1e-2 tras tran pss psn slrake ramp xtsf*1e-2 ytsf*1e-2 ztsf*1e-2 xtnf*1e-2 ytnf*1e-2 ztnf*1e-2];
% uncomment line below if searching for slip, TA, PS, rake and ramp angle
%x0= [sls*1e-2 sln*1e-2 tras tran pss psn slrake ramp];
% uncomment line below if searching for slip, TA, PS and rake
%x0= [sls*1e-2 sln*1e-2 tras tran pss psn slrake];
% Uncomment line below if searching for Slip TA PS
x0= [sls*1e-2 sln*1e-2 tras tran pss psn];
% -------------------------------------------------------------------------------
% Set lower lb and upper ub limits of the parameters
% REMEMBER TO USE SCALED VALUES OF THE PARAMETERS
% -------------------------------------------------------------------------------
% uncomment two lines below if searching for slip, TA, PS, rake, ramp angle, and
% fault tips
%lb = [0. 0. 20.*pi/180. 20.*pi/180. 1. 1. 60.*pi/180 20.*pi/180 2. 4. 0. 3. -1. 0.];
%ub = [2. 2. 100.*pi/180. 100.*pi/180. 3. 3. 120.*pi/180 40.*pi/180 4. 6. 2. 5. 1. 2.];
% uncomment two lines below if searching for slip, TA, PS, rake and ramp angle
%lb = [0. 0. 20.*pi/180. 20.*pi/180. 1. 1. 60.*pi/180 20.*pi/180];
%ub = [2. 2. 100.*pi/180. 100.*pi/180. 3. 3. 120.*pi/180 40.*pi/180];
% uncomment two lines below if searching for slip, TA, PS and rake
%lb = [0. 0. 20.*pi/180. 20.*pi/180. 1. 1. 60.*pi/180];
%ub = [2. 2. 100.*pi/180. 100.*pi/180. 3. 3. 120.*pi/180];
% uncomment two lines below if searching for slip, TA, PS
lb = [0. 0. 20.*pi/180. 20.*pi/180. 1. 1.];
ub = [2. 2. 100.*pi/180. 100.*pi/180. 3. 3.];

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
        % Nultiply slip, and x, y, and z tip location for 1e2
        % ------------------------------------------------------------------    
        sls = x(1)*1e2;
        sln = x(2)*1e2; 
        tras = x(3);
        tran = x(4);
        pss = x(5);
        psn = x(6);
         % Uncomment lines below if searching for slip rake
        %slrake = x(7);
        % Uncomment line below if searching for ramp angle
        %ramp = x(8);
        % Uncomment lines below if searching for final tips
        %xtsf = x(9)*1e2;
        %ytsf = x(10)*1e2;
        %ztsf = x(11)*1e2;
        %xtnf = x(12)*1e2;
        %ytnf = x(13)*1e2;
        %ztnf = x(14)*1e2;
        
        % objective function
        f = backtrishear(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake);   
    end
end