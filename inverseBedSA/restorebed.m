% Copyright: Nestor Cardozo 2009. All rights reserved

function [xpr,ypr,xpf,ypf,xtfr,ytfr] = restorebed(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc) 

% 2D
% Use this function to restore a bed and obtain straight line (xpf,ypf) that best 
% fits restored bed.

% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE AND NEGATIVE FOR NORMAL FAULTS

% INPUT: beds (xp,yp), fault tip (xtf, ytf), ramp angle (ramp), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc).

% OUTPUT: restored bed (xpr,ypr)
% Best fit line to restored bed (xpf,ypf)
% Restored location of fault tip (xtfr,ytfr)

% 2D Retrodeform bed for the griven model parameters

% RAMP ANGLE
a11=cos(ramp);
a12=cos(pi/2.-ramp);
a21=cos(pi/2.+ramp);
a22=a11;

% P/S (restore)
psr = ps * -1.;

% Tangent of half trishear angle
m=tan(tra/2.);

% SLIP INCREMENT (RESTORE): AS A CONTINUOUS PARAMETER
ninc=round(slip/sinc);
sincd = (slip/ninc - sinc);
sinc = sinc + sincd;
slipr = slip * -1.;
sincr = sinc * -1.;

% TRANSFORM TO CURRENT FAULT TIP
fx=(xp-xtf)*a11+(yp-ytf)*a12;
fy=(xp-xtf)*a21+(yp-ytf)*a22;

% RETRODEFORM TRISHEAR MODEL
for i=1:ninc
   	for j=1:size(fx,2)
         % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
         % FAULT TIP. NOTE: FIRST DEFORM AND THEN MOVE TIP BACK
         xx=fx(j)-(psr*(i-1)*abs(sincr));
         yy=fy(j);
         % compute velocity
         [vx,vy]=veltrishear(xx,yy,sincr,m);
         % UPDATE fx, fy coordinates
         fx(j)=fx(j)+vx;
         fy(j)=fy(j)+vy;
    end
end

% TRANSFORM BACK TO PLOTTING COORDINATE SYSTEM
xpr=(fx*a11+fy*a21)+xtf;
ypr=(fx*a12+fy*a22)+ytf;

% Find coefficients b of best fit line
XX = [ones(size(xpr,2),1) xpr'];
YY = ypr';
b = regress(YY,XX);
% Compute best fit line
maxxpr = max(xpr);
minxpr = min(xpr);
intxpr = (maxxpr-minxpr)/500.;
xpf = minxpr:intxpr:maxxpr;
ypf = b(1)+b(2)*xpf;
% Compute restored tip
xtfr=xtf+(psr*abs(slipr))*cos(ramp);
ytfr=ytf+(psr*abs(slipr))*sin(ramp);
