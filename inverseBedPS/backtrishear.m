% Copyright: Nestor Cardozo 2009. All rights reserved

function chisq = backtrishear(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc) 

% 2D
% Retrodeform bed for the given model parameters
% Compute best fit line
% And return sum of square of residuals (chisq) 
% between plane and best fit line

% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE AND NEGATIVE FOR NORMAL FAULTS

% INPUT: bed (xp,yp), fault tip (xtf, ytf), ramp angle (ramp), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc).

% OUTPUT: chisq which is the objective function

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

% Compute coefficients of best fit line to restored bed
XX = [ones(size(fx,2),1) fx'];
YY = fy';
b = regress(YY,XX);
% Sum of square of residuals (objective function)
chisq = sum((fy-b(1)-b(2)*fx).^2.);