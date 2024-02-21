% Copyright: Nestor Cardozo 2009. All rights reserved

function [xpd,ypd] = deformbed(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc) 

% 2D
% Deform bed for the given model parameters
% Often used to deform best fit line to restored bed
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE + FOR REVERSE AND - FOR NORMAL FAULT

% INPUT: bed (xp,yp), fault tip (xtf, ytf), ramp angle (ramp), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc).

% OUTPUT: deformed bed (xpd,ypd)

% RAMP ANGLE
a11=cos(ramp);
a12=cos(pi/2.-ramp);
a21=cos(pi/2.+ramp);
a22=a11;

% Tangent of half trishear angle
m=tan(tra/2.);

% SLIP INCREMENT (RESTORE): AS A CONTINUOUS PARAMETER
ninc=round(slip/sinc);
sincd = (slip/ninc - sinc);
sinc = sinc + sincd;

% TRANSFORM TO CURRENT FAULT TIP
fx=(xp-xtf)*a11+(yp-ytf)*a12;
fy=(xp-xtf)*a21+(yp-ytf)*a22;

% RUN TRISHEAR MODEL
for i=1:ninc
   	for j=1:size(fx,2)
         % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
         % FAULT TIP. NOTE: FIRST move tip and then deform
         xx=fx(j)-(ps*i*abs(sinc));
         yy=fy(j);
         % compute velocity
         [vx,vy]=veltrishear(xx,yy,sinc,m);
         % UPDATE fx, fy coordinates
         fx(j)=fx(j)+vx;
         fy(j)=fy(j)+vy;
    end
end

% TRANSFORM BACK TO PLOTTING COORDINATE SYSTEM
xpd=(fx*a11+fy*a21)+xtf;
ypd=(fx*a12+fy*a22)+ytf;