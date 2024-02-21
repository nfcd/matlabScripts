% Copyright: Nestor Cardozo 2009. All rights reserved

function [xpr,ypr,zpr,xpf,ypf,zpf,xtsfr,ytsfr,ztsfr,xtnfr,ytnfr,ztnfr] = restorebed(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake)

% Retrodeform bed for the griven model parameters
% Pseudo 3D Trishear algorithm of Cristallini and Allmendinger 2001
% (Journal of Structural Geology, 23, 1883-1899)

% INPUT: beds (xp,yp,zp), current tip: south (xtsf, ytsf, ztxf), north (xtnf, ytnf, ztnf),
% ramp angle (ramp), P/S south (pss), north (psn),trishear angle south (tras), north(tran),
% slip south (sls), north (sln), number of slip increments (ninc), slip rake (slrake), 
% maximum number of iterations (maxit).

% NOTICE: 1. ENTER ANGLES IN RADIANS
%         2. SLIP SHOULD BE POSTIVE FOR REVERSE FAULTS AND NEGATIVE FOR NORMAL FAULTS

% xp, yp, zp are vectors with size (1, numberOfPoints)

% The first coordinate system is the XP (east), YP(north), ZP(up) system,
% which is used for plotting. The direction cosines of this system (North, East, Down) are:
% For XP (east) axis
dcosxp2=1.0; %dcosxp1, and dcosxp3 are zero
% For YP (north) axis
dcosyp1=1.0; %dcosyp2, and dcosyp3 are zero
% For ZP (up) axis
dcoszp3=-1.0; %dcoszp1, and dcoszp2 are zero

% Find strike (in degrees) from tips and ramp angle
htlLength = sqrt((xtnf-xtsf)*(xtnf-xtsf)+(ytnf-ytsf)*(ytnf-ytsf));
if (abs(ytnf-ytsf) > 0.0) 
	tlAngle = atan((xtnf-xtsf)/(ytnf-ytsf))*180/pi;
	if (ytsf > ytnf)
		tlAngle = tlAngle + 180.0;
	end
	if (tlAngle >= 360.0)
		tlAngle = tlAngle - 360.0;
	end	
else
    if (xtsf < xtnf)
		tlAngle = 90.0;
	else
		tlAngle = 270.0;
    end
end
if ((ztsf == ztnf) || ramp == pi/2.0) 
	striked = tlAngle;
else 
	dacStrike = abs(ztsf - ztnf)/(tan(ramp));
	diffAngle = asin(dacStrike/htlLength)*180./pi;
	if (ztsf < ztnf)
		striked = tlAngle + diffAngle;
	end
	if (ztsf > ztnf)
		striked = tlAngle - diffAngle;
	end
end
if (striked >= 360.0)
	striked = striked - 360.0;
end
if (striked < 0.0)
	striked = striked + 360.0;
end

% The slip vector trend and plunge are
if (slrake <= pi/2.)
	slrakef = slrake;
else
	slrakef = pi-slrake;
end
angle1 = atan(tan(slrakef)*cos(ramp))*180.0/pi;
if (slrake <= pi/2.)
	oppstriked = striked + 180.0;
	if (oppstriked >= 360.0)
		oppstriked = oppstriked - 360.0;
	end
	sltrend = oppstriked + angle1;
	if (sltrend >= 360.0)
		sltrend = sltrend - 360.0;
	end
else
	sltrend = striked - angle1;
	if (sltrend < 0.0)
		sltrend = sltrend + 360.0;
	end
end
sltrend = sltrend*pi/180.0;
slplunge=asin(sin(slrakef)*sin(ramp));


% The second coordinate system is the FX (parallel to fault slip), 
% FY (perpendicular to fault plane) and FZ (perpendicular to fault
% slip) system. The FX, FY, FZ system has the following direction 
% cosines (North, east, down):

% For the FX (parallel to fault slip) axis
dcosfx1=cos(sltrend)*cos(slplunge);
dcosfx2=sin(sltrend)*cos(slplunge);
dcosfx3=-sin(slplunge);

% For the FY (perpendicular to fault plane) axis
angle1 = striked + 90.0;
if (angle1 >= 360.0)
	angle1 = angle1 - 360.0;
end
angle1 = angle1*pi/180.0;
dcosfy1=cos(angle1)*cos(pi/2.0-ramp);
dcosfy2=sin(angle1)*cos(pi/2.0-ramp);
dcosfy3=-sin(pi/2.-ramp);

% For the FZ (perpendicular to fault slip) axis, use the cross product
% between FX and FY
dcosfz1=dcosfx2*dcosfy3-dcosfx3*dcosfy2;
dcosfz2=dcosfx3*dcosfy1-dcosfx1*dcosfy3;
dcosfz3=dcosfx1*dcosfy2-dcosfx2*dcosfy1;

% The transformation matrix between the XP-YP-ZP system and the FX-FY-FZ system
% is therefore:
% cosine of the angle between FX and XP 
gp11=dcosfx2*dcosxp2;
% cosine of the angle between FX and YP
gp12=dcosfx1*dcosyp1;
% cosine of the angle between FX and ZP
gp13=dcosfx3*dcoszp3;
% cosine of the angle between FY and XP 
gp21=dcosfy2*dcosxp2;
% cosine of the angle between FY and YP
gp22=dcosfy1*dcosyp1;
% cosine of the angle between FY and ZP
gp23=dcosfy3*dcoszp3;
% cosine of the angle between FZ and XP
gp31=dcosfz2*dcosxp2;
% cosine of the angle between FZ and YP
gp32=dcosfz1*dcosyp1;
% cosine of the angle between FZ and ZP
gp33=dcosfz3*dcoszp3;

% BED, TRANSFORM xp-yp-zp TO fx-fy-fz COORDINATE SYSTEM
% WITH ORIGIN AT THE SOUTHERN FAULT TIP
fx= (xp-xtsf)*gp11+(yp-ytsf)*gp12+(zp-ztsf)*gp13;
fy= (xp-xtsf)*gp21+(yp-ytsf)*gp22+(zp-ztsf)*gp23;
fz= (xp-xtsf)*gp31+(yp-ytsf)*gp32+(zp-ztsf)*gp33;

% transform northern tip to FX, FY, FZ coordinate system
% WITH ORIGIN AT THE SOUTHERN FAULT TIP
northtipfx = (xtnf-xtsf)*gp11+(ytnf-ytsf)*gp12+(ztnf-ztsf)*gp13;
northtipfy = (xtnf-xtsf)*gp21+(ytnf-ytsf)*gp22+(ztnf-ztsf)*gp23;
northtipfz = (xtnf-xtsf)*gp31+(ytnf-ytsf)*gp32+(ztnf-ztsf)*gp33;

% Variation of trishear parameters along T1
T1 = northtipfz;
% P/S (restore) and variation
psnr = psn * -1.0;
pssr = pss * -1.0;
Aps = (psnr-pssr)/T1;
% Half trishear angle and variation
htran = tran/2.0;
htras = tras/2.0;
Ata = (htran-htras)/T1;
% Slip increment (restore) and variation
slnr = sln * -1.0;
slsr = sls * -1.0;
sincsr = slsr/ninc;
sincnr = slnr/ninc;
Av = (sincnr-sincsr)/T1;

%RETRODEFORM TRISHEAR MODEL
for i=1:ninc
    for j=1:size(fx,2)
        % P/S (restore)
        psr = Aps*fz(j)+pssr;
        % half trishear angle
        htra = Ata*fz(j)+htras;
        m = tan(htra);
        % slip increment (restore)
        sincr = Av*fz(j)+sincsr;
        % Notice that in the case of a slip vector not
        % perpendicular to the tip line, I have to correct
        % for the distance between XP = 0 and the tip line
        % I do this by substracting (fz(j)/northtipfz)*northtipfx to xx
        % NOTE: FIRST DEFORM AND THEN MOVE TIP BACK
        xx=fx(j)- (fz(j)/northtipfz)*northtipfx - psr*(i-1)*(abs(sincr));
        yy=fy(j);
        [vx,vy]=veltrishear(xx,yy,sincr,m);
        %update fx, fy coordinates
        fx(j)= fx(j) + vx;
        fy(j) = fy(j) + vy;
    end
end

% TRANSFORM FX-FY-FZ TO XP-YP-ZP COORDINATE SYSTEM FOR PLOTTING
xpr= fx*gp11+fy*gp21+fz*gp31+xtsf;
ypr= fx*gp12+fy*gp22+fz*gp32+ytsf;
zpr= fx*gp13+fy*gp23+fz*gp33+ztsf;

% Best fitting plane
XX = [ones(size(xpr,2),1) xpr' ypr'];
YY = zpr';
b = regress(YY,XX);
maxxpr = max(xpr);
minxpr = min(xpr);
maxypr = max(ypr);
minypr = min(ypr);
intxpr = (maxxpr-minxpr)/100.0;
intypr = (maxypr-minypr)/100.0;
xpfl = minxpr:intxpr:maxxpr;
ypfl = minypr:intypr:maxypr;
xpf = zeros(1,size(xpfl,2)*size(ypfl,2));
ypf = zeros(1,size(xpfl,2)*size(ypfl,2));
count = 1;
for i=1:size(ypfl,2)
    for j=1:size(xpfl,2)
        xpf(count) = xpfl(j);
        ypf (count) = ypfl(i);
        count = count + 1;
    end
end
zpf = b(1)+b(2)*xpf+b(3)*ypf;

% fault tips
xtsfr = xtsf+pssr*abs(slsr)*gp11;
xtnfr = xtnf+psnr*abs(slnr)*gp11;
ytsfr = ytsf+pssr*abs(slsr)*gp12;
ytnfr = ytnf+psnr*abs(slnr)*gp12;
ztsfr = ztsf+pssr*abs(slsr)*gp13;
ztnfr = ztnf+psnr*abs(slnr)*gp13;