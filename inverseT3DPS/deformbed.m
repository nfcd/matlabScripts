% Copyright: Nestor Cardozo 2009. All rights reserved

function [xpd,ypd,zpd,xtsfd,ytsfd,ztsfd,xtnfd,ytnfd,ztnfd] = deformbed(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake)

% Deform surface using the Trishear 3D algorithm 
% of Cristallini et al. GSA Bulletin v. 116, pp. 938-952 (2004)

% INPUT: beds (xp,yp,zp), fault tip: south (xtsf, ytsf, ztxf), north (xtnf, ytnf, ztnf),
% ramp angle (ramp), P/S south (pss), north (psn),trishear angle south (tras), north(tran),
% slip south (sls), north (sln), number of slip increments (ninc), slip rake (slrake).

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


% The second coordinate system is the FX (perpendicular to fault tip line), 
% FY (perpendicular to fault plane) and FZ (paralllel to fault
% tip line) system

% The FX, FY, FZ system has therefore the following direction cosines (North, east, down)

% For the FX (perpendicular to initial fault tip line) axis
angle1 = striked - 90.0;
if (angle1 < 0.0)
	angle1 = angle1 + 360.0;
end
angle1 = angle1*pi/180.0;
dcosfx1 = cos(angle1)*cos(ramp);
dcosfx2 = sin(angle1)*cos(ramp);
dcosfx3 = -sin(ramp);

% For the FY (perpendicular to fault plane) axis
angle1 = striked + 90.0;
if (angle1 >= 360.0)
	angle1 = angle1 - 360.0;
end
angle1 = angle1*pi/180.0;
dcosfy1=cos(angle1)*cos(pi/2.0-ramp);
dcosfy2=sin(angle1)*cos(pi/2.0-ramp);
dcosfy3=-sin(pi/2.0-ramp);

% For the FZ (parallel to initial fault tip line) axis, use the cross product
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

% BEDS, TRANSFORM XP,YP, ZP TO FX,FY,FZ COORDINATE SYSTEM
% WITH ORIGIN AT THE SOUTHERN FAULT TIP
fx= (xp-xtsf)*gp11+(yp-ytsf)*gp12+(zp-ztsf)*gp13;
fy= (xp-xtsf)*gp21+(yp-ytsf)*gp22+(zp-ztsf)*gp23;
fz= (xp-xtsf)*gp31+(yp-ytsf)*gp32+(zp-ztsf)*gp33;

% transform northern tip to FX, FY, FZ coordinate system
% WITH ORIGIN AT THE SOUTHERN FAULT TIP
northtipfx = (xtnf-xtsf)*gp11+(ytnf-ytsf)*gp12+(ztnf-ztsf)*gp13;
northtipfy = (xtnf-xtsf)*gp21+(ytnf-ytsf)*gp22+(ztnf-ztsf)*gp23;
northtipfz = (xtnf-xtsf)*gp31+(ytnf-ytsf)*gp32+(ztnf-ztsf)*gp33;

% Half trishear angle
htran = tran/2.0;
htras = tras/2.0;
% Slip increment
sincs = sls/ninc;
sincn = sln/ninc;

% Some constants to speed up computation
trueC1 = pss*abs(sincs)*sin(slrake);
trueC2 = pss*abs(sincs)*cos(slrake);
trueC3 = psn*abs(sincn)*sin(slrake);
trueC4 = psn*abs(sincn)*cos(slrake);

% RUN TRISHEAR
for i=1:ninc
    % calculate the angle between the tip line and the north-south (ZP)
    % direction. This is the angle alpha of Cristallini et al. 2004
    % NOTE: MOVE TIP FORWARD AND THEN DEFORM
    newxs = trueC1*i;
	newzs = trueC2*i;
	newxn = northtipfx + trueC3*i;
	newzn = northtipfz + trueC4*i;
    alpha = atan((newxs-newxn)/(newzn-newzs));
    cosAlpha = cos(alpha);
	cos90MinusAlpha = cos(pi/2.0 - alpha);
	cos90PlusAlpha = cos(pi/2.0 + alpha);
    % Length along the fault
    T1 = sqrt((newzn - newzs)*(newzn-newzs)+(newxn - newxs)*(newxn-newxs));
    % variation of half trishear angle
    At = (htran-htras)/T1;
	% variation of slip increment
	Av = (sincn - sincs)/T1;
	% velocities perpendicular and parallel to tip line
	% Equation 4 of Cristallini et al
	voxp = sin(slrake + alpha);
	vozp = cos(slrake + alpha);
	% Equation 14 of Cristallini et al. 2004
	Azp = Av * vozp/2.0;
    for j=1:size(fx,2) 
        % transform coordinates of points to a system attached to the current fault tip line
		% and with origin at current tip 1
		xx = (fx(j) - newxs)*cosAlpha + (fz(j) - newzs)*cos90MinusAlpha;
		zz = (fx(j) - newxs)*cos90PlusAlpha + (fz(j) - newzs)*cosAlpha;
		yy = fy(j);
		% half trishear angle
		htra = At * zz + htras;
        m = tan(htra);
		% slip increment
		sinc = Av * zz + sincs;
        %compute velocity field
        vox = sinc * voxp;
		voz = sinc * vozp;
		[vx,vy,vz]=veltrishear(xx,yy,vox,voz,Azp,At,htra,m);
        % update FX, FY coordinates
        fx(j) = fx(j) + vx*cosAlpha + vz*cos90PlusAlpha;
		fy(j) = fy(j) + vy;
		fz(j) = fz(j) + vx*cos90MinusAlpha + vz*cosAlpha;
    end
end

% Find deformed coordinates
xpd= fx*gp11+fy*gp21+fz*gp31+xtsf;
ypd= fx*gp12+fy*gp22+fz*gp32+ytsf;
zpd= fx*gp13+fy*gp23+fz*gp33+ztsf;

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
% Direction cosines of slip axis-and XP-YP-ZP system
% cosine of the angle between Slip and XP 
slip11=sin(sltrend)*cos(slplunge);
% cosine of the angle between Slip and YP
slip12=cos(sltrend)*cos(slplunge);
% cosine of the angle between Slip and ZP
slip13=sin(slplunge);

% Update coordinates of fault tips

xtsfd = xtsf+pss*abs(sls)*slip11;
xtnfd = xtnf+psn*abs(sln)*slip11;
ytsfd = ytsf+pss*abs(sls)*slip12;
ytnfd = ytnf+psn*abs(sln)*slip12;
ztsfd = ztsf+pss*abs(sls)*slip13;
ztnfd = ztnf+psn*abs(sln)*slip13;