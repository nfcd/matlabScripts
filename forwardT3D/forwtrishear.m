% Copyright: Nestor Cardozo 2009. All rights reserved

%program forwtrishear

% Forward trishear modeling using the Trishear 3D algorithm
% of Cristallini et al. GSA Bulletin v. 116, pp. 938-952 (2004)

% The first coordinate system is the XP (east), YP(north), ZP(up) system,
% which is used for plotting. The direction cosines of this system (North, East, Down) are:
% For XP (east) axis
dcosxp2=1.0; %dcosxp1, and dcosxp3 are zero
% For YP (north) axis
dcosyp1=1.0; %dcosyp2, and dcosyp3 are zero
% For ZP (up) axis
dcoszp3=-1.0; %dcoszp1, and dcoszp2 are zero

% The second coordinate system is the FX (perpendicular to fault tip line), 
% FY (perpendicular to fault plane) and FZ (paralllel to fault
% tip line) system

% The FX, FY, FZ system has therefore the following direction cosines (North, east, down)

% For the FX (perpendicular to initial fault tip line) axis
striked = strike*180.0/pi;
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
FX= (XP-xts)*gp11+(YP-yts)*gp12+(ZP-zts)*gp13;
FY= (XP-xts)*gp21+(YP-yts)*gp22+(ZP-zts)*gp23;
FZ= (XP-xts)*gp31+(YP-yts)*gp32+(ZP-zts)*gp33;

% transform northern tip to FX, FY, FZ coordinate system
% WITH ORIGIN AT THE SOUTHERN FAULT TIP
northtipfx = (xtn-xts)*gp11+(ytn-yts)*gp12+(ztn-zts)*gp13;
northtipfy = (xtn-xts)*gp21+(ytn-yts)*gp22+(ztn-zts)*gp23;
northtipfz = (xtn-xts)*gp31+(ytn-yts)*gp32+(ztn-zts)*gp33;

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
    for j=1:size(FX,1) 
        for k=1:size(FX,2)  
            for l=1:size(FX,3) 
           	    % transform coordinates of points to a system attached to the current fault tip line
			    % and with origin at current tip 1
			    xx = (FX(j,k,l) - newxs)*cosAlpha + (FZ(j,k,l) - newzs)*cos90MinusAlpha;
			    zz = (FX(j,k,l) - newxs)*cos90PlusAlpha + (FZ(j,k,l) - newzs)*cosAlpha;
			    yy = FY(j,k,l);
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
                FX(j,k,l) = FX(j,k,l) + vx*cosAlpha + vz*cos90PlusAlpha;
				FY(j,k,l) = FY(j,k,l) + vy;
				FZ(j,k,l) = FZ(j,k,l) + vx*cos90MinusAlpha + vz*cosAlpha;  
            end
        end
    end
end

% TRANSFORM FX-FY-FZ TO XP-YP-ZP COORDINATE SYSTEM FOR PLOTTING
XP= FX*gp11+FY*gp21+FZ*gp31+xts;
YP= FX*gp12+FY*gp22+FZ*gp32+yts;
ZP= FX*gp13+FY*gp23+FZ*gp33+zts;

% ***************************************************
% PLOTTING
% ***************************************************
% 3D Plot
%S = [48 0];
%K = [0.3, 0.0, 0.5, 10000];
%V = 100:10:250;
%V = 200:10:340;
%VV = 200:50:350;

for i=1:size(XP,3) % vary in ZP
    xplot = zeros(size(XP,1),size(XP,2));
    yplot = zeros(size(XP,1),size(XP,2));
    zplot = zeros(size(XP,1),size(XP,2));
    
    for j=1:size(XP,1) % vary in YP
        for k=1:size(XP,2) % vary in XP
            xplot(j,k)=XP(j,k,i);
            yplot(j,k)=YP(j,k,i);
            zplot(j,k)=ZP(j,k,i);
        end
    end
    mesh(xplot,yplot,zplot);
    %surfl(xplot,yplot,zplot,S,K);
    %shading interp
    %colormap(gray);
    %[CS,H] = contour3(xplot,yplot,zplot,V);
    %clabel(CS,H,VV,'fontsize',10);
    hold on;
    
    
    
end

% draw fault
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
xtsf = xts+pss*abs(sls)*slip11;
xtnf = xtn+psn*abs(sln)*slip11;
xfault=[xts,xtsf,xtnf,xtn,xts];
ytsf = yts+pss*abs(sls)*slip12;
ytnf = ytn+psn*abs(sln)*slip12;
yfault=[yts,ytsf,ytnf,ytn,yts];
ztsf = zts+pss*abs(sls)*slip13;
ztnf = ztn+psn*abs(sln)*slip13;
zfault=[zts,ztsf,ztnf,ztn,zts];
plot3(xfault,yfault,zfault,'b');

hold off;
axis equal;
axis off;
%view(33,9);
view(-31,29);
%xlabel('East');
%ylabel('North');
%zlabel('Up');

