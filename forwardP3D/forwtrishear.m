% Copyright: Nestor Cardozo 2009. All rights reserved

%program forwtrishear

% PSEUDO-3D
% Forward trishear modeling using the pseudo 3D Trishear algorithm 
% of Cristallini and Allmendinger 2001 (Journal of Structural Geology, 23, 1883-1899)

% The first coordinate system is the XP (east), YP(north), ZP(up) system,
% which is used for plotting. The direction cosines of this system (North, East, Down) are:
% For XP (east) axis
dcosxp2=1.0; %dcosxp1, and dcosxp3 are zero
% For YP (north) axis
dcosyp1=1.0; %dcosyp2, and dcosyp3 are zero
% For ZP (up) axis
dcoszp3=-1.0; %dcoszp1, and dcoszp2 are zero

% The slip vector trend and plunge are
striked = strike*180.0/pi;
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
dcosfy3=-sin(pi/2.0-ramp);

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

% Variation of trishear parameters along T1
T1 = northtipfz;
% Variation of P/S
Aps = (psn-pss)/T1;
% Half trishear angle and variation
htran = tran/2.0;
htras = tras/2.0;
Ata = (htran-htras)/T1;
% incremental slip and variation
sincs = sls/ninc;
sincn = sln/ninc;
Av = (sincn-sincs)/T1;

% RUN TRISHEAR
for i=1:ninc
    for j=1:size(FX,1) 
        for k=1:size(FX,2)  
            for l=1:size(FX,3) 
                % P/S
                ps = Aps*FZ(j,k,l)+pss;
                % half trishear angle
                htra = Ata*FZ(j,k,l)+htras;
                m = tan(htra);
                % slip increment
                sinc = Av*FZ(j,k,l)+sincs;
                % Notice that in the case of a slip vector not
                % perpendicular to the tip line, I have to correct
                % for the distance between XP = 0 and the tip line
                % I do this by substracting (FZ(j,k,l)/northtipfz)*northtipfx to xx
                % NOTE: MOVE TIP FORWARD AND THEN DEFORM
                xx=FX(j,k,l)- (FZ(j,k,l)/northtipfz)*northtipfx - ps*i*(abs(sinc));
                yy=FY(j,k,l);
                [vx,vy]=veltrishear(xx,yy,sinc,m);
                % update FX, FY coordinates
                FX(j,k,l)= FX(j,k,l) + vx;
                FY(j,k,l) = FY(j,k,l) + vy;  
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
    hold on;
end

% draw fault
xtsf = xts+pss*abs(sls)*gp11;
xtnf = xtn+psn*abs(sln)*gp11;
xfault=[xts,xtsf,xtnf,xtn,xts];
ytsf = yts+pss*abs(sls)*gp12;
ytnf = ytn+psn*abs(sln)*gp12;
yfault=[yts,ytsf,ytnf,ytn,yts];
ztsf = zts+pss*abs(sls)*gp13;
ztnf = ztn+psn*abs(sln)*gp13;
zfault=[zts,ztsf,ztnf,ztn,zts];
plot3(xfault,yfault,zfault,'b');

hold off;
axis equal;
xlabel('East');
ylabel('North');
zlabel('Up');

