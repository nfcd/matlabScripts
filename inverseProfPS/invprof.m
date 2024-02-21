% Copyright: Nestor Cardozo 2009. All rights reserved

% Read profile
load profile.txt;
xprof = profile(:,1)';
yprof = profile(:,2)';
% Read intersections 
% dips in degrees (positive dips to the left, negative dips to the right)
load intersections.txt;
xo = intersections(:,1)';
yo = intersections(:,2)';
TETA = intersections(:,3)';
teta = TETA*pi/180.0;
% Read undeformed stratigraphy
% Regional dip in degrees
load undstrat.txt;
xs=undstrat(1);
yu=undstrat(2:1:size(undstrat,1)-2)';
TETAU=undstrat(size(undstrat,1)-1);
tetau=TETAU*pi/180.0;
ftblock=undstrat(size(undstrat,1));
% ----------------------------------------------------------
% Initial guess
% Edit this section to change the initial estimate or guess
% -----------------------------------------------------------
% Fault tip
xtf = 400.;
ytf = 200.;
% Ramp angle in degrees
RA = 35.;
ramp = RA*pi/180.0;
% P by S
ps = 2.0;
% Trishear angle in degrees
TA = 70.;
tra = TA*pi/180.0;
% Slip and Slip increment
% Positive for reverse fault
% Negative for normal fault
slip = 250.;
sinc = 1.;
maxit = 500;

% Inversion
[history] = gradcongen(xprof,yprof,xo,yo,teta,xs,yu,tetau,ftblock,xtf,ytf,ramp,ps,tra,slip,sinc,maxit);

% Best fit parameters
slipf = history.x(size(history.x,1),1)*1e2;
traf=history.x(size(history.x,1),2);
TAf = traf*180./pi;
psf = history.x(size(history.x,1),3);
rampf=history.x(size(history.x,1),4);
RAf = rampf*180./pi;
xtff = history.x(size(history.x,1),5)*1e2;
ytff = history.x(size(history.x,1),6)*1e2;

% initial geometry
xtfr = xtff - abs(slipf)*psf*cos(rampf);
ytfr = ytff - abs(slipf)*psf*sin(rampf);
extprof = xprof(size(xprof,2)) - xprof(1);
if ftblock == 1
	xsc = xs - slipf*cos(rampf);
	yuc = yu - slipf*sin(rampf);
else
	xsc = xs;
	yuc = yu;
end
xsl = xprof(1)-extprof*0.5;
yul = yuc - (xsc - xsl)*tan(tetau);
xsr = xprof(size(xprof,2))+extprof*0.5;
xint = (xsr - xsl)/1000.0;
xpf = xsl:xint:xsr;
ypf=zeros(size(xpf));

% plot profile
plot(xprof,yprof);
hold on;
% plot modeled beds
for i=1:size(yul,2)
	for j=1:size(xpf,2)
		ypf(j)=yul(i)+(xpf(j)-xpf(1))*tan(tetau);
	end
	% deform
	[xpd,ypd] = deformbed(xpf,ypf,xtfr,ytfr,rampf,psf,traf,slipf,sinc);
	plot(xpd,ypd,'b-');
end

% plot intersections
tickl = 3.0*xint;
for i=1:size(xo,2)
	plot(xo(i),yo(i),'ro');
	if abs(teta(i)) > 0.0
        xob = xo(i)-tickl*cos(teta(i))*sign(sin(teta(i)));
    else
        xob = xo(i)-tickl*cos(teta(i));
    end
	yob = yo(i)-abs(tickl*sin(teta(i)));
	tickx = [xob xo(i)];
	ticky = [yob yo(i)];
	plot(tickx,ticky,'r-');
end

% plot fault
xfault = [xtfr xtff];
yfault = [ytfr ytff];
plot(xfault,yfault,'ro-','LineWidth',2);
hold off;
axis equal;