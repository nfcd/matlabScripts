% Copyright: Nestor Cardozo 2009. All rights reserved

% Inversion of bed 7
% Shijia trench site, Chelungpu fault
% west central Taiwan. Lin et al. 2007
% Journal of Structural Geology 29, 1267-1280

load bed7.txt;

% Bed
xp = bed7(:,1)';
yp = bed7(:,2)';
% ----------------------------------------------------------
% Initial guess
% Edit this section to change the initial estimate or guess
% -----------------------------------------------------------
% Fault tip
xtf = 12.;
ytf = 6.;
% Ramp angle in degrees
RA = 35.;
ramp = RA*pi/180.0;
% P by S
ps = 2.;
% Trishear angle in degrees
TA = 60.;
tra = TA*pi/180.0;
% Slip and Slip increment
% Positive for reverse fault
% Negative for normal fault
slip = 6.0;
sinc = 0.01;
maxit = 200;

% Inversion
[history] = gradcongen(xp,yp,xtf,ytf,ramp,ps,tra,slip,sinc,maxit);

% Best fit parameters
% REMEMBER TO RETURN THE PARAMETERS TO THEIR ACTUAL VALUE
% BEFORE COMPUTING THE OBJECTIVE FUNCTION
% For Lin et al. (2007) example, multiply slip, xtf, and ytf by 1e1
slipf = history.x(size(history.x,1),1)*1e1;
traf = history.x(size(history.x,1),2);
TAf = traf*180./pi;
psf = history.x(size(history.x,1),3);
rampf = history.x(size(history.x,1),4);
RAf = rampf*180./pi;
xtff = history.x(size(history.x,1),5)*1e1;
ytff = history.x(size(history.x,1),6)*1e1;

% restore and find best-fit line
[xpr,ypr,xpf,ypf,xtfr,ytfr] = restorebed(xp,yp,xtff,ytff,rampf,psf,traf,slipf,sinc);

% deform best-fit line
[xpd,ypd,xtfd,ytfd] = deformbed(xpf,ypf,xtfr,ytfr,rampf,psf,traf,slipf,sinc);

% Compare
figure('Name','Inversion of Bed 7','NumberTitle','off');
% Data in red
plot(xp,yp,'r.');
hold on;
% Model in blue
plot(xpd,ypd,'b-','LineWidth',2);
% fault in red
xfault = [xtfr xtfd];
yfault = [ytfr ytfd];
plot(xfault,yfault,'ro-','LineWidth',2);
hold off;
axis equal;