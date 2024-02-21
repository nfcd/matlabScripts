% Copyright: Nestor Cardozo 2009. All rights reserved

% Inversion of bed surface

% load bed

load bed.txt;

% x is east-west
% y is north-sout
% z is elevation

xp = bed(:,1)';
yp = bed(:,2)';
zp = bed(:,3)';

% -----------------------
% Known parameters
% -----------------------

% Increments of slip: let's try two hundred
ninc = 200;
% Maximum number of iterations for inversion: Let's try two hundred
maxit = 200;
% Southern fault tip
xtsf = 286.6025;
ytsf = 500.0;
ztsf = 50.0;
% Northern fault tip
xtnf = 373.2051;
ytnf = 0.0;
ztnf = 100.0;
% Ramp angle
ramp = 30.0*pi/180.0;
% Slip rake in degrees
slrake = 90.0*pi/180.0;

% ---------------------
% Inversion: Searching for P/S, Trishear angle, and fault slip at both tips
% ---------------------

% ----------------------------------------------------------
% Initial guess
% Edit this section to change the initial estimate or guess
% -----------------------------------------------------------
% P/S
pss = 1.5;
psn = 1.5;
% Trishear angle
tras = 60.0*pi/180.0;
tran = 60.0*pi/180.0;
% Slip
sls = 75.0;
sln = 75.0;

% Invert
[history] = gradcongen(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake,maxit);

% Results from inversion
sls=history.x(size(history.x,1),1)*1e2;
sln=history.x(size(history.x,1),2)*1e2;
tras=history.x(size(history.x,1),3);
tran=history.x(size(history.x,1),4);
pss = history.x(size(history.x,1),5);
psn = history.x(size(history.x,1),6);

% restore bed
[xpr,ypr,zpr,xpf,ypf,zpf,xtsfr,ytsfr,ztsfr,xtnfr,ytnfr,ztnfr] = restorebed(xp,yp,zp,xtsf,ytsf,ztsf,xtnf,ytnf,ztnf,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake);

% deform best fit plane
[xdef,ydef,zdef,xtsfd,ytsfd,ztsfd,xtnfd,ytnfd,ztnfd] = deformbed(xpf,ypf,zpf,xtsfr,ytsfr,ztsfr,xtnfr,ytnfr,ztnfr,ramp,pss,psn,tras,tran,sls,sln,ninc,slrake);

% plot
figure('Name','Deformed bed and trishear model','NumberTitle','off');
% data
plot3(xp,yp,zp,'k.');
hold on;
% best fit trishear model
plot3(xdef,ydef,zdef,'r.');
hold off;
axis equal;
