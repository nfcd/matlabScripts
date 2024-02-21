% Copyright: Nestor Cardozo 2012: Circular fault and inclined shear in the
% backlimb

% Inversion of bed 4 (synthetic model)

load bed4.txt;

% Bed
xp = bed4(:,1)';
yp = bed4(:,2)';
% ----------------------------------------------------------
% Initial guess
% Edit this section to change the initial estimate or guess
% -----------------------------------------------------------
% Center of curvature
ccx = 300.;
ccy = 300.;
% radius of curvature
ccr = 300.;
% maximum central angle of fault
maxarc = 50.0 * pi/180.0;
% P by S
ps = 1.5;
% Trishear angle in degrees
TA = 50.;
tra = TA*pi/180.0;
% Slip and Slip increment
% Positive for reverse fault
% Negative for normal fault
slip = 150.0;
sinc = 1;
% Shear angle in backlimb
ASHEAR = 30.0;
ashear = ASHEAR*pi/180.0;
% Maximum number of iterations in the search
maxit = 1000;

% Inversion
[history] = gradcongen(xp,yp,ccx,ccy,ccr,maxarc,ps,tra,slip,ashear,sinc,maxit);

% Best fit parameters: Use minimum of objective function
% REMEMBER TO RETURN THE PARAMETERS TO THEIR ACTUAL VALUE
% BEFORE COMPUTING THE OBJECTIVE FUNCTION
% For this example, multiply slip, ccx, ccy, ccr by 1e2
[minfval,imin]=min(history.fval);
slipf = history.x(imin,1)*1e2;
traf = history.x(imin,2);
psf = history.x(imin,3);
ccxf = history.x(imin,4)*1e2;
ccyf = history.x(imin,5)*1e2;
ccrf = history.x(imin,6)*1e2;
ashearf = history.x(imin,7);

% restore and find best-fit gaussian line
[xpr,ypr,xpf,ypf] = restorebed(xp,yp,ccxf,ccyf,ccrf,maxarc,psf,traf,slipf,ashearf,sinc);

% deform best-fit line
[xpd,ypd] = deformbed(xpf,ypf,ccxf,ccyf,ccrf,maxarc,psf,traf,slipf,ashearf,sinc);

% Plot
figure('Name','Inversion of Bed 4','NumberTitle','off');
% Data in red
plot(xp,yp,'r.');
hold on;
% Restored data in green
plot(xpr,ypr,'g.');
% Model in blue
plot(xpd,ypd,'b-','LineWidth',2);
% Best fit to restored in gray
plot(xpf,ypf,'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
% Fault
[kl] = makeFault(ccxf,ccyf,ccrf,maxarc,psf,slipf,ashearf,sinc);
fain = [0.0 ccyf-ccrf];
fatr = [fain;kl];
plot(fatr(:,1),fatr(:,2),'r-','LineWidth',2);
%Center of curvature
plot(ccxf,ccyf,'ro');

hold off;
axis equal;