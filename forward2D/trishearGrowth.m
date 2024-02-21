% Copyright: Nestor Cardozo 2009. All rights reserved

% Program trishearGrowth: 2D forward trishear model
% with growth strata

% Algorithm from Zehnder and Allmendinger (2000)
% Journal of Structural Geology 22, 1009-1014

% INPUT INITIAL GEOMETRY
prompt = {'Datums'};
def = {'50.0 80.0 110.0 140.0 170.0'};
%def = {'150.0 180.0 210.0 240.0 270.0'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yp = str2num(answer{1});
% Top of layers
top = yp(size(yp,2));

prompt = {'Extent of section','Number of points for digitizing beds'};
def = {'1000.0','500'};
title = 'EXTENT, POINTS, AND ORIGIN OF BEDS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extent = str2double(answer{1});
npoint = str2double(answer{2});

% INITIAL COORDINATES OF FAULT TIP
prompt = {'x coordinate','y coordinate'};
def = {'300.0','50.0'};
title = 'FAULT TIP';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
xt = str2double(answer{1});
yt = str2double(answer{2});

% RAMP ANGLE
prompt = {'Ramp Angle'};
def = {'30.0'};
title = 'RAMP ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ramp = str2double(answer{1})*(pi/180);
a11=cos(ramp);
a12=cos(pi/2-ramp);
a21=cos(pi/2+ramp);
a22=a11;

% PROPAGATION TO SLIP RATIO
prompt = {'P/S'};
def = {'1.5'};
title = 'PROPAGATION TO SLIP RATIO';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ps = str2double(answer{1});

% TRISHEAR ANGLE
prompt = {'Trishear Angle'};
def = {'60.0'};
title = 'TRISHEAR ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
tra = str2double(answer{1})*(pi/180);
% Tangent of half trishear angle
m = tan(tra/2);

% SLIP
prompt = {'Slip + if reverse, - if normal','slip increment + if reverse, - if normal'};
def = {'100.0','1.0'};
title = 'SLIP';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slip = str2double(answer{1});
sinc = str2double(answer{2});
ninc = round(slip/sinc);

prompt = {'Regional subsidence/crestal uplift'};
def = {'2.0'};
title = 'Regional subsidence vs crestal uplift';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
subUpl = str2double(answer{1});
% Make ten growth strata
nincG=round(ninc/10);
% Initialize count of growth strata to 1
countG = 1;

% MAKE BED GEOMETRY
xp=0:extent/npoint:extent;
[XP,YP]=meshgrid(xp,yp);

% TRANSFORM TO COORDINATE AXES PARALLEL AND PERPENDICULAR TO THE FAULT, AND
% WITH ORIGIN AT INITIAL FAULT TIP LOCATION
FX=(XP-xt)*a11+(YP-yt)*a12;
FY=(XP-xt)*a21+(YP-yt)*a22;

% RUN TRISHEAR MODEL
for i=1:ninc
    for j=1:size(FX,1)
        for k=1:size(FX,2)
            % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
            % FAULT TIP. NOTE: MOVE TIP FORWARD AND THEN DEFORM
            xx=FX(j,k)-(ps*i*abs(sinc));
            yy=FY(j,k);
            % compute velocity
            [vx,vy]=veltrishear(xx,yy,sinc,m);
            % UPDATE FX, FY COORDINATES
            FX(j,k)=FX(j,k)+vx;
            FY(j,k)=FY(j,k)+vy;
        end
    end
   
   % TRANSFORM BACK TO XP, YP COORDINATE SYSTEM FOR PLOTTING
   XP=(FX*a11+FY*a21)+xt;
   YP=(FX*a12+FY*a22)+yt;
   
   % MAKE FAULT ARRAY
   xtf=xt+(ps*i*abs(sinc))*cos(ramp);
   ytf=yt+(ps*i*abs(sinc))*sin(ramp);
   XF=[xt xtf];
   YF=[yt ytf];
   % MAKE TRISHEAR BOUNDARIES ARRAYS
   axlo=0:10:300;
   htz=axlo*m;
   ftz=-axlo*m;
   XHTZ=(axlo*a11+htz*a21)+xtf;
   YHTZ=(axlo*a12+htz*a22)+ytf;
   XFTZ=(axlo*a11+ftz*a21)+xtf;
   YFTZ=(axlo*a12+ftz*a22)+ytf;
   
   % PLOT INCREMENT
   % Fault
   plot(XF,YF,'r-','LineWidth',2);
   hold on;
   % Hanging wall trishear boundary
   plot(XHTZ,YHTZ,'b-');
   % Footwall trishear boundary
   plot(XFTZ,YFTZ,'b-');
   % Beds
   hw = zeros(1,size(XP,2));
   fw = zeros(1,size(XP,2));
   xhb = zeros(size(XP,1),size(XP,2));
   yhb = zeros(size(XP,1),size(XP,2));
   xfb = zeros(size(XP,1),size(XP,2));
   yfb = zeros(size(XP,1),size(XP,2));
   for j=1:size(XP,1)
      hw(j)=0.0;
      fw(j)=0.0;
      for k=1:size(XP,2)
         if XP(j,k)<=xt+(YP(j,k)-yt)/tan(ramp),
            hw(j)=hw(j)+1;
            xhb(j,hw(j))=XP(j,k);
            yhb(j,hw(j))=YP(j,k);
         else
            fw(j)=fw(j)+1;
            xfb(j,fw(j))=XP(j,k);
            yfb(j,fw(j))=YP(j,k);
         end
      end
      if (j <= size(yp,2))
          % Pregrowth strata
          plot(xhb(j,1:1:hw(j)),yhb(j,1:1:hw(j)),'k-');
          plot(xfb(j,1:1:fw(j)),yfb(j,1:1:fw(j)),'k-');
      else
          % Growth strata
          plot(xhb(j,1:1:hw(j)),yhb(j,1:1:hw(j)),'g-');
          plot(xfb(j,1:1:fw(j)),yfb(j,1:1:fw(j)),'g-');
      end
   end
   axis equal;
   axis([0 extent 0 0.5*extent]);
   hold off;
   frame(i) = getframe;
   
   % ADD GROWTH STRATA
    % CAREFUL: INTERSECTIONS PREGROWTH-GROWTH STRATA ARE NOT
    % CALCULATED. GROWTH STRATA WILL NOT LOOK RIGHT FOR
    % SUBSIDENCE RATE LOWER THAN UPLIFT RATE
    if (i == countG*nincG)
        % Make growth strata
        % Update top
        top = top + nincG*sinc*sin(ramp)*subUpl;
        % Make bed geometry
        xm=i*sinc:extent/npoint:extent+i*sinc;
        [GXP,GYP]=meshgrid(xm,top);
        % TRANSFORM TO COORDINATE AXES PARALLEL AND PERPENDICULAR TO THE FAULT, AND
        % WITH ORIGIN AT INITIAL FAULT TIP LOCATION
        GFX=(GXP-xt)*a11+(GYP-yt)*a12;
        GFY=(GXP-xt)*a21+(GYP-yt)*a22;
        % Add to beds
        XP = [XP; GXP];
        YP = [YP; GYP];
        FX = [FX; GFX];
        FY = [FY; GFY];
        % update count of growth strata
        countG = countG + 1;
    end
   
end