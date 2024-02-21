% Copyright: Nestor Cardozo 2009. All rights reserved

% Program fixedAxisGrowth: Simple step, fixed Axis fault propagation folding
% plus growth strata

% Algorithm from Hardy and Poblet (1995)
% Marine and Petroleum Geology 12, 165-176

% INPUT INITIAL GEOMETRY

prompt = {'Datums'};
def = {'50.0 100.0 150.0 200.0 250'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yp = str2num(answer{1});
% Base of layers
base = yp(1);
% Top of layers
top = yp(size(yp,2));

prompt = {'Extent of section','Number of points for digitizing beds'};
def = {'1000.0','500'};
title = 'EXTENT OF SECTION';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extent = str2double(answer{1});
npoint = str2double(answer{2});

% INPUT FAULT GEOMETRY, SLIP, AND SUBSIDENCE RATE

prompt = {'x coordinate of the lowest footwall cutoff','Ramp Angle'};
def = {'400.0','29.0'};
title = 'RAMP ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
xramp = str2double(answer{1});
ramp = str2double(answer{2})*(pi/180);

prompt = {'Slip','Slip increment'};
def = {'100.0','0.5'};
title = 'SLIP';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slip = str2double(answer{1});
sinc = str2double(answer{2});
ninc=round(slip/sinc);

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

% SOLVE MAIN PARAMETERS

% See Figure 1 B of Hardy and Poblet (1995)
% Eq. 7 of Hardy and Poblet (1995)
gam1=(pi-ramp)/2.;
% Eq. 5 of Hardy and Poblet (1995). Assume no excess shear Sp = 0
gamestar = acot((-1.*(2.*cos(ramp)-3.)/sin(ramp))/2.);
% Eq. 6 of Hardy and Poblet (1995)
gamistar=gam1-gamestar;
% Eq. 4 of Hardy and Poblet (1995)
game=acot(cot(gamestar)-2.*cot(gam1));
% Eq. 3 of Hardy and Poblet (1995)
gami = asin((sin(gamistar)*sin(game))/sin(gamestar));
% Eq. 8 of Hardy and Poblet (1995)
% Ratio of backlimb length to total slip
a1=cot(gamestar)-cot(gam1);
a2=1./sin(ramp)-(sin(gami)/sin(game))/sin(game+gami-ramp);
a3=sin(gam1+ramp)/sin(gam1);
lbrat=a1/a2 + a3;
% Change in slip across boundary between backlimb and forelimb regions
% Eq. 23 of Hardy and Poblet (1995) is wrong. Use instead
% Eq. 7 of Hardy and Poblet (2005): Basin Research 17, 417-424.
R=sin(gam1+ramp)/sin(gam1+game);
% Incremental Crestal uplift
% Eq. 15 of Hardy and Poblet (2005)
inccrupl = (sin(gam1)/sin(gam1+game))*sin(game);

% MAKE BED GEOMETRY

xp=0.0:extent/npoint:extent;
[XP,YP]=meshgrid(xp,yp);

% PREPARE PLOTTING
% Identify hangingwall region

hwid = zeros(ninc,size(yp,2));
for i=1:ninc
    uplift = lbrat*i*sinc*sin(ramp);
    for j=1:size(yp,2)
        if yp(j)-base<=uplift
            hwid(i,j)=0;
            for k=1:size(xp,2)
                if xp(k) <= xramp + (yp(j)-base)/tan(ramp)
                    hwid(i,j)=hwid(i,j)+1;
                end
            end
        else
            hwid(i,j)=size(xp,2);
        end
    end
end

% DEFORM BEDS
% Apply velocity field
% Eqs. 17-22 of Hardy and Poblet (1995)

for i=1:ninc
    % Compute uplift
    lb = lbrat*i*sinc;
    uplift = lb*sin(ramp);
    lbh = lb*cos(ramp);
    % Compute point at fault tip
    xt = xramp + lbh;
    yt = base + uplift;
    % For each bed
    for j=1:size(XP,1)
        if j<=size(yp,2)
            points=hwid(i,j);
        else
            points=size(XP,2);
        end
        % For every hanging wall point
        for k=1:points
            if XP(j,k) < xramp - (YP(j,k)-base)/tan(gam1)
                %Region 1, Eqs. 17 and 18 of Hardy and Poblet (1995)
                XP(j,k) = XP(j,k) + sinc;
            else
                if XP(j,k) < xt - (YP(j,k)-yt)/tan(gam1)
                     %Region 2, Eqs. 19 and 20 of Hardy and Poblet (1995)
                     XP(j,k) = XP(j,k) + sinc*cos(ramp);
                     YP(j,k) = YP(j,k) + sinc*sin(ramp);
                else
                    if XP(j,k) < xt + (YP(j,k)-yt)/tan(game)
                        %Region 3, Eqs. 21 and 22 of Hardy and Poblet(1995)
                        XP(j,k) = XP(j,k) + sinc*R*cos(game);
                        YP(j,k) = YP(j,k) + sinc*R*sin(game);
                    end
                end
            end
        end
    end
    % PLOT INCREMENT
    % Fault
    xf=[0 xramp xramp+lbh];
    yf=[base base uplift+base];
    plot(xf,yf,'r-','LineWidth',2);
    hold on;
    % Pregrowth beds
    for j=1:size(yp,2)
        if yp(j)-base <= uplift
            % Beds cut by the fault
            plot(XP(j,1:1:hwid(i,j)),YP(j,1:1:hwid(i,j)),'k-');
            plot(XP(j,hwid(i,j)+1:1:size(xp,2)),YP(j,hwid(i,j)+1:1:size(xp,2)),'k-');
        else
            % Beds not cut by the fault
            plot(XP(j,:),YP(j,:),'k-');
        end
    end
    % Growth beds
    for j=size(yp,2)+1:size(XP,1)
        plot(XP(j,:),YP(j,:),'g-');
    end
    
    axis equal;
    axis([0 1.5*extent 0 3.*max(yp)]);
    hold off;
    frame(i) = getframe;
    
    % ADD GROWTH STRATA
    % CAREFUL: INTERSECTIONS PREGROWTH-GROWTH STRATA ARE NOT
    % CALCULATED. GROWTH STRATA WILL NOT LOOK RIGHT FOR
    % SUBSIDENCE RATE LOWER THAN UPLIFT RATE
    if (i == countG*nincG)
        % Make growth strata
        % Update top
        top = top + nincG*sinc*inccrupl*subUpl;
        % Make bed geometry
        xp=i*sinc:extent/npoint:extent+i*sinc;
        [GXP,GYP]=meshgrid(xp,top);
        % Add to beds
        XP = [XP; GXP];
        YP = [YP; GYP];
        % update count of growth strata
        countG = countG + 1;
    end
end