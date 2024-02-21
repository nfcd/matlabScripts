% Copyright: Nestor Cardozo 2009. All rights reserved

% Program constantThickness: Simple step, constant thickness 
% fault propagation folding

% Algorithm from Hardy (1997)
% Journal of Structural Geology 19, 893-896

% INPUT INITIAL GEOMETRY

prompt = {'Datums'};
def = {'50.0 100.0 150.0 200.0 250'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yp = str2num(answer{1});
% Base of layers
base = yp(1);

prompt = {'Extent of section','Number of points for digitizing beds'};
def = {'1000.0','500'};
title = 'EXTENT OF SECTION';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extent = str2double(answer{1});
npoint = str2double(answer{2});

% INPUT FAULT GEOMETRY AND SLIP

prompt = {'x coordinate of lowest footwall cutoff','Ramp Angle'};
def = {'400.0','29.0'};
title = 'RAMP ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
xramp = str2double(answer{1});
ramp = str2double(answer{2})*(pi/180);

prompt = {'Slip','Slip increment'};
def = {'100','0.5'};
title = 'SLIP';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slip = str2double(answer{1});
sinc = str2double(answer{2});
ninc=round(slip/sinc);

% SOLVE MAIN PARAMETERS

% See figure 1 of Hardy (1997)
% Solve eq. 1 of Hardy (1997)
options=optimset('display','off');
gamstar = fzero('suppequ',0.5,options,ramp);
% Eq. 2 of Hardy (1997)
gam1 = pi/2. - ramp/2.;
% Eq. 3 of Hardy (1997)
gam = pi/2.+gamstar-gam1;
% Eq. 4 of Hardy (1997)
bet2 = pi - 2.*gamstar;
% other angle for computation
kap = pi - bet2 + ramp;
% Eq. 8 of Hardy (1997)
% Propagation to slip ratio
lbrat = 1./(1.-sin(ramp)/sin(2.*gam-ramp));
% Eqs. 17 and 18 of Hardy (1997)
% Note there are typos in Hardy (1997). Use instead
% Eqs. 13 and 14 of Hardy and Poblet (2005): : Basin Research 17, 417-424.
R1=sin(gam1+ramp)/sin(gam1+gam);
R2=sin(bet2)/sin(bet2-ramp+gam);

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
% Eqs. 9-16 of Hardy (1997)

for i=1:ninc
    % Compute uplift
    lb = lbrat*i*sinc;
    uplift = lb*sin(ramp);
    lbh = lb*cos(ramp);
    % Compute distance ef. Eq. 6 of Hardy (1997)
    ef=uplift/sin(2.*gamstar);
    % Compute fault tip
    xt=xramp+lbh;
    yt=base+uplift;
    % Compute location e in fig. 1 of Hardy (1997)
    xe=xt+ef*cos(kap);
    ye=yt+ef*sin(kap);
    % For each bed
    for j=1:size(XP,1)
        % For every hanging wall point
        for k=1:hwid(i,j)
            if XP(j,k) < xramp - (YP(j,k)-base)/tan(gam1)
                %Region 1, Eqs. 9 an 10 of Hardy (1997)
                XP(j,k) = XP(j,k) + sinc;
            else
                if YP(j,k) < ye
                    % if y lower than y at e
                    if XP(j,k) < xt + (YP(j,k)-yt)/tan(kap)
                        %Region 2, Eqs. 11 and 12 of Hardy (1997)
                        XP(j,k) = XP(j,k) + sinc*cos(ramp);
                        YP(j,k) = YP(j,k) + sinc*sin(ramp);
                    else
                        if XP(j,k) < xt + (YP(j,k)-yt)/tan(gam)
                            %Region 4, Eqs. 15 and 16 of Hardy (1997)
                            XP(j,k) = XP(j,k) + sinc*R2*cos(gam);
                            YP(j,k) = YP(j,k) + sinc*R2*sin(gam);
                        end
                    end
                else
                    % if y higher than y at e
                    if XP(j,k) < xe - (YP(j,k)-ye)/tan(gam1)
                        %Region 2, Eqs. 11 and 12 of Hardy (1997)
                        XP(j,k) = XP(j,k) + sinc*cos(ramp);
                        YP(j,k) = YP(j,k) + sinc*sin(ramp);
                    else
                        if XP(j,k) < xe + (YP(j,k)-ye)/tan(gam)
                            %Region 3, Eqs. 13 and 14 of Hardy (1997)
                            XP(j,k) = XP(j,k) + sinc*R1*cos(gam);
                            YP(j,k) = YP(j,k) + sinc*R1*sin(gam);
                        else
                            if XP(j,k) < xt + (YP(j,k)-yt)/tan(gam)
                                %Region 4, Eqs. 15 and 16 of Hardy (1997)
                                XP(j,k) = XP(j,k) + sinc*R2*cos(gam);
                                YP(j,k) = YP(j,k) + sinc*R2*sin(gam);
                            end
                            
                        end
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
    % Beds
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
    axis equal;
    axis([0 1.5*extent 0 3.*max(yp)]);
    hold off;
    frame(i) = getframe;
end