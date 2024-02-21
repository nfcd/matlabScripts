% Copyright: Nestor Cardozo 2012

%Program trishear: 2D forward trishear model with a circular fault,
%parallel shear in the backlimb, and adding growth strata

% Trishear velocity from Zehnder and Allmendinger (2000)
% Journal of Structural Geology 22, 1009-1014

% Input initial geometry
prompt = {'Datums','Beds origin in x'};
def = {'0.0 30.0 60.0 90.0 120.0 150.0 180.0','0.0'};
title = 'BEDS DATUMS AND ORIGIN';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yp = str2num(answer{1});
borig = str2double(answer{2}); 

prompt = {'Extent of section','Number of points for digitizing beds'};
def = {'1000.0','500'};
title = 'EXTENT, POINTS, AND ORIGIN OF BEDS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extent = str2double(answer{1});
npoint = str2double(answer{2});

% Fault
prompt = {'Center of curvature x','Center of curvature y','Radius of curvature','Maximum central angle'};
def = {'300.0','300.0','300.0','50'};
title = 'CENTER AND RADIUS OF CURVATURE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ccx = str2double(answer{1});
ccy = str2double(answer{2});
ccr = str2double(answer{3});
maxarc = str2double(answer{4})*pi/180.;

% Propagation to slip ratio
prompt = {'P/S'};
def = {'1.5'};
title = 'PROPAGATION TO SLIP RATIO';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ps = str2double(answer{1});

% Trishear angle
prompt = {'Trishear Angle'};
def = {'60.0'};
title = 'TRISHEAR ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
tra = str2double(answer{1})*(pi/180);
% Tangent of half trishear angle
m = tan(tra/2);

% Slip
prompt = {'Slip + if reverse, - if normal','slip increment + if reverse, - if normal'};
def = {'200.0','1.0'};
title = 'SLIP';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slip = str2double(answer{1});
sinc = str2double(answer{2});

% Growth strata
prompt = {'Regional subsidence/crestal uplift'};
def = {'2.0'};
title = 'Regional subsidence vs crestal uplift';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
subUpl = str2double(answer{1});

% Make bed geometry
xp=0:extent/npoint:extent;
xp = xp + borig;
[XP,YP]=meshgrid(xp,yp);
% Top of layers
top = max(YP(size(YP,1),:));

% Make fault geometry
[kl,ka,fd,nincs,sincs] = makeFault(ccx,ccy,ccr,maxarc,ps,slip,sinc);
fs = size(fd,1)-1;

% total number of slip increments
ninc = 0; 
for i=1:fs
    ninc = ninc + nincs(i);
end

% Make nine growth strata
nincG=round(ninc/10);
% Initialize count of growth strata to 1
countG = 1;

%Initial point of fault trajectory
fain = [0.0 ccy-ccr];
%initialize count of points in hanging wall
hwid=zeros(size(yp,2),1); 

% Run trishear model
ninc = 0; % number of increments
cslip = 0.0; % cummulative slip 
% for all fault segments
for ii=1:fs
    ramp = fd(ii+1);
    xt = kl(ii,1);
    yt = kl(ii,2);
    a11=cos(ramp);
    a12=cos(pi/2-ramp);
    a21=cos(pi/2+ramp);
    a22=a11;
    % For slip increments in fault segment
    for i=1:nincs(ii)
        for j=1:size(XP,1)
            %Find y of fault at x of bed points
            fatr = [fain;kl(1:ii,:)];
            yfi = interp1(fatr(:,1),fatr(:,2),XP(j,:),'linear','extrap');
            hwid(j) = 0;
            for k=1:size(XP,2)
                %point
                px = XP(j,k);
                py = YP(j,k);
                %point in trishear coordinate system
                fx=(px-xt)*a11+(py-yt)*a12;
                fy=(px-xt)*a21+(py-yt)*a22;
                % Backlimb (kinem = 0) or trishear (kinem = 1) kinematics?
                kinem = -1; %neither parallel nor trishear
                % If in hanging wall
                if py > yfi(k)
                    hwid(j)=hwid(j)+1;
                    xkp = xt - (py-yt)/tan(ka(ii));
                    % If behind most frontal kink = parallel
                    if XP(j,k) < xkp
                        kinem = 0;
                    % Else trishear
                    else
                        kinem = 1;
                    end    
                % If in footwall
                else
                    % If inside trishear zone
                    xx=fx-(ps*i*abs(sincs(ii)));
                    if xx >= 0.0 && fy >-xx*m
                        kinem = 1;
                    end
                end
                % Backlimb parallel kinematics
                if kinem == 0
                    idk = 1;
                    for kk=1:ii
                        xkp = kl(kk,1) - (py-kl(kk,2))/tan(ka(kk));
                        % If ahead of kink
                        if px >= xkp
                            idk = kk+1;
                        end
                    end
                    XP(j,k)=px+sincs(ii)*cos(fd(idk));
                    YP(j,k)=py+sincs(ii)*sin(fd(idk));
                end
                % Trishear kinematics
                if kinem == 1
                    % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
                    % FAULT TIP. NOTE: MOVE TIP FORWARD AND THEN DEFORM
                    xx=fx-(ps*i*abs(sincs(ii)));
                    yy=fy;
                    % compute velocity
                    [vx,vy]=veltrishear(xx,yy,m,sincs(ii));
                    % UPDATE fx, fy COORDINATES
                    fx=fx+vx;
                    fy=fy+vy;
                    % TRANSFORM BACK TO XP, YP COORDINATE SYSTEM FOR PLOTTING
                    XP(j,k)=(fx*a11+fy*a21)+xt;
                    YP(j,k)=(fx*a12+fy*a22)+yt;
                end 
            end
        end
        % MAKE FAULT ARRAY
        xtf=xt+(ps*i*abs(sincs(ii)))*cos(ramp);
        ytf=yt+(ps*i*abs(sincs(ii)))*sin(ramp);
        XF=[fain(1) kl(1:1:ii,1)' xtf];
        YF=[fain(2) kl(1:1:ii,2)' ytf];
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
        % Kinks
        for j=1:ii
            ykink=kl(j,2):5.0:kl(j,2)+300.0;
            xkink=kl(j,1)- (ykink-kl(j,2))/tan(ka(j));
            plot(xkink,ykink,'Color',[1.0 0.5 0.0]);
        end
        % Hanging wall trishear boundary
        plot(XHTZ,YHTZ,'b-');
        % Footwall trishear boundary
        plot(XFTZ,YFTZ,'b-');
        % Center of curvature
        plot(ccx,ccy,'ko');
        % Beds
        for j=1:size(XP,1)
            % Pregrowth strata
            if j <= size(yp,2)
                plot(XP(j,1:1:hwid(j)),YP(j,1:1:hwid(j)),'k-');
                plot(XP(j,hwid(j)+1:1:size(XP,2)),YP(j,hwid(j)+1:1:size(XP,2)),'k-');
            %Growth strata
            else
                plot(XP(j,1:1:hwid(j)),YP(j,1:1:hwid(j)),'g-');
                plot(XP(j,hwid(j)+1:1:size(XP,2)),YP(j,hwid(j)+1:1:size(XP,2)),'g-');
            end
        end
        axis equal;
        axis([0 extent -0.1*extent 0.5*extent]);
        hold off;
        frame(i) = getframe;
        %update number of increments and cummulative slip
        ninc = ninc + 1;
        cslip = cslip + sincs(ii);
        % ADD GROWTH STRATA
        % CAREFUL: INTERSECTIONS PREGROWTH-GROWTH STRATA ARE NOT
        % CALCULATED. GROWTH STRATA WILL NOT LOOK RIGHT FOR
        % SUBSIDENCE RATE LOWER THAN UPLIFT RATE
        if (ninc == countG*nincG)
            % Make growth strata
            % Update top: Notice that uplift rate increases and therefore
            % vertical spacing between growth beds increases
            top = top + (max(YP(size(YP,1),:))-top)*subUpl;
            % Make bed geometry
            xm=cslip:extent/npoint:extent+cslip;
            [GXP,GYP]=meshgrid(xm,top);
            % Add to beds
            XP = [XP; GXP];
            YP = [YP; GYP];
            % update count of growth strata
            countG = countG + 1;
        end   
    end
end