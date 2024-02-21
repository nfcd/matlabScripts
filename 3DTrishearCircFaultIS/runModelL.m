% Copyright: Nestor Cardozo 2018. All rights reserved

%program runModelL, one listrict thrust from the left

% Pseudo 3D: Listric fault: Linear variation of trishear parameters from
% southern to northern fault tip

%clear variables from workspace
clear variables;

% INITIAL GEOMETRY
% Origin, extent and cell size
prompt = {'Origin in the X (East) direction','Extent of section in the X (East) direction','Origin in the Y (North) direction','Extent of section in the Y (North) direction','Cell size'};
def = {'0.0','1000.0','0.0','1000.0','5.0'};
title = 'X-Y EXTENT AND CELL SIZE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
origx = str2double(answer{1});
extentx = str2double(answer{2});
origy = str2double(answer{3});
extenty = str2double(answer{4});
cellsize = str2double(answer{5});

% Datums (Z increases upwards)
prompt = {'Datums, Z increases upwards'};
def = {'25.0 50.0 75.0 100.0 125.0 150.0'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
datums = str2num(answer{1});

% Make beds (surfaces)
[XP,YP,ZP]=meshgrid(origx:cellsize:extentx+origx,origy:cellsize:extenty+origy,datums);

% Listric Fault
prompt = {'Center of curvature in X (East)','Center of curvature in Z (Up)','Radius of curvature','Maximum central angle'};
def = {'300.0','300.0','300.0','50'};
title = 'LISTRIC FAULT';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ccx = str2double(answer{1});
ccz = str2double(answer{2});
ccr = str2double(answer{3});
maxarc = str2double(answer{4})*(pi/180);

% Angle of shear in backlimb
prompt = {'Shear angle in backlimb. Positive angle is for antithetic shear'};
def = {'30.0'};
title = 'Angle of shear';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ashear = str2double(answer{1})*(pi/180);

% Propagation to Slip ratio (P/S)
prompt = {'P/S at southern tip','P/S at northern tip'};
def = {'2.5','2.5'};
title = 'PROPAGATION TO SLIP RATIO';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
pss = str2double(answer{1});
psn = str2double(answer{2});

% Trishear angle (angle of triangular zone of deformation)
prompt = {'Trishear angle at southern tip','Trishear angle at northern tip'};
def = {'60.0','60.0'};
title = 'TRISHEAR ANGLE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
tras = str2double(answer{1})*pi/180;
tran = str2double(answer{2})*pi/180;

% Fault Slip
% Slip should be positive for reverse faults and negative for normal faults
prompt = {'Fault slip at southern tip','Fault slip at northern tip','Slip increment'};
def = {'100.0','0.0','1'};
title = 'FAULT SLIP, + IF REVERSE, - IF NORMAL';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
sls = str2double(answer{1});
sln = str2double(answer{2});
sinc = str2double(answer{3});

% fault points and count fault points
faultp = zeros(50000,3);
countfp = 0;

% Run model
% For each section parallel to slip (east)
for hh=1:size(XP,1)
    
    % trishear parameters at section
    yc = YP(hh,1,1); % y (north of section)
    psc = pss + (psn - pss)/extenty*yc; % ps at section
    trac = tras + (tran - tras)/extenty*yc; % trishear angle at section
    mc = tan(trac/2); % tangent of half trishear angle
    slc = sls + (sln - sls)/extenty*yc; % slip at section
    
    % Make fault geometry at section
    [kl,ka,fd,R,nincs,sincs] = makeFault(ccx,ccz,ccr,maxarc,psc,slc,ashear,sinc);
    fs = size(fd,1)-1;
    
    % first 20 fault points
    fain = [0.0 ccz-ccr];
    finc = (kl(1,1)-fain(1))/20.0;
    fains = fain(1):finc:kl(1,1)-finc;
    for ii=1:size(fains,2)
        countfp = countfp + 1;
        faultp(countfp,:) = [fains(ii) yc fain(2)];
    end 
    
    %initialize count of points in hanging wall
    hwid=zeros(size(XP,3),1); 
    
    % Run trishear model
    % for all fault segments
    for ii=1:fs
        ramp = fd(ii+1);
        xt = kl(ii,1);
        zt = kl(ii,2);
        a11=cos(ramp);
        a12=cos(pi/2-ramp);
        a21=cos(pi/2+ramp);
        a22=a11;
        % For slip increments in fault segment
        for i=1:nincs(ii)
            % For all layers
            for j=1:size(XP,3)
                % x coordinate starting at left end
                xpr = XP(hh,:,j); 
                %Find z of fault at x of bed points
                fatr = [fain;kl(1:ii,:)];
                zfi = interp1(fatr(:,1),fatr(:,2),xpr,'linear','extrap');
                hwid(j) = 0;
                % For all points in x
                for k=1:size(XP,2)
                    %point
                    px = xpr(k);
                    pz = ZP(hh,k,j);
                    %point in trishear coordinate system
                    fx=(px-xt)*a11+(pz-zt)*a12;
                    fz=(px-xt)*a21+(pz-zt)*a22;
                    % Backlimb (kinem = 0) or trishear (kinem = 1) kinematics?
                    kinem = -1; %neither parallel nor trishear
                    % If in hanging wall
                    if pz > zfi(k)
                        hwid(j)=hwid(j)+1;
                        xkp = xt - (pz-zt)/tan(ka(ii));
                        % If behind most frontal kink = parallel
                        if xpr(k) < xkp
                            kinem = 0;
                        % Else trishear
                        else
                            kinem = 1;
                        end    
                    % If in footwall
                    else
                        % If inside trishear zone
                        xx=fx-(psc*i*abs(sincs(ii))*R(ii+1));
                        if xx >= 0.0 && fz >-xx*mc
                            kinem = 1;
                        end
                    end
                    % Backlimb kinematics
                    if kinem == 0
                        idk = 1;
                        for kk=1:ii
                            xkp = kl(kk,1) - (pz-kl(kk,2))/tan(ka(kk));
                            % If ahead of kink
                            if px >= xkp
                                idk = kk+1;
                            end
                        end
                        % Slip is reduced across each kink
                        xpr(k)=px+sincs(ii)*R(idk)*cos(fd(idk)); 
                        ZP(hh,k,j)=pz+sincs(ii)*R(idk)*sin(fd(idk)); 
                    end
                    % Trishear kinematics
                    if kinem == 1
                        % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
                        % FAULT TIP. NOTE: MOVE TIP FORWARD AND THEN DEFORM
                        % Slip is reduced across each kink
                        xx=fx-(psc*i*abs(sincs(ii))*R(ii+1));
                        zz=fz;
                        % compute velocity
                        % Slip is reduced across each kink
                        [vx,vz]=veltrishear(xx,zz,mc,sincs(ii)*R(ii+1)); 
                        % UPDATE fx, fz COORDINATES
                        fx=fx+vx;
                        fz=fz+vz;
                        % TRANSFORM BACK TO XP, ZP COORDINATE SYSTEM FOR PLOTTING
                        xpr(k)=(fx*a11+fz*a21)+xt;
                        ZP(hh,k,j)=(fx*a12+fz*a22)+zt;
                    end 
                end
                XP(hh,:,j) = xpr; 
            end
        end
        % intermediate fault points
        countfp = countfp+1;
        faultp(countfp,:) = [kl(ii,1) yc kl(ii,2)];
    end
    % last fault point
    if fs > 0
        xtf=xt+(psc*nincs(fs)*abs(sincs(fs)))*cos(ramp);
        ztf=zt+(psc*nincs(fs)*abs(sincs(fs)))*sin(ramp);
        countfp = countfp+1;
        faultp(countfp,:) = [xtf yc ztf];
    end
end

% Remove empty fault points
countfp = countfp+1;
faultp(countfp:1:size(faultp,1),:)= [];

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

% draw fault points
faultp(any(isnan(faultp),2),:)=[]; % remove any NAN values
plot3(faultp(:,1),faultp(:,2),faultp(:,3),'r.');

hold off;
axis equal;
xlabel('East');
ylabel('North');
zlabel('Up');