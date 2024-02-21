% Copyright: Nestor Cardozo 2018. All rights reserved

%program runModelLRG, two listric thrusts, one from the left and another from the
%right, and growth strata

% Pseudo 3D: Listric fault: Linear variation of trishear parameters from
% southern to northern fault tip and four growth strata, lowermost growth bed 1
% experiences 80%, growth bed 2 60%, growth bed 3 40%, and uppermost growth bed 4
% 20% of incoming fault slip

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

% Datums (Z increases upwards), 4 uppermost beds are growth
prompt = {'Datums, Z increases upwards and 4 uppermost beds are growth'};
def = {'25.0 50.0 60.0 70.0 80.0 90.0'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
datums = str2num(answer{1});

% Make beds (surfaces)
[XP,YP,ZP]=meshgrid(origx:cellsize:extentx+origx,origy:cellsize:extenty+origy,datums);

% Listric Faults
% Left fault
prompt = {'Center of curvature in X (East)','Center of curvature in Z (Up)','Radius of curvature','Maximum central angle'};
def = {'400.0','300.0','300.0','50'};
title = 'LEFT LISTRIC FAULT';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ccxl = str2double(answer{1});
cczl = str2double(answer{2});
ccrl = str2double(answer{3});
maxarcl = str2double(answer{4})*(pi/180);
% Right fault
prompt = {'Center of curvature in X (East)','Center of curvature in Z (Up)','Radius of curvature','Maximum central angle'};
def = {'600.0','300.0','300.0','50'};
title = 'RIGHT LISTRIC FAULT';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ccxr = extentx - str2double(answer{1});
cczr = str2double(answer{2});
ccrr = str2double(answer{3});
maxarcr = str2double(answer{4})*(pi/180);

% Angle of shear in backlimb:
% Left fault
prompt = {'Shear angle in backlimb. Positive angle is for antithetic shear'};
def = {'30.0'};
title = 'Angle of shear left fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ashearl = str2double(answer{1})*(pi/180);
% Right fault
prompt = {'Shear angle in backlimb. Positive angle is for antithetic shear'};
def = {'30.0'};
title = 'Angle of shear right fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
ashearr = str2double(answer{1})*(pi/180);

% Propagation to Slip ratio (P/S)
% Left fault
prompt = {'P/S at southern tip','P/S at northern tip'};
def = {'2.5','2.5'};
title = 'PROPAGATION TO SLIP RATIO left fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
pssl = str2double(answer{1});
psnl = str2double(answer{2});
% Right fault
prompt = {'P/S at southern tip','P/S at northern tip'};
def = {'2.5','2.5'};
title = 'PROPAGATION TO SLIP RATIO right fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
pssr = str2double(answer{1});
psnr = str2double(answer{2});

% Trishear angle (angle of triangular zone of deformation)
% Left fault
prompt = {'Trishear angle at southern tip','Trishear angle at northern tip'};
def = {'60.0','60.0'};
title = 'TRISHEAR ANGLE left fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
trasl = str2double(answer{1})*pi/180;
tranl = str2double(answer{2})*pi/180;
% Right fault
prompt = {'Trishear angle at southern tip','Trishear angle at northern tip'};
def = {'60.0','60.0'};
title = 'TRISHEAR ANGLE right fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
trasr = str2double(answer{1})*pi/180;
tranr = str2double(answer{2})*pi/180;

% Fault Slip
% Slip should be positive for reverse faults and negative for normal faults
% Left fault
prompt = {'Fault slip at southern tip','Fault slip at northern tip','Slip increment'};
def = {'100.0','0.0','1'};
title = 'FAULT SLIP, + IF REVERSE, - IF NORMAL left fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slsl = str2double(answer{1});
slnl = str2double(answer{2});
sincl = str2double(answer{3});
% Right fault
prompt = {'Fault slip at southern tip','Fault slip at northern tip','Slip increment'};
def = {'0.0','100.0','1'};
title = 'FAULT SLIP, + IF REVERSE, - IF NORMAL right fault';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
slsr = str2double(answer{1});
slnr = str2double(answer{2});
sincr = str2double(answer{3});

% fault points and count fault points
% Left fault
faultpl = zeros(50000,3);
countfpl = 0;
% Right fault
faultpr = zeros(50000,3);
countfpr = 0;

% Run model
% For each section parallel to slip (east)
for hh=1:size(XP,1)
    
    % trishear parameters at section
    yc = YP(hh,1,1); % y (north of section)
    % Left fault
    pscl = pssl + (psnl - pssl)/extenty*yc; % ps at section
    tracl = trasl + (tranl - trasl)/extenty*yc; % trishear angle at section
    mcl = tan(tracl/2); % tangent of half trishear angle
    slcl = slsl + (slnl - slsl)/extenty*yc; % slip at section
    % Right fault
    pscr = pssr + (psnr - pssr)/extenty*yc; % ps at section
    tracr = trasr + (tranr - trasr)/extenty*yc; % trishear angle at section
    mcr = tan(tracr/2); % tangent of half trishear angle
    slcr = slsr + (slnr - slsr)/extenty*yc; % slip at section
    
    % Make fault geometry at section
    % Left fault
    [kll,kal,fdl,Rl,nincsl,sincsl] = makeFault(ccxl,cczl,ccrl,maxarcl,pscl,slcl,ashearl,sincl);
    fsl = size(fdl,1)-1;
    % Right fault
    [klr,kar,fdr,Rr,nincsr,sincsr] = makeFault(ccxr,cczr,ccrr,maxarcr,pscr,slcr,ashearr,sincr);
    fsr = size(fdr,1)-1;
    
    % first 20 fault points
    % Left fault
    fainl = [0.0 cczl-ccrl];
    fincl = (kll(1,1)-fainl(1))/20.0;
    fainsl = fainl(1):fincl:kll(1,1)-fincl;
    for ii=1:size(fainsl,2)
        countfpl = countfpl + 1;
        faultpl(countfpl,:) = [fainsl(ii) yc fainl(2)];
    end
    % Right fault
    fainr = [0.0 cczr-ccrr];
    fincr = (klr(1,1)-fainr(1))/20.0;
    fainsr = fainr(1):fincr:klr(1,1)-fincr;
    for ii=1:size(fainsr,2)
        countfpr = countfpr + 1;
        faultpr(countfpr,:) = [fainsr(ii) yc fainr(2)];
    end
    
    %initialize count of points in hanging wall
    hwid=zeros(size(XP,3),1); 
    
    % Run trishear model
    % initialize cummulative slip for section
    % Left fault
    cumslipl = 0.0;
    % Right fault
    cumslipr = 0.0;
    % for all fault segments
    maxfs = fsl;
    if fsl < fsr
        maxfs = fsr;
    end
    for ii=1:maxfs
        % Left fault
        if ii <= fsl
            rampl = fdl(ii+1);
            xtl = kll(ii,1);
            ztl = kll(ii,2);
            a11l=cos(rampl);
            a12l=cos(pi/2-rampl);
            a21l=cos(pi/2+rampl);
            a22l=a11l;
        end
        % Right fault
        if ii <= fsr
            rampr = fdr(ii+1);
            xtr = klr(ii,1);
            ztr = klr(ii,2);
            a11r=cos(rampr);
            a12r=cos(pi/2-rampr);
            a21r=cos(pi/2+rampr);
            a22r=a11r;
        end
        % For slip increments in fault segment
        maxnincs = 0;
        if ii <= fsl
            maxnincs = nincsl(ii);
        end
        if ii <= fsr
            if maxnincs < nincsr(ii)
                maxnincs = nincsr(ii);
            end
        end
        for i=1:maxnincs
            % uppermost active layer - left fault
            uacll = size(XP,3);
            if cumslipl < 0.8*slcl
                uacll = uacll - 1;
            end
            if cumslipl < 0.6*slcl
                uacll = uacll - 1;
            end
            if cumslipl < 0.4*slcl
                uacll = uacll - 1;
            end
            if cumslipl < 0.2*slcl
                uacll = uacll - 1;
            end
            % uppermost active layer - right fault
            uaclr = size(XP,3);
            if cumslipr < 0.8*slcr
                uaclr = uaclr - 1;
            end
            if cumslipr < 0.6*slcr
                uaclr = uaclr - 1;
            end
            if cumslipr < 0.4*slcr
                uaclr = uaclr - 1;
            end
            if cumslipr < 0.2*slcr
                uaclr = uaclr - 1;
            end
            % Left fault
            if ii <= fsl
                if i <= nincsl(ii)
                    % For all active layers
                    for j=1:uacll
                        % x coordinate starting at left end
                        xpr = XP(hh,:,j); 
                        %Find z of fault at x of bed points
                        fatr = [fainl;kll(1:ii,:)];
                        zfi = interp1(fatr(:,1),fatr(:,2),xpr,'linear','extrap');
                        hwid(j) = 0;
                        % For all points in x
                        for k=1:size(XP,2)
                            %point
                            px = xpr(k);
                            pz = ZP(hh,k,j);
                            %point in trishear coordinate system
                            fx=(px-xtl)*a11l+(pz-ztl)*a12l;
                            fz=(px-xtl)*a21l+(pz-ztl)*a22l;
                            % Backlimb (kinem = 0) or trishear (kinem = 1) kinematics?
                            kinem = -1; %neither parallel nor trishear
                            % If in hanging wall
                            if pz > zfi(k)
                                hwid(j)=hwid(j)+1;
                                xkp = xtl - (pz-ztl)/tan(kal(ii));
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
                                xx=fx-(pscl*i*abs(sincsl(ii))*Rl(ii+1));
                                if xx >= 0.0 && fz >-xx*mcl
                                    kinem = 1;
                                end
                            end
                            % Backlimb kinematics
                            if kinem == 0
                                idk = 1;
                                for kk=1:ii
                                    xkp = kll(kk,1) - (pz-kll(kk,2))/tan(kal(kk));
                                    % If ahead of kink
                                    if px >= xkp
                                        idk = kk+1;
                                    end
                                end
                                % Slip is reduced across each kink
                                xpr(k)=px+sincsl(ii)*Rl(idk)*cos(fdl(idk)); 
                                ZP(hh,k,j)=pz+sincsl(ii)*Rl(idk)*sin(fdl(idk)); 
                            end
                            % Trishear kinematics
                            if kinem == 1
                                % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
                                % FAULT TIP. NOTE: MOVE TIP FORWARD AND THEN DEFORM
                                % Slip is reduced across each kink
                                xx=fx-(pscl*i*abs(sincsl(ii))*Rl(ii+1));
                                zz=fz;
                                % compute velocity
                                % Slip is reduced across each kink
                                [vx,vz]=veltrishear(xx,zz,mcl,sincsl(ii)*Rl(ii+1)); 
                                % UPDATE fx, fz COORDINATES
                                fx=fx+vx;
                                fz=fz+vz;
                                % TRANSFORM BACK TO XP, ZP COORDINATE SYSTEM FOR PLOTTING
                                xpr(k)=(fx*a11l+fz*a21l)+xtl;
                                ZP(hh,k,j)=(fx*a12l+fz*a22l)+ztl;
                            end 
                        end
                        XP(hh,:,j) = xpr; 
                    end
                    % update cumulative slip of section
                    cumslipl = cumslipl + sincsl(ii);
                end
            end
            % Right fault
            if ii <= fsr
                if i <= nincsr(ii)
                    % For all active layers
                    for j=1:uaclr
                        % x coordinate starting at right end
                        xpr = extentx - XP(hh,:,j); 
                        %Find z of fault at x of bed points
                        fatr = [fainr;klr(1:ii,:)];
                        zfi = interp1(fatr(:,1),fatr(:,2),xpr,'linear','extrap');
                        hwid(j) = 0;
                        % For all points in x
                        for k=1:size(XP,2)
                            %point
                            px = xpr(k);
                            pz = ZP(hh,k,j);
                            %point in trishear coordinate system
                            fx=(px-xtr)*a11r+(pz-ztr)*a12r;
                            fz=(px-xtr)*a21r+(pz-ztr)*a22r;
                            % Backlimb (kinem = 0) or trishear (kinem = 1) kinematics?
                            kinem = -1; %neither parallel nor trishear
                            % If in hanging wall
                            if pz > zfi(k)
                                hwid(j)=hwid(j)+1;
                                xkp = xtr - (pz-ztr)/tan(kar(ii));
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
                                xx=fx-(pscr*i*abs(sincsr(ii))*Rr(ii+1));
                                if xx >= 0.0 && fz >-xx*mcr
                                    kinem = 1;
                                end
                            end
                            % Backlimb kinematics
                            if kinem == 0
                                idk = 1;
                                for kk=1:ii
                                    xkp = klr(kk,1) - (pz-klr(kk,2))/tan(kar(kk));
                                    % If ahead of kink
                                    if px >= xkp
                                        idk = kk+1;
                                    end
                                end
                                % Slip is reduced across each kink
                                xpr(k)=px+sincsr(ii)*Rr(idk)*cos(fdr(idk)); 
                                ZP(hh,k,j)=pz+sincsr(ii)*Rr(idk)*sin(fdr(idk)); 
                            end
                            % Trishear kinematics
                            if kinem == 1
                                % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
                                % FAULT TIP. NOTE: MOVE TIP FORWARD AND THEN DEFORM
                                % Slip is reduced across each kink
                                xx=fx-(pscr*i*abs(sincsr(ii))*Rr(ii+1));
                                zz=fz;
                                % compute velocity
                                % Slip is reduced across each kink
                                [vx,vz]=veltrishear(xx,zz,mcr,sincsr(ii)*Rr(ii+1)); 
                                % UPDATE fx, fz COORDINATES
                                fx=fx+vx;
                                fz=fz+vz;
                                % TRANSFORM BACK TO XP, ZP COORDINATE SYSTEM FOR PLOTTING
                                xpr(k)=(fx*a11r+fz*a21r)+xtr;
                                ZP(hh,k,j)=(fx*a12r+fz*a22r)+ztr;
                            end 
                        end
                        XP(hh,:,j) = extentx - xpr; 
                    end
                    % update cumulative slip of section - right fault
                    cumslipr = cumslipr + sincsr(ii);
                end
            end
        end
        % intermediate fault points
        % Left fault
        if ii <= fsl
            countfpl = countfpl+1;
            faultpl(countfpl,:) = [kll(ii,1) yc kll(ii,2)];
        end
        % Right fault
        if ii <= fsr
            countfpr = countfpr+1;
            faultpr(countfpr,:) = [klr(ii,1) yc klr(ii,2)];
        end
    end
    % last fault point
    % Left fault
    if fsl > 0
        xtf=xtl+(pscl*nincsl(fsl)*abs(sincsl(fsl)))*cos(rampl);
        ztf=ztl+(pscl*nincsl(fsl)*abs(sincsl(fsl)))*sin(rampl);
        countfpl = countfpl+1;
        faultpl(countfpl,:) = [xtf yc ztf];
    end
    % Right fault
    if fsr > 0
        xtf=xtr+(pscr*nincsr(fsr)*abs(sincsr(fsr)))*cos(rampr);
        ztf=ztr+(pscr*nincsr(fsr)*abs(sincsr(fsr)))*sin(rampr);
        countfpr = countfpr+1;
        faultpr(countfpr,:) = [xtf yc ztf];
    end
end

% Remove empty fault points
% Left fault
countfpl = countfpl+1;
faultpl(countfpl:1:size(faultpl,1),:)= [];
% Right fault
countfpr = countfpr+1;
faultpr(countfpr:1:size(faultpr,1),:)= [];

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
% Left fault
faultpl(any(isnan(faultpl),2),:)=[]; % remove any NAN values
plot3(faultpl(:,1),faultpl(:,2),faultpl(:,3),'r.');
% Right fault
faultpr(any(isnan(faultpr),2),:)=[]; % remove any NAN values
faultpr(:,1)= extentx - faultpr(:,1);
plot3(faultpr(:,1),faultpr(:,2),faultpr(:,3),'r.');

hold off;
axis equal;
xlabel('East');
ylabel('North');
zlabel('Up');