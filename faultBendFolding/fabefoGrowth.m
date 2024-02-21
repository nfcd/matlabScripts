% Copyright: Nestor Cardozo 2009. All rights reserved

% Program fabefo

% Simple step, Mode I fault bend folding
% Algorithm by Stuart Hardy (1995)
% Journal of Structural Geology 17, 1785-1788

% INPUT INITIAL GEOMETRY

prompt = {'Datums'};
def = {'50.0 100.0 150.0 200.0 250'};
title = 'DATUMS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
yp = str2num(answer{1});
% Top of layers
top = yp(size(yp,2));

prompt = {'Extent of section','Number of points for digitizing beds'};
def = {'1000.0','500'};
title = 'EXTENT OF SECTION';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
extent = str2double(answer{1});
npoint = str2double(answer{2});

% INPUT FAULT GEOMETRY AND SLIP

prompt = {'x coordinate of lowest footwall cutoff','Ramp Angle','Ramp Height'};
def = {'400','25','100'};
title = 'RAMP GEOMETRY';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
xramp = str2double(answer{1});
ramp = str2double(answer{2})*(pi/180);
height = str2double(answer{3});

prompt = {'Slip','Slip increment'};
def = {'300','1'};
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

% Avoid convergence problems when ramp
% is greater or equal to 30 degrees

if ramp > 30*pi/180
    error('ramp angle cannot be more than 30');
elseif ramp == 30*pi/180
    ramp=29.9*pi/180;
end

% SOLVE MAIN PARAMETERS
% Equations 1 to 3 of Hardy (1995)

options=optimset('display','off');
gama = fzero('suppequ',1.5,options,ramp);
beta = pi - 2*gama;
R = sin(gama - ramp)/sin(gama);

% MAKE BED GEOMETRY

xp=0.0:extent/npoint:extent;
[XP,YP]=meshgrid(xp,yp);

% PREPARE PLOTTING

% Make fault
xf=[0 xramp xramp+height/tan(ramp) 1.5*extent];
yf=[0 0 height height];
    
% Identify hangingwall region
hwid = zeros(size(yp,2));
for i=1:size(yp,2)
   if yp(i) <= height
     hwid(i)=0;
   	 for j=1:size(xp,2)
      	if xp(j) <= xramp + yp(i)/tan(ramp)
         	hwid(i)= hwid(i)+1;
      	end
   	 end
   else
      hwid(i)=size(xp,2);
   end
end

% DEFORM BEDS
% Apply velocity field
% Eqs. 8 to 13 of Hardy (1995)

for i=1:ninc
    % For each bed
    for j=1:size(XP,1)
        if j<=size(yp,2)
            points=hwid(j);
        else
            points=size(XP,2);
        end
        % For every hanging wall point
        for k=1:points
            if XP(j,k) < xramp - YP(j,k)*tan(ramp/2)
                XP(j,k) = XP(j,k) + sinc;
                YP(j,k) = YP(j,k);
            else
                if YP(j,k) < height
                    XP(j,k) = XP(j,k) + sinc*cos(ramp);
                    YP(j,k) = YP(j,k) + sinc*sin(ramp);
                else  
                    if i*sinc*sin(ramp) < height
                        if XP(j,k) < xramp + height/tan(ramp) + (YP(j,k)-height)*tan(pi/2-gama)
            		        XP(j,k) = XP(j,k) + sinc*cos(ramp);
                            YP(j,k) = YP(j,k) + sinc*sin(ramp);	
                        else
                            XP(j,k)= XP(j,k) + sinc*R;
            	            YP(j,k)= YP(j,k);
                        end    
                    else      
                        if XP(j,k) < xramp + height/tan(ramp)-(YP(j,k)-height)*tan(ramp/2)
            	            XP(j,k)= XP(j,k) + sinc*cos(ramp);
            		        YP(j,k)= YP(j,k) + sinc*sin(ramp);
         		        else
            		        XP(j,k) = XP(j,k) + sinc*R;
                            YP(j,k) = YP(j,k);
                        end
                    end
                end
            end
        end
    end
    % PLOT INCREMENT
    % Fault
    plot(xf,yf,'r-','LineWidth',2);
    hold on;
    % Pregrowth beds
    for j=1:size(yp,2)
      if yp(j) <= height
         plot(XP(j,1:1:hwid(j)),YP(j,1:1:hwid(j)),'k-');
         plot(XP(j,hwid(j)+1:1:size(xp,2)),YP(j,hwid(j)+1:1:size(xp,2)),'k-');
      else
         plot(XP(j,:),YP(j,:),'k-');
      end
    end
    % Growth beds
    for j=size(yp,2)+1:size(XP,1)
        plot(XP(j,:),YP(j,:),'g-');
    end
    
    axis equal
    axis([0 1.5*extent 0 3.0*max(yp)]);
    hold off
    frame(i) = getframe;
    % ADD GROWTH STRATA
    % CAREFUL: INTERSECTIONS PREGROWTH-GROWTH STRATA ARE NOT
    % CALCULATED. GROWTH STRATA WILL NOT LOOK RIGHT FOR
    % SUBSIDENCE RATE LOWER THAN UPLIFT RATE
    if (i == countG*nincG)
        % Make growth strata
        % Update top
        totUpl = nincG*sinc*sin(ramp);
        if (totUpl <= height)
            top = top + totUpl*subUpl;
        end
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