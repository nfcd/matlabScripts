% Copyright: Nestor Cardozo 2012: Circular fault model and parallel shear
% in the backlimb

function [xpd,ypd] = deformbed(xp,yp,ccx,ccy,ccr,maxarc,ps,tra,slip,sinc) 

% 2D
% Deform bed for the given model parameters
% Use this function to deform best fit line to restored bed in order to get 
% best fit trishear model. See invBed5.m, invBed6.m or invBed7.m
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE AND NEGATIVE FOR NORMAL FAULTS

% INPUT: bed (xp,yp), x and y of center of curvature (ccx, ccy), radius of
% curvature (ccr), maximum central angle of fault (maxarc), P/S (ps),
% trishear angle (tra), fault slip, slip increment (sinc)

% OUTPUT: deformed bed (xpd,ypd)

% SLIP INCREMENT (RESTORE): AS A CONTINUOUS PARAMETER
ninc=round(slip/sinc);
sincd = (slip/ninc - sinc);
sinc = sinc + sincd;

% Make fault geometry
[kl,ka,fd,nincs,sincs] = makeFault(ccx,ccy,ccr,maxarc,ps,slip,sinc);
fs = size(fd,1)-1;

%Initial point of fault trajectory
fain = [0.0 ccy-ccr];

% Tangent of half trishear angle
m=tan(tra/2.);

% Run trishear model
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
        %Find y of fault at x of bed points
        fatr = [fain;kl(1:ii,:)];
        yfi = interp1(fatr(:,1),fatr(:,2),xp,'linear','extrap');
        % For all points
        for j=1:size(xp,2)
            %point
            px = xp(j);
            py = yp(j);
            %point in trishear coordinate system
            fx=(px-xt)*a11+(py-yt)*a12;
            fy=(px-xt)*a21+(py-yt)*a22;
            % Backlimb (kinem = 0) or trishear (kinem = 1) kinematics?
            kinem = -1; %neither parallel nor trishear
            % If in hanging wall
            if py > yfi(j)
                xkp = xt - (py-yt)/tan(ka(ii));
                % If behind most frontal kink = parallel
                if xp(j) < xkp
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
                xp(j)=px+sincs(ii)*cos(fd(idk));
                yp(j)=py+sincs(ii)*sin(fd(idk));
            end
            % Trishear kinematics
            if kinem == 1
                % SOLVE TRISHEAR IN A COORDINATE SYSTEM ATTACHED TO CURRENT 
                % FAULT TIP. NOTE: FIRST move tip and then deform
                xx=fx-(ps*i*abs(sincs(ii)));
                yy=fy;
                % compute velocity
                [vx,vy]=veltrishear(xx,yy,m,sincs(ii));
                % UPDATE fx, fy COORDINATES
                fx=fx+vx;
                fy=fy+vy;
                % TRANSFORM BACK TO xp, yp COORDINATES
                xp(j)=(fx*a11+fy*a21)+xt;
                yp(j)=(fx*a12+fy*a22)+yt;
            end
        end
    end
end

% DEFORMED COORDINATES
xpd=xp;
ypd=yp;
