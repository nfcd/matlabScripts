% Copyright: Nestor Cardozo 2012: Circular fault model and parallel shear
% in the backlimb

function chisq = backtrishear(xp,yp,ccx,ccy,ccr,maxarc,ps,tra,slip,sinc) 

% Retrodeform bed for the given model parameters
% Compute best fit line (gaussian function)
% And return sum of square of residuals (chisq) 
% between bed and best fit line

% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE POSITIVE FOR REVERSE AND NEGATIVE FOR NORMAL FAULTS

% INPUT: bed (xp,yp), x and y of center of curvature (ccx, ccy), radius of
% curvature (ccr), maximum central angle of fault (maxarc), P/S (ps),trishear angle 
%(tra), slip, slip increment (sinc)

% OUTPUT: chisq which is the objective function

% SLIP INCREMENT (RESTORE): AS A CONTINUOUS PARAMETER
ninc=round(slip/sinc);
sincd = (slip/ninc - sinc);
sinc = sinc + sincd;
sincr = sinc * -1.;

% Make fault geometry
[kl,ka,fd,nincs,sincs] = makeFault(ccx,ccy,ccr,maxarc,ps,slip,sincr);
fs = size(fd,1)-1;

% P/S (restore)
psr = ps * -1.;

% Tangent of half trishear angle
m=tan(tra/2.);

%Initial point of fault trajectory
fain = [0.0 ccy-ccr];

% Retrodeform trishear model
% for all fault segments
for ii=fs:-1:1  %Notice the change here from forward
    ramp = fd(ii+1);
    xt = kl(ii+1,1); %Notice the change here from forward
    yt = kl(ii+1,2); %Notice the change here from forward
    a11=cos(ramp);
    a12=cos(pi/2-ramp);
    a21=cos(pi/2+ramp);
    a22=a11;
    % For slip increments in fault segment
    for i=1:nincs(ii)
        %Find y of fault at x of bed points
        fatr = [fain;kl(1:ii+1,:)]; % Notice the change here from forward
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
                xx=fx-(psr*(i-1)*abs(sincs(ii)));
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
                % FAULT TIP. NOTE: FIRST DEFORM AND THEN MOVE TIP BACK
                xx=fx-(psr*(i-1)*abs(sincs(ii)));
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

% Compute coefficients of best fit line to restored bed
XX = [ones(size(xp,2),1) xp'];
YY = yp';
b = regress(YY,XX);
% Sum of square of residuals (objective function)
yfit = b(1)+b(2)*xp;
chisq = sum((yp-yfit).^2.0);