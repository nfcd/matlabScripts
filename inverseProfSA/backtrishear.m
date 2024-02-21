% Copyright: Nestor Cardozo 2009. All rights reserved

function chisq = backtrishear(xprof,yprof,xo,yo,teta,xs,yu,tetau,ftblock,xtf,ytf,ramp,ps,tra,slip,sinc) 

% 2D
% COMPUTE OBJECTIVE FUNCTION: Error between model and data
% NOTE: 1. ENTER ANGLES IN RADIANS
%       2. SLIP AND SLIP INCREMENT SHOULD BE + FOR REVERSE AND - FOR NORMAL FAULT

% INPUT: Topographic profile (xprof,yprof), location of bed intersections
% (xo,yo), dips of intersections (teta), x location of undeformed
% stratigraphy (xs), tops of undeformed stratigraphy (yu),
% dip of undeformed stratigraphy (tetau), side of fault where
% undeformed stratigraphy is given (ftblock): footwall (0), hangingwall (1),
% current tip (xtf, ytf), ramp angle (ramp), P/S (ps),
% trishear angle (tra), slip, slip increment (sinc), max number of
% iterations(maxit).

% OUTPUT: chisq which is the objective function

% initial geometry

% Compute initial location of fault tip
xtfr = xtf - abs(slip)*ps*cos(ramp);
ytfr = ytf - abs(slip)*ps*sin(ramp);
% Construc initial layer template
extprof = xprof(size(xprof,2)) - xprof(1);
% If undeformed stratigraphy is in the hangingwall, correct for fault slip
if ftblock == 1
	xsc = xs - slip*cos(ramp);
	yuc = yu - slip*sin(ramp);
else
	xsc = xs;
	yuc = yu;
end
% Layer template extend outside the topographic profile, 
% half the profile on both sides of the profile
xsl = xprof(1)-extprof*0.5;
yul = yuc - (xsc - xsl)*tan(tetau);
xsr = xprof(size(xprof,2))+extprof*0.5;
xint = (xsr - xsl)/1000.0;
xpf = xsl:xint:xsr;
ypf = zeros(size(xpf));

% Sort intersections in x
INT = [xo' yo' teta'];
SORTINT = sortrows(INT);

% Initialize modeled intersections
MODINT = zeros(100,3);
% Initialize x and y of intersections to 2*xsr and
% dips to 720 degrees: This will create big
% errors if intersections are not found
dumval1 = 2.0*xsr;
dumval2 = 720.0*pi/180.0;
for i=1:100
    MODINT(i,1)=dumval1;
    MODINT(i,2)=dumval1;
    MODINT(i,3)=dumval2;
end

% count of intersections
citsc = 0;

% compute trishear fold

% for each bed
for i=1:size(yul,2)
    % make bed y's
	for j=1:size(xpf,2)
		ypf(j)=yul(i)+(xpf(j)-xpf(1))*tan(tetau);
	end
	% deform bed
	[xpd,ypd] = deformbed(xpf,ypf,xtfr,ytfr,ramp,ps,tra,slip,sinc);
    % Compute intersections of deformed bed with profile
    % for each profile segment
    for j=1:size(xprof,2)-1
        x1 = xprof(j);
        x2 = xprof(j+1);
        y1 = yprof(j);
        y2 = yprof(j+1);
        % find x limits where profile segment is
        if x1 <= x2
            left = x1;
            right = x2;
        else
            left = x2;
            right = x1;
        end
        % find y limits where profile segment is
         if y1 <= y2
            bottom = y1;
            top = y2;
        else
            bottom = y2;
            top = y1;
        end
        % for each bed segment
        for k=1:size(xpd,2)-1
            % Find if the bed segment is partially or totally
            % within the box defined by the profile segment
            cins = 0;
            if (xpd(k) >= left && xpd(k) <= right)
                if (ypd(k) >= bottom && ypd(k) <= top)
                    cins = cins + 1;
                end
            end
            if (xpd(k+1) >= left && xpd(k+1) <= right)
                if (ypd(k+1) >= bottom && ypd(k+1) <= top)
                    cins = cins + 1;
                end
            end
            % if bed segment is within box defined by profile segment,
            % it is worth estimating intersection
            if (cins > 0)
                % compute intersection: Paul Bourke's equation
                x3 = xpd(k);
                x4 = xpd(k+1);
                y3 = ypd(k);
                y4 = ypd(k+1);
                denom = (y4-y3)*(x2-x1)-(x4-x3)*(y2-y1);
                if (abs(denom) > 0.0)
                    ua = ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))/denom;
                    if (ua >= 0.0 && ua <= 1.0)
                        ub = ((x2-x1)*(y1-y3)-(y2-y1)*(x1-x3))/denom;
                        % if the lines intersect
                        if (ub >= 0.0 && ub <= 1.0)
                        	% compute intersection
                            xitsc = x1 + ua*(x2-x1);
                            yitsc = y1 + ua*(y2-y1);
                            dipitsc = atan((y4-y3)/(x4-x3));
                            % increase intersection count
                            citsc = citsc + 1;
                            % update modeled intersections
                            MODINT(citsc,1)=xitsc;
                            MODINT(citsc,2)=yitsc;
                            MODINT(citsc,3)=dipitsc;
                        end
                    end
                end
            end
        end
    end
end


% sort modeled intersections in x
SORTMODINT = sortrows(MODINT);

% Compute error
% Note: If error in location or dip is lower than 0.01, make error equal to zero
chisq = 0.0;
for i=1:size(SORTINT,1)
    if abs(SORTINT(i,1)-SORTMODINT(i,1)) < 0.01
		errx = 0.0;
   	else
   		errx = (SORTINT(i,1)-SORTMODINT(i,1))^2/(SORTINT(i,1)+SORTMODINT(i,1));
    end
    if abs(SORTINT(i,2)-SORTMODINT(i,2)) < 0.01
    	erry = 0.0;
   	else
   		erry = (SORTINT(i,2)-SORTMODINT(i,2))^2/(SORTINT(i,2)+SORTMODINT(i,2));
    end
    % Add pi to actual or modeled intersections that are dipping to the right (negative dip)
    if SORTINT(i,3) < 0.0
        SORTINT(i,3) = pi + SORTINT(i,3);
    end
    if SORTMODINT(i,3) < 0.0
        SORTMODINT(i,3) = pi + SORTMODINT(i,3);
    end
    % compute difference in dip in degrees
    if abs((SORTINT(i,3)-SORTMODINT(i,3))*180.0/pi) < 0.01
    	errdip = 0.0;
    else
    	errdip = ((SORTINT(i,3)-SORTMODINT(i,3))*180.0/pi)^2/((SORTINT(i,3)+SORTMODINT(i,3))*180.0/pi);
    end
    % error
    chisq = chisq + errx + erry + errdip;
end