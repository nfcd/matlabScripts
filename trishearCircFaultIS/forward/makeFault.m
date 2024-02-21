function [kl,ka,fd,R,nincs,sincs] = makeFault(ccx,ccy,ccr,maxarc,ps,slip,ashear,sinc)
% Compute fault geometry and return necessary parameters
% 
% Input: ccx, ccy, ccr = location and radius of curvature of fault; 
% maxarc = maximum central angle of fault, ps = propagation to slip ratio, 
% slip = fault slip, ashear = angle of shear in backlimb
% and sinc = fault slip increment,
%
% Output: kl = location of fault kinks, ka = angle of fault kinks with 
% the horizontal, fd = dip of fault segment, R = change of slip across 
% kink, nincs = number of increments on fault segment, sincs = slip 
% increment on fault segment

% Initialize arrays
kl = zeros(100,2); % x and y location of fault kinks
ka = zeros(100,1); % fault kink orientation
fd = zeros(100,1); % dip of fault segments
R = zeros(100,1); % Change of slip across fault kinks

% Increment central angle by 2.0 deg
iarc = 2.0*pi/180.0;

% location of first kink and dip of first fault segment
kl(1,:)=[ccx ccy-ccr];
fd(1)=0.0;
% In the lowest fault segment, fault slip is entering slip
R(1) = 1.0; 

% next kinks and fault segments
carc = iarc; % central angle of first fault segment
i = 1; % counter
while carc <= maxarc
    % increase count
    i = i + 1;
    % x location of kink
    kl(i,1)=ccx+ccr*sin(carc);
    % y location of kink
    kl(i,2)=ccy-ccr*cos(carc);
    % fault dip 
    fd(i)=atan((kl(i,2)-kl(i-1,2))/(kl(i,1)-kl(i-1,1)));
    % orientation of kink (inclined shear)
    ka(i-1)=pi/2.0-ashear; 
    % Eq. 7 of Hardy (1995) applied to his Fig. 4b
    % R = 1.0/(sin(teta)/(df/dx)+cos(teta))
    R(i)=1.0/(sin(fd(i))/tan(ka(i-1))+cos(fd(i)));
    % increase central angle
    carc = carc + iarc;
end

% Get rid of empty elements
if i < 100
    kl(i+1:100,:)=[];
    ka(i+1:100)=[];
    fd(i+1:100)=[];
    R(i+1:100)=[];
end

% Compute ninc and sinc in each fault segment
fs = size(fd,1)-1;
nincs = zeros(fs,1);
sincs = zeros(fs,1);
rslip = abs(slip); % remaining fault slip
count = 0; % a counter

for i=1:fs
    % length of fault segment
    lfs = sqrt((kl(i+1,1)-kl(i,1))*(kl(i+1,1)-kl(i,1))+(kl(i+1,2)-kl(i,2))*(kl(i+1,2)-kl(i,2)));
    % slip in fault segment
    sfs = lfs/(ps*R(i+1));
    % If remaning slip
    if rslip > 0.0
        % Increase count
        count = count + 1;
        % If remaining slip less than slip on fault segment
        % or last fault segment
        if rslip < sfs || i == fs
            sfs = rslip;
            % Update location of fault tip
            kl(i+1,1) = kl(i,1)+sfs*ps*R(i+1)*cos(fd(i+1));
            kl(i+1,2) = kl(i,2)+sfs*ps*R(i+1)*sin(fd(i+1));
        end
        % number of slip increments for segment
        nincs(i) = round(sfs/abs(sinc));
        % slip increment for segment
        sincs(i) = sfs/nincs(i)*sign(sinc);
        % Update fault slip
        rslip = rslip - sfs;
    end
end

% Get rid of not wanted elements
if count < fs
    kl(count+2:fs+1,:)=[];
    ka(count+2:fs+1)=[];
    fd(count+2:fs+1)=[];
    R(count+2:fs+1)=[];
    nincs(count+1:fs)=[];
    sincs(count+1:fs)=[];
end


end

