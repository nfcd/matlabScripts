% Copyright: Nestor Cardozo 2009. All rights reserved

function [xp,yp,zp] = makesmallp(XP,YP,ZP,bedi)

% 3D
% make xp, yp, zp vectors out of XP, YP, ZP matrices


xp = zeros(1,(size(XP,1)*size(XP,2)));
yp = zeros(1,(size(XP,1)*size(XP,2)));
zp = zeros(1,(size(XP,1)*size(XP,2)));
count = 1;

for i=1:1:size(XP,1)
    for j=1:1:size(XP,2)
        xp(count) = XP(i,j,bedi);
        yp(count) = YP(i,j,bedi);
        zp(count) = ZP(i,j,bedi);
        count = count + 1;
    end
end