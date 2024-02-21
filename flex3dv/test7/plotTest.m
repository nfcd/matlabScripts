% Copyright Nestor Cardozo 2009

% plot deflections produced by test7.m

% number of points in x and y
pointsx = 201;
pointsy = 201;
% grid size in kilometers
delta = 5.0;

XG=zeros(pointsy,pointsx);
YG=zeros(pointsy,pointsx);
WG=zeros(pointsy,pointsx);

load deflection.txt;

count = 1;
for i=1:pointsy
    for j=1:pointsx
        XG(i,j) = (j-1)*delta-250;
        YG(i,j) = (i-1)*delta-250;
        WG(i,j) = deflection(count);
        count = count + 1;
    end
end

figure('Name','Downward displacement in meters','NumberTitle','off');
[cs,h] = contour(XG,YG,WG);
clabel(cs,h);
axis equal;
axis([0 500 0 500]);
grid on;
xlabel('km');
ylabel('km');