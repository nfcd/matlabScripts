% Copyright: Nestor Cardozo 2012

% Plot beds: Synthetic model

% Bed 4
load bed4.txt;
x4 = bed4(:,1)';
y4 = bed4(:,2)';
% Bed 5
load bed5.txt;
x5 = bed5(:,1)';
y5 = bed5(:,2)';
% Bed 6
load bed6.txt;
x6 = bed6(:,1)';
y6 = bed6(:,2)';
% Bed 7
load bed7.txt;
x7 = bed7(:,1)';
y7 = bed7(:,2)';

% Beds in black
plot(x4,y4,'k.');
hold on;
plot(x5,y5,'k.');
plot(x6,y6,'k.');
plot(x7,y7,'k.');

% Area where center of curvature will be searched
axcc = [290 290 310 310 290];
aycc = [290 310 310 290 290];
plot(axcc,aycc,'r');
hold off;
axis equal;