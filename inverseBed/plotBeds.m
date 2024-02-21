% Copyright: Nestor Cardozo 2009. All rights reserved

% Plot beds: Shijia trench site, Chelungpu fault
% west central Taiwan. Lin et al. 2007
% Journal of Structural Geology 29, 1267-1280

% Pregrowth
% Bed 7
load bed7.txt;
x7 = bed7(:,1)';
y7 = bed7(:,2)';
% Bed 6
load bed6.txt;
x6 = bed6(:,1)';
y6 = bed6(:,2)';
% Bed 5
load bed5.txt;
x5 = bed5(:,1)';
y5 = bed5(:,2)';
% Growth
% Bed 4
load bedGrowth4.txt;
xG4 = bedGrowth4(:,1)';
yG4 = bedGrowth4(:,2)';
% Bed 3
load bedGrowth3.txt;
xG3 = bedGrowth3(:,1)';
yG3 = bedGrowth3(:,2)';
% Bed 2
load bedGrowth2.txt;
xG2 = bedGrowth2(:,1)';
yG2 = bedGrowth2(:,2)';
% Bed 1
load bedGrowth1.txt;
xG1 = bedGrowth1(:,1)';
yG1 = bedGrowth1(:,2)';

% Plot
% Pregrowth in black
plot(x7,y7,'k.');
hold on;
plot(x6,y6,'k.');
plot(x5,y5,'k.');
% Growth in blue
plot(xG4,yG4,'b.');
plot(xG3,yG3,'b.');
plot(xG2,yG2,'b.');
plot(xG1,yG1,'b.');
% Area where fault tip will be searched
axft = [7.5 7.5 16. 16. 7.5];
ayft = [2.5 10. 10. 2.5 2.5];
plot(axft,ayft,'r');

hold off;
axis equal;