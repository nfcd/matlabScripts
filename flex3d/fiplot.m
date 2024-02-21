function fiplot
% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

%plots final geometry of basin and orogen

global canvas
global WIDTH LENGTH
global dence denre denor
global nel nnel nnode
global np
global xc
global yc
global disp
global hel
global w
global flagsol

if flagsol==0
   
   msgbox('Problem is not solved');

else
   
set(canvas,'Pointer','watch');


%getting the deformed topography
w = zeros(nnode,1);
divfac = zeros(nnode,1);
for i=1:nel
    for j=1:nnel
        j1=np(i,j);
        if hel(i)>0
            w(j1)=w(j1)+(hel(i)-disp(j1));
            divfac(j1)=divfac(j1)+1;
        end
    end
end

for i=1:nnode
    if divfac(i)==0
        w(i)=-disp(i);
        divfac(i)=1;
    end
end

w = w./divfac;

%contour the data

xx = zeros((denre+2*denor+1),(dence+2*denor+1));
yy = zeros((denre+2*denor+1),(dence+2*denor+1));
dd = zeros((denre+2*denor+1),(dence+2*denor+1));
ww = zeros((denre+2*denor+1),(dence+2*denor+1));

for i=1:(denre+2*denor+1)
   for j=1:(dence+2*denor+1)
      xx(i,j)=xc((i-1)*(dence+2*denor+1)+j);
      yy(i,j)=yc((i-1)*(dence+2*denor+1)+j);
      dd(i,j)=-disp((i-1)*(dence+2*denor+1)+j);
      ww(i,j)=w((i-1)*(dence+2*denor+1)+j);
   end
end

surfl(xx,yy,ww); %surf plot
shading interp;
colormap pink;

if max(hel)==0.0
   axis([0.0 WIDTH 0.0 LENGTH -5 5]);
elseif max(hel)~=0.0
   axis([0.0 WIDTH 0.0 LENGTH -max(hel)*3/2 max(hel)*3/2]);
end

xlabel('km');
ylabel('km');
zlabel('km');
title('DEFORMED TOPOGRAPHY');

figure;

surfl(xx,yy,dd); %surf plot
shading interp;
colormap pink;

if max(hel)==0.0
   axis([0.0 WIDTH 0.0 LENGTH -5 5]);
elseif max(hel)~=0.0
   axis([0.0 WIDTH 0.0 LENGTH -max(hel)*3/2 max(hel)*3/2]);
end

xlabel('km');
ylabel('km');
zlabel('km');
title('DISPLACEMENT');

set(canvas,'Pointer','arrow');

end