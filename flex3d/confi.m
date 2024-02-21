% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function confi

%plots contours of deformed surface

global canvas
global dence denre denor
global nel nnel nnode
global np
global xc
global yc
global hel
global disp
global w
global flagsol
global vp


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

wc=w*1e3;   

xx = zeros((denre+2*denor+1),(dence+2*denor+1));
yy = zeros((denre+2*denor+1),(dence+2*denor+1));
ww = zeros((denre+2*denor+1),(dence+2*denor+1));

for i=1:(denre+2*denor+1)
   for j=1:(dence+2*denor+1)
      xx(i,j)=xc((i-1)*(dence+2*denor+1)+j);
      yy(i,j)=yc((i-1)*(dence+2*denor+1)+j);
      ww(i,j)=wc((i-1)*(dence+2*denor+1)+j);
   end
end

prompt = {'Enter preferred value of contours in meters, from minimum to maximum'};
def = {'-100 -50 0 50 100'};
title = 'CONTOUR VALUES';
lineNo=10;
answer=inputdlg(prompt,title,lineNo,def);
vp = str2num(answer{1});

prompt = {'Enter minimum and maximum values of x axis in km','Enter minimum and maximum values of y axis in km'};
def = {'75 125','75 125'};
title = 'AXIS VALUES';
lineNo=2;
answer=inputdlg(prompt,title,lineNo,def);
xaxis = str2num(answer{1});
yaxis = str2num(answer{2});

[C,h]=contour(xx,yy,ww,vp','k');

axis([xaxis(1) xaxis(2) yaxis(1) yaxis(2)]);
xlabel('West-East in km');
ylabel('South-North in km');

clabel(C,h,'manual');

set(canvas,'Pointer','arrow');

end