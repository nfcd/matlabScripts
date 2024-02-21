% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function plotin

%plots undeformed load topography

global canvas
global WIDTH LENGTH
global nel nnel 
global xc
global yc
global np
global hel
global dens
global flagload


if flagload==0
   
   msgbox('NO LOADS');
   
else

% plotting the mesh

set(canvas,'Pointer','watch');

a = zeros(nnel,1);
b = zeros(nnel,1);
h = zeros(nnel,1);
bh = zeros(nnel,1);
ah = zeros(nnel,1);
bb = zeros(nnel,1);
ab = zeros(nnel,1);
hh = zeros(nnel,1);
hb = zeros(nnel,1);
al = zeros(nnel,1);
bl = zeros(nnel,1);
ar = zeros(nnel,1);
br = zeros(nnel,1);
hl = zeros(nnel,1);
hr = zeros(nnel,1);

for i=1:1:nel
    for j = 1:1:nnel
        j1 = np(i,j);
        a(j) = xc(j1);
        b(j) = yc(j1);
        h(j) = hel(i);
    end
    
    if hel(i)>0.0
        for j=1:1:nnel
            bh(j)=b(1);
            ah(j)=a(j);
            bb(j)=b(3);
            ab(j)=a(j);
            if j<=2
               hh(j)=0.0;
               hb(j)=0.0;
            else
               hh(j)=hel(i);
               hb(j)=hel(i);
            end
            al(j)=a(1);
            bl(j)=b(j);
            ar(j)=a(2);
            br(j)=b(j);        
            if j==1 || j==3
               hl(j)=0.0;
               hr(j)=0.0;
            else
               hl(j)=hel(i);
               hr(j)=hel(i);
            end
        end
        
        if dens(i)>2400
         	fill3(a,b,h,[0 1 0]);
      		fill3(ah,bh,hh,[0 1 0]);
      		fill3(ab,bb,hb,[0 1 0]);
      		fill3(al,bl,hl,[0 1 0]);
         	fill3(ar,br,hr,[0 1 0]);
      	elseif dens(i)<=2400
         	fill3(a,b,h,[1 1 0]);
      		fill3(ah,bh,hh,[1 1 0]);
      		fill3(ab,bb,hb,[1 1 0]);
      		fill3(al,bl,hl,[1 1 0]);
            fill3(ar,br,hr,[1 1 0]);
        end
    else
        fill3(a,b,h,'w');
    end
      
      hold on   
   
end

hold off

if max(hel)==0.0
   axis([0.0 WIDTH 0.0 LENGTH -5 5]);
elseif max(hel)~=0.0
   axis([0.0 WIDTH 0.0 LENGTH 0.0 max(hel)*3]);
end
   
title(' LOAD TOPOGRAPHY');
xlabel('km');
ylabel('km');
zlabel('km');


set(canvas,'Pointer','arrow');



end