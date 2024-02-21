% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function immesh

global canvas
global ISO
global flagfig
global x_pos y_pos
global kmw kml
global x_pix y_pix
global nel nnel
global nnode
global WIDTH LENGTH
global xc yc
global np
global elcol
global dens
global flagima
global flagscal
global flagdisc

if flagima==0 
   
   msgbox('No Image');
   
elseif flagscal==0
   
   msgbox('Image is not scaled');
   
elseif flagdisc==0
   
   msgbox('Geometry is not defined');
   
else
   
%plot figure

set(canvas,'Pointer','watch');
AW = (WIDTH - kmw)/2;
AL = (LENGTH - kml)/2;
left = x_pos(1)-AW*x_pix/kmw;
rigth = x_pos(1)+kmw*x_pix/kmw+AW*x_pix/kmw;
top = y_pos(1)-AL*y_pix/kml;
bottom = y_pos(1)+kml*y_pix/kml+AL*y_pix/kml;


%transforming coordinates from kilometers to pixels

xcpix = zeros(nnode,1);
ycpix = zeros(nnode,1);

for i=1:nnode
   
   xcpix(i) = (xc(i)-AW)*x_pix/kmw + x_pos(1);
   ycpix(i) = (LENGTH-AL-yc(i))*y_pix/kml + y_pos(1);
   
end
   
   
%plotting the image
   
image(ISO)
hold on
  
   
%plotting the mesh

a = zeros(nnel+1,1);
b = zeros(nnel+1,1);
XV = zeros(nnel,1);
YV = zeros(nnel,1);

for i=1:nel
    for j = 1:nnel
        j1 = np(i,j);
        a(j) = xcpix(j1);
        b(j) = ycpix(j1);
    end
    j1 = nnel+1;
    a(j1)= a(1);
    b(j1)= b(1);
    plot(a,b);
end 

for i=1:nel
    for j=1:nnel
        id=np(i,j);
        XV(j)=xcpix(id);
        YV(j)=ycpix(id);
    end
    if elcol(i)==1
        if dens(i)>2400
            fill(XV,YV,'g');
        elseif dens(i)<=2400
            fill(XV,YV,'y');
        end
    end
end

hold off

axis([left rigth top bottom]);


title('PLAN VIEW OF THE PLATE')

set(canvas,'Pointer','arrow');

axis equal;
axis off;

flagfig=1;

end