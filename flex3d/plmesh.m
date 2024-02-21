% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function plmesh

global canvas
global flagfig
global nel nnel
global xc yc
global np
global elcol
global dens
global flagdisc


if flagdisc==0
   
   msgbox('GEOMETRY IS NOT DEFINED');
   
elseif flagdisc==1

%plotting the mesh

set(canvas,'Pointer','watch');

a = zeros(nnel+1,1);
b = zeros(nnel+1,1);

for i=1:1:nel
    for j = 1:1:nnel
      j1 = np(i,j);
      a(j) = xc(j1);
      b(j) = yc(j1);
    end

	 j1 = nnel+1;
     a(j1)= a(1);
     b(j1)= b(1);


   plot(a,b);
   hold on

end 

XV = zeros(nnel,1);
YV = zeros(nnel,1);

for i=1:nel
   for j=1:nnel
      id=np(i,j);
      XV(j)=xc(id);
      YV(j)=yc(id);
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

title('PLAN VIEW OF THE PLATE')
xlabel('East-West in km');
ylabel('North-South in km');
axis equal;

set(canvas,'Pointer','arrow');

flagfig=1;

end