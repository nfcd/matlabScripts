% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function discre

%input of plate dimensions, discretize and plot the mesh

global canvas
global WIDTH LENGTH 
global denor dence denre
global nel nnel nnode  
global xc 
global yc
global np
global elcol
global hel
global dens
global ff
global flagdat
global flagdisc

prompt = {'Enter the width of the plate in km','Enter the length of the plate in km','Enter the number of columns and rows in the external regions','Enter the number of columns in the central region','Enter the number of rows in the central region'};
def = {'200','200','5','20','20'};
title = 'PLATE DIMENSIONS';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
WIDTH = str2double(answer{1});
LENGTH = str2double(answer{2});
denor = str2num(answer{3});
dence = str2num(answer{4});
denre = str2num(answer{5});

%initialize basic mesh parameters and plot the mesh.

set(canvas,'Pointer','watch');


%set mesh parameters

a=0:WIDTH/(4*denor):WIDTH/4;
b=(WIDTH/4+WIDTH/(2*dence)):WIDTH/(2*dence):3*WIDTH/4;
c=(3*WIDTH/4+WIDTH/(4*denor)):WIDTH/(4*denor):WIDTH;
x=[a b c];

d=0:LENGTH/(4*denor):LENGTH/4;
e=(LENGTH/4+LENGTH/(2*denre)):LENGTH/(2*denre):3*LENGTH/4;
f=(3*LENGTH/4+LENGTH/(4*denor)):LENGTH/(4*denor):LENGTH;
y=[d e f];

columns = (denor+dence+denor);
rows = (denor+denre+denor);
nel = columns*rows;
nnel = 4;
nnode = (columns+1)*(rows+1);

%initialize geometrical arrays

xc=zeros(nnode,1);
yc=zeros(nnode,1);
np=zeros(nel,nnel);

% initialize filled elements to zero

elcol=zeros(nel,1);

% initialize heights and densities to zero

hel=zeros(nel,1);
dens=zeros(nel,1);

%initialize loads to zero

ff=zeros(nnode,1);


%set coordinates of nodes

for i=1:rows+1
   for j=1:columns+1
      xc((i-1)*(columns+1)+j)=x(j);
      yc((i-1)*(columns+1)+j)=y(i);
   end
end

%set nodal connectivity matrix

for i=1:rows
   for j=1:columns
      np((i-1)*columns+j,1)=(i-1)*(columns+1)+j;
      np((i-1)*columns+j,2)=(i-1)*(columns+1)+j+1;
      np((i-1)*columns+j,3)=(i)*(columns+1)+j+1;
      np((i-1)*columns+j,4)=(i)*(columns+1)+j;
   end
end

flagdat=1;
flagdisc=1;


set(canvas,'Pointer','arrow');