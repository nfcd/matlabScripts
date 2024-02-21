% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function isomap

global canvas
global ISO
global flagima
global flagdat

[filein,pathin]=uigetfile('*.jpg');
set(canvas,'Pointer','watch');
ISO = imread([pathin,filein]);
image(ISO)
axis equal
xlabel('Pixels');
ylabel('Pixels');
set(canvas,'Pointer','arrow');

flagdat=1;
flagima=1;