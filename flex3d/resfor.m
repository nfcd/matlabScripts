% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function resfor

global densup grav
global flagrav
global flagdat

prompt = {'Enter the density of underlying media in kg/m^3','Enter Gravity in m/s^2'};
def = {'3300','9.8'};
title = 'RESTORING FORCE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
densup = str2double(answer{1});
grav = str2double(answer{2})*1e-3;

flagdat=1;
flagrav=1;