% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function elas

global emodule poisson t
global flagelas
global flagdat

prompt = {'Enter Young''s Modulus in Pa','Enter Poisson''s Ratio','Enter Elastic thickness in km'};
def = {'7e10','0.25','30'};
title = 'ELASTIC PROPERTIES';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
emodule = str2double(answer{1})*1e-6;
poisson = str2double(answer{2});
t =       str2double(answer{3});

flagelas=1;
flagdat=1;