% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function lomodres

id = inputdlg('File Name');

a = str2mat(id);

eval(['load ', a]);