% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function modresa

if flagdat==0
   
   msgbox('No Data');
   
elseif flagdat==1

   id = inputdlg('File Name');

	a = str2mat(id);

	eval(['save ', a]);
   
end