function inunits

global zt zb
global flagunits

prompt = {'Depth to tops in kilometers','Depth to bases in kilometers'};
def = {'0.0','0.0'};
title = 'N SEDIMENTARY UNITS, FROM OLDEST TO YOUNGEST';
lineNo = 10;
answer=inputdlg(prompt,title,lineNo,def);
zt  = str2num(answer{1});
zb = str2num(answer{2});
   
flagunits= 1;	
   





	
   



