function agein

global agetop agebas
global flagage

prompt = {'Age of tops in Ma','Age of bases in Ma'};
def = {'0.0','0.0'};
title = 'N SEDIMENTARY UNITS, FROM OLDEST TO YOUNGEST';
lineNo = 10;
answer=inputdlg(prompt,title,lineNo,def);
agetop  = str2num(answer{1});
agebas = str2num(answer{2});
  
flagage = 1;	
   





	
   



