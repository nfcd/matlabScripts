function porpar

global c surpor
global flagporo

prompt = {'Compaction constant of each unit in 1/km','Surface porosity of each unit'};
def = {'0.5','0.5'};
title = 'POROSITY PARAMETERS, FROM OLDEST TO YOUNGEST';
lineNo = 10;
answer=inputdlg(prompt,title,lineNo,def);
c  = str2num(answer{1});
surpor = str2num(answer{2});

flagporo = 1;
	
   





	
   



