function densit

global psed pm pw
global flagdens


prompt = {'Water density in kg/m^3','Mantle density in kg/m^3'};
def = {'1000.0','3300.0'};
title = 'DENSITY OF WATER AND MANTLE';
lineNo = 1;
answer=inputdlg(prompt,title,lineNo,def);
pw  = str2double(answer{1});
pm = str2double(answer{2});


prompt = {'Grain density of units in kg/m^3'};
def = {'2650.0'};
title = 'N SEDIMENTARY UNITS, FROM OLDEST TO YOUNGEST';
lineNo = 10;
answer=inputdlg(prompt,title,lineNo,def);
psed  = str2num(answer{1});

flagdens = 1;
   
	
   





	
   



