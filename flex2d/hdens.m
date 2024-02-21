function hdens

global flaggeom
global flaghede
global noint
global height pload

if flaggeom==0
   warndlg('Geometry of the loads is not defined','WARNING');
elseif flaggeom==1
	prompt = {'Enter the height of the load in meters, from most distal to most proximal columns to the basin','Enter the average density of the tectonic or sedimentary columns in kg/m^3'};
	def = {'1000','2400'};
	title = 'LOAD COLUMNS';
	lineNo = noint;
	answer=inputdlg(prompt,title,lineNo,def);
	height  = str2num(answer{1});
    pload = str2num(answer{2});
	flaghede=1;
end