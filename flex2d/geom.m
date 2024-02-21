function geom

global flagtyp
global flaggeom
global beamtype
global loadl lmax hmax
global noint
global lprox 

if flagtyp==0
   warndlg('Type of analysis is not specified','WARNING');
elseif flagtyp==1
   if beamtype==1
		prompt = {'Enter the extent of the loaded region in kilometers','Enter the extent of the analysis in kilometers','Enter the maximum height of the loads in meters','Enter the number of intervals used to discretize the load'};
		def = {'240','500','10000','120'};
		title = 'INPUT LOAD GEOMETRY';
		lineNo=1;
		answer=inputdlg(prompt,title,lineNo,def);
		loadl = str2double(answer{1})*1e3;
   	    lmax = str2double(answer{2})*1e3;
   	    hmax =str2double(answer{3});
   	    noint = str2double(answer{4});
	elseif beamtype==2
        prompt = {'Enter the extent of the loaded region in kilometers','Enter the distance from the broken end of the beam to the load, in kilometers','Enter the maximum height of the loads in meters','Enter the number of intervals used to discretize the load'};
		def = {'30','500','1000','10'};
		title = 'INPUT LOAD GEOMETRY';
		lineNo=1;
		answer=inputdlg(prompt,title,lineNo,def);
		loadl = str2double(answer{1})*1e3;
   	    lprox = str2double(answer{2})*1e3;
   	    hmax =str2double(answer{3});
   	    noint = str2double(answer{4});
   end
	flaggeom=1;
end




   
      