% Copyright: Nestor Cardozo, 2009 ....
% These scripts are intended for academic use
% To use the scripts for commercial purposes, please contact the author

function isoscal

global x_pos y_pos
global x_pix y_pix
global kmw kml
global flagima
global flagscal

if flagima==0
   
   msgbox('No Image');
   
elseif flagima==1
    questdlg('Click the upper left and lower rigtht ends of the image','Continue Operation','OK','OK');
   
    [x_pos,y_pos] = ginput(2);
	x_pix = x_pos(2)-x_pos(1);
	y_pix = y_pos(2)-y_pos(1);

	prompt = {'Enter the horizontal extension of image in km:','Enter vertical extension of image in km'};
	def = {'45','50'};
	lineNo=1;
	title = 'SCALE IMAGE';
	answer=inputdlg(prompt,title,lineNo,def);

	kmw = str2double(answer{1});
	kml = str2double(answer{2});

	axis off;

	flagscal=1;

end