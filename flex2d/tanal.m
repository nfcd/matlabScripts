function tanal

global beamtype densup
global flagtyp


prompt = {'Enter 1 for an infinite or 2 for a semiinfinite plate','Enter the density of the underlying astenosphere in kg/m^3'};
def = {'1','3300'};
title = 'TYPE OF PLATE';
lineNo=1;
answer=inputdlg(prompt,title,lineNo,def);
beamtype = str2double(answer{1});
densup = str2double(answer{2});

flagtyp=1;

