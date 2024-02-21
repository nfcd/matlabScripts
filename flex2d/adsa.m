function adsa

global flagfig

if flagfig==0
   warndlg('No Figure on display','WARNING');
elseif flagfig==1
    id = inputdlg('File Name');
    a = str2mat(id);
    eval(['print -dill ', a]);
end
