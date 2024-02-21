function water

global pm pw
global zt
global Ym
global Yplot
global Ywmax Ywmin Ywplot
global e
global flagsedwe
global flagwater

if flagsedwe ~= 1
   warndlg('Decompaction or-and sediment load correction are missing','WARNING');
else
   
a=size(zt,1);

for i=1:a+1
   if i==1
      Yplot(i)=zt(a);
   else
      Yplot(i)=Ym(i-1)+zt(a);
   end
end

%correcting for changes in sea level and water depth

prompt = {'Enter sea level in m'};
def = {'100'};
title = 'SEA LEVEL, FOR NUMBER OF UNITS + 1';
lineNo = a+1;
answer=inputdlg(prompt,title,lineNo,def);
sea  = str2num(answer{1})/1e3;

prompt = {'Enter water depth in m','Enter the error in your estimate in m'};
def = {'100','10'};
title = 'WATER DEPTH, FOR NUMBER OF UNITS + 1';
lineNo = a+1;
answer=inputdlg(prompt,title,lineNo,def);
watdep  = str2num(answer{1})/1e3;
errwat = str2num(answer{2})/1e3;

deltasea=zeros(1,a+1);
deltawatd=zeros(1,a+1);
deltamaxw=zeros(1,a+1);
deltaminw=zeros(1,a+1);

for i=1:a+1
   if i==1
      deltasea (i)=0.0;
      deltawatd (i) = 0.0;
      deltamaxw (i) = 0.0;
      deltaminw (i) = 0.0;
   else
      deltasea(i) = sea(i)-sea(i-1);
      deltawatd(i) = watdep(i) - watdep(i-1);
      deltamaxw (i) = (watdep(i) + errwat(i))- (watdep(i-1) + errwat(i-1));
      deltaminw (i) = (watdep(i) - errwat(i))- (watdep(i-1) - errwat(i-1));
   end
end

Ywmax =  Yplot - deltasea*(pw/(pm-pw)) + (deltamaxw-deltasea);     
Ywmin =  Yplot - deltasea*(pw/(pm-pw)) + (deltaminw-deltasea); 
Ywplot = Yplot - deltasea*(pw/(pm-pw)) + (deltawatd-deltasea);
e = abs(Ywmax-Ywmin);

flagwater = 1;

end






