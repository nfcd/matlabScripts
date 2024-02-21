function ptecwow

global zt zb
global S Ym
global agetop agebas
global flagsedwe

if flagsedwe==0
   warndlg('Decompaction or-and weight of sediments correction are missing','WARNING');
else

a=size(agetop,1);
Splot=zeros(1,a+1);
Yplot=zeros(1,a+1);
ageplot=zeros(1,a+1);
Dplot=zeros(1,a+1);

for i=1:a+1
   if i==1
      Splot(i)=zt(a);
      Yplot(i)=zt(a);
      ageplot(i)=-agebas(i);
   else
      Splot(i)=S(i-1)+zt(a);
      Yplot(i)=Ym(i-1)+zt(a);
      ageplot(i)=-agetop(i-1);
   end
end

for i=1:a+1
   if i==a+1
      Dplot(i)=zb(1);
   else
      Dplot(i)=zb(1)-zb(i);
   end
end

plot(ageplot, Dplot, 'kx-',ageplot,Splot,'r+-',ageplot,Yplot,'b.-');
legend('Compacted','Decompacted','After sediment load correction');
axis([-agebas(1) -agetop(a) Splot(1) Splot(a+1)]);
axis ij;

xlabel('Age (Ma)');
ylabel('Decompacted depth(km)');

end
