function ptecsub

global zt zb
global agetop agebas
global S
global Yplot
global Ywplot
global e
global flagwater

if flagwater==0
   warndlg('Missing water depth correction','WARNING');
else

a=size(zt,1);

Splot=zeros(1,a+1);
ageplot=zeros(1,a+1);
Dplot=zeros(1,a+1);

for i=1:a+1
   if i==1
      Splot(i)=zt(a);
      ageplot(i)=-agebas(i);
   else
      Splot(i)=S(i-1)+zt(a);
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

plot(ageplot,Dplot,'kx-',ageplot,Splot,'r.-',ageplot,Yplot,'g*-');
hold on
errorbar(ageplot,Ywplot,e,'b');
legend('Compacted','Decompacted','Sediment Load corrected', 'Tectonic subsidence');
axis([-agebas(1) -agetop(a) Splot(1) Splot(a+1)]);
axis ij;
xlabel('Age (Ma)');
ylabel('Decompacted depth(km)');

hold off

end





