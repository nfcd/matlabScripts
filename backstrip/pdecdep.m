function pdecdep

%plot decompacted depth versus time

global zt zb
global agetop agebas
global Y2 
global flagdecom

if flagdecom==0
   warndlg('Sediments are not decompacted','WARNING');
else

a = size(agebas,1);
ageplot=zeros(1,a+1);
depstrat=zeros(1,a+1);

for i=1:a+1
   if i==a+1
      ageplot(i)=-agetop(i-1);
   else
      ageplot(i)=-agebas(i);
   end
end

for j=1:a
   
   for i=1:a+1
      if i==1
         depstrat(j,i)=zt(a);
      else
         depstrat(j,i)=Y2(j,i-1)+zt(a);
      end
   end
   
   plot(ageplot,depstrat(j,:),'.-');
   hold on;
   
end

axis([-agebas(1) -agetop(a) zt(a) zb(1)]);
axis ij;
xlabel('Age (Ma)');
ylabel('Decompacted depth(km)');

hold off
   
end

         
      
      
   
   
   
      
   