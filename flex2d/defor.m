function defor

global flagsol
global beamtype
global loadl lmax hmax
global noint
global lprox
global height
global x
global w
global wplot he xplotd

if flagsol==0
   warndlg('Problem has not been solved','WARNING');
elseif flagsol==1
    if beamtype==1
        wi=w;
    elseif beamtype==2
        wi=fliplr(w);
    end
    xplotd=x/1e3;
    he =zeros(1,length(xplotd));
    for i=1:noint
        he(i)=(height(i)-wi(i));
    end
    for i=noint+1:1:size(he,2)
        he(i)=-wi(i);
    end
    wplot = -wi;
    plot(xplotd,wplot,'k');
    hold on
    stairs(xplotd,he,'k');
    title('DEFORMED TOPOGRAPHY')
    xlabel('horizontal distance, km')
    ylabel('height m')
    if beamtype==1
        axis([0.0 lmax/1e3 -hmax hmax]);
    elseif beamtype==2
        ldist=lprox+loadl; 
        axis([0.0 ldist/1e3 -hmax hmax]);
    end
    hold off
end