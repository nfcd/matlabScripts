function undef

global flaghede
global beamtype
global loadl lmax hmax
global noint
global lprox 
global height
global hu 
global xplotu


if flaghede==0
   warndlg('Loads are not defined','WARNING');
elseif flaghede==1
    int = loadl/(noint*1e3);
    %set limits in the plot
    if beamtype==1
        xplotu=0:int:lmax/1e3;
    elseif beamtype==2
        ldist=lprox+loadl;
        xplotu=0:int:ldist/1e3;
    end
    hu=zeros(1,length(xplotu));
    for i=1:noint
        hu(i) = height(i);
    end
    stairs(xplotu,hu,'k')
    if beamtype==1
        axis([0.0 lmax/1e3 -hmax hmax]);
    elseif beamtype==2
        axis([0.0 ldist/1e3 -hmax hmax]);
    end
    title('UNDEFORMED TOPOGRAPHY');
    xlabel('horizontal distance, km');
    ylabel('load height, m');
end 
   