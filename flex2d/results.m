function results

global flagsol
global beamtype 
global x
global w

if flagsol==0
   warndlg('Analysis has not been solved','WARNING');
elseif flagsol==1
    if beamtype==1
        wi=w;
    elseif beamtype==2
        wi=fliplr(w);
    end
    amp=max(wi);
    [fbh,id]=min(wi);
    foreb=-fbh;
    xb=x(id)/1000;
    %displaying numerical values in a graphical user interface
    n=sprintf('The maximum amplitude of the vertical deformation is %.4g meters. The maximum amplitude of the forebulge is %.4g meters. The distance to the maximum amplitude of the forebulge is %.0f kilometers ',amp,foreb,xb);
    msgbox(n,'MODEL RESULTS');
end