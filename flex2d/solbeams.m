function solbeams

%solves for the deflection of an infinite
%beam on an elastic foundation

global flagelas
global flaghede
global flagsol
global beamtype
global rigid
global densup
global loadl lmax
global lprox
global noint
global height pload
global x
global w


if flagelas==0
   warndlg('Elastic properties of the lithosphere are not specified','WARNING');
elseif flaghede==0
   warndlg('Loads are not defined','WARNING'); 
elseif flaghede==1
    q=zeros(1,noint);
    if beamtype==1
        int = loadl/noint;
        x=0:int:lmax;
        %Estimate the flexural parameter 
        g = 9.81;
        k=densup*g;
        alfa = ((4*rigid)/(k))^(1/4);
        %Estimate the deflection 
        %set the deflection initially to zero
        w =zeros(1,length(x));
        %compute the deflections
        for i=1:noint
            q(i)=height(i)*pload(i)*g;
            for j=1:length(x)
                if x(j)<=x(i)
                    a=x(i)-x(j);
                    b=x(i+1)-x(j);
                    Da = exp(-a/alfa)*cos(a/alfa);
                    Db = exp(-b/alfa)*cos(b/alfa);
                    w(j)=w(j)+(q(i)/(2*k))*(Da-Db);
                elseif x(j)>=x(i+1)
                    a=x(j)-x(i);
                    b=x(j)-x(i+1);
                    Da = exp(-a/alfa)*cos(a/alfa);
                    Db = exp(-b/alfa)*cos(b/alfa);
                    w(j)=w(j)+(-q(i)/(2*k))*(Da-Db);
                end
      
            end
   
        end
    elseif beamtype==2
        int = loadl/noint;
        ldist=lprox+loadl;
        x=0:int:ldist;
        %Estimate the flexural parameter 
        g = 9.81;
        k=densup*g;
        alfa = ((4*rigid)/(k))^(1/4);
        %Estimate the deflection 
        %set the deflection initially to zero
        w = zeros(1,length(x));
        %compute the deflections
        for i=1:noint
            q(i)=height(i)*pload(i)*g;
            asi=ldist-(i*int);
            bsi=ldist-(i-1)*int;
            Basi=exp(-asi/alfa)*sin(asi/alfa);
            Bbsi=exp(-bsi/alfa)*sin(bsi/alfa);
            Casi=exp(-asi/alfa)*(cos(asi/alfa)-sin(asi/alfa));
            Cbsi=exp(-bsi/alfa)*(cos(bsi/alfa)-sin(bsi/alfa));
            Po=(q(i)*alfa)*((Casi-Cbsi)-(Basi-Bbsi));
            Mo=(q(i)*alfa*alfa)*((Basi-Bbsi)-0.5*(Casi-Cbsi));
            for j=1:length(x)
                if x(j)<=(ldist-i*int)
                    a=(ldist-i*int)-x(j);
                    b=(ldist-(i-1)*int)-x(j);
                    amq=x(j);
                    Da = exp(-a/alfa)*cos(a/alfa);
                    Db = exp(-b/alfa)*cos(b/alfa);
                    Aamq = exp(-amq/alfa)*(cos(amq/alfa)+sin(amq/alfa));
                    Bamq = exp(-amq/alfa)*sin(amq/alfa);
                    w(j)=w(j)+(q(i)/(2*k))*(Da-Db)+ (Po/(2*k*alfa))*Aamq + (Mo/(k*alfa*alfa))*Bamq;         
                elseif x(j)>=(ldist-(i-1)*int)
                    a=x(j)-(ldist-i*int);
                    b=x(j)-(ldist-(i-1)*int);
                    amq=x(j);
                    Da = exp(-a/alfa)*cos(a/alfa);
                    Db = exp(-b/alfa)*cos(b/alfa);
                    Aamq = exp(-amq/alfa)*(cos(amq/alfa)+sin(amq/alfa));
                    Bamq = exp(-amq/alfa)*sin(amq/alfa);
                    w(j)=w(j)+(-q(i)/(2*k))*(Da-Db)+(Po/(2*k*alfa))*Aamq + (Mo/(k*alfa*alfa))*Bamq; 
                end
      
            end
   
        end
    end
    flagsol=1;     
end