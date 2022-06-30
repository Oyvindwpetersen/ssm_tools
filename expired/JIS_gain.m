function [gain trace_Rk trace_Pkk trace_Ppkk] = JIS_gain(A,B,G,J,R,Q,P01,varargin)

p=inputParser;
addParameter(p,'inverse','conv',@ischar)
addParameter(p,'nt',10*1e3,@isnumeric)
addParameter(p,'plottrace','no',@ischar)

parse(p,varargin{:});
inverse = p.Results.inverse;
nt = p.Results.nt;
plottrace = p.Results.plottrace;


np=size(B,2);

Pkk_=P01;


for k=1:nt;

    %input estimation
    if strcmp(inverse,'conv')
        
        Rk =G*Pkk_*G'+R;  Rk=forcesym(Rk);
        Mk = (J'/Rk*J)\J'/Rk;  
        Ppkk=eye(np)/(J'/Rk*J);  Ppkk=forcesym(Ppkk);

    elseif  strcmp(inverse,'svd')
        Rk =G*Pkk_*G'+R;  Rk=forcesym(Rk);
        a2=svdinv(Rk,1e-12);
        a1=svdinv(J'*a2*J,1e-12);
        Mk = a1*J'*a2;
        Ppkk=a1; Ppkk=forcesym(Ppkk);
    end
    
    %measurement update
    Lk=Pkk_*G'/Rk;
    Pkk = Pkk_ - Lk*(Rk-J*Ppkk*J')*Lk'; Pkk=forcesym(Pkk);
    Pxpkk=-Lk*J*Ppkk;
    Ppxkk=Pxpkk';

    %time update
    Pkk_ = [A B]*[ Pkk Pxpkk;Ppxkk Ppkk]*[A';B']+Q; Pkk_=forcesym(Pkk_);
    
    trace_Rk(k)=trace(Rk);
    trace_Pkk(k)=trace(Pkk);
    trace_Ppkk(k)=trace(Ppkk);

end

gain.Rk=Rk;
gain.Lk=Lk;
gain.Mk=Mk;
gain.Ppkk=Ppkk;
gain.Pkk=Pkk;
gain.Pkk_=Pkk_;

%% plot trace

if strcmpi(plottrace,'yes')

figure(); 
ha = tight_subplot( 3,1,[.08 .05],[.05 .05],[.05 .05]);

axes(ha(1));
plot(trace_Rk);
grid on;
set(gca,'Yscale','log');

ti=title('Trace Rk');
set(ti, 'FontSize', 10);

axes(ha(2));
plot(trace_Pkk);
grid on;
set(gca,'Yscale','log');

ti=title('Trace Pkk');
set(ti, 'FontSize', 10);

axes(ha(3));
plot(trace_Ppkk);
grid on;
set(gca,'Yscale','log');

ti=title('Trace Ppkk');
set(ti, 'FontSize', 10);

set(gcf,'Units', 'Normalized', 'OuterPosition', [0.05 0.05 0.6 0.8]);

end
