function p_all=mc_simple_wind(fstart,fend,A,t,n);

freq=linspace(fstart,fend,1000);

p_all=zeros(n,length(t));

for j=1:n
p=0;
for k=1:length(freq)
    Ak=A*0.001/freq(k);
    p=p+Ak*cos(2*pi*freq(k)*t+rand*2*pi);
end

p_all(j,:)=p;
end