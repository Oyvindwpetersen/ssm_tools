function p=trianglepulse(t_start,t_peak,t_end,t)

%% Triangle pulse

% Inputs:
% t_start: vector with start times
% t_peak: vector with peak times
% t_end: vector with end times
% t: time vector

% Outputs:
% p: signal with peaks

%% Fill

p=zeros(size(t));

for k=1:length(t_start)

[~,ind_min1]=min(abs(t_start(k)-t));
[~,ind_min2]=min(abs(t_peak(k)-t));
[~,ind_min3]=min(abs(t_end(k)-t));

t1=t(ind_min1);
t2=t(ind_min2);
t3=t(ind_min3);

p(ind_min1:ind_min2)=linspace(0,1,ind_min2-ind_min1+1);
p(ind_min2:ind_min3)=linspace(1,0,ind_min3-ind_min2+1);

end