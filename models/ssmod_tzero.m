function tz=ssmod_tzero(A,B,G,J)

%% Transmission zeros for state space model
%
% Inputs:
% A: state matrix (disc)
% B: input matrix (disc)
% G: output matrix (disc)
% J: direct transmission matrix (disc)
%
% Outputs:
% tz: transmission zeros
%
%%

tz=tzero(A,B,G,J);

tz_unstable=sum(abs(tz)>1);
tz_marginalstable=sum(abs(tz)==1);
tz_stable=sum(abs(tz)<1);

tz_max=max(abs(tz));

if isempty(tz_max);
disp(['Transmission zero max: NA']);
else
disp(['Transmission zero max: ' num2str(tz_max,'%3.5f') ' (' num2str(tz_stable) ' stable) (' num2str(tz_marginalstable) ' marginal stable) (' num2str(tz_unstable) ' unstable) ']);
end

if ~isempty(tz)
figure();
h_fake=compass(0,max(abs(tz)),'--y');hold on;

ind=find(abs(tz)<1); compass(real(tz(ind)),imag(tz(ind)),'b');
ind=find(abs(tz)==1); compass(real(tz(ind)),imag(tz(ind)),'g');
ind=find(abs(tz)>1); compass(real(tz(ind)),imag(tz(ind)),'r');

set(h_fake, 'Visible', 'Off');
title('Transmission zeros');
end



disp(['rank(J) = ' num2str(rank(J)) ', cond(J) = ' num2str(cond(J),'%0.3e')]);
