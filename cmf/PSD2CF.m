function [tau,R]=PSD2CF(w_axis,S,n,interpol_type,varargin)

%%

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'input','onesided',@ischar)
addParameter(p,'output','onesided',@ischar)
addParameter(p,'diagonly',false,@islogical)

parse(p,varargin{:});
input_type=p.Results.input;
output_type=p.Results.output;
diagonly=p.Results.diagonly;

%%

if isempty(interpol_type)
    interpol_type='spline';
end

if any(w_axis<=0)
   error('Positive range only');
end

%%

if strcmpi(input_type,'twosided')
   error('Not implemented yet') 
end

%% Frequency axis

if w_axis(1)~=0
    w_axis=[0 w_axis];
    S=cat(3,zeros(size(S,1)),S);
end
    
w_min=0; w_max=max(w_axis);
w_axis_int=linspace(w_min,w_max,n);

w=[-flip(w_axis_int) w_axis_int(2:end) ];
f=w/(2*pi);

N=length(f);

Fs=f(end)*2;
df=Fs/(N-1);
dt=1/Fs;


%%

if diagonly==false
    
    S_int=interp1z(w_axis,S,w_axis_int,interpol_type,0);
    S_int_neg=flip(conj(permute(S_int,[2 1 3])),3);
    Sw=cat(3,S_int_neg,S_int(:,:,2:end))/2;
    Sf=Sw*(2*pi);

    R=ifft(ifftshift(Sf,3),[],3)*N*df;
    R=fftshift(R,3);

end

%%

if diagonly==true
        
    S_2d=diag3d(S);
    
    S_int=interp1(w_axis,S_2d.',w_axis_int,interpol_type).';
    S_int_neg=flip(conj(S_int),2);
    Sw=cat(2,S_int_neg,S_int(:,2:end))/2;
    Sf=Sw*(2*pi);

    R_2d=ifft(ifftshift(Sf,2),[],2)*N*df;
    R_2d=fftshift(R_2d,2);
    
	R=diag3d(R_2d);

end

%%

tau=dt*[-(N-1)/2:(N-1)/2];

if strcmpi(output_type,'onesided')
   ind_neg=tau<0;
   tau(ind_neg)=[];
   R(:,:,ind_neg)=[];
end

%%

r1=sum(sum(sum(abs(real(R)))))/numel(R);
r2=sum(sum(sum(abs(imag(R)))))/numel(R);

if r2/r1>1e-12 
   warning('Some imaginary part here')
end

R=real(R);




