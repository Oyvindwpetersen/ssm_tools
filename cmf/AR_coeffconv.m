function A_out=AR_coeffconv(A_in,outputtype)

%% Convert AR-coefficients between matrix/cell form
%
% Inputs:
% A_in: cell or A_in={A1,A2,A3,...} or matrix A_in(:,:,1)=A1, A_in(:,:,2)=A2, A_in(:,:,3)=A3,...
% outputtype: 'cell' or 'matrix'
%
% Outputs:
% A_out: AR-coefficients in cell or matrix form 
%
%%

if strcmpi(outputtype,'cell')
    
	if iscell(A_in)
        A_out=A_in;
    else
        
        A_out=cell(1,size(A_in,3));
        
        for k=1:size(A_in,3)
           A_out{1,k}=A_in(:,:,k); 
        end
        
    end
    
end

if strcmpi(outputtype,'matrix')
    
	if ~iscell(A_in)
        A_out=A_in;
    else
        
        A_out=zeros(size(A_in{1},1),size(A_in{1},2),length(A_in));
        
        for k=1:size(A_in,3)
           A_out(:,:,k)=A_in{k};
        end
        
    end 
    
end