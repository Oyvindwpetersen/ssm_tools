function D=diag3d(M,n)

%%

if nargin==1
    n=[];
end


if isempty(n)
    
    if size(M,3)>1
        
        %     D=zeros(size(M,2),size(M,3));
        %     for k=1:size(M,3)
        %     D(:,k)=diag(M(:,:,k));
        %     end
        
        D=zeros(size(M,2),size(M,3));
        for k=1:size(M,1)
            D(k,:)=M(k,k,:);
        end
        
    else
        
        D=zeros(size(M,1),size(M,1),size(M,2));
        for k=1:size(M,1)
            D(k,k,:)=M(k,:);
        end
        
    end
    
end


%% M is a vector, D is M diagonal repeated n times on the third dim

if ~isempty(n)
    
    if ~isvector(M)
        error('M must be vector');
    end
    
	D=zeros(length(M),length(M),n);
	for k=1:length(M)
            D(k,k,:)=M(k);
	end    
    
    
end

