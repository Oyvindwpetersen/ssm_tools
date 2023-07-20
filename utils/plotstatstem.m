function hbar=stat(data_stem,data_point,varargin)

%% Plot variance of 


%% Input

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'gap',[0.1 0.1],@isnumeric)
addParameter(p,'marg_h',[0.15 0.1],@isnumeric)
addParameter(p,'marg_w',[0.1 0.05],@isnumeric)
addParameter(p,'weight',[],@isnumeric)
addParameter(p,'marker',{'x' 'x' 'x' 'x' 'x'},@iscell)
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'color',[[0 0 1] ; [1 0 0] ; [0 0 0] ;  [0 0.4 0] ; [1 0 1] ])
addParameter(p,'xlabel','',@ischar)
addParameter(p,'ylabel','Error variance',@ischar)
addParameter(p,'xticklabel',{},@iscell)
addParameter(p,'displayname',{},@iscell)
addParameter(p,'displayname2',{},@iscell)
addParameter(p,'split',{},@iscell)
addParameter(p,'space',2,@isnumeric)

parse(p,varargin{1:end});

gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;
weight=p.Results.weight;
marker=p.Results.marker;
markersize=p.Results.markersize;
color=p.Results.color;
xl=p.Results.xlabel;
yl=p.Results.ylabel;
xticklabel=p.Results.xticklabel;
displayname=p.Results.displayname;
displayname2=p.Results.displayname2;
split=p.Results.split;
space=p.Results.space;

%%

% Matrix nx*npoints

data_stem_matrix=[];
for k=1:length(data_stem)
    data_stem_matrix(:,k)=data_stem{k};
end

[nx,npoints]=size(data_stem_matrix);


data_point_matrix=[];
for k=1:length(data_point)
    data_point_matrix(:,k)=data_point{k};
end

% [nx,npoints]=size(data_stem_matrix);

%%

% if isempty(marker)
%     marker
% end

if ~iscell(color)
    color_matrix=color; clear color;
    for k=1:size(color_matrix,1)
        color{k}=color_matrix(k,:);
    end
end

if isempty(xticklabel)
    for k=1:length(data_stem{1})
        xticklabel{k}=['x_' '{' num2str(k) '}'];
    end
end

if isempty(displayname)
    for k=1:length(data_stem)
        displayname{k}=['Set ' num2str(k)];
    end
end

if isempty(displayname2)
    for k=1:length(data_stem)
        displayname2{k}=['Set ' num2str(k)];
    end
end

if length(markersize)==1
    markersize=markersize*ones(1,length(data_stem));
end

if isempty(split)
    split={[1:size(data_stem_matrix,1)]};
end

if isempty(weight)
    for k=1:length(split)
        weight(k)=length(split{k});
    end
end


%%

n_sub=length(split);

% close all
figure();
ha=tight_subplot(1,n_sub,[gap],[marg_h],[marg_w],[],weight);

hbar={};
for j=1:n_sub

    idx_plot=split{j};
    nx_sub=length(idx_plot);

    axes(ha(j));
    hold on; grid on; ylog;

    x_plot_all=[1:(npoints+space)*nx_sub]; x_plot_all=x_plot_all(1:end-space)

    for k=1:npoints
        
        x_plot_k{k}=x_plot_all(k:(npoints+space):end);

        stem(x_plot_k{k},data_stem_matrix(idx_plot,k),...
            'Marker',marker{k},...
            'Color',color{k},...
            'displayname',displayname{k}...
            );

        plot(x_plot_k{k},data_point_matrix(idx_plot,k),...
            'LineStyle','none',...
            'Marker','o',...
            'Color',color{k},...
            'displayname',displayname2{k}...
            );


    end

    axistight(gca,[ 0.1],'ylog');
    xlim([x_plot_all(1)-1 x_plot_all(end)+1]);

    set(gca,'XTick',x_plot_k{1}*0.5+x_plot_k{end}*0.5,'XTickLabel',xticklabel(idx_plot));

end

tilefigs

%%

axes(ha(1));
xlabel(xl);
ylabel(yl);

legend show; 


% set(get(gca,'Xaxis'),'TickLabelInterpreter','latex');



