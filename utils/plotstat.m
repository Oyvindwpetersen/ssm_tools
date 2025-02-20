function plotstat(varargin)

%% Plot variance of

ind_data=nargin;
for k=1:nargin

    if ischar(varargin{k}) | isstruct(varargin{k})
        ind_data=k-1; break;
    else
        data_point{k}=varargin{k};
    end

end

%% Input

p=inputParser;

addParameter(p,'gap',[0.1 0.1],@isnumeric)
addParameter(p,'marg_h',[0.15 0.1],@isnumeric)
addParameter(p,'marg_w',[0.1 0.05],@isnumeric)
addParameter(p,'weight',[],@isnumeric)
addParameter(p,'marker',{'x' 'o' 'd' 's' '*' '.'},@iscell)
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'color',[ [0 0 1] ; [1 0 0] ; [0 0 0] ;  [0 0.4 0] ; [1 0 1] ; [ 0 0.5 0.5]])
addParameter(p,'xlabel','',@ischar)
addParameter(p,'ylabel','Error deviation',@ischar)
addParameter(p,'xticklabel',{},@iscell)
addParameter(p,'displayname',{},@iscell)
addParameter(p,'split',{},@iscell)
addParameter(p,'space',2,@isnumeric)
addParameter(p,'spread',true,@islogical)
addParameter(p,'normalize',false,@islogical)
addParameter(p,'interpreter','latex',@ischar)

parse(p,varargin{ind_data+1:end});

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
split=p.Results.split;
space=p.Results.space;
spread=p.Results.spread;
normalize=p.Results.normalize;
interpreter=p.Results.interpreter;

%%

[nsources]=length(data_point);
[nx]=length(data_point{1});

%%

if ~iscell(color)
    color_matrix=color; clear color;
    for k=1:size(color_matrix,1)
        color{k}=color_matrix(k,:);
    end
end

if isempty(xticklabel)
    for k=1:nx
        xticklabel{k}=['x_' '{' num2str(k) '}'];
    end
end

if isempty(displayname)
    for k=1:nsources
        displayname{k}=['Set ' num2str(k)];
    end
end

if length(markersize)==1
    markersize=markersize*ones(1,nsources);
end

if isempty(split)
    split={[1:nx]};
end

if isempty(weight)
    for k=1:length(split)
        weight(k)=length(split{k});
    end
end


%%

n_sub=length(split);

figure();
ha=tight_subplot(1,n_sub,gap,marg_h,marg_w,[],weight);

% hbar={};
for j=1:n_sub

    idx_plot=split{j};
    nx_sub=length(idx_plot);

    axes(ha(j));
    hold on; grid on; ylog;

    if spread
        x_plot_all=[1:(nsources+space)*nx_sub]; x_plot_all=x_plot_all(1:end-space);
    else
        x_plot_all=[1:(space)*nx_sub]; %x_plot_all=x_plot_all()
    end

    for k=1:nsources

        if spread
            x_plot_k{k}=x_plot_all(k:(nsources+space):end);
        else
            x_plot_k{k}=x_plot_all(1:(space):end);
        end

        if normalize
            norm_scale=data_point{1}(idx_plot);
        else
            norm_scale=1;
        end

        plot(x_plot_k{k},data_point{k}(idx_plot)./norm_scale,...
            'LineStyle','none',...
            'Marker',marker{k},...
            'Color',color{k},...
            'displayname',displayname{k}...
            );

    end

    axistight(gca,[ 0.1],'ylog');
    xlim([ min(x_plot_k{1})-1 max(x_plot_k{end})+1]);

    set(ha(j),'XTick',x_plot_k{1}*0.5+x_plot_k{end}*0.5,'XTickLabel',xticklabel(idx_plot));

    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex';
end

% tilefigs

%%

axes(ha(1));
xlabel(xl);
ylabel(yl);

legend show;


% set(get(gca,'Xaxis'),'TickLabelInterpreter','latex');



