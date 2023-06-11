function plotvariance(plot_data,varargin)

%% Plot variance of 


%% Input

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'gap',[0.1 0.1],@isnumeric)
addParameter(p,'marg_h',[0.15 0.1],@isnumeric)
addParameter(p,'marg_w',[0.1 0.05],@isnumeric)
addParameter(p,'marker',{'o' 'x' 'd' 's' 'v'},@iscell)
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'color',[[0 0 1] ; [1 0 0] ; [0 0 0] ;  [0 0.4 0] ; [1 0 1] ])
addParameter(p,'xlabel','',@ischar)
addParameter(p,'ylabel','Error variance',@ischar)
addParameter(p,'xticklabel',{},@iscell)
addParameter(p,'displayname',{},@iscell)

parse(p,varargin{1:end});

gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;
marker=p.Results.marker;
markersize=p.Results.markersize;
color=p.Results.color;
xl=p.Results.xlabel;
yl=p.Results.ylabel;
xticklabel=p.Results.xticklabel;
displayname=p.Results.displayname;

%%

if ~iscell(color)
    color_matrix=color; clear color;
    for k=1:size(color_matrix,1)
        color{k}=color_matrix(k,:);
    end
end


if isempty(xticklabel)
    for k=1:length(plot_data{1})
        xticklabel{k}=['x_' num2str(k)];
    end
end

if isempty(displayname)
    for k=1:length(plot_data)
        displayname{k}=['Set ' num2str(k)];
    end
end

if length(markersize)==1
    markersize=markersize*ones(1,length(plot_data));
end


%%

% if ~iscell(plot_data)
%     plot_data={plot_data};
% end

% n_sub=length(split);

% plot_data_all=[];
% for k=1:lenght(
%     plot_data_all=[plot_data_all plot_data{k}];
% end

% if isempty(split)
%    split=cell();
%    split{1}=1:length(plot_data_all);
% end


n_sub=1;

figure();
ha=tight_subplot(1,n_sub,[gap],[marg_h],[marg_w]);

for j=1:n_sub
    
    axes(ha(j));
    hold on; grid on; ylog;
    
    for k=1:length(plot_data)

        plot(plot_data{k},'LineStyle','None','Color',color{k},'Marker',marker{k},'MarkerSize',markersize(k),'DisplayName',displayname{k});

    end

end

xlabel(xl);
ylabel(yl);

legend show; 
axistight(gca,[0.05 0.05],'x','ylog');

set(gca,'XTick',[1:length(xticklabel)],'XTickLabel',xticklabel);

% set(get(gca,'Xaxis'),'TickLabelInterpreter','latex');



