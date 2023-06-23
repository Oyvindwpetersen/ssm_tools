function plotsingval(ss,varargin)

%% Plot singular value spectrum


%% Input

p=inputParser;
p.KeepUnmatched=true;

addParameter(p,'gap',[0.1 0.1],@isnumeric)
addParameter(p,'marg_h',[0.15 0.1],@isnumeric)
addParameter(p,'marg_w',[0.1 0.05],@isnumeric)

addParameter(p,'linestyle','None')
addParameter(p,'marker','o')
addParameter(p,'markersize',3,@isnumeric)
addParameter(p,'color','k')
addParameter(p,'xlabel','Model order',@ischar)
addParameter(p,'ylabel','Singular value',@ischar)
addParameter(p,'xlimsub',[0 12],@iscell)
addParameter(p,'xlim',[0 100],@isnumeric)

parse(p,varargin{1:end});

gap=p.Results.gap;
marg_h=p.Results.marg_h;
marg_w=p.Results.marg_w;

linestyle=p.Results.linestyle;
marker=p.Results.marker;
markersize=p.Results.markersize;
color=p.Results.color;
xl=p.Results.xlabel;
yl=p.Results.ylabel;
xlimsub=p.Results.xlimsub;
xlimit=p.Results.xlim;

%%

close all

figure();
ha=tight_subplot(1,1,[gap],[marg_h],[marg_w]);

hold on; grid on; ylog;

plotopt=struct();
plotopt.linestyle=linestyle;
plotopt.marker=marker;
plotopt.markersize=markersize;
plotopt.color=color;

plot(ss,plotopt);

xlim(xlimit);
xlabel(xl);
ylabel(yl);

axistight(ha,[0 0.05],'x','ylog2');

% Small window
pos_ha=get(ha,'Position');
pos_small=[pos_ha(1)+pos_ha(3)-0.2-0.1 pos_ha(2)+pos_ha(4)-0.2-0.1 0.2 0.2];
ha_small=axes('Units','normalized','Position',pos_small);

copyaxescontent(ha,ha_small);
xlim(xlimsub);
xlabel('');
ylabel('');

axistight(ha_small,[0 0.05],'x','ylog2');


set(ha_small,'box','on');

% set(gca,'XTick',[1:length(xticklabel)],'XTickLabel',xticklabel);

% set(get(gca,'Xaxis'),'TickLabelInterpreter','latex');



