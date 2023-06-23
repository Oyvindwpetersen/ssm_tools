function plotcovmatrix(C,XTickLabel,YTickLabel,xlab,ylab)


%%

[n1,n2]=size(C);

if isempty(YTickLabel)
    for k=1:n1
        YTickLabel{k}=['X' num2str(k)];
    end
end

if isempty(XTickLabel)
    for k=1:n2
        XTickLabel{k}=['X' num2str(k)];
    end
end

figure();
ha = tight_subplot(1,1,[.1 .1],[.05 .05],[.05 .05]);
bar3z(C,...
'xlabel',xlab,'XTickLabel',XTickLabel,...
'ylabel',ylab,'YTickLabel',YTickLabel);
view([0 90]);

% set(gca,'TickLabelInterpreter','latex'); 
xtickangle(90);
cmap=colormap(brewermap(20,'PuBu')); %'OrRd'
axis image;
colorbar

for k1=1:n1
for k2=1:n2
text(k2,k1,max(max(C)),num2str(C(k1,k2),'%0.2f'),'Color',[1 0.5 0],'HorizontalAlignment','center','FontSize',6,'BackGroundColor','none','Margin',1,'Interpreter','latex');
end
end

% plotname=['cross_mode' num2str(modeNo(end)) ]; %'_modenorm'

% plotscriptmain('h',8,'w',10,'name',plotname,'path','fig2','labelsize',6,'ticksize',6,'legendsize',6,'titlesize',4,'box','on','format',{'pdf' 'jpg'});
