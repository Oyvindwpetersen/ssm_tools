function antires_sens(omega,xi,phi,idx_d,idx_p)

if ~isvector(omega); omega=diag(omega); end
if ~isvector(xi); xi=diag(xi); end

[anti_omega]=antires(omega,xi,phi,idx_d,idx_p);

%%

scale=1.01;

anti_omega_p_all=[];

for n=1:2
    for idx_mode=1:length(omega)

        omega_p=omega;
        xi_p=xi;
        phi_p=phi;

        if n==1
            omega_p(idx_mode)=omega_p(idx_mode)*scale;

            theta=omega(idx_mode);
            delta_theta=omega_p(idx_mode)-omega(idx_mode);
        elseif n==2
            xi_p(idx_mode)=xi_p(idx_mode)*scale;

            theta=xi(idx_mode);
            delta_theta=xi_p(idx_mode)-xi(idx_mode);
        % elseif n==3
            % phi_p(:,idx_mode)=phi_p(:,idx_mode)*sqrt(scale);
        end

        [anti_omega_p]=antires(omega_p,xi_p,phi_p,idx_d,idx_p);

        delta_anti_omega=anti_omega_p-anti_omega;

        anti_omega_sens{n}(idx_mode,:)=delta_anti_omega./delta_theta;

        anti_omega_sens_norm{n}(idx_mode,:)=delta_anti_omega./delta_theta .*theta./anti_omega;
        
    end
end


% dwa/dtheta * theta/onega

figure(); sizefig('m');

% ha=tight_subplot(1,2,[0.2 0.1],[0.15 0.15],[0.05 0.05]);

% plotopt=optlib1([1 2]);

plotopt.gap=[0.2];
plotopt.marg_h=[0.15 0.15];
plotopt.marg_w=[0.1 0.15];



ha=tight_subplot(1,2,plotopt.gap,plotopt.marg_h,plotopt.marg_w);

for n=1:2
    axes(ha(n)); hold on; grid on;

    for idx=1:length(omega); xtl{idx}=['m=' num2str(idx)]; end
    for idx=1:length(anti_omega); ytl{idx}=['j=' num2str(idx)]; end

    bar3z(abs(anti_omega_sens_norm{n}).','xticklabel',xtl,'yticklabel',ytl);
    view([0 90]);

    xlabel('Mode number');
    ylabel('Anti-resonant frequency');
    
    if n==1
        var_str='\omega_m';
    elseif n==2
        var_str='\xi_m';
    elseif n==3
        var_str='\tilde{m}';
    end
    
    tit=[ '$' ... 
        '\frac{ d \omega_{a,j} }{ d ' var_str  '}' ...
        '\cdot' ...
        '\frac{ ' var_str   '}{ \omega_{a,j} }' ...
        '$' ];
    

        tit=[ '$' ... 
        '\frac{ d \omega_{a,j} /  \omega_{a,j} }{ d ' var_str  '/' var_str '}' ...
        '$' ];
    


    title(tit,'Fontsize',8,'FontWeight','normal','Interpreter','latex');

    hc=colorbar();

    hc_lim=[0 max(max(abs(anti_omega_sens_norm{n}))) ];

    set(hc,'Limits',hc_lim);


    % xtickangle(90);

    % set(hc,'TicksMode','manual');
    % set(hc,'Ticks',linspace(0,hc_lim(2),4));

    colorbarpos(hc,1,1.5,0.05);

    % caxis([0 max(max(abs(anti_omega_sens_norm{n}))) ]);

    colormap(brewermap(100,'GnBu'));

end

