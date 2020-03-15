clear;clc

% Load saved results
load results
chain_zs = x_zs';
chain_kzs = x_kzs';
x_zs =  permute(reshape(chain_zs',[Nx,N,T]),[3,1,2]);
x_kzs =  permute(reshape(chain_kzs',[Nx,N,T]),[3,1,2]);
t_1 = 1:1:size(chain_zs,1);
t_2 = 1:1:size(chain_zs,1);
x_ti = 0:size(chain_zs,1)/3:size(chain_zs,1);
x_ti_label = {'0','8,000','16,000','24,000'};
para = {'\itI\rm_{max}','\itS\rm_{max}','\itQ\rm_{max}','\it\alpha\rm_E',...
    '\it\alpha\rm_F','\itK\rm_F','\itK\rm_S','Error model parameter 1',...
    'Error model parameter 2'};
la_str = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)'};

% Trace plots of model parameters obtained by dream_zs
figure('Color',[1 1 1]);
for i = 1:Nx
    subplot(3,3,i,'FontWeight','norm','FontSize',14);
    plot(t_2,chain_zs(t_2,i),'linestyle','none','marker','.','markersize',16);
    hold on;
    plot(size(chain_zs,1),xreal(i),'linewidth',3,'linestyle','none',...
        'marker','x','markersize',15,'color','r');
    xlabel('Number of model evaluations','FontWeight','norm','FontSize',14);
    ylabel(para{i},'FontWeight','norm','FontSize',14)
    set(gca,'xtick',x_ti);
    set(gca,'xticklabel',x_ti_label);
    axis([0 size(chain_zs,1) range(i,1) range(i,2)]);
    text(size(chain_zs,1)*0.8,range(i,1)+(range(i,2)-range(i,1))*0.8,...
        la_str{i},'FontWeight','norm','FontSize',15);
    if i == 2
        title('DREAM_{(ZS)}');
    end
end
legend('Parameter samples','True value')

% Trace plots of model parameters obtained by dream_kzs
figure('Color',[1 1 1]);
for i = 1:Nx
    subplot(3,3,i,'FontWeight','norm','FontSize',14);
    plot(t_1,chain_kzs(t_1,i),'linestyle','none','marker','.','markersize',16);
    hold on;
    plot(size(chain_kzs,1),xreal(i),'linewidth',3,'linestyle','none',...
        'marker','x','markersize',15,'color','r');
    xlabel('Number of model evaluations','FontWeight','norm','FontSize',14);
    ylabel(para{i},'FontWeight','norm','FontSize',14)
    set(gca,'xtick',x_ti);
    set(gca,'xticklabel',x_ti_label);
    axis([0 size(chain_kzs,1) range(i,1) range(i,2)]);
    text(size(chain_kzs,1)*0.8,range(i,1)+(range(i,2)-range(i,1))*0.8,...
        la_str{i},'FontWeight','norm','FontSize',15);
    if i == 2
        title('DREAM_{(KZS)}');
    end
end
legend('Parameter samples','True value')

% Marginal PPDFs obtained by dream_kzs and dream_zs
figure('Color',[1 1 1]);
for  i = 1:Nx
    subplot(3,3,i,'FontWeight','demi','FontSize',13);
    [xp,xx] = ksdensity(chain_zs(end-8000:end,i));
    [xp1,xx1] = ksdensity(chain_kzs(end-8000:end,i));
    plot(xx,xp,'linewidth',2,'color','r');hold on;
    plot(xx1,xp1,'linestyle','--','linewidth',2,'color','b');
    x_min = min([min(xx) min(xx1)]);
    x_max = max([max(xx) max(xx1)]);
    p_max = max([max(xp) max(xp1)]);
    plot([xreal(i) xreal(i)],[0 p_max*1.1],'linewidth',2,'color','k');
    axis([x_min - (x_max-x_min)*0.1 x_max + (x_max-x_min)*0.1 0 p_max*1.1]);
    text(x_max - (x_max-x_min)*0.2,p_max*0.9,la_str{i},'FontWeight','norm','FontSize',12);
    xlabel(para{i},'FontWeight','norm','FontSize',10);
    ylabel('Marginal PPDF','FontWeight','norm','FontSize',10);
end
legend('DREAM_{(ZS)}','DREAM_{(KZS)}','True value')

% Convergence analysis
figure('color',[1 1 1]);
subplot(1,2,1,'FontWeight','demi','FontSize',13)
[nn,R_stat,~] = Convergence(x_zs);
plot(nn,R_stat,'linewidth',2); hold on;
xlabel('Number of model evaluations','interpreter','latex','fontsize',13);
ylabel('$\bf\it\hat{R}\bf\rm-statistic$','interpreter','latex','fontsize',13);
plot([0 N*T],[1.2 1.2],'k--','linewidth',1.5,'color','r');
axis([0 N*T 0.8 3.5]);
set(gca,'xtick',x_ti);
set(gca,'xticklabel',x_ti_label);
text(19000,3.2,'(a)','fontsize',13)
title('DREAM_{(ZS)}','fontsize',10)
subplot(1,2,2,'FontWeight','demi','FontSize',13);
[nn,R_stat,~] = Convergence(x_kzs);
plot(nn,R_stat,'linewidth',2); hold on;
xlabel('Number of model evaluations','interpreter','latex','fontsize',13);
ylabel('$\bf\it\hat{R}\bf\rm-statistic$','interpreter','latex','fontsize',13);
plot([0 N*T],[1.2 1.2],'k--','linewidth',1.5,'color','r');
axis([0 N*T 0.8 3.5]);
set(gca,'xtick',x_ti);
set(gca,'xticklabel',x_ti_label);
text(19000,3.2,'(b)','fontsize',13)
title('DREAM_{(KZS)}','fontsize',10)
