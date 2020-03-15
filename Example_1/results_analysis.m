clear;clc

% Load saved results
load results
chain_kzs = x_kzs';
chain_zs = x_zs';
x_zs =  permute(reshape(chain_zs',[Nx,N,T]),[3,1,2]);
x_kzs =  permute(reshape(chain_kzs',[Nx,N,T]),[3,1,2]);
para = {'\itC\rm_{max}','\itb\rm_{exp}','\it\alpha','\itR\rm_s','\itR\rm_q',...
    'Error model parameter 1','Error model parameter 2'};
la_str = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)'};
x_ti = [0 2000 4000 6000 8000];
x_ti_label = {'0','2,000','4,000','6,000','8,000'};
if isempty(sd)
    I = 7;a1 = 3;b1 = 3;
else
    I = 5;a1 = 2;b1 = 3;
end

% Trace plots of model parameters obtained by dream_zs
figure('Color',[1 1 1]);
for i = 1:I
    subplot(a1,b1,i,'FontWeight','norm','FontSize',14);
    plot(1:N*T,chain_zs(:,i),'linestyle','none','marker','.','markersize',16);
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
for i = 1:I
    subplot(a1,b1,i,'FontWeight','norm','FontSize',14);
    plot(1:N*T,chain_kzs(:,i),'linestyle','none','marker','.','markersize',16);
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
for  i = 1:I
    subplot(a1,b1,i,'FontWeight','demi','FontSize',13);
    [xp,xx] = ksdensity(chain_zs(end-4000:end,i));
    [xp1,xx1] = ksdensity(chain_kzs(end-4000:end,i));
    plot(xx,xp,'linewidth',2,'color','r');hold on;
    plot(xx1,xp1,'linestyle','--','linewidth',2,'color','b');
    x_min = min([min(xx) min(xx1)]);
    x_max = max([max(xx) max(xx1)]);
    p_max = max([max(xp) max(xp1)]);
    plot([xreal(i) xreal(i)],[0 p_max*1.1],'linewidth',2,'color','k');
    axis([x_min - (x_max-x_min)*0.1 x_max + (x_max-x_min)*0.1 0 p_max*1.1]);
    text(x_max - (x_max-x_min)*0.15,p_max*0.9,la_str{i},'FontWeight','norm','FontSize',12);
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
ylabel('$\bf\it\hat{R}\bf\rm-statistics$','interpreter','latex','fontsize',13);
plot([0 N*T],[1.2 1.2],'k--','linewidth',1.5,'color','r');
axis([0 N*T 0.8 3.5]);
set(gca,'xtick',x_ti);
set(gca,'xticklabel',x_ti_label);
set(gca,'ytick',[1 1.5 2 2.5 3.0 3.5]);
set(gca,'yticklabel',{'1.0','1.5','2.0','2.5','3.0','3.5'});
text(7000,3.3,'(a)','fontsize',13)
title('DREAM_{(ZS)}','fontsize',10)
subplot(1,2,2,'FontWeight','demi','FontSize',13);
[nn,R_stat,~] = Convergence(x_kzs);
plot(nn,R_stat,'linewidth',2); hold on;
xlabel('Number of model evaluations','interpreter','latex','fontsize',13);
ylabel('$\bf\it\hat{R}\bf\rm-statistics$','interpreter','latex','fontsize',13);
plot([0 N*T],[1.2 1.2],'k--','linewidth',1.5,'color','r');
axis([0 N*T 0.8 3.5]);
set(gca,'xtick',x_ti);
set(gca,'xticklabel',x_ti_label);
set(gca,'ytick',[1 1.5 2 2.5 3.0 3.5]);
set(gca,'yticklabel',{'1.0','1.5','2.0','2.5','3.0','3.5'});
text(7000,3.3,'(b)','fontsize',13)
title('DREAM_{(KZS)}','fontsize',10)
