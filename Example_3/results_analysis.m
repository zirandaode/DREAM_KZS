clear;clc

% Load saved results
load results
chain_zs = x_zs';
chain_kzs = x_kzs';
x_zs =  permute(reshape(chain_zs',[Nx,N,T]),[3,1,2]);
x_kzs =  permute(reshape(chain_kzs',[Nx,N,size(chain_kzs,1)/N]),[3,1,2]);
para = {'\itx_s','\ity_s','\itS\rm_1','\itS\rm_2','\itS\rm_3',...
    '\itS\rm_4','\itS\rm_5','\itS\rm_6'};
la_str = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','i'};

% Trace plots of the contaminant source parameters obtained by dream_kzs
figure('Color',[1 1 1]);
t_1 = 1:1:size(chain_kzs,1);
x_ti = [0 50000 100000];
x_ti_label = {'0','50,000','100,000'};
for i = 1:8
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
end
legend('Parameter samples','True value')

% Trace plots of the contaminant source parameters obtained by dream_zs
hf = figure('Color',[1 1 1]);
t_1 = 1:50:size(chain_zs,1);
x_ti = [0 500000 1000000];
x_ti_label = {'0','500,000','1,000,000'};
for i = 1:8
    subplot(3,3,i,'FontWeight','norm','FontSize',14);
    plot(t_1,chain_zs(t_1,i),'linestyle','none','marker','.','markersize',16);
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
end
legend('Parameter samples','True value')

% The estimated Y fields by dream_zs and dream_kzs
load fun.mat
field_real = 2 + fun(:,1:100)*xreal(9:end,1);
field_zs = mean(2 + fun(:,1:100)*chain_zs(end-5000:10:end,9:end)',2);
field_kzs = mean(2 + fun(:,1:100)*chain_kzs(end-5000:10:end,9:end)',2);
figure1 = figure('PaperSize',[20.98 29.68],'color',[1 1 1]);
axes('Parent',figure1,'Position',[0.3552 0.5838 0.3347 0.3412],...
    'FontWeight','norm','FontSize',12);
contourf(0:0.25:20,0:0.25:10,reshape(field_real,81,41)',25,'LineStyle','none');hold on;
title('Reference')
xlabel({'\itx\rm [\itL\rm]'},'FontSize',12);
ylabel({'\ity\rm [\itL\rm]'},'FontSize',12);
text(1,9,'(a)','FontWeight','bold',...
    'FontSize',12,'LineStyle','none','Color','k');
subplot(2,7,[8 9 10],'Parent',figure1,'FontWeight','norm','FontSize',12);
contourf(0:0.25:20,0:0.25:10,reshape(field_zs,81,41)',25,'LineStyle','none');hold on;
title('DREAM_{(ZS)}')
xlabel({'\itx\rm [\itL\rm]'},'FontSize',12);
ylabel({'\ity\rm [\itL\rm]'},'FontSize',12);
text(1,9,'(b)','FontWeight','bold',...
    'FontSize',12,'LineStyle','none','Color','k');
subplot(2,7,[12 13 14],'Parent',figure1,'FontWeight','norm','FontSize',12);
contourf(0:0.25:20,0:0.25:10,reshape(field_kzs,81,41)',25,'LineStyle','none');hold on;
title('DREAM_{(KZS)}')
xlabel({'\itx\rm [\itL\rm]'},'FontSize',12);
ylabel({'\ity\rm [\itL\rm]'},'FontSize',12);
text(1,9,'(c)','FontWeight','bold',...
    'FontSize',12,'LineStyle','none','Color','k');
