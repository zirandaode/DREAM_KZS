clear;clc

% Load saved results
load results.mat
load data.mat
chain_zs = x_zs'; y_zs = y_zs'; 
chain_kzs = x_kzs'; y_kzs = y_kzs';
x_zs =  permute(reshape(chain_zs',[Nx,N,T]),[3,1,2]);
x_kzs =  permute(reshape(chain_kzs',[Nx,N,T]),[3,1,2]);
if size(xreal,1) == 1, xreal = xreal'; end
x_ti = [0 20000 40000 60000 80000];
x_ti_label = {'0','20,000','40,000','60,000','80,000'};
y_ti = [0 0.1 0.2 0.3];
y_ti_label = {'0','0.1','0.2','0.3'};

% The reference Y fields, the dream_zs estimations and the dream_kzs estimations
ref_1 = MeanY + Fi(:,1:kl_num)*xreal(1:kl_num,1);
ref_2 = MeanY + Fi(:,1:kl_num)*xreal(kl_num+1:kl_num*2,1);
ref_3 = MeanY + Fi(:,1:kl_num)*xreal(kl_num*2+1:kl_num*3,1);
zs_1 = mean(MeanY + Fi(:,1:kl_num)*chain_zs(end-1000:end,1:kl_num)',2);
zs_2 = mean(MeanY + Fi(:,1:kl_num)*chain_zs(end-1000:end,1+kl_num:kl_num*2)',2);
zs_3 = mean(MeanY + Fi(:,1:kl_num)*chain_zs(end-1000:end,kl_num*2+1:kl_num*3)',2);
kzs_1 = mean(MeanY + Fi(:,1:kl_num)*chain_kzs(end-1000:end,1:kl_num)',2);
kzs_2 = mean(MeanY + Fi(:,1:kl_num)*chain_kzs(end-1000:end,1+kl_num:kl_num*2)',2);
kzs_3 = mean(MeanY + Fi(:,1:kl_num)*chain_kzs(end-1000:end,kl_num*2+1:kl_num*3)',2);
fg = figure('Color',[1 1 1]);
subplot(9,9,[1 2 3 10 11 12 19 20 21],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(ref_1,35,40)',20,'LineStyle','none');
title({'Reference'},'FontSize',10);
text(5,35,'Layer 1','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
ylabel('Row','FontSize',10);
set(gca,'xtick',[])
subplot(9,9,[28 29 30 37 38 39 46 47 48],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(ref_2,35,40)',20,'LineStyle','none');
text(5,35,'Layer 2','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
ylabel('Row','FontSize',10);
set(gca,'xtick',[])
subplot(9,9,[55 56 57 64 65 66 73 74 75],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(ref_3,35,40)',20,'LineStyle','none');
text(5,35,'Layer 3','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
ylabel('Row','FontSize',10);
xlabel('Column','FontSize',10);
subplot(9,9,[4 5 6 13 14 15 22 23 24],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(zs_1,35,40)',20,'LineStyle','none')
title({'DREAM_{(ZS)}'},'FontSize',10);
text(5,35,'Layer 1','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(9,9,[31 32 33 40 41 42 49 50 51],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(zs_2,35,40)',20,'LineStyle','none')
text(5,35,'Layer 2','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(9,9,[58 59 60 67 68 69 76 77 78],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(zs_3,35,40)',20,'LineStyle','none')
text(5,35,'Layer 3','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
xlabel('Column','FontSize',10);
set(gca,'ytick',[])
subplot(9,9,[7 8 9 16 17 18 25 26 27],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(kzs_1,35,40)',20,'LineStyle','none');
title({'DREAM_{(KZS)}'},'FontSize',10);
text(5,35,'Layer 1','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(9,9,[34 35 36 43 44 45 52 53 54],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(kzs_2,35,40)',20,'LineStyle','none');
text(5,35,'Layer 2','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
set(gca,'xtick',[])
set(gca,'ytick',[])
subplot(9,9,[61 62 63 70 71 72 79 80 81],'Parent',fg,'Layer','top','FontSize',10)
contourf(reshape(kzs_3,35,40)',20,'LineStyle','none');
text(5,35,'Layer 3','FontWeight','demi','FontSize',10,'LineStyle','none','Color','k');
xlabel('Column','FontSize',10);
set(gca,'ytick',[])

% RMSEs between the simulated model outputs in the chains and the measurements
rmse_zs = sqrt( mean( (y_zs-repmat(Obs',size(y_zs,1),1)).^2,2 ) );
rmse_kzs = sqrt( mean( (y_kzs-repmat(Obs',size(y_kzs,1),1)).^2,2 ) );
thin_zs = 1:5:size(chain_zs,1);
fig = figure('color',[1 1 1]);
axes1 = axes('Parent',fig,'FontSize',12);
plot(thin_zs,(rmse_zs(thin_zs)),'marker','.','color','b',...
    'linestyle','none','markersize',16);
hold on;
thin_kzs = 1:5:size(chain_kzs,1);
plot(thin_kzs,(rmse_kzs(thin_kzs)),'marker','.','color','r',...
    'linestyle','none','markersize',16);
hold off
legend('DREAM_{(ZS)}','DREAM_{(KZS)}')
xlabel('Number of model evaluations','FontSize',12);
ylabel('RMSE','FontSize',12)
set(gca,'xtick',x_ti);
set(gca,'xticklabel',x_ti_label);
set(gca,'ytick',y_ti);
set(gca,'yticklabel',y_ti_label);
axis([0 N*T 0 0.3])

if isempty(sd)
    figure('color',[1 1 1]);
    subplot(1,2,1,'FontWeight','norm','FontSize',13)
    plot(chain_zs(:,end),'marker','.','markersize',16,'linestyle','none');hold on
    plot(size(chain_zs,1),0.01,'linewidth',3,'linestyle','none',...
        'marker','x','markersize',16,'color','r');
    xlabel('Number of model evaluations','FontSize',13);
    ylabel('Error model parameter','FontSize',13);
    title('DREAM_{(ZS)}','FontSize',13)
    set(gca,'xtick',x_ti);
    set(gca,'xticklabel',x_ti_label);
    axis([0 N*T 0 0.02])
    subplot(1,2,2,'FontWeight','norm','FontSize',13)
    plot(chain_kzs(:,end),'marker','.','markersize',16,'linestyle','none');hold on
    plot(size(chain_zs,1),0.01,'linewidth',3,'linestyle','none',...
        'marker','x','markersize',16,'color','r');
    xlabel('Number of model evaluations','FontSize',13);
    ylabel('Error model parameter','FontSize',13);
    title('DREAM_{(KZS)}','FontSize',13)
    set(gca,'xtick',x_ti);
    set(gca,'xticklabel',x_ti_label);
    axis([0 N*T 0 0.02])
end

