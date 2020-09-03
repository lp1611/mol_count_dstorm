%load sims_3state_slow.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_3state_slow.cum_mode],'Normalization','probability');
hist([sims_3state_slow.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([70,130])
hold off 
xlabel({'M'},'FontSize',18);
ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_3state_slow.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_3state_slow.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_3state_slow.pdf');
%export_fig post_dens_3state_slow.png -transparent

%load sims_3state_med.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_3state_med.cum_mode],'Normalization','probability');
hist([sims_3state_med.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([60,140])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_3state_med.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_3state_med.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_3state_med.pdf');
%export_fig post_dens_3state_med.png -transparent

%load sims_3state_fast.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_3state_fast.cum_mode],'Normalization','probability');
hist([sims_3state_fast.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([60,140])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_3state_fast.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_3state_fast.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_3state_fast.pdf');
%export_fig post_dens_3state_fast.png -transparent

%load sims_4state_slow.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_4state_slow.cum_mode],'Normalization','probability');
hist([sims_4state_slow(1:9.7e3).cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([80,120])
hold off 
xlabel({'M'},'FontSize',18);
ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_4state_slow.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_4state_slow.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_4state_slow.pdf');
%export_fig post_dens_4state_slow.png -transparent

%load sims_4state_med.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_4state_med.cum_mode],'Normalization','probability');
hist([sims_4state_med.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([60,140])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_4state_med.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_4state_med.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_4state_med.pdf');
%export_fig post_dens_4state_med.png -transparent

%load sims_4state_fast.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_4state_fast.cum_mode],'Normalization','probability');
hist([sims_4state_fast.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([60,140])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_4state_fast.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_4state_fast.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_4state_fast.pdf');
%export_fig post_dens_4state_fast.png -transparent

%load sims_5state_slow.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_5state_slow.cum_mode],'Normalization','probability');
hist([sims_5state_slow.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([60,140])
hold off 
xlabel({'M'},'FontSize',18);
ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_5state_slow.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_5state_slow.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_5state_slow.pdf');
%export_fig post_dens_5state_slow.png -transparent

%load sims_5state_med.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_5state_med.cum_mode],'Normalization','probability');
hist([sims_5state_med.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([70,130])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_5state_med.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_5state_med.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_5state_med.pdf');
%export_fig post_dens_5state_med.png -transparent

%load sims_5state_fast.mat 

figure 
hold on 
box on 
set(gca,'FontSize',18);
%h = histogram([sims_5state_fast.cum_mode],'Normalization','probability');
hist([sims_5state_fast.cum_mode])
%title('Posterior of K');
xline(100,'LineWidth',3,'Color','red')
xticks([80 100 120]); xlim([70,130])
hold off 
xlabel({'M'},'FontSize',18);
%ylabel({'Counts'},'FontSize',18);

filename = 'post_dens_5state_fast.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
%export_fig post_dens_5state_fast.pdf -transparent
%export_fig(filename, '-eps', '-opengl')
%print(gcf, '-dpdf', '-bestfit');
%saveas(gcf,'post_dens_5state_fast.pdf');
%export_fig post_dens_5state_fast.png -transparent