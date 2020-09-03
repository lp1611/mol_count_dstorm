load p_sum_of_locs_sims.mat
load p_sum_of_locs_true.mat

figure 
hold on 
box on 
set(gca,'FontSize',18);
stairs(p_sum_of_locs_sims(1,:),'LineWidth',1.5,'Color','black')
stairs(p_sum_of_locs_true(1,:),'LineWidth',3.5,'Color','black','LineStyle','--')
stairs(p_sum_of_locs_sims(2,:),'LineWidth',1.5,'Color','blue')
stairs(p_sum_of_locs_true(2,:),'LineWidth',3.5,'Color','blue','LineStyle','--')
stairs(p_sum_of_locs_sims(3,:),'LineWidth',1.5,'Color','cyan')
stairs(p_sum_of_locs_true(3,:),'LineWidth',3.5,'Color','cyan','LineStyle','--')
stairs(p_sum_of_locs_sims(4,:),'LineWidth',1.5,'Color','magenta')
stairs(p_sum_of_locs_true(4,:),'LineWidth',3.5,'Color','magenta','LineStyle','--')
stairs(p_sum_of_locs_sims(5,:),'LineWidth',1.5,'Color','red')
stairs(p_sum_of_locs_true(5,:),'LineWidth',3.5,'Color','red','LineStyle','--')
xlim([0, 500])
xticks([0,250,500])
hold off 
xlabel({'N_F'},'FontSize',18);
ylabel({'Probability'},'FontSize',18);
savefig('sum_locs_mus.fig')
filename = 'sum_locs_mus.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')

%%%%% 
figure 
hold on 
box on
set(gca,'FontSize',18);
%h = histogram([sims_4state_med.cum_mode],'Normalization','probability');
stairs(p_sum_of_locs_sims(6,:),'LineWidth',1.5,'Color','black')
stairs(p_sum_of_locs_true(6,:),'LineWidth',3.5,'Color','black','LineStyle','--') %true
stairs(p_sum_of_locs_sims(2,:),'LineWidth',1.5,'Color','blue','LineStyle','--')
stairs(p_sum_of_locs_true(2,:),'LineWidth',3.5,'Color','blue','LineStyle','--') %true 
stairs(p_sum_of_locs_sims(7,:),'LineWidth',1.5,'Color','red')
stairs(p_sum_of_locs_true(7,:),'LineWidth',3.5,'Color','red','LineStyle','--') %true
xlim([0, 600])
xticks([0,250,500])
hold off 
xlabel({'N_F'},'FontSize',18);
ylabel({'Probability (x 0.001)'},'FontSize',18);
filename = 'sum_locs_d.eps'; 
savefig('sum_locs_d.fig')
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')

figure 
hold on 
box on 
set(gca,'FontSize',18);
xlabel({'Number of localisations'},'FontSize',18);
ylabel({'Probability (\times 10^{-3})'},'FontSize',18);
xlim([0 600])
xticks([0 200 400 600])
ylim([0 0.005])
kk=2; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','red')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','red')

kk=7; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','green')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','green')

kk=6; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','blue')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','blue')
filename = 'cum_sims_1.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 20, 'width',20,'Color','rgb')
hold off 


figure 
hold on 
box on 
set(gca,'FontSize',18);
xlabel({'Number of localisations'},'FontSize',18);
ylabel({'Probability'},'FontSize',18);
xlim([0 600])
xticks([0 200 400 600])
ylim([0 0.0155])
kk=6; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','green')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','green')

kk=3; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','red')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','red')

kk=4; 
hist_counts = []; 
for ii=0:1e3
    if p_sum_of_locs_sims(1,ii+1) ~= 0 
        hist_counts = [hist_counts ii*ones(1,round(p_sum_of_locs_sims(kk,ii+1)*1e6))];
    end 
end 
histogram(hist_counts,'Normalization','pdf','FaceAlpha',0.1,'FaceColor','blue')
stairs(p_sum_of_locs_true(kk,:),'LineWidth',2.5,'Color','blue')
filename = 'cum_sims_2.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 20, 'width',20,'Color','rgb')
hold off 


