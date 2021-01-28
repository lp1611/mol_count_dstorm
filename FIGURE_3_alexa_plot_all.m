load testing_number_of_locs_alexa_2.mat
load testing_number_of_locs_alexa_negbin.mat

M = [testing_number_of_locs_alexa_2.cum_mode];
M_negbin = [testing_number_of_locs_alexa_negbin.cum_mode];
y = 1:27;
% 
% yneg = [1 3 5 3 5 3 6 4 3 3];
% ypos = [2 5 3 5 2 5 2 2 5 5];
xneg = 1:27; 
xpos = 1:27; 

for i=1:27 
    ci = testing_number_of_locs_alexa_2(i).cum_CI; 
    mval = M(i); 
    xneg(i) = mval-ci(1); 
    xpos(i) = ci(2)-mval; 
end 

xneg_negbin = 1:27; 
xpos_negbin = 1:27; 

for i=1:27 
    ci = testing_number_of_locs_alexa_negbin(i).cum_CI; 
    mval = M_negbin(i); 
    xneg_negbin(i) = mval-ci(1); 
    xpos_negbin(i) = ci(2)-mval; 
end 

figure 
hold on 
box on 
set(gca,'FontSize',18);
xlabel({'M'},'FontSize',18);
xlim([40, 440])

yyaxis right 
errorbar(M,y,xneg,xpos,'horizontal','o','LineWidth',2,'MarkerSize',8,'Color','blue')
errorbar(M_negbin,y,xneg_negbin,xpos_negbin,'horizontal','o','LineWidth',2,'MarkerSize',8,'Color','magenta')
for i=1:27 
    scatter(testing_number_of_locs_alexa_2(i).N_true,i,70,'xr','LineWidth',2); 
end

xticks([60 120 180 240 300 360 420])

plot([40, 440], [1.5, 1.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [2.5, 2.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [4.5, 4.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [8.5, 8.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [12.5, 12.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [17.5, 17.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], [22.5, 22.5],'--','Color','black','LineWidth',1.5)
ylim([0, 28])
ylabel({'Dataset'},'FontSize',18);

yyaxis left
yticks([1 2 3.5 6.5 10.5 15 20 25]./30)
yticklabels({'1.0','1.9','3.9','7.8','16','31','62','97'})
ylabel({'Laser Intensity (kW/cm^2)'},'FontSize',18);

export_fig alexa_all_both.pdf -transparent
% filename = 'alexa_all.eps'; 
% exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
% 



  
