load testing_number_of_locs_alexa_2.mat
load testing_number_of_locs_alexa_3state.mat
load testing_number_of_locs_alexa_4state.mat
load testing_number_of_locs_alexa_negbin.mat 

M_5state = [testing_number_of_locs_alexa_2.cum_mode];
M_4state = [testing_number_of_locs_alexa_4state.cum_mode];
M_3state = [testing_number_of_locs_alexa_3state.cum_mode];
M_negbin = [testing_number_of_locs_alexa_negbin.cum_mode];
y = 4:6:81*2;
% 
% yneg = [1 3 5 3 5 3 6 4 3 3];
% ypos = [2 5 3 5 2 5 2 2 5 5];
xneg_5state = 1:27; 
xpos_5state = 1:27; 
xneg_4state = 1:27; 
xpos_4state = 1:27; 
xneg_3state = 1:27; 
xpos_3state = 1:27; 

for i=1:27 
    ci_5state = testing_number_of_locs_alexa_2(i).cum_CI; 
    ci_4state = testing_number_of_locs_alexa_4state(i).cum_CI; 
    ci_3state = testing_number_of_locs_alexa_3state(i).cum_CI; 
    mval_5state = M_5state(i); 
    mval_4state = M_4state(i); 
    mval_3state = M_3state(i); 
    xneg_5state(i) = mval_5state-ci_5state(1); 
    xpos_5state(i) = ci_5state(2)-mval_5state; 
    xneg_4state(i) = mval_4state-ci_4state(1); 
    xpos_4state(i) = ci_4state(2)-mval_4state; 
    xneg_3state(i) = mval_3state-ci_3state(1); 
    xpos_3state(i) = ci_3state(2)-mval_3state; 
end 

xneg_negbin = 1:27; 
xpos_negbin = 1:27; 

for i=1:27 
    ci_5state = testing_number_of_locs_alexa_negbin(i).cum_CI; 
    mval_5state = M_negbin(i); 
    xneg_negbin(i) = mval_5state-ci_5state(1); 
    xpos_negbin(i) = ci_5state(2)-mval_5state; 
end 

figure 
hold on 
box on 
set(gca,'FontSize',18);
xlabel({'M'},'FontSize',18);
xlim([40, 240])

yyaxis right 
errorbar(M_5state,2:6:162,xneg_5state,xpos_5state,'horizontal','x',...
    'LineWidth',2,'MarkerSize',8,'Color','black')
d1 = errorbar(M_4state,4:6:162,xneg_4state,xpos_4state,'horizontal',...
    'x','LineWidth',2,'MarkerSize',8,'Color','blue');
d1.Bar.LineStyle = 'dotted'; 
d2 = errorbar(M_3state,6:6:162,xneg_3state,xpos_3state,'horizontal',...
    'x','LineWidth',2,'MarkerSize',8,'Color',[0.4 0.4 0.4]);
d2.Bar.LineStyle = 'dashdot'; 
% errorbar(M_negbin,y,xneg_negbin,xpos_negbin,'horizontal','o','LineWidth',2,'MarkerSize',8,'Color','magenta')

xs = [testing_number_of_locs_alexa_2.N_true];
ys = 4:6:162;
errorbar(xs,ys,2.5*ones(size(ys)),'vertical','.','LineWidth',3,...
    'MarkerSize',8,'Color',[.85 0 0])

% val = 1; 
% for i=1:27
%     for j=1:3
%         scatter(testing_number_of_locs_alexa_2(i).N_true,j+val-1,70,'xr','LineWidth',2); 
%     end 
%     val = j+val;
% end

% xticks([60 120 180 240 300 360 420])
xticks([60 90 120 150 180 210 240])

plot([40, 440], 2*[3.5, 3.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[9.5, 9.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[12.5, 12.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[24.5, 24.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[36.5, 36.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[51.5, 51.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[66.5, 66.5],'--','Color','black','LineWidth',1.5)
ylim([0, 2*3*27+1])
yticks((6:12:162)-1)
yticklabels((6:12:162)/6)
ylabel({'Dataset'},'FontSize',18);

yyaxis left
% yticks([2 4 5.5 7.5 10.5 15 20 25]./30)
yticks(2*[2 6.5 11.5 18.5 30.5 44 59 74]./165)
yticklabels({'1.0','1.9','3.9','7.8','16','31','62','97'})
ylabel({'Laser Intensity (kW/cm^2)'},'FontSize',18);

% export_fig alexa_all_models_2.pdf -transparent
filename = 'alexa_all_models_2.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18,...
    'height', 10, 'width',10,'Color','rgb')


figure 
hold on 
box on 
set(gca,'FontSize',18);
xlabel({'M'},'FontSize',18);
xlim([40, 440])

yyaxis right 
errorbar(M_5state,2:6:162,xneg_5state,xpos_5state,'horizontal','x',...
    'LineWidth',2,'MarkerSize',8,'Color','black')
d1 = errorbar(M_4state,4:6:162,xneg_4state,xpos_4state,'horizontal',...
    'x','LineWidth',2,'MarkerSize',8,'Color','blue');
d1.Bar.LineStyle = 'dotted'; 
d2 = errorbar(M_3state,6:6:162,xneg_3state,xpos_3state,'horizontal',...
    'x','LineWidth',2,'MarkerSize',8,'Color',[0.4 0.4 0.4]);
d2.Bar.LineStyle = 'dashdot'; 
errorbar(M_negbin,y,xneg_negbin,xpos_negbin,'horizontal','x',...
    'LineWidth',2,'MarkerSize',8,'Color','magenta')

xs = [testing_number_of_locs_alexa_2.N_true];
ys = 4:6:162;
errorbar(xs,ys,2.5*ones(size(ys)),'vertical','.','LineWidth',2,...
    'MarkerSize',8,'Color',[.85 0 0])

% val = 1; 
% for i=1:27
%     for j=1:3
%         scatter(testing_number_of_locs_alexa_2(i).N_true,j+val-1,70,'xr','LineWidth',2); 
%     end 
%     val = j+val;
% end

xticks([60 120 180 240 300 360 420])
% xticks([60 90 120 150 180 210 240])

plot([40, 440], 2*[3.5, 3.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[9.5, 9.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[12.5, 12.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[24.5, 24.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[36.5, 36.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[51.5, 51.5],'--','Color','black','LineWidth',1.5)
plot([40, 440], 2*[66.5, 66.5],'--','Color','black','LineWidth',1.5)
ylim([0, 2*3*27+1])
yticks((6:12:162)-1)
yticklabels((6:12:162)/6)
ylabel({'Dataset'},'FontSize',18);

yyaxis left
% yticks([2 4 5.5 7.5 10.5 15 20 25]./30)
yticks(2*[2 6.5 11.5 18.5 30.5 44 59 74]./165)
yticklabels({'1.0','1.9','3.9','7.8','16','31','62','97'})
ylabel({'Laser Intensity (kW/cm^2)'},'FontSize',18);

export_fig alexa_all_models_all_2.pdf -transparent
filename = 'alexa_all_models_all_2.eps'; 
exportfig(gcf, filename, 'FontMode', 'fixed','FontSize', 18, ...
    'height', 10, 'width',10,'Color','rgb')

params_4state = zeros(27,10); 
for i=1:27 
    params = testing_number_of_locs_alexa_4state(i).raw_params * ...
        testing_number_of_locs_alexa_4state(i).Delta;
    params_4state(i,1:4) = params(1:4) .* [10 1 1e4 1]; 
    params_4state(i,5) = params(7) * 100; 
    params_4state(i,6) = params(8)/...
        (testing_number_of_locs_alexa_4state(i).Delta)^2; 
    params_4state(i,7) = params(9) * 1e5 / ...
        (testing_number_of_locs_alexa_4state(i).Delta); 
    nu = testing_number_of_locs_alexa_4state(i).nu_X; 
    params_4state(i,8:10) = nu(1:end-1); 
end 

params_3state = zeros(27,7); 
for i=1:27 
    params = testing_number_of_locs_alexa_3state(i).raw_params * ...
        testing_number_of_locs_alexa_3state(i).Delta;
    params_3state(i,1:2) = params(1:2) .* [1e3 1]; 
    params_3state(i,3) = params(4) * 10; 
    params_3state(i,4) = params(5)/ ...
        (testing_number_of_locs_alexa_3state(i).Delta)^2; 
    params_3state(i,5) = params(6) * 1e5/ ...
        (testing_number_of_locs_alexa_4state(i).Delta); 
    nu = testing_number_of_locs_alexa_3state(i).nu_X; 
    params_3state(i,6:7) = nu(1:end-1); 
end 


  