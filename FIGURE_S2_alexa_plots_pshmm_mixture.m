load testing_number_of_locs_alexa_2.mat

Dataset = 'Dataset'; 
file_name = 'post_k_alexa_';

for i=1:27 
    s = num2str(i); 
    p_no_locs_blinks_alexa = testing_number_of_locs_alexa_2(i).p_no_locs_blinks_alexa;
    N_val = max(p_no_locs_blinks_alexa); 
    ks = testing_number_of_locs_alexa_2(i).cum_CI; 
    N_mode = testing_number_of_locs_alexa_2(i).cum_mode;
    N_true = testing_number_of_locs_alexa_2(i).N_true; 
    N_mode_nb = testing_number_of_locs_alexa_negbin(i).cum_mode; 
    
    figure 
    box on 
    hold on
    set(gca,'FontSize',18);
    stairs(p_no_locs_blinks_alexa,'LineWidth',2,'color','black');
    line([N_mode N_mode], [0 N_val+1e-3],'color','cyan','LineWidth',2);
    line([N_mode_nb N_mode_nb], [0 N_val+1e-3],'color','red','LineWidth',2);
    line([N_true N_true], [0 N_val+1e-3],'color','magenta','LineStyle','--','LineWidth',2);
    line([ks(1) ks(1)], [0 N_val+1e-3],'LineWidth',2,'LineStyle','--','color','black');
    line([ks(2) ks(2)], [0 N_val+1e-3],'LineWidth',2,'LineStyle','--','color','black');
    xlim([ks(1)-20, ks(2)+20])
    ylim([0 N_val+1e-3])
    title(strcat({'Dataset '},s));
    xlabel({'M'},'FontSize',18);
    mid = roundn(N_mode,1); 
    xlim([min(N_mode_nb-20,N_mode-20),max(N_mode_nb+20,N_mode+20)])
    xticks([mid-60 mid-30 mid mid+30 mid+60])
    if mod(i,3) == 1
        ylabel({'Posterior probability'},'FontSize',18);
    end
    hold off 
    % title('Dataset 4');
    % xline(N_mode,'LineWidth',3,'Color','red')
    % xline(N_true,'LineWidth',3,'Color','green')
    % xlim([N_mode-40,N_mode+40])
    % hold off 
    % xlabel({'K'},'FontSize',18);
    % ylabel({'Probability'},'FontSize',18);
    FIGURE_NAME = strcat(file_name,s,'.eps'); 
%     
%     if ischar(FIGURE_NAME)
%         filename = sprintf('%s', FIGURE_NAME);
%     else
%         filename = sprintf('%d', FIGURE_NAME);
%     end
%     FIGURE_NAME = 'post_m_alexa_1.eps';
    exportfig(gcf, FIGURE_NAME, 'FontMode', 'fixed','FontSize', 18,'height', 10, 'width',10,'Color','rgb')
    %saveas(gcf,FIGURE_NAME);
    %export_fig(filename, '-pdf', '-transparent')
end     
