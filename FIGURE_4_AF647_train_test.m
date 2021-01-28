mex -O analytic_prs_3.c
mex -O fwdbkwd_norm_c2.c
filename_start = 'AF647_2pc405_';
filename_end = '_by_frame.csv';

Data_1 = zeros(50000,22); 
for i=1:22
    s = num2str(i);
    Data_1(:,i) = csvread(strcat(filename_start,s,filename_end));
end 

Y_1 = Data_1'; 

cd ..

Delta = 0.02; %Frame acqusition time in seconds for training and testing
% [nu_X_all,lambda_all,mu_all, delta_all, alpha_all, B0_all,B1_all] = PSHMM_model_select(Y_3,Delta);
% [nu_X_no,lambda_no,mu_no, delta_no, alpha_no, B0_no,B1_no] = PSHMM_model_select(Y_2,Delta);
[nu_X,lambda,mu, delta, alpha, B0,B1] = PSHMM_model_select(Y_1,Delta);

AF647_est_params = {}; 
AF647_est_params.nu_X = nu_X; 
AF647_est_params.lambda = lambda;
AF647_est_params.mu = mu;
AF647_est_params.delta = delta;
AF647_est_params.FPR = alpha;
AF647_est_params.B0 = B0;
AF647_est_params.B1 = B1;

save AF647_est_params.mat AF647_est_params 

% filename_start = 'AF647DriCorMrgFil_'; %With no drift correction
filename_start = 'AF647DriCor_'; %With drift correction
filename_end = '.csv';
data_files = cell(4); 
NFs = 1:4; 
for i=1:4
    s = num2str(i);
    data_files{i} = csvread(strcat(filename_start,s,filename_end),1,0);
    NFs(i) = max(data_files{i}(:,2)); 
end 

cd .. 
%First file contains fewer frames than files 2:4 
[probs_1, expected_1, variance_1] = no_locs_pmf(nu_X,B0,B1,NFs(1));
[probs_2, expected_2, variance_2] = no_locs_pmf(nu_X,B0,B1,40000); %N_F = 40000 

figure 
hold on 
box on 
set(gca,'FontSize',18);
title('AF647 distribution of localizations');
xlabel({'Localizations'},'FontSize',18);
ylabel({'Probability'},'FontSize',18);
stairs(0:1000,probs_2(1:1001),'LineWidth',2,'Color','blue')
% plot(0:1000,exp_pdf,'LineWidth',2,'Color','green')
xline(expected_2,'--r','LineWidth',2)
xline(expected_2+sqrt(variance_2),'--','LineWidth',2)
legend('Probability mass function',strcat('Mean =',num2str(round(expected_2,2))),...
    strcat('Standard deviation =',num2str(round(sqrt(variance_2),2))))
hold off 
%saveas(gcf,'asilomar_ks_4.pdf');
export_fig dist_of_locs_af647.pdf -transparent

probs = cell(4); 
expected = 1:4; 
variance = 1:4; 
probs{1} = probs_1; 
expected(1) = expected_1; 
variance(1) = variance_1; 
for i=2:4 
    probs{i} = probs_2; 
    expected(i) = expected_2; 
    variance(i) = variance_2; 
end 

%Ok now use Dave's txt file
[cell_name, cell_number , file_number, x_cent, y_cent] = textread('coords.txt', ...
'%s %f %f %d %d');
Nl = cell_number; 
for i=1:length(Nl)
    file = file_number(i); 
    vec = (data_files{file}(:,3) > x_cent(i)-1500) & ...
        (data_files{file}(:,3) < x_cent(i)+1500) & ...
        (data_files{file}(:,4) > y_cent(i)-1500) & ...
        (data_files{file}(:,4) < y_cent(i)+1500); %find all localisations within the specified grid
    Nl(i) = sum(vec); 
end 

% p_no_locs = cell(length(Nl)); 
testing_number_of_locs_AF647 = struct; 

for i=1:length(Nl) 
    testing_number_of_locs_AF647(i).Nl = Nl(i); 
    testing_number_of_locs_AF647(i).cell_number = cell_number(i);
    testing_number_of_locs_AF647(i).file_number = file_number(i);
    testing_number_of_locs_AF647(i).x_cent = x_cent(i);
    testing_number_of_locs_AF647(i).y_cent = y_cent(i);
end 

for i=1:length(Nl) 
    fprintf("i is %d \n",i)
    file = file_number(i); 
    khat = ceil(Nl(i)/expected(file));
    kmin = max(ceil(Nl(i)/NFs(file)), floor(khat - 4*sqrt(variance(file))));
    kmax = ceil(khat + 3*sqrt(variance(file)));
    p_no_locs = zeros(1,kmax+1);
    prior = 1:kmax+1;
    
    parfor j=kmin:kmax
        vec_probs = abs(ifft(fft([probs{file} ...
            zeros(1,j*NFs(file)+1-length(probs{file}))]).^j));
        probv = vec_probs(Nl(i)+1);
        p_no_locs(j) = probv;
    end

    clear vec_probs
    p_no_locs_grid = p_no_locs/sum(p_no_locs);
    testing_number_of_locs_AF647(i).p_no_locs_grid = p_no_locs_grid;
    [N_val,N_mode] = max(p_no_locs_grid);
    v = cumsum(p_no_locs_grid) < 0.975 & cumsum(p_no_locs_grid) > 0.025; %95% of the variance
    ks = prior(diff(v)~=0);
    testing_number_of_locs_AF647(i).cum_mode = N_mode;
    testing_number_of_locs_AF647(i).cum_CI = ks;
end 

save('testing_number_of_locs_AF647.mat','testing_number_of_locs_AF647','-v7.3')
save -v7.3 testing_number_of_locs_AF647.mat testing_number_of_locs_AF647
