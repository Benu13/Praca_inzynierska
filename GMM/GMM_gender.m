%% Gaussian Mixture Model 
% gaussian_form is used for estimating normal distribution applayed to
% vector.
%% DATA FOR ALGORITHM
% Input data in form of tables, specified somewhere idk
clear all; close all;
train_table = readtable('CCOF_TABLE_TRAIN_AGE32m.txt');
test_table = readtable('CCOF_TABLE_TEST_AGE32m.txt');
% fmm = readtable('Fem_table.txt');
% mff = readtable('Men_table.txt');

% train_table = [fmm(1:90,:);mff(1:90,:)];
% test_table = [fmm(91:end,:);mff(91:end,:)];
answ_num = 2; % number of answers we can predict eg. gender = 2 (male/female); not implemented rn;
max_iter = 400; %max iterations
gauss_models = 5; % Number of gaussians mixtures for model
cepst_num = 32; %16: 16-3,9,16 %322: 9,7,16

% iidx = 1;
% for cepst_num = 20
%     gauss_models = 19;
%     while gauss_models <= cepst_num
%% Calculate data for modeling Gaussian Mixture models
% We need to separate male and female samples to model both female and male
% gaussian mixtures. It best when there is same ammount of male and female
% samples, so when samples for one sex exceeds samples for other we need to
% take the same amount.
fem_samples_num = sum(cell2mat(strfind(train_table.Sex,'F')));
male_samples_num = size(train_table.Sex,1)-fem_samples_num;
samp_num = min([male_samples_num fem_samples_num]);
t_data{1} = table2array(train_table(1:samp_num,2:cepst_num+1));
t_data{2} = table2array(train_table(male_samples_num+1:male_samples_num+samp_num, 2:cepst_num+1));
% t_data{1} = table2array(train_table(1:male_samples_num,2:cepst_num+1));
% t_data{2} = table2array(train_table(male_samples_num+1:fem_samples_num, 2:cepst_num+1));
labels = train_table(:,1);
%% Initialization of parameters
for model_num = 1:answ_num
    x = cell2mat(t_data(model_num)); % Read data matrice for current model
    %% Get starting mean and covariance form k-means
    [pi,muu,sigma] = kmean_param(x,gauss_models);
    %% Fit data
    loss = 0;
    best_loss = -Inf;
    for i = 1:max_iter
        gamma = e_step(x,gauss_models,pi,muu,sigma);
        [pi,muu,sigma] = m_step(x,gamma);
        loss_old = loss;
        loss = cl_fun(x,pi,muu,sigma,gamma);
        lo(model_num, i) = loss;
        if mod(i,10) == 0
            sprintf(['Iteration: ' num2str(i) ' Loss: ' num2str(loss)])
        end
        if loss == loss_old
            break
        end
        if loss > best_loss
            best_pi = pi; best_muu = muu; best_sigma = sigma; best_gamma = gamma;
            best_loss = loss;
        end
    end
    %Save created model
    model(model_num).pi = best_pi; model(model_num).muu = best_muu;
    model(model_num).sigma = best_sigma;  model(model_num).gamma = best_gamma;
end
%% Testing created models
test_mat = table2array(test_table(:,2:cepst_num+1));
test_len = size(test_mat,1); male_ind = sum(cell2mat(strfind(test_table.Sex,'M')));
true_labels = [ones(male_ind,1); 2*ones(test_len-male_ind,1)];

predicted_labes = zeros(size(test_mat,1),1);
for i =1:size(test_mat,1)
    logL_m = cl_fun(test_mat(i,:),model(1).pi,model(1).muu,model(1).sigma,model(1).gamma);
    logL_f = cl_fun(test_mat(i,:),model(2).pi,model(2).muu,model(2).sigma,model(2).gamma);
    if logL_m > logL_f
        predicted_labes(i) = 1;
    else
        predicted_labes(i) = 2;
    end
end
%% confusion matrix and accuracy for windows of signal
conf_mat_windows = confusionmat(true_labels,predicted_labes)
accuracy_windows = trace(conf_mat_windows)/sum(conf_mat_windows,'All')
% anww(iidx,:) = [cepst_num gauss_models accuracy_windows];
% iidx = iidx + 1;
% gauss_models = gauss_models+1;
%     end
% end
%% confusion matrix and accuracy for speaker probe
num_male_spk = max(test_table.ManNum(1:male_ind));
num_female_spk = max(test_table.ManNum(male_ind+1:end));
true_gender = [ones(num_male_spk,1); 2*ones(num_female_spk,1)];

pred_man_g = zeros(num_male_spk,1);
pred_fem_g = zeros(num_female_spk,1);

man_test_answ= test_table.ManNum(1:male_ind);
female_test_answ = test_table.ManNum(male_ind+1:end);
%%
for i = 1:num_male_spk
    speaker_probes = find(man_test_answ == i);
    if sum(ismember(predicted_labes(speaker_probes),1)) > fix(length(speaker_probes)/2)
        pred_man_g(i) = 1;
    else
        pred_man_g(i) = 2;
    end
end
for i = 1:num_female_spk
    speaker_probes = find(female_test_answ == i)+male_ind;
    if sum(ismember(predicted_labes(speaker_probes),2)) > fix(length(speaker_probes)/2)
        pred_fem_g(i) = 2;
    else
        pred_fem_g(i) = 1;
    end
end
pred_sp_gender = [pred_man_g;pred_fem_g];

conf_mat_speaker = confusionmat(true_gender,pred_sp_gender)
accuracy_s = trace(conf_mat_speaker)/sum(conf_mat_speaker,'All')
%%
% subplot(1,2,1); 
% confusionchart(conf_mat_windows)
% title(['Dokładność = ' num2str(round(accuracy_windows*100,2,'significant')) '%'])
% subplot(1,2,2);
% confusionchart(conf_mat_speaker)
% title(['Dokładność = ' num2str(round(accuracy_s*100,2,'significant')) '%'])
