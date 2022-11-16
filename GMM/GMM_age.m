%% Gaussian Mixture Model 
% gaussian_form is used for estimating normal distribution applayed to
% vector.
%% DATA FOR ALGORITHM
% Input data in form of tables, specified somewhere idk
clear all; close all;
train_table = readtable('CCOF_TABLE_TRAIN_AGE_32.txt');
test_table = readtable('CCOF_TABLE_TEST_AGE_32.txt');
answ_num = 3; % number of answers we can predict eg. gender = 2 (male/female); not implemented rn;
max_iter = 400; %max iterations
gauss_models =5; % Number of gaussians mixtures for model
cepst_num =32; %% 4/15 ftw, 17/5 32/5 ftw
%% Calculate data for modeling Gaussian Mixture models
% We need to separate male and female samples to model both female and male
% gaussian mixtures. It best when there is same ammount of male and female
% samples, so when samples for one sex exceeds samples for other we need to
% take the same amount.
% train_table(isnan(train_table.Delta_1),:) = [];
% test_table(isnan(test_table.Delta_1),:) = [];
% 
% train_table(train_table.Delta_1 == -Inf,:) = [];
% test_table(train_table.Delta_1 == -Inf,:) = [];
% 
% train_table(train_table.Delta_1 == Inf,:) = [];
% test_table(train_table.Delta_1 == Inf,:) = [];
% 
% train_table(train_table.Delta_1 == 0,:) = [];
% test_table(train_table.Delta_1 == 0,:) = [];

twenties_samples_num = length((find(train_table.Age == 20)));
thirties_samples_num = length((find(train_table.Age == 30)));
fourties_samples_num = length((find(train_table.Age == 40)));

minnum = min([twenties_samples_num,thirties_samples_num,fourties_samples_num]);

% t_data{1} = table2array(train_table(1:twenties_samples_num,2:cepst_num+1));
% t_data{2} = table2array(train_table(twenties_samples_num+1:twenties_samples_num+thirties_samples_num , 2:cepst_num+1));
% t_data{3} = table2array(train_table(twenties_samples_num+thirties_samples_num+1:twenties_samples_num+thirties_samples_num+fourties_samples_num, 2:cepst_num+1));


t_data{1} = table2array(train_table(1:minnum,2:cepst_num+1));
t_data{2} = table2array(train_table(twenties_samples_num+1:twenties_samples_num+minnum, 2:cepst_num+1));
t_data{3} = table2array(train_table(twenties_samples_num+thirties_samples_num+1:twenties_samples_num+thirties_samples_num+minnum, 2:cepst_num+1));

labels = train_table(:,1);
%% Initialization of parameters
for model_num = 1:answ_num
    x = cell2mat(t_data(model_num)); % Read data matrice for current model
    %% Get starting mean and covariance form k-means
    [pi,muu,sigma] = kmean_param(x,gauss_models);
    %% Fit data
    loss = 0;
    gamma = [];
    best_loss = -Inf;
    for i = 1:max_iter
        gamma = e_step(x,gauss_models,pi,muu,sigma);
        [pi,muu,sigma] = m_step(x,gamma);
        loss_old = loss;
        lo(model_num,i)=loss;
        loss = cl_fun(x,pi,muu,sigma,gamma);
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
    model(model_num).pi = pi; model(model_num).muu = muu;
    model(model_num).sigma = sigma;  model(model_num).gamma = gamma;
end
%% Testing created models

test_mat = table2array(test_table(:,2:cepst_num+1));
twenties_stl = length((find(test_table.Age == 20)));
thirties_stl= length((find(test_table.Age == 30)));
fourties_stl= length((find(test_table.Age == 40)));

test_len = size(test_mat,1);
true_labels = [ones(twenties_stl,1); 2*ones(thirties_stl,1); 3*ones(fourties_stl,1)];

predicted_labes = zeros(size(test_mat,1),1);
for i =1:size(test_mat,1)
    logL_20 = cl_fun(test_mat(i,:),model(1).pi,model(1).muu,model(1).sigma,model(1).gamma);
    logL_30 = cl_fun(test_mat(i,:),model(2).pi,model(2).muu,model(2).sigma,model(2).gamma);
    logL_40 = cl_fun(test_mat(i,:),model(3).pi,model(3).muu,model(3).sigma,model(3).gamma);
    answ = [logL_20 logL_30 logL_40];
    [~, predicted_labes(i)] = max(answ);
end
%% confusion matrix and accuracy for windows of signal
conf_mat_windows = confusionmat(true_labels,predicted_labes)
accuracy_windows = trace(conf_mat_windows)/sum(conf_mat_windows,'All')

%% confusion matrix and accuracy for speaker probe
test_table.ManNum = test_table.ManNum ;

num_tw_spk = max(test_table.ManNum(1:twenties_stl));
num_th_spk = max(test_table.ManNum(twenties_stl+1:twenties_stl+thirties_stl));
num_fr_spk = max(test_table.ManNum(twenties_stl+thirties_stl+1:end));

true_age = [ones(num_tw_spk,1); 2*ones(num_th_spk,1); 3*ones(num_fr_spk,1) ];

pred_20_g = zeros(num_tw_spk,1);
pred_30_g = zeros(num_th_spk,1);
pred_40_g = zeros(num_fr_spk,1);

tw_test_answ= test_table.ManNum(1:twenties_stl);
th_test_answ = test_table.ManNum(twenties_stl+1:twenties_stl+thirties_stl);
fr_test_answ = test_table.ManNum(twenties_stl+thirties_stl+1:end);
%%
for i = 1:num_tw_spk
    speaker_probes = find(tw_test_answ == i);
    pred_20_g(i) = mode(predicted_labes(speaker_probes));
end

for i = 1:num_th_spk
    speaker_probes = find(th_test_answ == i)+twenties_stl;
    pred_30_g(i) = mode(predicted_labes(speaker_probes));
end

for i = 1:num_fr_spk
    speaker_probes = find(fr_test_answ == i)+twenties_stl+thirties_stl;
    pred_40_g(i) = mode(predicted_labes(speaker_probes));
end

pred_sp_age = [pred_20_g;pred_30_g; pred_40_g];

conf_mat_speaker = confusionmat(true_age,pred_sp_age)
accuracy_s = trace(conf_mat_speaker)/sum(conf_mat_speaker,'All')
%%
subplot(1,2,1); 
confusionchart(conf_mat_windows)
title(['Dokładność = ' num2str(round(accuracy_windows*100,2,'significant')) '%'])
subplot(1,2,2);
confusionchart(conf_mat_speaker)
title(['Dokładność = ' num2str(round(accuracy_s*100,2,'significant')) '%'])
