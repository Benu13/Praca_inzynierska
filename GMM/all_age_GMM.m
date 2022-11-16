%% Gaussian Mixture Model 
% gaussian_form is used for estimating normal distribution applayed to
% vector.
%% DATA FOR ALGORITHM
% Input data in form of tables, specified somewhere idk
clear all; close all;
train_table = readtable('all_TABLE_TRAIN_AGE.txt');
test_table = readtable('all_TABLE_TEST_AGE.txt');
answ_num = 4; % number of answers we can predict eg. gender = 2 (male/female); not implemented rn;
max_iter = 400; %max iterations
gauss_models =11; % Number of gaussians mixtures for model
cepst_num =20; %% 4/15 ftw, 17/5
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
teens_samples_num = length((find(train_table.Age == 13)));
twenties_samples_num = length((find(train_table.Age == 20)));
thirties_samples_num = length((find(train_table.Age == 30)));
fourties_samples_num = length((find(train_table.Age == 40)));

minnum = min([teens_samples_num, twenties_samples_num,thirties_samples_num,fourties_samples_num]);
t_data{1} = table2array(train_table(1:minnum,2:cepst_num+1));
t_data{2} = table2array(train_table(teens_samples_num+1:teens_samples_num+minnum,2:cepst_num+1));
t_data{3} = table2array(train_table(teens_samples_num+twenties_samples_num+1:teens_samples_num +twenties_samples_num+minnum, 2:cepst_num+1));
t_data{4} = table2array(train_table(teens_samples_num+twenties_samples_num+thirties_samples_num+1:teens_samples_num+twenties_samples_num+thirties_samples_num+minnum, 2:cepst_num+1));

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
teens_stl = length((find(test_table.Age == 13)));
twenties_stl = length((find(test_table.Age == 20)));
thirties_stl= length((find(test_table.Age == 30)));
fourties_stl= length((find(test_table.Age == 40)));

test_len = size(test_mat,1);
true_labels = [ones(teens_stl,1); 2*ones(twenties_stl,1); 3*ones(thirties_stl,1); 4*ones(fourties_stl,1)];

predicted_labes = zeros(size(test_mat,1),1);
for i =1:size(test_mat,1)
    logL_10 = cl_fun(test_mat(i,:),model(1).pi,model(1).muu,model(1).sigma,model(1).gamma);
    logL_20 = cl_fun(test_mat(i,:),model(2).pi,model(2).muu,model(2).sigma,model(2).gamma);
    logL_30 = cl_fun(test_mat(i,:),model(3).pi,model(2).muu,model(3).sigma,model(3).gamma);
    logL_40 = cl_fun(test_mat(i,:),model(4).pi,model(4).muu,model(4).sigma,model(4).gamma);
    answ = [logL_10 logL_20 logL_30 logL_40];
    [~, predicted_labes(i)] = max(answ);
end
%% confusion matrix and accuracy for windows of signal
conf_mat_windows = confusionmat(true_labels,predicted_labes)
accuracy_windows = trace(conf_mat_windows)/sum(conf_mat_windows,'All')

%% confusion matrix and accuracy for speaker probe
test_table.ManNum = test_table.ManNum +20;

num_te_spk = max(test_table.ManNum(1:teens_stl));
num_tw_spk = max(test_table.ManNum(teens_stl+1:teens_stl+twenties_stl));
num_th_spk = max(test_table.ManNum(teens_stl+twenties_stl+1:teens_stl+twenties_stl+thirties_stl));
num_fr_spk = max(test_table.ManNum(teens_stl+twenties_stl+thirties_stl+1:end));

true_age = [ones(num_tw_spk,1); 2*ones(num_tw_spk,1); 3*ones(num_th_spk,1); 4*ones(num_fr_spk,1) ];

pred_10_g = zeros(num_te_spk,1);
pred_20_g = zeros(num_tw_spk,1);
pred_30_g = zeros(num_th_spk,1);
pred_40_g = zeros(num_fr_spk,1);

te_test_answ= test_table.ManNum(1:teens_stl);
tw_test_answ= test_table.ManNum(teens_stl+1:teens_stl+twenties_stl);
th_test_answ = test_table.ManNum(teens_stl+twenties_stl+1:teens_stl+twenties_stl+thirties_stl);
fr_test_answ = test_table.ManNum(teens_stl+twenties_stl+thirties_stl+1:end);
%%
for i = 1:num_tw_spk
    speaker_probes = find(te_test_answ == i);
    pred_10_g(i) = mode(predicted_labes(speaker_probes));
end

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

pred_sp_age = [pred_10_g; pred_20_g;pred_30_g; pred_40_g];

conf_mat_speaker = confusionmat(true_age,pred_sp_age)
accuracy_s = trace(conf_mat_speaker)/sum(conf_mat_speaker,'All')
%%
figure()
subplot(1,2,1)
plotConfMat(conf_mat_windows',['13-19';'20-29';'30-39'])
subplot(1,2,2)
plotConfMat(conf_mat_speaker',['13-19'; '20-29';'30-39'])
%%
function plotConfMat(varargin)
%PLOTCONFMAT plots the confusion matrix with colorscale, absolute numbers
%   and precision normalized percentages
%
%   usage: 
%   PLOTCONFMAT(confmat) plots the confmat with integers 1 to n as class labels
%   PLOTCONFMAT(confmat, labels) plots the confmat with the specified labels
%
%   Vahe Tshitoyan
%   20/08/2017
%
%   Arguments
%   confmat:            a square confusion matrix
%   labels (optional):  vector of class labels
% number of arguments
switch (nargin)
    case 0
       confmat = 1;
       labels = {'1'};
    case 1
       confmat = varargin{1};
       labels = 1:size(confmat, 1);
    otherwise
       confmat = varargin{1};
       labels = varargin{2};
end
confmat(isnan(confmat))=0; % in case there are NaN elements
numlabels = size(confmat, 1); % number of labels
% calculate the percentage accuracies
confpercent = 100*confmat./repmat(sum(confmat, 1),numlabels,1);
% plotting the colors
imagesc(confpercent);
title(sprintf('Dokładność: %.2f%%', 100*trace(confmat)/sum(confmat(:))));
ylabel('Wynik klasyfikacji'); xlabel('Prawdziwa wartość');
% set the colormap
colormap(flipud(gray));
% Create strings from the matrix values and remove spaces
textStrings = num2str([confpercent(:), confmat(:)], '%.1f%%\n%d\n');
textStrings = strtrim(cellstr(textStrings));
% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:numlabels);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center');
% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));
% Choose white or black for the text color of the strings so
% they can be easily seen over the background color
textColors = repmat(confpercent(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'TickLength',[0 0]);
end