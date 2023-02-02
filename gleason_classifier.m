function gleason_classifer
addpath /Volumes/Hydra/Analyses/SD_Recurrence/roi_tile_output/

%% Load and partition data
load feature_tbl.mat
load feature_tbl2.mat
%%
% train/test = unbalanced data
% train2 = balanaced
% train_hln/test_hln = train2 split into high, low, normal (semi-balanced lol)
% train_ng/test_ng = train2 split into normal vs all other grades

training = []; % training subject labels
testing = []; % testing subject labels

% Z-score first
train = [];
for i = 1:numel(training)
    trainingIdx = find(tbl.Subject == training(i));
    train = cat(1,train,tbl(trainingIdx,[1,2,7:end]));
end
train2 = train;
for i = 1:40
    train2{:,i+2} = zscore(train{:,i+2});
end


test = [];
for i = 1:numel(testing)
    testingIdx = find(tbl.Subject == testing(i));
    test = cat(1,test,tbl(testingIdx,[1,2,7:end]));
end
test2 = test;
for i = 1:24
    ave = mean(train{:,i+2});
    stdev = std(train{:,i+2},'omitnan');
    
    for j = 1:height(test2)
        test2{j,i+2} = (test2{j,i+2} - ave)/stdev;
    end
end


% Z-score first
trainingIdx = find(tbl.Category == "Training");
train = tbl(trainingIdx,[2,7:end-2]);
train2 = train;
for i = 1:40
    train2{:,i+1} = zscore(train{:,i+1});
end

testingIdx = find(tbl.Category == "Testing");
test = tbl(testingIdx,[2,7:end-2]);

test2 = test;
for i = 1:24
    ave = mean(train{:,i+1});
    stdev = std(train{:,i+1},'omitnan');
    
    for j = 1:height(test2)
        test2{j,i+1} = (test2{j,i+1} - ave)/stdev;
    end
end

% Z-score whole mounts used for demonstrations
load feature_tbl_2235.mat
load feature_tbl_2181.mat
load feature_tbl_1163.mat

for i = 1:24
    ave = mean(train{:,i+2});
    stdev = std(train{:,i+2},'omitnan');
    
    for j = 1:height(feature_tbl_2235)
        feature_tbl_2235{j,i+7} = (feature_tbl_2235{j,i+7} - ave)/stdev;
    end
end

for i = 1:24
    ave = mean(train{:,i+2});
    stdev = std(train{:,i+2},'omitnan');
    
    for j = 1:height(feature_tbl_2181)
        feature_tbl_2181{j,i+7} = (feature_tbl_2181{j,i+7} - ave)/stdev;
    end
end

for i = 1:24
    ave = mean(train{:,i+2});
    stdev = std(train{:,i+2},'omitnan');
    
    for j = 1:height(feature_tbl_1163)
        feature_tbl_1163{j,i+5} = (feature_tbl_1163{j,i+5} - ave)/stdev;
    end
end

%% H/L/N
train_hln = train2;
loc = find(train_hln.Anot_Class == 'Atrophy' | train_hln.Anot_Class == 'HGPIN' ...
    | train_hln.Anot_Class == 'Seminal_Vesicles' | train_hln.Anot_Class == 'Tissue');
train_hln.Anot_Class(loc) = 'Normal';

loc = find(train_hln.Anot_Class == 'G3');
train_hln.Anot_Class(loc) = 'Low';

loc = find(train_hln.Anot_Class == 'G4FG' | train_hln.Anot_Class == 'G4CG' | train_hln.Anot_Class == 'G5');
train_hln.Anot_Class(loc) = 'High';
%
test_hln = test2;
loc = find(test_hln.Anot_Class == 'Atrophy' | test_hln.Anot_Class == 'HGPIN' ...
    | test_hln.Anot_Class == 'Seminal_Vesicles' | test_hln.Anot_Class == 'Tissue');
test_hln.Anot_Class(loc) = 'Normal';

loc = find(test_hln.Anot_Class == 'G3');
test_hln.Anot_Class(loc) = 'Low';

loc = find(test_hln.Anot_Class == 'G4FG' | test_hln.Anot_Class == 'G4CG' | test_hln.Anot_Class == 'G5');
test_hln.Anot_Class(loc) = 'High';
%
loc = find(feature_tbl_2235.Anot_Class == 'Atrophy' | feature_tbl_2235.Anot_Class == 'HGPIN' ...
    | feature_tbl_2235.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2235.Anot_Class == 'Tissue');
feature_tbl_2235.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_2235.Anot_Class == 'G3');
feature_tbl_2235.Anot_Class(loc) = 'Low';

loc = find(feature_tbl_2235.Anot_Class == 'G4FG' | feature_tbl_2235.Anot_Class == 'G4CG' | feature_tbl_2235.Anot_Class == 'G5');
feature_tbl_2235.Anot_Class(loc) = 'High';
%
loc = find(feature_tbl_2181.Anot_Class == 'Atrophy' | feature_tbl_2181.Anot_Class == 'HGPIN' ...
    | feature_tbl_2181.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2181.Anot_Class == 'Tissue');
feature_tbl_2181.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_2181.Anot_Class == 'G3');
feature_tbl_2181.Anot_Class(loc) = 'Low';

loc = find(feature_tbl_2181.Anot_Class == 'G4FG' | feature_tbl_2181.Anot_Class == 'G4CG' | feature_tbl_2181.Anot_Class == 'G5');
feature_tbl_2181.Anot_Class(loc) = 'High';
%
loc = find(feature_tbl_1163.Anot_Class == 'Atrophy' | feature_tbl_1163.Anot_Class == 'HGPIN' ...
    | feature_tbl_1163.Anot_Class == 'Seminal_Vesicles' | feature_tbl_1163.Anot_Class == 'Tissue');
feature_tbl_1163.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_1163.Anot_Class == 'G3');
feature_tbl_1163.Anot_Class(loc) = 'Low';

loc = find(feature_tbl_1163.Anot_Class == 'G4FG' | feature_tbl_1163.Anot_Class == 'G4CG' | feature_tbl_1163.Anot_Class == 'G5');
feature_tbl_1163.Anot_Class(loc) = 'High';

mdl_hln_nb = fitcnb(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','all');
mdl_hln_knn = fitcknn(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','all');
mdl_hln_tree = fitctree(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','all');
mdl_hln_ecoc = fitcecoc(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','auto'); %
mdl_hln_ensemble = fitcensemble(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','all');
mdl_hln_net = fitcnet(train_hln(:,2:25),train_hln.Anot_Class,'OptimizeHyperparameters','all');

%% N/G 
train_ng = train2;
loc = find(train_ng.Anot_Class == 'Atrophy' | train_ng.Anot_Class == 'HGPIN' ...
    | train_ng.Anot_Class == 'Seminal_Vesicles' | train_ng.Anot_Class == 'Tissue');
train_ng.Anot_Class(loc) = 'Normal';
%
test_ng = test2;
loc = find(test_ng.Anot_Class == 'Atrophy' | test_ng.Anot_Class == 'HGPIN' ...
    | test_ng.Anot_Class == 'Seminal_Vesicles' | test_ng.Anot_Class == 'Tissue');
test_ng.Anot_Class(loc) = 'Normal';
%
loc = find(feature_tbl_2235.Anot_Class == 'Atrophy' | feature_tbl_2235.Anot_Class == 'HGPIN' ...
    | feature_tbl_2235.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2235.Anot_Class == 'Tissue');
feature_tbl_2235.Anot_Class(loc) = 'Normal';
%
loc = find(feature_tbl_2181.Anot_Class == 'Atrophy' | feature_tbl_2181.Anot_Class == 'HGPIN' ...
    | feature_tbl_2181.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2181.Anot_Class == 'Tissue');
feature_tbl_2181.Anot_Class(loc) = 'Normal';
%
loc = find(feature_tbl_1163.Anot_Class == 'Atrophy' | feature_tbl_1163.Anot_Class == 'HGPIN' ...
    | feature_tbl_1163.Anot_Class == 'Seminal_Vesicles' | feature_tbl_1163.Anot_Class == 'Tissue');
feature_tbl_1163.Anot_Class(loc) = 'Normal';

mdl_ng_nb = fitcnb(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','all');
mdl_ng_knn = fitcknn(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','all');
mdl_ng_tree = fitctree(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','all');
mdl_ng_ecoc = fitcecoc(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','auto');
mdl_ng_ensemble = fitcensemble(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','all');
mdl_ng_net = fitcnet(train_ng(:,2:25),train_ng.Anot_Class,'OptimizeHyperparameters','all');



%% Cancer v Noncancer
train_nc = train2;
loc = find(train_nc.Anot_Class == 'Atrophy' | train_nc.Anot_Class == 'HGPIN' ...
    | train_nc.Anot_Class == 'Seminal_Vesicles' | train_nc.Anot_Class == 'Tissue');
train_nc.Anot_Class(loc) = 'Normal';

loc = find(train2.Anot_Class == 'G3' | train2.Anot_Class == 'G4FG' ...
    | train2.Anot_Class == 'G4CG' | train2.Anot_Class == 'G5');
train_nc.Anot_Class(loc) = 'Cancer';
%
test_nc = test2;
loc = find(test_nc.Anot_Class == 'Atrophy' | test_nc.Anot_Class == 'HGPIN' ...
    | test_nc.Anot_Class == 'Seminal_Vesicles' | test_nc.Anot_Class == 'Tissue');
test_nc.Anot_Class(loc) = 'Normal';

loc = find(test2.Anot_Class == 'G3' | test2.Anot_Class == 'G4FG' ...
    | test2.Anot_Class == 'G4CG' | test2.Anot_Class == 'G5');
test_nc.Anot_Class(loc) = 'Cancer';
%
loc = find(feature_tbl_2235.Anot_Class == 'Atrophy' | feature_tbl_2235.Anot_Class == 'HGPIN' ...
    | feature_tbl_2235.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2235.Anot_Class == 'Tissue');
feature_tbl_2235.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_2235.Anot_Class == 'G3' | feature_tbl_2235.Anot_Class == 'G4FG' ...
    | feature_tbl_2235.Anot_Class == 'G4CG' | feature_tbl_2235.Anot_Class == 'G5');
feature_tbl_2235.Anot_Class(loc) = 'Cancer';
%
loc = find(feature_tbl_2181.Anot_Class == 'Atrophy' | feature_tbl_2181.Anot_Class == 'HGPIN' ...
    | feature_tbl_2181.Anot_Class == 'Seminal_Vesicles' | feature_tbl_2181.Anot_Class == 'Tissue');
feature_tbl_2181.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_2181.Anot_Class == 'G3' | feature_tbl_2181.Anot_Class == 'G4FG' ...
    | feature_tbl_2181.Anot_Class == 'G4CG' | feature_tbl_2181.Anot_Class == 'G5');
feature_tbl_2181.Anot_Class(loc) = 'Cancer';
%
loc = find(feature_tbl_1163.Anot_Class == 'Atrophy' | feature_tbl_1163.Anot_Class == 'HGPIN' ...
    | feature_tbl_1163.Anot_Class == 'Seminal_Vesicles' | feature_tbl_1163.Anot_Class == 'Tissue');
feature_tbl_1163.Anot_Class(loc) = 'Normal';

loc = find(feature_tbl_1163.Anot_Class == 'G3' | feature_tbl_1163.Anot_Class == 'G4FG' ...
    | feature_tbl_1163.Anot_Class == 'G4CG' | feature_tbl_1163.Anot_Class == 'G5');
feature_tbl_1163.Anot_Class(loc) = 'Cancer';

mdl_nc_nb = fitcnb(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','all');
mdl_nc_knn = fitcknn(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','all');
mdl_nc_tree = fitctree(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','all');
mdl_nc_ecoc = fitcecoc(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','auto');
mdl_nc_ensemble = fitcensemble(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','all');
mdl_nc_net = fitcnet(train_nc(:,2:25),train_nc.Anot_Class,'OptimizeHyperparameters','all');

%% All Gleason
mdl_nb = fitcnb(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','all');
mdl_knn = fitcknn(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','all');
mdl_tree = fitctree(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','all');
mdl_ecoc = fitcecoc(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','auto');
mdl_ensemble = fitcensemble(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','all');
mdl_net = fitcnet(train2(:,2:25),train2.Anot_Class,'OptimizeHyperparameters','all');

%% Check Model Performance
mdl_nb_pred = predict(mdl_nb, test2(:,2:25));
mdl_knn_pred = predict(mdl_knn, test2(:,2:25));
mdl_tree_pred = predict(mdl_tree, test2(:,2:25));
mdl_ecoc_pred = predict(mdl_ecoc, test2(:,2:25));
mdl_ensemble_pred = predict(mdl_ensemble, test2(:,2:25));
mdl_net_pred = predict(mdl_net, test2(:,2:25));

mdl_hln_nb_pred = predict(mdl_hln_nb, test_hln(:,2:25));
mdl_hln_knn_pred = predict(mdl_hln_knn, test_hln(:,2:25));
mdl_hln_tree_pred = predict(mdl_hln_tree, test_hln(:,2:25));
mdl_hln_ecoc_pred = predict(mdl_hln_ecoc, test_hln(:,2:25));
mdl_hln_ensemble_pred = predict(mdl_hln_ensemble, test_hln(:,2:25));
mdl_hln_net_pred = predict(mdl_hln_net, test_hln(:,2:25));

mdl_nc_nb_pred = predict(mdl_nc_nb, test_nc(:,2:25));
mdl_nc_knn_pred = predict(mdl_nc_knn, test_nc(:,2:25));
mdl_nc_tree_pred = predict(mdl_nc_tree, test_nc(:,2:25));
mdl_nc_ecoc_pred = predict(mdl_nc_ecoc, test_nc(:,2:25));
mdl_nc_ensemble_pred = predict(mdl_nc_ensemble, test_nc(:,2:25));
mdl_nc_net_pred = predict(mdl_nc_net, test_nc(:,2:25));

mdl_ng_nb_pred = predict(mdl_ng_nb, test_ng(:,2:25));
mdl_ng_knn_pred = predict(mdl_ng_knn, test_ng(:,2:25));
mdl_ng_tree_pred = predict(mdl_ng_tree, test_ng(:,2:25));
mdl_ng_ecoc_pred = predict(mdl_ng_ecoc, test_ng(:,2:25));
mdl_ng_ensemble_pred = predict(mdl_ng_ensemble, test_ng(:,2:25));
mdl_ng_net_pred = predict(mdl_ng_net, test_ng(:,2:25));
%
tbl_mdl_nb = calculate_error(test2.Anot_Class,mdl_nb_pred);
tbl_mdl_knn = calculate_error(test2.Anot_Class,mdl_knn_pred);
tbl_mdl_tree = calculate_error(test2.Anot_Class,mdl_tree_pred);
tbl_mdl_ecoc = calculate_error(test2.Anot_Class,mdl_ecoc_pred);
tbl_mdl_ensemble = calculate_error(test2.Anot_Class,mdl_ensemble_pred);
tbl_mdl_net = calculate_error(test2.Anot_Class,mdl_net_pred);

tbl_mdl_hln_nb = calculate_error(test_hln.Anot_Class,mdl_hln_nb_pred);
tbl_mdl_hln_knn = calculate_error(test_hln.Anot_Class,mdl_hln_knn_pred);
tbl_mdl_hln_tree = calculate_error(test_hln.Anot_Class,mdl_hln_tree_pred);
tbl_mdl_hln_ecoc = calculate_error(test_hln.Anot_Class,mdl_hln_ecoc_pred);
tbl_mdl_hln_ensemble = calculate_error(test_hln.Anot_Class,mdl_hln_ensemble_pred);
tbl_mdl_hln_net = calculate_error(test_hln.Anot_Class,mdl_hln_net_pred);

tbl_mdl_nc_nb = calculate_error(test_nc.Anot_Class,mdl_nc_nb_pred);
tbl_mdl_nc_knn = calculate_error(test_nc.Anot_Class,mdl_nc_knn_pred);
tbl_mdl_nc_tree = calculate_error(test_nc.Anot_Class,mdl_nc_tree_pred);
tbl_mdl_nc_ecoc = calculate_error(test_nc.Anot_Class,mdl_nc_ecoc_pred);
tbl_mdl_nc_ensemble = calculate_error(test_nc.Anot_Class,mdl_nc_ensemble_pred);
tbl_mdl_nc_net = calculate_error(test_nc.Anot_Class,mdl_nc_net_pred);

tbl_mdl_ng_nb = calculate_error(test_ng.Anot_Class,mdl_nc_nb_pred);
tbl_mdl_ng_knn = calculate_error(test_ng.Anot_Class,mdl_nc_knn_pred);
tbl_mdl_ng_tree = calculate_error(test_ng.Anot_Class,mdl_nc_tree_pred);
tbl_mdl_ng_ecoc = calculate_error(test_ng.Anot_Class,mdl_nc_ecoc_pred);
tbl_mdl_ng_ensemble = calculate_error(test_ng.Anot_Class,mdl_nc_ensemble_pred);
tbl_mdl_ng_net = calculate_error(test_ng.Anot_Class,mdl_nc_net_pred);
%
prec_mdl_nb = sum(tbl_mdl_nb.TP) / (sum(tbl_mdl_nb.TP) + sum(tbl_mdl_nb.FP));
prec_mdl_knn = sum(tbl_mdl_knn.TP) / (sum(tbl_mdl_knn.TP) + sum(tbl_mdl_knn.FP));
prec_mdl_tree = sum(tbl_mdl_tree.TP) / (sum(tbl_mdl_tree.TP) + sum(tbl_mdl_tree.FP));
prec_mdl_ecoc = sum(tbl_mdl_ecoc.TP) / (sum(tbl_mdl_ecoc.TP) + sum(tbl_mdl_ecoc.FP));
prec_mdl_ensemble = sum(tbl_mdl_ensemble.TP) / (sum(tbl_mdl_ensemble.TP) + sum(tbl_mdl_ensemble.FP));
prec_mdl_net = sum(tbl_mdl_net.TP) / (sum(tbl_mdl_net.TP) + sum(tbl_mdl_net.FP));

prec_mdl_hln_nb = sum(tbl_mdl_hln_nb.TP) / (sum(tbl_mdl_hln_nb.TP) + sum(tbl_mdl_hln_nb.FP));
prec_mdl_hln_knn = sum(tbl_mdl_hln_knn.TP) / (sum(tbl_mdl_hln_knn.TP) + sum(tbl_mdl_hln_knn.FP));
prec_mdl_hln_tree = sum(tbl_mdl_hln_tree.TP) / (sum(tbl_mdl_hln_tree.TP) + sum(tbl_mdl_hln_tree.FP));
prec_mdl_hln_ecoc = sum(tbl_mdl_hln_ecoc.TP) / (sum(tbl_mdl_hln_ecoc.TP) + sum(tbl_mdl_hln_ecoc.FP));
prec_mdl_hln_ensemble = sum(tbl_mdl_hln_ensemble.TP) / (sum(tbl_mdl_hln_ensemble.TP) + sum(tbl_mdl_hln_ensemble.FP));
prec_mdl_hln_net = sum(tbl_mdl_hln_net.TP) / (sum(tbl_mdl_hln_net.TP) + sum(tbl_mdl_hln_net.FP));

prec_mdl_nc_nb = sum(tbl_mdl_nc_nb.TP) / (sum(tbl_mdl_nc_nb.TP) + sum(tbl_mdl_nc_nb.FP));
prec_mdl_nc_knn = sum(tbl_mdl_nc_knn.TP) / (sum(tbl_mdl_nc_knn.TP) + sum(tbl_mdl_nc_knn.FP));
prec_mdl_nc_tree = sum(tbl_mdl_nc_tree.TP) / (sum(tbl_mdl_nc_tree.TP) + sum(tbl_mdl_nc_tree.FP));
prec_mdl_nc_ecoc = sum(tbl_mdl_nc_ecoc.TP) / (sum(tbl_mdl_nc_ecoc.TP) + sum(tbl_mdl_nc_ecoc.FP));
prec_mdl_nc_ensemble = sum(tbl_mdl_nc_ensemble.TP) / (sum(tbl_mdl_nc_ensemble.TP) + sum(tbl_mdl_nc_ensemble.FP));
prec_mdl_nc_net = sum(tbl_mdl_nc_net.TP) / (sum(tbl_mdl_nc_net.TP) + sum(tbl_mdl_nc_net.FP));

prec_mdl_ng_nb = sum(tbl_mdl_ng_nb.TP) / (sum(tbl_mdl_ng_nb.TP) + sum(tbl_mdl_ng_nb.FP));
prec_mdl_ng_knn = sum(tbl_mdl_ng_knn.TP) / (sum(tbl_mdl_ng_knn.TP) + sum(tbl_mdl_ng_knn.FP));
prec_mdl_ng_tree = sum(tbl_mdl_ng_tree.TP) / (sum(tbl_mdl_ng_tree.TP) + sum(tbl_mdl_ng_tree.FP));
prec_mdl_ng_ecoc = sum(tbl_mdl_ng_ecoc.TP) / (sum(tbl_mdl_ng_ecoc.TP) + sum(tbl_mdl_ng_ecoc.FP));
prec_mdl_ng_ensemble = sum(tbl_mdl_ng_ensemble.TP) / (sum(tbl_mdl_ng_ensemble.TP) + sum(tbl_mdl_ng_ensemble.FP));
prec_mdl_ng_net = sum(tbl_mdl_ng_net.TP) / (sum(tbl_mdl_ng_net.TP) + sum(tbl_mdl_ng_net.FP));
%
rec_mdl_nb = sum(tbl_mdl_nb.TP) / (sum(tbl_mdl_nb.TP) + sum(tbl_mdl_nb.FN));
rec_mdl_knn = sum(tbl_mdl_knn.TP) / (sum(tbl_mdl_knn.TP) + sum(tbl_mdl_knn.FN));
rec_mdl_tree = sum(tbl_mdl_tree.TP) / (sum(tbl_mdl_tree.TP) + sum(tbl_mdl_tree.FN));
rec_mdl_ecoc = sum(tbl_mdl_ecoc.TP) / (sum(tbl_mdl_ecoc.TP) + sum(tbl_mdl_ecoc.FN));
rec_mdl_ensemble = sum(tbl_mdl_ensemble.TP) / (sum(tbl_mdl_ensemble.TP) + sum(tbl_mdl_ensemble.FN));
rec_mdl_net = sum(tbl_mdl_net.TP) / (sum(tbl_mdl_net.TP) + sum(tbl_mdl_net.FN));

rec_mdl_hln_nb = sum(tbl_mdl_hln_nb.TP) / (sum(tbl_mdl_hln_nb.TP) + sum(tbl_mdl_hln_nb.FN));
rec_mdl_hln_knn = sum(tbl_mdl_hln_knn.TP) / (sum(tbl_mdl_hln_knn.TP) + sum(tbl_mdl_hln_knn.FN));
rec_mdl_hln_tree = sum(tbl_mdl_hln_tree.TP) / (sum(tbl_mdl_hln_tree.TP) + sum(tbl_mdl_hln_tree.FN));
rec_mdl_hln_ecoc = sum(tbl_mdl_hln_ecoc.TP) / (sum(tbl_mdl_hln_ecoc.TP) + sum(tbl_mdl_hln_ecoc.FN));
rec_mdl_hln_ensemble = sum(tbl_mdl_hln_ensemble.TP) / (sum(tbl_mdl_hln_ensemble.TP) + sum(tbl_mdl_hln_ensemble.FN));
rec_mdl_hln_net = sum(tbl_mdl_hln_net.TP) / (sum(tbl_mdl_hln_net.TP) + sum(tbl_mdl_hln_net.FN));

rec_mdl_nc_nb = sum(tbl_mdl_nc_nb.TP) / (sum(tbl_mdl_nc_nb.TP) + sum(tbl_mdl_nc_nb.FN));
rec_mdl_nc_knn = sum(tbl_mdl_nc_knn.TP) / (sum(tbl_mdl_nc_knn.TP) + sum(tbl_mdl_nc_knn.FN));
rec_mdl_nc_tree = sum(tbl_mdl_nc_tree.TP) / (sum(tbl_mdl_nc_tree.TP) + sum(tbl_mdl_nc_tree.FN));
rec_mdl_nc_ecoc = sum(tbl_mdl_nc_ecoc.TP) / (sum(tbl_mdl_nc_ecoc.TP) + sum(tbl_mdl_nc_ecoc.FN));
rec_mdl_nc_ensemble = sum(tbl_mdl_nc_ensemble.TP) / (sum(tbl_mdl_nc_ensemble.TP) + sum(tbl_mdl_nc_ensemble.FN));
rec_mdl_nc_net = sum(tbl_mdl_nc_net.TP) / (sum(tbl_mdl_nc_net.TP) + sum(tbl_mdl_nc_net.FN));

rec_mdl_ng_nb = sum(tbl_mdl_ng_nb.TP) / (sum(tbl_mdl_ng_nb.TP) + sum(tbl_mdl_ng_nb.FN));
rec_mdl_ng_knn = sum(tbl_mdl_ng_knn.TP) / (sum(tbl_mdl_ng_knn.TP) + sum(tbl_mdl_ng_knn.FN));
rec_mdl_ng_tree = sum(tbl_mdl_ng_tree.TP) / (sum(tbl_mdl_ng_tree.TP) + sum(tbl_mdl_ng_tree.FN));
rec_mdl_ng_ecoc = sum(tbl_mdl_ng_ecoc.TP) / (sum(tbl_mdl_ng_ecoc.TP) + sum(tbl_mdl_ng_ecoc.FN));
rec_mdl_ng_ensemble = sum(tbl_mdl_ng_ensemble.TP) / (sum(tbl_mdl_ng_ensemble.TP) + sum(tbl_mdl_ng_ensemble.FN));
rec_mdl_ng_net = sum(tbl_mdl_ng_net.TP) / (sum(tbl_mdl_ng_net.TP) + sum(tbl_mdl_ng_net.FN));
%
F_mdl_nb = 2*(prec_mdl_nb * rec_mdl_nb) / (prec_mdl_nb + rec_mdl_nb);
F_mdl_knn = 2*(prec_mdl_knn * rec_mdl_knn) / (prec_mdl_knn + rec_mdl_knn);
F_mdl_tree = 2*(prec_mdl_tree * rec_mdl_tree) / (prec_mdl_tree + rec_mdl_tree);
F_mdl_ecoc = 2*(prec_mdl_ecoc * rec_mdl_ecoc) / (prec_mdl_ecoc + rec_mdl_ecoc);
F_mdl_ensemble = 2*(prec_mdl_ensemble * rec_mdl_ensemble) / (prec_mdl_ensemble + rec_mdl_ensemble);
F_mdl_net = 2*(prec_mdl_net * rec_mdl_net) / (prec_mdl_net + rec_mdl_net);

F_mdl_hln_nb = 2*(prec_mdl_hln_nb * rec_mdl_hln_nb) / (prec_mdl_hln_nb + rec_mdl_hln_nb);
F_mdl_hln_knn = 2*(prec_mdl_hln_knn * rec_mdl_hln_knn) / (prec_mdl_hln_knn + rec_mdl_hln_knn);
F_mdl_hln_tree = 2*(prec_mdl_hln_tree * rec_mdl_hln_tree) / (prec_mdl_hln_tree + rec_mdl_hln_tree);
F_mdl_hln_ecoc = 2*(prec_mdl_hln_ecoc * rec_mdl_hln_ecoc) / (prec_mdl_hln_ecoc + rec_mdl_hln_ecoc);
F_mdl_hln_ensemble = 2*(prec_mdl_hln_ensemble * rec_mdl_hln_ensemble) / (prec_mdl_hln_ensemble + rec_mdl_hln_ensemble);
F_mdl_hln_net = 2*(prec_mdl_hln_net * rec_mdl_hln_net) / (prec_mdl_hln_net + rec_mdl_hln_net);

F_mdl_nc_nb = 2*(prec_mdl_nc_nb * rec_mdl_nc_nb) / (prec_mdl_nc_nb + rec_mdl_nc_nb);
F_mdl_nc_knn = 2*(prec_mdl_nc_knn * rec_mdl_nc_knn) / (prec_mdl_nc_knn + rec_mdl_nc_knn);
F_mdl_nc_tree = 2*(prec_mdl_nc_tree * rec_mdl_nc_tree) / (prec_mdl_nc_tree + rec_mdl_nc_tree);
F_mdl_nc_ecoc = 2*(prec_mdl_nc_ecoc * rec_mdl_nc_ecoc) / (prec_mdl_nc_ecoc + rec_mdl_nc_ecoc);
F_mdl_nc_ensemble = 2*(prec_mdl_nc_ensemble * rec_mdl_nc_ensemble) / (prec_mdl_nc_ensemble + rec_mdl_nc_ensemble);
F_mdl_nc_net = 2*(prec_mdl_nc_net * rec_mdl_nc_net) / (prec_mdl_nc_net + rec_mdl_nc_net);

F_mdl_ng_nb = 2*(prec_mdl_ng_nb * rec_mdl_ng_nb) / (prec_mdl_ng_nb + rec_mdl_ng_nb);
F_mdl_ng_knn = 2*(prec_mdl_ng_knn * rec_mdl_ng_knn) / (prec_mdl_ng_knn + rec_mdl_ng_knn);
F_mdl_ng_tree = 2*(prec_mdl_ng_tree * rec_mdl_ng_tree) / (prec_mdl_ng_tree + rec_mdl_ng_tree);
F_mdl_ng_ecoc = 2*(prec_mdl_ng_ecoc * rec_mdl_ng_ecoc) / (prec_mdl_ng_ecoc + rec_mdl_ng_ecoc);
F_mdl_ng_ensemble = 2*(prec_mdl_ng_ensemble * rec_mdl_ng_ensemble) / (prec_mdl_ng_ensemble + rec_mdl_ng_ensemble);
F_mdl_ng_net = 2*(prec_mdl_ng_net * rec_mdl_ng_net) / (prec_mdl_ng_net + rec_mdl_ng_net);
%
acc_mdl_nb = (sum(tbl_mdl_nb.TP )+ sum(tbl_mdl_nb.TN)) / (sum(tbl_mdl_nb.TP) + sum(tbl_mdl_nb.FP) + sum(tbl_mdl_nb.TN) + sum(tbl_mdl_nb.FN));
acc_mdl_knn = (sum(tbl_mdl_knn.TP )+ sum(tbl_mdl_knn.TN)) / (sum(tbl_mdl_knn.TP) + sum(tbl_mdl_knn.FP) + sum(tbl_mdl_knn.TN) + sum(tbl_mdl_knn.FN));
acc_mdl_tree = (sum(tbl_mdl_tree.TP )+ sum(tbl_mdl_tree.TN)) / (sum(tbl_mdl_tree.TP) + sum(tbl_mdl_tree.FP) + sum(tbl_mdl_tree.TN) + sum(tbl_mdl_tree.FN));
acc_mdl_ecoc = (sum(tbl_mdl_ecoc.TP )+ sum(tbl_mdl_ecoc.TN)) / (sum(tbl_mdl_ecoc.TP) + sum(tbl_mdl_ecoc.FP) + sum(tbl_mdl_ecoc.TN) + sum(tbl_mdl_ecoc.FN));
acc_mdl_ensemble = (sum(tbl_mdl_ensemble.TP )+ sum(tbl_mdl_ensemble.TN)) / (sum(tbl_mdl_ensemble.TP) + sum(tbl_mdl_ensemble.FP) + sum(tbl_mdl_ensemble.TN) + sum(tbl_mdl_ensemble.FN));
acc_mdl_net = (sum(tbl_mdl_net.TP )+ sum(tbl_mdl_net.TN)) / (sum(tbl_mdl_net.TP) + sum(tbl_mdl_net.FP) + sum(tbl_mdl_net.TN) + sum(tbl_mdl_net.FN));

acc_mdl_hln_nb = (sum(tbl_mdl_hln_nb.TP )+ sum(tbl_mdl_hln_nb.TN)) / (sum(tbl_mdl_hln_nb.TP) + sum(tbl_mdl_hln_nb.FP) + sum(tbl_mdl_hln_nb.TN) + sum(tbl_mdl_hln_nb.FN));
acc_mdl_hln_knn = (sum(tbl_mdl_hln_knn.TP )+ sum(tbl_mdl_hln_knn.TN)) / (sum(tbl_mdl_hln_knn.TP) + sum(tbl_mdl_hln_knn.FP) + sum(tbl_mdl_hln_knn.TN) + sum(tbl_mdl_hln_knn.FN));
acc_mdl_hln_tree = (sum(tbl_mdl_hln_tree.TP )+ sum(tbl_mdl_hln_tree.TN)) / (sum(tbl_mdl_hln_tree.TP) + sum(tbl_mdl_hln_tree.FP) + sum(tbl_mdl_hln_tree.TN) + sum(tbl_mdl_hln_tree.FN));
acc_mdl_hln_ecoc = (sum(tbl_mdl_hln_ecoc.TP )+ sum(tbl_mdl_hln_ecoc.TN)) / (sum(tbl_mdl_hln_ecoc.TP) + sum(tbl_mdl_hln_ecoc.FP) + sum(tbl_mdl_hln_ecoc.TN) + sum(tbl_mdl_hln_ecoc.FN));
acc_mdl_hln_ensemble = (sum(tbl_mdl_hln_ensemble.TP )+ sum(tbl_mdl_hln_ensemble.TN)) / (sum(tbl_mdl_hln_ensemble.TP) + sum(tbl_mdl_hln_ensemble.FP) + sum(tbl_mdl_hln_ensemble.TN) + sum(tbl_mdl_hln_ensemble.FN));
acc_mdl_hln_net = (sum(tbl_mdl_hln_net.TP )+ sum(tbl_mdl_hln_net.TN)) / (sum(tbl_mdl_hln_net.TP) + sum(tbl_mdl_hln_net.FP) + sum(tbl_mdl_hln_net.TN) + sum(tbl_mdl_hln_net.FN));

acc_mdl_nc_nb = (sum(tbl_mdl_nc_nb.TP )+ sum(tbl_mdl_nc_nb.TN)) / (sum(tbl_mdl_nc_nb.TP) + sum(tbl_mdl_nc_nb.FP) + sum(tbl_mdl_nc_nb.TN) + sum(tbl_mdl_nc_nb.FN));
acc_mdl_nc_knn = (sum(tbl_mdl_nc_knn.TP )+ sum(tbl_mdl_nc_knn.TN)) / (sum(tbl_mdl_nc_knn.TP) + sum(tbl_mdl_nc_knn.FP) + sum(tbl_mdl_nc_knn.TN) + sum(tbl_mdl_nc_knn.FN));
acc_mdl_nc_tree = (sum(tbl_mdl_nc_tree.TP )+ sum(tbl_mdl_nc_tree.TN)) / (sum(tbl_mdl_nc_tree.TP) + sum(tbl_mdl_nc_tree.FP) + sum(tbl_mdl_nc_tree.TN) + sum(tbl_mdl_nc_tree.FN));
acc_mdl_nc_ecoc = (sum(tbl_mdl_nc_ecoc.TP )+ sum(tbl_mdl_nc_ecoc.TN)) / (sum(tbl_mdl_nc_ecoc.TP) + sum(tbl_mdl_nc_ecoc.FP) + sum(tbl_mdl_nc_ecoc.TN) + sum(tbl_mdl_nc_ecoc.FN));
acc_mdl_nc_ensemble = (sum(tbl_mdl_nc_ensemble.TP )+ sum(tbl_mdl_nc_ensemble.TN)) / (sum(tbl_mdl_nc_ensemble.TP) + sum(tbl_mdl_nc_ensemble.FP) + sum(tbl_mdl_nc_ensemble.TN) + sum(tbl_mdl_nc_ensemble.FN));
acc_mdl_nc_net = (sum(tbl_mdl_nc_net.TP )+ sum(tbl_mdl_nc_net.TN)) / (sum(tbl_mdl_nc_net.TP) + sum(tbl_mdl_nc_net.FP) + sum(tbl_mdl_nc_net.TN) + sum(tbl_mdl_nc_net.FN));

acc_mdl_ng_nb = (sum(tbl_mdl_ng_nb.TP )+ sum(tbl_mdl_ng_nb.TN)) / (sum(tbl_mdl_ng_nb.TP) + sum(tbl_mdl_ng_nb.FP) + sum(tbl_mdl_ng_nb.TN) + sum(tbl_mdl_ng_nb.FN));
acc_mdl_ng_knn = (sum(tbl_mdl_ng_knn.TP )+ sum(tbl_mdl_ng_knn.TN)) / (sum(tbl_mdl_ng_knn.TP) + sum(tbl_mdl_ng_knn.FP) + sum(tbl_mdl_ng_knn.TN) + sum(tbl_mdl_ng_knn.FN));
acc_mdl_ng_tree = (sum(tbl_mdl_ng_tree.TP )+ sum(tbl_mdl_ng_tree.TN)) / (sum(tbl_mdl_ng_tree.TP) + sum(tbl_mdl_ng_tree.FP) + sum(tbl_mdl_ng_tree.TN) + sum(tbl_mdl_ng_tree.FN));
acc_mdl_ng_ecoc = (sum(tbl_mdl_ng_ecoc.TP )+ sum(tbl_mdl_ng_ecoc.TN)) / (sum(tbl_mdl_ng_ecoc.TP) + sum(tbl_mdl_ng_ecoc.FP) + sum(tbl_mdl_ng_ecoc.TN) + sum(tbl_mdl_ng_ecoc.FN));
acc_mdl_ng_ensemble = (sum(tbl_mdl_ng_ensemble.TP )+ sum(tbl_mdl_ng_ensemble.TN)) / (sum(tbl_mdl_ng_ensemble.TP) + sum(tbl_mdl_ng_ensemble.FP) + sum(tbl_mdl_ng_ensemble.TN) + sum(tbl_mdl_ng_ensemble.FN));
acc_mdl_ng_net = (sum(tbl_mdl_ng_net.TP )+ sum(tbl_mdl_ng_net.TN)) / (sum(tbl_mdl_ng_net.TP) + sum(tbl_mdl_ng_net.FP) + sum(tbl_mdl_ng_net.TN) + sum(tbl_mdl_ng_net.FN));
%
spec_mdl_nb = sum(tbl_mdl_nb.TN) / (sum(tbl_mdl_nb.TN) + sum(tbl_mdl_nb.FP));
spec_mdl_knn = sum(tbl_mdl_knn.TN) / (sum(tbl_mdl_knn.TN) + sum(tbl_mdl_knn.FP));
spec_mdl_tree = sum(tbl_mdl_tree.TN) / (sum(tbl_mdl_tree.TN) + sum(tbl_mdl_tree.FP));
spec_mdl_ecoc = sum(tbl_mdl_ecoc.TN) / (sum(tbl_mdl_ecoc.TN) + sum(tbl_mdl_ecoc.FP));
spec_mdl_ensemble = sum(tbl_mdl_ensemble.TN) / (sum(tbl_mdl_ensemble.TN) + sum(tbl_mdl_ensemble.FP));
spec_mdl_net = sum(tbl_mdl_net.TN) / (sum(tbl_mdl_net.TN) + sum(tbl_mdl_net.FP));

spec_mdl_hln_nb = sum(tbl_mdl_hln_nb.TN) / (sum(tbl_mdl_hln_nb.TN) + sum(tbl_mdl_hln_nb.FP));
spec_mdl_hln_knn = sum(tbl_mdl_hln_knn.TN) / (sum(tbl_mdl_hln_knn.TN) + sum(tbl_mdl_hln_knn.FP));
spec_mdl_hln_tree = sum(tbl_mdl_hln_tree.TN) / (sum(tbl_mdl_hln_tree.TN) + sum(tbl_mdl_hln_tree.FP));
spec_mdl_hln_ecoc = sum(tbl_mdl_hln_ecoc.TN) / (sum(tbl_mdl_hln_ecoc.TN) + sum(tbl_mdl_hln_ecoc.FP));
spec_mdl_hln_ensemble = sum(tbl_mdl_hln_ensemble.TN) / (sum(tbl_mdl_hln_ensemble.TN) + sum(tbl_mdl_hln_ensemble.FP));
spec_mdl_hln_net = sum(tbl_mdl_hln_net.TN) / (sum(tbl_mdl_hln_net.TN) + sum(tbl_mdl_hln_net.FP));

spec_mdl_nc_nb = sum(tbl_mdl_nc_nb.TN) / (sum(tbl_mdl_nc_nb.TN) + sum(tbl_mdl_nc_nb.FP));
spec_mdl_nc_knn = sum(tbl_mdl_nc_knn.TN) / (sum(tbl_mdl_nc_knn.TN) + sum(tbl_mdl_nc_knn.FP));
spec_mdl_nc_tree = sum(tbl_mdl_nc_tree.TN) / (sum(tbl_mdl_nc_tree.TN) + sum(tbl_mdl_nc_tree.FP));
spec_mdl_nc_ecoc = sum(tbl_mdl_nc_ecoc.TN) / (sum(tbl_mdl_nc_ecoc.TN) + sum(tbl_mdl_nc_ecoc.FP));
spec_mdl_nc_ensemble = sum(tbl_mdl_nc_ensemble.TN) / (sum(tbl_mdl_nc_ensemble.TN) + sum(tbl_mdl_nc_ensemble.FP));
spec_mdl_nc_net = sum(tbl_mdl_nc_net.TN) / (sum(tbl_mdl_nc_net.TN) + sum(tbl_mdl_nc_net.FP));

spec_mdl_ng_nb = sum(tbl_mdl_ng_nb.TN) / (sum(tbl_mdl_ng_nb.TN) + sum(tbl_mdl_ng_nb.FP));
spec_mdl_ng_knn = sum(tbl_mdl_ng_knn.TN) / (sum(tbl_mdl_ng_knn.TN) + sum(tbl_mdl_ng_knn.FP));
spec_mdl_ng_tree = sum(tbl_mdl_ng_tree.TN) / (sum(tbl_mdl_ng_tree.TN) + sum(tbl_mdl_ng_tree.FP));
spec_mdl_ng_ecoc = sum(tbl_mdl_ng_ecoc.TN) / (sum(tbl_mdl_ng_ecoc.TN) + sum(tbl_mdl_ng_ecoc.FP));
spec_mdl_ng_ensemble = sum(tbl_mdl_ng_ensemble.TN) / (sum(tbl_mdl_ng_ensemble.TN) + sum(tbl_mdl_ng_ensemble.FP));
spec_mdl_ng_net = sum(tbl_mdl_ng_net.TN) / (sum(tbl_mdl_ng_net.TN) + sum(tbl_mdl_ng_net.FP));

%% Make Final Tables
tbl_mdl = [tbl_mdl_nb;tbl_mdl_knn;tbl_mdl_tree;tbl_mdl_ecoc;tbl_mdl_ensemble;tbl_mdl_net;...
    tbl_mdl_hln_nb;tbl_mdl_hln_knn;tbl_mdl_hln_tree;tbl_mdl_hln_ecoc;tbl_mdl_hln_ensemble;tbl_mdl_hln_net;...
    tbl_mdl_nc_nb;tbl_mdl_nc_knn;tbl_mdl_nc_tree;tbl_mdl_nc_ecoc;tbl_mdl_nc_ensemble;tbl_mdl_nc_net;...
    tbl_mdl_ng_nb;tbl_mdl_ng_knn;tbl_mdl_ng_tree;tbl_mdl_ng_ecoc;tbl_mdl_ng_ensemble;tbl_mdl_ng_net];

overall_tbl = [prec_mdl_nb, rec_mdl_nb, spec_mdl_nb, F_mdl_nb, acc_mdl_nb;...
    prec_mdl_knn, rec_mdl_knn, spec_mdl_knn, F_mdl_knn, acc_mdl_knn;...
    prec_mdl_tree, rec_mdl_tree, spec_mdl_tree, F_mdl_tree, acc_mdl_tree;...
    prec_mdl_ecoc, rec_mdl_ecoc, spec_mdl_ecoc, F_mdl_ecoc, acc_mdl_ecoc;...
    prec_mdl_ensemble, rec_mdl_ensemble, spec_mdl_ensemble, F_mdl_ensemble, acc_mdl_ensemble;...
    prec_mdl_net, rec_mdl_net, spec_mdl_net, F_mdl_net, acc_mdl_net;...
    
    prec_mdl_hln_nb, rec_mdl_hln_nb, spec_mdl_hln_nb, F_mdl_hln_nb, acc_mdl_hln_nb;...
    prec_mdl_hln_knn, rec_mdl_hln_knn, spec_mdl_hln_knn, F_mdl_hln_knn, acc_mdl_hln_knn;...
    prec_mdl_hln_tree, rec_mdl_hln_tree, spec_mdl_hln_tree, F_mdl_hln_tree, acc_mdl_hln_tree;...
    prec_mdl_hln_ecoc, rec_mdl_hln_ecoc, spec_mdl_hln_ecoc, F_mdl_hln_ecoc, acc_mdl_hln_ecoc;...
    prec_mdl_hln_ensemble, rec_mdl_hln_ensemble, spec_mdl_hln_ensemble, F_mdl_hln_ensemble, acc_mdl_hln_ensemble;...
    prec_mdl_hln_net, rec_mdl_hln_net, spec_mdl_hln_net, F_mdl_hln_net, acc_mdl_hln_net;...
    
    prec_mdl_nc_nb, rec_mdl_nc_nb, spec_mdl_nc_nb, F_mdl_nc_nb, acc_mdl_nc_nb;...
    prec_mdl_nc_knn, rec_mdl_nc_knn, spec_mdl_nc_knn, F_mdl_nc_knn, acc_mdl_nc_knn;...
    prec_mdl_nc_tree, rec_mdl_nc_tree, spec_mdl_nc_tree, F_mdl_nc_tree, acc_mdl_nc_tree;...
    prec_mdl_nc_ecoc, rec_mdl_nc_ecoc, spec_mdl_nc_ecoc, F_mdl_nc_ecoc, acc_mdl_nc_ecoc;...
    prec_mdl_nc_ensemble, rec_mdl_nc_ensemble, spec_mdl_nc_ensemble, F_mdl_nc_ensemble, acc_mdl_nc_ensemble;...
    prec_mdl_nc_net, rec_mdl_nc_net, spec_mdl_nc_net, F_mdl_nc_net, acc_mdl_nc_net;...
    
    prec_mdl_ng_nb, rec_mdl_ng_nb, spec_mdl_ng_nb, F_mdl_ng_nb, acc_mdl_ng_nb;...
    prec_mdl_ng_knn, rec_mdl_ng_knn, spec_mdl_ng_knn, F_mdl_ng_knn, acc_mdl_ng_knn;...
    prec_mdl_ng_tree, rec_mdl_ng_tree, spec_mdl_ng_tree, F_mdl_ng_tree, acc_mdl_ng_tree;...
    prec_mdl_ng_ecoc, rec_mdl_ng_ecoc, spec_mdl_ng_ecoc, F_mdl_ng_ecoc, acc_mdl_ng_ecoc;...
    prec_mdl_ng_ensemble, rec_mdl_ng_ensemble, spec_mdl_ng_ensemble, F_mdl_ng_ensemble, acc_mdl_ng_ensemble;...
    prec_mdl_ng_net, rec_mdl_ng_net, spec_mdl_ng_net, F_mdl_ng_net, acc_mdl_ng_net];

%
figure(); confusionchart(char(test2.Anot_Class), char(mdl_nb_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test2.Anot_Class), char(mdl_knn_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test2.Anot_Class), char(mdl_tree_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test2.Anot_Class), char(mdl_ecoc_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test2.Anot_Class), char(mdl_ensemble_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test2.Anot_Class), char(mdl_net_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');

figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_nb_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_knn_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_tree_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_ecoc_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_ensemble_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_hln.Anot_Class), char(mdl_hln_net_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');

figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_nb_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_knn_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_tree_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_ecoc_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_ensemble_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_nc.Anot_Class), char(mdl_nc_net_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');

figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_nb_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_knn_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_tree_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_ecoc_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_ensemble_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure(); confusionchart(char(test_ng.Anot_Class), char(mdl_ng_net_pred), 'RowSummary','row-normalized','ColumnSummary','column-normalized');

%% Demonstration models
figure
confusionchart(feature_tbl_2235.Anot_Class,predict(mdl_auto_nc_unbal,feature_tbl_2235(:,(8:end))),'RowSummary','row-normalized','ColumnSummary','column-normalized')
title('confusion matrix 2235')
figure
confusionchart(feature_tbl_2181.Anot_Class,predict(mdl_auto_nc_unbal,feature_tbl_2181(:,(8:end))),'RowSummary','row-normalized','ColumnSummary','column-normalized')
title('confusion matrix 2181')


mdl_2181 = 1 - loss(mdl_auto_nc_unbal,feature_tbl_2181(:,[2,8:end]),'Anot_Class');
mdl_2181_lbs = predict(mdl_auto_nc_unbal,feature_tbl_2181(:,[8:end]));

mdl_2235 = 1 - loss(mdl_auto_nc_unbal,feature_tbl_2235(:,[4,8:end]),'Anot_Class');
mdl_2235_lbs = predict(mdl_auto_nc_unbal,feature_tbl_2235(:,[8:end]));

mdl_1163 = 1 - loss(mdl_auto_nc_unbal,feature_tbl_1163(:,[2,6:end]),'Anot_Class');
mdl_1163_lbs = predict(mdl_auto_nc_unbal,feature_tbl_1163(:,[6:end]));

%%%%%%
mdl_2181 = 1 - loss(mdl_auto_unbal_ng,feature_tbl_2181(:,[2,8:end]),'Anot_Class');
mdl_2181_lbs = predict(new_mdl,table2array(feature_tbl_2181(:,8:end)));

mdl_2235 = 1 - loss(mdl_auto_unbal_ng,feature_tbl_2235(:,4,8:end),'Anot_Class');
mdl_2235_lbs = predict(new_mdl,table2array(feature_tbl_2235(:,8:end)));

mdl_1163 = 1 - loss(mdl_auto_unbal_ng,feature_tbl_1163(:,2,6:end),'Anot_Class');
mdl_1163_lbs = predict(new_mdl,table2array(feature_tbl_1163(:,6:end)));

figure
confusionchart(train_ng.Anot_Class,predict(mdl_auto_unbal_ng,table2array(train_ng(:,2:25))),'RowSummary','row-normalized','ColumnSummary','column-normalized')
title('confusion matrix test')


%% Create Whole Mounts
% anot_names = {'Atrophy','Seminal_Vesicles','HGPIN','Tissue','G3','G4CG','G4FG','G5','Background'};
% map = [0 0 0 % atrophy
%     0 0 255
%     255 120 0
%     220 140 200 %%
%     50 255 0
%     255 20 255
%     255 250 20
%     0 255 255
%     220 222 221]; 
% 
anot_names = {'Normal','G3','G4CG','G4FG','G5','Background'};
map = [220 140 200
    50 255 0
    255 20 255
    255 250 20
    0 255 255
    220 222 221];

% anot_names = {'Normal','Cancer','Background'};
% map = [220 140 200
%     50 255 0
%     255 20 255
%     255 250 20
%     0 255 255
%     255 255 255];
% %     220 222 221];

anot_names = {'Normal','Low','High','Background'};
map = [220 140 200
    255 0 0
    0 0 255
    220 222 221];
map = uint8(map);

feature_tbl_2181 = [feature_tbl_2181 table(mdl_2181_lbs)];
feature_tbl_2181.Properties.VariableNames{32} = 'prediction';
feature_tbl_2181.prediction = categorical(feature_tbl_2181.prediction);

for i = 1:numel(anot_names)
    loc = find(feature_tbl_2181.prediction == anot_names{i});
    feature_tbl_2181(loc,33) = {i};
end
feature_tbl_2181.Properties.VariableNames{33} = 'pred_indx';

wm2181 = zeros(34,41);
for h = 1:height(feature_tbl_2181)
    tile = feature_tbl_2181.Tile_Name{h};
    name = strsplit(tile,{'_','.'});
    i = str2num(name{4}(2:end));
    j = str2num(name{5}(2:end));
    wm2181(j,i) = feature_tbl_2181.pred_indx(h);
end
wm2181(find(wm2181==0)) = numel(anot_names);

figure();
imagesc(wm2181); colormap(map); title('classifier');


feature_tbl_2235 = [feature_tbl_2235 table(mdl_2235_lbs)];
feature_tbl_2235.Tile_Name = cellstr(feature_tbl_2235.Tile_Name);
feature_tbl_2235.Properties.VariableNames{32} = 'prediction';
feature_tbl_2235.prediction = categorical(feature_tbl_2235.prediction);
for i = 1:numel(anot_names)
    loc = find(feature_tbl_2235.prediction == anot_names{i});
    feature_tbl_2235(loc,33) = {i};
end
feature_tbl_2235.Properties.VariableNames{33} = 'pred_indx';


wm2235 = zeros(24,31);
for h = 1:height(feature_tbl_2235)
    tile = feature_tbl_2235.Tile_Name{h};
    name = strsplit(tile,{'_','.'});
    i = str2num(name{4}(2:end));
    j = str2num(name{5}(2:end));
    wm2235(j,i) = feature_tbl_2235.pred_indx(h);
end

wm2235(find(wm2235==0)) = numel(anot_names);

figure();
imagesc(wm2235); colormap(map); title('classifier');


feature_tbl_1163 = [feature_tbl_1163 table(mdl_1163_lbs)];
feature_tbl_1163.Tile_Name = cellstr(feature_tbl_1163.Tile_Name);
feature_tbl_1163.Properties.VariableNames{30} = 'prediction';
feature_tbl_1163.prediction = categorical(feature_tbl_1163.prediction);
for i = 1:numel(anot_names)
    loc = find(feature_tbl_1163.prediction == anot_names{i});
    feature_tbl_1163(loc,31) = {i};
end
feature_tbl_1163.Properties.VariableNames{31} = 'pred_indx';


wm1163 = zeros(39,33);
for h = 1:height(feature_tbl_1163)
    tile = feature_tbl_1163.Tile_Name{h};
    name = strsplit(tile,{'_','.'});
    i = str2num(name{4}(2:end));
    j = str2num(name{5}(2:end));
    wm1163(j,i) = feature_tbl_1163.pred_indx(h);
end

wm1163(find(wm1163==0)) = numel(anot_names);

figure();
imagesc(wm1163); colormap(map); title('classifier');


anot_2181 = imread('/Volumes/Hera/Prostate_Data/2181/Hist/9/Olympus/large_recon_8_ANOT.tif');
anot_2181 = anot_2181(:,:,1:3);
whole_image2181 = imread('/Volumes/Hera/Prostate_Data/2181/Hist/9/Olympus/large_recon_8_nowhite.tiff');

[Seminal_Vesicles,Atrophy,HGPIN,G3,G4FG,G4CG,G5,tissue_mask] = color_thresh_masks(anot_2181,whole_image2181,0,0);

% Normal vs Cancer
ken2181 = (Atrophy+HGPIN+Seminal_Vesicles+tissue_mask) + 2*G3 + 3*G4CG + 4*G4FG + 5*G5;
ken2181(find(ken2181==0)) = numel(anot_names);
% figure(); imagesc(ken2181); colormap(map);

atari_2181 = imresize(ken2181,size(wm2181),'nearest');
% figure(); imagesc(atari_2181); colormap(map)

wm2181(find(atari_2181 == numel(anot_names))) = numel(anot_names);
atari_2181(find(wm2181 == numel(anot_names))) = numel(anot_names);
figure(); 
subplot(121); imagesc(wm2181); colormap(map); title('classifier'); axis image
subplot(122); imagesc(atari_2181); colormap(map); title('ken'); axis image

for i = 1:numel(anot_names)
    comp_table2(i,1) = {anot_names{i}};
    comp_table2(i,2) = {numel(find(atari_2181 == i))}; % Ken volume
    comp_table2(i,3) = {numel(find(wm2181 == i))}; % Classifier volume
    comp_table2(i,4) = {dice((atari_2181==i),(wm2181==i))};
end

figure
b = bar(categorical(anot_names),cell2mat(comp_table2(:,2:3)));
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('2181 Volume');



anot_1163 = imread('/Volumes/Hera/Prostate_Data/1163/Hist/7/Olympus/large_recon_8_ANOT.tif');
anot_1163= anot_1163(:,:,1:3);
whole_image1163 = imread('/Volumes/Hera/Prostate_Data/1163/Hist/7/Olympus/large_recon_8_nowhite.tiff');

[Seminal_Vesicles,Atrophy,HGPIN,G3,G4FG,G4CG,G5,tissue_mask] = color_thresh_masks(anot_1163,whole_image1163,0,0);

% % % All Annotation Classes
% ken2181 = Atrophy + 2*Seminal_Vesicles + 3*HGPIN + 4*tissue_mask + 5*G3 + 6*G4CG + 7*G4FG + 8*G5;
% ken2181(find(ken2181==0)) = numel(anot_names);
% figure(); imagesc(ken2181); colormap(map)
    
% Normal vs Cancer
ken1163 = (Atrophy+HGPIN+Seminal_Vesicles+tissue_mask) + 2*G3 + 3*G4CG + 4*G4FG + 5*G5;
ken1163(find(ken1163==0)) = numel(anot_names);
% figure(); imagesc(ken2181); colormap(map);

atari_1163 = imresize(ken1163,size(wm1163),'nearest');
% figure(); imagesc(atari_2181); colormap(map)

wm1163(find(atari_1163 == numel(anot_names))) = numel(anot_names);
atari_1163(find(wm1163 == numel(anot_names))) = numel(anot_names);
figure(); 
subplot(121); imagesc(wm1163); colormap(map); title('classifier'); axis image
subplot(122); imagesc(atari_1163); colormap(map); title('ken'); axis image

for i = 1:numel(anot_names)
    comp_table2(i,1) = {anot_names{i}};
    comp_table2(i,2) = {numel(find(atari_2181 == i))}; % Ken volume
    comp_table2(i,3) = {numel(find(wm2181 == i))}; % Classifier volume
    comp_table2(i,4) = {dice((atari_2181==i),(wm2181==i))};
end

figure
b = bar(categorical(anot_names),cell2mat(comp_table2(:,2:3)));
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
title('2181 Volume');
