function [t_pred, score, y_target] = evaluateMBIMDS(network,imds)

%%
% To be run on the following networks:
% network = '1024x1024_34_true_5_class_prostate.onnx'
% network = '1024x1024_50_true_5_class_prostate.onnx'
% network = '1024x1024_152_true_5_class_prostate.onnx'
%%
% Use michaels original test set for 1024 images
% imds = imageDatastore('ABSOLUTE_PATH','IncludeSubfolders',true,'labelsource','foldernames');
% Absolute path should be the top level folder containing folders that
% seperate classes
%%
% For ease of use, set code directory as available path and run each
% network in its own network folder. This will ensure the saves are
% appropriately placed.

warning off
net = importONNXNetwork(network,'OutputLayerType','classification');
warning on

imdsT = transform(imds,@(x) preNetNorm(x));

[t_pred, score] = classify(net,imdsT);

y_target = imdsT.UnderlyingDatastores{1,1}.Labels;

csvwrite(sprintf('%s_prediction.csv',network(1:end-5)),grp2idx(t_pred));
csvwrite(sprintf('%s_score.csv',network(1:end-5)),score);
csvwrite(sprintf('%s_target.csv',network(1:end-5)),grp2idx(y_target));
