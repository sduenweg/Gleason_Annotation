function tbl = calculate_error(true,pred)

% Calculate TP, FP, TN, FN, and class accuracy for matrices of true and
% predicted labels

class = unique(true);

for i = 1:numel(class)

    TP(i,1) = sum(pred==class(i) & true==class(i)) / sum(pred==class(i)) * 100;
    FP(i,1) = sum(pred==class(i) & true~=class(i)) / sum(pred==class(i)) * 100;

    TN(i,1) = sum(pred~=class(i) & true~=class(i)) / sum(pred~=class(i)) * 100;
    FN(i,1) = sum((pred~=class(i)) & (true==class(i))) / sum(pred~=class(i))  * 100;

    class_acc(i,1) = sum(pred==class(i) & true==class(i))/sum(true==class(i));
end

TP(isnan(TP)) = 0;
FP(isnan(FP)) = 0;
TN(isnan(TN)) = 0;
FN(isnan(FN)) = 0;


tbl = array2table([TP, FP, TN, FN,class_acc]);
tbl.Class = class;
tbl.Properties.VariableNames{1} = 'TP';
tbl.Properties.VariableNames{2} = 'FP';
tbl.Properties.VariableNames{3} = 'TN';
tbl.Properties.VariableNames{4} = 'FN';
tbl.Properties.VariableNames{5} = 'Class_Acc';
tbl = movevars(tbl, 'Class', 'Before', 'TP');

end
    


