function [tpr,auc,acc,cm] = cross_val_ens_tree_confmat(xdata,ydata,folds,niter,posclass,gridlen)
% Common grid for differente ROC curves
% gridlen = 1000;
fpr = linspace(0,1,gridlen);

tpr = zeros(gridlen,niter);
auc = zeros(niter,1);
acc = zeros(niter,1);
cm = zeros(2,2,niter);

% t = templateTree('MaxNumSplits',10,'MergeLeaves','off','MinLeafSize',1,...
%     'MinParentSize',2,'PredictorSelection','allsplits','Prune','off',...
%     'Surrogate','off','type','classification','SplitCriterion','twoing');
t = templateTree('Surrogate','all');

parfor i=1:niter
    cv = cvpartition(ydata,'Kfold',folds);
    aux_tpr2 =zeros(1,gridlen);
    aux_auc2 =0;
    aux_acc = 0;
    aux_conf  = zeros(2,2);
        
    for f=1:folds
        test_id = test(cv,f);
        train_id = training(cv,f);
        
        x_new = xdata(:,train_id);
        x_test = xdata(:,test_id);
        
        y_new = ydata(train_id);
        y_test = ydata(test_id);
        
%         Train classifier        
%         class_struct = fitcensemble(x_new',y_new); % default parameters yield best results
        class_struct = fitcensemble(x_new',y_new,'Learners',t); % 
%         class_struct = fitcensemble(x_new',y_new,'NumLearningCycles',100,...
%             'LearnRate',1,'Learners',t,'Method','LogitBoost');

        
        [pred_lab,pred_score] = predict(class_struct,x_test');
        [aux_fpr,aux_tpr,~,aux_auc] = perfcurve(y_test,pred_score(:,1),posclass);
%         aux_acc = aux_acc + sum(strcmp(pred_lab,y_test))./numel(y_test);
        aux_acc = aux_acc + sum(pred_lab & y_test)./numel(y_test);
        
        % Confusion matrix
        aux_conf = aux_conf + confusionmat(y_test,pred_lab);
        
        [~, ind] = unique(aux_fpr);
        aux_tpr2 = aux_tpr2 + interp1(aux_fpr(ind),aux_tpr(ind),fpr);
        
        aux_auc2 = aux_auc2+aux_auc;

    end
    tpr(:,i) = aux_tpr2./folds; %average interpolated tpr
    auc(i) = aux_auc2./folds;%average auc
    acc(i) =  aux_acc./folds; % average accuracy
    cm(:,:,i) = aux_conf./folds; % average confusion matrix  
    
end
