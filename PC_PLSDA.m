%% 
% *PC-PLSDA SCRIPT, LOAD THE WORKDATASET WORKSPACE FROM DIRECTORY*

%vector to be used later to determine test order (necessary for both
%training and validation sets)
%contvars = TCTable.FrozenCurdPercent;
categories = cellstr(TCTable.categories);
%%
%RIMOZIONE DEGLI OUTLIERS
pattern = ["TC1_DX_2_" "TC3_DX_2_"];
idx = contains(TCTable.Sample,pattern);
TCTable(idx,:) = [];
categories(idx,:) = [];
%contvars(idx,:) = [];
%%
% Principal Components Analysis of the dataset used for pls training
curves = zeros(1129,height(TCTable));
for k = 1:height(TCTable)
    curves(:,k) = TCTable.Data{k}; 
end
curves = curves';
curves = normalize(curves);
[n,m] = size(curves);

[coeff,score,~,~,explained] = pca(curves);

%PC ranking based on discriminant power
pvalues = zeros(length(score)-1,1);
for k = 1:length(score)-1
    [p,~,~] = anova1(score(:,k)',TCTable{:,"categories"}');
    pvalues(k,1) = p;
end

pvalues = array2table(pvalues);
PCorder = cell(length(score)-1,1);
form = 'PC%d';
for k = 1:length(score)-1
    PCorder{k} = sprintf(form,k);
end
PCorder = cell2table(PCorder);
PCorder = [PCorder pvalues];
PCorder = sortrows(PCorder,'pvalues','ascend'); %sorted table of PC's p-values 

pcorderidx = zeros(height(PCorder),1); %row index to sort scores matrix in descending discrimination power
for k = 1:height(PCorder)
    pcorderidx(k,1) = sscanf(PCorder.PCorder{k},strcat("PC","%d"));
end
%%
%PCA scores plot
colourmap = zeros(length(TCTable.Sample),3); %colourmap
for k = 1:length(TCTable.Sample)
    if ismember(TCTable.categories(k),"TC1") == 1
        colourmap(k,:) = [1 0 0];
    elseif ismember(TCTable.categories(k),"TC2") == 1
          colourmap(k,:) = [0 1 0];
    elseif ismember(TCTable.categories(k),"TC3") == 1
          colourmap(k,:) = [0 0 1];
    else 
          colourmap(k,:) = [0 0 0];
    end
end

scatter3(score(:,5),score(:,1),score(:,9),20,colourmap)
title("PCA");
xlabel("PC5 (0.84%)");
ylabel("PC1 (0.02%)");
zlabel("PC9 (0.01%)");
%%
%plot of the PC loadings vs. the relaxation time, good to look at it
%anytime the dataset used in the PCA is modified
cla reset
plot(TCTable.Time{1},coeff(:,[1 3 4 2 5]),'-');
xlabel('Time (ms)');
ylabel('PC coeff');
legend("PC1","PC3","PC4","PC2","PC5","location","bestoutside");

yyaxis right
plot(TCTable.Time{1},TCTable.Data{1},'-');
ylabel('Intensity (A.U.)');
%%
%Partial Least Square Discriminant Analysis
Y = zeros(length(TCTable.Sample),4); %Assign values of froze curd percent to categorical variables
for k = 1:length(TCTable.Sample)     %categorical variables are passed to the program as binary codes
    if ismember(TCTable.categories(k),"TC1") == 1   %the length of the binary code corresponds to the 
        Y(k,:) = [1 0 0 0];                         %number of categories
    elseif ismember(TCTable.categories(k),"TC2") == 1
            Y(k,:) = [0 1 0 0];
    elseif ismember(TCTable.categories(k),"TC3") == 1
            Y(k,:) = [0 0 1 0];
    else 
            Y(k,:) = [0 0 0 1];
    end
end

score2 = score(:, pcorderidx); %score matrix sorted by descending disciminant power of pc
cvpart = cvpartition(TCTable.categories, 'LeaveOut'); %dataset partition for cross-validation training

[XL,YL,XS,YS,BETA,PCTVAR,MSPE,stats] = plsregress(score2(:,[1 2 4 5 6 7 8]),Y,5,"cv",cvpart); %PLS
%remember THIS is the MODEL, the stuff in the following section is the crossvalidation of
%this model
DAmdl = fitcdiscr(XS(:,1:2),categories);
%%
cla reset                                 %percent of variance explained by the LVs showed in a plot
plot(1:5,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of LV');
ylabel('Percent Variance Explained in y');

[~,plsdacomp] = maxk(PCTVAR(2,:),3) ;      %combination of the LVs explaining the highest variance
percyvar = sum(PCTVAR(2,plsdacomp));

cla reset                   %MSEP plot vs. n° of LV to have a look on optimal LV number
plot(0:5,MSPE(2,:),'-or');   %repeating this plot with different PC number allows to optimize
xlabel('Number of LV');      %PC number used in the model as well
ylabel('Estimated Mean Squared Prediction Error');
%%
%CrossValidation of the PC_PLSDA training model
cvscores = zeros(height(TCTable),5);
for k = 1:length(cvscores)          %remove risk of overfitting by the calculation of the
    idxTrain = training(cvpart,k);  %score by multiplying the PCA score by the weigths matrix
    idxTest = test(cvpart,k);        %of the PLS model without it with this aproach is safe to obtain a 
    [~,~,~,~,bet,~,~,stat] = plsregress(score2(idxTrain,[1 2 4 5 6 7 8]),Y(idxTrain,:),5);%true prediction without risk
    cvscores(k,:) = score2(idxTest,[1 2 4 5 6 7 8])*stat.W; %of overparametrisation
end                                             
                                                                                                                                  
orderedY = cell(height(TCTable),1);                                                        
for k = 1:length(orderedY)
    orderedY(k) = categories(test(cvpart,k));
end

mdl = fitcdiscr(cvscores,orderedY,'CrossVal',"on","CVPartition",cvpart); %the scores calculated above are
label = kfoldPredict(mdl);               %then classified by a linear discriminant classifier, a model 
confmtx = confusionmat(orderedY,label,'Order',{'TC1' 'TC2' 'TC3' 'TC4'}); %which computes classification
confusionchart(orderedY,label)      %based on the possibility to spatially divide scores of different categories
NER = (confmtx(1,1)+confmtx(2,2)+confmtx(3,3)+confmtx(4,4))/length(TCTable.Sample);
ER = 1-NER
sensitivityTC1 = confmtx(1,1)/sum(confmtx(1,:),2)
sensitivityTC2 = confmtx(2,2)/sum(confmtx(2,:),2)
sensitivityTC3 = confmtx(3,3)/sum(confmtx(3,:),2)
sensitivityTC4 = confmtx(4,4)/sum(confmtx(4,:),2)
%%
%Validation of the model with a Validation Set
%Prediction of the validation set
curvesValidationSet = zeros(1129,height(OldTable)); %extraction of the curves from the validation set
for k = 1:height(OldTable)
    curvesValidationSet(:,k) = OldTable.Data{k}; 
end
curvesValidationSet = curvesValidationSet';
curvesValidationSet = normalize(curvesValidationSet);
curvesValidationSet(13:end,:) = [];
coeff2 = coeff(:,pcorderidx);
PCscoresValiSet = curvesValidationSet*coeff2; %calculation of the scores from the PCA model multiplying
                                                   %the curves with PC loadings
PLSscoresValiSet = PCscoresValiSet(:,[1 2 4 5 6 7 8])*stats.W;

labelValidation = predict(DAmdl,PLSscoresValiSet); %Something goes wrong, prediction of validation set is somehow 
%not correct and weird as well, as the PLS is defentley not overparametrized
%maybe something is messing with the Discriminant Classifier, another
%question is wether predicting sample classes on PLS predictors instead of responses is correct or not
%%
gscatter(cvscores(:,2), cvscores(:,5),categorical(orderedY))