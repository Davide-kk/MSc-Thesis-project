%Al momento il modello UNEQ predice male perchè lo spazio delle PC in cui
%sono proiettati gli scores è formato da quelle PC che spiegano la maggior
%parte della varianza, la distinzione tra diversi contenuti di cagliata
%congelata in sti dati da come si  visto anche nelle PLS non è lungo la
%direzione di massima varianza. Quindi quello che manca in sto modello è un
%modo per selezionare le PC che contengono l'informazione che mi interessa.
load('WorkDataSets.mat'); %Loading all data workspace
category = zeros(height(All0Percent),1);
for k = 1:height(All0Percent)
    if contains(All0Percent.Sample{k,1},"_T0.001") == 1
        category(k,1) = 1;
    elseif contains(All0Percent.Sample{k,1},"_T1.001") == 1
        category(k,1) = 2;
    elseif contains(All0Percent.Sample{k,1},"_T2.001") == 1
        category(k,1) = 3;
    elseif contains(All0Percent.Sample{k,1},"MozzT2_") == 1
        category(k,1) = 4;
    elseif contains(All0Percent.Sample{k,1},"MozzT3_") == 1
        category(k,1) = 5;
    elseif contains(All0Percent.Sample{k,1},"MozzT4_") == 1
        category(k,1) = 6;
    else
        category(k,1) = 7;
    end
end
All0Percent = [All0Percent array2table(category)];
cats = {'TC_T0' 'TC_T1' 'TC_T2' 'OLD_T2' 'OLD_T3' 'OLD_T4' 'MARKET'};
All0Percent.category = categorical(category,[1 2 3 4 5 6 7],cats);
%%
%PCA for UNEQ 
%control whether there are outlier in the training set
All0Percent(43:end,:) = [];
curves = zeros(1129,height(All0Percent));
for k = 1:height(All0Percent)
    curves(:,k) = All0Percent.Data{k}; 
end
curves = curves';
curves = normalize(curves);
time = All0Percent.Time{1,1};
curves(:,1:5) = []; %strumental noisy data
time(1:5) = [];
curves(:,826:end) = []; %high noisy data
time(826:end) = [];

ncomp = 5;
[loadings,scores,eigenVal,tsquared,explained] = pca(curves,'NumComponents',ncomp);
normscores = scores./eigenVal(1:ncomp,1)'; %normalisation of PC scores by their Eigen Value

gscatter(normscores(:,4),normscores(:,5),All0Percent.category)
title("PCA");
xlabel("PC1 (93.4%)");
ylabel("PC2 (5.8%)");
%Calculate the critical T^2 value with 95% confidence
n = length(normscores);
p = 2;
F95 = 2.8387; %F value for alpha=0.05, df1 = p and df2 = 40 (nearest number to n-p)
F90 = 2.22609; %F value for alpha=0.10, df1 = p and df2 = 40 (nearest number to n-p)
t2crit95 = ((p*(n-1)*(n+1))/(n*(n-p)))*F95; %critical t squares value rapresents the surface
t2crit90 = ((p*(n-1)*(n+1))/(n*(n-p)))*F90; %in which the population is included
t2reduced = mahal(normscores(:,4:5),normscores(:,4:5)); %equal to the mahalanobis distance between not normalized scores

idxOut95 = zeros(n,1);
for k= 1:n
    if t2reduced(k,1)>t2crit95
       idxOut95(k,1) = 1;
    else
        idxOut95(k,1) = 0;
    end
end

idxOut90 = zeros(n,1);
for k= 1:n
    if t2reduced(k,1)>t2crit90
       idxOut90(k,1) = 1;
    else
        idxOut90(k,1) = 0;
    end
end
gscatter(1:n,t2reduced,All0Percent.category) %scatterplot of scores distance from population centroid
xlabel("Campioni");
ylabel("Hotelling's T^{2}");
%yline(t2crit95,'--',"T^{2}crit. (\alpha=0.05)");
yline(t2crit90,'--',"T^{2}crit. (\alpha=0.10)");
%Outliers removal from training set
curves(logical(idxOut90),:) = [];
scores(logical(idxOut90),:) = [];
normscores(logical(idxOut90),:) = [];
%%
%trying to get the most informative PCs by computing the MSPE

%%
%plottin loadings
subplot(3,2,1)
plot(time,loadings(:,1),'-');
subplot(3,2,2)
plot(time,loadings(:,2),'-');
subplot(3,2,3)
plot(time,loadings(:,3),'-');
subplot(3,2,4)
plot(time,loadings(:,4),'-');
subplot(3,2,5)
plot(time,loadings(:,5),'-');
subplot(3,2,6)
plot(time,loadings(:,6),'-');
%%
%Creating test set for UNEQ validation (checking outliers)
TC = AllTCSamples(25:end,:);
TC = removevars(TC,'Category');
OLD = ShelfLifeTable(19:end,:);
testset = [TC ; OLD];
testset = sortrows(testset,'FrozenCurdPercent','ascend');
explorer = zeros(1129,height(testset));
for k = 1:height(testset)
    explorer(:,k) = testset.Data{k}; 
end
testset.FrozenCurdPercent = categorical(testset.FrozenCurdPercent);
explorer = explorer';
explorer = normalize(explorer);
explorer(:,1:5) = [];
explorer(:,826:end) = [];

[testload,testscor,~,testT2,testexpl] = pca(explorer,"NumComponents",5);

h = scatter3(testscor(1:24,3),testscor(1:24,4),testscor(1:24,5),20,"green");
hold on
h1 = scatter3(testscor(25:42,3),testscor(25:42,4),testscor(25:42,5),20,"blue");
h2 = scatter3(testscor(43:66,3),testscor(43:66,4),testscor(43:66,5),20,"magenta");
h3 = scatter3(testscor(67:91,3),testscor(67:91,4),testscor(67:91,5),20,"black");
h4 = scatter3(testscor(92:end,3),testscor(92:end,4),testscor(92:end,5),20,"red");
hold off
xlabel("PC3 (93.52%)");
ylabel("PC4 (5.63%)");
zlabel("PC5 (0.60%)");
legend([h h1 h2 h3 h4],{'15' '16' '30' '50' '100'});
%%
%merge explorer with 4 0% samples which will be used for validation instead
%of training
idx0fc = [7 14 28 32];  %4 random 0% samples
curves(idx0fc,:) = [];   %remove these samples from the training set
All0Percent(logical(idxOut90),:) = []; %romeving outliers from training table
zerofc = All0Percent(idx0fc,:);  %getting raw data of 0% destined to validation
ValiSet = [zerofc ; testset];  %merging 0% with all the other categories

valid = zeros(1129,height(ValiSet));
for k = 1:height(ValiSet)
    valid(:,k) = ValiSet.Data{k}; %getting raw data for validation
end
valid = valid'; %transposing them
valid = normalize(valid);  %centering and scaling all together
%%
%UNEQ training
uneqscores = normscores;
[n2,p2] = size(uneqscores);
uneqmd = mahal(uneqscores,uneqscores);
criticalmd95 = ((p2*(n2-1)*(n2+1))/(n2*(n2-p2)))*2.45;
criticalmd90 = ((p2*(n2-1)*(n2+1))/(n2*(n2-p2)))*2;
%UNEQ test
valiscores = explorer*loadings;
valiscores = valiscores./eigenVal(1:ncomp,1)';
mdtest = mahal(valiscores,uneqscores);

gscatter(1:109,mdtest,testset.FrozenCurdPercent)
yline(criticalmd95);
scatter3(uneqscores(:,1),uneqscores(:,4),uneqscores(:,5),20,"blue")
hold on
scatter3(valiscores(:,1),valiscores(:,4), valiscores(:,5),20,'red')
hold off