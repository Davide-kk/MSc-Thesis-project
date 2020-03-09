load('WorkDataSets.mat') %Loading the workspace
%%
%Plot CPMG curves del TCset
for k = 1:6
    p1 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'r');
    hold on
end
for k = 7:14
    p2 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'g');
end

for k = 15:20
    p3 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'b');
end

for k = 21:28
    p4 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'k');
end
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
legend([p1 p2 p3 p4],{'0%' '15%' '30%' '50%'});
xlim([0 2330]);
hold off
%%
%plot of 0% OLD and TC samples
cla reset
subplot(1,2,1);
for k = 1:8
    p1 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'r');
    hold on
end
for k = 1:6
    p52 = plot(OldTable.Time{k,1},OldTable.Data{k,1},'m');
end
xlim([0 2330]);
ylim([0 60]);
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
legend([p1 p52], {'0% New' '0% Old'});
hold off
%plot of 15% OLD and TC samples
subplot(1,2,2);
for k = 7:14
    p2 = plot(TCTable.Time{k,1},TCTable.Data{k,1},'b');
    hold on
end
for k = 7:12
    p6 = plot(OldTable.Time{k,1},OldTable.Data{k,1},'c');
end
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
legend([p2 p6],{'15% New' '16% Old'});
xlim([0 2330]);
ylim([0 60]);
hold off
%%
%Plot CPMG curves del OldSet
for k = 1:6
    p5 = plot(OldTable.Time{k,1},OldTable.Data{k,1},'r');
    hold on
end
for k = 7:12
    p6 = plot(OldTable.Time{k,1},OldTable.Data{k,1},'g');
end

for k = 13:18
    p7 = plot(OldTable.Time{k,1},OldTable.Data{k,1},'b');
end
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
legend([p5 p6 p7],{'0%' '16%' '100%'});
xlim([0 2330]);
ylim([0 60]);
hold off
%%
x0 = [10 10 5 5 0 10 30 100 700];
% 4-exp regression on TC samples
TCfit = cell(1,height(TCTable));
for k=1:height(TCTable)
    time = TCTable.Time{1};
    data = TCTable.Data{k,1};
    ft = fittype('offset + abs(I1*exp(-time./t21)) + abs(I2*exp(-time./t22)) + (I3*exp(-time./t33)) + (I4*exp(-time./t44))','independent', 'time', 'dependent', 'data');
    TCfit{k} = fit(time, data, ft, 'StartPoint', x0,'Algorithm','Levenberg-Marquardt');
end

% 4-exp regression on OLD samples
Oldfit = cell(1,height(OldTable));
for k=1:height(OldTable)
    time = OldTable.Time{1};
    data = OldTable.Data{k,1};
    ft = fittype('offset + abs(I1*exp(-time./t21)) + abs(I2*exp(-time./t22)) + (I3*exp(-time./t33)) + (I4*exp(-time./t44))','independent', 'time', 'dependent', 'data');
    Oldfit{k} = fit(time, data, ft, 'StartPoint', x0,'Algorithm','Levenberg-Marquardt');
end
%%
%Bontà delle regressioni (Goodness of fit)
%TC dataset
TCgof = cell(1,height(TCTable));
for k=1:height(TCTable)
    time = TCTable.Time{1};
    data = TCTable.Data{k,1};
    ft = fittype('offset + abs(I1*exp(-time./t21)) + abs(I2*exp(-time./t22)) + (I3*exp(-time./t33)) + (I4*exp(-time./t44))','independent', 'time', 'dependent', 'data');
    [~,TCgof{k}] = fit(time, data, ft, 'StartPoint', x0,'Algorithm','Levenberg-Marquardt');
end

TCrmse = zeros(height(TCTable),1);
for k = 1:height(TCTable)
TCrmse(k,1) = TCgof{1,k}.rmse;
end

TCsse = zeros(height(TCTable),1);
for k = 1:height(TCTable)
TCsse(k,1) = TCgof{1,k}.sse;
end

%Old Dataset
OLDgof = cell(1,height(OldTable));
for k=1:height(OldTable)
    time = OldTable.Time{1};
    data = OldTable.Data{k,1};
    ft = fittype('offset + abs(I1*exp(-time./t21)) + abs(I2*exp(-time./t22)) + (I3*exp(-time./t33)) + (I4*exp(-time./t44))','independent', 'time', 'dependent', 'data');
    [~,OLDgof{k}] = fit(time, data, ft, 'StartPoint', x0,'Algorithm','Levenberg-Marquardt');
end

OLDrmse = zeros(height(OldTable),1);
for k = 1:height(OldTable)
OLDrmse(k,1) = OLDgof{1,k}.rmse;
end

OLDsse = zeros(height(OldTable),1);
for k = 1:height(OldTable)
OLDsse(k,1) = OLDgof{1,k}.sse;
end
%%
% Average statistics for goodness of fit of the regressions
AvgRmseTC = mean(TCrmse)
StdRmseTC = std(TCrmse)
AvgSseTC = mean(TCsse)
StdSseTC = std(TCsse)
AvgRmseOLD = mean(OLDrmse)
StdRmseOLD = std(OLDrmse)
AvgSseOLD = mean(OLDsse)
StdSseOLD = std(OLDsse)
%%
%Grafico dei residui dei TC
cla reset
subplot(1,2,1)
plot(TCfit{14},TCTable.Time{14,1},TCTable.Data{14,1},'residuals')
xlim([0 2330]);
ylim([-0.3 0.31]);
ylabel("Residui");
xlabel("Tempo (ms)");
%Grafico dei residui degli OLD
subplot(1,2,2)
plot(Oldfit{4},OldTable.Time{4,1},OldTable.Data{4,1},'residuals')
xlim([0 2330]);
ylim([-0.3 0.31]);
ylabel("Residui");
xlabel("Tempo (ms)");
%%
%Confidence and Prediction bounds
cla reset
PredBounds = predint(TCfit{1},TCTable.Time{1,1},0.99,'observation','on');
plot(TCfit{1}, TCTable.Time{1,1},TCTable.Data{1,1})
hold on
plot(TCTable.Time{1,1},PredBounds,'m--')
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
xlim([0 2330]);
hold off
%magnification of previous plot to insert in it 
cla reset
plot(TCfit{1}, TCTable.Time{1,1},TCTable.Data{1,1})
hold on
plot(TCTable.Time{1,1},PredBounds,'m--')
xlabel("Tempo (ms)");
ylabel("Intensità (A.U.)");
xlim([350 500]);
ylim([3 6.5]);
hold off
%%
% Tables of coefficents per i TC samples
TCCoefTable = zeros(height(TCTable),9);
for k = 1:length(TCCoefTable)
    TCCoefTable(k,:) = coeffvalues(TCfit{k});
end
SampleNames = TCTable.Sample;                       
TCCoefTable = array2table(TCCoefTable,"VariableNames",{'Ijz' 'Idw' 'Ifat' 'Icw' 'offset' 'T2jz' 'T2dw' 'T2fat' 'T2cw'},"RowNames",SampleNames);

FrozenCurd = zeros(height(TCTable),1);
for k = 1:length(FrozenCurd)
    if contains(TCCoefTable.Row(k),"TC1") == 1
        FrozenCurd(k,1) = 0;
    elseif contains(TCCoefTable.Row(k),"TC2") == 1
         FrozenCurd(k,1) = 15;
    elseif contains(TCCoefTable.Row(k),"TC3") == 1
        FrozenCurd(k,1) = 30;
    else 
        FrozenCurd(k,1) = 50;
    end
end

TCCoefTable = [TCCoefTable array2table(FrozenCurd)];
TCCoefTable = removevars(TCCoefTable, 'offset');

ITot = zeros(height(TCCoefTable),1);   
for k = 1:height(TCCoefTable)          
    ITot(k,1) = sum(TCCoefTable{k,1:4},2);   % I(%) = (I(a.u.)/Itot)*100
end
TCCoefTable(:,1:4) = array2table((TCCoefTable{:,1:4}./ITot)*100); %Intensities normalisation on Itot 'cos samples ain't weighted
%%
% Tables of coefficents of OLD samples regressions
OLDCoefTable = zeros(height(OldTable),9);
for k = 1:length(OLDCoefTable)
    OLDCoefTable(k,:) = coeffvalues(Oldfit{k});
end
SampleNames = OldTable.Sample;                       
OLDCoefTable = array2table(OLDCoefTable,"VariableNames",{'Ijz' 'Idw' 'Ifat' 'Icw' 'offset' 'T2jz' 'T2dw' 'T2fat' 'T2cw'},"RowNames",SampleNames);

FrozenCurd2 = zeros(height(OldTable),1);
for k = 1:length(FrozenCurd2)
    if contains(OLDCoefTable.Row(k),"_0F_") == 1
        FrozenCurd2(k,1) = 0;
    elseif contains(OLDCoefTable.Row(k),"_16_") == 1
         FrozenCurd2(k,1) = 15;
    else 
        FrozenCurd2(k,1) = 100;
    end
end

OLDCoefTable = [OLDCoefTable array2table(FrozenCurd2)];
OLDCoefTable = removevars(OLDCoefTable, 'offset');

ITot = zeros(height(OLDCoefTable),1);   
for k = 1:height(OLDCoefTable)          
    ITot(k,1) = sum(OLDCoefTable{k,1:4},2);   % I(%) = (I(a.u.)/Itot)*100
end
OLDCoefTable(:,1:4) = array2table((OLDCoefTable{:,1:4}./ITot)*100); %Intensities normalisation on Itot 'cos samples ain't weighted
%%
%Boxplots of coefficents against categories
%TC Intensities
cla reset
subplot(2,2,1);
boxplot(TCCoefTable.Ijz,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)","FontSize",14);
title("I_{jz}","FontSize",18);
subplot(2,2,2);
boxplot(TCCoefTable.Idw,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)","FontSize",14);
title("I_{dw}","FontSize",18);
subplot(2,2,3);
boxplot(TCCoefTable.Ifat,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)","FontSize",14);
title("I_{fat}","FontSize",18);
subplot(2,2,4);
boxplot(TCCoefTable.Icw,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)","FontSize",14);
title("I_{cw}","FontSize",18);
%%
%TC T2 times
cla reset
subplot(2,2,1);
boxplot(TCCoefTable.T2jz,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{jz}","FontSize",18);
subplot(2,2,2);
boxplot(TCCoefTable.T2dw,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{dw}","FontSize",18);
subplot(2,2,3);
boxplot(TCCoefTable.T2fat,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{fat}","FontSize",18);
subplot(2,2,4);
boxplot(TCCoefTable.T2cw,TCCoefTable.FrozenCurd);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{cw}","FontSize",18);
%%
%OLD Intensities
cla reset
subplot(2,2,1);
boxplot(OLDCoefTable.Ijz,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)");
title("I_{jz}");
subplot(2,2,2);
boxplot(OLDCoefTable.Idw,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)");
title("I_{dw}");
subplot(2,2,3);
boxplot(OLDCoefTable.Ifat,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)");
title("I_{fat}");
subplot(2,2,4);
boxplot(OLDCoefTable.Icw,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Intensità (A.U.)");
title("I_{cw}");
%%
% OLD T2 times
cla reset
subplot(2,2,1);
boxplot(OLDCoefTable.T2jz,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{jz}","FontSize",18);
subplot(2,2,2);
boxplot(OLDCoefTable.T2dw,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{dw}","FontSize",18);
subplot(2,2,3);
boxplot(OLDCoefTable.T2fat,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{fat}","FontSize",18);
subplot(2,2,4);
boxplot(OLDCoefTable.T2cw,OLDCoefTable.FrozenCurd2);
xlabel("Cagliata Congelata (%)","FontSize",14);
ylabel("Tempo (ms)","FontSize",14);
title("T2_{cw}","FontSize",18);
%%
%OUTLIERS REMOVAL
pattern = ["TC1_DX_2_" "TC3_DX_2_"];
idx = contains(TCTable.Sample,pattern);
TCCoefTable(idx,:) = [];
TCTable(idx,:) = [];
%%
%Mann-Whitney U-test significant differences for coefficents between
%categories
%TCSAMPLES
bau = [0 15 30 50];  %Test's results interpretation:
UtestIjz = zeros(4,4); %'1' means rejection of null hypothesis (h0) 
for k = 1:4           %with 95% confidence = significant diff between distributions
    for i = 1:4       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,UtestIjz(i,k)] = ranksum(TCCoefTable.Ijz(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.Ijz(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for Idw
UtestIdw = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestIdw(i,k)] = ranksum(TCCoefTable.Idw(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.Idw(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for Ifat
UtestIfat = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestIfat(i,k)] = ranksum(TCCoefTable.Ifat(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.Ifat(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for Icw
UtestIcw = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestIcw(i,k)] = ranksum(TCCoefTable.Icw(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.Icw(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for T2jz
UtestT2jz = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestT2jz(i,k)] = ranksum(TCCoefTable.T2jz(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.T2jz(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for T2dw
UtestT2dw = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestT2dw(i,k)] = ranksum(TCCoefTable.T2dw(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.T2dw(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for T2fat
UtestT2fat = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestT2fat(i,k)] = ranksum(TCCoefTable.T2fat(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.T2fat(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%test for T2cw
UtestT2cw = zeros(4,4);
for k = 1:4
    for i = 1:4
[~,UtestT2cw(i,k)] = ranksum(TCCoefTable.T2cw(ismember(TCCoefTable.FrozenCurd,bau(k))),TCCoefTable.T2cw(ismember(TCCoefTable.FrozenCurd,bau(i))));
    end
end
%%
%Mann-Whitney U-test significant differences for coefficents between
%categories
%OLD SAMPLES
bau = [0 15 100];  %Test's results interpretation:
OLDUtestIjz = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestIjz(i,k)] = ranksum(OLDCoefTable.Ijz(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.Ijz(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for Idw
OLDUtestIdw = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestIdw(i,k)] = ranksum(OLDCoefTable.Idw(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.Idw(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for Ifat
OLDUtestIfat = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestIfat(i,k)] = ranksum(OLDCoefTable.Ifat(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.Ifat(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for Icw
OLDUtestIcw = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestIcw(i,k)] = ranksum(OLDCoefTable.Icw(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.Icw(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for T2jz
OLDUtestT2jz = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestT2jz(i,k)] = ranksum(OLDCoefTable.T2jz(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.T2jz(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for T2dw
OLDUtestT2dw = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestT2dw(i,k)] = ranksum(OLDCoefTable.T2dw(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.T2dw(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for T2fat
OLDUtestT2fat = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestT2fat(i,k)] = ranksum(OLDCoefTable.T2fat(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.T2fat(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%diff for T2cw
OLDUtestT2cw = zeros(3,3); %'1' means rejection of null hypothesis (h0) 
for k = 1:3           %with 95% confidence = significant diff between distributions
    for i = 1:3       %'0' means h0 is confirmed (95% prob) aka no significant difference
[~,OLDUtestT2cw(i,k)] = ranksum(OLDCoefTable.T2cw(ismember(OLDCoefTable.FrozenCurd2,bau(k))),OLDCoefTable.T2cw(ismember(OLDCoefTable.FrozenCurd2,bau(i))));
    end
end
%%
%Mann-Whitney between 0 and 15% of the OLD and TC samples
p0 = zeros(8,1);
p0(1,1) =  ranksum(TCCoefTable.Ijz(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.Ijz(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(2,1) =  ranksum(TCCoefTable.Idw(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.Idw(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(3,1) =  ranksum(TCCoefTable.Ifat(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.Ifat(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(4,1) =  ranksum(TCCoefTable.Icw(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.Icw(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(5,1) =  ranksum(TCCoefTable.T2jz(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.T2jz(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(6,1) =  ranksum(TCCoefTable.T2dw(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.T2dw(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(7,1) =  ranksum(TCCoefTable.T2fat(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.T2fat(ismember(OLDCoefTable.FrozenCurd2,0)));
p0(8,1) =  ranksum(TCCoefTable.T2cw(ismember(TCCoefTable.FrozenCurd,0)),OLDCoefTable.T2cw(ismember(OLDCoefTable.FrozenCurd2,0)));

p15 = zeros(8,1);
p15(1,1) =  ranksum(TCCoefTable.Ijz(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.Ijz(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(2,1) =  ranksum(TCCoefTable.Idw(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.Idw(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(3,1) =  ranksum(TCCoefTable.Ifat(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.Ifat(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(4,1) =  ranksum(TCCoefTable.Icw(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.Icw(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(5,1) =  ranksum(TCCoefTable.T2jz(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.T2jz(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(6,1) =  ranksum(TCCoefTable.T2dw(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.T2dw(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(7,1) =  ranksum(TCCoefTable.T2fat(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.T2fat(ismember(OLDCoefTable.FrozenCurd2,15)));
p15(8,1) =  ranksum(TCCoefTable.T2cw(ismember(TCCoefTable.FrozenCurd,15)),OLDCoefTable.T2cw(ismember(OLDCoefTable.FrozenCurd2,15)));
%%
%mean and std of fit coefficents
%TC SAMPLE
a = mean(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,0)),:});
b = mean(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,15)),:});
c = mean(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,30)),:});
d = mean(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,50)),:});
avgTCcoef = [a ; b; c; d]
a = std(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,0)),:});
b = std(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,15)),:});
c = std(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,30)),:});
d = std(TCCoefTable{(ismember(TCCoefTable.FrozenCurd,50)),:});
stdTCcoef = [a ; b; c; d]
%%
%OLD SAMPLES
a = mean(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,0)),:});
b = mean(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,15)),:});
c = mean(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,100)),:});
avgOLDcoeff = [a;b;c]
a = std(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,0)),:});
b = std(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,15)),:});
c = std(OLDCoefTable{(ismember(OLDCoefTable.FrozenCurd2,100)),:});
stdOLDcoeff = [a;b;c]