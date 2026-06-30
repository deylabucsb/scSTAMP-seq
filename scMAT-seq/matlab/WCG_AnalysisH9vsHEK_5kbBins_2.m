clc
clear
minReadsNeeded = 10000;
PercentofCellsNeeded = 0.2;
%%
%read in data
a = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib3\Lib3_AllFaba\P39_5mC_Access_mRNA_Fresh_H9_Lib3_24_Cells_A_AllSeqRun_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
aa = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib4\Lib4_AllFaba\P39_5mC_Access_mRNA_Fresh_H9_Lib4_AllSeqRun_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

b = table2array(a(1:606218,2:97)); %remove Y chromosome
bb = table2array(aa(1:606218,2:97));
RownamesTemp = table2array(a(1:606218,1));


Celltypetemp = [b,bb];
TotperCellTypeH9 = sum(Celltypetemp);
KeepCellsType = find(TotperCellTypeH9>minReadsNeeded);
CellTypetempH9 = Celltypetemp(:,KeepCellsType);
% BinaryCellTypetempH9 = CellTypetempH9 > 0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts
% BinaryCellTypetempH9 = CellTypetemp2;

% KeepRowsH9 = find(sum(BinaryCellTypetempH9,2) > ceil(PercentofCellsNeeded*length(KeepCellsType)));

clear a aa b bb Celltypetemp KeepCellsType
%%
a = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\Plate4_FrozenHEK293T_mRNA_5mC_Accessibilty\AllCombined_5mC\P4_AllCondition_CellNumberCorrected_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
aa = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\Plate3_McviPI\P3_McvipI_1C-3_HEK293T_48Cells2-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

b = table2array(a(1:606218,2:97));
bb = table2array(aa(1:606218,2:97));
RownamesTemp = table2array(a(1:606218,1));


Buffertable = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\Plate4_FrozenHEK293T_mRNA_5mC_Accessibilty\AllCombined_5mC\P4_ConditionUsed_CellNumberCorrected.txt','headerlines',0,'ReadVariableNames',false);
Buffer = table2array(Buffertable);
KeepBufferRows1 = find(contains(Buffer,'FSB-GCB.'));
KeepBufferRows2 = find(contains(Buffer,'200nM'));
% KeepBufferRows2 = find(contains(Buffer,'GCB.')); %200nM

KeepBufferRows = intersect(KeepBufferRows1,KeepBufferRows2);

% Celltypetemp = [b,bb];
Celltypetemp = [b(:,KeepBufferRows), bb];

BufferHEKtemp(1:size(Celltypetemp,2),1) = "FSB-GCB.125nM";
BufferHEKtemp(1:size(Buffer(KeepBufferRows),1),1) = Buffer(KeepBufferRows);

TotperCellTypeHEK = sum(Celltypetemp);
KeepCellsType = find(TotperCellTypeHEK>minReadsNeeded);
CellTypetempHEK = Celltypetemp(:,KeepCellsType);
BinaryCellTypetempHEK = CellTypetempHEK > 0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts
% BinaryCellTypetempHEK = CellTypetemp2;

BufferHEK = BufferHEKtemp(KeepCellsType);

% KeepRowsHEK = find(sum(BinaryCellTypetempHEK,2) > ceil(PercentofCellsNeeded*length(KeepCellsType)));

clear a aa b bb Celltypetemp KeepCellsType CellTypetemp2

%%
%%Lets remove bins that are clearly way to high (i.e. blacklist)
combinedtemp = [CellTypetempH9,CellTypetempHEK];
KeepNonBlacklist = find(sum(combinedtemp,2)<size(combinedtemp,2)*3);
combinedtempNonBL = combinedtemp(KeepNonBlacklist,:);
RownamesTemp2 = RownamesTemp(KeepNonBlacklist,:);
CountPerCelltemp2 = sum(combinedtempNonBL,1);

combinedtemp2 = combinedtempNonBL >0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts
CountPerCell = sum(combinedtemp);

%%
%Using the Cell type information we know
BinaryCellTypetempH9 = combinedtemp2(:,1:size(CellTypetempH9,2));
BinaryCellTypetempHEK = combinedtemp2(:,1+size(CellTypetempH9,2):end);

pcrtilekeep = 98;
KeepRowsH9 = find(sum(BinaryCellTypetempH9,2) > prctile(sum(BinaryCellTypetempH9,2),pcrtilekeep));
KeepRowsHEK = find(sum(BinaryCellTypetempHEK,2) > prctile(sum(BinaryCellTypetempHEK,2),pcrtilekeep));
KeepRows = unique([KeepRowsH9;KeepRowsHEK]);

AllBinaryTemp = [BinaryCellTypetempH9,BinaryCellTypetempHEK];

AllBinary = double(AllBinaryTemp(KeepRows,:));
AllRownames = RownamesTemp2(KeepRows,:);


% BinaryH9 = BinaryCellTypetempH9(KeepRows,:);
% H9Rownames = RownamesTemp(KeepRows,:);
% 
% BinaryHEK = BinaryCellTypetempHEK(KeepRows,:);
% HEKRownames = RownamesTemp(KeepRows,:);

%%
%Which Bins are used
BinarySTD = std(combinedtempNonBL,0,2);
BinaryMean = mean(combinedtempNonBL,2);
BinaryCV = BinarySTD./BinaryMean;

sz1=20;
sz2=10;
figure
scatter(BinaryMean,log10(BinaryCV),sz1,'sk','filled')
hold on
scatter(BinaryMean(KeepRowsH9),log10(BinaryCV(KeepRowsH9))+1,sz2,'or')
hold on
scatter(BinaryMean(KeepRowsHEK),log10(BinaryCV(KeepRowsHEK))+2,sz2,'ob')
xlabel("Mean Bin Val")
ylabel("Log10 CV")

%%
%addpath('C:\Users\alexc\Box Sync\Research\mRNA_Mapping\mRNA_Analysis\files')

figure;
M = AllBinary;

C = corr(M);

D = pdist(C);
Z = linkage(D,'ward');
[h,t,perm] = dendrogram(Z,0);

clear MOrg
for i = 1:1:length(perm)
    MOrg(:,i) = M(:,perm(i));
end
COrg = corr(MOrg);

% figure;
% imagesc(C);colorbar;

figure;
imagesc(COrg);
colormap jet; 
cb = colorbar;
set(gca,'XTick',1:length(perm),'XTickLabel',perm,'YTick',1:length(perm),'YTickLabel',perm,'FontSize',9,'FontName','Arial');
cb.Label.String = 'Pearson {\itr}';
cb.FontSize = 16;
xlabel('Single Cell','FontSize',18,'FontName','Arial');
ylabel('Single Cell','FontSize',18,'FontName','Arial');

%%
[PCAcoeff,score,latent,tsquared,explained,mu] = pca(AllBinary);
clear CellTypes
CellTypes(1:size(BinaryCellTypetempH9,2),1) = 1; % 1 is H9
CellTypes(size(BinaryCellTypetempH9,2)+1:size(AllBinary,2)) = 2; % 2 is HEK

BuffersUsedtemp(1:size(BinaryCellTypetempH9,2),1) = "FSB-GCB.200nM";
BuffersUsed = [BuffersUsedtemp;BufferHEK];

PCAcoeffMeanH9 = mean(PCAcoeff(1:size(CellTypetempH9,2),:));
PCAcoeffMeanHEK = mean(PCAcoeff(1+size(CellTypetempH9,2):end,:));
PCAcoeffMean = [PCAcoeffMeanH9;PCAcoeffMeanHEK];

figure
gscatter(PCAcoeff(:,1),PCAcoeff(:,2),CellTypes)
hold on
scatter(PCAcoeffMean(:,1),PCAcoeffMean(:,2),'ok','filled')
hold on
xlabel('PCA 1')
ylabel('PCA 2')

figure
gscatter(PCAcoeff(:,1),PCAcoeff(:,3),CellTypes)
xlabel('PCA 1')
ylabel('PCA 3')

figure
gscatter(PCAcoeff(:,1),PCAcoeff(:,2),BuffersUsed)
xlabel('PCA 1')
ylabel('PCA 2')

figure
bar(explained(1:6))

%%
%Correlations between PCA and the sequencing depth per cell (Thinking that
%PCA 1 is highly correlated to the sequencing depth
rng(100);
rndPermCells = randperm(size(CountPerCell,2));
RandomizedCountPerCell = CountPerCell(rndPermCells);
for i=1:1:20
   PCA_To_Depth(i) = corr(PCAcoeff(:,i),CountPerCell','Type','Pearson');
   PCA_To_DepthRand(i) = corr(PCAcoeff(:,i),RandomizedCountPerCell','Type','Pearson');
%    PCA_To_Depth(i) = corr(PCAcoeff(:,i),CountPerCell','Type','Spearman');
   PCA_To_DepthLookedAt(i) = i;
   Explained_To_DepthLookedAt(i) = explained(i);
end

figure
scatter(PCA_To_DepthLookedAt,PCA_To_Depth,'ok','filled')
% hold on
% scatter(PCA_To_DepthLookedAt,PCA_To_DepthRand,'ob','filled')
hold on
scatter(PCA_To_DepthLookedAt,abs(PCAcoeffMean(1,1:20)-PCAcoeffMean(2,1:20)),'or','filled')
legend("r PCA to Depth","Mean Dist between cell types")
% print(gcf, 'WCG_H9_HEK_PCA_vs_Counts_CurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')


figure
scatter(PCAcoeff(:,2),CountPerCell','ok','filled')
% hold on
% plot([-0.2,0.2],[abc(2)+abc(1)*-0.2,abc(2)+abc(1)*0.2])
xlabel('PCA 2')
ylabel('Detected Marks')

%%
%Cluster on PCA 1 based on elbow plot and corr with sequencing depth
figure
D2 = pdist(PCAcoeff(:,1));
Z2 = linkage(D2,'ward');
[h2,t2,perm2] = dendrogram(Z2,0);
ClustNumberDendo2temp = cluster(Z2,'MaxClust',2);

%Rearragne numbers to match cell types numbers
clear ClustNumberDendo2
ClustNumberDendo2(find(ClustNumberDendo2temp ==1),1) = 2;
ClustNumberDendo2(find(ClustNumberDendo2temp ==2),1) = 1;

%See if cell type and cluster asignment match
CorrectAsignment = CellTypes == ClustNumberDendo2;
CorrectAsignment2 = find(CorrectAsignment==0);

%%
figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial')
gscat1 = gscatter(PCAcoeff(:,1),PCAcoeff(:,2),CellTypes,[165/255,25/255,165/255;1,120/255,0/255])
hold on
gscatter(PCAcoeff(CorrectAsignment2,1),PCAcoeff(CorrectAsignment2,2),CorrectAsignment(CorrectAsignment2),'k','o')
xlabel('PCA 1')
ylabel('PCA 2')
title("Endogneous CpG Methylation")
% ylim([0,0.225])
% yticks([0.1,0.2,0.3,0.4])
xlim([-0.2 0.25])
%xticks([0,10,20,30,510,520])
leg1 =legend(gscat1,{'H9','HEK293T'})
set(leg1,'Position',[0.75, 0.20, .15, .15])
box off
grid off
legend boxoff
set(gcf,'color','w')


% print(gcf, 'WCG_H9_HEK_CurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')
% saveas(gcf,'WCG_H9_HEK_CurratedCellType_5kbBins.pdf')
%%

%%---%%
%Without using cell type information
KeepRowsNonCellType = find(sum(combinedtemp2,2) > prctile(sum(combinedtemp2,2),pcrtilekeep));

AllBinary2 = double(combinedtemp2(KeepRowsNonCellType,:));
AllRownames2 = RownamesTemp2(KeepRowsNonCellType,:);


%%
%addpath('C:\Users\alexc\Box Sync\Research\mRNA_Mapping\mRNA_Analysis\files')

figure;
M = AllBinary2;

C = corr(M);

D = pdist(C);
Z = linkage(D,'ward');
[h,t,perm] = dendrogram(Z,0);

clear MOrg
for i = 1:1:length(perm)
    MOrg(:,i) = M(:,perm(i));
end
COrg = corr(MOrg);

% figure;
% imagesc(C);colorbar;

figure;
imagesc(COrg);
colormap jet; 
cb = colorbar;
set(gca,'XTick',1:length(perm),'XTickLabel',perm,'YTick',1:length(perm),'YTickLabel',perm,'FontSize',9,'FontName','Arial');
cb.Label.String = 'Pearson {\itr}';
cb.FontSize = 16;
xlabel('Single Cell','FontSize',18,'FontName','Arial');
ylabel('Single Cell','FontSize',18,'FontName','Arial');

%%
[PCAcoeff2,score2,latent2,tsquared2,explained2,mu2] = pca(AllBinary2);
CellTypes(1:size(BinaryCellTypetempH9,2),1) = 1; % 1 is H9
CellTypes(size(BinaryCellTypetempH9,2)+1:size(AllBinary2,2)) = 2; % 2 is HEK

PCAcoeffMeanH92 = mean(PCAcoeff2(1:size(CellTypetempH9,2),:));
PCAcoeffMeanHEK2 = mean(PCAcoeff2(1+size(CellTypetempH9,2):end,:));
PCAcoeffMean2 = [PCAcoeffMeanH92;PCAcoeffMeanHEK2];

figure
gscatter(PCAcoeff2(:,1),PCAcoeff2(:,2),CellTypes)
hold on
scatter(PCAcoeffMean2(:,1),PCAcoeffMean2(:,2),'ok','filled')
xlabel('PCA 1')
ylabel('PCA 2')

figure
bar(explained2(1:6))
%%
%Correlations between PCA and the sequencing depth per cell (Thinking that
%PCA 1 is highly correlated to the sequencing depth
for i=1:1:20
   PCA_To_Depth2(i) = corr(PCAcoeff2(:,i),CountPerCell','Type','Pearson');
%    PCA_To_Depth2(i) = corr(PCAcoeff(:,i),CountPerCell','Type','Spearman');
   PCA_To_DepthLookedAt2(i) = i;
   Explained_To_DepthLookedAt2(i) = explained2(i);
end
figure
scatter(PCA_To_DepthLookedAt2,PCA_To_Depth2,'ok','filled')
hold on
scatter(PCA_To_DepthLookedAt2,abs(PCAcoeffMean2(1,1:20)-PCAcoeffMean2(2,1:20)),'or','filled')
legend("r PCA to Depth","Mean Dist between cell types")
% print(gcf, 'WCG_H9_HEK_PCA_vs_Counts_nonCurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

figure
scatter(PCAcoeff2(:,1),CountPerCell,'ok','filled')
xlabel('PCA 1')
ylabel('Detected Marks')

%%
%Cluster on PCA 1& 2 based on elbow plot and corr with sequencing depth
figure
D3 = pdist(PCAcoeff2(:,1:2));
Z3 = linkage(D3,'ward');
[h3,t3,perm3] = dendrogram(Z3,0);
ClustNumberDendo3temp = cluster(Z3,'MaxClust',2);

%Rearragne numbers to match cell types numbers
clear ClustNumberDendo3
ClustNumberDendo3(find(ClustNumberDendo3temp ==1),1) = 1;
ClustNumberDendo3(find(ClustNumberDendo3temp ==2),1) = 2;


%See if cell type and cluster asignment match
CorrectAsignmentnon = CellTypes == ClustNumberDendo3;
CorrectAsignmentnon2 = find(CorrectAsignmentnon==0);

%%
figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Ariel','DefaultAxesFontName','Ariel')
gscat1 = gscatter(PCAcoeff2(:,1),PCAcoeff2(:,2),CellTypes)
hold on
gscatter(PCAcoeff2(CorrectAsignmentnon2,1),PCAcoeff2(CorrectAsignmentnon2,2),CorrectAsignmentnon(CorrectAsignmentnon2),'k','o')
xlabel('PCA 1')
ylabel('PCA 2')
ylim([-0.13,0.24])
yticks([-0.1,0.0,0.1,0.2])
% xlim([-0.15 0.18])
%xticks([0,10,20,30,510,520])
leg1 =legend(gscat1,{'H9','HEK293T'},'location','north')
set(leg1,'Position',[0.46, 0.78, .15, .15])
box off
grid off
legend boxoff
set(gcf,'color','w')


% print(gcf, 'WCG_H9_HEK_NonCurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')


%%

%%---%%
%Miss assigning cell type information then trying it
RealCellTypes = CellTypes(rndPermCells);
RealCountPerCell = CountPerCell(rndPermCells);

BinaryCellTypetempH92 = combinedtemp2(:,rndPermCells(1:size(CellTypetempH9,2)));
BinaryCellTypetempHEK2 = combinedtemp2(:,rndPermCells(1+size(CellTypetempH9,2):end));

pcrtilekeep = 98;
KeepRowsH92 = find(sum(BinaryCellTypetempH92,2) > prctile(sum(BinaryCellTypetempH92,2),pcrtilekeep));
KeepRowsHEK2 = find(sum(BinaryCellTypetempHEK2,2) > prctile(sum(BinaryCellTypetempHEK2,2),pcrtilekeep));
KeepRows2 = unique([KeepRowsH92;KeepRowsHEK2]);

AllBinaryTemp2 = [BinaryCellTypetempH92,BinaryCellTypetempHEK2];

AllBinary3 = double(AllBinaryTemp2(KeepRows2,:));
AllRownames3 = RownamesTemp2(KeepRows2,:);


%%
%addpath('C:\Users\alexc\Box Sync\Research\mRNA_Mapping\mRNA_Analysis\files')

figure;
M = AllBinary3;

C = corr(M);

D = pdist(C);
Z = linkage(D,'ward');
[h,t,perm] = dendrogram(Z,0);

clear MOrg
for i = 1:1:length(perm)
    MOrg(:,i) = M(:,perm(i));
end
COrg = corr(MOrg);

% figure;
% imagesc(C);colorbar;

figure;
imagesc(COrg);
colormap jet; 
cb = colorbar;
set(gca,'XTick',1:length(perm),'XTickLabel',perm,'YTick',1:length(perm),'YTickLabel',perm,'FontSize',9,'FontName','Arial');
cb.Label.String = 'Pearson {\itr}';
cb.FontSize = 16;
xlabel('Single Cell','FontSize',18,'FontName','Arial');
ylabel('Single Cell','FontSize',18,'FontName','Arial');

%%
[PCAcoeff3,score3,latent3,tsquared3,explained3,mu3] = pca(AllBinary3);


PCAcoeffMeanH93 = mean(PCAcoeff3(find(RealCellTypes==1),:));
PCAcoeffMeanHEK3 = mean(PCAcoeff3(find(RealCellTypes==2),:));
PCAcoeffMean3 = [PCAcoeffMeanH93;PCAcoeffMeanHEK3];

figure
gscatter(PCAcoeff3(:,1),PCAcoeff3(:,2),RealCellTypes)
hold on
scatter(PCAcoeffMean3(:,1),PCAcoeffMean3(:,2),'ok','filled')
xlabel('PCA 1')
ylabel('PCA 2')

figure
bar(explained3(1:6))
%%
%Correlations between PCA and the sequencing depth per cell (Thinking that
%PCA 1 is highly correlated to the sequencing depth
for i=1:1:20
   PCA_To_Depth3(i) = corr(PCAcoeff3(:,i),RealCountPerCell','Type','Pearson');
%    PCA_To_Depth2(i) = corr(PCAcoeff(:,i),CountPerCell','Type','Spearman');
   PCA_To_DepthLookedAt3(i) = i;
   Explained_To_DepthLookedAt3(i) = explained3(i);
end
figure
scatter(PCA_To_DepthLookedAt3,PCA_To_Depth3,'ok','filled')
hold on
scatter(PCA_To_DepthLookedAt3,abs(PCAcoeffMean3(1,1:20)-PCAcoeffMean3(2,1:20)),'or','filled')
legend("r PCA to Depth","Mean Dist between cell types")
% print(gcf, 'WCG_H9_HEK_PCA_vs_Counts_InitialRandomizedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

figure
scatter(PCAcoeff3(:,1),RealCountPerCell,'ok','filled')
xlabel('PCA 1')
ylabel('Detected Marks')

%%
%Cluster on PCA 1,2,3 based on elbow plot and corr with sequencing depth
figure
D3 = pdist(PCAcoeff3(:,1:3));
Z3 = linkage(D3,'ward');
[h3,t3,perm3] = dendrogram(Z3,0);
ClustNumberDendo3temp = cluster(Z3,'MaxClust',2);

%Rearragne numbers to match cell types numbers
clear ClustNumberDendo3
ClustNumberDendo3(find(ClustNumberDendo3temp ==1),1) = 2;
ClustNumberDendo3(find(ClustNumberDendo3temp ==2),1) = 1;


%See if cell type and cluster asignment match
CorrectAsignmentnon = RealCellTypes == ClustNumberDendo3;
CorrectAsignmentnon2 = find(CorrectAsignmentnon==0);

%%
figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Ariel','DefaultAxesFontName','Ariel')
gscat1 = gscatter(PCAcoeff3(:,1),PCAcoeff3(:,2),RealCellTypes)
hold on
gscatter(PCAcoeff3(CorrectAsignmentnon2,1),PCAcoeff3(CorrectAsignmentnon2,2),CorrectAsignmentnon(CorrectAsignmentnon2),'k','o')
xlabel('PCA 1')
ylabel('PCA 2')
ylim([-0.15,0.32])
% yticks([0.1,0.2,0.3,0.4])
xlim([-0.22 0.22])
%xticks([0,10,20,30,510,520])
leg1 =legend(gscat1,{'H9','HEK293T'},'location','northeast')
% set(leg1,'Position',[0.35, 0.78, .15, .15])
box off
grid off
legend boxoff
set(gcf,'color','w')


% print(gcf, 'WCG_H9_HEK_InitialRandomizedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

%%
% figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Ariel','DefaultAxesFontName','Ariel')
% gscat1 = gscatter(PCAcoeff3(:,1),PCAcoeff3(:,3),RealCellTypes)
% hold on
% gscatter(PCAcoeff3(CorrectAsignmentnon2,1),PCAcoeff3(CorrectAsignmentnon2,3),CorrectAsignmentnon(CorrectAsignmentnon2),'k','o')
% xlabel('PCA 1')
% ylabel('PCA 3')
% ylim([-0.02,0.22])
% % yticks([0.1,0.2,0.3,0.4])
% xlim([-0.15 0.15])
% %xticks([0,10,20,30,510,520])
% leg1 =legend(gscat1,{'H9','HEK293T'},'location','north')
% % set(leg1,'Position',[0.35, 0.78, .15, .15])
% box off
% grid off
% legend boxoff
% set(gcf,'color','w')
% 
% 
% % print(gcf, 'GC_H9_HEK_InitialRandomizedCellType_5kbBins-2.png', '-dpng', '-r600', '-painters')
