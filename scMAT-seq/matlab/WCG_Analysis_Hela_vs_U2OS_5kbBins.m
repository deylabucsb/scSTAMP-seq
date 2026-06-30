clc
clear
minReadsNeeded = 10000; %may need to change this for running endogenous WCG files
PercentofCellsNeeded = 0.2;
%% Cell Type 2 (Hela)
%read in data
%a = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib3\Lib3_AllFaba\P39_5mC_Access_mRNA_Fresh_H9_Lib3_24_Cells_A_AllSeqRun_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
a = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

%aa = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib4\Lib4_AllFaba\P39_5mC_Access_mRNA_Fresh_H9_Lib4_AllSeqRun_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
aa = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

% MULTIPLEX DETERMINED CELL TYPES
%b = table2array(a(1:606218,2:97)); %remove Y chromosome
cellType_2_mp_indx = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/cellType_clustering_in_R/scMATseq_TGS_P4L1_multiplex_celltype_1_indxNumbs.csv');
cellType_2_mp_indx = table2array(cellType_2_mp_indx); %convert table to array
cellType_2_mp_indx = cellType_2_mp_indx + 1; %add 1 to each row of array
cellType_2_mp_indx = cellType_2_mp_indx'; %transpose array
b = table2array(a(1:606218,cellType_2_mp_indx)); %remove Y chromosome (Replace 2:97 with the cell numbers I want (+1) for cluster 1)

% SEURAT CLUSTERING DETERMINGED CELL TYPES
%bb = table2array(aa(1:606218,2:97));
cellType_2_sc_indx = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/cellType_clustering_in_R/scMATseq_TGS_P4L1_seuratCluster_celltype_1_indxNumbs.csv');
cellType_2_sc_indx = table2array(cellType_2_sc_indx); %convert table to array
cellType_2_sc_indx = cellType_2_sc_indx + 1; %add 1 to each row of array
cellType_2_sc_indx = cellType_2_sc_indx'; %transpose array
bb = table2array(aa(1:606218,cellType_2_sc_indx)); %(Replace 2:97 with the cell numbers I want (+1) for cluster 1)


RownamesTemp = table2array(a(1:606218,1));


% UNCOMMENT BELOW FOR MULTIPLEX CLUSTERING-----
%Celltypetemp = b; % Multiplex Assigned Cell Types

% UNCOMMENT BELOW FOR SEURAT CLUSTERING-----
Celltypetemp = bb; % Seurat Clustering Assigned Cell Types


%Celltypetemp = [b,bb];
TotperCellTypeHela = sum(Celltypetemp);
KeepCellsType = find(TotperCellTypeHela>minReadsNeeded);
CellTypetempHela = Celltypetemp(:,KeepCellsType);
BinaryCellTypetempHela = CellTypetempHela > 0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts (doesn't appear to make my data look different)
% BinaryCellTypetempHela = CellTypetemp2;

% KeepRowsHela = find(sum(BinaryCellTypetempHela,2) > ceil(PercentofCellsNeeded*length(KeepCellsType)));

% Use these variables to add spatial ratio information (chunk after figure 10)
KeepHelaCellsIndxNum_multiplex = (cellType_2_mp_indx(KeepCellsType) - 1)'; % subtract one to get original index number
KeepHelaCellsIndxNum_SeuratCluster = (cellType_2_sc_indx(KeepCellsType) -1)'; % subtract one to get original index number

clear a aa b bb Celltypetemp KeepCellsType

%% Cell Type 1 (U2OS)

%a = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\Plate4_FrozenHEK293T_mRNA_5mC_Accessibilty\AllCombined_5mC\P4_AllCondition_CellNumberCorrected_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
a = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

%aa = readtable('C:\Users\alexc\Box Sync\Research\scMAT-seq\Plate3_McviPI\P3_McvipI_1C-3_HEK293T_48Cells2-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);
aa = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_TotalDetection_5kb_Bins.txt','headerlines',1);

% MULTIPLEX DETERMINED CELL TYPES
%b = table2array(a(1:606218,2:97)); %(Replace 2:97 with the cell numbers I want (+1) for cluster 0)
cellType_1_mp_indx = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/cellType_clustering_in_R/scMATseq_TGS_P4L1_multiplex_celltype_0_indxNumbs.csv');
cellType_1_mp_indx = table2array(cellType_1_mp_indx); %convert table to array
cellType_1_mp_indx = cellType_1_mp_indx + 1; %add 1 to each row of array
cellType_1_mp_indx = cellType_1_mp_indx'; %transpose array
b = table2array(a(1:606218,cellType_1_mp_indx)); %remove Y chromosome (Replace 2:97 with the cell numbers I want (+1) for cluster 0)

% SEURAT CLUSTERING DETERMINGED CELL TYPES
%bb = table2array(aa(1:606218,2:97)); %(Replace 2:97 with the cell numbers I want (+1) for cluster 0)
cellType_1_sc_indx = readtable('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/cellType_clustering_in_R/scMATseq_TGS_P4L1_seuratCluster_celltype_0_indxNumbs.csv');
cellType_1_sc_indx = table2array(cellType_1_sc_indx); %convert table to array
cellType_1_sc_indx = cellType_1_sc_indx + 1; %add 1 to each row of array
cellType_1_sc_indx = cellType_1_sc_indx'; %transpose array
bb = table2array(aa(1:606218,cellType_1_sc_indx)); %(Replace 2:97 with the cell numbers I want (+1) for cluster 0)


RownamesTemp = table2array(a(1:606218,1));


% UNCOMMENT BELOW FOR MULTIPLEX CLUSTERING-----
%Celltypetemp = b; % Multiplex Assigned Cell Types

% UNCOMMENT BELOW FOR SEURAT CLUSTERING-----
Celltypetemp = bb; % Seurat Clustering Assigned Cell Types


%Celltypetemp = [b,bb];
TotperCellTypeU2OS = sum(Celltypetemp);
KeepCellsType = find(TotperCellTypeU2OS>minReadsNeeded);
CellTypetempU2OS = Celltypetemp(:,KeepCellsType);
BinaryCellTypetempU2OS = CellTypetempU2OS > 0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts (doesn't appear to make my data look different)
% BinaryCellTypetempU2OS = CellTypetemp2;

% KeepRowsU2OS = find(sum(BinaryCellTypetempU2OS,2) > ceil(PercentofCellsNeeded*length(KeepCellsType)));

% Use these variables to add spatial ratio information (chunk after figure 10)
KeepU2OSCellsIndxNum_multiplex = (cellType_1_mp_indx(KeepCellsType) - 1)'; % subtract one to get original index number
KeepU2OSCellsIndxNum_SeuratCluster = (cellType_1_sc_indx(KeepCellsType) -1)'; % subtract one to get original index number

clear a aa b bb Celltypetemp KeepCellsType CellTypetemp2

%%
%%Lets remove bins that are clearly way to high (i.e. blacklist)
combinedtemp = [CellTypetempHela,CellTypetempU2OS];
KeepNonBlacklist = find(sum(combinedtemp,2)<size(combinedtemp,2)*3);
combinedtempNonBL = combinedtemp(KeepNonBlacklist,:);
RownamesTemp2 = RownamesTemp(KeepNonBlacklist,:);
CountPerCelltemp2 = sum(combinedtempNonBL,1);

combinedtemp2 = combinedtempNonBL >0; %Makes it binary, shown to be useful to help remove sequencing depth artifacts
CountPerCell = sum(combinedtemp);

%%
%Using the Cell type information we know
BinaryCellTypetempHela = combinedtemp2(:,1:size(CellTypetempHela,2));
BinaryCellTypetempU2OS = combinedtemp2(:,1+size(CellTypetempHela,2):end);

pcrtilekeep = 98;
KeepRowsHela = find(sum(BinaryCellTypetempHela,2) > prctile(sum(BinaryCellTypetempHela,2),pcrtilekeep));
KeepRowsU2OS = find(sum(BinaryCellTypetempU2OS,2) > prctile(sum(BinaryCellTypetempU2OS,2),pcrtilekeep));
KeepRows = unique([KeepRowsHela;KeepRowsU2OS]);

AllBinaryTemp = [BinaryCellTypetempHela,BinaryCellTypetempU2OS];

AllBinary = double(AllBinaryTemp(KeepRows,:));
AllRownames = RownamesTemp2(KeepRows,:);


% BinaryHela = BinaryCellTypetempHela(KeepRows,:);
% HelaRownames = RownamesTemp(KeepRows,:);
% 
% BinaryU2OS = BinaryCellTypetempU2OS(KeepRows,:);
% U2OSRownames = RownamesTemp(KeepRows,:);

%% Figure 1
%Which Bins are used
BinarySTD = std(combinedtempNonBL,0,2);
BinaryMean = mean(combinedtempNonBL,2);
BinaryCV = BinarySTD./BinaryMean;

sz1=20;
sz2=10;
figure
scatter(BinaryMean,log10(BinaryCV),sz1,'sk','filled')
hold on
scatter(BinaryMean(KeepRowsHela),log10(BinaryCV(KeepRowsHela))+1,sz2,'or')
hold on
scatter(BinaryMean(KeepRowsU2OS),log10(BinaryCV(KeepRowsU2OS))+2,sz2,'ob')
xlabel("Mean Bin Val")
ylabel("Log10 CV")

%% Figures 2-3
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

%% Figures 4-6
[PCAcoeff,score,latent,tsquared,explained,mu] = pca(AllBinary);
clear CellTypes
CellTypes(1:size(BinaryCellTypetempHela,2),1) = 1; % 1 is Hela
CellTypes(size(BinaryCellTypetempHela,2)+1:size(AllBinary,2)) = 2; % 2 is U2OS

% BuffersUsedtemp(1:size(BinaryCellTypetempHela,2),1) = "FSB-GCB.200nM";
% BuffersUsed = [BuffersUsedtemp;BufferU2OS];

PCAcoeffMeanHela = mean(PCAcoeff(1:size(CellTypetempHela,2),:));
PCAcoeffMeanU2OS = mean(PCAcoeff(1+size(CellTypetempHela,2):end,:));
PCAcoeffMean = [PCAcoeffMeanHela;PCAcoeffMeanU2OS];

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

% figure
% gscatter(PCAcoeff(:,1),PCAcoeff(:,2),BuffersUsed)
% xlabel('PCA 1')
% ylabel('PCA 2')

figure
bar(explained(1:6))

%% Figures 7-8
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
% print(gcf, 'WCG_Hela_U2OS_PCA_vs_Counts_CurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')


figure
scatter(PCAcoeff(:,2),CountPerCell','ok','filled')
% hold on
% plot([-0.2,0.2],[abc(2)+abc(1)*-0.2,abc(2)+abc(1)*0.2])
xlabel('PCA 2')
ylabel('Detected Marks')

%% Figure 9
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

%% Figure 10 THIS THE ONE
%%% DNA Endogenous 5mCpG PCA1vPCA2 Plot
figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial')
gscat1 = gscatter(PCAcoeff(:,1),PCAcoeff(:,2),CellTypes,[165/255,25/255,165/255;1,120/255,0/255])
xlabel('PCA 1')
ylabel('PCA 2')
title("Endogeneous CpG Methylation")
leg1 =legend(gscat1,{'Hela','U2OS'})
box off
grid off
legend boxoff
set(gcf,'color','w')

%%% Correct Assignment Overlay Plot
% figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial')
% gscatter(PCAcoeff(CorrectAsignment2,1),PCAcoeff(CorrectAsignment2,2),CorrectAsignment(CorrectAsignment2),'k','o')
% xlabel('PCA 1')
% ylabel('PCA 2')
% title("Endogeneous CpG Methylation")
% leg1 =legend(gscat1,{'Hela','U2OS'})
% box off
% grid off
% legend boxoff
% set(gcf,'color','w')

%%% Both Plots Overlayed
% figure('DefaultAxesFontSize',16,'DefaultTextFontSize',16,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial')
% gscat1 = gscatter(PCAcoeff(:,1),PCAcoeff(:,2),CellTypes,[165/255,25/255,165/255;1,120/255,0/255])
% hold on
% gscatter(PCAcoeff(CorrectAsignment2,1),PCAcoeff(CorrectAsignment2,2),CorrectAsignment(CorrectAsignment2),'k','o')
% xlabel('PCA 1')
% ylabel('PCA 2')
% title("Endogeneous CpG Methylation")
% % ylim([0,0.225])
% % yticks([0.1,0.2,0.3,0.4])
% % xlim([-0.2 0.25])
% % xticks([0,10,20,30,510,520])
% leg1 =legend(gscat1,{'Hela','U2OS'})
% % % set(leg1,'Position',[0.75, 0.20, .15, .15])
% box off
% grid off
% legend boxoff
% set(gcf,'color','w')


% print(gcf, 'WCG_Hela_U2OS_CurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')
% saveas(gcf,'WCG_Hela_U2OS_CurratedCellType_5kbBins.pdf')

%% To add spatial ratio information, extract cell index information below into csv

% cd("/Users/piscopio/Downloads");
% mkdir(strcat("WCG_Analysis_Hela_vs_U2OS_5kbBins_", date));
% cd(strcat("WCG_Analysis_Hela_vs_U2OS_5kbBins_", date));
% 
% writematrix(CellTypes, 'CellTypes_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'); % add this arrary to 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'
% 
% % writematrix(KeepHelaCellsIndxNum_multiplex, 'IndxCellNums_Hela_multiplex_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv');         % if multiplex_cluster used:    add this arrary to 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'
% writematrix(KeepHelaCellsIndxNum_SeuratCluster, 'IndxCellNums_Hela_SeuratCluster_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'); % if seurat_cluster used:       add this arrary to 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'
% % writematrix(KeepU2OSCellsIndxNum_multiplex, 'IndxCellNums_U2OS_multiplex_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv');         % if multiplex_cluster used:    add this arrary to 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'
% writematrix(KeepU2OSCellsIndxNum_SeuratCluster, 'IndxCellNums_U2OS_SeuratCluster_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'); % if seurat_cluster used:       add this arrary to 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv'
% 
% writematrix(PCAcoeff, 'PCAcoeff_(WCG_Analysis_Hela_vs_U2OS_5kbBins).csv');
% 
% cd("/Users/piscopio/Desktop");

%%

%%---%%
%Without using cell type information
KeepRowsNonCellType = find(sum(combinedtemp2,2) > prctile(sum(combinedtemp2,2),pcrtilekeep));

AllBinary2 = double(combinedtemp2(KeepRowsNonCellType,:));
AllRownames2 = RownamesTemp2(KeepRowsNonCellType,:);


%% Figures 11-12
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

%% Figures 13-14
[PCAcoeff2,score2,latent2,tsquared2,explained2,mu2] = pca(AllBinary2);
CellTypes(1:size(BinaryCellTypetempHela,2),1) = 1; % 1 is Hela
CellTypes(size(BinaryCellTypetempHela,2)+1:size(AllBinary2,2)) = 2; % 2 is U2OS

PCAcoeffMeanHela2 = mean(PCAcoeff2(1:size(CellTypetempHela,2),:));
PCAcoeffMeanU2OS2 = mean(PCAcoeff2(1+size(CellTypetempHela,2):end,:));
PCAcoeffMean2 = [PCAcoeffMeanHela2;PCAcoeffMeanU2OS2];

figure
gscatter(PCAcoeff2(:,1),PCAcoeff2(:,2),CellTypes)
hold on
scatter(PCAcoeffMean2(:,1),PCAcoeffMean2(:,2),'ok','filled')
xlabel('PCA 1')
ylabel('PCA 2')

figure
bar(explained2(1:6))
%% Figures 15-16
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
% print(gcf, 'WCG_Hela_U2OS_PCA_vs_Counts_nonCurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

figure
scatter(PCAcoeff2(:,1),CountPerCell,'ok','filled')
xlabel('PCA 1')
ylabel('Detected Marks')

%% Figure 17
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

%% Figure 18
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
leg1 =legend(gscat1,{'Hela','U2OS'},'location','north')
set(leg1,'Position',[0.46, 0.78, .15, .15])
box off
grid off
legend boxoff
set(gcf,'color','w')


% print(gcf, 'WCG_Hela_U2OS_NonCurratedCellType_5kbBins.png', '-dpng', '-r600', '-painters')


%%

%%---%%
%Miss assigning cell type information then trying it
RealCellTypes = CellTypes(rndPermCells);
RealCountPerCell = CountPerCell(rndPermCells);

BinaryCellTypetempHela2 = combinedtemp2(:,rndPermCells(1:size(CellTypetempHela,2)));
BinaryCellTypetempU2OS2 = combinedtemp2(:,rndPermCells(1+size(CellTypetempHela,2):end));

pcrtilekeep = 98;
KeepRowsHela2 = find(sum(BinaryCellTypetempHela2,2) > prctile(sum(BinaryCellTypetempHela2,2),pcrtilekeep));
KeepRowsU2OS2 = find(sum(BinaryCellTypetempU2OS2,2) > prctile(sum(BinaryCellTypetempU2OS2,2),pcrtilekeep));
KeepRows2 = unique([KeepRowsHela2;KeepRowsU2OS2]);

AllBinaryTemp2 = [BinaryCellTypetempHela2,BinaryCellTypetempU2OS2];

AllBinary3 = double(AllBinaryTemp2(KeepRows2,:));
AllRownames3 = RownamesTemp2(KeepRows2,:);


%% Figure 19-20
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

%% Figures 21-22
[PCAcoeff3,score3,latent3,tsquared3,explained3,mu3] = pca(AllBinary3);


PCAcoeffMeanHela3 = mean(PCAcoeff3(find(RealCellTypes==1),:));
PCAcoeffMeanU2OS3 = mean(PCAcoeff3(find(RealCellTypes==2),:));
PCAcoeffMean3 = [PCAcoeffMeanHela3;PCAcoeffMeanU2OS3];

figure
gscatter(PCAcoeff3(:,1),PCAcoeff3(:,2),RealCellTypes)
hold on
scatter(PCAcoeffMean3(:,1),PCAcoeffMean3(:,2),'ok','filled')
xlabel('PCA 1')
ylabel('PCA 2')

figure
bar(explained3(1:6))
%% Figures 23-24
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
% print(gcf, 'WCG_Hela_U2OS_PCA_vs_Counts_InitialRandomizedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

figure
scatter(PCAcoeff3(:,1),RealCountPerCell,'ok','filled')
xlabel('PCA 1')
ylabel('Detected Marks')

%% Figure 25
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

%% Figure 26
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
leg1 =legend(gscat1,{'Hela','U2OS'},'location','northeast')
% set(leg1,'Position',[0.35, 0.78, .15, .15])
box off
grid off
legend boxoff
set(gcf,'color','w')


% print(gcf, 'WCG_Hela_U2OS_InitialRandomizedCellType_5kbBins.png', '-dpng', '-r600', '-painters')

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
% leg1 =legend(gscat1,{'Hela','U2OS'},'location','north')
% % set(leg1,'Position',[0.35, 0.78, .15, .15])
% box off
% grid off
% legend boxoff
% set(gcf,'color','w')
% 
% 
% % print(gcf, 'GC_Hela_U2OS_InitialRandomizedCellType_5kbBins-2.png', '-dpng', '-r600', '-painters')
