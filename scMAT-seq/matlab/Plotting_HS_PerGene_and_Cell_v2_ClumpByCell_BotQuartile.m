%Plotting UCSC TSS to TSE
clear
clc

%load('wgEncodeAwgDnaseDukeH9esUniPk_Reformated_BottomQuartile.mat');

% REFERENCE FILES PATH
addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/ReferenceData/ReferenceMatrixFiles');

% HS DATA PATH
addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/HypersensitivityMatrixFiles/20230130_CellClusters_(10kThreshold)');
%addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61/Genomic/HypersensitivityData/HypersensitivityMatrixFiles/20230130_CellClusters_(10kThreshold)');

%addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Exp02_(Run56_Run61)/Genomic/ReferenceData/MatrixFilesSaved/')
%addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Exp02_(Run56_Run61)/Genomic/HypersensitivityData/MatrixFilesGenerated/20230130_CellClusters_(10kThreshold)');

%addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/Hela/HSMatrixFiles');
%addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/ReferenceData/ReferenceMatrixFiles')
%addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/MAT_FABA_ACCESS/HSMatrixFiles/20230509_10kThreshold/Cluster0');
%addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/MAT_FABA_ACCESS/HSMatrixFiles/20230509_10kThreshold/Cluster1');

load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_BottomQuartile.mat');


Sex = 23; %23 if female, 24 if male

% MovingSize=10;
% WindowSize=100;
% 
% MinimumNeeded_GC =10000;

%% 
% Accessibiilty (GC)

%Background Model based on H9 Open DNA (proteased) and Full GC
%load('Bulk_GC_2020_1_21_H9_Hypersensitivity_Bottom_Quartile_Per_Cell_25BP_Bin_plusMinus2000.mat')


%Background Model based on H9 Open DNA (proteased) and Full GC
%load('Bulk_GC_2022-07-19_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat');

%Background Model based on Hela Open DNA (proteased) and Full GC
%load('Bulk_GC_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_PlusMinus3000.mat')
load('Bulk_GC_2022-07-19_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat');


BulkBackGroundPerGene = CellGeneDistStruct;

%Load Genome WCG
%load('../WCG_Genome/Genome_WCG_2020_4_24_H9_Hypersensitivity_Bottom_Quartile_Per_Cell_25BP_Bin_plusMinus2000.mat')

%Load Genome WCG
%load('Genome_WCG_2022-07-19_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat');
%load('Genome_WCG_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_PlusMinus3000.mat');
load('Genome_WCG_2022-07-19_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat');


GenomeWCGBackGroundPerGene = CellGeneDistStruct;


%% ACCESSIBILITY GC old

%---%
%Looking at Plate 39 Lib 3
%load('P39_Lib3_GC_2020_1_21_H9_Hypersenstivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat')

%load('P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat');

%load('R_Reformatted_Cluster_0_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat');
load('R_Reformatted_Cluster_1_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat');

%% ACCESSIBILITY GC

        %%%%% Cluster 0 %%%%%

%load('Run71_Cluster_0_P3L1_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L2_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L3_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L4_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');

        %%%%% Cluster 1 %%%%%

%load('Run71_Cluster_1_P3L1_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L2_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L3_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L4_ACCESS_GC_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');

%%

%CellGeneDistStruct_P39_GC_Lib3 = CellGeneDistStruct;
%CellGeneDistStruct_P4L1_GC = CellGeneDistStruct;
CellGeneDistStruct_GC = CellGeneDistStruct;

fn = fieldnames(TSSstruct);
%TSSVecByGene_P39_GC_Lib3 = struct([]);
%TSSVecByGene_P4L1_GC = struct([]);
TSSVecByGene_GC = struct([]);

% for k=1:Sex
%   for i=1:1:length(WorkingCells)
%        j=WorkingCells(i);
%        
%        
%        if i==1
%         TSSVecByGene_P39_GC_Lib3(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
%      
%        else
%        TSSVecByGene_P39_GC_Lib3(1).(fn{k}) = TSSVecByGene_P39_GC_Lib3(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
% %        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
%        end
% %     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
%   end 
%   TSSVecByChr_P39_GC_Lib3(k,:)=sum(TSSVecByGene_P39_GC_Lib3(1).(fn{k}));
% end

for i=1:1:length(WorkingCells)
       j=WorkingCells(i);
       for k=1:Sex
        %temp(k,:)= sum(CellGeneDistStruct_P39_GC_Lib3(j).(fn{k}),1);
        %temp(k,:)= sum(CellGeneDistStruct_P4L1_GC(j).(fn{k}),1);
        temp(k,:)= sum(CellGeneDistStruct_GC(j).(fn{k}),1);
       end
     %TSSVecByCell_P39_GC_Lib3(i,:)=sum(temp);
     %TSSVecByCell_P4L1_GC(i,:)=sum(temp);
     TSSVecByCell_GC(i,:)=sum(temp);
     %CountsPerCell_P39_GC_Lib3(i,:)=CountsPerCell(j);
     %CountsPerCell_P4L1_GC(i,:)=CountsPerCell(j);
     CountsPerCell_GC(i,:)=CountsPerCell(j);
end

%%
% 
% %---%
% % addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\scMAT-seq_Figures\H9_HypersenstivitySites\HEK293T\P4')
% %Looking at Plate 39 Lib 4
% load('P39_Lib4_GC_2020_1_21_H9_Hypersenstivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat')
% 
% CellGeneDistStruct_P39_GC_Lib4 = CellGeneDistStruct;
% 
% fn = fieldnames(TSSstruct);
% TSSVecByGene_P39_GC_Lib4 = struct([]);
% 
% % for k=1:Sex
% %   for i=1:1:length(WorkingCells)
% %        j=WorkingCells(i);
% %        if i==1
% %         TSSVecByGene_P39_GC_Lib4(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
% %      
% %        else
% %        TSSVecByGene_P39_GC_Lib4(1).(fn{k}) = TSSVecByGene_P39_GC_Lib4(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
% % %        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
% %        end
% % %     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
% %   end 
% %   TSSVecByChr_P39_GC_Lib4(k,:)=sum(TSSVecByGene_P39_GC_Lib4(1).(fn{k}));
% % end
% 
% for i=1:1:length(WorkingCells)
%        j=WorkingCells(i);
%        for k=1:Sex
%         temp(k,:)= sum(CellGeneDistStruct_P39_GC_Lib4(j).(fn{k}),1);
%        end
%      TSSVecByCell_P39_GC_Lib4(i,:)=sum(temp);
%      CountsPerCell_P39_GC_Lib4(i,:)=CountsPerCell(j);
% end

%%
clear TSSVecByChr_NormalizedByBulkBackground
BinSize2 = 25; %25 means do not adjust bin size to make it larger bins (previous scripts used a bin size of 25)
EdgesOfBins2=-Distance:BinSize2:Distance;
[RightBin NewBinIndex] = discretize(BulkBackGroundPerGene(2).(fn{1}),EdgesOfBins2);
RightBin_Useful = RightBin(1:length(RightBin)-1);
BulkBackGroundPerGene_NewBins = struct([]);
% TSSVecByGene_P39_GC_Total = struct([]);
% TSSVecByGene_P39_GC_Total_NewBins = struct([]);
% TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground = struct([]);
% TSSVecByGene_oneVec_P39_GC_Total_NormalizedByBulkBackground = [];
TSSGeneNameVect = [];

GenomeWCGBackGroundPerGene_NewBins = struct([]);

for k=1:Sex
%     TSSVecByGene_P39_GC_Total(1).(fn{k}) = TSSVecByGene_P39_GC_Lib3.(fn{k})+TSSVecByGene_P39_GC_Lib4.(fn{k});
        
    for i = 1:1:length(NewBinIndex)-1
        tempIndex = find(RightBin_Useful==i);
        BulkBackGroundPerGene_NewBins(1).(fn{k})(:,i) = sum(BulkBackGroundPerGene(1).(fn{k})(:,tempIndex),2);
%         TSSVecByGene_P39_GC_Total_NewBins(1).(fn{k})(:,i) = sum(TSSVecByGene_P39_GC_Total(1).(fn{k})(:,tempIndex),2);
        GenomeWCGBackGroundPerGene_NewBins(1).(fn{k})(:,i) = sum(GenomeWCGBackGroundPerGene(1).(fn{k})(:,tempIndex),2);

    end
    BulkBackGroundPerGene_NewBins(2).(fn{k}) = NewBinIndex;
%     TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground(1).(fn{k}) = TSSVecByGene_P39_GC_Total_NewBins.(fn{k}) ./ (BulkBackGroundPerGene_NewBins(1).(fn{k}) + TSSVecByGene_P39_GC_Total_NewBins.(fn{k})); %Normalizes by GC detected in bulk and found in single cells
%     TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground(1).(fn{k}) = (TSSVecByGene_P39_GC_Total_NewBins.(fn{k})) ./ (BulkBackGroundPerGene_NewBins(1).(fn{k})+1); %added 1 to each bin but could add a factor that relates to the number of GC (incase it is found in sc but not in bulk (which happens))
%     TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground(1).(fn{k}) = (TSSVecByGene_P39_GC_Total_NewBins.(fn{k}));
            GenomeWCGBackGroundPerGene_NewBins(2).(fn{k}) = NewBinIndex;

    
%     TSSVecByChr_NormalizedByBulkBackground(k,:) = nansum(TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground(1).(fn{k}));
%     TSSVecByGene_oneVec_P39_GC_Total_NormalizedByBulkBackground = [TSSVecByGene_oneVec_P39_GC_Total_NormalizedByBulkBackground;TSSVecByGene_P39_GC_Total_NormalizedByBulkBackground(1).(fn{k})];
    TSSGeneNameVect = [TSSGeneNameVect;TSSstruct(6).(fn{k})];
end

%%
bin_mdpt2=(EdgesOfBins2(2:end)+EdgesOfBins2(1:(end-1)))/2;
% TSSVecTotal_NormalizedByBulkBackground = sum(TSSVecByChr_NormalizedByBulkBackground);
% figure
% plot(bin_mdpt2,TSSVecTotal_NormalizedByBulkBackground)

%%

% AllEnhancerBins_Norm = nansum(TSSVecByGene_oneVec_P39_GC_Total_NormalizedByBulkBackground)/length(TSSGeneNameVect);
% MovingAvg_GC_Norm = movmean(AllEnhancerBins_Norm,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],'Endpoints','discard');


%%
% figure
% plot(bin_mdpt2MovingAvg,MovingAvg_GC_Norm,'k')
% % ylim([0.02 0.08])
% xlim([-2600 2600])
% xticks([-2000  0  2000])
% % xlabs = xticklabels; xlabs(2) = {'0'}; xticklabels(xlabs);
% ylabel("Normalized GC Detection")
% title("Active H9 Enhancers")
% % legend(["Top";"Middle";"Bottom";"Undetected"])
% box off
% set(gcf,'color','w')
% print(gcf, 'H9_Normalized_ToBulk_Accessibility_GC_Active_Enhancers.png', '-dpng', '-r600', '-painters')


%%
%Look at bulk fully open chromatin detection
BulkBackGround_TSSVecByGene_oneVec = [];
for k=1:Sex
    temp=sum(BulkBackGroundPerGene_NewBins(1).(fn{k}));
    BulkBackGround_TSSVecByGene_oneVec = [BulkBackGround_TSSVecByGene_oneVec;temp];
end
BulkBackGround = sum(BulkBackGround_TSSVecByGene_oneVec);

%WCG Genome background
GenomeWCGBackGroundPerGene_oneVec = [];
for k=1:Sex
    temp=sum(GenomeWCGBackGroundPerGene_NewBins(1).(fn{k}));
    GenomeWCGBackGroundPerGene_oneVec = [GenomeWCGBackGroundPerGene_oneVec;temp];
end
WCGGenomeBackGround = sum(GenomeWCGBackGroundPerGene_oneVec);
%%
%By Cell Moving Average
%GC_Norm_ByCell = [TSSVecByCell_P39_GC_Lib3;TSSVecByCell_P39_GC_Lib4]./BulkBackGround.*median(BulkBackGround);
%GC_Norm_ByCell = (TSSVecByCell_P4L1_GC)./BulkBackGround.*median(BulkBackGround);
GC_Norm_ByCell = (TSSVecByCell_GC)./BulkBackGround.*median(BulkBackGround);

%CountsPerCell_GC = [CountsPerCell_P39_GC_Lib3;CountsPerCell_P39_GC_Lib4];
%CountsPerCell_GC = CountsPerCell_P4L1_GC;
CountsPerCell_GC = CountsPerCell_GC;

%% ENDOGENOUS WCG 

%Now Look at WCG Sites

%Looking at Plate 39 Lib 3
%load('P39_Lib3_WCG_2020_1_21_H9_Hypersenstivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat')

%load('P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat')

%load('R_Reformatted_Cluster_0_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat');
load('R_Reformatted_Cluster_1_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat');

        %%%%% Cluster 0 %%%%%

%load('Run71_Cluster_0_P3L1_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L2_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L3_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_0_P3L4_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');

        %%%%% Cluster 1 %%%%%

%load('Run71_Cluster_1_P3L1_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L2_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L3_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');
%load('Run71_Cluster_1_P3L4_ENDOG_WCG_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_PlusMinus3000.mat');

%CellGeneDistStruct_P39_WCG_Lib3 = CellGeneDistStruct;
%CellGeneDistStruct_P4L1_WCG = CellGeneDistStruct;
CellGeneDistStruct_WCG = CellGeneDistStruct;

fn = fieldnames(TSSstruct);
%TSSVecByGene_P39_WCG_Lib3 = struct([]);
%TSSVecByGene_P4L1_WCG = struct([]);
TSSVecByGene_WCG = struct([]);

% for k=1:Sex
%   for i=1:1:length(WorkingCells)
%        j=WorkingCells(i);
%        if i==1
%         TSSVecByGene_P39_WCG_Lib3(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
%      
%        else
%        TSSVecByGene_P39_WCG_Lib3(1).(fn{k}) = TSSVecByGene_P39_WCG_Lib3(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
% %        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
%        end
% %     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
%   end 
%   TSSVecByChr_P39_WCG_Lib3(k,:)=sum(TSSVecByGene_P39_WCG_Lib3(1).(fn{k}));
% end

for i=1:1:length(WorkingCells)
    j=WorkingCells(i);
    for k=1:Sex
        %temp(k,:)= sum(CellGeneDistStruct_P39_WCG_Lib3(j).(fn{k}),1);
        %temp(k,:)= sum(CellGeneDistStruct_P4L1_WCG(j).(fn{k}),1);
        temp(k,:)= sum(CellGeneDistStruct_WCG(j).(fn{k}),1);
        %CellGeneDistStruct_P4L1_WCG = CellGeneDistStruct;
    end
    %TSSVecByCell_P39_WCG_Lib3(i,:)=sum(temp);
    %TSSVecByCell_P4L1_WCG(i,:)=sum(temp);
    TSSVecByCell_WCG(i,:)=sum(temp);
    %CountsPerCell_P39_WCG_Lib3(i,:)=CountsPerCell(j);
    %CountsPerCell_P4L1_WCG(i,:)=CountsPerCell(j);
    CountsPerCell_WCG(i,:)=CountsPerCell(j);
end

%%
% 
% %---%
% %Looking at Plate 39 Lib 4
% load('P39_Lib4_WCG_2020_1_21_H9_Hypersenstivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat')
% 
% CellGeneDistStruct_P39_WCG_Lib4 = CellGeneDistStruct;
% 
% fn = fieldnames(TSSstruct);
% TSSVecByGene_P39_WCG_Lib4 = struct([]);
% % 
% % for k=1:Sex
% %   for i=1:1:length(WorkingCells)
% %        j=WorkingCells(i);
% %        if i==1
% %         TSSVecByGene_P39_WCG_Lib4(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
% %      
% %        else
% %        TSSVecByGene_P39_WCG_Lib4(1).(fn{k}) = TSSVecByGene_P39_WCG_Lib4(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
% % %        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
% %        end
% % %     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
% %   end 
% %   TSSVecByChr_P39_WCG_Lib4(k,:)=sum(TSSVecByGene_P39_WCG_Lib4(1).(fn{k}));
% % end
% 
% for i=1:1:length(WorkingCells)
%        j=WorkingCells(i);
%        for k=1:Sex
%         temp(k,:)= sum(CellGeneDistStruct_P39_WCG_Lib4(j).(fn{k}),1);
%        end
%      TSSVecByCell_P39_WCG_Lib4(i,:)=sum(temp);
%      CountsPerCell_P39_WCG_Lib4(i,:)=CountsPerCell(j);
% end

%%
%Combine all cells for each gene
% TSSVecByGene_P39_WCG_Total = struct([]);
% TSSVecByGene_P39_WCG_Total_NewBins = struct([]);
% TSSVecByGene_P39_WCG = struct([]);
% TSSVecByGene_oneVec_P39_WCG_Total = [];
% 
% for k=1:Sex
%     TSSVecByGene_P39_WCG_Total(1).(fn{k}) = TSSVecByGene_P39_WCG_Lib3.(fn{k})+TSSVecByGene_P39_WCG_Lib4.(fn{k});
%         
%     for i = 1:1:length(NewBinIndex)-1
%         tempIndex = find(RightBin_Useful==i);
%         TSSVecByGene_P39_WCG_Total_NewBins(1).(fn{k})(:,i) = sum(TSSVecByGene_P39_WCG_Total(1).(fn{k})(:,tempIndex),2);
%     end
%     TSSVecByGene_P39_WCG(1).(fn{k}) = TSSVecByGene_P39_WCG_Total_NewBins.(fn{k}) ; %No Normalization for CG
%     TSSVecByChr_SumCells(k,:) = nansum(TSSVecByGene_P39_WCG(1).(fn{k}));
%     
%     TSSVecByGene_oneVec_P39_WCG_Total = [TSSVecByGene_oneVec_P39_WCG_Total;TSSVecByGene_P39_WCG(1).(fn{k})];
% end

%%
% TSSVecTotal_WCG = sum(TSSVecByChr_SumCells);
% figure
% plot(bin_mdpt2,TSSVecTotal_WCG)
% %%
% 
% 
% AllEnhancerBins_WCG = nansum(TSSVecByGene_oneVec_P39_WCG_Total);
% MovingAvg_WCG = movmean(AllEnhancerBins_WCG,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],'Endpoints','discard');

%%
%By Cell Moving Average
BinsToAverageBeforeandAfter = 1;

bin_mdpt2MovingAvg = movmean(bin_mdpt2,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],'Endpoints','discard');

BulkBackGround_MovingAvg = movmean(BulkBackGround,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],'Endpoints','discard');

MovingAvg_GC_Norm_ByCell = movmean(GC_Norm_ByCell,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],2,'Endpoints','discard');
MovingAvg_GC_Norm = sum(MovingAvg_GC_Norm_ByCell);

%WCG_ByCell = [TSSVecByCell_P39_WCG_Lib3;TSSVecByCell_P39_WCG_Lib4]./WCGGenomeBackGround.*median(WCGGenomeBackGround);
%WCG_ByCell = (TSSVecByCell_P4L1_WCG)./WCGGenomeBackGround.*median(WCGGenomeBackGround);
WCG_ByCell = (TSSVecByCell_WCG)./WCGGenomeBackGround.*median(WCGGenomeBackGround);

%CountsPerCell_WCG = [CountsPerCell_P39_WCG_Lib3;CountsPerCell_P39_WCG_Lib4];
%CountsPerCell_WCG = CountsPerCell_P4L1_WCG;
CountsPerCell_WCG = CountsPerCell_WCG;

MovingAvg_WCG_ByCell = movmean(WCG_ByCell,[BinsToAverageBeforeandAfter BinsToAverageBeforeandAfter],2,'Endpoints','discard');
MovingAvg_WCG = sum(MovingAvg_WCG_ByCell);
%%
% figure
% plot(bin_mdpt2MovingAvg,MovingAvg_WCG,'k')
% % ylim([0.02 0.08])
% xlim([-2600 2600])
% xticks([-2000  0  2000])
% xlabs = xticklabels; xlabs(2) = {'0'}; xticklabels(xlabs);
% ylabel("WCG Detection")
% title("Active H9 Enhancers")
% % legend(["Top";"Middle";"Bottom";"Undetected"],'Location','SouthEast')
% box off
% set(gcf,'color','w')
% % print(gcf, 'H9_Endogenous_WCG_Active_Enhancers.png', '-dpng', '-r600', '-painters')
%%
figure
plot(bin_mdpt2MovingAvg,BulkBackGround_MovingAvg,'k')
xlim([-2600 2600])
xticks([-2000  0  2000])
xlabs = xticklabels; xlabs(2) = {'Active H9 Enhancers'}; xticklabels(xlabs);
ylabel("GC Detection")
title("Bulk Fully Open Chromatin");
%legend(["Bulk Fully Open Chromatin"])
box off
set(gcf,'color','w')
% print(gcf, 'H9_Bulk_Fully_OpenChromatin_GC_Detection_atActiveEnhancers.png', '-dpng', '-r600', '-painters')

%%
%Plot Access and Endogenous Together
LNWidth = 1;
figure
plot(bin_mdpt2MovingAvg,MovingAvg_GC_Norm,'-k','LineWidth',LNWidth)
ylabel("GC Detection")
hold on
plot(bin_mdpt2MovingAvg,MovingAvg_WCG,'-b','LineWidth',LNWidth)
xlim([-1000 1000])
xticks([-500  0  500])
xlabs = xticklabels; xlabs(2) = {'0'}; xticklabels(xlabs);
ylabel("Detection")
%title("HEK293T Hypersenstivity Sites")
title("HeLa Hypersenstivity Sites")
legend(["Accessibility";"Endogenous"],'location','southoutside','orientation','horizontal')
box off
set(gcf,'color','w')


LNWidth = 1;
figure
plot(bin_mdpt2MovingAvg,MovingAvg_GC_Norm,'-k','LineWidth',LNWidth)
ylabel("GC Detection")
hold on
yyaxis right
plot(bin_mdpt2MovingAvg,MovingAvg_WCG,'-b','LineWidth',LNWidth)
xlim([-1000 1000])
xticks([-500  0  500])
xlabs = xticklabels; xlabs(2) = {'0'}; xticklabels(xlabs);
ylabel("WCG Detection")
set(gca,'ycolor','b')
%title("HEK293T Hypersenstivity Sites")
title("HeLa Hypersenstivity Sites")
legend(["Accessibility";"Endogenous"],'location','southoutside','orientation','horizontal')
box off
set(gcf,'color','w')
% print(gcf, 'H9_Normalized_ToBulk_Access_GC_and_Endogenous_WCG_Active_Enhancers_1k_150bp_MoveAvg.png', '-dpng', '-r600', '-painters')

%%
%addpath('D:\Box Sync\MATLAB')
addpath('/Users/piscopio/Library/CloudStorage/Box-Box/MATLAB/functions/')

LNWidth = 1;

%Doing a RPM Type thing
Mean_MovingAvg_WCG_ByCell=mean(MovingAvg_WCG_ByCell./CountsPerCell_WCG.*10^6./size(TSS,1));
std_MovingAvg_WCG_ByCell=std(MovingAvg_WCG_ByCell./CountsPerCell_WCG.*10^6./size(TSS,1));

Mean_MovingAvg_GC_Norm_ByCell=mean(MovingAvg_GC_Norm_ByCell./CountsPerCell_GC.*10^6./size(TSS,1));
std_MovingAvg_GC_Norm_ByCell=std(MovingAvg_GC_Norm_ByCell./CountsPerCell_GC.*10^6./size(TSS,1));
shadeval = 0; %up to 256 (RGB scale)
transparentval = 0.25; %from 0 to 1
ymaxval=0.025; %0.034


%%

figure('DefaultAxesFontSize',18,'DefaultTextFontSize',18,'DefaultTextFontName','Arial','DefaultAxesFontName','Arial')
shadedplot(bin_mdpt2MovingAvg,Mean_MovingAvg_GC_Norm_ByCell-std_MovingAvg_GC_Norm_ByCell,Mean_MovingAvg_GC_Norm_ByCell+std_MovingAvg_GC_Norm_ByCell,[shadeval/256, shadeval/256, 1],'none');
alpha(transparentval)
hold on
p1=plot(bin_mdpt2MovingAvg,Mean_MovingAvg_GC_Norm_ByCell,'-b','LineWidth',LNWidth);
hold on
shadedplot(bin_mdpt2MovingAvg,Mean_MovingAvg_WCG_ByCell-std_MovingAvg_WCG_ByCell,Mean_MovingAvg_WCG_ByCell+std_MovingAvg_WCG_ByCell,[1, shadeval/256, shadeval/256],'none');
alpha(transparentval)
hold on
p2=plot(bin_mdpt2MovingAvg,Mean_MovingAvg_WCG_ByCell,'-r','LineWidth',LNWidth);
xlim([-1500 1500])
ylim([0 ymaxval])
yticks([0 ymaxval/4*1 ymaxval/4*2 ymaxval/4*3 ymaxval/4*4])
xticks([-1000 -500  0  500 1000])
xlabs = xticklabels; xlabs(2) = {''}; xlabs(4) = {''}; xticklabels(xlabs);
ylabs = yticklabels; ylabs(2) = {''}; ylabs(4) = {''}; yticklabels(ylabs);
%ylabel("Detection Per Million Per H9 DNase Hypersenstivity Site")
%title("H9 DNaseI hypersensitivity sites")
title("HeLa DNaseI hypersensitivity sites")
%legend([p1,p2],["Accessibility";"Endogenous"],'location','southoutside','orientation','horizontal')
box off
grid off
legend boxoff
set(gcf,'color','w')
ylabel("Normalized reads")
xlabel("Distance from center of region (bp)")
legend([p1,p2],["Chromatin accessibility";"DNA methylation"],'location','northwest','orientation','vertical','FontSize',12)
legend boxoff
text(1000,ymaxval-0.004,{"\bfLower Quartile"},'FontName','Arial','FontSize',12,'HorizontalAlignment','center')



% print(gcf, 'H9_FreshCells_Normalized_ToBulk_Access_GC_and_Endogenous_WCG_Normalized_to_WCGGenome_H9_HypersensitvitySites_BottomQuartile_25BpBin_1bin_BothSides_MoveAvg_RPM.png', '-dpng', '-r600', '-painters')
% saveas(gcf, 'H9_FreshCells_Normalized_ToBulk_Access_GC_and_Endogenous_WCG_Normalized_to_WCGGenome_H9_HypersensitvitySites_BottomQuartile_25BpBin_1bin_BothSides_MoveAvg_RPM.pdf')