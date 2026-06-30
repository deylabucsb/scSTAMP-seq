clc
clear



addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Exp_02/Run56_Run61_Run62/P4L1/Genomic/ReferenceData/ReferenceMatrixFiles')
%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_BottomQuartile.mat');
%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleBottomQuartile.mat');
%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleTopQuartile.mat');
load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_TopQuartile.mat');

% (20230113) (20230116) (20230117)
%addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/HypersensitivityData/');
%CountsPerCell = dlmread('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_CountsPerCellFaba.txt');
%CountsPerCell = dlmread('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_CountsPerCellFaba.txt');

% (20230130)
%__________R REFORMATTED FABA.TXT FILES TO GET RID OF LOW COUNT CELLS AND U2OS
addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Exp02_(Run56_Run61)/Genomic/HypersensitivityData/');

%CountsPerCell = dlmread('Accessibility_GC_Final_Rndup_dnaseI-HS_CountsPerCellFaba_R_Reformatted_Cluster_0.txt');
%CountsPerCell = dlmread('Accessibility_GC_Final_Rndup_dnaseI-HS_CountsPerCellFaba_R_Reformatted_Cluster_1.txt');

%CountsPerCell = dlmread('Endogenous_WCG_Final_Rndup_dnaseI-HS_CountsPerCellFaba_R_Reformatted_Cluster_0.txt');
CountsPerCell = dlmread('Endogenous_WCG_Final_Rndup_dnaseI-HS_CountsPerCellFaba_R_Reformatted_Cluster_1.txt');

%__________SAVE ACCESSIBILITY.MAT FILES TO DIRECTORY
% Cluster 0
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_Bottom_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_BottomMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_MiddleBottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_TopMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_MiddleTopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_Top_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_TopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';

% Cluster 1
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_Bottom_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_BottomMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_MiddleBottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_TopMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_MiddleTopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_Top_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ACCESS_GC_2022-07-19_HeLa_Hypersenstivity_TopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';

%____________________SAVE ENDOGENOUS.MAT FILES TO DIRECTORY
% Cluster 0
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_Bottom_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_BottomMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_MiddleBottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_TopMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_MiddleTopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_Top_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_0_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_TopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';

% Cluster 1
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_Bottom_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_BottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_BottomMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_MiddleBottomQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
%fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_TopMiddle_quartile-3000-Upstream_3000-Downstream.txt');
%SaveName='R_Reformatted_Cluster_1_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_MiddleTopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';
fileID= fopen('P4L1_G-MSPJI-se_Full5mC_Endogenous_WCG_Final_Rmdup_dnaseI-HS_Top_quartile-3000-Upstream_3000-Downstream.txt');
SaveName='R_Reformatted_Cluster_1_P4L1_ENDOG_WCG_2022-07-19_HeLa_Hypersenstivity_TopQuartile_Per_Cell_25BP-Bin_plusMinus2000.mat';


%WorkingCells = find(CountsPerCell>1000);
WorkingCells = find(CountsPerCell>10000);
%WorkingCells = find(CountsPerCell>30000);
%WorkingCells = find(CountsPerCell>50000);

%%

% addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\ReferenceDataUsedForPaper\H9ActiveEnhancer')
% load('wgEncodeAwgDnaseDukeH9esUniPk_Reformated_TopQuartile.mat');
% load('wgEncodeAwgDnaseDukeH9esUniPk_Reformated_MiddleTopQuartile.mat');
% load('wgEncodeAwgDnaseDukeH9esUniPk_Reformated_MiddleBottomQuartile.mat');
% load('wgEncodeAwgDnaseDukeH9esUniPk_Reformated_BottomQuartile.mat');

% addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib3\Lib3_AllFaba');
% addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\P39_H9_scMATseq\Lib4\Lib4_AllFaba');
% addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\P38_Frozen_H9_scMATseq\P38_Lib1\Faba_P38Lib1');
% addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\P38_Frozen_H9_scMATseq\P38_Lib2\Faba_P38Lib2');


%addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\scMAT-seq_Figures\H9_HypersenstivitySites\H9\FrozenP38');
% CountsPerCell = dlmread('P39_5mC_Access_mRNA_Fresh_H9_Lib4_AllSeqRun_Full5mC_Accessibilty_GC_Fin_CountsPerCellFaba.txt');
% CountsPerCell = dlmread('P39_5mC_Access_mRNA_Fresh_H9_Lib4_AllSeqRun_Full5mC_Endogenous_WCG_Fin_CountsPerCellFaba.txt');
% CountsPerCell = dlmread('P39_5mC_Access_mRNA_Fresh_H9_Lib3_24_Cells_A_AllSeqRun_Full5mC_Endogenous_WCG_Fin_CountsPerCellFaba.txt');
% CountsPerCell = dlmread('P39_5mC_Access_mRNA_Fresh_H9_Lib3_24_Cells_A_AllSeqRun_Full5mC_Accessibilty_GC_Fin_CountsPerCellFaba.txt');
% CountsPerCell = dlmread('P38_Frozen_H9_scMAT_Lib2_AllSeqRun_Full5mC_Accessibilty_GC_Fin_CountsPerCellFaba.txt');

% WorkingCells = find(CountsPerCell>10000);
% WorkingCells = find(CountsPerCell>30000);


%fileID= fopen('P38_Frozen_H9_scMAT_Lib2_AllSeqRun_Full5mC_Accessibilty_GC_Final_Rmdup_H9HyperSensitivitySites-2000-Upstream_2000-Downstream.txt');
%SaveName='P38_Lib2_GC_2020_6_06_H9_Hypersenstivity_BottomQuartile_Per_Cell_25BP_Bin_plusMinus2000.mat';


C = textscan(fileID,'%f %s %f %f','HeaderLines',0);
%C = textscan(fileID,'%f %s %f %s %f','HeaderLines',0);
fclose(fileID);

SingleCellChr = C{1,2};
SingleCellLoci = C{1,3};
SingleCellCell= C{1,1};
SingleCellStrand= C{1,4};



%%
%Filter Out Loci that are not in range and then cells that were poorly sequenced
SingleCellCell_temp = SingleCellCell;
SingleCellChr_temp = SingleCellChr;
SingleCellLoci_temp = SingleCellLoci;
SingleCellStrand_temp = SingleCellStrand;

keepRows2=zeros(length(SingleCellCell_temp),1);

for i=1:1:length(SingleCellCell_temp)
    tempsum = sum(SingleCellCell_temp(i) == WorkingCells);
    if tempsum == 1
        keepRows2(i,1)=i;
    end
end
keepRows2(keepRows2==0)=[];


SingleCellCell_Filtered_InRightLocation = SingleCellCell_temp(keepRows2);
SingleCellChr_Filtered_InRightLocation = SingleCellChr_temp(keepRows2);
SingleCellLoci_Filtered_InRightLocation = SingleCellLoci_temp(keepRows2);
SingleCellStrand_Filtered_InRightLocation = SingleCellStrand_temp(keepRows2);


clear SingleCellCell SingleCellChr SingleCellLoci SingleCellStrand
clear SingleCellCell_temp SingleCellChr_temp SingleCellLoci_temp SingleCellStrand_temp
clear keepRows keepRows2
%%
%Create Downstream of TSS by Distance Value and store in Slot 8 of TSSstruct
%Create a new structure that will contain binned Data around each gene for
%each cell.  Chr by Cell structure, inside each is a Bin by gene format

fn = fieldnames(TSSstruct);
CellGeneDistStruct = struct([]);
BinSize = 25;
EdgesOfBins=-Distance:BinSize:Distance;

for k=1:numel(fn)
    TSSstruct(8).(fn{k}) = TSSstruct(1).(fn{k}) + Distance*TSSstruct(2).(fn{k});
    
    for i=1:1:96
       if(find(WorkingCells == i))
        CellGeneDistStruct(i).(fn{k}) = zeros(size(TSSstruct(1).(fn{k}),1),(Distance*2)/BinSize);
       end 
    end
    
    CellGeneDistStruct(97).(fn{k}) = EdgesOfBins;
end

%%
% Want to report % of the way to each element aka Promoter start to TSS,
% TSS to TSE, and TSE to Downstream TSE
Prom_to_TSS = zeros(length(SingleCellLoci_Filtered_InRightLocation),1);
Prom_to_TSS(Prom_to_TSS==0)=-10000;
TSS_to_DownStream = zeros(length(SingleCellLoci_Filtered_InRightLocation),1);
TSS_to_DownStream(TSS_to_DownStream==0)=-10000;


% PlusStrandTSSstruct
% 
% find(TSSstruct(2).(char(SingleCellChr_Filtered_InRightLocation(i)))==-1)
%1 is TSS, 2 is Strand, 3 is Promoter, 4 is Downstream TSE, 7 is TSE
P_TSS_Counter = 1;
TSS_DS_Counter = 1;

tic
for i=1:1:length(SingleCellLoci_Filtered_InRightLocation)
    %
temploc =SingleCellLoci_Filtered_InRightLocation(i);
tempchr = SingleCellChr_Filtered_InRightLocation(i);
tempCell =SingleCellCell_Filtered_InRightLocation(i);

%promoter to TSS for plus strand
b=temploc>=TSSstruct(3).(char(tempchr)) & temploc<TSSstruct(1).(char(tempchr)) & TSSstruct(2).(char(tempchr)) == 1;
b1 = find(b >0);
if ~isempty(b1)
   for j=1:1:length(b1)
        TempDistToTSS = temploc - TSSstruct(1).(char(tempchr))(b1(j)); %Find DistanceToTSS for each gene
        TempBin = discretize(TempDistToTSS,EdgesOfBins); %Find bin it will be in for the Distance To TSS for each gene
        if ~isnan(TempBin)
        CellGeneDistStruct(tempCell).(char(tempchr))(b1(j),TempBin)=CellGeneDistStruct(tempCell).(char(tempchr))(b1(j),TempBin)+1; %+1 to that bin value for that gene
        
        Prom_to_TSS(P_TSS_Counter:P_TSS_Counter+1,1) = TempDistToTSS;
        P_TSS_Counter=P_TSS_Counter+1;        
        end
    end
end

%TSS to Downstream for plus strand
c=temploc>=TSSstruct(1).(char(tempchr)) & temploc<=TSSstruct(8).(char(tempchr)) & TSSstruct(2).(char(tempchr)) == 1;
c1 = find(c >0);
if ~isempty(c1)
    for j=1:1:length(c1)
        TempDistToTSS = temploc - TSSstruct(1).(char(tempchr))(c1(j));
        TempBin = discretize(TempDistToTSS,EdgesOfBins); 
        if ~isnan(TempBin)
        CellGeneDistStruct(tempCell).(char(tempchr))(c1(j),TempBin)=CellGeneDistStruct(tempCell).(char(tempchr))(c1(j),TempBin)+1; 
        
        TSS_to_DownStream(TSS_DS_Counter:TSS_DS_Counter+1,1) = TempDistToTSS;
        TSS_DS_Counter=TSS_DS_Counter+1;
        end
    end
end


%---%

%promoter to TSS for minus strand
bm=temploc<=TSSstruct(3).(char(tempchr)) & temploc>TSSstruct(1).(char(tempchr)) & TSSstruct(2).(char(tempchr)) == -1;
bm1 = find(bm >0);
if ~isempty(bm1)
   for j=1:1:length(bm1)
        TempDistToTSS = TSSstruct(1).(char(tempchr))(bm1(j)) - temploc;
        TempBin = discretize(TempDistToTSS,EdgesOfBins); 
        if ~isnan(TempBin)
        CellGeneDistStruct(tempCell).(char(tempchr))(bm1(j),TempBin)=CellGeneDistStruct(tempCell).(char(tempchr))(bm1(j),TempBin)+1; 
        
        
        Prom_to_TSS(P_TSS_Counter:P_TSS_Counter+1,1) = TempDistToTSS;
        P_TSS_Counter=P_TSS_Counter+1;
        end
    end
end

%TSS to Downstream for minus strand
cm=temploc<=TSSstruct(1).(char(tempchr)) & temploc>=TSSstruct(8).(char(tempchr)) & TSSstruct(2).(char(tempchr)) == -1;
cm1 = find(cm >0);
if ~isempty(cm1)
    for j=1:1:length(cm1)
        TempDistToTSS = TSSstruct(1).(char(tempchr))(cm1(j)) - temploc;
        TempBin = discretize(TempDistToTSS,EdgesOfBins); 
        if ~isnan(TempBin)
        CellGeneDistStruct(tempCell).(char(tempchr))(cm1(j),TempBin)=CellGeneDistStruct(tempCell).(char(tempchr))(cm1(j),TempBin)+1; 
        
        
        TSS_to_DownStream(TSS_DS_Counter:TSS_DS_Counter+1,1) = TempDistToTSS;
        TSS_DS_Counter=TSS_DS_Counter+1;
        end
    end
end

end
toc


Prom_to_TSS(Prom_to_TSS==-10000)=[];
TSS_to_DownStream(TSS_to_DownStream==-10000)=[];

%%
% save('StrippedDNA_H9_GC_ACCESS_Bulk_AllData_2019_11_23_Genes.mat','keepRows','keepRows2','Prom_to_TSS','TSS_to_TSE','TSE_to_Downstream');
%  save('Human_GC_Mappable_2019_11_11.mat','keepRows','keepRows2','Prom_to_TSS','TSS_to_TSE','TSE_to_Downstream');

%%
%I will create a psuedogene of 3000 length, need to normalize by number of
%bp looked at  
clear TSS_ForPloting_Pseudogene xaxis
MovingSize=25;
WindowSize=25;
counts=1;


for i=-Distance:MovingSize:Distance
    xaxis(counts,1)=i;
     
    TSS_ForPloting_Pseudogene(counts,1) = sum(Prom_to_TSS<i+WindowSize/2 & Prom_to_TSS>=i-WindowSize/2)+sum(TSS_to_DownStream<i+WindowSize/2 & TSS_to_DownStream>=i-WindowSize/2);
 
    counts = counts+1;    
end


%% COMMENT IF COMBINED PLOTS DESIRED
for k=1:numel(fn)
  for i=1:1:length(WorkingCells)
       j=WorkingCells(i);
       if i==1
        TSSVecByGene(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
     
       else
       TSSVecByGene(1).(fn{k}) = TSSVecByGene(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
%        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
       end
%     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
  end 
  TSSVecByChr(k,:)=sum(TSSVecByGene(1).(fn{k}));
end

bin_mdpt=(EdgesOfBins(2:end)+EdgesOfBins(1:(end-1)))/2;

%% UNCOMMENT IF COMBINED PLOTS DESIRED
% figure
% plot(xaxis,TSS_ForPloting_Pseudogene/length(WorkingCells))
% % xlim([-Distance+MovingSize 2*Distance-MovingSize])
% % xticklabels({'-2k';'';'TSS';'';'';'TES';'';'2k'})
% xlabel('Relative Position')
% ylabel('Normalized Detection Per Cell')
% title('Bulk Fully Proteased HeLa Genome GC Sites')
% % print(gcf, 'Bulk_Fully_Proteased_H9_Genome_GC_Sitesat_TSS-TSE.png', '-dpng', '-r600', '-painters')
% 
% %To make psudeo genes each gene should be plus or minus a certain percentage of the gene.
% 
% TSSVecByGene = struct([]);
% 
% % for k=1:numel(fn)
% %   for i=1:1:length(WorkingCells)
% %        j=WorkingCells(i);
% %        if i==1
% %         TSSVecByGene(1).(fn{k}) = CellGeneDistStruct(j).(fn{k});
% %      
% %        else
% %        TSSVecByGene(1).(fn{k}) = TSSVecByGene(1).(fn{k}) + CellGeneDistStruct(j).(fn{k});
% % %        TSSVecByCell(i,:)=sum(CellGeneDistStruct(j).(fn{k}));  
% %        end
% % %     TSSVecByChr(k,:)=sum(CellGeneDistStruct(i).(fn{k}));
% %   end 
% %   TSSVecByChr(k,:)=sum(TSSVecByGene(1).(fn{k}));
% % end
% 
% bin_mdpt=(EdgesOfBins(2:end)+EdgesOfBins(1:(end-1)))/2;
% 
% figure
% plot(bin_mdpt,sum(TSSVecByChr))

%% COMMENT OUT IF INDIVIDUAL PLOTS DESIRED

figure
tiledlayout(1,2)
nexttile;
plot(xaxis,TSS_ForPloting_Pseudogene/length(WorkingCells))
% xlim([-Distance+MovingSize 2*Distance-MovingSize])
% xticklabels({'-2k';'';'TSS';'';'';'TES';'';'2k'})
xlabel('Relative Position')
ylabel('Normalized Detection Per Cell')
title('Bulk Fully Proteased HeLa Genome GC Sites')
nexttile;
plot(bin_mdpt,sum(TSSVecByChr))
title('')
xlabel('')
ylabel('')
set(gcf,'Position',[0, 0, 800, 400])


%%
save(SaveName,'-v7.3','BinSize','CellGeneDistStruct','EdgesOfBins','bin_mdpt','WorkingCells','CountsPerCell')

