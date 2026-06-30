clc
clear

addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/ReferenceData');
addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/Hela');
addpath('/Users/piscopio/Desktop/STAMP/Experiments/scMATseq/Run67/genome/Run71/Hela/HSMatrixFiles');


%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_BottomQuartile.mat');
%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleBottomQuartile.mat');
%load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleTopQuartile.mat');
load('wgEncodeOpenChromDnaseHelas3Pk_Reformated_TopQuartile.mat');

%fileID= fopen('Hela-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_BottomQuartile-3000-Upstream_3000-Downstream.txt');
%SaveName='Bulk_GC_HeLa_Hypersensitivity_BottomQuartile_Per_Cell_25BP_Bin_PlusMinus3000.mat';
%fileID= fopen('Hela-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_MidBotQuartile-3000-Upstream_3000-Downstream.txt');
%SaveName='Bulk_GC_HeLa_Hypersensitivity_MiddleBottomQuartile_Per_Cell_25BP_Bin__PlusMinus3000.mat';
%fileID= fopen('Hela-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_MidTopQuartile-3000-Upstream_3000-Downstream.txt');
%SaveName='Bulk_GC_HeLa_Hypersensitivity_MiddleTopQuartile_Per_Cell_25BP_Bin__PlusMinus3000.mat';
fileID= fopen('Hela-MSPJI-se_Full5mC_Accessibilty_GC_Final_Rmdup_dnaseI-HS_TopQuartile-3000-Upstream_3000-Downstream.txt');
SaveName='Bulk_GC_HeLa_Hypersensitivity_TopQuartile_Per_Cell_25BP_Bin__PlusMinus3000.mat';

%%

WorkingCells = 1:96;


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
    
    CellGeneDistStruct(1).(fn{k}) = zeros(size(TSSstruct(1).(fn{k}),1),(Distance*2)/BinSize);
    
    CellGeneDistStruct(2).(fn{k}) = EdgesOfBins;
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
tempCell =1;

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


%%
figure
plot(xaxis,TSS_ForPloting_Pseudogene/length(WorkingCells))
% xlim([-Distance+MovingSize 2*Distance-MovingSize])
% xticklabels({'-2k';'';'TSS';'';'';'TES';'';'2k'})
xlabel('Relative Position')
ylabel('Normalized Detection Per Cell')
title('Bulk Fully Proteased H9 Genome GC Sites')
% print(gcf, 'Bulk_Fully_Proteased_H9_Genome_GC_Sitesat_TSS-TSE.png', '-dpng', '-r600', '-painters')

%To make psudeo genes each gene should be plus or minus a certain percentage of the gene.

TSSVecByGene = struct([]);

for k=1:numel(fn)
  
    TSSVecByChr(k,:)=sum(CellGeneDistStruct(1).(fn{k}));
    
end

bin_mdpt=(EdgesOfBins(2:end)+EdgesOfBins(1:(end-1)))/2;

figure
plot(bin_mdpt,sum(TSSVecByChr))

%%
save(SaveName,'-v7.3','BinSize','CellGeneDistStruct','EdgesOfBins','bin_mdpt','WorkingCells')

