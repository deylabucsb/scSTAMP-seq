
clc
clear all

%%

% READ IN ACTIVE ENHANCER SITES OR JUST LOAD IT FROM FIRST TIME I READ IT IN

%addpath('C:\Users\alexc\Box Sync\Research\scMAT-seq\ReferenceDataUsedForPaper\CTCF\Ensemble75');
%fileID = fopen('Ensembel75_CTCF_BindingMotifs_GRCh37.p13_Reformat_TopQuartile.txt');

%addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/ReferenceData/');
addpath('/Users/piscopio/Library/CloudStorage/Box-Box/scMATseq/Run56_Run61_Run62/P4L1/Genomic/ReferenceData');

fileID = fopen('wgEncodeOpenChromDnaseHelas3Pk_Reformated_BottomQuartile.txt');
%fileID = fopen('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleBottomQuartile.txt');
%fileID = fopen('wgEncodeOpenChromDnaseHelas3Pk_Reformated_MiddleTopQuartile.txt');
%fileID = fopen('wgEncodeOpenChromDnaseHelas3Pk_Reformated_TopQuartile.txt');




%wgEncodeAwgDnaseDukeH9esUniPk_Reformated_

% C = textscan(fileID,'%f %f %f %f %s','HeaderLines',0);
C = textscan(fileID,'%f %f %f %f %s %f','HeaderLines',0);
fclose(fileID);

ListGenes = C{1,5};

%TSS goes Chr (just the number), Loci, Strand
    TSS(:,1) = C{1,1};
    TSS(:,2) = C{1,2};
    TSS(:,3) = C{1,4};
% load('H9_Active_Enhancer_Sites.mat'); %Load in active enhancer sites and nameclear




%%
%plus strand subtract to go upstream.  Minus strand add to go upstream
Distance=2000;
MAvDist=100;
Distance1=Distance+MAvDist/2;
%create structure
field1 = 'chr1';  value1 = TSS(find(TSS(:,1)==1),2); strand1 = TSS(find(TSS(:,1)==1),3); numb1= zeros(length(strand1),1); promot1=value1-Distance1*strand1; DownStream1=value1+Distance1*strand1;
field2 = 'chr2';  value2 = TSS(find(TSS(:,1)==2),2); strand2 = TSS(find(TSS(:,1)==2),3); numb2= zeros(length(strand2),1); promot2=value2-Distance1*strand2; DownStream2=value2+Distance1*strand2;
field3 = 'chr3';  value3 = TSS(find(TSS(:,1)==3),2); strand3 = TSS(find(TSS(:,1)==3),3); numb3= zeros(length(strand3),1); promot3=value3-Distance1*strand3; DownStream3=value3+Distance1*strand3;
field4 = 'chr4';  value4 = TSS(find(TSS(:,1)==4),2); strand4 = TSS(find(TSS(:,1)==4),3); numb4= zeros(length(strand4),1); promot4=value4-Distance1*strand4; DownStream4=value4+Distance1*strand4;
field5 = 'chr5';  value5 = TSS(find(TSS(:,1)==5),2); strand5 = TSS(find(TSS(:,1)==5),3); numb5= zeros(length(strand5),1); promot5=value5-Distance1*strand5; DownStream5=value5+Distance1*strand5;
field6 = 'chr6';  value6 = TSS(find(TSS(:,1)==6),2); strand6 = TSS(find(TSS(:,1)==6),3); numb6= zeros(length(strand6),1); promot6=value6-Distance1*strand6; DownStream6=value6+Distance1*strand6;
field7 = 'chr7';  value7 = TSS(find(TSS(:,1)==7),2); strand7 = TSS(find(TSS(:,1)==7),3); numb7= zeros(length(strand7),1); promot7=value7-Distance1*strand7; DownStream7=value7+Distance1*strand7;
field8 = 'chr8';  value8 = TSS(find(TSS(:,1)==8),2); strand8 = TSS(find(TSS(:,1)==8),3); numb8= zeros(length(strand8),1); promot8=value8-Distance1*strand8; DownStream8=value8+Distance1*strand8;
field9 = 'chr9';  value9 = TSS(find(TSS(:,1)==9),2); strand9 = TSS(find(TSS(:,1)==9),3); numb9= zeros(length(strand9),1); promot9=value9-Distance1*strand9; DownStream9=value9+Distance1*strand9;
field10 = 'chr10';  value10 = TSS(find(TSS(:,1)==10),2); strand10 = TSS(find(TSS(:,1)==10),3); numb10= zeros(length(strand10),1); promot10=value10-Distance1*strand10; DownStream10=value10+Distance1*strand10;
field11 = 'chr11';  value11 = TSS(find(TSS(:,1)==11),2); strand11 = TSS(find(TSS(:,1)==11),3); numb11= zeros(length(strand11),1); promot11=value11-Distance1*strand11; DownStream11=value11+Distance1*strand11;
field12 = 'chr12';  value12 = TSS(find(TSS(:,1)==12),2); strand12 = TSS(find(TSS(:,1)==12),3); numb12= zeros(length(strand12),1); promot12=value12-Distance1*strand12; DownStream12=value12+Distance1*strand12;
field13 = 'chr13';  value13 = TSS(find(TSS(:,1)==13),2); strand13 = TSS(find(TSS(:,1)==13),3); numb13= zeros(length(strand13),1); promot13=value13-Distance1*strand13; DownStream13=value13+Distance1*strand13;
field14 = 'chr14';  value14 = TSS(find(TSS(:,1)==14),2); strand14 = TSS(find(TSS(:,1)==14),3); numb14= zeros(length(strand14),1); promot14=value14-Distance1*strand14; DownStream14=value14+Distance1*strand14;
field15 = 'chr15';  value15 = TSS(find(TSS(:,1)==15),2); strand15 = TSS(find(TSS(:,1)==15),3); numb15= zeros(length(strand15),1); promot15=value15-Distance1*strand15; DownStream15=value15+Distance1*strand15;
field16 = 'chr16';  value16 = TSS(find(TSS(:,1)==16),2); strand16 = TSS(find(TSS(:,1)==16),3); numb16= zeros(length(strand16),1); promot16=value16-Distance1*strand16; DownStream16=value16+Distance1*strand16;
field17 = 'chr17';  value17 = TSS(find(TSS(:,1)==17),2); strand17 = TSS(find(TSS(:,1)==17),3); numb17= zeros(length(strand17),1); promot17=value17-Distance1*strand17; DownStream17=value17+Distance1*strand17;
field18 = 'chr18';  value18 = TSS(find(TSS(:,1)==18),2); strand18 = TSS(find(TSS(:,1)==18),3); numb18= zeros(length(strand18),1); promot18=value18-Distance1*strand18; DownStream18=value18+Distance1*strand18;
field19 = 'chr19';  value19 = TSS(find(TSS(:,1)==19),2); strand19 = TSS(find(TSS(:,1)==19),3); numb19= zeros(length(strand19),1); promot19=value19-Distance1*strand19; DownStream19=value19+Distance1*strand19;
field20 = 'chr20';  value20 = TSS(find(TSS(:,1)==20),2); strand20 = TSS(find(TSS(:,1)==20),3); numb20= zeros(length(strand20),1); promot20=value20-Distance1*strand20; DownStream20=value20+Distance1*strand20;
field21 = 'chr21';  value21 = TSS(find(TSS(:,1)==21),2); strand21 = TSS(find(TSS(:,1)==21),3); numb21= zeros(length(strand21),1); promot21=value21-Distance1*strand21; DownStream21=value21+Distance1*strand21;
field22 = 'chr22';  value22 = TSS(find(TSS(:,1)==22),2); strand22 = TSS(find(TSS(:,1)==22),3); numb22= zeros(length(strand22),1); promot22=value22-Distance1*strand22; DownStream22=value22+Distance1*strand22;
field23 = 'chr23';  value23 = TSS(find(TSS(:,1)==23),2); strand23 = TSS(find(TSS(:,1)==23),3); numb23= zeros(length(strand23),1); promot23=value23-Distance1*strand23; DownStream23=value23+Distance1*strand23;
field24 = 'chr24';  value24 = TSS(find(TSS(:,1)==24),2); strand24 = TSS(find(TSS(:,1)==24),3); numb24= zeros(length(strand24),1); promot24=value24-Distance1*strand24; DownStream24=value24+Distance1*strand24;
field25 = 'chr25';  value25 = TSS(find(TSS(:,1)==25),2); strand25 = TSS(find(TSS(:,1)==25),3); numb25= zeros(length(strand25),1); promot25=value25-Distance1*strand25; DownStream25=value25+Distance1*strand25;

% TSSstruct = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
%     field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,...
%     field10,value10,field11,value11,field12,value12,field13,value13,field14,value14,...
%     field15,value15,field16,value16,field17,value17,field18,value18,field19,value19,...
%     field20,value20,field21,value21,field22,value22,field23,value23,field24,value24,...
%     field25,value25);
GenesNameValue1 = ListGenes(find(TSS(:,1)==1),1); 
GenesNameValue2 = ListGenes(find(TSS(:,1)==2),1);
GenesNameValue3 = ListGenes(find(TSS(:,1)==3),1);
GenesNameValue4 = ListGenes(find(TSS(:,1)==4),1);
GenesNameValue5 = ListGenes(find(TSS(:,1)==5),1);
GenesNameValue6 = ListGenes(find(TSS(:,1)==6),1);
GenesNameValue7 = ListGenes(find(TSS(:,1)==7),1);
GenesNameValue8 = ListGenes(find(TSS(:,1)==8),1);
GenesNameValue9 = ListGenes(find(TSS(:,1)==9),1);
GenesNameValue10 = ListGenes(find(TSS(:,1)==10),1); 
GenesNameValue11 = ListGenes(find(TSS(:,1)==11),1); 
GenesNameValue12 = ListGenes(find(TSS(:,1)==12),1); 
GenesNameValue13 = ListGenes(find(TSS(:,1)==13),1); 
GenesNameValue14 = ListGenes(find(TSS(:,1)==14),1); 
GenesNameValue15 = ListGenes(find(TSS(:,1)==15),1); 
GenesNameValue16 = ListGenes(find(TSS(:,1)==16),1); 
GenesNameValue17 = ListGenes(find(TSS(:,1)==17),1); 
GenesNameValue18 = ListGenes(find(TSS(:,1)==18),1); 
GenesNameValue19 = ListGenes(find(TSS(:,1)==19),1); 
GenesNameValue20 = ListGenes(find(TSS(:,1)==20),1); 
GenesNameValue21 = ListGenes(find(TSS(:,1)==21),1); 
GenesNameValue22 = ListGenes(find(TSS(:,1)==22),1); 
GenesNameValue23 = ListGenes(find(TSS(:,1)==23),1); 
GenesNameValue24 = ListGenes(find(TSS(:,1)==24),1); 
GenesNameValue25 = ListGenes(find(TSS(:,1)==25),1); 

TSSstruct = struct(...
    field1,{value1;strand1;promot1;DownStream1;numb1;GenesNameValue1},...
    field2,{value2;strand2;promot2;DownStream2;numb2;GenesNameValue2},...
    field3,{value3;strand3;promot3;DownStream3;numb3;GenesNameValue3},...
    field4,{value4;strand4;promot4;DownStream4;numb4;GenesNameValue4},...
    field5,{value5;strand5;promot5;DownStream5;numb5;GenesNameValue5},...
    field6,{value6;strand6;promot6;DownStream6;numb6;GenesNameValue6},...
    field7,{value7;strand7;promot7;DownStream7;numb7;GenesNameValue7},...
    field8,{value8;strand8;promot8;DownStream8;numb8;GenesNameValue8},...
    field9,{value9;strand9;promot9;DownStream9;numb9;GenesNameValue9},...
    field10,{value10;strand10;promot10;DownStream10;numb10;GenesNameValue10},...
    field11,{value11;strand11;promot11;DownStream11;numb11;GenesNameValue11},...
    field12,{value12;strand12;promot12;DownStream12;numb12;GenesNameValue12},...
    field13,{value13;strand13;promot13;DownStream13;numb13;GenesNameValue13},...
    field14,{value14;strand14;promot14;DownStream14;numb14;GenesNameValue14},...
    field15,{value15;strand15;promot15;DownStream15;numb15;GenesNameValue15},...
    field16,{value16;strand16;promot16;DownStream16;numb16;GenesNameValue16},...
    field17,{value17;strand17;promot17;DownStream17;numb17;GenesNameValue17},...
    field18,{value18;strand18;promot18;DownStream18;numb18;GenesNameValue18},...
    field19,{value19;strand19;promot19;DownStream19;numb19;GenesNameValue19},...
    field20,{value20;strand20;promot20;DownStream20;numb20;GenesNameValue20},...
    field21,{value21;strand21;promot21;DownStream21;numb21;GenesNameValue21},...
    field22,{value22;strand22;promot22;DownStream22;numb22;GenesNameValue22},...
    field23,{value23;strand23;promot23;DownStream23;numb23;GenesNameValue23},...
    field24,{value24;strand24;promot24;DownStream24;numb24;GenesNameValue24},...
    field25,{value25;strand25;promot25;DownStream25;numb25;GenesNameValue25});


% 
clear field* value* strand* promot* DownStream* numb* GenesNameValue*
% TSSstt = struct(field1,{value1;strand1},field2,{value2;strand2}
%%
%read in possible locations
% % fileID = fopen('AluI_human.txt');
% % D = textscan(fileID,'%s %s %f','HeaderLines',0);
% % fclose(fileID);
% % AGCTChr=D{1,2};
% % AGCT(:,1) = D{1,3};
%%
ChrNames=fieldnames(TSSstruct);
