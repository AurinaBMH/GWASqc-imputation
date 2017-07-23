
%check the number of subjects in each HapMap group in order to fix MDS
%plotting and legends
DATA = readtable('rawGWAdata_HM3_mds.txt');
group = DATA.FID; 


nASW = sum(strcmp(group, 'ASW')); 
nCEU = sum(strcmp(group, 'CEU')); 
nCHB = sum(strcmp(group, 'CHB')); 
nCHD = sum(strcmp(group, 'CHD')); 
nGIH = sum(strcmp(group, 'GIH')); 
nJPT = sum(strcmp(group, 'JPT')); 
nLWK = sum(strcmp(group, 'LWK')); 
nMEX = sum(strcmp(group, 'MEX')); 
nMKK = sum(strcmp(group, 'MKK')); 
nTSI = sum(strcmp(group, 'TSI')); 
nYRI = sum(strcmp(group, 'YRI')); 

