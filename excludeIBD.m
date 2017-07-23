%-------------------------------------------------------------------------------
% Aurina Arnatkeviciute 01-06-2017
% This script finds family ID and individual ID for subjects that don't
% pass relatedness test in GWAS QC. One subject out of a pair is selected
% by taking one column of IDs. 
%-------------------------------------------------------------------------------
% create FID and IID list for subjects that have high IBD rates
%open a list of IDs for selected related individuals
fid = fopen('rawGWAdata.genome.PIHat.sorted.txt');
fseek(fid, 0, 'eof');
chunksize = ftell(fid);
fseek(fid, 0, 'bof');
ch = fread(fid, chunksize, '*uchar');
nol = sum(ch == sprintf('\n')); % number of lines
fid = fclose(fid);

fid = fopen('rawGWAdata.genome.PIHat.sorted.txt'); 
data = textscan(fid, '%d %d %*[^\n]', nol-1);
fid = fclose(fid);

% select subjects in the first column as subjects to remove from the
% analysis
IID = data{1}; 

% open raw .fam file where all subject IDs and family IDs are stored.
fid = fopen('rawGWAdata.fam'); 
allIDs = textscan(fid, '%d %d %*[^\n]');
fid = fclose(fid);

IIDall = allIDs{2}; 
FIDall = allIDs{1}; 

% find the indexes for subjects that we want to exclude and save their
% family IDs. 

[~, ind] = intersect(IIDall, IID); 
% save thei FIDs
FID = FIDall(ind); 

% save results into the table
T = table(FID, IID); 
writetable(T,'fail-BID-QC.txt', 'Delimiter','\t'); 
