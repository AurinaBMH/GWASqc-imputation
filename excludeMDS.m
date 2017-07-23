% choose a threshold for subject exclusion based on standart deviations on
% their first two components
threshold = 2; % +-2DS
% exclude subjects with +-3SD on 1 and 2 components
mixedDATA = readtable('rawGWAdata_HM3_mds.txt');
% keep only those subjects from our studies
ourDATA = readtable('rawGWAdata.txt');

% get IDs to filter
IIDall = mixedDATA.IID;
IIDour = ourDATA.IID;
% select C1 and C2 just for subjects in our studies (exclude HM3 samples)
[~, indkeep] = intersect(IIDall, IIDour);
C1 = mixedDATA.C1(indkeep);
C2 = mixedDATA.C2(indkeep);

% select ind based on the criteria
excludeC11 = find(C1>(mean(C1)+threshold*(std(C1))));
excludeC12 = find(C1<(mean(C1)-threshold*(std(C1))));

excludeC21 = find(C2>(mean(C2)+threshold*(std(C2))));
excludeC22 = find(C2<(mean(C2)-threshold*(std(C2))));

% combine all indexes and choose the unique ones
excludeall = unique([excludeC11; excludeC12; excludeC21; excludeC22]);

% get their famili IDs and individual IDs from the main file
FID = mixedDATA.FID(excludeall);
IID = mixedDATA.IID(excludeall);

% put them into a table
T = table(FID, IID);
% save to file
writetable(T,'fail-mds-qc.txt', 'Delimiter','\t', ...
    'WriteVariableNames', false, 'WriteRowNames', false);
