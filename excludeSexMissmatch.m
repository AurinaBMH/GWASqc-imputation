% read txt file for sexcheck
DATA = readtable('rawGWAdata.txt');
k=1; 
% select subjects that have missmatching information and also have initial
% sex information
for i=1:size(DATA,1)
    phenoS = DATA.PEDSEX(i); 
    geneS = DATA.SNPSEX(i); 
    if (phenoS==1 && geneS==2) || (phenoS==2 && geneS==1)
        FID(k) = DATA.FID(i); 
        IID(k) = DATA.IID(i); 
        k=k+1; 
    end
end

% save their family ID (FID) and individual ID (IID) into a txt file to
% exclude
T = table(FID',IID');
writetable(T,'fail-sex-qc.txt', 'Delimiter','\t', ...
    'WriteVariableNames', false, 'WriteRowNames', false); 
