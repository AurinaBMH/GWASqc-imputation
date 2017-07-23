% load .het data

DATA = readtable('rawGWAdata_rsex_filteredIDcontrols.txt');

%threshold = 0.03;
Nnm = DATA.N_NM_;
Ohom = DATA.O_HOM_;
%Calculate the observed heterozygosity rate per individual
Hzrate = (Nnm-Ohom)./Nnm;
% select to remove individuals with missingness rate higher than a threshold
%ind1 = (F>=threshold);
%select to remove individuals with +-3SD from the mean on heterozigosity
%rate
ind2 = (Hzrate>mean(Hzrate)+3*std(Hzrate));
ind3 = (Hzrate<mean(Hzrate)-3*std(Hzrate));
% combine indivicuals to remove
ind = ind2+ind3;

% select FID and IID for individuals intended to remove
FID = DATA.FID(ind>0);
IID = DATA.IID(ind>0);

%save data into a table
T = table(FID,IID);
writetable(T,'fail-imisshet-qc.txt', 'Delimiter','\t'); 
