fid = fopen('rawGWAdataRM.sexprobs'); 
data = textscan(fid, '%d %d %d %*[^\n]');
fid = fclose(fid);

% find subjects with 0 in the third column and exclude them from the list,
% therefore keep only those subjects that had initial information about sex
% and it doesn't match with the original.
keep = find(data{3}>0); 
FID = data{1}(keep); 
IID = data{2}(keep); 

T = table(FID, IID); 
writetable(T,'fail-sexcheck-qc.txt', 'Delimiter','\t'); 
