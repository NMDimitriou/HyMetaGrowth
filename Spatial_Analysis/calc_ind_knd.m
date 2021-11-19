function varargout = calc_ind_knd(coord,name)
%% Computes the inter-nucleic distances and the nearest-neighbor distances.
% input: the coordinates of the nuclei centroids
% output: saved inter-nucleic, nearest-neighbor distances, and their
% distributions
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

dirINDist='Distances/INDist'; % directory to save the inter-nucleic distances
dirKNDist='Distances/KNDist'; % directory to save the nearest-neighbor distances
dirDistri='Distances/Distributions/'; % directory to save their distance distributions

if ~exist(dirINDist, 'dir')
    mkdir(dirINDist)
end

if ~exist(dirKNDist, 'dir')
    mkdir(dirKNDist)
end

if ~exist(dirDistri, 'dir')
    mkdir(dirDistri)
end

INDist = pdist(coord, 'euclidean'); % Inter nucleic distance 
[f1,xi1] = ksdensity(INDist);
fileID   = fopen([dirINDist '/' 'INDist_' name '.bin'],'w');
fwrite(fileID,INDist,'double');
fclose(fileID);
writematrix([f1;xi1],[dirDistri 'IND_' name]);


[~,NNDist] = knnsearch(coord, coord,'K',2); % nearest neighbor search
[f2,xi2]     = ksdensity(NNDist(:,2));
fileID       = fopen([dirKNDist '/' 'KNDist_' name '.bin'],'w');
fwrite(fileID,NNDist(:,2),'double');
fclose(fileID);
writematrix([f2;xi2],[dirDistri '/' 'KND_' name]);


