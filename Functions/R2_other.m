function [R2_other_dir,SIstd,SIall]=R2_other(score,idx_dir,idx_duration,ndim,normalised_t)
%% R2_other calculates the coefficient of determination between trajectories
%% of different directions
%
% INPUTS
%
% score: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score
% matrix
% 
% idx_duration: array containing the duration bin of each row in the score
% matrix
%
% ndim: number of dimensions of the trajectories
%
% normalised_t: vector of timebins to normalise trajectories 
% 
% OUTPUTS
%
% R2_other_dir: mean coefficient of determination (R2) across all trajectories
%
% SIstd: standard error of the mean R2 across all trajectories
%
% SIall: array of coefficient of determination (R2) across all trajectories
%
% 24/05/2023
% Andrea Colins Rodriguez


Ndir=max(idx_dir);
Nbins=max(idx_duration);
R2_other_dir=zeros(Ndir,Nbins-1);

for i_dir=1:Ndir
    counter=1;
    for j_dir=i_dir+1:Ndir
        
        for i_bin=1:Nbins
            idx=idx_dir==i_dir & idx_duration==i_bin;
            
            for j_bin=i_bin+1:Nbins
                idx_j=idx_dir==j_dir & idx_duration==j_bin;
                
                R2_other_dir(i_dir,counter)=scaled_TR_R2(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
                counter=counter+1;
            end
        end
    end
end

SIall=R2_other_dir(:);
SIstd=std(R2_other_dir(:))/sqrt(numel(R2_other_dir));
R2_other_dir=nanmean(R2_other_dir(:));
end