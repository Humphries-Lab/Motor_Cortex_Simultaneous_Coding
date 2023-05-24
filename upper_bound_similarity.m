function upper_bound_similarity(Sessions,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin)
%% upper_bound_similarity calculates the upper bound of the similarity between trajectories for all recordings denoted in Sessions and Area
%
% INPUTS
%
% Sessions: cell array containing the names of the sessions to be analysed.
% e.g {'MC_S1_raw.mat','MC_S2_raw.mat'}
%
% Areas: cell array containing the names of the areas to be analysed in
% each session. Areas and Sessions must have the same number of elements.
% e.g {'M1','M1'}
% 
% threshold: percentage of the variance to be explained by the first nPCs
%
% Ndir: number of direction to bin the movements
%
% Nbins: number of durations to bin the movements
%
% t_from: start time of the neural activity relative to movement onset for all recordings[S]
% e.g t_from=[-0.5 -0.5]
%
% t_upto: end time of the neural activity relative to movement end for all recordings [S]
% e.g t_from=[0.3 0.3]
%
% edges_dur_bin= array containing the edges of each duration bin [S]
%
%
% 24/01/2023
% Andrea Colins Rodriguez

for isession=1:size(Sessions,2)
    r=embedding_dimensions_upper_bound(Sessions{isession}, Area{isession},Ndir,Nbins,threshold,t_from(isession),t_upto(isession),edges_dur_bin);
end
end