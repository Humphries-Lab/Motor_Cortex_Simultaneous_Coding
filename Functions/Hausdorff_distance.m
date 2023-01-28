function [fraction_above,p,Hdist_vel,Hdist]=Hausdorff_distance(score_vel,idx_dir,idx_vel,Ndim,colourArea,nplot)
%% Hausdorff_distance calculates the Hausdorff distance between trjectories
%% in score_vel. This function assumes there are only two speed bins. 
%
% INPUTS
%
% score_vel: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score_vel
% matrix
% 
% idx_vel: array containing the duration bin of each row in the
% score_vel matrix
%
% Ndim: number of dimensions of the trajectories
%
% ColourArea: Colour of the dots to plot the distances.
%
% nplot: number of the subplot in the current figure where the results will
% be plotted
%
% OUTPUTS
%
% fraction_above: fraction of trajectories from same speed/distance that are
% closer than trajectories of adjacent direction bins
%
% p: p-value related to the comparison of distance between trajectories of adjacent direction bins and
% trajectories of same speed/distance. 
%
% Hdist_vel: Distance between trajectories of same direction and
% different speeds
%
% Hdist: Distance between trajectories of same speed and
% adjacent direction bins
%
% 27/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
Hdist_vel_dir=zeros(Ndir,Ndir);
Hdist_vel_dir2=zeros(Ndir,Ndir);
Hdist_vel=zeros(Ndir,1);
Ncomp=4;
Hdist=zeros(1,Ndir*Ncomp);
all_points=zeros(Ndir*Ncomp,2);

for i_dir=1:Ndir
    
    idx1=find(idx_dir==i_dir & idx_vel==1);
    idx2=find(idx_dir==i_dir & idx_vel==2);
    %% distance between trajectories
    recurrence=pdist2(score_vel(idx1,1:Ndim),score_vel(idx2,1:Ndim));
    Hdist_vel(i_dir)=max(min([recurrence,recurrence']));
    
    for j_dir=1:Ndir
        idx1_j=idx_dir==j_dir & idx_vel==1;
        idx2_j=idx_dir==j_dir & idx_vel==2;
        recurrence_1=pdist2(score_vel(idx1,1:Ndim),score_vel(idx1_j,1:Ndim));
        Hdist_vel_dir(i_dir,j_dir)=max(min([recurrence_1,recurrence_1']));
        
        recurrence_2=pdist2(score_vel(idx2,1:Ndim),score_vel(idx2_j,1:Ndim));
        Hdist_vel_dir2(i_dir,j_dir)=max(min([recurrence_2,recurrence_2']));
    end
    
    distance_slow=circshift(Hdist_vel_dir(i_dir,:),Ndir/2-i_dir);
    distance_fast=circshift(Hdist_vel_dir2(i_dir,:),Ndir/2-i_dir);
     
    all_points((i_dir-1)*Ncomp+1:i_dir*Ncomp,:)=[Hdist_vel(i_dir)*ones(4,1),[distance_fast(Ndir/2-1);distance_fast(Ndir/2+1);distance_slow(Ndir/2+1);distance_slow(Ndir/2-1)]];
    
    subplot(2,3,nplot)
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2-1),'.','Color',colourArea)
    hold on
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2-1),'.','Color',colourArea)
    
    Hdist((i_dir-1)*Ncomp+1:i_dir*Ncomp)=[distance_fast(Ndir/2-1),distance_fast(Ndir/2+1),distance_slow(Ndir/2+1),distance_slow(Ndir/2-1)];

end
fraction_above=sum(all_points(:,1)<all_points(:,2))/size(all_points,1);
[~,p]=ttest2(all_points(:,1),all_points(:,2));
box off
end