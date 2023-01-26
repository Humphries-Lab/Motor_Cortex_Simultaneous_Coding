function [fraction_above_vel,p_vel,fraction_above_dist,p_dist,average_speed_bin,average_dist_bin,Distance_vel,Distance_v_dir,Distance_dist,Distance_d_dir]=speed_distance_dimensions(session, Area,threshold,Ndir,t_1,t_2)
% aim: normalise time of individual trajectories and then do PCA
event=2;
%consider reaches from lenght "from" to "to"
from=0.2;
to=from+1;

% final number of points for the normalised trajectories
final_length=600;

average_cond_1=[];
idx_dir=[];
averages_vel=[];
averages_dist=[];
idx_dir_2=[];
idx_vel_2=[];
threshold_speed=20;
nsamples_vel=zeros(Ndir,2);


if strcmp(Area,'M1')
    colourArea=[85 30 116]./256;
else
    colourArea=[89 156 153]./256;
end

from_to=[from to];
if strcmp(Area,'PMd')
    load(session,'PMd','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(PMd.units,startt,endt);
else
    load(session,'M1','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(M1.units,startt,endt);
end
sigma_filter=round(median(ISI));
[Neural_info,Mov_params]=neural_data_per_duration_normalised(session,Area,sigma_filter,t_1,t_2,event,from_to);

%discretise dir, vel and distance
%better do this by percentil
direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));

limit_vel=median(Mov_params.max_speed(abs(Mov_params.distance-5)<1));
speed_discrete=ones(size(Mov_params.max_speed));
speed_discrete(Mov_params.max_speed>limit_vel)=2;


limit_dist=median(Mov_params.distance(abs(Mov_params.max_speed-threshold_speed)<2))+0.001;
dist_discrete=ones(size(Mov_params.max_speed));
dist_discrete(Mov_params.distance>=limit_dist)=2;
average_speed_bin=[mean(Mov_params.max_speed(speed_discrete==1 & abs(Mov_params.distance-5)<1)) mean(Mov_params.max_speed(speed_discrete==2 & abs(Mov_params.distance-5)<1))];
average_dist_bin=[mean(Mov_params.distance(dist_discrete==1 & abs(Mov_params.max_speed-threshold_speed)<2)) mean(Mov_params.distance(dist_discrete==2 & abs(Mov_params.max_speed-threshold_speed)<2))];

%% Histograms of the speed and distance of the selected movements
%if do_plot
%     subplot(2,3,2)
%     histogram(max_speed(speed_discrete==1 & abs(dist_mov_dir-5)<1),0:1:max(max_speed),'Normalization','Probability')
%     hold on
%     histogram(max_speed(speed_discrete==2 & abs(dist_mov_dir-5)<1),0:1:max(max_speed),'Normalization','Probability')
%     box off
%     xlabel('Speed [cm/s]')
%     ylabel('Normalised Frequency')
%     subplot(2,5,5+2)
%     histogram(dist_mov_dir(dist_discrete==1 & abs(max_speed-threshold_speed)<2),0:0.5:max(dist_mov_dir),'Normalization','Probability')
%     hold on
%     histogram(dist_mov_dir(dist_discrete==2 & abs(max_speed-threshold_speed)<2),0:0.5:max(dist_mov_dir),'Normalization','Probability')
%     box off
%     xlabel('Distance [cm]')
%     ylabel('Normalised Frequency')
% end

for cond=1:Ndir
    average_cond_1=[average_cond_1,mean(Neural_info.FR(:,:,direction1==cond),3)];
    idx_dir=[idx_dir;zeros(final_length,1)+cond];
    for i_vel=1:2
        
        %averages_vel=[averages_vel,mean(condition_matrix(:,:,speed_discrete==i_vel & direction1==cond),3)];
        % to select only the ones with the same distance
        averages_vel=[averages_vel,mean(Neural_info.FR(:,:,speed_discrete==i_vel & direction1==cond & abs(Mov_params.distance-5)<1),3)];
        idx_dir_2=[idx_dir_2;zeros(final_length,1)+cond];
        idx_vel_2=[idx_vel_2;zeros(final_length,1)+i_vel];
        
        averages_dist=[averages_dist,mean(Neural_info.FR(:,:,dist_discrete==i_vel & direction1==cond & abs(Mov_params.max_speed-threshold_speed)<2),3)];
        nsamples_vel(cond,i_vel)=[sum(speed_discrete==i_vel & direction1==cond & abs(Mov_params.distance-5)<1)];
        
    end
end


%% Project same directions coloring different times
%[coeff,score,~,~,exp_all]=pca(average_cond_1);
% [~,score_vel,~,~,exp_1]=pca(averages_vel);
% [~,score_dist,~,~,exp_dist_1]=pca(averages_dist);
load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'coeffs','variance','delete_units','normalisation')
   
%load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'coeffs','score','idx_dir','idx_duration','variance','t_1','t_2','from','delete_units','normalisation')
exp_1=variance;
exp_dist_1=variance;
% soft normalization
averages_vel(delete_units,:)=[];
averages_dist(delete_units,:)=[];
%average_cond_1=average_cond_1'./repmat(range(average_cond_1')+5,size(average_cond_1,2),1);
averages_vel=averages_vel'./repmat(normalisation,size(averages_vel,2),1);
averages_dist=averages_dist'./repmat(normalisation,size(averages_dist,2),1);


Ndim=find(cumsum(exp_1)>threshold,1,'First');

%projecting trajectories onto the subspace defined from all durations and
%all dir
score_vel=averages_vel*coeffs(:,1:Ndim);
score_dist=averages_dist*coeffs(:,1:Ndim);

%% Hausdorff distance for different speeds
Ndim=find(cumsum(exp_1)>threshold,1,'First');
[fraction_above_vel,p_vel,Distance_vel,Distance_v_dir]=Hausdorff_distance(score_vel,idx_dir_2,idx_vel_2,Ndim,colourArea,5);


%%
%% now for distance
Ndim=find(cumsum(exp_dist_1)>threshold,1,'First');
[fraction_above_dist,p_dist,Distance_dist,Distance_d_dir]=Hausdorff_distance(score_dist,idx_dir_2,idx_vel_2,Ndim,colourArea,6);


end

function [fraction_above,p,Hdist_vel,Hdist]=Hausdorff_distance(score_vel,idx_dir_2,idx_vel_2,Ndim,colourArea,nplot)
Ndir=max(idx_dir_2);
Hdist_vel_dir=zeros(Ndir,Ndir);
Hdist_vel_dir2=zeros(Ndir,Ndir);
Hdist_vel=zeros(Ndir,1);
Hdist=[];
all_points=[];
for i_dir=1:Ndir
    
    idx1=find(idx_dir_2==i_dir & idx_vel_2==1);
    idx2=find(idx_dir_2==i_dir & idx_vel_2==2);
    %% distance between trajectories
    recurrence=pdist2(score_vel(idx1,1:Ndim),score_vel(idx2,1:Ndim));
    Hdist_vel(i_dir)=max([min(recurrence),min(recurrence')]);
    
    for j_dir=1:Ndir
        idx1_j=idx_dir_2==j_dir & idx_vel_2==1;
        idx2_j=idx_dir_2==j_dir & idx_vel_2==2;
        recurrence_1=pdist2(score_vel(idx1,1:Ndim),score_vel(idx1_j,1:Ndim));
        Hdist_vel_dir(i_dir,j_dir)=max([min(recurrence_1),min(recurrence_1')]);
        
        recurrence_2=pdist2(score_vel(idx2,1:Ndim),score_vel(idx2_j,1:Ndim));
        Hdist_vel_dir2(i_dir,j_dir)=max([min(recurrence_2),min(recurrence_2')]);
    end
    
    distance_slow=circshift(Hdist_vel_dir(i_dir,:),Ndir/2-i_dir);
    distance_fast=circshift(Hdist_vel_dir2(i_dir,:),Ndir/2-i_dir);
    all_points=[all_points;Hdist_vel(i_dir)*ones(4,1),[distance_fast(Ndir/2-1);distance_fast(Ndir/2+1);distance_slow(Ndir/2+1);distance_slow(Ndir/2-1)]];
    
    subplot(2,3,nplot)
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2-1),'.','Color',colourArea)
    hold on
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2-1),'.','Color',colourArea)
    
    Hdist=[Hdist,distance_fast(Ndir/2-1),distance_fast(Ndir/2+1),distance_slow(Ndir/2+1),distance_slow(Ndir/2-1)];
end
fraction_above=sum(all_points(:,1)<all_points(:,2))/size(all_points,1);
[~,p]=ttest2(all_points(:,1),all_points(:,2));
box off
end