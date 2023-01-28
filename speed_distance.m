function [fraction_above_speed,p_speed,fraction_above_dist,p_dist,average_speed_bin,average_dist_bin,Distance_speed,Distance_dir_speed,Distance_dist,Distance_dir_dist]=speed_distance(Session, Area,threshold,Ndir,t_from,t_upto)
%% speed_distance compares the distance between trajectories corresponding
%% to different speeds (& distances) vs adjacent direction bins
%
% INPUTS
%
% Session: Name of the session to be analysed.
% e.g 'MC_S1_raw.mat'
%
% Area: Name of the area to be analysed in the Session.
% e.g 'M1'
%
% threshold: percentage of the variance to be explained by the first nPCs
%
% Ndir: number of directions to bin the movements
%
% t_from: start time of the neural activity relative to movement onset [S]
% e.g t_from=-0.5
%
% t_upto: end time of the neural activity relative to movement end [S]
% e.g t_from=0.3
%
%
% OUTPUTS
%
% fraction_above_speed(dist): fraction of trajectories from same speed (distance) that are
% closer than trajectories of adjacent direction bins
%
% p_speed(dist): related to the comparison of distance between trajectories of adjacent direction bins and
% trajectories of same speed(distance). 
%
% average_speed(dist)_bin: average speed(distance) for slow and fast (short
% and long) movements 
%
% Distance_speed(dist): Distance between trajectories of same direction and
% different speeds (distance)
%
% Distance_dir_speed(dist): Distance between trajectories of same speed (distance) and
% adjacent direction bins
%
% 27/05/2023
% Andrea Colins Rodriguez

% aim: normalise time of individual trajectories and then do PCA
event=2;
%consider reaches of duration between 0.2 and 1.2 s
from_to=[0.2 1.2]; 
% final number of points for the normalised trajectories
final_length=600;

averages_vel=[];
averages_dist=[];
idx_dir_2=[];
idx_vel_2=[];
threshold_speed=20;

if strcmp(Area,'M1')
    colourArea=[85 30 116]./256;
else
    colourArea=[89 156 153]./256;
end

load(['../Output_files/PCA_' Session(1:end-4) '_' Area '.mat'],'coeffs','variance','delete_units','normalisation','sigma_filter')


[Neural_info,Mov_params]=neural_data_per_duration_normalised(Session,Area,sigma_filter,t_from,t_upto,event,from_to);

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

for cond=1:Ndir
     for i_vel=1:2
        
        % to select only the ones with the same distance
        averages_vel=[averages_vel,mean(Neural_info.FR(:,:,speed_discrete==i_vel & direction1==cond & abs(Mov_params.distance-5)<1),3)];
        idx_vel_2=[idx_vel_2;zeros(final_length,1)+i_vel];
        
        averages_dist=[averages_dist,mean(Neural_info.FR(:,:,dist_discrete==i_vel & direction1==cond & abs(Mov_params.max_speed-threshold_speed)<2),3)];
        idx_dir_2=[idx_dir_2;zeros(final_length,1)+cond];
         
    end
end


%% Project same directions coloring different times
%soft normalization
% averages_vel(delete_units,:)=[];
% averages_dist(delete_units,:)=[];
% averages_vel=averages_vel'./repmat(range(averages_vel')+5,size(averages_vel,2),1);
% averages_dist=averages_dist'./repmat(range(averages_dist')+5,size(averages_dist,2),1);
% 
% [~,score_vel,~,~,exp_1]=pca(averages_vel);
% [~,score_dist,~,~,exp_dist_1]=pca(averages_dist);


% soft normalization
averages_vel(delete_units,:)=[];
averages_dist(delete_units,:)=[];
averages_vel=averages_vel'./repmat(normalisation,size(averages_vel,2),1);
averages_dist=averages_dist'./repmat(normalisation,size(averages_dist,2),1);
Ndim=find(cumsum(variance)>threshold,1,'First');
%projecting trajectories onto the subspace defined from all durations and
%all dir
score_vel=averages_vel*coeffs(:,1:Ndim);
score_dist=averages_dist*coeffs(:,1:Ndim);

%% Hausdorff distance for different speeds

[fraction_above_speed,p_speed,Distance_speed,Distance_dir_speed]=Hausdorff_distance(score_vel,idx_dir_2,idx_vel_2,Ndim,colourArea,5);

%% now for distance
[fraction_above_dist,p_dist,Distance_dist,Distance_dir_dist]=Hausdorff_distance(score_dist,idx_dir_2,idx_vel_2,Ndim,colourArea,6);

end

