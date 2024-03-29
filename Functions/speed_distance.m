function [fraction_above_speed,p_speed,fraction_above_dist,p_dist,average_speed_bin,average_dist_bin,Distance_speed,Distance_dir_speed,Distance_dist,Distance_dir_dist,R2,power_sample]=speed_distance(Session, Area,threshold,Ndir,t_from,t_upto)
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
%
% fraction_above_speed(dist): fraction of trajectories from the same speed (distance) that are
% closer than trajectories of adjacent direction bins
%
% p_speed(dist): related to the comparison of the distance between trajectories of adjacent direction bins and
% trajectories of same speed(distance). 
%
% average_speed(dist)_bin: average speed(distance) for slow and fast (short
% and long) movements 
%
% Distance_speed(dist): Distance between trajectories of the same direction and
% different speeds (distance)
%
% Distance_dir_speed(dist): Distance between trajectories of the same speed (distance) and
% adjacent direction bins
%
%
% 27/05/2023
% Andrea Colins Rodriguez

% aim: normalise time of individual trajectories and then do PCA

%consider reaches of duration between 0.2 and 1.2 s
from_to=[0.2 1.2]; 
% final number of points for the normalised trajectories
final_length=600;
% 
% averages_vel=[];
% idx_dir_2=[];
% idx_vel_2=[];
threshold_speed=20;
threshold_dist=5;
if strcmp(Area,'M1')
    colourArea=[85 30 116]./256;
else
    colourArea=[89 156 153]./256;
end

load(['../Output_files/PCA_' Session(1:end-4) '_' Area '.mat'],'coeffs','variance','delete_units','normalisation','sigma_filter','nsamples_condition')


[Neural_info,Mov_params]=neural_data_per_duration_normalised(Session,Area,sigma_filter,t_from,t_upto,from_to);

%discretise dir, vel and distance
%better do this by percentil
direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));

limit_vel=median(Mov_params.max_speed);
speed_discrete=ones(size(Mov_params.max_speed));
speed_discrete(Mov_params.max_speed>limit_vel)=2;


limit_dist=median(Mov_params.distance(abs(Mov_params.max_speed-threshold_speed)<2))+0.001;
dist_discrete=ones(size(Mov_params.max_speed));
dist_discrete(Mov_params.distance>=limit_dist)=2;
average_speed_bin=[mean(Mov_params.max_speed(speed_discrete==1 & abs(Mov_params.distance-threshold_dist)<1)) mean(Mov_params.max_speed(speed_discrete==2 & abs(Mov_params.distance-threshold_dist)<1))];
average_dist_bin=[mean(Mov_params.distance(dist_discrete==1 & abs(Mov_params.max_speed-threshold_speed)<2)) mean(Mov_params.distance(dist_discrete==2 & abs(Mov_params.max_speed-threshold_speed)<2))];

averages_vel=nan(size(Neural_info.FR,1),Ndir*2*final_length);
averages_dist=nan(size(Neural_info.FR,1),Ndir*2*final_length);
idx_vel_2=nan(Ndir*2*final_length,1);
idx_dir_2=nan(Ndir*2*final_length,1);
counter=1;
Nsamples=nan(Ndir*2,1);
for i_dir=1:Ndir
     for i_vel=1:2
        
        % to select only the ones with the same distance
        constant_dist=(speed_discrete==i_vel & direction1==i_dir & abs(Mov_params.distance-threshold_dist)<1);
        Nsamples(counter)=sum(constant_dist);
        averages_vel=[averages_vel,mean(Neural_info.FR(:,:,constant_dist),3)];
        idx_vel_2=[idx_vel_2;zeros(final_length,1)+i_vel];
        
        averages_dist=[averages_dist,mean(Neural_info.FR(:,:,dist_discrete==i_vel & direction1==i_dir & abs(Mov_params.max_speed-threshold_speed)<2),3)];
        idx_dir_2=[idx_dir_2;zeros(final_length,1)+i_dir];
         
        counter=counter+1;
        
      
 
    end
end

power_sample=[mean(Nsamples) mean(nsamples_condition,'all') mean(Nsamples) mean(nsamples_condition,'all') median(Nsamples) median(nsamples_condition,'all')];
power_sample=[100*power_sample(1:2)./power_sample(1) power_sample(3:end)];
%% Project same directions coloring different times
%soft normalization
% averages_vel(delete_units,:)=[];
% averages_dist(delete_units,:)=[];
% averages_vel=averages_vel'./repmat(range(averages_vel')+threshold_dist,size(averages_vel,2),1);
% averages_dist=averages_dist'./repmat(range(averages_dist')+threshold_dist,size(averages_dist,2),1);
% 
% [~,score_vel,~,~,exp_1]=pca(averages_vel);
% [~,score_dist,~,~,exp_dist_1]=pca(averages_dist);
%colour_dir=hsv(Ndir);


% % soft normalization
averages_vel(delete_units,:)=[];
averages_dist(delete_units,:)=[];
averages_vel=averages_vel'./repmat(normalisation,size(averages_vel,2),1);
averages_dist=averages_dist'./repmat(normalisation,size(averages_dist,2),1);
 Ndim=find(cumsum(variance)>threshold,1,'First');
%projecting trajectories onto the subspace defined from all durations and
%all dir
score_vel=averages_vel*coeffs(:,1:Ndim);
score_dist=averages_dist*coeffs(:,1:Ndim);

% figure
% hold on
%     for i_dir=1:Ndir
%         for i_speed=1:2
%             idx=idx_dir_2==i_dir & idx_vel_2==i_speed;
%             plot3(score_dist(idx,1),score_dist(idx,2),score_dist(idx,3),'Color',colour_dir(i_dir,:)/i_speed)    
%         end
%     end

%% Hausdorff distance for different speeds

[fraction_above_speed,p_speed,Distance_speed,Distance_dir_speed]=Hausdorff_distance(score_vel,idx_dir_2,idx_vel_2,Ndim,colourArea,3);

%% now for distance
%[fraction_above_dist,p_dist,Distance_dist,Distance_dir_dist]=Hausdorff_distance(score_dist,idx_dir_2,idx_vel_2,Ndim,colourArea,6);

%% R^2
fraction_above_dist=fraction_above_speed;
p_dist=p_speed;
Distance_dist=Distance_speed;
Distance_dir_dist=Distance_dir_speed;
%% R^2
final_length=600;
normalised_t=linspace(0,1,final_length);
[~,~,R2.vel_same_dir]=R_same(score_vel,idx_dir_2,idx_vel_2,Ndim,normalised_t);
[~,~,R2.vel_other_dir]=R2_other(score_vel,idx_dir_2,idx_vel_2,Ndim,normalised_t);

[~,~,R2.dist_same_dir]=R_same(score_dist,idx_dir_2,idx_vel_2,Ndim,normalised_t);
[~,~,R2.dist_other_dir]=R2_other(score_dist,idx_dir_2,idx_vel_2,Ndim,normalised_t);
end

