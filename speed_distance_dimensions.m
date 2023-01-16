function [fraction_above_vel,p_vel,fraction_above_dist,p_dist,average_speed_bin,average_dist_bin,Distance_vel,Distance_v_dir,Distance_dist,Distance_d_dir]=speed_distance_dimensions(session, Area,threshold,Ndir,t_1,t_2,do_plot)
% aim: normalise time of individual trajectories and then do PCA

%% Average across bins of segments durations (200,300,400,500 ms)
colour_dir=hsv(Ndir);
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
    colourArea='m';
else
    colourArea='b';
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
ms=nanmedian(ISI);

[condition_matrix,direction,~,~,~,reach_number,dist_mov_dir,mov_duration,max_speed,prep_duration]=neural_data_per_duration_normalised(session,Area,ms,t_1,t_2,event,from_to);


%discretise dir, vel and distance
%better do this by percentil
direction1=ceil(Ndir*(direction+pi)/(2*pi));

limit_vel=median(max_speed(abs(dist_mov_dir-5)<1));
speed_discrete=ones(size(max_speed));
speed_discrete(max_speed>limit_vel)=2;


limit_dist=median(dist_mov_dir(abs(max_speed-threshold_speed)<2))+0.001;
dist_discrete=ones(size(max_speed));
dist_discrete(dist_mov_dir>=limit_dist)=2;
average_speed_bin=[mean(max_speed(speed_discrete==1 & abs(dist_mov_dir-5)<1)) mean(max_speed(speed_discrete==2 & abs(dist_mov_dir-5)<1))];
average_dist_bin=[mean(dist_mov_dir(dist_discrete==1 & abs(max_speed-threshold_speed)<2)) mean(dist_mov_dir(dist_discrete==2 & abs(max_speed-threshold_speed)<2))];
if do_plot
    subplot(2,5,2)
    histogram(max_speed(speed_discrete==1 & abs(dist_mov_dir-5)<1),0:1:max(max_speed),'Normalization','Probability')
    hold on
    histogram(max_speed(speed_discrete==2 & abs(dist_mov_dir-5)<1),0:1:max(max_speed),'Normalization','Probability')
    box off
    xlabel('Speed [cm/s]')
    ylabel('Normalised Frequency')
    subplot(2,5,5+2)
    histogram(dist_mov_dir(dist_discrete==1 & abs(max_speed-threshold_speed)<2),0:0.5:max(dist_mov_dir),'Normalization','Probability')
    hold on
    histogram(dist_mov_dir(dist_discrete==2 & abs(max_speed-threshold_speed)<2),0:0.5:max(dist_mov_dir),'Normalization','Probability')
    box off
    xlabel('Distance [cm]')
    ylabel('Normalised Frequency')
end

for cond=1:Ndir
    average_cond_1=[average_cond_1,mean(condition_matrix(:,:,direction1==cond),3)];
    idx_dir=[idx_dir;zeros(final_length,1)+cond];
    for i_vel=1:2
        
        %averages_vel=[averages_vel,mean(condition_matrix(:,:,speed_discrete==i_vel & direction1==cond),3)];
        % to select only the ones with the same distance
        averages_vel=[averages_vel,mean(condition_matrix(:,:,speed_discrete==i_vel & direction1==cond & abs(dist_mov_dir-5)<1),3)];
        idx_dir_2=[idx_dir_2;zeros(final_length,1)+cond];
        idx_vel_2=[idx_vel_2;zeros(final_length,1)+i_vel];
        
        averages_dist=[averages_dist,mean(condition_matrix(:,:,dist_discrete==i_vel & direction1==cond & abs(max_speed-threshold_speed)<2),3)];
        nsamples_vel(cond,i_vel)=[sum(speed_discrete==i_vel & direction1==cond & abs(dist_mov_dir-5)<1)];
        
    end
end

%% Compute a common subspace

%% Project same directions coloring different times
%[coeff,score,~,~,exp_all]=pca(average_cond_1);
% [~,score_vel,~,~,exp_1]=pca(averages_vel);
% [~,score_dist,~,~,exp_dist_1]=pca(averages_dist);
load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'coeffs','score','idx_dir','idx_duration','variance','t_1','t_2','from','delete_units','normalisation')
exp_1=variance;
exp_dist_1=variance;
% soft normalization
averages_vel(delete_units,:)=[];
averages_dist(delete_units,:)=[];
%average_cond_1=average_cond_1'./repmat(range(average_cond_1')+5,size(average_cond_1,2),1);
averages_vel=averages_vel'./repmat(normalisation,size(averages_vel,2),1);
averages_dist=averages_dist'./repmat(normalisation,size(averages_dist,2),1);
Ndim=find(cumsum(exp_1)>threshold,1,'First');
score_vel=averages_vel*coeffs(:,1:Ndim);
score_dist=averages_dist*coeffs(:,1:Ndim);
optional_plot=0;
if optional_plot
    
    for i_dir=1:Ndir
        [~,score2]=pca(score_vel(idx_dir_2==i_dir,1:Ndim));
        for i_vel=1:2
            idx=idx_vel_2==i_vel & idx_dir_2==i_dir;
            idx2=idx_vel_2(idx_dir_2==i_dir)==i_vel;
            subplot(2,5,4)
            plot(score2(idx2,1)+1,score2(idx2,2)+1,'Color',colour_dir(i_dir,:)./i_vel)
           
            Area_curve_speed(i_dir,i_vel)=abs(area_closed_curve(score2(idx2,1)+1,score2(idx2,2)+1));
            %plot3(score_vel(idx,1),score_vel(idx,2),score_vel(idx,3),'Color',colour_dir(i_dir,:)./i_vel,'LineWidth',i_vel)
            hold on
        end
        
    end

    subplot(2,5,5)
        plot(1:2,mean(Area_curve_speed)./mean(Area_curve_speed(:,1)),'.-','Color',colourArea)
        hold on
    title('speed')
    
    
    for i_dir=1:Ndir
        [~,score2]=pca(score_dist(idx_dir_2==i_dir,1:Ndim));
        for i_vel=1:2
            idx=idx_vel_2==i_vel & idx_dir_2==i_dir;
            idx2=idx_vel_2(idx_dir_2==i_dir)==i_vel;
            subplot(2,5,5+4)
            %plot3(score_dist(idx,1),score_dist(idx,2),score_dist(idx,3),'Color',colour_dir(i_dir,:)./i_vel,'LineWidth',i_vel)
            plot(score2(idx2,1)+1,score2(idx2,2)+1,'Color',colour_dir(i_dir,:)./i_vel)
            hold on
            Area_curve_dist(i_dir,i_vel)=abs(area_closed_curve(score2(idx2,1)+1,score2(idx2,2)+1));
        end
        title('Distance')
        
    end
    subplot(2,5,5+5)
        plot(1:2,mean(Area_curve_dist)./mean(Area_curve_dist(:,1)),'.-','Color',colourArea)
        hold on
    
end

%% Hausdorff distance for different speeds
Ndim=find(cumsum(exp_1)>threshold,1,'First');
[fraction_above_vel,p_vel,Distance_vel,Distance_v_dir]=Hausdorff_distance(score_vel,idx_dir_2,idx_vel_2,Ndim,colourArea,1);



%%
%% now for distance
Ndim=find(cumsum(exp_dist_1)>threshold,1,'First');
[fraction_above_dist,p_dist,Distance_dist,Distance_d_dir]=Hausdorff_distance(score_dist,idx_dir_2,idx_vel_2,Ndim,colourArea,6);
p_dist


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
    
    subplot(2,5,nplot)
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2-1),'.','Color',colourArea)
    hold on
    plot(Hdist_vel(i_dir),distance_fast(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2+1),'.','Color',colourArea)
    plot(Hdist_vel(i_dir),distance_slow(Ndir/2-1),'.','Color',colourArea)
    
    Hdist=[Hdist,distance_fast(Ndir/2-1),distance_fast(Ndir/2+1),distance_slow(Ndir/2+1),distance_slow(Ndir/2-1)];
end
fraction_above=sum(all_points(:,1)<all_points(:,2))/size(all_points,1);
[h,p]=ttest2(all_points(:,1),all_points(:,2));
box off
end

function A=area_closed_curve(x,y)
N=numel(x);
I=zeros(N,1);
for i=1:N-1
    I(i)=(y(i+1)+y(i))*(x(i+1)-x(i))/2;
end
I(N)=(y(1)+y(N))*(x(1)-x(N))/2;
A=sum(I);
end