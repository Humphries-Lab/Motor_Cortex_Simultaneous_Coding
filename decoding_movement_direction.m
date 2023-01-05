function [mean_total_error,error_duration]=decoding_movement_direction(session,Area,ref_bin,Ndim)
%% approach to decode direction
%% 1) select half of the trials (we need examples of all durations and directions)
%% 2) train the decoder
%% 3) test in the other half
do_plot=1;
% Load data from normalised trajectories
%% Load PCA output
% Load data only from segments 300-400 ms length
%load(['scores_LDS_diff_duration_' session(1:5) '_' Area '.mat'],'score','idx_dir','idx_duration')
%load(['scores_LDS_diff_duration_ms_original_' session '_' Area '.mat'],'score','idx_dir','idx_duration')
load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance')

Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
% session='MM_S1_raw.mat';
% Area='M1';

% from=[0.2 0.3 0.4 0.5];
% to=from+0.1;
% t_1=-0.3;
% t_2=to+0.2;

test_bins=1:4;
test_bins(ref_bin)=[];
traj_ref=sum((idx_dir==1)& (idx_duration==ref_bin));
Ntest_bins=numel(test_bins);
for i_bin=1:Ntest_bins
    bin_i=test_bins(i_bin);
    traj_length(i_bin)=sum((idx_dir==1)& (idx_duration==bin_i));
end
fraction=zeros(Ndir,Ntest_bins);
estimated_duration=zeros(Ndir,Ntest_bins);



idx=find(idx_duration==ref_bin);
idx_dir_ref=idx_dir(idx);
%reference/ longest/slowest
x1=score(idx,1:Ndim);
clear idx

angle_bin=linspace(-pi,pi,Ndir+1);
angle_bin(1)=[];

 
 Mdl=train_bayes(x1,idx_dir_ref,do_plot)

for i_dir=1:Ndir
    
    
    
    for i_bin=1:Ntest_bins
        bin_i=test_bins(i_bin);
        idx=find((idx_dir==i_dir) & (idx_duration==bin_i));
        x2=score(idx,1:Ndim);
        %% shuffle
        %x2=x2(randperm(size(x2,1)),:);
        %fraction(i_dir,i_bin)=compute_fraction_speed_general_v3(x1,x2);
        sample=x2(250:300,:);
        angle=predict_from_bayes(sample,Mdl,angle_bin);
              
        %angle=estimate_direction(x1,x2,idx_dir_ref);

        if do_plot
            subplot(2,2,1)
            %from bayes
            plot(angle,angle_bin(i_dir),'.','Color',colour_dir(i_dir,:))
            %[angle,angle_bin(i_dir)]
            %plot(angle_bin(angle),angle_bin(i_dir),'.','Color',colour_dir(i_dir,:))
            hold on
            subplot(2,2,2)
            plot3(sample(:,1),sample(:,2),sample(:,3),'.','Color',colour_dir(i_dir,:))
            pause
        end
        %from Bayes
        error_angle(i_dir,i_bin)=abs(angdiff(mean(angle),angle_bin(i_dir))*180/pi)
        %error_angle(i_dir,i_bin)=abs(angdiff(angle_bin(angle),angle_bin(i_dir))*180/pi);
        clear idx
    end
    
    if do_plot
        
        subplot(2,2,3)
        plot(error_angle(i_dir,:),'.','Color',colour_dir(i_dir,:))
        hold on
    end
end
mean_error=mean(error_angle);
mean_total_error=mean(error_angle(:));

if do_plot
    subplot(2,2,1)
    plot(1:10,1:10,'.k-','Linewidth',1.5)
    xlabel('Prediction')
    ylabel('Angle')
    box off
    %
    % subplot(2,1,2)
    %
    % plot(traj_length,mean_error,'.k-','Linewidth',1.5)
    % xlabel('Trajectory length')
    % ylabel('Decoding error')
    % box off
end


end

function  angle=estimate_direction(x1,x2,idx_dir_ref)
Ndir=max(idx_dir_ref);
colour_dir=hsv(Ndir);
subplot(2,1,1)
for i_dir=1:Ndir
    idx=find(idx_dir_ref==i_dir);
    %      plot3(x1(idx,1),x1(idx,2),x1(idx,3),'Color',colour_dir(i_dir,:))
    %      hold on
    dist_matrix=pdist2(x1(idx,:),x2);
    mean_dist(i_dir)=mean(dist_matrix(:));
end
plot3(x2(:,1),x2(:,2),x2(:,3),'k')
%
% subplot(2,1,2)
% plot(mean_dist)
% pause
[~,angle]=min(mean_dist);
end

function Mdl=train_bayes(score,idx_dir,debugging)
Ndir=max(idx_dir);
from=250;
to=300;
if debugging
    colour_dir=hsv(Ndir);
end
score2=[];
idx_dir2=[];
for i=1:Ndir
    idx=find(idx_dir==i);
    if debugging
        subplot(2,2,2)
        plot3(score(idx(from:to),1),score(idx(from:to),2),score(idx(from:to),3),'.','Color',colour_dir(i,:))
        hold on
    end
    score2=[score2;score(idx(from:to),:)];
    idx_dir2=[idx_dir2;idx_dir(idx(from:to))];
end

Mdl = fitcnb(score2,idx_dir2);
end

function angle=predict_from_bayes(sample,Mdl,xangle)
debugging=0;
[label,Posterior] = predict(Mdl,sample)
if debugging
figure
polarplot(xangle,ones(numel(xangle),1))
hold on 
Posterior.*xangle

sum(Posterior.*xangle)
polarplot(xangle,Posterior)
polarplot([0 sum(Posterior.*xangle)],[0 1.5])
%polarplot([0 sum(Posterior.*xangle)],[0])
end
angle=sum(Posterior.*xangle,2)
end