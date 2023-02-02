function plot_trajectories_prep_exec_after(score,idx_dir,t_from,duration_range)
%% plot_trajectories_prep_exec_after plots the trajectories corresponding to different directions (of the same duration range) during movement preparation,
%% movement execution and after movement end
%
% INPUTS 
%
% score: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score_vel
% matrix
%
% t_from: start time of the neural activity relative to movement onset [S]
%
% duration_range= range of the durations of the selected movements [S] 
%
% 02/02/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
mov_onset=round(abs(t_from)*1000);
mov_end=round((duration_range(2)+abs(t_from))*1000);


for i_dir=1:Ndir
        
        idx=find(idx_dir==i_dir);
       
        subplot(2,6,1)
        plot(score(idx(1:mov_onset),1),score(idx(1:mov_onset),2),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(mov_onset),1),score(idx(mov_onset),2),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
        
        subplot(2,6,2)
        plot(score(idx(mov_onset+1:mov_end),1),score(idx(mov_onset+1:mov_end),2),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(mov_end),1),score(idx(mov_end),2),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
      
        subplot(2,6,3)
        plot(score(idx(mov_end+1:end),1),score(idx(mov_end+1:end),2),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(end),1),score(idx(end),2),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
        
        subplot(2,6,7)
        plot(score(idx(1:mov_onset),3),score(idx(1:mov_onset),4),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(mov_onset),3),score(idx(mov_onset),4),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
      
        subplot(2,6,8)
        plot(score(idx(mov_onset+1:mov_end),3),score(idx(mov_onset+1:mov_end),4),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(mov_end),3),score(idx(mov_end),4),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
       
        subplot(2,6,9)
        plot(score(idx(mov_end+1:end),3),score(idx(mov_end+1:end),4),'Color',colour_dir(i_dir,:),'LineWidth',2)
        hold on
        plot(score(idx(end),3),score(idx(end),4),'^','Color',colour_dir(i_dir,:),'LineWidth',2)
       
        
end
subplot(2,6,1)
xlabel('PC 1')
ylabel('PC 2')
box off
xlim([min(score(:,1)) max(score(:,1))])
ylim([min(score(:,2)) max(score(:,2))])
title('Preparation')

subplot(2,6,2)
xlabel('PC 1')
ylabel('PC 2')
box off
xlim([min(score(:,1))-0.001 max(score(:,1))+0.001])
ylim([min(score(:,2))-0.001 max(score(:,2))+0.001])
title('Execution')

subplot(2,6,3)
xlabel('PC 1')
ylabel('PC 2')
box off
xlim([min(score(:,1))-0.001 max(score(:,1))+0.001])
ylim([min(score(:,2))-0.001 max(score(:,2))+0.001])
title('After Movement')

subplot(2,6,7)
xlabel('PC 3')
ylabel('PC 4')
box off
xlim([min(score(:,3))-0.001 max(score(:,3))+0.001])
ylim([min(score(:,4))-0.001 max(score(:,4))+0.001])

subplot(2,6,8)
xlabel('PC 3')
ylabel('PC 4')
box off
xlim([min(score(:,3))-0.001 max(score(:,3))+0.001])
ylim([min(score(:,4))-0.001 max(score(:,4))+0.001])

subplot(2,6,9)
xlabel('PC 3')
ylabel('PC 4')
box off
xlim([min(score(:,3))-0.001 max(score(:,3))+0.001])
ylim([min(score(:,4))-0.001 max(score(:,4))+0.001])


%% legend of colour directions
colour_dir=hsv(Ndir);
polaraxes('Position',[.8 .4 .15 .15],'Box','off');
thetacolor_axis=(-pi+0.0001:2*pi/Ndir:pi)+pi/Ndir;
coloraxis=ceil(Ndir*(thetacolor_axis+pi)/(2*pi));
hold on
for i=1:numel(thetacolor_axis)
    
    polarplot([thetacolor_axis(i) thetacolor_axis(i)],[0 1],'.-','Color',colour_dir(coloraxis(i),:),'LineWidth', 2)
    hold on
end
rticks([0 180])
thetaticks([0 90 180 270])
end