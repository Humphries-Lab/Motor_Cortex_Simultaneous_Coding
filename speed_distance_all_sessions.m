function speed_distance_all_sessions(session,Area,threshold,Ndir,t_from,t_upto)
%% speed_distance_all_sessions compares the distance between
%% trajectories corresponding to movements with different speeds (same
%% direction and distance) and movemements of adjacent bin directions (same speed)
%% It then repeats the process for movement distance (same speed) vs directions
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
% Ndir: number of directions to bin the movements
%
% t_from: start time of the neural activity relative to movement onset for all recordings[S]
% e.g t_from=[-0.5 -0.5]
%
% t_upto: end time of the neural activity relative to movement end for all recordings [S]
% e.g t_from=[0.3 0.3]
%
% 27/05/2023
% Andrea Colins Rodriguez

Nsessions=size(session,2);
fraction_above_vel=zeros(Nsessions,1);
p_vel=zeros(Nsessions,1);
fraction_above_dist=zeros(Nsessions,1);
Dur_var_exp=zeros(Nsessions,3);
p_dist=zeros(Nsessions,1);
average_speed_bin=zeros(Nsessions,2);
average_dist_bin=average_speed_bin;
Distance_vel=nan(Nsessions,Ndir);
Distance_dist=nan(Nsessions,Ndir);
Distance_v_dir=nan(Nsessions,32);
Distance_d_dir=nan(Nsessions,32);
power_sample=nan(Nsessions,6);
do_plot_r2=0;


for isession=1:Nsessions
    
    [fraction_above_vel(isession),p_vel(isession),fraction_above_dist(isession),p_dist(isession),average_speed_bin(isession,:),average_dist_bin(isession,:),Distance_vel(isession,:),Distance_v_dir(isession,:),Distance_dist(isession,:),Distance_d_dir(isession,:),R2,power_sample(isession,:)]=speed_distance(session{isession}, Area{isession},threshold,Ndir,t_from(isession),t_upto(isession));
    
    if do_plot_r2
        % compute coefficient of determination between trajectories of diff
        % speed or distances
        plot_R2(R2,Area{isession},isession)
        [~,p_vel_R(isession)]=ttest2(R2.vel_same_dir,R2.vel_other_dir);
        [~,p_dist_R(isession)]=ttest2(R2.dist_same_dir,R2.dist_other_dir);
    end
    
end

disp(['Mean power sample speed = ' num2str(mean(power_sample))])
if do_plot_r2
    subplot(2,3,1)
    text(1, 0, [' max p-val = ' num2str(max(p_vel_R),2)])
    
    subplot(2,3,2)
    text(1, 0, [' max p-val = ' num2str(max(p_dist_R),2)])
end

subplot(2,3,3)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance speed')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])
text(0.005,0.02,[' max p-value = ' num2str(max(p_vel))],'FontSize',8)
tmp=reshape(repmat(Distance_vel,4,1),13,32);
diag_hist(tmp(:)-Distance_v_dir(:),0.025)
axis square

% subplot(2,3,6)
% plot([0 0.025],[0 0.025],'k')
% xlabel('Distance distance')
% ylabel('Distance closest directions')
% xlim([0 0.025])
% ylim([0 0.025])
% text(0.005,0.02,[' max p-value = ' num2str(max(p_dist))],'FontSize',8)
% disp('--------------------------------')
% disp('Recordings which p-values > 0.05 (Distance vs direction, FIG 5F, right)')
% Session_name=session(p_dist>0.05)';
% Area_name=Area(p_dist>0.05)';
% pvalue_dist=p_dist(p_dist>0.05);
% pvalues=table(Session_name,Area_name,pvalue_dist)
% 
% disp('--------------------------------')


% axes('Position',[0.4116    0.47    0.21   0.041]);
% histogram(Distance_vel(:),'Normalization','probability')
% xlim([0 0.025])
% xlabel('Distance Speed')
% box off
% 
% axes('Position',[ 0.6331    0.1116    0.0164    0.34]);
% histogram(Distance_v_dir(:),'Normalization','probability','orientation','horizontal')
% ylim([0 0.025])
% ylabel('Distance Closest direction')
% box off

% axes('Position',[ 0.6936    0.47    0.2100    0.0410]);
% histogram(Distance_dist(:),'Normalization','probability')
% xlim([0 0.025])
% box off

% axes('Position',[ 0.9152    0.1187    0.0164    0.3380]);
% histogram(Distance_d_dir(:),'Normalization','probability','orientation','horizontal')
% ylim([0 0.025])
% box off

disp('--------------------------------')
disp(['Mean lower speed = ' num2str(mean(average_speed_bin(:,1))) ' [cm/s]'])
disp(['Mean higher speed = ' num2str(mean(average_speed_bin(:,2))) ' [cm/s]'])
disp('--------------------------------')
%disp(['Mean shorter distance = ' num2str(mean(average_dist_bin(:,1))) ' [cm]'])
%disp(['Mean longer distance = ' num2str(mean(average_dist_bin(:,2))) ' [cm]'])
%disp('--------------------------------')
disp(['Variance explained by speed = ',num2str(mean(Dur_var_exp)),' +- ', num2str(std(Dur_var_exp))])

% figure
% plot(Dur_var_exp','.-')
% ylabel('Variance explained')
% box off
% xlim([0.8 3.2])
% ylim([0 100])
% 
% mean(Dur_var_exp)
% std(Dur_var_exp)
end

function plot_R2(R2,Area,isession)

if strcmp(Area,'M1')
    colourArea=[85 30 116]./256;
else
    colourArea=[89 156 153]./256;
end


subplot(2,3,3)

hold on
errorbar(isession,mean(R2.vel_same_dir),std(R2.vel_same_dir)/sqrt(numel(R2.vel_same_dir)),'.','Color',colourArea)
errorbar(isession,mean(R2.vel_other_dir),std(R2.vel_other_dir)/sqrt(numel(R2.vel_other_dir)),'.','Color',[0.5 0.5 0.5])

% subplot(2,3,3)
% hold on
% errorbar(isession,mean(R2.dist_same_dir),std(R2.dist_same_dir)/sqrt(numel(R2.dist_same_dir)),'.','Color',colourArea)
% errorbar(isession,mean(R2.dist_other_dir),std(R2.dist_other_dir)/sqrt(numel(R2.dist_other_dir)),'.','Color',[0.5 0.5 0.5])

end