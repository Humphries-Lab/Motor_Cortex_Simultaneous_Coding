function speed_distance_all_sessions(session,Area,threshold,Ndir,t_from,t_upto)
%% speed_distance_all_sessions compares the distance between
%% trajectories corresponding to movement with different speed (same
%% direction and distance) and movemements of adjacent bin directions (same speed)
%% It then repeats the process for distance (same speed) vs directions 
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
p_dist=zeros(Nsessions,1);
average_speed_bin=zeros(Nsessions,2);
average_dist_bin=average_speed_bin;
Distance_vel=nan(Nsessions,Ndir);
Distance_dist=nan(Nsessions,Ndir);
Distance_v_dir=nan(Nsessions,32);
Distance_d_dir=nan(Nsessions,32);

for isession=1:Nsessions

    [fraction_above_vel(isession),p_vel(isession),fraction_above_dist(isession),p_dist(isession),average_speed_bin(isession,:),average_dist_bin(isession,:),Distance_vel(isession,:),Distance_v_dir(isession,:),Distance_dist(isession,:),Distance_d_dir(isession,:)]=speed_distance(session{isession}, Area{isession},threshold,Ndir,t_from(isession),t_upto(isession));

end
subplot(2,3,5)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance speed')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])
text(0.005,0.02,[' max p-value = ' num2str(max(p_vel))],'FontSize',8)


subplot(2,3,6)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance distance')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])
text(0.005,0.02,[' max p-value = ' num2str(max(p_dist))],'FontSize',8)
disp('--------------------------------')
disp('Recordings which p-values > 0.05 (Distance vs direction, FIG 5F, right)')
Session_name=session(p_dist>0.05)';
Area_name=Area(p_dist>0.05)';
pvalue_dist=p_dist(p_dist>0.05);
pvalues=table(Session_name,Area_name,pvalue_dist)

disp('--------------------------------')


axes('Position',[0.4116    0.47    0.21   0.041]);
histogram(Distance_vel(:),'Normalization','probability')
xlim([0 0.025])
xlabel('Distance Speed')
box off

axes('Position',[ 0.6331    0.1116    0.0164    0.34]);
histogram(Distance_v_dir(:),'Normalization','probability','orientation','horizontal')
ylim([0 0.025])
ylabel('Distance Closest direction')
box off

axes('Position',[ 0.6936    0.47    0.2100    0.0410]);
histogram(Distance_dist(:),'Normalization','probability')
xlim([0 0.025])
box off

axes('Position',[ 0.9152    0.1187    0.0164    0.3380]);
histogram(Distance_d_dir(:),'Normalization','probability','orientation','horizontal')
ylim([0 0.025])
box off

 disp('--------------------------------')
 disp(['Mean lower speed = ' num2str(mean(average_speed_bin(:,1))) ' [cm/s]'])
 disp(['Mean higher speed = ' num2str(mean(average_speed_bin(:,2))) ' [cm/s]'])
 disp('--------------------------------')
 disp(['Mean shorter distance = ' num2str(mean(average_dist_bin(:,1))) ' [cm]'])
 disp(['Mean longer distance = ' num2str(mean(average_dist_bin(:,2))) ' [cm]'])
 disp('--------------------------------')


end