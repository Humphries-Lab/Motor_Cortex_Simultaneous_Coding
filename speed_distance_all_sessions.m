function speed_distance_all_sessions(session,Area,threshold,Ndir,t_1,t_2b)
Nsessions=size(session,2);
fraction_above_vel=zeros(Nsessions,1);
p_vel=zeros(Nsessions,1);
fraction_above_dist=zeros(Nsessions,1);
p_dist=zeros(Nsessions,1);
average_speed_bin=zeros(Nsessions,2);
average_dist_bin=average_speed_bin;
Distance_vel=nan(Nsessions,8);
Distance_dist=nan(Nsessions,8);
Distance_v_dir=nan(Nsessions,32);
Distance_d_dir=nan(Nsessions,32);

for isession=1:Nsessions

    [fraction_above_vel(isession),p_vel(isession),fraction_above_dist(isession),p_dist(isession),average_speed_bin(isession,:),average_dist_bin(isession,:),Distance_vel(isession,:),Distance_v_dir(isession,:),Distance_dist(isession,:),Distance_d_dir(isession,:)]=speed_distance_dimensions(session{isession}, Area{isession},threshold,Ndir,t_1(isession),t_2b(isession));

end
subplot(2,3,5)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance speed')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])

subplot(2,3,6)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance distance')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])

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


% p_vel
% p_dist
% p_vel_max=max(p_vel)
% p_dist_max=max(p_dist)
% average_speed_bin
% mean(average_speed_bin)
% mean(average_dist_bin)
% average_dist_bin
end