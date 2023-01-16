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

figure
for i=1:Nsessions
    if i==1
        do_plot=1;
    else
        do_plot=0;
    end
    [fraction_above_vel(i),p_vel(i),fraction_above_dist(i),p_dist(i),average_speed_bin(i,:),average_dist_bin(i,:),Distance_vel(i,:),Distance_v_dir(i,:),Distance_dist(i,:),Distance_d_dir(i,:)]=speed_distance_dimensions(session{i}, Area{i},threshold,Ndir,t_1(i),t_2b(i),do_plot);
    if strcmp(Area{i},'M1')
        colourArea='m';
    else
        colourArea='b';
    end

end
subplot(2,5,1)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance speed')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])

subplot(2,5,5+1)
plot([0 0.025],[0 0.025],'k')
xlabel('Distance distance')
ylabel('Distance closest directions')
xlim([0 0.025])
ylim([0 0.025])

subplot(2,5,3)
histogram(Distance_vel(:),'Normalization','probability')
xlim([0 0.025])
xlabel('Distance Speed')
box off

subplot(2,5,4)
histogram(Distance_v_dir(:),'Normalization','probability','orientation','horizontal')
ylim([0 0.02])
ylabel('Distance Closest direction')
box off


subplot(2,5,3+5)
histogram(Distance_dist(:),'Normalization','probability')
xlim([0 0.025])
box off
xlabel('Distance distance')

subplot(2,5,4+5)
histogram(Distance_d_dir(:),'Normalization','probability','orientation','horizontal')
ylim([0 0.025])
box off
ylabel('Distance closest directions')

p_vel
p_dist
p_vel_max=max(p_vel)
p_dist_max=max(p_dist)
average_speed_bin
mean(average_speed_bin)
mean(average_dist_bin)
average_dist_bin
end