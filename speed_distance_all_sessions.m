function speed_distance_all_sessions(session,Area,threshold,Ndir,t_1,t_2b)
Nsessions=size(session,2);
fraction_above_vel=zeros(Nsessions,1);
p_vel=zeros(Nsessions,1);
fraction_above_dist=zeros(Nsessions,1);
p_dist=zeros(Nsessions,1);
average_speed_bin=zeros(Nsessions,2);
average_dist_bin=average_speed_bin;
figure
for i=1:Nsessions
    if i==1
        do_plot=1;
    else
        do_plot=0;
    end
    [fraction_above_vel(i),p_vel(i),fraction_above_dist(i),p_dist(i),average_speed_bin(i,:),average_dist_bin(i,:)]=speed_distance_dimensions(session{i}, Area{i},threshold,Ndir,t_1(i),t_2b(i),do_plot);
    if strcmp(Area{i},'M1')
        colourArea='m';
    else
        colourArea='b';
    end
    
    subplot(2,5,3)
    hold on
    plot(i,fraction_above_vel(i),'o','Color',colourArea)
    
    subplot(2,5,5+3)
    hold on
    plot(i,fraction_above_dist(i),'o','Color',colourArea)
    pause(0.001)
end

subplot(2,5,3)
xlim([0 Nsessions+1])
ylim([0 1])
box off
xlabel('Recording Number')
ylabel('Fraction above (Speed)')
subplot(2,5,5+3)
xlim([0 Nsessions+1])
ylim([0 1])
box off
xlabel('Recording Number')
ylabel('Fraction above (Speed)')
p_vel
p_dist
p_vel_max=max(p_vel)
p_dist_max=max(p_dist)
average_speed_bin
mean(average_speed_bin)
mean(average_dist_bin)
average_dist_bin
end