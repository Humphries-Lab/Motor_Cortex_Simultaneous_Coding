function plot_kinematics(dist_mov_dir,mov_duration,max_speed,Nbins)
colour_dur=plasma(Nbins);
for i=1:Nbins
    
    subplot(2,3,1)
    hold on
    plot(dist_mov_dir{i},mov_duration{i},'.','Color',colour_dur(i,:))
    
    
    subplot(2,3,2)
    hold on
    plot(dist_mov_dir{i},max_speed{i},'.','Color',colour_dur(i,:))
    
    
    subplot(2,3,3)
    hold on
    plot(max_speed{i},mov_duration{i},'.','Color',colour_dur(i,:))
    
    
end
subplot(2,3,1)
R_dist_dur=corr([dist_mov_dir{1:Nbins}]',[mov_duration{1:Nbins}]');
text(12,0.25,['Mean = ' num2str(mean(R_dist_dur),'%.2f')])
xlabel('Distance [cm]')
ylabel('Duration [s]')

subplot(2,3,2)
R_dist_speed=corr([dist_mov_dir{1:Nbins}]',[max_speed{1:Nbins}]');
text(12,5,['Mean = ' num2str(mean(R_dist_speed),'%.2f')])
xlabel('Distance [cm]')
ylabel('Max speed [cm/s]')

subplot(2,3,3)
R_speed_dur=corr([max_speed{1:Nbins}]',[mov_duration{1:Nbins}]');
text(60,0.22,['Mean = ' num2str(mean(R_speed_dur),'%.2f')])
xlabel('Max speed [cm/s]')
ylabel('Duration [s]')
end