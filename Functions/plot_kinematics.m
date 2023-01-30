function plot_kinematics(mov_dist,mov_duration,max_speed,Nbins)
%% plot_kinematics plots the correlation between movement distance,duration and distance
%
% INPUTS
% mov_dist: a cell array containing Nbins cells. Each cell contains a vector with the
% distance of each movement within a duration bin.
%
% mov_duration: a cell array containing Nbins cells. Each cell contains a vector with the
% duration of each movement within a duration bin.
%
% max_speed: a cell array containing Nbins cells. Each cell contains a vector with the
% maximum speed reached in each movement within a duration bin.
% 
% Nbins: number of movement duration bins. 

% 23/01/2023
% Andrea Colins Rodriguez

colour_dur=plasma(Nbins);

for i=1:Nbins
    
    subplot(2,3,1)
    hold on
    plot(mov_dist{i},mov_duration{i},'.','Color',colour_dur(i,:))
    
    
    subplot(2,3,2)
    hold on
    plot(mov_dist{i},max_speed{i},'.','Color',colour_dur(i,:))
    
    
    subplot(2,3,3)
    hold on
    plot(max_speed{i},mov_duration{i},'.','Color',colour_dur(i,:))
    
    
end
subplot(2,3,1)
R_dist_dur=corr([mov_dist{1:Nbins}]',[mov_duration{1:Nbins}]');
text(12,0.25,['Corr = ' num2str(mean(R_dist_dur),'%.2f')])
xlabel('Distance [cm]')
ylabel('Duration [s]')

subplot(2,3,2)
R_dist_speed=corr([mov_dist{1:Nbins}]',[max_speed{1:Nbins}]');
text(12,5,['Corr = ' num2str(mean(R_dist_speed),'%.2f')])
xlabel('Distance [cm]')
ylabel('Max speed [cm/s]')

subplot(2,3,3)
R_speed_dur=corr([max_speed{1:Nbins}]',[mov_duration{1:Nbins}]');
text(60,0.22,['Corr = ' num2str(mean(R_speed_dur),'%.2f')])
xlabel('Max speed [cm/s]')
ylabel('Duration [s]')
end