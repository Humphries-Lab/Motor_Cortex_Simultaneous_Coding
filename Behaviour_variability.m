function Behaviour_variability(Session,Area,Ndir,edges_dur_bin)
%% Behaviour_variability plots the distribution of each movement parameter
% (Direction, duration, speed and distance) in the Session
% Movements are binned into Ndir directions and Nbins durations
%
% INPUTS 
%
% Session: Name of the session to be analysed.
% e.g 'MC_S1_raw.mat'
%
% Area: Name of the area to be analysed in the Session.
% e.g 'M1'
%
% Ndir: number of directions to bin the movements
%
% edges_dur_bin= array containing the edges of each duration bin [S]
%
%
% Example
% Behaviour_variability('MC_S1_raw.mat','M1',8,[200 300 400 500 600])
%
% 23/01/2023
% Andrea Colins Rodriguez

Nbins=numel(edges_dur_bin)-1;
colour_dir=My_hsv(Ndir);
colour_dur=plasma(Nbins);
ms=1000; %to convert from s to ms
%these parameters are related to the neural activity. They don't really
%matter for the purpose of the analysis of the behaviour
sigma_filter=20;
t_from=-0.5;
t_upto=0.2;

duration_range=[0.05 2]; %select movements between 50 ms and 2s

% extract all movements and their parameters from the session
load(Session,Area,'trial_table2','cont')
if strcmp(Area,'PMd')
    neural_data=PMd.units;
else
    neural_data=M1.units;
end

[~,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto,duration_range);

%% Plot histogram of movement direction
subplot(4,4,13)
h=polarhistogram(Mov_params.direction,Ndir);
Edges=h.BinEdges;
direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));
% colour bars by direction
for i=1:Ndir
    polarhistogram(Mov_params.direction(direction1==i),Edges,'FaceColor',colour_dir(i,:));
    hold on
end
    
%% Plot histogram of movement duration  
subplot(4,4,14)
hold on
dur_binsize=(edges_dur_bin(2)-edges_dur_bin(1))*ms;

h=histogram(Mov_params.duration*ms,0:dur_binsize:max(Mov_params.duration)*ms,'normalization','probability','EdgeAlpha',0,'FaceColor','k');

middlebin=edges_dur_bin(1:Nbins)*ms+dur_binsize/2;

for i=1:numel(middlebin)
bar(middlebin(i),h.Values(i+round(edges_dur_bin(1)*ms/dur_binsize)),dur_binsize,'FaceColor',colour_dur(i,:))
end

box off
xlabel('Movement duration [ms]')
ylabel('Fraction of movements')
ylim([0 0.6])

selected_mov=Mov_params.duration<=0.6 & Mov_params.duration>=0.2;

%% Plot histogram of movement speed
subplot(4,4,15)
hold on
plot(Mov_params.distance(selected_mov),Mov_params.max_speed(selected_mov),'.k')
corrtex=corr(Mov_params.distance(selected_mov)',Mov_params.max_speed(selected_mov)');
text(0.7,0.1,['R = ' num2str(corrtex)],'Unit','Normalized')
box off
ylabel('Max speed [cm/s]')
xlabel('Distance [cm]')


%% Plot histogram of movement distance
subplot(4,4,16)
hold on
plot(Mov_params.distance(selected_mov),Mov_params.duration(selected_mov),'.k')
corrtex=corr(Mov_params.distance(selected_mov)',Mov_params.duration(selected_mov)');
text(0.7,0.1,['R = ' num2str(corrtex)],'Unit','Normalized')
xlabel('Distance [cm]')
ylabel('Duration [s]')
box off

end