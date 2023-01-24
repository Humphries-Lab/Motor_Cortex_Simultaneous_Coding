function Behaviour_variability(Session,Area,Ndir,Nbins)
%% Behaviour_variability plots the distribution of each movement parameter
% (Direction, duration, speed and distance)

% Behaviour_variability(Session,Area,Ndir,Nbins) plots the distribution of each movement parameter
% of the session Session and area Area. 
% Movements are binned into Ndir directions and Nbins durations

% Example
% Behaviour_variability('MC_S1_raw.mat','M1',8,4)
%
% 23/01/2023
% Andrea Colins Rodriguez

colour_dir=hsv(Ndir);
colour_dur=plasma(Nbins);

%these parameters are related to the neural activity. They don't really
%matter for the purpose of the analyses of the behaviour
event=2; %align neural data to movement onser
ms=20;% this number doesn't matter because its related to the neural activity, not the behaviour
t_1=-0.5;
t_2=0.5;

duration_range=[0.05 2]; %select movements between 50 ms and 2s

% extract all movements and their parameters from the session
[~,direction,~,~,~,~,dist_mov,mov_duration,max_speed]=neural_data_per_duration(Session,Area,ms,t_1,t_2,event,duration_range);


%% Plot histogram of movement direction
subplot(4,4,13)
h=polarhistogram(direction,Ndir);
Edges=h.BinEdges;
direction1=ceil(Ndir*(direction+pi)/(2*pi));
% colour bars by direction
for i=1:Ndir
    polarhistogram(direction(direction1==i),Edges,'FaceColor',colour_dir(i,:));
    hold on
end
    
%% Plot histogram of movement duration  
subplot(4,4,14)
hold on
h=histogram(mov_duration*1000,0:100:max(mov_duration)*1000,'normalization','probability','EdgeAlpha',0,'FaceColor','k');
middlebin=[250 350 450 550];
for i=1:numel(middlebin)
bar(middlebin(i),h.Values(i+2),100,'FaceColor',colour_dur(i,:))
end

box off
xlabel('Movement duration [ms]')
ylabel('Fraction of movements')
ylim([0 0.6])

%% Plot histogram of movement speed
subplot(4,4,15)
histogram(max_speed,0:1:max(max_speed),'normalization','probability','EdgeAlpha',0,'FaceColor','k')
box off
xlabel('Max speed [cm/s]')
ylim([0 0.1])


%% Plot histogram of movement distance
subplot(4,4,16)
histogram(dist_mov,0:1:max(dist_mov),'normalization','probability','EdgeAlpha',0,'FaceColor','k')
xlabel('Distance [cm]')
ylim([0 0.4])
box off

end