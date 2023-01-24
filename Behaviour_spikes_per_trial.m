function Behaviour_spikes_per_trial(Session,Area,Ndir,Ntrial)
%% Behaviour_spikes_per_trial plots the timing, movement kinematic and neural activity of a selected trial

% Behaviour_spikes_per_trial(session,Area,Ndir,Ntrial) plots:
% 1) The hand speed of the trial Ntrial of the beavioural session Session.
% The times of the 4 movements are highlighted in the colours corresponding
% to their directions. Directions are binned into Ndir bins.
% 2) The raster plot of the population from the area Area (M1 or PMd)
% 3) The average firing rate across the population

% Red lines indicate movement onset, grey lines indicate target onset
% Example
% 
% Behaviour_spikes_per_trial('MT_S3_raw.mat','PMd',8,4)
%
% 23/01/2023
% Andrea Colins Rodriguez

load(Session,'trial_table2','cont',Area)

if strcmp(Area,'M1')
    %if selecting area M1, create a copy of the neural activity info called
    %PMd
    PMd=M1;
    clear M1
end

colour_dir=hsv(Ndir);
startt=trial_table2(Ntrial,1);
endt=trial_table2(Ntrial,22);


%calculate targets' position
ntarget_time=trial_table2(Ntrial,[2 7 12 17]);
idx_target=round(1000*(ntarget_time-1))+1;
xfirst=cont.pos(idx_target(1),1);
yfirst=cont.pos(idx_target(1),2);
ntargetx=[xfirst trial_table2(Ntrial,[5 10 15 20])];
ntargety=[yfirst trial_table2(Ntrial,[6 11 16 21])];

%speed
idx_start=round(1000*(startt-1))+1;
idx_end=round(1000*(endt-1))+1;
vel=sqrt(cont.vel(idx_start:idx_end,1).^2+cont.vel(idx_start:idx_end,2).^2);

% Define movement onset and movement end
nonset_time=round((trial_table2(Ntrial,[3 8 13 18])-startt)*1000);
nend_time=round((trial_table2(Ntrial,[3 8 13 18]+1)-startt)*1000);


figure

subplot(4,1,1)
plot((1:numel(vel))/1000,vel,'k')
hold on
ylabel('Hand speed [cm/s]')
title(['Trial number=' num2str(Ntrial)])
box off


matrix=spikest2vector(PMd.units,startt,endt);
filter=[zeros(1,20), ones(1,20)];
sub_psth=conv(mean(matrix)*1000,filter,'same');

subplot(4,1,3)
plot((1:numel(sub_psth))/1000,sub_psth,'k')% mean(matrix)*1000 in Hz
hold on
box off

for i=2:5
    direction=find_upcoming_direction([ntargetx(i-1) ntargety(i-1)],[ntargetx(i) ntargety(i)]);
    dir_color=ceil(Ndir*(direction+pi)/(2*pi));
    
    subplot(4,1,1)
    plot((nonset_time(i-1):nend_time(i-1))/1000,vel(nonset_time(i-1):nend_time(i-1)),'Color',colour_dir(dir_color,:))
    
    subplot(4,1,3)
    plot((nonset_time(i-1):nend_time(i-1))/1000,sub_psth(nonset_time(i-1):nend_time(i-1)),'Color',colour_dir(dir_color,:))
    
end

xlabel('Time [s]')

% Rasterplot
subplot(4,1,2)
title(Area)
for unit=1:size(PMd.units,2)
    
    idxsp=find(((PMd.units(unit).ts>startt) & (PMd.units(unit).ts<endt)));
    
    tsp=PMd.units(unit).ts(idxsp)-startt;
    
    if ~isempty(idxsp)
        plot(tsp,repmat(unit,1,numel(idxsp)),'.k')
        hold on
    end
end

ntarget_time=trial_table2(Ntrial,[2 7 12 17])-startt;
nonset_time=trial_table2(Ntrial,[3 8 13 18])-startt;

% draw target and onset times
for i=1:4
    plot([ntarget_time(i) ntarget_time(i)],[0 size(PMd.units,2)],'Color',[0.5 0.5 0.5])
    plot([nonset_time(i) nonset_time(i)],[0 size(PMd.units,2)],'r')
    
end
box off
xlabel('Time [s]')
ylabel('Units')


end
