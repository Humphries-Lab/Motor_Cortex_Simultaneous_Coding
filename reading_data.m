% add the path to matnwb and generate the core classes
addpath('C:\Users\controlmotor\Desktop\Andrea\codes_from_papers\matnwb-master')

% Reminder: YOU DO NOT NORMALLY NEED TO CALL THIS FUNCTION. Only attempt this method if you
% encounter read errors.
% generateCore(util.getSchemaVersion('sub-M_ses-CO-20140203_behavior+ecephys.nwb'));

% ignorecache informs the `nwbRead` call to not generate files by default. Since we have already
% done this, we can skip this automated step when reading. If you are reading the file before
% generating, you can omit this argument flag.
nwb = nwbRead('sub-M_ses-RT-20140224_behavior+ecephys.nwb', 'ignorecache');


%% Neural activity
Nunits = nwb.units.id.data.dims;

for iunit=1:Nunits
unit_spikes = nwb.units.getRow(iunit, 'columns', {'spike_times'}).spike_times{1};
units(iunit).ts=unit_spikes;
end
PMd.units=units;

% results = cell(1, length(stimulus_times));
% 
% for itime = 1:length(stimulus_times)
%     stimulus_time = stimulus_times(itime);
%     spikes = unit_spikes - stimulus_time;
%     spikes = spikes(spikes > -before);
%     spikes = spikes(spikes < after);
%     results{itime} = spikes;
% end

% figure();
% hold on
% for i = 1:length(results)
%     spikes = results{i};
%     yy = ones(length(spikes)) * i;
%  
%     plot(spikes, yy, 'k.');
% end
% hold off
% ylabel('trial');
% xlabel('time (s)');
% axis('tight')

%% behaviour
Pos=nwb.processing.get('behavior').nwbdatainterface.get('Position').spatialseries.get('cursor_pos').data;

%X=Pos(1,:)';
%Y=Pos(2,:)';

% Velocity
Vel=nwb.processing.get('behavior').nwbdatainterface.get('Velocity').timeseries.get('cursor_vel').data;
%Vx=Speed(1,:);
%Vy=Speed(2,:);


% Acceleration
Acc=nwb.processing.get('behavior').nwbdatainterface.get('Acceleration').timeseries.get('cursor_acc').data;
%Vx=Speed(1,:);
%Vy=Speed(2,:);
t=nwb.processing.get('behavior').nwbdatainterface.get('Position').spatialseries.get('cursor_pos').timestamps;


cont.t=t(1:end);
cont.pos=Pos(1:2,:)';
cont.vel=Vel(1:2,:)';
cont.acc=Acc(1:2,:)';


%% trial info
trial_start = nwb.intervals_trials.start_time.data;
trial_start=trial_start(1:end);

trial_end = nwb.intervals_trials.stop_time.data;
trial_end=trial_end(1:end);

Ntrials=length(trial_start);

%this table contain the trial info following the same format than in Lawlor
%2018(?)
trial_table2=nan(Ntrials,22);
trial_table2(:,1)=trial_start(:);
trial_table2(:,end)=trial_end(:);

% not positive about this one. Check!
target_onset=nwb.intervals_trials.vectordata.get('go_cue_time_array').data';

trial_table2(:,[2 7 12 17])=target_onset;



