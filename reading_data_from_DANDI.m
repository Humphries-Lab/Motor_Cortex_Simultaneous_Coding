function reading_data_from_DANDI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reading_data_from_DANDI is a function that converts the data
% from the DANDI repository into Matlab structures to perform subsequent
% analyses
%
% 30/01/2025
% Andrea Colins Rodriguez
%
% Please define the path where the MATNWB toolbox is stored
MATNWB_path='C:\Users\controlmotor\Desktop\Andrea\codes_from_papers\matnwb-master';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add the path to matnwb and generate the core classes
addpath(MATNWB_path)
[Dandi_name,rec_name]=index_files;
counter=1;
thr_dist=1;

for i_file=1:size(Dandi_name,2)
    counter_nan=0;

    nwb = nwbRead(Dandi_name{i_file});


    %% Neural activity
    %Nunits = nwb.units.id.data.dims;

    Loc=nwb.general_extracellular_ephys_electrodes.vectordata.get('location').data(:);
    Loc_M1=find(strcmp(Loc,'Primary Motor Cortex'));
    Loc_PMd=find(strcmp(Loc,'Dorsal Premotor Cortex'));

    El=nwb.units.electrodes.data(:);
    %% asign M1 units
    idx_M1 = find(ismember(El,Loc_M1));
    if ~isempty(idx_M1)
        NM1=numel(idx_M1);

        for iunit=1:NM1
            unit_spikes = nwb.units.getRow(idx_M1(iunit), 'columns', {'spike_times'}).spike_times{1};
            units(iunit).ts=unit_spikes;
        end
        M1.units=units;
        clear units
    end

    %% asign PMd units
    idx_PMd = find(ismember(El,Loc_PMd));

    if ~isempty(idx_PMd)
        NPMd=numel(idx_PMd);

        for iunit=1:NPMd
            unit_spikes = nwb.units.getRow(idx_PMd(iunit), 'columns', {'spike_times'}).spike_times{1};
            units(iunit).ts=unit_spikes;
        end
        PMd.units=units;
        clear units
    end



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

    % Velocity
    Vel=nwb.processing.get('behavior').nwbdatainterface.get('Velocity').timeseries.get('cursor_vel').data;

    % Acceleration
    Acc=nwb.processing.get('behavior').nwbdatainterface.get('Acceleration').timeseries.get('cursor_acc').data;
    t=nwb.processing.get('behavior').nwbdatainterface.get('Position').spatialseries.get('cursor_pos').timestamps;

    %% interpolate because for some reason the kinematics is formatted in 10 ms bin.
    %this conversion works just fine for Mihili but not for Chewie
    Pos=interp1(t(1:end),Pos(1:2,:)',0:0.001:t(end));
    Vel=interp1(t(1:end),Vel(1:2,:)',0:0.001:t(end));
    Acc=interp1(t(1:end),Acc(1:2,:)',0:0.001:t(end));
    t2=interp1(t(1:end),t(1:end),0:0.001:t(end));

    startt=find(t2>=1,1,'first');

    cont.t=t2(startt:end);
    cont.pos=Pos(startt:end,:);
    cont.vel=Vel(startt:end,:);
    cont.acc=Acc(startt:end,:);


    %% trial info
    trial_start = nwb.intervals_trials.start_time.data;
    trial_start=trial_start(1:end);

    trial_end = nwb.intervals_trials.stop_time.data;
    trial_end=trial_end(1:end);

    Ntrials=length(trial_start);

    %this table contains the trial info following the same format than in Lawlor
    %2018(?)
    trial_table2=nan(Ntrials,22);
    trial_table2(:,1)=trial_start(:);
    trial_table2(:,end)=trial_end(:);

    % not positive about this one. Check!
    target_onset=nwb.intervals_trials.vectordata.get('go_cue_time_array').data;

    trial_table2(:,[2 7 12 17])=target_onset(:,:)';

    %[3 8 13 18] mov onset
    %[4 9 14 19] time of peak velocity
    %[5 10 15 20] X position of targets 1 to 4
    %[6 11 16 21] Y position of targets 1 to 4

    %% target positions
    % 1
    X_1=nwb.intervals_trials.vectordata.get('target_1_x_position').data;
    X_2=nwb.intervals_trials.vectordata.get('target_2_x_position').data;
    X_3=nwb.intervals_trials.vectordata.get('target_3_x_position').data;
    X_4=nwb.intervals_trials.vectordata.get('target_4_x_position').data;

    trial_table2(:,[5 10 15 20])=[X_1(:) X_2(:) X_3(:) X_4(:)];


    Y_1=nwb.intervals_trials.vectordata.get('target_1_y_position').data;
    Y_2=nwb.intervals_trials.vectordata.get('target_2_y_position').data;
    Y_3=nwb.intervals_trials.vectordata.get('target_3_y_position').data;
    Y_4=nwb.intervals_trials.vectordata.get('target_4_y_position').data;

    trial_table2(:,[6 11 16 21])=[Y_1(:) Y_2(:) Y_3(:) Y_4(:)];
    idx_tar=[2 7 12 17];
    target_pos=abs(trial_table2(:,idx_tar+3))+abs(trial_table2(:,idx_tar+4));

    % some of the targets' position have impossible number like 1e23? ignore
    % those trials

    %% mov onset
    speed=sqrt(Vel(:,1).^2+Vel(:,2).^2);

    % plot(cont.t,speed)
    % hold on
    idx_next=[7 12 17 22];

    for i_trial=1:Ntrials
        for i_target=1:4
            T_app=trial_table2(i_trial,idx_tar(i_target));
            T_next=trial_table2(i_trial,idx_next(i_target));

            if ~isnan(T_app)
                idx_T=find(t2>=T_app,1,'first');
                idx_T2=find(t2>=T_next,1,'first');

                % define mov onset
                idx_onset=find(speed(idx_T:idx_T2)>=8,1,'first');

                if idx_onset>1
                    t_onset=t2(idx_T+idx_onset-1);

                    %% find the time when cursor reach target area

                    xtarget=trial_table2(i_trial,idx_tar(i_target)+3);
                    ytarget=trial_table2(i_trial,idx_tar(i_target)+4);


                    tmp_x=(Pos(idx_T:idx_T2,1)>=xtarget-thr_dist) &(Pos(idx_T:idx_T2,1)<=xtarget+thr_dist);
                    tmp_y=(Pos(idx_T:idx_T2,2)>=ytarget-thr_dist)&(Pos(idx_T:idx_T2,2)<=ytarget+thr_dist);


                    idx_offset=find((tmp_x+tmp_y)>1,1,'first');

                    % check that speed it's ok
                    % find time of max speed
                    if ~isnan(idx_T2)
                        MaxS=max(speed(idx_T+idx_onset:idx_T2));

                    end

                    if ~isempty(idx_offset) && MaxS<200
                        t_offset=t2(idx_T+idx_offset);
                        trial_table2(i_trial,idx_tar(i_target)+2)=t_offset;
                    end


                else
                    t_onset=nan;
                    counter_nan=counter_nan+1;
                end

                trial_table2(i_trial,idx_tar(i_target)+1)=t_onset;

                %plot([T_app T_app],[0 16],'k')
                %plot([t_onset t_onset],[0 8],'r')
                %plot([t_maxspeed t_maxspeed],[0 25],'b')
            end

        end
    end

    if ~isempty(idx_PMd) && ~isempty(idx_M1)
        save([rec_name{i_file} '_raw.mat'],'trial_table2','cont','PMd','M1')
    elseif ~isempty(idx_PMd)
        save([rec_name{i_file} '_raw.mat'],'trial_table2','cont','PMd')
    else
        save([rec_name{i_file} '_raw.mat'],'trial_table2','cont','M1')
    end

    counter=counter+1;
    %counter_nan

    clear PMD M1 units trial_table2 cont
    toc
end

end

function [Dandi_name,rec_name]=index_files
Dandi_name={'sub-M_ses-RT-20140116_behavior+ecephys.nwb',...
    'sub-M_ses-RT-20140214_behavior+ecephys.nwb',...
    'sub-M_ses-RT-20140221_behavior+ecephys.nwb',...
    'sub-C_ses-RT-20150318_behavior+ecephys.nwb',...
    'sub-C_ses-RT-20150320_behavior+ecephys.nwb',...
    'sub-C_ses-RT-20150317_behavior+ecephys.nwb',...
    'sub-C_ses-RT-20150316_behavior+ecephys.nwb',...
    'sub-T_ses-RT-20130904_behavior+ecephys.nwb',...
    'sub-T_ses-RT-20130906_behavior+ecephys.nwb',...
    'sub-T_ses-RT-20130910_behavior+ecephys.nwb'};

rec_name={'MM_S1',...
    'MM_S2',...
    'MM_S3',...
    'MC_S1',...
    'MC_S2',...
    'MC_S3',...
    'MC_S4',...
    'MT_S1',...
    'MT_S2',...
    'MT_S3'};
end
