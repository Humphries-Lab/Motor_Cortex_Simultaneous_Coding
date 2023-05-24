%% Script to reproduce the results of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: A few sessions contain simultaneous recordings from M1 and PMd.
%
% We define:
%
% Session: Behavioural session. A file that contain one set of movements
% and their kinematic parameters. e.g the file MM_S1_raw.mat
%
% Area: Neural recording corresponding to one area. e.g. 'M1' or 'PMd'
%
% Recording: A recording is defined by a session and an area.

% This way ('MM_S1_raw.mat','M1') and ('MM_S1_raw.mat','PMd') are two
% recordings from the same behavioural session.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clearvars
addpath(genpath('../'))
mkdir('../Output_files')
%% Define recordings to analyse

session={'MC_S1_raw.mat','MC_S2_raw.mat','MC_S3_raw.mat','MC_S4_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat','MT_S3_raw.mat','MT_S2_raw.mat','MT_S1_raw.mat'};
Area={'M1','M1','M1','M1','M1','M1','M1','PMd','PMd','PMd','PMd','PMd','PMd'};
session_N=[1 2 3 4 5 6 7 5 6 7 8 9 10]; %Number of behavioural sessions

%% Parameters
Ndir=8; % Number of movement directions
Nbins=4; % Number of movement durations
dur_bin_size=0.1; % duration bin size [s]
start_dur_bin=0.2; % minimum duration for the first duration bin [s]
threshold=80; % Percentage of variance explained by PCs
threshold_dist=10; %threshold of the distance (theta) for recurrence analyses

% decoding direction parameters
k_fold=6; %Number of folds for cross-validation for predicting movement direction
Nrep=10; %Number of repetitions for cross-validation for predicting movement direction

edges_dur_bin=(0:Nbins)*dur_bin_size+start_dur_bin; %Movement duration for each bin [S]

%% Run analyses pipeline

%% Figure 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example trial
trial=4;
Behaviour_spikes_per_trial(session{11}, Area{11},Ndir,trial)
% Example session stats
Behaviour_variability(session{11}, Area{11},Ndir,edges_dur_bin)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2
% before selecting the end and start of each area, select neural activity from 500 ms before the
% movement onset upto 300 ms after the movement end

t_from=-0.5*ones(1,size(session_N,2));% from 500 ms before mov onset
t_upto=0.3*ones(1,size(session_N,2)); %movement duration+ 300 ms

plot_traj_all_rec=0; % plot trajectories for all recordings
plot_supp=0; %Don't plot supplementary yet

embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin,session_N,plot_traj_all_rec,plot_supp);
[new_t_1,new_t_2]=Trajectories_differ_by_dir_all_sessions(session,Area,threshold,Ndir,Nbins);

%% select the neural activity from the selected segments for each area

t_from=[-0.25*ones(1,sum(strcmp(Area,'M1'))),-0.45*ones(1,sum(strcmp(Area,'PMd')))];
t_upto=[0.2*ones(sum(strcmp(Area,'M1')),1);0.05*ones(sum(strcmp(Area,'PMd')),1)];%movement duration+ 200 ms for M1 and 50 ms for PMd

plot_supp=1; %plot supplementary

embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin,session_N,plot_traj_all_rec,plot_supp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3
do_plot_supp=1;
i_bin=1; % plot results for this bin only
recurrence_all_sessions(session,Area,threshold,Ndir,i_bin,threshold_dist,do_plot_supp)
recurrence_region_all_sessions(session,Area,threshold,Ndir,Nbins,threshold_dist)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 4
warning('off','stats:pca:ColRankDefX') % do not show the warning "Columns of X are linearly dependent to within machine precision."
decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep)
warning('on','stats:pca:ColRankDefX') % back on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 5
upper_bound_similarity(session,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin)
do_plot_supp=1;
%fig2=figure;
distance_duration_vs_direction_all_sessions(session,Area,threshold,Ndir,Nbins,do_plot_supp);
speed_distance_all_sessions(session,Area,threshold,Ndir,t_from,t_upto)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 6
decoding_movement_duration_all_sessions(session,Area,threshold,Ndir,Nbins)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 7
InputRange=[-1,1];
RunDifferentTau(InputRange)

% for supp figure
InputRange=[0,1];
RunDifferentTau(InputRange)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%