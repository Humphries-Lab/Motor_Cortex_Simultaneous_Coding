%% aim: this script should call all functions to reproduce the results of the paper
close all 
clearvars
addpath(genpath('../'))  
%% Define recordings to analyse

session={'MC_S1_raw.mat','MC_S2_raw.mat','MC_S3_raw.mat','MC_S4_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat','MT_S3_raw.mat','MT_S2_raw.mat','MT_S1_raw.mat'};
Area={'M1','M1','M1','M1','M1','M1','M1','PMd','PMd','PMd','PMd','PMd','PMd'};
session_N=[1 2 3 4 5 6 7 5 6 7 8 9 10]; %Number of behavioural sessions

%% Parameters
Ndir=8; % Number of movement directions
Nbins=4; % Number of movement durations (100 ms bins)
threshold=80; % Percentage of variance explained by PCs
threshold_dist=10; %threshold of the distance (theta) for recurrence analyses


k_fold=6; %Number of folds for cross-validation for predicting movement direction  
Nrep=10; %Nunmber of repetitions for cross-validation for predicting movement direction
shuffle=0; %Perform shuffle of direction labels 
do_extra_plot=1;

mov_durS=(1:Nbins)/10+0.1; %Movement duration for each bin [S]

t_2b=[0.2*ones(sum(strcmp(Area,'M1')),1);0.05*ones(sum(strcmp(Area,'PMd')),1)];

% Run analyses pipeline

%% Figure 1
% % Example trial
%  trial=4;
% Behaviour_spikes_per_trial(session{11}, Area{11},Ndir,trial)
% % Example session stats
% Behaviour_variability(session{11}, Area{11},Ndir,Nbins)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2

% before selecting the end and start of each area, select neural activity from 500 ms before the
% movement onset upto 300 ms after the movement end

t_1=[-0.5*ones(1,sum(strcmp(Area,'M1'))),-0.5*ones(1,sum(strcmp(Area,'PMd')))];
t_2=[repmat(mov_durS+0.1+0.3,sum(strcmp(Area,'M1')),1);repmat(mov_durS+0.1+0.3,sum(strcmp(Area,'PMd')),1)];%movement duration+ 300 ms
do_plot=0;
plot_supp=0; %Don't plot supplementary

embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,do_plot,t_1,t_2,mov_durS,session_N,plot_supp);

[new_t_1,new_t_2]=Trajectories_differ_by_dir_all_sessions(session,Area,threshold,Ndir,Nbins);

% select the neural activity from the selected segments for each area
 t_1=[-0.25*ones(1,sum(strcmp(Area,'M1'))),-0.45*ones(1,sum(strcmp(Area,'PMd')))];
 t_2=[repmat(mov_durS+0.1+0.2,sum(strcmp(Area,'M1')),1);repmat(mov_durS+0.1+0.05,sum(strcmp(Area,'PMd')),1)];%movement duration+ 200 ms for M1 and 50 ms for PMd

do_plot=0;
plot_supp=1; %plot supplementary
embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,do_plot,t_1,t_2,mov_durS,session_N,plot_supp);

%% Figure 3

%recurrence_all_sessions(session,Area,threshold,Ndir,threshold_dist)

%recurrence_region_all_sessions(session,Area,threshold,Ndir)
 
%% Figure 4

% decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep,shuffle,do_extra_plot)
% shuffle=1;
% decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep,shuffle,do_extra_plot)


%% Figure 5
% [h,p]=distance_duration_vs_direction_all_sessions(session,Area,threshold,Ndir,Nbins,do_extra_plot);
% 
% speed_distance_all_sessions(session,Area,threshold,Ndir,t_1,t_2b)

%% Figure 6
%decoding_movement_duration_all_sessions(session,Area,Ndir,Nbins)

%% Figure 7
% RunDifferentTau