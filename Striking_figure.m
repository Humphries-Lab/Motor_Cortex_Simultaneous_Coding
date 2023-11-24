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
%mkdir('../Output_files')
%% Define recordings to analyse

session={'MC_S1_raw.mat'};
Area={'M1'};
session_N=[1]; %Number of behavioural sessions

%% Parameters
Ndir=8; % Number of movement directions
Nbins=2; % Number of movement durations
dur_bin_size=0.4; % duration bin size [s]
start_dur_bin=0.1; % minimum duration for the first duration bin [s]
threshold=80; % Percentage of variance explained by PCs
threshold_dist=10; %threshold of the distance (theta) for recurrence analyses

% decoding direction parameters
k_fold=6; %Number of folds for cross-validation for predicting movement direction
Nrep=10; %Number of repetitions for cross-validation for predicting movement direction

edges_dur_bin=(0:Nbins)*dur_bin_size+start_dur_bin; %Movement duration for each bin [S]

%% Run analyses pipeline
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2
% before selecting the end and start of each area, select neural activity from 500 ms before the
% movement onset upto 300 ms after the movement end

t_from=[-0.25*ones(1,sum(strcmp(Area,'M1'))),-0.45*ones(1,sum(strcmp(Area,'PMd')))];
t_upto=[0.2*ones(sum(strcmp(Area,'M1')),1);0.05*ones(sum(strcmp(Area,'PMd')),1)];%movement duration+ 200 ms for M1 and 50 ms for PMd

plot_supp=1; %plot supplementary

%embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin,session_N,plot_traj_all_rec,plot_supp);
[variance,mov_distance,mov_duration,max_speed]=striking_subfigure(session{1},Area{1},Ndir,Nbins,t_from,t_upto,edges_dur_bin,1)