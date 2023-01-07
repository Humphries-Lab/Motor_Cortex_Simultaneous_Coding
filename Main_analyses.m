%% aim: this script should call all functions to reproduce the results of the paper
close all 
clear all
session={'MC_S1_raw.mat','MC_S2_raw.mat','MC_S3_raw.mat','MC_S4_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat',...
    'MM_S1_raw.mat','MM_S2_raw.mat','MM_S3_raw.mat','MT_S3_raw.mat','MT_S2_raw.mat','MT_S1_raw.mat'};
Area={'M1','M1','M1','M1','M1','M1','M1','PMd','PMd','PMd','PMd','PMd','PMd'};

%parameters
do_extra_plot=1;
threshold=80;
Ndir=8;
Nbins=4;
threshold_dist=10;
k_fold=6;
Nrep=10;
shuffle=0;

from=(1:Nbins)/10+0.1;
t_1=[-0.25*ones(1,sum(strcmp(Area,'M1'))),-0.45*ones(1,sum(strcmp(Area,'PMd')))];
t_2=[repmat(from+0.1+0.2,sum(strcmp(Area,'M1')),1);repmat(from+0.1+0.05,sum(strcmp(Area,'PMd')),1)];%movement duration+ 200 ms for M1 and 50 ms for PMd

%t_1=[-0.5*ones(1,sum(strcmp(Area,'M1'))),-0.5*ones(1,sum(strcmp(Area,'PMd')))];
%t_2=[repmat(from+0.1+0.3,sum(strcmp(Area,'M1')),1);repmat(from+0.1+0.3,sum(strcmp(Area,'PMd')),1)];%movement duration+ 300 ms



t_2b=[0.2*ones(sum(strcmp(Area,'M1')),1);0.05*ones(sum(strcmp(Area,'PMd')),1)];

%Run analyses

% Figure 1

%Figure 2

%embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,0,t_1,t_2,from);

% Figure 3
%%%%%%%%distance_position_all_sessions(session,Area,threshold,Ndir,Nbins,1,t_1,t_2,from);
%[t_1,t_2]=Trajectories_differ_by_dir_all_sessions(session,Area,threshold,Ndir,Nbins);

%embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,1,t_1,t_2,from)
%variance_trajectories_v2_all_sessions(session,Area,threshold,Ndir,Nbins,1,t_1,t_2,from);
recurrence_all_sessions(session,Area,threshold,Ndir,threshold_dist)
 
recurrence_region_all_sessions(session,Area,threshold,Ndir,do_extra_plot)
  
[h,p]=distance_duration_vs_direction_all_sessions(session,Area,threshold,Ndir,Nbins,do_extra_plot);
%speed_distance_all_sessions(session,Area,threshold,Ndir,t_1,t_2b)
%speed_distance_all_sessions_pos(session,Area,threshold,Ndir,t_1,t_2b)
% 
%decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep,shuffle,do_extra_plot)
%shuffle=1;
%decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep,shuffle,do_extra_plot)
%decoding_instantaneuos_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,Nrep,shuffle,1)

%decoding_movement_duration_all_sessions(session,Area,Ndir,Nbins)