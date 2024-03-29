function decoding_movement_direction_all_sessions_v2(Sessions,Areas,threshold,Ndir,k_fold,Nrep)
%% decoding_movement_direction_all_sessions predicts the movement direction bin and angle from the preparatory neural activity
%
% INPUTS
%
% Sessions: cell array containing the names of the sessions to be analysed.
% e.g {'MC_S1_raw.mat','MC_S2_raw.mat'}
%
% Areas: cell array containing the names of the areas to be analysed in
% each session. Areas and Sessions must have the same number of elements.
% e.g {'M1','M1'}
%
% threshold: percentage of the variance to be explained by the first nPCs
%
% Ndir: number of directions to bin the movements
%
% k_fold: number of folds for Cross-Validation (CV)
%
% Nrep: Number of repetitions of CV
%
% 24/05/2023
% Andrea Colins Rodriguez

%% Overview
%   step 1: extract the information from all trials
%   step 2: each movement has an ID. Divide the total number of
%   trials into Ndir bins. Select a random sample of the movements so that all bins have the same number
%   of samples.
%   step 3: for those selected movements, split the data to do CV
%   step 4: train the naive bayes model inside CV
% step 5: test using 1 class prediction accuracy and the estimation of the
% angle

colour_dir=hsv(Ndir);
median_error=zeros(size(Sessions,2),1);
median_error_angle_bin=zeros(size(Sessions,2),1);
std_error_bin=zeros(size(Sessions,2),1);
std_error=zeros(size(Sessions,2),1);
acc_error_mean=zeros(size(Sessions,2),1);
acc_error_std=zeros(size(Sessions,2),1);
Colour_M1=[85 30 116]./256;
Colour_PMd=[89 156 153]./256;
M1=strcmp(Areas,'M1');
Colour_Area(M1,:)=repmat(Colour_M1,sum(M1),1);
Colour_Area(~M1,:)=repmat(Colour_PMd,sum(~M1),1);
colour_shuffle=[0.5 0.5 0.5];

figure

for isession=1:size(Sessions,2)
    disp(['Starting decoding direction for recording ' num2str(isession)])
    %% Step 0 
    load(['../Output_files/PCA_' Sessions{isession}(1:end-4) '_' Areas{isession} '.mat'],'coeffs','variance','delete_units','normalisation')
    ndims=find(cumsum(variance)>threshold,1,'First');
    coeffs=coeffs(:,1:ndims);

    %% Step 1
    [angle_dir,idx_dir,neural_mov]=extract_movement_for_decoding(Sessions{isession}, Areas{isession},Ndir);
    
    %% Step 2
    Nsample=zeros(Ndir,1);
    
    for i_dir=1:Ndir
        Nsample(i_dir)=sum(idx_dir==i_dir);
    end
    
    New_N=min(Nsample);
    selected_angle=nan(New_N*Ndir,1);
    selected_class=nan(New_N*Ndir,1);
    selected_neural=nan(New_N*Ndir,size(coeffs,1));
    
    % transform FR to a matrix [n mov, units]
    %neural_mov=squeeze(mean(neural_mov,2))';
    neural_mov=squeeze(neural_mov(:,end,:))';
    % preprocess data for PCA
    % soft normalization
    neural_mov(:,delete_units)=[];
    neural_mov= neural_mov./repmat(normalisation,size(neural_mov,1),1);
    xangle=zeros(1,Ndir);
    
    for i_dir=1:Ndir
        idx_i=find(idx_dir==i_dir);
        %randon permutation to break temporal patterns throughout the session
        p_tmp = randperm(numel(idx_i));
        selected_idx=idx_i(p_tmp(1:New_N));
        
        selected_angle((i_dir-1)*New_N+1:New_N*i_dir)=angle_dir(selected_idx)';
        selected_class((i_dir-1)*New_N+1:New_N*i_dir)=idx_dir(selected_idx)';
        selected_neural((i_dir-1)*New_N+1:New_N*i_dir,:)=neural_mov(selected_idx,:);
        
        xangle(i_dir)=mean(angle_dir(selected_idx));
        
    end
    
    %% Step 3
    %%%% repeat CV Nrep times to make a robust statistic analysis
    shuffle=0;
    do_plot=isession==1;

    [Acc_per_cv,error_angle,error_angle_bin]=decoding_movement_direction(coeffs,selected_class,selected_angle,selected_neural,colour_dir,k_fold,Nrep,shuffle,do_plot);
    
    acc_error_mean(isession)=mean(Acc_per_cv(:));
    acc_error_std(isession)=std(Acc_per_cv(:));
    median_error_angle_bin(isession)=median(error_angle_bin);
    std_error_bin(isession)=std(error_angle_bin);
    median_error(isession)=median(error_angle);
    std_error(isession)=std(error_angle);
    
    subplot(2,3,6)
    hold on
    errorbar(isession+0.1,median_error(isession),std_error(isession),'.','Color',Colour_Area(isession,:))
    subplot(2,3,5)
    hold on
    errorbar(isession,acc_error_mean(isession),acc_error_std(isession),'.','Color',Colour_Area(isession,:))
    
    %% repeat above for shuffle test
    shuffle=1;
    do_plot=0;
    [Acc_per_cv_sh,error_angle_sh]=decoding_movement_direction(coeffs,selected_class,selected_angle,selected_neural,colour_dir,k_fold,Nrep,shuffle,do_plot);
    subplot(2,3,6)
    hold on
    errorbar(isession-0.1,median(error_angle_sh),std(error_angle_sh),'.','Color',colour_shuffle)
    
    subplot(2,3,5)
    hold on
    errorbar(isession,mean(Acc_per_cv_sh(:)),std(Acc_per_cv_sh(:)),'Color',colour_shuffle)
    pause(0.01)
end


%% Summary for all sessions

subplot(2,3,6)
hold on
box off
xlabel('Session')
ylabel('Absolute angle error [^o]')


disp('---------------------------------------')
disp(['Median error predicted bin = ' num2str(mean(median_error_angle_bin),2)])
disp(['Median error predicted angle = ' num2str(mean(median_error),2)])
disp('---------------------------------------')


subplot(2,3,5)
box off
hold on
xlabel('Recording')
ylabel('Average Accuracy per fold')
ylim([0 1])
plot([1 size(Sessions,2)],[1/8 1/8],'k')

text(1,0.4,['Mean M1 =' num2str(mean(acc_error_mean(M1)),3)],'FontSize',8)
text(10,0.4,['Mean PMd =' num2str(mean(acc_error_mean(~M1)),3)],'FontSize',8)



end

function [angle_dir,class_dir,neural_mov]=extract_movement_for_decoding(Session, Area,Ndir)
%% extract_movement_for_decoding extracts all movements and their respective
%% *preparatory* neural activity.
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
% OUTPUTS
%
% angle_dir: angle of each movment [rad]
%
% class_dir: number of class/ direction bin of each movement
%
% neural_mov: Preparatory activity (upto 200 ms before mov onset) of each
% movement
%
% 27/05/2023
% Andrea Colins Rodriguez

load(Session,Area,'trial_table2','cont')
startt=trial_table2(1,1);
endt=trial_table2(size(trial_table2,1),22);

if strcmp(Area,'PMd')
    neural_data=PMd.units;
else
    neural_data=M1.units;
end

ISI=compute_ISI(neural_data,startt,endt);
sigma_filter=round(median(ISI));

t_from=-0.2;
t_upto=-0.1;
from_to=[0.2 1.2];

[Neural_info,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto,from_to);
neural_mov=Neural_info.FR;
angle_dir=Mov_params.direction;
class_dir=ceil(Ndir*(angle_dir+pi)/(2*pi));

end
