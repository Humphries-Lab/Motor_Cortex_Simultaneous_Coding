function [variance,mov_distance,mov_duration,max_speed]=embedding_dimensions(Session,Area,Ndir,Nbins,t_from,t_upto,edges_dur_bin,do_plot)
%% embedding_dimensions performs PCA on the neural data from Session and
%% Area binning the movements into Ndir directions and Nbins durations
% This function also saves the PCA results on a file for later analyses
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
% Nbins: number of durations to bin the movements
%
% t_from: start time of the neural activity relative to movement onset [S]
% e.g t_from=-0.5
%
% t_upto: end time of the neural activity relative to movement end [S]
% e.g t_from=0.3
%
% edges_dur_bin= array containing the edges of each duration bin [S]
%
% do_plot= 1- plots the average FR of the population for each duration bin
%             and the PCA trajectories for the first duration bin
%
%          0- omit the above plot
% OUTPUTS
%
% variance: array containing the variance explained by each PC
%
% mov_distance= cell array containing the movement distance of all
% movements selected for PCA. Each cell contains an array with the distance of each movement within a duration bin
%
% mov_duration= cell array containing the movement duration of all
% movements selected for PCA. Each cell contains an array with the duration of each movement within a duration bin
%
% max_speed= cell array containing the maximum speed reached of all
% movements selected for PCA. Each cell contains an array with the maximum speed of each movement within a duration bin
%
% 24/05/2023
% Andrea Colins Rodriguez


if do_plot
    figure
    colour_dir=hsv(Ndir);
    colour_plasma=plasma(Nbins);
end
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

nsamples_condition=zeros(Ndir,Nbins);
mov_distance=cell(Nbins,1);
mov_duration=cell(Nbins,1);
max_speed=cell(Nbins,1);
dur_binsize=(edges_dur_bin(2)-edges_dur_bin(1));
t_upto=edges_dur_bin(1:Nbins)+dur_binsize+t_upto;
total_time_bins=sum(round((t_upto-t_from)*1000))*Ndir;

average_cond=nan(size(neural_data,2),total_time_bins);
idx_dir=nan(total_time_bins,1);
idx_duration=nan(total_time_bins,1);
counter=0;
% Extract movements and related neural activity for eac duration bin
for i_dur=1:Nbins
    
    current_dur_bin=[edges_dur_bin(i_dur) edges_dur_bin(i_dur+1)]; % select movements in this range of duration only
    
    [Neural_info,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);
    
    %Save Mov params for later analyses
    mov_distance{i_dur}=Mov_params.distance;
    mov_duration{i_dur}=Mov_params.duration;
    max_speed{i_dur}=Mov_params.max_speed;
    
    
    % bin by direction
    direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));
    
    if do_plot
        xtime=round(t_from*1000):1:round(t_upto(i_dur)*1000-1);
        FR=mean(mean(Neural_info.FR,3));
        
        subplot(4,5,1:3)
        plot(xtime,FR,'Color',colour_plasma(i_dur,:))
        hold on
        ylabel('Firing rate')
        title([ Area ' ' Session ])
        
        
    end
    
    for i_dir=1:Ndir
        nsamples_condition(i_dir,i_dur)=sum(direction1==i_dir);
        if nsamples_condition(i_dir,i_dur)>=1
            timebins=round((t_upto(i_dur)-t_from)*1000);
            average_cond(:,counter+1:counter+timebins)=mean(Neural_info.FR(:,:,direction1==i_dir),3);
            idx_dir(counter+1:counter+timebins)=zeros(timebins,1)+i_dir;
            idx_duration(counter+1:counter+timebins)=zeros(timebins,1)+i_dur;

            counter=counter+timebins;
        else
            disp([Session ' ' Area ' direction = ' num2str(i_dir) ' has less than 2 samples'])
            keyboard
        end
        
    end
    
end


%% Compute a common subspace
% soft normalization
%delete neurons that have minimal number of spikes
delete_units=sum(average_cond,2)<5;
average_cond(delete_units,:)=[];
normalisation=range(average_cond')+5;
average_cond=average_cond'./repmat(normalisation,size(average_cond,2),1);
mean_norm=mean(average_cond);
%% Project same directions coloring different times
[coeffs,score,~,~,variance]=pca(average_cond);

if do_plot
    
    for i_dir=1:Ndir
        for i_bin=1
            
            idx=find(idx_dir==i_dir & idx_duration==i_bin);
            xtime=round(t_from*1000):1:round(t_upto(i_bin)*1000-1);
            
            subplot(4,5,[4 5 9 10])
            plot3(score(idx,1),score(idx,2),score(idx,3),'Color',colour_dir(i_dir,:),'LineWidth',2)
            hold on
            plot3(score(idx(1),1),score(idx(1),2),score(idx(1),3),'o','MarkerFaceColor',colour_dir(i_dir,:),'MarkerEdgeColor',colour_dir(i_dir,:))
            subplot(4,5,6:8)
            plot(xtime,score(idx,1),'Color',colour_dir(i_dir,:),'LineWidth',2)
            hold on
            subplot(4,5,11:13)
            plot(xtime,score(idx,2),'Color',colour_dir(i_dir,:),'LineWidth',2)
            hold on
            subplot(4,5,16:18)
            plot(xtime,score(idx,3),'Color',colour_dir(i_dir,:),'LineWidth',2)
            hold on
            
        end
    end
    subplot(4,5,[4 5 9 10])
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3')
    subplot(4,5,16:18)
    xlabel('Time to movement onset')
    ylabel('PC 3')
    subplot(4,5,11:13)
    ylabel('PC 2')
    subplot(4,5,6:8)
    ylabel('PC 1')
    
    
end
save(['../Output_files/PCA_' Session(1:end-4) '_' Area '.mat'],'coeffs','score','idx_dir','idx_duration','variance','edges_dur_bin','sigma_filter','nsamples_condition','delete_units','normalisation','mean_norm','t_from','t_upto')
end