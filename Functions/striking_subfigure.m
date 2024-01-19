function [variance,mov_distance,mov_duration,max_speed]=striking_subfigure(Session,Area,Ndir,Nbins,t_from,t_upto,edges_dur_bin,do_plot)
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

ms=1000;

if do_plot
    figure
    colour_dir=My_colour_palette;
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
%sigma_filter=50;

nsamples_condition=zeros(Ndir,Nbins);
mov_distance=cell(Nbins,1);
mov_duration=cell(Nbins,1);
max_speed=cell(Nbins,1);
dur_binsize=(edges_dur_bin(2)-edges_dur_bin(1));
t_upto=edges_dur_bin(1:Nbins)+dur_binsize+t_upto;
total_time_bins=sum(round((t_upto(1:Nbins)-t_from)*ms))*Ndir;

average_cond=nan(size(neural_data,2),total_time_bins);
idx_dir=nan(total_time_bins,1);
idx_duration=nan(total_time_bins,1);
counter=0;
Neural_all=cell(Nbins,1);
% Extract movements and related neural activity for eac duration bin
for i_dur=1:Nbins
    
    current_dur_bin=[edges_dur_bin(i_dur) edges_dur_bin(i_dur+1)] % select movements in this range of duration only
    
    [Neural_info,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);
    %[Neural_info,Mov_params]=neural_data_per_duration_normalised(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);
    [Mov_params.duration,iddir]=sort(Mov_params.duration);
    FR=squeeze(mean(Neural_info.FR(:,:,iddir),1));
    FR=FR(:,abs(Mov_params.direction(iddir))<0.3927);
    current_dur=round(ms*Mov_params.duration(abs(Mov_params.direction(iddir))<0.3927)-t_from*ms);
    
    
    figure
    set(gca,'Color',[0 0 0])
   hold on
   Ntrials=size(FR,2);
   colour_dur=plasma(Ntrials);
    for trial=1:4:Ntrials-20
        plot3(zeros(1,current_dur(trial))+Ntrials-trial,1:current_dur(trial),mean(FR(1:current_dur(trial),trial:trial+20),2),'Color',colour_dur(trial,:))
    end

    direction1=ceil(Ndir*(Mov_params.direction(iddir)+pi)/(2*pi));
    duration=Mov_params.duration(iddir);
    FR2=squeeze(mean(Neural_info.FR(:,:,iddir),1));
    colour_dur=flipud(plasma(Ntrials*1.1));
    
    %% circle
    figure
    %set(gca,'Color',[0 0 0])
    
    for i_dir=1:Ndir
    FRtmp=FR2(:,direction1==i_dir);
    current_dur=round(ms*duration(direction1==i_dir)-t_from*ms);
    Ntrials=size(FRtmp,2)
    
    for trial=1:4:Ntrials-20
    polarplot(linspace(pi+(i_dir-1)*2*pi/Ndir,5*pi/4+(i_dir-1)*2*pi/Ndir,current_dur(trial)),mean(FRtmp(1:current_dur(trial),trial:trial+20),2)+trial*0.0002,'Color',colour_dur(trial,:))
    hold on
    
    end
    
    end
    rlim([0.006 0.07])
    axis off
    ax=gca;
    ax.Color=[0 0 0];
    % s=surface(1:size(FR,2),1:size(FR,1),FR);
    % s.EdgeColor = 'none';
   %imagesc(squeeze(mean(Neural_info.FR(:,:,iddir),1)))
    %Save Mov params for later analyses
    mov_distance{i_dur}=Mov_params.distance;
    mov_duration{i_dur}=Mov_params.duration;
    max_speed{i_dur}=Mov_params.max_speed;
    
    Neural_all{i_dur}=Neural_info;
    % bin by direction
    direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));
    
    
    
    
    for i_dir=1:Ndir
        nsamples_condition(i_dir,i_dur)=sum(direction1==i_dir);
        %close_dir{i_dur,i_dir}=direction1==i_dir;
        if nsamples_condition(i_dir,i_dur)>=2
            
            
            timebins=round((t_upto(i_dur)-t_from)*ms);
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
        for i_bin=1:Nbins
            
            idx=find(idx_dir==i_dir & idx_duration==i_bin);
            xtime=round(t_from*ms):1:round(t_upto(i_bin)*ms-1);
            
            subplot(4,5,[4 5 9 10])
            plot3(score(idx,1),score(idx,2),score(idx,3),'Color',colour_dir(i_dir,:)./sqrt(i_bin),'LineWidth',2)
            hold on
            plot3(score(idx(1),1),score(idx(1),2),score(idx(1),3),'o','MarkerFaceColor',colour_dir(i_dir,:)./sqrt(i_bin),'MarkerEdgeColor',colour_dir(i_dir,:)./sqrt(i_bin))
           
            subplot(4,5,6:8)

            %[size(xtime) size(idx)]
            plot(xtime,score(idx,1),'Color',colour_dir(i_dir,:)./sqrt(i_bin),'LineWidth',2)
            hold on
            subplot(4,5,11:13)
            plot(xtime,score(idx,2),'Color',colour_dir(i_dir,:)./sqrt(i_bin),'LineWidth',2)
            hold on
            subplot(4,5,16:18)
            plot(xtime,score(idx,3),'Color',colour_dir(i_dir,:)./sqrt(i_bin),'LineWidth',2)
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


end