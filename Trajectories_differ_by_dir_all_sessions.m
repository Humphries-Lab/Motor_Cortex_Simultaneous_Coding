function [t_from_new,t_upto_new]=Trajectories_differ_by_dir_all_sessions(Sessions,Areas,threshold,Ndir,Nbins)
%% Trajectories_differ_by_dir_all_sessions calculates the distance between
%% trajectories corresponding to different movement directions
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
% Nbins: number of durations to bin the movements
%
% OUTPUTS
%
% t_from_new= array containing the time at which the trajectories are
% closest to each other before movement onset. Each element corresponds to
% the value of each Session
%
% t_upto_new= array containing the time at which the trajectories are
% closest to each other after the movement ends. Each element corresponds to
% the value of each Session.
%
% 24/01/2023
% Andrea Colins Rodriguez

fig1=figure;
counter=1;
start=nan(Ndir*Nbins,1);
endt=nan(Ndir*Nbins,1);
max_sep=nan(Ndir*Nbins,1);
M1=nan(Ndir*Nbins,1);
t_from_new=nan(numel(Areas),1);
t_upto_new=nan(numel(Areas),Nbins);
ms=1000; %to convert from s to ms


figure(fig1)
for i_area=1:numel(Areas)
    Area=Areas{i_area};
    session=Sessions{i_area};
    
    load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto','edges_dur_bin')
    
    if i_area==1
        %plot_trajectories for duration bin =1
        plot_trajectories_prep_exec_after(score(idx_duration==1,:),idx_dir(idx_duration==1),t_from,edges_dur_bin(1:2))
    end
    
    if strcmp(Area,'M1')
        colourArea=[85 30 116]./256;
    else
        colourArea=[89 156 153]./256;
    end
    
    t1_tmp=nan(Nbins,1);
    t2_tmp=nan(Nbins,1);
    for i_bin=1:Nbins
        
        ndim=find(cumsum(variance)>threshold,1,'First');
        score2=score(idx_duration==i_bin,:);
        
        %% Compute distace between trajectories
        segment_length=round((t_upto(i_bin)-t_from)*ms);
        distances=zeros(Ndir,Ndir,segment_length);
        for i=1:Ndir
            for j=1:Ndir
                traj_i=score2(idx_dir(idx_duration==i_bin)==i,1:ndim);
                traj_j=score2(idx_dir(idx_duration==i_bin)==j,1:ndim);
                distances(i,j,:)=sqrt(sum((traj_i-traj_j).^2,2));
            end
        end
        max_distance=max(distances(:));
        
        angle_diff=linspace(0,360,Ndir+1)-180;
        angle_diff(end)=[];
        angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); %center at zero
        all_distances_time=zeros(Ndir,size(distances(i,:,:),3));
        all_distances=zeros(Ndir,Ndir);
        
        for i=1:Ndir
            
            all_distances(i,:)=circshift(mean(distances(i,:,:),3),Ndir/2-i)./max_distance;
            dirs=1:Ndir;
            dirs(i)=[];
            all_distances_time(i,:)=squeeze(mean(distances(i,dirs,:),2))./max_distance;
        end
        
        [~,idx_prep]=min(mean(all_distances_time(:,1:round(-t_from*ms))));
        [~,idx_max]=max(mean(all_distances_time(:,1:round(-t_from*ms))));
        [~,idx_mov]=min(mean(all_distances_time(:,round(-t_from*ms)+1:end)));
        
        binsize=edges_dur_bin(2)-edges_dur_bin(1);
        
        start(counter)=(idx_prep+t_from*ms)/ms;
        endt(counter)=(idx_mov-edges_dur_bin(i_bin)*ms-binsize*ms)/ms;
        max_sep(counter)=(idx_max+t_from*ms)/ms;
        t1_tmp(i_bin)=(idx_prep+t_from*ms)/ms;
        t2_tmp(i_bin)=(idx_mov-edges_dur_bin(i_bin)*ms-binsize)/ms;
        
        M1(counter)=strcmp(Areas{i_area},'M1');
        %% Plots!!
        % distance across directions
        subplot(2,6,5:6)
        plot(t_from*ms:round(t_upto(i_bin)*ms-1),mean(all_distances_time),'Color',colourArea)
        hold on
        
        subplot(2,6,10)
        plot(~M1(counter)+1,t1_tmp(i_bin) ,'.','Color',colourArea)
        hold on
        
        subplot(2,6,11)
        plot(~M1(counter)+1,max_sep(counter) ,'.','Color',colourArea)
        hold on
        
        subplot(2,6,12)
        plot(~M1(counter)+1,t2_tmp(i_bin) ,'.','Color',colourArea)
        hold on
        
        counter=counter+1;
        
        if  i_bin==2
            subplot(2,6,4)
            hold on
            errorbar(angle_diff,mean(all_distances),std(all_distances),'Color',colourArea)
            
            
        end
        clear all_distances all_distances_time
    end
    
    t_from_new(i_area)=round(mean(t1_tmp),3);
    t_upto_new(i_area,:)=max(round(mean(t2_tmp),3),0)+edges_dur_bin(1:Nbins)+binsize;
    
end
subplot(2,6,4)
hold on
xlabel('Angle difference [Degrees]')
ylabel('Normalised distance')
box off


subplot(2,6,5:6)
box off
xlabel('Time to movement onset [ms]')
ylabel('Average distace to other trajectories')


% minimum separation before movement onset
subplot(2,6,10)
hold on
errorbar(1,mean(start(M1==1)),std(start(M1==1)),'Color',[85 30 116]./256)
errorbar(2,mean(start(~M1==1)),std(start(M1~=1)),'Color',[89 156 153]./256)
box off
xlim([0.8  2.2])
xticks([1 2])
xticklabels({'M1','PMd'})

title('Minimum sep before movement onset')
ylabel('Time to movement onset [s]')
[~,p_start]=ttest2(start(M1==1),start(M1~=1));
text(0.9,-0.5,['p-value = ' num2str(p_start,1)],'FontSize',8)
text(0.9,-0.4,['Mean M1 = ' num2str(mean(start(M1==1)),2)],'FontSize',8)
text(1.4,-0.2,['Mean PMd = ' num2str(mean(start(~M1==1)),2)],'FontSize',8)

% maximum separation before movement onset
subplot(2,6,11)
errorbar(1,mean(max_sep(M1==1)),std(max_sep(M1==1)),'Color',[85 30 116]./256)
errorbar(2,mean(max_sep(~M1==1)),std(max_sep(M1~=1)),'Color',[89 156 153]./256)
box off
ylabel('Time to movement onset')
title('Time of maximum separation')
xlim([0.8  2.2])
xticks([1 2])
xticklabels({'M1','PMd'})
[~,p_max]=ttest2(max_sep(M1==1),max_sep(M1~=1));
text(0.9,-0.15,['p-value = ' num2str(p_max,1)],'FontSize',8)
text(0.9,-0.12,['Mean M1 = ' num2str(mean(max_sep(M1==1)),2)],'FontSize',8)
text(1.4,-0.02,['Mean PMd = ' num2str(mean(max_sep(M1~=1)),2)],'FontSize',8)


% minimum separation after movement end
subplot(2,6,12)
errorbar(1,mean(endt(M1==1)),std(endt(M1==1)),'Color',[85 30 116]./256)
errorbar(2,mean(endt(M1~=1)),std(endt(M1~=1)),'Color',[89 156 153]./256)
box off
xlim([0.8  2.2])
xticks([1 2])
xticklabels({'M1','PMd'})
title('Time of minimum separation after movement end')
ylabel('Time to movement end [s]')
[~,p_end]=ttest2(endt(M1==1),endt(M1~=1));
text(0.9,-0.31,['p-value = ' num2str(p_end,1)],'FontSize',8)
text(0.9,-0.2,['Mean M1 = ' num2str(mean(endt(M1==1)),2)],'FontSize',8)
text(1.4,0.2,['Mean PMd = ' num2str(mean(endt(M1~=1)),2)],'FontSize',8)

end