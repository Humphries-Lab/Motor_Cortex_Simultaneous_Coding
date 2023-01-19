function [t_1_new,t_2_new]=Trajectories_differ_by_dir_all_sessions(sessions,Areas,threshold,Ndir,Nbins)
fig1=figure;
counter=1;
start=nan(Ndir*Nbins,1);
endt=nan(Ndir*Nbins,1);
max_sep=nan(Ndir*Nbins,1);
M1=nan(Ndir*Nbins,1);
t_1_new=nan(numel(Areas),1);
t_2_new=nan(numel(Areas),Nbins);

figure(fig1)
for i_area=1:numel(Areas)
    Area=Areas{i_area};
    session=sessions{i_area};
    
    load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','from')
    
    if i_area==1
    plot_trajectories_prep_exec_after(score(idx_duration==1,:),idx_dir(idx_duration==1),t_1,from)
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
        segment_length=round((t_2(i_bin)-t_1)*1000);
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

        [~,idx_prep]=min(mean(all_distances_time(:,1:round(-t_1*1000))));
        [~,idx_max]=max(mean(all_distances_time(:,1:round(-t_1*1000))));
        [~,idx_mov]=min(mean(all_distances_time(:,round(-t_1*1000)+1:end)));
        
        start(counter)=(idx_prep+t_1*1000)/1000;
        endt(counter)=(idx_mov-from(i_bin)*1000-100)/1000;
        max_sep(counter)=(idx_max+t_1*1000)/1000;  
        t1_tmp(i_bin)=(idx_prep+t_1*1000)/1000;
        t2_tmp(i_bin)=(idx_mov-from(i_bin)*1000-100)/1000;
        
        M1(counter)=strcmp(Areas{i_area},'M1');
        %% Plots!!
        % distance across directions
        subplot(2,6,5:6)
        plot(t_1*1000:round(t_2(i_bin)*1000-1),mean(all_distances_time),'Color',colourArea)
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
        
        if  i_bin==1
            subplot(2,6,4)
            errorbar(angle_diff,mean(all_distances),std(all_distances),'Color',colourArea)
            xlabel('Angle difference [Degrees]')
            ylabel('Normalised distance')
            hold on
            box off
            clear all_distances all_distances_time
        end
    end
    
   t_1_new(i_area)=round(mean(t1_tmp),3);
   t_2_new(i_area,:)=max(round(mean(t2_tmp),3),0)+from+0.1;
   
end

    subplot(2,6,5:6)
        box off
        xlabel('Time to movement onset [ms]')
        ylabel('Average distace to other trajectories')
        
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
    text(0.9,-0.14,['p-value = ' num2str(p_max,1)],'FontSize',8)
    
    subplot(2,6,10)
    hold on
    errorbar(1,mean(start(M1==1)),std(start(M1==1)),'Color',[85 30 116]./256)
    errorbar(2,mean(start(~M1==1)),std(start(M1~=1)),'Color',[89 156 153]./256)
    box off
    xlim([0.8  2.2])
    xticks([1 2])
    xticklabels({'M1','PMd'})
    
    title('Time of minimum separation before movement onset')
    ylabel('Time to movement onset [s]')
    %% ttest
    [~,p_start]=ttest2(start(M1==1),start(M1~=1));
    text(0.9,-0.5,['p-value = ' num2str(p_start,1)],'FontSize',8)
    
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
end