function [t_1_new,t_2_new]=Trajectories_differ_by_dir_all_sessions(sessions,Areas,threshold,Ndir,Nbins)
fig1=figure;
colour_dir=hsv(Ndir);
%colour_Area=jet(numel(Areas));
counter=1;
start=nan(Ndir*Nbins,1);
endt=nan(Ndir*Nbins,1);
max_sep=nan(Ndir*Nbins,1);
M1=nan(Ndir*Nbins,1);
counter2=0;
for i_area=1:numel(Areas)
    Area=Areas{i_area};
    session=sessions{i_area};
    %load(['scores_LDS_diff_duration_ms_original_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance')
    load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','from')
    
    % plot_trajectories_prep_exec_after(score(idx_duration==1,:),idx_dir(idx_duration==1),t_1,from)
    
    if strcmp(Area,'M1')
        colourArea=[85 30 116]./256;
    else
        colourArea=[89 156 153]./256;
    end
    
    %figure
    for i_bin=1:Nbins
        t1_tmp=nan(Nbins,1);
        t2_tmp=nan(Nbins,1);
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
        %% Plots!!
        %% find time when are maximally separated
        %         [~,idx_max]=max(distances,[],3);
        %         idx_max=idx_max+diag(ones(Ndir,1)*nan);
        % imagesc(idx_max>300)
        % pause
        % distance across directions
        angle_diff=linspace(0,360,Ndir+1)-180;
        angle_diff(end)=[];
        angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); %center at zero
        all_distances_time=zeros(Ndir,size(distances(i,:,:),3));
        all_distances=zeros(Ndir,Ndir);
        for i=1:Ndir
            
            all_distances(i,:)=circshift(mean(distances(i,:,:),3),Ndir/2-i)./max_distance;
            subplot(3,2,1)
            plot(angle_diff,circshift(mean(distances(i,:,:),3),Ndir/2-i)./max_distance,'Color',colour_dir(i,:))
            hold on
            dirs=1:Ndir;
            dirs(i)=[];
            all_distances_time(i,:)=squeeze(mean(distances(i,dirs,:),2))./max_distance;
        end
        subplot(3,2,2)
        plot(t_1*1000:round(t_2(i_bin)*1000-1),mean(all_distances_time),'Color',colourArea)
        hold on
        counter2=counter2+1;
        box off
        xlabel('Time to movement onset [ms]')
        ylabel('Average distace to other trajectories')
        
        [~,idx_prep]=min(mean(all_distances_time(:,1:round(-t_1*1000))));
        [~,idx_max]=max(mean(all_distances_time(:,1:round(-t_1*1000))));
        [~,idx_mov]=min(mean(all_distances_time(:,round(-t_1*1000)+1:end)));
        start(counter)=(idx_prep+t_1*1000)/1000;
        endt(counter)=(idx_mov-from(i_bin)*1000-100)/1000;
        max_sep(counter)=(idx_max+t_1*1000)/1000;
        M1(counter)=strcmp(Areas{i_area},'M1');
        t1_tmp(i_bin)=(idx_prep+t_1*1000)/1000;
        t2_tmp(i_bin)=(idx_mov-from(i_bin)*1000-100)/1000;
        
        subplot(3,2,1)
        box off
        xlabel('Angle difference [\circ]')
        ylabel('Mean euclidean distance')
               
        subplot(3,2,4)
        plot(~M1(counter)+1,max_sep(counter) ,'.','Color',colourArea)
        hold on
        
        
        subplot(3,2,5)
        plot(~M1(counter)+1,t1_tmp(i_bin) ,'.','Color',colourArea)
        hold on
        
        subplot(3,2,6)
        plot(~M1(counter)+1,t2_tmp(i_bin) ,'.','Color',colourArea)
        hold on
        
        counter=counter+1;
        if  i_bin==1
            subplot(3,2,3)
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
    counter2
    
    subplot(3,2,4)
    errorbar(1,mean(max_sep(M1==1)),std(max_sep(M1==1)),'Color',[85 30 116]./256)
    errorbar(2,mean(max_sep(~M1==1)),std(max_sep(M1~=1)),'Color',[89 156 153]./256)
    box off
    ylabel('Time to movement onset') 
    title('Time of maximum separation')
    xlim([0.8  2.2])
    [h,p_max]=ttest2(max_sep(M1==1),max_sep(M1~=1))
    
    subplot(3,2,5)
    hold on
    errorbar(1,mean(start(M1==1)),std(start(M1==1)),'Color',[85 30 116]./256)
    errorbar(2,mean(start(~M1==1)),std(start(M1~=1)),'Color',[89 156 153]./256)
    box off
    xlim([0.8  2.2])
    title('Time of minimum separation before movement onset')
    ylabel('Time to movement onset')
    
    %% ttest
    [h,p_start]=ttest2(start(M1==1),start(M1~=1))
    
    subplot(3,2,6)
    errorbar(1,mean(endt(M1==1)),std(endt(M1==1)),'Color',[85 30 116]./256)
    errorbar(2,mean(endt(M1~=1)),std(endt(M1~=1)),'Color',[89 156 153]./256)  
    box off
    xlim([0.8  2.2])
    title('Time of minimum separation after movement end')
    ylabel('Time to movement end')
    [h,p_end]=ttest2(endt(M1==1),endt(M1~=1))
end