function embedding_dimensions_all_sessions(session,Area,threshold,Ndir,Nbins,do_plot,t_1,t_2,from,session_N)

fig1=figure;
fig1sup=figure;
embedding_dim=zeros(size(session));
R(max(session_N))=struct();
for i=1:size(session,2)
    disp(['Recording ' num2str(i)])
    %      [~,t_1_new,t_2_new]=embedding_dimensions(session{i}, Area{i},0,Ndir,t_1(i),t_2(i,:),from);
    %      t_1(i)=t_1_new;
    %      t_2(i,:)=t_2_new;
    
    [variance,~,~,dist_mov_dir,mov_duration,max_speed]=embedding_dimensions(session{i}, Area{i},do_plot,Ndir,Nbins,t_1(i),t_2(i,:),from);
    %% kinematics
    if i==1
        figure(fig1sup)
        plot_kinematics(dist_mov_dir,mov_duration,max_speed,Nbins)
    end
    % Here the sessions are repeated because of the shared regions
    % Andrea needs to fix this
    R(session_N(i)).R_dist_dur=corr([dist_mov_dir{1:Nbins}]',[mov_duration{1:Nbins}]');
    R(session_N(i)).R_dist_speed=corr([dist_mov_dir{1:Nbins}]',[max_speed{1:Nbins}]');
    R(session_N(i)).R_speed_dur=corr([max_speed{1:Nbins}]',[mov_duration{1:Nbins}]');
    
    % Dimensionality of the subspace
    embedding_dim(i)=find(cumsum(variance)>threshold,1,'First');
    
    if strcmp(Area{i},'M1')
        colourArea=[85 30 116]./256;
    else
        colourArea=[89 156 153]./256;
    end
    
    figure(fig1)
    subplot(2,1,1)
    plot(cumsum(variance),'Color',colourArea)
    hold on
    
    subplot(2,1,2)
    hold on
    plot(i,embedding_dim(i),'o','Color',colourArea)
    
    clear variance dist_mov_dir mov_duration max_speed
    
end

figure(fig1)
subplot(2,1,1)
plot([0 100],[threshold threshold])
xlabel('N units')
ylabel('Variance explained [%]')
box off

subplot(2,1,2)
hold on
box off
xlabel('N session')
ylabel('Embedding dimensions')
ylim([0 13])
%%
figure(fig1sup)

subplot(2,3,4)
histogram([R.R_dist_dur])
box off
xlabel('Corr Distance-Duration')
ylabel('Number of sessions')
xlim([0 1])
title(['Mean = ' num2str(mean([R.R_dist_dur]),'%.2f')])

subplot(2,3,5)
histogram([R.R_dist_speed])
box off
xlabel('Corr Distance-Speed')
xlim([0 1])
title(['Mean = ' num2str(mean([R.R_dist_speed]),'%.2f')])

subplot(2,3,6)
histogram([R.R_speed_dur])
box off
xlabel('Corr Distance-Speed')
xlim([0 1])
title(['Mean = ' num2str(mean([R.R_speed_dur]),'%.2f')])
end


function [variance,t_1_new,t_2_new,dist_mov_dir,mov_duration,max_speed]=embedding_dimensions(session, Area,do_plot,Ndir,Nbins,t_1,t_2,from)
%% Average across bins of segments durations (200,300,400,500 ms)
event=2;
if do_plot
    figure
    colour_dir=hsv(Ndir);
end

if strcmp(Area,'PMd')
    load(session,'PMd','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(PMd.units,startt,endt);
    %t_1=-0.5;
else
    load(session,'M1','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(M1.units,startt,endt);
    %t_1=-0.5;
end

ms=nanmedian(ISI); %round(ms*sqrt(12))

% from=[0.2 0.3 0.4 0.5];
%
% t_2=from+0.1+0.5;%movement duration+ 200 ms
average_cond_1=[];
idx_dir=[];
idx_duration=[];
colour_plasma=plasma(Nbins);
counter=1;
nsamples_condition=zeros(Ndir,Nbins);
for i=1:Nbins
    from_to=[from(i) from(i)+0.1];
    [condition_matrix{i},direction,~,~,~,reach_number{i},dist_mov_dir{i},mov_duration{i},max_speed{i}]=neural_data_per_duration(session,Area,ms,t_1,t_2(i),event,from_to,0);
    xtime=round(t_1*1000):1:round(t_2(i)*1000-1);
    %take average by direction
    direction1=ceil(Ndir*(direction+pi)/(2*pi));
    if do_plot
        subplot(4,5,1:3)
        FR=mean(mean(condition_matrix{i},3));
        plot(xtime,FR,'Color',colour_plasma(i,:))
        hold on
        ylabel('Firing rate')
        title([ Area ' ' session ])
        
        
    end
    
    for cond=1:Ndir
        nsamples_condition(cond,i)=sum(direction1==cond);
        if nsamples_condition(cond,i)>=2
            average_cond_1=[average_cond_1,mean(condition_matrix{i}(:,:,direction1==cond),3)];
            idx_dir=[idx_dir;zeros(round((t_2(i)-t_1)*1000),1)+cond];
            idx_duration=[idx_duration;zeros(round((t_2(i)-t_1)*1000),1)+i];
            counter=counter+1;
            %         if do_plot
            %             subplot(4,5,20)
            %             polarplot(direction(direction1==cond),1,'.','Color',colour_dir(cond,:))
            %             hold on
            %         end
            
        else
            cond
            [session, Area]
            pause
        end
        
    end
    
end


%% Compute a common subspace
% soft normalization
%delete neurons that have minimal number of spikes
delete_units=sum(average_cond_1,2)<5;
average_cond_1(delete_units,:)=[];
normalisation=range(average_cond_1')+5;
average_cond_1=average_cond_1'./repmat(normalisation,size(average_cond_1,2),1);
mean_norm=mean(average_cond_1);
%% Project same directions coloring different times
[coeffs,score,~,~,variance]=pca(average_cond_1);
threshold=80;
embedding_dim=find(cumsum(variance)>threshold,1,'First');
for i_bin=1:Nbins
    idx=find(idx_duration==i_bin);
    occupancy=subspace_occupancy(score(idx,1:embedding_dim),idx_dir(idx));
    [~,idx_prep]=min(occupancy(1:round(-t_1*1000)));
    t1_tmp(i_bin)=(idx_prep+t_1*1000)/1000;
    [~,idx_mov]=min(occupancy(round(-t_1*1000)+1:end));
    t2_tmp(i_bin)=(idx_mov-from(i_bin)*1000-100)/1000;
    
    
    if do_plot
        xtime=round(t_1*1000):1:round(t_2(i_bin)*1000-1);
        subplot(4,5,20)
        hold on
        plot(xtime,occupancy,'Color',colour_plasma(i_bin,:))
    end
    
end
t_1_new=round(mean(t1_tmp),3);
t_2_new=max(round(mean(t2_tmp),3),0)+from+0.1;

if do_plot
    
    for i_dir=1:Ndir
        for i_bin=1
            
            idx=find(idx_dir==i_dir & idx_duration==i_bin);
            xtime=round(t_1*1000):1:round(t_2(i_bin)*1000-1);
            
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
save(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'coeffs','score','idx_dir','idx_duration','variance','t_1','t_2','from','ms','nsamples_condition','delete_units','normalisation','mean_norm')
end

function occupancy=subspace_occupancy(score,idx_dir)
Ndir=max(idx_dir);
%reshape data to compute variance across directions for each point in time
score2=zeros(sum(idx_dir==1),size(score,2),Ndir);
for i=1:Ndir
    score2(:,:,i)=score(idx_dir==i,:);
end
occupancy=sum(var(score2,[],3),2);
end
