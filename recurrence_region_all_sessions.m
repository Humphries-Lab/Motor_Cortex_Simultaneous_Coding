function recurrence_region_all_sessions(sessions,Areas,threshold,Ndir,do_plot)

fig1=figure;
fig2=figure;
colormap(flipud(colormap('gray')))
upto=200;% number of points to analyse from the beginning of prep
not_coming_back=zeros(numel(Areas),1);
not_starting_there=zeros(numel(Areas),1);
Nbins=4;
colour_dur=plasma(Nbins);
for j=1:numel(Areas)
    Area=Areas{j};
    session=sessions{j};
    for i_bin=[2:Nbins 1]
        
        
        %load(['scores_LDS_diff_duration_ms_original_longer_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance')
        load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','from')
        ndim=find(cumsum(variance)>threshold,1,'First');
        if i_bin==2
            %region=define_recurrence_region_v2(score(:,1:ndim),idx_dir(idx_duration==i_bin),upto);
            [limit,center]=define_recurrence_region_v3(score(idx_duration==i_bin,1:ndim),idx_dir(idx_duration==i_bin));
        end
        if do_plot
            figure(j+1) %per bin
        end
        %% compute recurrence time
        segment_length=round((t_2(i_bin)-t_1)*1000);
        Its=zeros(segment_length,Ndir);
        for i_dir=1:Ndir
            %Its(:,i_dir)=IsInRecRegion(score(idx_dir(idx_duration==i_bin)==i_dir,1:ndim),region);
            Its(:,i_dir)=IsInRecRegion2(score(idx_dir==i_dir & idx_duration==i_bin,1:ndim),[limit,center]);
            Start_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Its(1:upto,i_dir),1)>0;
            End_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Its(end-upto:end,i_dir),1)>0;
        end
        if do_plot
            subplot(4,1,i_bin)
            imagesc(Its')
            caxis([0 1])
        end
        
        figure(fig2)
        if strcmp(Areas{j},'M1')
            subplot(2,2,1)
            plot((1:size(Its,1))+t_1*1000,sum(Its,2),'Color',colour_dur(i_bin,:))
            hold on
        else
            subplot(2,2,3)
            plot((1:size(Its,1))+t_1*1000,sum(Its,2),'Color',colour_dur(i_bin,:))
            hold on
        end
        if  j==numel(Areas)
            xlabel('Time to movement onset [ms]')
            ylabel('Number of trajectories in the recurrence region')
            box off
        end
        clear Its
    end
    %% Compute number of points at the beginning and end. Or just if it gets there? both?
    figure(fig1)
    subplot(1,numel(Areas),j)
    imagesc([Start_recurrence;End_recurrence]')
    caxis([0 1])
    
    not_coming_back(j,1)=sum(End_recurrence<1)/numel(End_recurrence);
    not_starting_there(j,1)=sum(Start_recurrence<1)/numel(End_recurrence);
    
    clear Start_recurrence End_recurrence
end

M1=strcmp(Areas,'M1');
figure(fig2)
subplot(2,2,1)
title('M1')
box off
subplot(2,2,3)
title('PMd')

subplot(2,2,2)
plot([not_coming_back(M1) 1-not_coming_back(M1)]','m')
hold on
plot([not_coming_back(~M1) 1-not_coming_back(~M1)]','b')
xlim([0.8 2.2])
xticks([1 2])
xticklabels({'Do not recur', 'Recur'})
ylim([0 1])
title('End')
box off

subplot(2,2,4)
plot([not_starting_there(M1) 1-not_starting_there(M1)]','m')
hold on
plot([not_starting_there(~M1) 1-not_starting_there(~M1)]','b')
[mean(1-not_starting_there(M1)) mean(1-not_starting_there(~M1))]
[mean(1-not_coming_back(M1)) mean(1-not_coming_back(~M1))]
xlim([0.8 2.2])
xticks([1 2])
xticklabels({'Do not start at recurrence region','Start at recurrence region'})
ylim([0 1])
title('Start')
box off
end


function region=define_recurrence_region(score,idx_dir,rec_length)
debugging=0;
%rec_length=50;
Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
all_points_beginning=zeros(rec_length*Ndir,size(score,2));
for i_dir=1:Ndir
    all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,:)=score(find(idx_dir==i_dir,rec_length,'First'),:);
    if debugging
        subplot(2,2,1)
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        %pause
    end
end

% define region
lower_edge=min(all_points_beginning,[],1);
upper_edge=max(all_points_beginning,[],1);
region=[lower_edge;upper_edge];
if debugging
    plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) lower_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    Its=zeros(sum(idx_dir==i_dir),Ndir);
    for i_dir=1:Ndir
        Its(:,i_dir)=IsInRecRegion(score(idx_dir==i_dir,:),region);
    end
    
    subplot(2,2,2)
    imagesc(Its')
    subplot(2,2,3)
    for i_dir=1:Ndir
        
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) lower_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        subplot(2,2,4)
        imagesc(Its(:,i_dir)')
        pause
        subplot(2,2,3)
        cla
    end
    
    
    
end


end

function region=define_recurrence_region_v2(score,idx_dir)
%this version uses the threshold to define a region
debugging=0;
rec_length=1;
Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
all_points_beginning=zeros(rec_length*Ndir,size(score,2));
threshold=5;
limit=prctile(pdist(score),threshold);
for i_dir=1:Ndir
    all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,:)=score(find(idx_dir==i_dir,rec_length,'First'),:);
    if debugging
        subplot(2,2,1)
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        %pause
    end
end

center=mean(all_points_beginning);

% define region
lower_edge=center-limit;
upper_edge=center+limit;
region=[lower_edge;upper_edge];
if debugging
    plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
    
    plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([lower_edge(1) lower_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    plot3([upper_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
    
    for i_dir=1:Ndir
        Its(:,i_dir)=IsInRecRegion(score(idx_dir==i_dir,:),region);
    end
    
    subplot(2,2,2)
    imagesc(Its')
    subplot(2,2,3)
    for i_dir=1:Ndir
        
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[lower_edge(3) lower_edge(3)],'k')
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) upper_edge(2)],[upper_edge(3) upper_edge(3)],'k')
        
        plot3([lower_edge(1) lower_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([lower_edge(1) lower_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[lower_edge(2) lower_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        plot3([upper_edge(1) upper_edge(1)],[upper_edge(2) upper_edge(2)],[lower_edge(3) upper_edge(3)],'k')
        subplot(2,2,4)
        imagesc(Its(:,i_dir)')
        pause
        subplot(2,2,3)
        cla
    end
    
    
    
end


end

function [limit,center]=define_recurrence_region_v3(score,idx_dir)
%this version uses the threshold to define a region
debugging=0;
rec_length=1;
Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
all_points_beginning=zeros(rec_length*Ndir,size(score,2));
threshold=10;
limit=prctile(pdist(score),threshold);
plot_axes=[2 3 4];
for i_dir=1:Ndir
    %all_points_beginning((i_dir-1)*2*rec_length+1:i_dir*rec_length*2,:)=score([find(idx_dir==i_dir,rec_length,'First');find(idx_dir==i_dir,rec_length,'Last')],:);
    
    all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,:)=score(find(idx_dir==i_dir,rec_length,'First'),:);
    if debugging
        subplot(2,2,1)
        plot3(score(idx_dir==i_dir,plot_axes(1)),score(idx_dir==i_dir,plot_axes(2)),score(idx_dir==i_dir,plot_axes(3)),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,plot_axes(1)),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,plot_axes(2)),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,plot_axes(3)),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        %pause
    end
end

center=mean(all_points_beginning);
if debugging
    subplot(2,2,1)
    colormap(colormap('white')*0.5)
    [X,Y,Z] = sphere;
    mesh(X.*limit+center(plot_axes(1)),Y.*limit+center(plot_axes(2)),Z.*limit+center(plot_axes(3)),'FaceAlpha','0.5')
    
    hold on
    min_all=min(min(score(:,plot_axes)));
    max_all=max(max(score(:,plot_axes)));
    xlim([min_all max_all])
    ylim([min_all max_all])
    zlim([min_all max_all])
end

% region is a sphere radius=threshold around the center
%compare distances of all points to the center

if debugging
    
    for i_dir=1:Ndir
        Its(:,i_dir)=pdist2(center,score(idx_dir==i_dir,:))<limit;
        %Its(:,i_dir)=IsInRecRegion(score(idx_dir==i_dir,:),region);
    end
    
    subplot(2,2,2)
    imagesc(Its')
    
    subplot(2,2,3)
    for i_dir=1:Ndir
        
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        
        subplot(2,2,4)
        imagesc(Its(:,i_dir)')
        pause
        subplot(2,2,3)
        cla
    end
    
end


end

function Its=IsInRecRegion(score,region)
Ndim=size(score,2);

%Its=zeros(size(score,1),1);
all_dim=zeros(size(score));
for i=1:Ndim
    all_dim(:,i)=score(:,i)>=region(1,i) & score(:,i)<=region(2,i);
    
end
Its=sum(all_dim,2)==Ndim;
end

function Its=IsInRecRegion2(score,region)
limit=region(1);
center=region(2:end);

Its=pdist2(center,score)<limit;
end