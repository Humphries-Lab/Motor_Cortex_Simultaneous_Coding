function recurrence_region_all_sessions(sessions,Areas,threshold,Ndir)

fig2=figure;
colormap(flipud(colormap('gray')))
upto=200;% number of points to analyse from the beginning of prep
not_coming_back=zeros(numel(Areas),1);
not_starting_there=zeros(numel(Areas),1);
Nbins=4;
colour_dur=plasma(Nbins);

for isession=1:numel(Areas)
    Area=Areas{isession};
    session=sessions{isession};
    End_recurrence=true(Nbins*Ndir,1);
    Start_recurrence=true(Nbins*Ndir,1);
    for i_bin=[2:Nbins 1]
        
        
        load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto')

       % load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2')
        ndim=find(cumsum(variance)>threshold,1,'First');
        if i_bin==2
             if isession==1
                 do_example_plot=1;
             else
                 do_example_plot=0;
             end
            [limit,center]=define_recurrence_region_v3(score(idx_duration==i_bin,1:ndim),idx_dir(idx_duration==i_bin),do_example_plot);

        end
        
        %% compute recurrence time
        segment_length=round((t_upto(i_bin)-t_from)*1000);
        Its=zeros(segment_length,Ndir);
        for i_dir=1:Ndir
            %Its(:,i_dir)=IsInRecRegion(score(idx_dir(idx_duration==i_bin)==i_dir,1:ndim),region);
            Its(:,i_dir)=IsInRecRegion2(score(idx_dir==i_dir & idx_duration==i_bin,1:ndim),[limit,center]);
            Start_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Its(1:upto,i_dir),1)>0;
            End_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Its(end-upto:end,i_dir),1)>0;
        end
        
        figure(fig2)
        if strcmp(Areas{isession},'M1')
            subplot(2,2,2)
            hold on
        else
            subplot(2,2,4)
            hold on
        end
        plot((1:size(Its,1))+t_from*1000,sum(Its,2),'Color',colour_dur(i_bin,:))
        
        
        clear Its
    end
    %% Compute number of points at the beginning and end. Or just if it gets there? both?
    
    not_coming_back(isession,1)=sum(End_recurrence<1)/numel(End_recurrence);
    not_starting_there(isession,1)=sum(Start_recurrence<1)/numel(End_recurrence);
    
    clear Start_recurrence End_recurrence
end

M1=strcmp(Areas,'M1');
figure(fig2)
subplot(2,2,2)
title('M1')
box off
subplot(2,2,4)
title('PMd')
xlabel('Time to movement onset [ms]')
ylabel('Number of trajectories in the recurrence region')
box off

[mean(1-not_starting_there(M1)) mean(1-not_starting_there(~M1))]
[mean(1-not_coming_back(M1)) mean(1-not_coming_back(~M1))]

end


function [limit,center]=define_recurrence_region_v3(score,idx_dir,do_plot)
%this version uses the threshold to define a region
rec_length=1;
Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
all_points_beginning=zeros(rec_length*Ndir,size(score,2));
threshold=10;
limit=prctile(pdist(score),threshold);


for i_dir=1:Ndir
      
    all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,:)=score(find(idx_dir==i_dir,rec_length,'First'),:);
    if do_plot
        subplot(2,2,1)
        plot3(score(idx_dir==i_dir,1),score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,1),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)
        
         subplot(2,2,3)
        plot3(score(idx_dir==i_dir,2),score(idx_dir==i_dir,3),score(idx_dir==i_dir,4),'Color',colour_dir(i_dir,:),'LineWidth',1)
        hold on
        plot3(all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,2),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,3),'Color',colour_dir(i_dir,:)/2,'LineWidth',3)

    end
end

center=mean(all_points_beginning);
if do_plot
     min_all=min(min(score(:,[1 2])));
    max_all=max(max(score(:,[1 2])));
    [X,Y,Z] = sphere;
    
    subplot(2,2,1)
    colormap(colormap('white')*0.5)
    mesh(X.*limit+center(1),Y.*limit+center(2),Z.*limit+center(3),'FaceAlpha','0.5')
    
    hold on
    xlim([min_all max_all])
    ylim([min_all max_all])
    zlim([min_all max_all])
    xlabel('PC 1')
    ylabel('PC 2')
    zlabel('PC 3')
    
    subplot(2,2,3)
    colormap(colormap('white')*0.5)
    mesh(X.*limit+center(2),Y.*limit+center(3),Z.*limit+center(4),'FaceAlpha','0.5')
    
    hold on
    xlim([min_all max_all])
    ylim([min_all max_all])
    zlim([min_all max_all])
    xlabel('PC 2')
    ylabel('PC 3')
    zlabel('PC 4')
end

end


function Its=IsInRecRegion2(score,region)
limit=region(1);
center=region(2:end);

Its=pdist2(center,score)<limit;
end