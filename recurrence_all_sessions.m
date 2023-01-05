function recurrence_all_sessions(sessions,Areas,threshold,Ndir,threshold_dist)
ndim=zeros(numel(Areas),1);
%thetacolor_axis=(-pi+0.0001:2*pi/Ndir:pi)+pi/Ndir;
colormap2=hsv(Ndir);

fig1=figure;
fig2=figure;
colormap(flipud(colormap('gray')))
upto=200;% number of points to analyse from the beginning of prep
nrows=max(sum(strcmp(Areas,'M1')),sum(strcmp(Areas,'PMd')));
j_m1=0;
j_PMd=0;

for j=1:numel(Areas)
    Area=Areas{j};
    session=sessions{j};
    time_rec=zeros(Ndir,4);
    if strcmp(Areas{j},'M1')
        colourArea='m';
        j_m1=j_m1+1;
        
    else
        colourArea='b';
        j_PMd=j_PMd+1;
    end
    
    for i_bin=1:4
        
        %load(['scores_LDS_diff_duration_' session(1:5) '_' Area '.mat'],'score','idx_dir','idx_duration')
        load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','from','ms')
        ndim(j)=find(cumsum(variance)>threshold,1,'First');
        score=score(idx_duration==i_bin,:);
        idx_dir=idx_dir(idx_duration==i_bin,:);
        recurrence=pdist2(score(:,1:ndim(j)),score(:,1:ndim(j)));
        limit=prctile(pdist(score(:,1:ndim(j))),threshold_dist);
        proportion=zeros(Ndir,Ndir);
        
        
        if j==1 && i_bin==1
            
            figure(fig2)
            subplot(1,3,1)
            plot3(score(idx_dir==1,1),score(idx_dir==1,2),score(idx_dir==1,3),'Color',colormap2(1,:,:))
            hold on
            plot3(score(find(idx_dir==1,1,'First'),1),score(find(idx_dir==1,1,'First'),2),score(find(idx_dir==1,1,'First'),3),'.','MarkerSize',12,'Color',colormap2(1,:,:))
            
            xlabel('PC 1')
            ylabel('PC 2')
            zlabel('PC 3')
            
            subplot(1,3,3)
            imagesc(flipud(recurrence<limit))
            hold on
            xlim([1 size(recurrence,1)])
            ylim([1 size(recurrence,1)])
            colormap(flipud(colormap('gray')))
            box off
            total_time=size(score,1);
            for i_line=1:Ndir
                time_line=i_line*sum(idx_dir==i_line);
                plot([time_line time_line],[0 total_time],'r')
                plot([0 total_time],[time_line time_line],'r')
            end
        end
        
        %figure(j+1) %per bin
        
        %% compute recurrence time
        segment_length=round((t_2(i_bin)-t_1)*1000);
        idx=zeros(Ndir,upto);
        hall=zeros(Ndir,numel(0:50:segment_length)-1);
        
        for i=1:Ndir
            
            recurrencebin=recurrence((i-1)*segment_length+1:i*segment_length,(i-1)*segment_length+1:i*segment_length)<limit;
            %             plot3(score(idx_dir==i,1),score(idx_dir==i,2),score(idx_dir==i,3),'Color',colormap2(i,:))
            %             pause
            %             cla
            idx(i,:)=recurrence_time(recurrencebin,upto,2*ms);
            figure(fig1)
            subplot(nrows,3,3)
            h=histogram(idx(i,:),0:50:segment_length,'Normalization','Probability','Visible','off');
            hall(i,:)=h.Values;
            %v1- median of all points
            %time_rec(i,i_bin)=nanmedian(idx(i,:));
            %v2- select bin with maximum number of recurrent points
            [val,bin]=max(hall(i,:),[],'includenan');
            
            if val==0
                time_rec(i,i_bin)=nan;
            else
                time_rec(i,i_bin)=h.BinEdges(bin)+h.BinWidth/2;
            end
             cla
            if j==1 && i_bin==1 && i==1
                
                figure(fig2)
                subplot(1,3,2)
                imagesc(flipud(recurrencebin))
                colormap(flipud(colormap('gray')))
                %pause
                box off
            end
            %figure(j+1)
        end
        
%         subplot(1,4,i_bin)
%         imagesc(0:50:segment_length,1:Ndir,hall)
%         caxis([0 1])
        
        %% Compute proportion of recurrence points across directions
        figure(fig1)
        
        for i=1:Ndir
            recurrence_matrix=sum(sum(recurrence((i-1)*segment_length+1:i*segment_length,(i-1)*segment_length+1:i*segment_length)<limit));
            
            
            for k=1:Ndir
                recurrence_matrix_other=recurrence((i-1)*segment_length+1:i*segment_length,(k-1)*segment_length+1:k*segment_length)<limit;
                if k~=i
                    tmp=diag(recurrence((i-1)*segment_length+1:i*segment_length,(k-1)*segment_length+1:k*segment_length));
                    [~,idx_max(k)]=max(tmp(1:-t_1*1000));
                    [~,idx_max_mov(k)]=max(tmp(-t_1*1000+1:end));
                else
                    idx_max(k)=nan;
                    idx_max_mov(k)=nan;
                    tmp=nan;
                end
                %plot(t_1*1000:(t_2*1000)-1,tmp,'Color',colormap2(k,:))
                proportion(i,k)=sum(recurrence_matrix_other(:))/recurrence_matrix;
            end
            
            max_dist_av(i)=nanmean(idx_max+t_1*1000);
            max_dist_av_mov(i)=nanmean(idx_max_mov+t_1*1000);
            clear idx_max
            
            
        end
        
        if i_bin==1
            if strcmp(colourArea,'m')
                subplot(nrows,3,(j_m1-1)*3+1)
            else
                subplot(nrows,3,(j_PMd-1)*3+2)
            end
            imagesc((0:50:segment_length)+t_1*1000,1:Ndir,hall)
            title([session(1:2) ' ' session(4:5) ' ' Area])
            caxis([0 1])
            colormap(flipud(colormap('gray')))
            box off
            pause
            if j_m1==sum(strcmp(Areas,'M1'))
            subplot(nrows,3,(j_m1-1)*3+1)
            xlabel('Time from movement onset [ms]')
            end
            if j_PMd==sum(strcmp(Areas,'PMd'))
            subplot(nrows,3,(j_PMd-1)*3+2)
            xlabel('Time from movement onset [ms]')
            end
        end
        
        
        
        
    end
    
    figure(fig1)
    
    subplot(4,3,12)
    %(t_2(1:i_bin)-t_1)*1000
    errorbar((t_2(1:i_bin)-t_1)*1000, nanmean(time_rec),nanstd(time_rec),'Color',colourArea)
    hold on
    plot([600 1100],[600 1100],'r')
    box off
    xlabel('Trajectory duration [ms]')
    ylabel('Recurrence time [ms]')
end
end

function times_idx=recurrence_time(recurrencebin,upto,ms)
idx=nan(1,upto);
%
% subplot(2,2,1)
% imagesc(recurrencebin(1:end,1:upto))

for point=1:upto
    %find last contiguous point
    
    tmp=find(diff(recurrencebin(point:end,point))>0,1,'last');
    if ~isempty(tmp)
        if tmp>=ms
            idx(point)=tmp;
        end
    end
    clear tmp
end

times_idx=idx;
end
