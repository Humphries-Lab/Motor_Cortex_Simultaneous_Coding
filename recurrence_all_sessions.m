function recurrence_all_sessions(sessions,Areas,threshold,Ndir,threshold_dist)
ndim=zeros(numel(Areas),1);
%thetacolor_axis=(-pi+0.0001:2*pi/Ndir:pi)+pi/Ndir;
colour_dir=hsv(Ndir);

fig1=figure;
fig2=figure;

colormap(flipud(colormap('gray')))
upto=200;% number of points to analyse from the beginning of prep
nrows=max(sum(strcmp(Areas,'M1')),sum(strcmp(Areas,'PMd')));
j_m1=0;
j_PMd=0;

for iarea=1:numel(Areas)
    Area=Areas{iarea};
    session=sessions{iarea};
    time_rec=zeros(Ndir,4);
    if strcmp(Areas{iarea},'M1')
        colourArea='m';
        j_m1=j_m1+1;
        
    else
        colourArea='b';
        j_PMd=j_PMd+1;
    end
    
    for i_bin=1:4
        
        %load(['scores_LDS_diff_duration_' session(1:5) '_' Area '.mat'],'score','idx_dir','idx_duration')
        load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','from','ms')
        ndim(iarea)=find(cumsum(variance)>threshold,1,'First');
        score=score(idx_duration==i_bin,:);
        idx_dir=idx_dir(idx_duration==i_bin,:);
        recurrence=pdist2(score(:,1:ndim(iarea)),score(:,1:ndim(iarea)));
        limit=prctile(pdist(score(:,1:ndim(iarea))),threshold_dist);
        proportion=zeros(Ndir,Ndir);
        
        
        if iarea==1 && i_bin==1
            
            figure(fig2)
            subplot(1,4,1)
            plot3(score(idx_dir==1,1),score(idx_dir==1,2),score(idx_dir==1,3),'Color',colour_dir(1,:,:))
            hold on
            plot3(score(find(idx_dir==1,1,'First'),1),score(find(idx_dir==1,1,'First'),2),score(find(idx_dir==1,1,'First'),3),'.','MarkerSize',12,'Color',colour_dir(1,:,:))
            
            xlabel('PC 1')
            ylabel('PC 2')
            zlabel('PC 3')
            
            subplot(1,4,3)
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
        
        for i_dir=1:Ndir
            
            recurrencebin=recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(i_dir-1)*segment_length+1:i_dir*segment_length)<limit;
            idx(i_dir,:)=recurrence_time(recurrencebin,upto,2*ms);
            figure(fig1)
            
            Edges=0:50:segment_length;
            midbin=Edges+50/2;

            hall(i_dir,:)=histcounts(idx(i_dir,:),Edges,'Normalization','probability');
            %v1- median of all points
            %time_rec(i,i_bin)=nanmedian(idx(i,:));
            %v2- select bin with maximum number of recurrent points
            [val,bin]=max(hall(i_dir,:),[],'includenan');
            
            if val==0
                time_rec(i_dir,i_bin)=nan;
            else
                time_rec(i_dir,i_bin)=midbin(bin);
            end
           
            
            if iarea==1 && i_bin==1 && i_dir==1
                
                figure(fig2)
                subplot(1,4,2)
                hold on
                imagesc([t_1*1000 t_2(i_bin)*1000],[t_2(i_bin)*1000 t_1*1000],flipud(recurrencebin))
                colormap(flipud(colormap('gray')))
                plot([0 0],[t_1*1000 t_2(i_bin)*1000],'--k')
                plot([t_1*1000 t_2(i_bin)*1000],[0 0],'--k')
                plot([250 250],[t_1*1000 t_2(i_bin)*1000],'--k')
                plot([t_1*1000 t_2(i_bin)*1000],[250 250],'--k')
                xlim([t_1*1000 t_2(i_bin)*1000])
                ylim([t_1*1000 t_2(i_bin)*1000])
                xlabel('Time from movement onset [ms]')
                ylabel('Time from movement onset [ms]')
                box off
            end
        end
        
        
        
        %% Compute proportion of recurrence points across directions
        figure(fig1)
        
        for i_dir=1:Ndir
            recurrence_matrix=sum(sum(recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(i_dir-1)*segment_length+1:i_dir*segment_length)<limit));
            
            
            for k=1:Ndir
                recurrence_matrix_other=recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(k-1)*segment_length+1:k*segment_length)<limit;
                if k~=i_dir
                    tmp=diag(recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(k-1)*segment_length+1:k*segment_length));
                    [~,idx_max(k)]=max(tmp(1:-t_1*1000));
                    [~,idx_max_mov(k)]=max(tmp(-t_1*1000+1:end));
                else
                    idx_max(k)=nan;
                    idx_max_mov(k)=nan;
                    tmp=nan;
                end
                %plot(t_1*1000:(t_2*1000)-1,tmp,'Color',colormap2(k,:))
                proportion(i_dir,k)=sum(recurrence_matrix_other(:))/recurrence_matrix;
            end
            
            max_dist_av(i_dir)=nanmean(idx_max+t_1*1000);
            max_dist_av_mov(i_dir)=nanmean(idx_max_mov+t_1*1000);
            clear idx_max
            
            
        end
        
        if i_bin==1
            
            
            figure(fig1)
            if strcmp(colourArea,'m')
                subplot(nrows,2,(j_m1-1)*2+1)
            else
                subplot(nrows,2,(j_PMd-1)*2+2)
            end
            
            imagesc((0:50:segment_length)+t_1*1000,1:Ndir,hall)
            title([session(1:2) ' ' session(4:5) ' ' Area])
            caxis([0 1])
            colormap(flipud(colormap('gray')))
            box off
            if j_m1==sum(strcmp(Areas,'M1'))
                subplot(nrows,2,(j_m1-1)*2+1)
                xlabel('Time from movement onset [ms]')
            end
            if j_PMd==sum(strcmp(Areas,'PMd'))
                subplot(nrows,2,(j_PMd-1)*2+2)
                xlabel('Time from movement onset [ms]')
            end
            
            % Examples figure 3
            % M1 (MC S1)
            if iarea==1
                figure(fig2)
                subplot(2,4,4)
                hold on
                imagesc((0:50:segment_length)+t_1*1000,1:Ndir,hall)
                plot([0 0],[0 8.5],'--k')
                plot([250 250],[0 8.5],'--k')
                caxis([0 1])
                ylim([0.5 8.5])
                title('M1')
                set(gca, 'YDir','reverse')
                colormap(flipud(colormap('gray')))
                box off
            end
            % PMd
            if iarea==11
                figure(fig2)
                subplot(2,4,8)
                hold on
                imagesc((0:50:segment_length)+t_1*1000,1:Ndir,hall)
                plot([0 0],[0 8.5],'--k')
                plot([250 250],[0 8.5],'--k')
                caxis([0 1])
                ylim([0.5 8.5])
                title('PMd')
                set(gca, 'YDir','reverse')
                colormap(flipud(colormap('gray')))
                box off
                xlabel('Time to movement onset [ms]')
            end
        end
        
        
        
        
    end
% plotting the recurrence time average per session    
%     figure(fig1)
%     
%     subplot(4,3,12)
%     errorbar((t_2(1:i_bin)-t_1)*1000, nanmean(time_rec),nanstd(time_rec),'Color',colourArea)
%     hold on
%     plot([600 1100],[600 1100],'r')
%     box off
%     xlabel('Trajectory duration [ms]')
%     ylabel('Recurrence time [ms]')
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
