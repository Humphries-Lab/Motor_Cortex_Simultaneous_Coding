function recurrence_all_sessions(sessions,Areas,threshold,Ndir,Nbins,threshold_dist,do_plot_supp)

ndim=zeros(numel(Areas),1);
colour_dir=hsv(Ndir);

upto=200;% number of points to analyse from the beginning of prep

nrows_plot=max(sum(strcmp(Areas,'M1')),sum(strcmp(Areas,'PMd')));

fig1=figure;

if do_plot_supp
    fig_supp=figure;
    j_m1=0;
    j_PMd=0;
end

colormap(flipud(colormap('gray')))

for iarea=1:numel(Areas)
    Area=Areas{iarea};
    session=sessions{iarea};
    %time_rec=nan(Ndir,Nbins);
    
    if do_plot_supp
        
    if strcmp(Area,'M1')
        j_m1=j_m1+1;
        
    else
        j_PMd=j_PMd+1;
    end
    
    end
    
    for i_bin=1:Nbins
        
        load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto','sigma_filter')
        
        ndim(iarea)=find(cumsum(variance)>threshold,1,'First');
        score=score(idx_duration==i_bin,:);
        idx_dir=idx_dir(idx_duration==i_bin,:);
        recurrence=pdist2(score(:,1:ndim(iarea)),score(:,1:ndim(iarea)));
        limit=prctile(pdist(score(:,1:ndim(iarea))),threshold_dist);
        
        
        
        if iarea==1 && i_bin==1
            
            figure(fig1)
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
        
        %% Compute recurrence
        segment_length=round((t_upto(i_bin)-t_from)*1000);
        idx=zeros(Ndir,upto);
        hall=zeros(Ndir,numel(0:50:segment_length)-1);
        Edges=0:50:segment_length;
%        midbin=Edges+50/2;
        
        for i_dir=1:Ndir
            
            recurrencebin=recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(i_dir-1)*segment_length+1:i_dir*segment_length)<limit;
            idx(i_dir,:)=recurrence_time(recurrencebin,upto,2*sigma_filter);
            
            hall(i_dir,:)=histcounts(idx(i_dir,:),Edges,'Normalization','probability');

            %To compute time of recurrene, select bin with maximum number of recurrent points
%             [val,bin]=max(hall(i_dir,:),[],'includenan');
%             
%             if val~=0
%                 time_rec(i_dir,i_bin)=midbin(bin);
%             end
            
            
            if iarea==1 && i_bin==1 && i_dir==1
                
                figure(fig1)
                subplot(1,4,2)
                hold on
                imagesc([t_from*1000 t_upto(i_bin)*1000],[t_upto(i_bin)*1000 t_from*1000],flipud(recurrencebin))
                colormap(flipud(colormap('gray')))
                plot([0 0],[t_from*1000 t_upto(i_bin)*1000],'--k')
                plot([t_from*1000 t_upto(i_bin)*1000],[0 0],'--k')
                plot([250 250],[t_from*1000 t_upto(i_bin)*1000],'--k')
                plot([t_from*1000 t_upto(i_bin)*1000],[250 250],'--k')
                xlim([t_from*1000 t_upto(i_bin)*1000])
                ylim([t_from*1000 t_upto(i_bin)*1000])
                xlabel('Time from movement onset [ms]')
                ylabel('Time from movement onset [ms]')
                box off
            end
        end
        
        
        
        %% Compute proportion of recurrence points across directions
        %         proportion=zeros(Ndir,Ndir);
        %         for i_dir=1:Ndir
        %             recurrence_matrix=sum(sum(recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(i_dir-1)*segment_length+1:i_dir*segment_length)<limit));
        %
        %             for k=1:Ndir
        %                 recurrence_matrix_other=recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(k-1)*segment_length+1:k*segment_length)<limit;
        %                 proportion(i_dir,k)=sum(recurrence_matrix_other(:))/recurrence_matrix;
        %             end
        %
        %         end
        
        %% Plots !!
        
        
        
        if i_bin==1
            if do_plot_supp
                figure(fig_supp)
                if strcmp(Area,'M1')
                    subplot(nrows_plot,2,(j_m1-1)*2+1)
                else
                    subplot(nrows_plot,2,(j_PMd-1)*2+2)
                end
                
                imagesc((0:50:segment_length)+t_from*1000,1:Ndir,hall)
                title([session(1:2) ' ' session(4:5) ' ' Area])
                caxis([0 1])
                colormap(flipud(colormap('gray')))
                box off
                
                if j_m1==sum(strcmp(Areas,'M1'))
                    subplot(nrows_plot,2,(j_m1-1)*2+1)
                    xlabel('Time from movement onset [ms]')
                end
                if j_PMd==sum(strcmp(Areas,'PMd'))
                    subplot(nrows_plot,2,(j_PMd-1)*2+2)
                    xlabel('Time from movement onset [ms]')
                end
            end
            % Examples figure 3
            % M1 (MC S1)
            if iarea==1
                figure(fig1)
                subplot(2,4,4)
                hold on
                imagesc((0:50:segment_length)+t_from*1000,1:Ndir,hall)
                plot([0 0],[0 8.5],'--k')
                plot([300 300],[0 8.5],'--k')
                caxis([0 1])
                ylim([0.5 8.5])
                title('M1')
                set(gca, 'YDir','reverse')
                colormap(flipud(colormap('gray')))
                box off
            end
            % PMd
            if iarea==11
                figure(fig1)
                subplot(2,4,8)
                hold on
                imagesc((0:50:segment_length)+t_from*1000,1:Ndir,hall)
                plot([0 0],[0 8.5],'--k')
                plot([300 300],[0 8.5],'--k')
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
    
end
end

function times_idx=recurrence_time(recurrencebin,upto,sigma_filter)
idx=nan(1,upto);

for point=1:upto
    %find last contiguous point
    
    tmp=find(diff(recurrencebin(point:end,point))>0,1,'last');
    if ~isempty(tmp)
        if tmp>=sigma_filter
            idx(point)=tmp;
        end
    end
    clear tmp
end

times_idx=idx;
end
