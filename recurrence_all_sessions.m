function recurrence_all_sessions(Sessions,Areas,threshold,Ndir,i_bin,threshold_dist,do_plot_supp)
%% recurrence_all_sessions plots the recurrence plots and the histograms of recurrence for the selected recordings and bin number
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
% i_bin: duration bin to be analysed
%
% threshold_dist: threshold of the percentil of distance(theta) for a point to be
% considered recurrent. threshold_dist ranges from 0 to 100.
%
% do_plot_supp: 1- plot histograms of recurrence for all recordings
%               0- omit the plot above
%
% 26/01/2023
% Andrea Colins Rodriguez
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
    session=Sessions{iarea};


    if do_plot_supp

        if strcmp(Area,'M1')
            j_m1=j_m1+1;

        else
            j_PMd=j_PMd+1;
        end

    end



    load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto','sigma_filter')


    if iarea==1
        figure(fig1)
       
        plot_initial_condition(score,idx_dir,idx_duration,round(-t_from*1000),round(1000*t_upto),colour_dir)

    end

    ndim(iarea)=find(cumsum(variance)>threshold,1,'First');
    score=score(idx_duration==i_bin,:);
    idx_dir=idx_dir(idx_duration==i_bin,:);
    recurrence=pdist2(score(:,1:ndim(iarea)),score(:,1:ndim(iarea)));
    limit=prctile(pdist(score(:,1:ndim(iarea))),threshold_dist);



    if iarea==1

        figure(fig1)
        subplot(2,4,1)
        plot3(score(idx_dir==1,1),score(idx_dir==1,2),score(idx_dir==1,3),'Color',colour_dir(1,:,:))
        hold on
        start=find(idx_dir==1,1,'First');
        plot3(score(start,1),score(start,2),score(start,3),'.','MarkerSize',12,'Color',colour_dir(1,:,:))

        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')

        subplot(2,4,3)
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

    for i_dir=1:Ndir

        recurrencebin=recurrence((i_dir-1)*segment_length+1:i_dir*segment_length,(i_dir-1)*segment_length+1:i_dir*segment_length)<limit;
        idx(i_dir,:)=recurrence_time(recurrencebin,upto,2*sigma_filter);

        hall(i_dir,:)=histcounts(idx(i_dir,:),Edges,'Normalization','probability');

        %% Plots !!

        if iarea==1 && i_dir==1

            figure(fig1)
            subplot(2,4,2)
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
        subplot(4,4,4)
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
        subplot(4,4,8)
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

function times_idx=recurrence_time(recurrencebin,first_points,sigma_filter)

%% recurrence_times finds the times in recurrencebin matrix that recur to the first_points
%
% INPUTS
% recurrencebin= recurrence matrix. i.e binarised distance between points
% in square format
%
% first_points= number of time bins to be tested from the start of the
% trajectories [ms]
%
% sigma_filter= time bins lower bound to be considered recurrence. If a
% time of recurrence is less than sigma_filer, the point is considered not
% recurrent.
%
% OUTPUTS
%
% times_idx= array containing the times at which the points recur. Nan
% values indicate non-recurrence.
%
% 26/01/2023
% Andrea Colins Rodriguez


times_idx=nan(1,first_points);

for point=1:first_points

    %find last contiguous point
    tmp=find(diff(recurrencebin(point:end,point))>0,1,'last');

    if ~isempty(tmp)
        if tmp>=sigma_filter
            times_idx(point)=tmp;
        end
    end
    clear tmp
end


end


function plot_initial_condition(score,idx_dir,idx_duration,t_from,t_upto,colour_dir)
Ndir=max(idx_dir);
Ndur=max(idx_duration);
colour_dur=plasma(Ndur);
for i_bin=1:Ndur
    for i_dir=[1 5 6 7]
        idx=find(idx_dir==i_dir & idx_duration==i_bin);
        subplot(2,8,[11 12 13])
        plot3(score(idx(t_from:t_upto(i_bin)),1),score(idx(t_from:t_upto(i_bin)),2),score(idx(t_from:t_upto(i_bin)),3),'Color',colour_dir(i_dir,:)./sqrt(i_bin))
        hold on
        %plot3(score(idx(t_from-50:t_from),1),score(idx(t_from-50:t_from),2),score(idx(t_from-50:t_from),3),'--','Color',colour_dir(i_dir,:)./sqrt(i_bin))
        plot3(score(idx(t_from),1),score(idx(t_from),2),score(idx(t_from),3),'s','MarkerSize',8,'Color',colour_dur(i_bin,:),'MarkerFaceColor',colour_dur(i_bin,:))
        
        subplot(2,8,[14 15 16])
        plot3(score(idx(t_from:t_upto(i_bin)),1),score(idx(t_from:t_upto(i_bin)),2),score(idx(t_from:t_upto(i_bin)),3),'Color',colour_dir(i_dir,:)./sqrt(i_bin))
        hold on
        %plot3(score(idx(t_from-50:t_from),1),score(idx(t_from-50:t_from),2),score(idx(t_from-50:t_from),3),'--','Color',colour_dir(i_dir,:)./sqrt(i_bin))
        plot3(score(idx(t_from),1),score(idx(t_from),2),score(idx(t_from),3),'s','MarkerSize',8,'Color',colour_dur(i_bin,:),'MarkerFaceColor',colour_dur(i_bin,:))
       
    end
end

subplot(2,8,[11 12 13])
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

subplot(2,8,[14 15 16])
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
view(-90.6429,73.5081)
end