function [h,p,pSI]=distance_duration_vs_direction_all_sessions(sessions,Areas,threshold,Ndir,Nbins,extra_plot)
ndim=zeros(numel(sessions),1);
colour_dur=plasma(4);
%colour_dir=hsv(Ndir);
final_length=600;
normalised_t=linspace(0,1,final_length);
Nsessions=numel(sessions);
h=zeros(1,Nsessions);
p=zeros(1,Nsessions);
SIj=zeros(1,Nsessions);
SI2j=zeros(1,Nsessions);
pR2=zeros(1,Nsessions);
pSI=zeros(1,Nsessions);
fig1=figure;
fig2=figure;
Hdist=[];
Hdur=[];
for isession=1:Nsessions
    
    Area=Areas{isession};
    session=sessions{isession};
    load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','nsamples_condition')

    %load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','ms','t_1','t_2','from','nsamples_condition')
    ndim(isession)=find(cumsum(variance)>threshold,1,'First');
    
    %% Plot example trajectories for same direction different durations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if extra_plot && isession==1
        figure(fig1)
        i_dir=1;
        
        for i_bin=1:Nbins
            idx=find(idx_dir==i_dir & idx_duration==i_bin);
            
            subplot(3,2,[1 3 5])
            hold on
            plot3(score(idx(1:600),1),score(idx(1:600),2),score(idx(1:600),3),'Color',colour_dur(i_bin,:),'LineWidth',2)
            plot3(score(idx(1),1),score(idx(1),2),score(idx(1),3),'o','MarkerFaceColor',colour_dur(i_bin,:),'MarkerEdgeColor',colour_dur(i_bin,:))
            
            %normalise time
            xtime2=linspace(0,1,numel(idx));
            % interpolate to match all vector lengths
            norm_score= interp1(xtime2,score(idx,1:ndim(isession)),normalised_t);
            
            subplot(3,2,2)
            hold on
            plot(normalised_t,norm_score(:,1),'Color',colour_dur(i_bin,:))
            
            subplot(3,2,4)
            hold on
            plot(normalised_t,norm_score(:,2),'Color',colour_dur(i_bin,:))
            
            subplot(3,2,6)
            hold on
            plot(normalised_t,norm_score(:,3),'Color',colour_dur(i_bin,:))
            
            
        end
        
        subplot(3,2,[1 3 5])
        xlim([min(score(:,1)) max(score(:,1))])
        ylim([min(score(:,2)) max(score(:,2))])
        zlim([min(score(:,3)) max(score(:,3))])
        xlabel('PC 1')
        ylabel('PC 2')
        zlabel('PC 3')
        
        subplot(3,2,2)
        ylabel('PC 1')
        box off
        xlim([0 1])
        
        subplot(3,2,4)
        ylabel('PC 2')
        box off
        xlim([0 1])
        
        subplot(3,2,6)
        box off
        xlim([0 1])
        xlabel('Normalised time')
        ylabel('PC 3')
    end
    
    
    %% Calculate R^2 for trajectories corresponding to same direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(Areas{isession},'M1')
        colourArea=[85 30 116]./256;
    else
        colourArea=[89 156 153]./256;
    end
    
    [SI,pSI(isession),~,SIstd,SIall]=scale_index(score,idx_dir,idx_duration,ndim(isession),normalised_t);
    [SI2,~,~,SI2std,SIallother]=scale_index_other(score,idx_dir,idx_duration,ndim(isession),normalised_t);
    % Stats
    SIj(isession)=SI;
    SI2j(isession)=SI2;
    [~,pR2(isession)]=ttest2(SIall,SIallother);
    
    figure(fig2)
    subplot(2,3,1)
    errorbar(isession,SI,SIstd,'.','Color',colourArea)
    hold on
    errorbar(isession,SI2,SI2std,'.','Color',[0.5 0.5 0.5])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test if Hausdorff distance for different durations is shorter than for different directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [fraction_above,Delta_distances,Hdist_tmp,Hdur_tmp]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim(isession),1,colourArea,nsamples_condition);
    %% add a plot with the histogram for each axis
    Hdist=[Hdist;Hdist_tmp];
    Hdur=[Hdur;Hdur_tmp];
    %% add a plot of the Delta distances for each session to show that they are over zero.
    [h(isession),p(isession)]=ttest(Delta_distances);
    
    subplot(2,3,3)
    errorbar(isession,mean(Delta_distances),std(Delta_distances),'.','Color',colourArea)
    hold on
    
    clear Delta_distances fraction_above idx_duration idx_dir
end
subplot(2,3,1)
plot([0 Nsessions],[0 0],'Color',[0.5 0.5 0.5])
box off
ylabel('R^2')
xlabel('Recording number')

text(1,0.9,['mean = ' num2str(mean(SIj),2)])

text(1,-1,['mean = ' num2str(mean(SI2j),2)])



subplot(2,3,3)
plot([0 Nsessions],[0 0],'Color',[0.5 0.5 0.5])
box off
ylabel('\Delta Distance')
xlabel('Recording')
text(1,-0.001,['max p-value = ' num2str(max(p),'%.e')])


subplot(2,3,2)
xlim([0 0.025])
ylim([0 0.025])

axes('Position',[0.41,0.93,0.21,0.041]);
histogram(Hdist,0:0.001:0.025,'Normalization','probability')
xlim([0 0.025])
box off

axes('Position',[0.63,0.58,0.028,0.338]);
histogram(Hdur,0:0.001:0.025,'Normalization','probability','orientation','horizontal')
ylim([0 0.025])
box off
end

function [total_fraction,Delta_distances,All_dist_dir,All_dist_dur]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim,do_plot,ColourArea,nsamples_condition)
Ndir=max(idx_dir);
Nbins=max(idx_duration);
Hdist_dir=zeros(Ndir,Ndir);
%Hdist_dur=nan(Ndir,sum(1:Nbins-1));
idx_others=1:Ndir;
idx_others(round(Ndir/2))=[];
if do_plot
    colour_dir=hsv(Ndir);
end

for i_dir=1:Ndir
    counter=1;
    %%trajectories with same direction but different durations
    for i_bin=1:Nbins-1
        idx1=find(idx_dir==i_dir & idx_duration==i_bin);
        for j_bin=i_bin+1:Nbins
            idx2=find(idx_dir==i_dir & idx_duration==j_bin);
            %% distance between trajectories
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            
            Hdist_dur(i_dir,counter)=max([min(recurrence),min(recurrence,[],2)']);
            nsamp(i_dir,counter)=min(nsamples_condition(i_dir,i_bin),nsamples_condition(i_dir,j_bin));
            counter=counter+1;
        end
    end
end
% if do_plot
%     subplot(3,3,1)
%     hold on
%     Delta_ms=[100 200 300];
%     Delta100=Hdist_dur(:,[1 4 6]);
%     Delta200=Hdist_dur(:,[2 5]);
%     Delta300=Hdist_dur(:,3);
%     HausdorffD_mean=[mean(Delta100(:)) mean(Delta200(:)) mean(Delta300(:))];
%     HausdorffD_SD=[std(Delta100(:)) std(Delta200(:)) std(Delta300(:))];
%     %[hdelta200,pdelta200]=ttest2(Delta100(:),Delta300(:),'Vartype','unequal');
%     [hdelta100,pdelta100]=ttest2(Delta100(:),Delta200(:),'alpha',0.025);
%     [hdelta100b,pdelta100]=ttest2(Delta200(:),Delta300(:),'alpha',0.025);
%     hdelta200=hdelta100 & hdelta100b;
%     errorbar(Delta_ms,HausdorffD_mean,HausdorffD_SD,'Color',ColourArea)
%     %plot(Delta_dur,mean(Hdist_dur)','Color',ColourArea)
%     box off
%     xlabel('Delta dur')
%     ylabel('H Distance')
% end
Delta_distances=[];
counter_points_above=[];
All_dist_dur=[];
All_dist_dir=[];
for i_bin=1:Nbins
    
    for i_dir=1:Ndir
        idx1=find(idx_dir==i_dir & idx_duration==i_bin);
        for j_dir=1:Ndir
            idx2=find(idx_dir==j_dir & idx_duration==i_bin);
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            
            Hdist_dir(i_dir,j_dir,i_bin)=max([min(recurrence),min(recurrence,[],2)']);
        end
        tmp=circshift(Hdist_dir(i_dir,:,i_bin),round(Ndir/2)-i_dir);
        %min_tmp=min(tmp(idx_others));
        %counter_points_below=[counter_points_below, min_tmp>Hdist_dur(i_dir,:)];
        
        counter_points_above=[counter_points_above, tmp(Ndir/2-1)>Hdist_dur(i_dir,:) ,tmp(Ndir/2+1)>Hdist_dur(i_dir,:)];
        Delta_distances=[Delta_distances,tmp(Ndir/2-1)-Hdist_dur(i_dir,:),tmp(Ndir/2+1)-Hdist_dur(i_dir,:)];
        All_dist_dur=[All_dist_dur;Hdist_dur(i_dir,:)';Hdist_dur(i_dir,:)'];
        All_dist_dir=[All_dist_dir;ones(size(Hdist_dur,2),1)*tmp(Ndir/2-1);ones(size(Hdist_dur,2),1)*tmp(Ndir/2+1)];
        if do_plot
            %             subplot(3,3,2)
            %
            %             plot(tmp,'.-','Color',colour_dir(i_dir,:))
            %             hold on
            %
            %             plot(ones(size(Hdist_dur,2))*Ndir/2,Hdist_dur(i_dir,:),'.','Color',colour_dir(i_dir,:))
            %
            subplot(2,3,2)
            hold on
            if strcmp(ColourArea,'')
                plot(Hdist_dur(i_dir,:),ones(size(Hdist_dur,2),1)*tmp(Ndir/2-1),'.','Color',colour_dir(i_dir,:))
                plot(Hdist_dur(i_dir,:),ones(size(Hdist_dur,2),1)*tmp(Ndir/2+1),'.','Color',colour_dir(i_dir,:))
            else
                transparency=nsamp./max(nsamp(:));
                for iplot=1:size(Hdist_dur,2)
                    s=scatter(Hdist_dur(i_dir,iplot),tmp(Ndir/2-1),8,ColourArea,'filled');
                    alpha(s,transparency(i_dir,iplot))
                    s=scatter(Hdist_dur(i_dir,iplot),tmp(Ndir/2+1),8,ColourArea,'filled');
                    alpha(s,transparency(i_dir,iplot))
                    
                end
                
                %                 s=scatter(Hdist_dur(i_dir,:),ones(size(Hdist_dur,2),1)*tmp(Ndir/2+1),ColourArea,'filled');
                %                 alpha(s,transparency)
                %plot(Hdist_dur(i_dir,:),ones(size(Hdist_dur,2),1)*tmp(Ndir/2-1),'.','Color',ColourArea,'FaceAlpha',0.5)
                %plot(Hdist_dur(i_dir,:),ones(size(Hdist_dur,2),1)*tmp(Ndir/2+1),'.','Color',ColourArea)
            end
        end
    end
    
    
end

total_fraction=sum(counter_points_above)/sum(~isnan(counter_points_above(:)));
if do_plot
    subplot(2,3,2)
    %text(0.015,0.015,['Fraction Above ' num2str(total_fraction)])
    plot([0 0.02],[0 0.02],'r')
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dir),'.k');
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dur),'horizontal','.k');
    
    box off
    xlabel('Distance durations')
    ylabel('Distance direction')
    
end

end


function [SI,p,SI2,SIstd,SIall]=scale_index(score,idx_dir,idx_duration,ndim,normalised_t)
% define scaling index as difference between non scaled
% version / distance scaled version
SI=zeros(max(idx_dir),max(idx_duration)-1);
for i_dir=1:max(idx_dir)
    counter=1;
    
    for i_bin=1:max(idx_duration)
        idx=find(idx_dir==i_dir & idx_duration==i_bin);
        for j_bin=i_bin+1:max(idx_duration)
            idx_j=find(idx_dir==i_dir & idx_duration==j_bin);
            [SI(i_dir,counter),SI2(i_dir,counter)]=scale_index_R2(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
            %SI(i_dir,counter)=scale_index_two_TR(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
            counter=counter+1;
        end
    end
    %% Example for the paper (Fig 3B)
    %     R2example=mean(SI(i_dir,:))
    %     pause
end

[~,p]=ttest(SI(:));
SIall=SI(:);
SIstd=std(SI(:))/sqrt(numel(SI));
SI=nanmean(SI(:));
SI2=nanmean(SI2(:));
end

function [SI,p,SI2,SIstd,SIall]=scale_index_other(score,idx_dir,idx_duration,ndim,normalised_t)
% define scaling index as difference between non scaled
% version / distance scaled version
SI=zeros(max(idx_dir),max(idx_duration)-1);
for i_dir=1:max(idx_dir)
    counter=1;
    for j_dir=i_dir+1:max(idx_dir)
        
        for i_bin=1:max(idx_duration)
            idx=find(idx_dir==i_dir & idx_duration==i_bin);
            
            for j_bin=i_bin+1:max(idx_duration)
                idx_j=find(idx_dir==j_dir & idx_duration==j_bin);
                %[i_dir j_dir i_bin j_bin]
                [SI(i_dir,counter),SI2(i_dir,counter)]=scale_index_R2(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
                %SI(i_dir,counter)=scale_index_two_TR(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
                counter=counter+1;
            end
        end
    end
end
[~,p]=ttest(SI(:));
SIall=SI(:);
SIstd=std(SI(:))/sqrt(numel(SI));
SI=nanmean(SI(:));
SI2=nanmean(SI2(:));
end


function [SI,SI2]=scale_index_R2(X,Y,normalised_t)
%normalise time
Xtime=linspace(0,1,size(X,1));
Ytime=linspace(0,1,size(Y,1));
% interpolate to match all vector lengths
Xnorm= interp1(Xtime,X,normalised_t);
Ynorm= interp1(Ytime,Y,normalised_t);
C=min(size(X,1),size(Y,1));

SI2=R2(X(1:C,:),Y(1:C,:));
scaled_distance=R2(Xnorm,Ynorm);
SI=scaled_distance;

end
function r=R2(X,Y)
X2=reshape(X,numel(X),1);
Y2=reshape(Y,numel(Y),1);
ymean=mean(Y2);
Sres=sum((Y2-X2).^2);
Stot=sum((Y2-ymean).^2);
r=1-Sres/Stot;
end