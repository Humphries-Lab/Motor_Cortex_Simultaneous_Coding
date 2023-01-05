function [h,p,pSI]=distance_duration_vs_direction_all_sessions(sessions,Areas,threshold,Ndir,Nbins,extra_plot)
ndim=zeros(numel(sessions),1);
do_all_dim=0;
colour_dur=plasma(4);
colour_dir=hsv(Ndir);
final_length=600;
normalised_t=linspace(0,1,final_length);
Nsessions=numel(sessions);
h=zeros(1,Nsessions);
p=zeros(1,Nsessions);
fig1=figure;
fig2=figure;
Hdist=[];
Hdur=[];
for j=1:numel(sessions)
    
    Area=Areas{j};
    session=sessions{j};
    %load(['scores_LDS_diff_duration_' session(1:5) '_' Area '.mat'],'score','idx_dir','idx_duration')
    load(['scores_LDS_diff_duration_newfilter_' session '_' Area '.mat'],'score','idx_dir','idx_duration','variance','ms','t_1','t_2','from','nsamples_condition')
    %video2(score,idx_dir,idx_duration,t_1,t_2)
    %pause
    if strcmp(Areas{j},'M1')
        colourArea=[85 30 116]./256;
    else
        colourArea=[89 156 153]./256;
    end
    ndim(j)=find(cumsum(variance)>threshold,1,'First');
    [SI,pSI(j),~,SIstd,SIall]=scale_index(score,idx_dir,idx_duration,ndim(j),normalised_t);
    [SI2,~,~,SI2std,SIallother]=scale_index_other(score,idx_dir,idx_duration,ndim(j),normalised_t);
    %% Stats
    SIj(j)=SI;
    [~,pR2(j)]=ttest2(SIall,SIallother);
    subplot(3,3,9)
    errorbar(j,SI,SIstd,'.','Color',colourArea)
    hold on 
    %plot(j,SIb,'.','Color',colourArea)
    errorbar(j,SI2,SI2std,'.','Color',[0.5 0.5 0.5]) 
    %plot(j,SIb2,'^','Color',colourArea)
    box off
    ylabel('Scaling index')
    xlabel('Session number')
    if extra_plot
    if j==1
        figure(fig1)
         %% video for neuromatch
                for t=1:750
                    for i_bin=1:Nbins
                         idx=find(idx_dir==1 & idx_duration==i_bin);
                plot3(score(idx(1:t),1),score(idx(1:t),2),score(idx(1:t),3),'Color',colour_dur(i_bin,:),'LineWidth',2)
                hold on
                plot3(score(idx(1),1),score(idx(1),2),score(idx(1),3),'o','MarkerFaceColor',colour_dur(i_bin,:),'MarkerEdgeColor',colour_dur(i_bin,:))
                    end
                xlim([min(score(:,1)) max(score(:,1))])
                ylim([min(score(:,2)) max(score(:,2))])
                zlim([min(score(:,3)) max(score(:,3))])
                xlabel('PC 1')
                ylabel('PC 2')
                zlabel('PC 3')
                plot3([-0.01 0],[0.005 0.005],[-0 -0],'Color',[0.3 0.3 0.3],'LineWidth',2)
                plot3([0 0 ],[-0.005 0.005],[0 0],'Color',[0.3 0.3 0.3],'LineWidth',2)
                plot3([0 0],[0.005 0.005],[0 0.01],'Color',[0.3 0.3 0.3],'LineWidth',2)
                view([-77,29])
                axis off
                pause(0.01)
                print(['.\neuromatch_video\traj_im' num2str(t)],'-dpng','-r0')
                cla
                end
                pause
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        Ndir=max(idx_dir);
        for i_dir=1
            
            [~,score2]=pca(score(idx_dir==i_dir,1:ndim(j)));
            
            
            for i_bin=1:Nbins
                idx=find(idx_dir==i_dir & idx_duration==i_bin);
                idx2=idx_duration(idx_dir==i_dir)==i_bin;
                subplot(3,2,3)
                plot(score2(idx2,1)+1,score2(idx2,2)+1,'Color',colour_dur(i_bin,:))
                hold on
                Area_curve_bin(i_dir,i_bin)=abs(area_closed_curve(score2(idx2,1)+1,score2(idx2,2)+1));
                
                
                
%                 subplot(3,2,1)
%                
%                 
%                 plot3(score(idx(1:600),1),score(idx(1:600),2),score(idx(1:600),3),'Color',colour_dur(i_bin,:),'LineWidth',2)
%                 hold on
%                 plot3(score(idx(1),1),score(idx(1),2),score(idx(1),3),'o','MarkerFaceColor',colour_dur(i_bin,:),'MarkerEdgeColor',colour_dur(i_bin,:))
%                 xlim([min(score(:,1)) max(score(:,1))])
%                 ylim([min(score(:,2)) max(score(:,2))])
%                 zlim([min(score(:,3)) max(score(:,3))])
%                 xlabel('PC 1')
%                 ylabel('PC 2')
%                 zlabel('PC 3')
                %normalise time
                xtime2=linspace(0,1,numel(idx));
                % interpolate to match all vector lengths
                norm_score= interp1(xtime2,score(idx,1:ndim(j)),normalised_t);
                subplot(3,2,2)
                plot(normalised_t,norm_score(:,1),'Color',colour_dur(i_bin,:))
                hold on
                box off
                xlim([0 1])
                subplot(3,2,4)
                plot(normalised_t,norm_score(:,2),'Color',colour_dur(i_bin,:))
                hold on
                box off
                xlim([0 1])
                subplot(3,2,6)
                plot(normalised_t,norm_score(:,3),'Color',colour_dur(i_bin,:))
                hold on
                box off
                xlim([0 1])
            end
            
        end
        subplot(3,2,5)
        plot(1:Nbins,mean(Area_curve_bin)./mean(Area_curve_bin(:,1)),'.-','Color',colourArea)
        hold on
        
        xlabel('Normalised time')
    end
    end
    pause
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test if Hausdorff distance for different durations is shorter than for different directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(fig2)
    if do_all_dim
        fraction_dim=zeros(1,ndim(j));
        Dmean=zeros(1,ndim(j));
        Dstd=zeros(1,ndim(j));
        hdim=zeros(1,ndim(j));
        pdim=zeros(1,ndim(j));
        for i=1:ndim(j)
            [fraction_dim(i),Delta_distances]=Test_distance_between_trajectories(score,idx_dir,idx_duration,i,0);
            Dmean(i)=mean(Delta_distances);
            Dstd(i)=std(Delta_distances);
            [hdim(i),pdim(i)]=ttest(Delta_distances);
            
            clear Delta_distances
        end
        
        subplot(3,3,4)
        plot(fraction_dim)
        hold on
        xlabel('Ndim')
        ylabel('Fraction over diagonal')
        subplot(3,3,5)
        errorbar(1:ndim(j),Dmean,Dstd)
        xlabel('Ndim')
        ylabel('Dist directions- Dist duration (\Delta dist)')
        hold on
        subplot(3,3,6)
        plot(1:ndim(j),pdim)
        hold on
        plot([0 ndim(j)],[0.05 0.05],'k')
        ylabel('p-value (\Delta dist >0)')
        hold on
        box off
        clear h p
    end
    
    [fraction_above,Delta_distances,Hdist_tmp,Hdur_tmp]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim(j),1,colourArea,nsamples_condition);
    %% add a plot with the histogram for each axis
    Hdist=[Hdist;Delta_distances];
    Hdur=[Hdur;Hdur_tmp];
    %% add a plot of the Delta distances for each session to show that they are over zero. 
    [h(j),p(j)]=ttest(Delta_distances);
    subplot(3,3,4)
    plot(j,fraction_above,['o' colourArea])
    hold on
    box off
    ylabel('Fraction above')
    xlabel('Recording')
    ylim([0 1])
    
    subplot(3,3,5)
    errorbar(j,mean(Delta_distances),std(Delta_distances),'Color',colourArea)
    hold on
    box off
    ylabel('Delta distance')
    xlabel('Recording')

    
    subplot(3,3,6)
    plot(j,p(j),['o' colourArea])
    hold on
    plot([0 ndim(j)],[0.05 0.05],'k')
    ylabel('p-value')
    xlabel('Recording')
    clear Delta_distances fraction_above idx_duration idx_dir
end
%mean R^2
subplot(3,3,7)
histogram(Hdist,-0.025:0.001:0.025,'Normalization','probability')
xlabel('Distance direction')
xlim([-0.025 0.025])
box off
    
subplot(3,3,8)
histogram(Hdur,0:0.001:0.025,'Normalization','probability')
xlabel('Distance duration')
xlim([0 0.025])
box off
meanR2=mean(SIj)
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
            %identify the time when this happen
            
%             [~,idxtmp1]=min(recurrence);
%             [~,idxtmp3]=max(min(recurrence));
%             [~,idxtmp2]=min(recurrence,[],2);
%             [~,idxtmp4]=max(min(recurrence,[],2));
%             [idxtmp1(idxtmp3) idxtmp3 idxtmp2(idxtmp4) idxtmp4]
%             pause
            Hdist_dur(i_dir,counter)=max([min(recurrence),min(recurrence,[],2)']);
            nsamp(i_dir,counter)=min(nsamples_condition(i_dir,i_bin),nsamples_condition(i_dir,j_bin));
            counter=counter+1;
        end
    end
end
if do_plot
    subplot(3,3,1)
    plot(Hdist_dur')
    box off
    xlabel('Combinations')
    ylabel('H Distance')
end
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
            %             [~,idx_tmp]=min(recurrence,[],2);
            %             [i_dir j_dir i_bin]
            %             plot(score(idx1,1))
            %             hold on
            %             plot(score(idx2,1))
            %             plot(score(idx_dir==i_dir & idx_duration==i_bin+1,1))
            %             pause
            %             hold off
            %             imagesc(recurrence)
            %             pause
            
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
            subplot(3,3,2)
            
            plot(tmp,'.-','Color',colour_dir(i_dir,:))
            hold on
            
            plot(ones(size(Hdist_dur,2))*Ndir/2,Hdist_dur(i_dir,:),'.','Color',colour_dir(i_dir,:))
            
            subplot(3,3,3)
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
    subplot(3,3,3)
    %text(0.015,0.015,['Fraction Above ' num2str(total_fraction)])
    plot([0 0.02],[0 0.02],'r')
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dir),'.k');
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dur),'horizontal','.k');
    
    box off
    xlabel('Distance durations')
    ylabel('Distance direction')
    
    subplot(3,3,2)
    box off
    xlabel('Angle difference')
    ylabel('Hausdorff distance')
end

end

function A=area_closed_curve(x,y)
N=numel(x);
I=zeros(N,1);
for i=1:N-1
    I(i)=(y(i+1)+y(i))*(x(i+1)-x(i))/2;
end
I(N)=(y(1)+y(N))*(x(1)-x(N))/2;
A=sum(I);
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

function SI=scale_index_two_TR(X,Y,normalised_t)
debugging=0;
%normalise time
Xtime=linspace(0,1,size(X,1));
Ytime=linspace(0,1,size(Y,1));
% interpolate to match all vector lengths
Xnorm= interp1(Xtime,X,normalised_t);
Ynorm= interp1(Ytime,Y,normalised_t);
C=min(size(X,1),size(Y,1));

non_scaled_distance=sum(sum(abs(X(1:C,:)-Y(1:C,:))));
scaled_distance=sum(sum(abs(Xnorm-Ynorm)));

SI=(non_scaled_distance-scaled_distance)./(non_scaled_distance+scaled_distance);
%debugging
if debugging
figure
subplot(2,2,1)
plot(X(1:C,1))
hold on
plot(Y(1:C,1))
% plot(X(1:C,1),Y(1:C,1))
% hold on
% plot(Xnorm(:,1),Ynorm(:,1))

subplot(2,2,2)
plot(Xnorm(:,1))
hold on
plot(Ynorm(:,1))

subplot(2,2,3)
plot(X(1:C,2))
hold on
plot(Y(1:C,2))

subplot(2,2,4)
plot(Xnorm(:,2))
hold on
plot(Ynorm(:,2))
end
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