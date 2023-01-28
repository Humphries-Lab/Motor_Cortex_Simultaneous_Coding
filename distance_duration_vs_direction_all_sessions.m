function p=distance_duration_vs_direction_all_sessions(Sessions,Areas,threshold,Nbins,do_plot_supp)
%% distance_duration_vs_direction_all_sessions compares the distance between
%% trajectories corresponding to different movement durations (same
%% direction) and movemements of adjacent bin directions (same duration)
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
% Nbins: number of durations to bin the movements
%
% do_plot_supp: 1- plots example trajectories of different durations and
% same direction for 1 recording 
%               0- omits the plot above
% 
% OUTPUTS
% p: p-values of the Delta distance (distances between trajectories of
% adjacent direction bins- distance between trajectories of same duration)
% for all recordings 
% 24/05/2023
% Andrea Colins Rodriguez


ndim=zeros(numel(Sessions),1);
colour_dur=plasma(Nbins);
final_length=600;
normalised_t=linspace(0,1,final_length);
Nsessions=numel(Sessions);
p=zeros(1,Nsessions);
pR2=zeros(1,Nsessions);
fig1=figure;
fig2=figure;
Hdist=[];
Hdur=[];
for isession=1:Nsessions
    
    Area=Areas{isession};
    session=Sessions{isession};
    
    load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance')
    
    ndim(isession)=find(cumsum(variance)>threshold,1,'First');
    
    %% Plot example trajectories for same direction different durations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_plot_supp && isession==1
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
    
    [R2_same,SIstd,SI_same]=R_same(score,idx_dir,idx_duration,ndim(isession),normalised_t);
    [R2_other_dir,SI2std,SIallother]=R2_other(score,idx_dir,idx_duration,ndim(isession),normalised_t);

    [~,pR2(isession)]=ttest2(SI_same,SIallother);
    
    figure(fig2)
    subplot(2,3,1)
    errorbar(isession,R2_same,SIstd,'.','Color',colourArea)
    hold on
    errorbar(isession,R2_other_dir,SI2std,'.','Color',[0.5 0.5 0.5])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test if Hausdorff distance for different durations is shorter than for different directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [~,Delta_distances,Hdist_tmp,Hdur_tmp]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim(isession),1,colourArea);
    %% add a plot with the histogram for each axis
    
    Hdist=[Hdist;Hdist_tmp];
    Hdur=[Hdur;Hdur_tmp];
    
    
   [~,p(isession)]=ttest(Delta_distances);
    
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

text(1,0.9,['mean = ' num2str(mean(R2_same),2)],'FontSize',8)
text(1,-1,['mean = ' num2str(mean(R2_other_dir),2)],'FontSize',8)
text(1, 0, ['Max p-value = ' num2str(max(pR2),2)],'FontSize',8)


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

function [total_fraction,Delta_distances,All_dist_dir,All_dist_dur]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim,do_plot,ColourArea)
%% Test_distance_between_trajectories calculates the Hausdorff distance between
%% trajectories of different durations and different directions
%
% INPUTS
%
% score: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score
% matrix
% 
% idx_duration: array containing the duration bin of each row in the score
% matrix
%
% ndim: number of dimensions of the trajectories
%
% do_plot: 1- plots distance between trajectories of different durations vs
%             distance between trajectories from adjacent direction bins (same
%             duration).
%          0- omits the plot above
%
% ColourArea: Colour of the dots for the plot above. 
% 
% OUTPUTS
%
% total_fraction: fraction of trajectories from same duration that are
% closer than trajectories of adjacent bins
%
% Delta_distances: difference between distance between trajectories of adjacent direction bins and
% trajectories of same duration.
%
% All_dist_dir: Distance between trajectories of same duration and
% adjacent direction bins
%
% All_dist_dur:Distance between trajectories of same direction and
% different duration
%
% 27/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
Nbins=max(idx_duration);
Hdist_dir=zeros(Ndir,Ndir,Nbins);
Hdist_dur=nan(Ndir,factorial(Nbins-1));
%nsamp=nan(Ndir,factorial(Nbins-1));
%idx_others=1:Ndir;
%idx_others(round(Ndir/2))=[];

%% calculate distance between trajectories with same direction but different durations
for i_dir=1:Ndir
    counter=1;
    
    for i_bin=1:Nbins-1
        idx1=idx_dir==i_dir & idx_duration==i_bin;
        for j_bin=i_bin+1:Nbins
            idx2=idx_dir==i_dir & idx_duration==j_bin;
            %% distance between trajectories
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            
            Hdist_dur(i_dir,counter)=max([min(recurrence),min(recurrence,[],2)']);
            %nsamp(i_dir,counter)=min(nsamples_condition(i_dir,i_bin),nsamples_condition(i_dir,j_bin));
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

%% Calculate distance between trajectories with duration but different directions
Delta_distances=[];
counter_points_above=[];
All_dist_dur=[];
All_dist_dir=[];

for i_bin=1:Nbins
    
    for i_dir=1:Ndir
        idx1=(idx_dir==i_dir & idx_duration==i_bin);
        for j_dir=1:Ndir
            idx2=(idx_dir==j_dir & idx_duration==i_bin);
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            Hdist_dir(i_dir,j_dir,i_bin)=max([min(recurrence),min(recurrence,[],2)']);
        end
        
        tmp=circshift(Hdist_dir(i_dir,:,i_bin),round(Ndir/2)-i_dir);
        
        counter_points_above=[counter_points_above, tmp(Ndir/2-1)>Hdist_dur(i_dir,:) ,tmp(Ndir/2+1)>Hdist_dur(i_dir,:)];
        Delta_distances=[Delta_distances,tmp(Ndir/2-1)-Hdist_dur(i_dir,:),tmp(Ndir/2+1)-Hdist_dur(i_dir,:)];
        All_dist_dur=[All_dist_dur;Hdist_dur(i_dir,:)';Hdist_dur(i_dir,:)'];
        All_dist_dir=[All_dist_dir;ones(size(Hdist_dur,2),1)*tmp(Ndir/2-1);ones(size(Hdist_dur,2),1)*tmp(Ndir/2+1)];
        
        if do_plot
            
            subplot(2,3,2)
            hold on
            
            transparency=nsamp./max(nsamp(:));
            for iplot=1:size(Hdist_dur,2)
                s=scatter(Hdist_dur(i_dir,iplot),tmp(Ndir/2-1),8,ColourArea,'filled');
                alpha(s,transparency(i_dir,iplot))
                s=scatter(Hdist_dur(i_dir,iplot),tmp(Ndir/2+1),8,ColourArea,'filled');
                alpha(s,transparency(i_dir,iplot))
                
            end
        end
    end
    
end

total_fraction=sum(counter_points_above)/sum(~isnan(counter_points_above(:)));

if do_plot
    subplot(2,3,2)
    plot([0 0.02],[0 0.02],'r')
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dir),'.k');
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dur),'horizontal','.k');
    
    box off
    xlabel('Distance durations')
    ylabel('Distance direction')
    
end

end


function [R2_same_dir,SIstd,SIall]=R_same(score,idx_dir,idx_duration,ndim,normalised_t)
%% R2_same calculates the coefficient of determination between trajectories
%% of same directions
%
% INPUTS
%
% score: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score
% matrix
% 
% idx_duration: array containing the duration bin of each row in the score
% matrix
%
% ndim: number of dimensions of the trajectories
%
% normalised_t: vector of timebins to normalise trajectories 
% 
% OUTPUTS
%
% R2_other_dir: mean coefficient of determination (R2) across all trajectories
%
% SIstd: standard error of the mean R2 across all trajectories
%
% SIall: array of coefficient of determination (R2) across all trajectories
%
% 24/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
Nbins=max(idx_duration);
R2_same_dir=zeros(Ndir,Nbins-1);

for i_dir=1:Ndir
    
    counter=1;
    for i_bin=1:Nbins
        idx=idx_dir==i_dir & idx_duration==i_bin;
        for j_bin=i_bin+1:max(idx_duration)
            idx_j=idx_dir==i_dir & idx_duration==j_bin;
            R2_same_dir(i_dir,counter)=scaled_TR_R2(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
               counter=counter+1;
        end
    end
end

SIall=R2_same_dir(:);
SIstd=std(R2_same_dir(:))/sqrt(numel(R2_same_dir));
R2_same_dir=nanmean(R2_same_dir(:));

end

function [R2_other_dir,SIstd,SIall]=R2_other(score,idx_dir,idx_duration,ndim,normalised_t)
%% R2_other calculates the coefficient of determination between trajectories
%% of different directions
%
% INPUTS
%
% score: Projection of the neural activity into the subspace. Rows are
% samples, columns are neurons
% 
% idx_dir: array containing the direction bin of each row in the score
% matrix
% 
% idx_duration: array containing the duration bin of each row in the score
% matrix
%
% ndim: number of dimensions of the trajectories
%
% normalised_t: vector of timebins to normalise trajectories 
% 
% OUTPUTS
%
% R2_other_dir: mean coefficient of determination (R2) across all trajectories
%
% SIstd: standard error of the mean R2 across all trajectories
%
% SIall: array of coefficient of determination (R2) across all trajectories
%
% 24/05/2023
% Andrea Colins Rodriguez


Ndir=max(idx_dir);
Nbins=max(idx_duration);
R2_other_dir=zeros(Ndir,Nbins-1);

for i_dir=1:Ndir
    counter=1;
    for j_dir=i_dir+1:Ndir
        
        for i_bin=1:Nbins
            idx=idx_dir==i_dir & idx_duration==i_bin;
            
            for j_bin=i_bin+1:Nbins
                idx_j=idx_dir==j_dir & idx_duration==j_bin;
                
                R2_other_dir(i_dir,counter)=scaled_TR_R2(score(idx,1:ndim),score(idx_j,1:ndim),normalised_t);
                counter=counter+1;
            end
        end
    end
end

SIall=R2_other_dir(:);
SIstd=std(R2_other_dir(:))/sqrt(numel(R2_other_dir));
R2_other_dir=nanmean(R2_other_dir(:));
end


function R2_scaled=scaled_TR_R2(X,Y,normalised_t)
%% scaled_TR_R2 scales trajectories X and Y to a normalised_t lenght and
%% computes their coefficient of determination
% 
% INPUTS
%
% X: n-dimensional trajectory. Rows are times, columns are dimensions 
%
% Y: n-dimensional trajectory. Rows are times, columns are dimensions 
%
% normalised_t: vector of timebins to normalise trajectories 
%
% OUTPUTS
%
% SI: Coefficient of determination between the scaled trajectories.
%
% 24/05/2023
% Andrea Colins Rodriguez

%normalise time
Xtime=linspace(0,1,size(X,1));
Ytime=linspace(0,1,size(Y,1));
% interpolate to match all vector lengths
Xnorm= interp1(Xtime,X,normalised_t);
Ynorm= interp1(Ytime,Y,normalised_t);

R2_scaled=R2(Xnorm,Ynorm);
end
function r=R2(X,Y)
%% R2 Computes the coefficient of determination between n-dimensional
%% trajectories X and Y
X2=reshape(X,numel(X),1);
Y2=reshape(Y,numel(Y),1);
ymean=mean(Y2);
Sres=sum((Y2-X2).^2);
Stot=sum((Y2-ymean).^2);
r=1-Sres/Stot;
end