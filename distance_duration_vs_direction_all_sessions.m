function distance_duration_vs_direction_all_sessions(Sessions,Areas,threshold,Ndir,Nbins,do_plot_supp)
%% distance_duration_vs_direction_all_sessions compares the distance between
%% trajectories corresponding to different movement durations (same
%% direction) and movements of adjacent bin directions (same duration)
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
% Nbins: number of durations to bin the movements
%
% do_plot_supp: 1- plots example trajectories of different durations and
% same direction for 1 recording 
%               0- omits the plot above
%
% 24/05/2023
% Andrea Colins Rodriguez


ndim=zeros(numel(Sessions),1);
colour_dur=plasma(Nbins);
final_length=600;
normalised_t=linspace(0,1,final_length);
Nsessions=numel(Sessions);
p=zeros(1,Nsessions);
pR_noise=nan(Nsessions,1);
p_delta_distance_same=nan(Nsessions,1);
pR2=zeros(1,Nsessions);
R2_same=zeros(1,Nsessions);
R2_other_dir=zeros(1,Nsessions);
mean_r=nan(Nsessions,1);
fraction_explained=zeros(1,Nsessions);
Delta_all=[];

if do_plot_supp
fig1=figure;
end

fig2=figure;
Ncomb_total=Nbins*(Nbins-1)*2*Ndir;
Hdist=nan(Nsessions*Ncomb_total,1);
Hdur=nan(Nsessions*Ncomb_total,1);

for isession=1:Nsessions
    
    Area=Areas{isession};
    session=Sessions{isession};
    
    load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','nsamples_condition','r','Hdist_same')
    %Hdist_same=nan(Nbins,Ndir);
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
    
    [R2_same(isession),SIstd,SI_same]=R_same(score,idx_dir,idx_duration,ndim(isession),normalised_t);
    [R2_other_dir(isession),SI2std,SIallother]=R2_other(score,idx_dir,idx_duration,ndim(isession),normalised_t);

    [~,pR2(isession)]=ttest2(SI_same,SIallother);
    [~,pR_noise(isession)]=ttest2(SI_same,r);
    
    figure(fig2)
    subplot(2,3,6)
    errorbar(isession,R2_same(isession),SIstd,'.','Color',colourArea)
    hold on
    errorbar(isession,R2_other_dir(isession),SI2std,'.','Color',[0.5 0.5 0.5])
    errorbar(isession,mean(r),std(r)/sqrt(numel(r)),'Color','k')
    mean_r(isession)=mean(r);
    fraction_explained(isession)=(R2_same(isession)-R2_other_dir(isession))/(mean(r)-R2_other_dir(isession));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Test if Hausdorff distance for different durations is shorter than for different directions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [~,Delta_distances,Hdist_tmp,Hdur_tmp,Delta_dist_all_same]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim(isession),1,colourArea,nsamples_condition,Hdist_same);
    %% add a plot with the histogram for each axis
    
    Hdist((isession-1)*Ncomb_total+1:isession*Ncomb_total)=Hdist_tmp;
    Hdur((isession-1)*Ncomb_total+1:isession*Ncomb_total)=Hdur_tmp;
    Delta_all=[Delta_all;Delta_dist_all_same];
    
   [~,p(isession)]=ttest(Delta_distances);
   [~,p_delta_distance_same(isession)]=ttest(Delta_dist_all_same,0,'Tail','right');

    subplot(2,3,4)
    errorbar(isession,mean(Delta_distances),std(Delta_distances),'.','Color',colourArea)
    hold on
    
       subplot(2,3,5)
       hold on
    errorbar(isession,mean(Delta_dist_all_same),std(Delta_dist_all_same),'.','Color',colourArea)
    
    clear Delta_distances fraction_above idx_duration idx_dir
end

subplot(2,3,6)
plot([0 Nsessions],[0 0],'Color',[0.5 0.5 0.5])
box off
ylabel('R^2')
xlabel('Recording number')
text(1,0.9,['mean = ' num2str(mean(R2_same),2)],'FontSize',8)
text(1,-1,['mean = ' num2str(mean(R2_other_dir),2)],'FontSize',8)
text(1, 0, ['Max p-value = ' num2str(max(pR2),2)],'FontSize',8)
text(6,0.8,['mean = ' num2str(mean(mean_r,'omitnan'),2)],'FontSize',8)
text(7,-1,['% explained = ' num2str(mean(fraction_explained),2)],'FontSize',8)

subplot(2,3,4)
plot([0 Nsessions],[0 0],'Color',[0.5 0.5 0.5])
box off
ylabel('\Delta Distance')
xlabel('Recording')
text(1,-0.001,['max p-value = ' num2str(max(p),'%.e')])
ylim([-0.005 0.025])

subplot(2,3,2)
xlabel('\Delta Duration')
ylabel('\Delta within bin')
xlim([0 0.025])
ylim([0 0.025])
diag_hist(Delta_all,0.025)
axis square

subplot(2,3,1)
xlim([0 0.025])
ylim([0 0.025])
axis square

diag_hist(Hdur-Hdist,0.025)
%axes('Position',[0.15,0.58,0.2,0.2]);
%h=histogram(Hdur,-2*a:0.001:a,'Normalization','probability');

%xlim([-a a])
%ylim([0 0.3])

%view(45,90)
%axis square
%box off

disp('----------------------------------')

pvalues=table(Sessions',Areas',pR_noise)

% axes('Position',[0.3,0.58,0.028,0.338]);
% histogram(Hdist,0:0.001:0.025,'Normalization','probability','orientation','horizontal')
% ylim([0 0.025])
% box off

subplot(2,3,2)
xlim([0 0.025])
ylim([0 0.025])
plot([0 0.025],[0 0.025],'r')
box off
%diag_hist(Hdur-Hdist,0.025)

subplot(2,3,5)
plot([0 13],[0 0],'Color',[0.5 0.5 0.5])
ylim([0 0.025])
ylim([-0.005 0.025])

pvalues_control=table(Sessions',Areas',p_delta_distance_same)
disp(['Average p-value delta control =', num2str(mean(p_delta_distance_same))])
sum(p_delta_distance_same<0.05)
end