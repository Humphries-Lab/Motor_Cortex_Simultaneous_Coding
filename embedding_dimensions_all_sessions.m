function embedding_dimensions_all_sessions(Sessions,Area,threshold,Ndir,Nbins,t_from,t_upto,edges_dur_bin,session_N,do_plot,plot_supp)
%% embedding_dimensions_all_sessions calculate the number of embedding dimensions for all recordings denoted in Sessions and Area
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
% Ndir: number of direction to bin the movements
%
% Nbins: number of durations to bin the movements
%
% t_from: start time of the neural activity relative to movement onset for all recordings[S]
% e.g t_from=[-0.5 -0.5]
%
% t_upto: end time of the neural activity relative to movement end for all recordings [S]
% e.g t_from=[0.3 0.3]
%
% edges_dur_bin= array containing the edges of each duration bin [S]
%
% session_N = Number of the behavioural session associated to Session and
% Area. There should be as many behavioural sessions as unique names in the
% array Session. If a Session contains simultaneous recordings from M1 and PMd, they
% belong to the same behavioural session. e.g
% session={'MC_S1_raw.mat','MC_S2_raw.mat','MM_S1_raw.mat','MM_S1_raw.mat'};
% Area={'M1','M1','M1','PMd'};
% session_N=[1 2 3 3]
%
% do_plot= 1- plots the average FR of the population for each duration bin
%             and the PCA trajectories for the first duration bin. This
%             generates one figure per recording
%
%          0- omit the above plot
%
% plot_supp= 1- plots 1) the variance explained and the number of embedding
%               dimensions for all recordings and 2) the correlation
%               between  the parameters of the movements (duration, distance, max speed)
%
%            0- omit the above plot
%
%
% 24/01/2023
% Andrea Colins Rodriguez

if plot_supp
    fig_sup_1=figure;
    fig_sup_2=figure;
end

embedding_dim=zeros(size(Sessions));
R(max(session_N))=struct();
for isession=1:size(Sessions,2)
    disp(['Starting PCA for recording ' num2str(isession)])
    
    [variance,dist_mov_dir,mov_duration,max_speed]=embedding_dimensions(Sessions{isession}, Area{isession},Ndir,Nbins,t_from(isession),t_upto(isession),edges_dur_bin,do_plot);
    
    %% kinematics
    if isession==1 && plot_supp
        figure(fig_sup_2)
        plot_kinematics(dist_mov_dir,mov_duration,max_speed,Nbins)
        
    end
    
    R(session_N(isession)).R_dist_dur=corr([dist_mov_dir{1:Nbins}]',[mov_duration{1:Nbins}]');
    R(session_N(isession)).R_dist_speed=corr([dist_mov_dir{1:Nbins}]',[max_speed{1:Nbins}]');
    R(session_N(isession)).R_speed_dur=corr([max_speed{1:Nbins}]',[mov_duration{1:Nbins}]');
    
    % Dimensionality of the subspace
    embedding_dim(isession)=find(cumsum(variance)>threshold,1,'First');
    
    if plot_supp
        if strcmp(Area{isession},'M1')
            colourArea=[85 30 116]./256;
        else
            colourArea=[89 156 153]./256;
        end
        
        figure(fig_sup_1)
        subplot(1,3,1)
        plot(cumsum(variance),'Color',colourArea)
        hold on
        
        subplot(1,3,2)
        hold on
        plot(isession,embedding_dim(isession),'o','Color',colourArea)


        
        subplot(1,3,3)
        hold on
        plot(isession,numel(variance),'o','Color',colourArea)
    end
    clear variance dist_mov_dir mov_duration max_speed
    
end
if plot_supp
    figure(fig_sup_1)
    subplot(1,3,1)
    plot([0 100],[threshold threshold])
    xlabel('N dimensions')
    ylabel('Variance explained [%]')
    box off
    
    subplot(1,3,2)
    hold on
    box off
    xlabel('N recording')
    ylabel('Embedding dimensions')
    ylim([0 13])
    legend(Area)
    
    subplot(1,3,3)
    hold on
    box off
    xlabel('N recording')
    ylabel('Number of neurons')
    ylim([0 100])
    legend(Area)
    %%
    figure(fig_sup_2)
    
    subplot(2,3,4)
    histogram([R.R_dist_dur])
    box off
    xlabel('Corr Distance-Duration')
    ylabel('Number of sessions')
    xlim([0 1])
    title(['Mean = ' num2str(mean([R.R_dist_dur]),'%.2f')])
    
    subplot(2,3,5)
    histogram([R.R_dist_speed])
    box off
    xlabel('Corr Distance-Speed')
    xlim([0 1])
    title(['Mean = ' num2str(mean([R.R_dist_speed]),'%.2f')])
    
    subplot(2,3,6)
    histogram([R.R_speed_dur])
    box off
    xlabel('Corr Duration-Speed')
    xlim([0 1])
    title(['Mean = ' num2str(mean([R.R_speed_dur]),'%.2f')])
end
end