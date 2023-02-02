function recurrence_region_all_sessions(Sessions,Areas,threshold,Ndir,Nbins,threshold_dist)
%% recurrence_region_all_sessions defines a region of recurrence 
%% for each recording and counts how many trajectories are inside this region at each moment in time
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
% threshold_dist: threshold of the percentil of distance(theta) for a point to be
% considered recurrent. threshold_dist ranges from 0 to 100.
%
% 26/01/2023
% Andrea Colins Rodriguez

ref_bin=2; % duration bin used as reference to define the region of recurrence
upto=200;% number of points to analyse from the beginning of prep

fig1=figure;
colormap(flipud(colormap('gray')))
colour_dur=plasma(Nbins);
colour_dir=hsv(Ndir);

not_coming_back=zeros(numel(Areas),1);
not_starting_there=zeros(numel(Areas),1);

for isession=1:numel(Areas)
    Area=Areas{isession};
    session=Sessions{isession};
    End_recurrence=true(Nbins*Ndir,1);
    Start_recurrence=true(Nbins*Ndir,1);
    
    for i_bin=[ref_bin:Nbins 1:ref_bin-1]
        
        
        load(['../Output_files/PCA_' session(1:end-4) '_' Area '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto')
        
        ndim=find(cumsum(variance)>threshold,1,'First');
        
        %% Define region of recurrence for the recording
        if i_bin==ref_bin

            [limit,center]=define_recurrence_region(score(idx_duration==i_bin,1:ndim),idx_dir(idx_duration==i_bin),threshold_dist);
            
            if isession==1
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
                
                for i_dir=1:Ndir
                    subplot(2,2,1)
                    plot3(score(idx_dir==i_dir & idx_duration==i_bin,1),score(idx_dir==i_dir & idx_duration==i_bin,2),score(idx_dir==i_dir & idx_duration==i_bin,3),'Color',colour_dir(i_dir,:),'LineWidth',1)
                    hold on
                    
                    subplot(2,2,3)
                    plot3(score(idx_dir==i_dir & idx_duration==i_bin,2),score(idx_dir==i_dir & idx_duration==i_bin,3),score(idx_dir==i_dir & idx_duration==i_bin,4),'Color',colour_dir(i_dir,:),'LineWidth',1)
                    hold on
                end
            end
        end
        
        %% compute recurrence time
        segment_length=round((t_upto(i_bin)-t_from)*1000);
        
        Is_in_Rec=zeros(segment_length,Ndir);
        
        for i_dir=1:Ndir
            % Check if points of traj(i_dir,i_bin) are in the region of
            % recurrence
            Is_in_Rec(:,i_dir)=IsInRecRegion(score(idx_dir==i_dir & idx_duration==i_bin,1:ndim),[limit,center]);
            
            %Check if the points at the start and end are in the region of
            %recurrence
            Start_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Is_in_Rec(1:upto,i_dir),1)>0;
            End_recurrence(Ndir*(i_bin-1)+i_dir)=sum(Is_in_Rec(end-upto:end,i_dir),1)>0;
        end
        
        figure(fig1)
        if strcmp(Areas{isession},'M1')
            subplot(2,2,2)
            hold on
        else
            subplot(2,2,4)
            hold on
        end
        
        plot((1:segment_length)+t_from*1000,sum(Is_in_Rec,2),'Color',colour_dur(i_bin,:))
        
        
    end
    %% Compute number of points at the beginning and end.
    
    not_coming_back(isession,1)=sum(End_recurrence<1)/numel(End_recurrence);
    not_starting_there(isession,1)=sum(Start_recurrence<1)/numel(End_recurrence);
    
    clear Start_recurrence End_recurrence
end

figure(fig1)
subplot(2,2,2)
title('M1')
box off
subplot(2,2,4)
title('PMd')
xlabel('Time to movement onset [ms]')
ylabel('Number of trajectories in the recurrence region')
box off

%[mean(1-not_starting_there(M1)) mean(1-not_starting_there(~M1))]
%[mean(1-not_coming_back(M1)) mean(1-not_coming_back(~M1))]

end


function [radius,center]=define_recurrence_region(score,idx_dir,threshold_dist)
%% define_recurrence_region defines a region of recurrence as a hypersphere around the startpoints of the trajectories
%
% INPUTS
% score: neural trajectories for the selected reference duration 
% 
% idx_dir: array indicating the direction bin corresponding to each point
% in score
%
% threshold_dist: threshold of the percentil of distance(theta) for a point to be
% considered recurrent. threshold_dist ranges from 0 to 100.
%
% OUTPUTS
% 
% radius: radius of the recurrence region
% 
% center: center of the recurrence region
%
% 26/01/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);

rec_length=1;

all_points_beginning=zeros(rec_length*Ndir,size(score,2));

radius=prctile(pdist(score),threshold_dist);

for i_dir=1:Ndir
    all_points_beginning((i_dir-1)*rec_length+1:i_dir*rec_length,:)=score(find(idx_dir==i_dir,rec_length,'First'),:);
end

center=mean(all_points_beginning);
end


function Its=IsInRecRegion(score,region)
limit=region(1);
center=region(2:end);

Its=pdist2(center,score)<limit;
end