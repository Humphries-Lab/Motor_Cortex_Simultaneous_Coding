function [total_fraction,Delta_distances,All_dist_dir,All_dist_dur,Delta_dist_all_same]=Test_distance_between_trajectories(score,idx_dir,idx_duration,ndim,do_plot,ColourArea,nsamples_condition,Hdist_same)
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
% nsamples_condition: number of movements selected for each specific
% condition
% 
% OUTPUTS
%
% total_fraction: fraction of trajectories from the same duration that are
% closer than trajectories of adjacent bins
%
% Delta_distances: difference between the distance between trajectories of adjacent direction bins and
% trajectories of the same duration.
%
% All_dist_dir: Distance between trajectories of the same duration and
% adjacent direction bins
%
% All_dist_dur: Distance between trajectories of the same direction and
% different duration
%
% 27/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
Nbins=max(idx_duration);
Hdist_dir=zeros(Ndir,Nbins,Ndir);
Hdist_dur=nan(Nbins,Nbins,Ndir);
nsamp=nan(Nbins,Nbins,Ndir);


%% calculate distance between trajectories with same direction but different durations
for i_dir=1:Ndir
    counter=1;
    
    for i_bin=1:Nbins-1
        idx1=idx_dir==i_dir & idx_duration==i_bin;
        for j_bin=i_bin+1:Nbins
            idx2=idx_dir==i_dir & idx_duration==j_bin;
            %% distance between trajectories
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            
            Hdist_dur(i_bin,j_bin,i_dir)=max([min(recurrence),min(recurrence,[],2)']);
            Hdist_dur(j_bin,i_bin,i_dir)=Hdist_dur(i_bin,j_bin,i_dir);
            
            nsamp(i_bin,j_bin,i_dir)=min(nsamples_condition(i_dir,i_bin),nsamples_condition(i_dir,j_bin));
            nsamp(j_bin,i_bin,i_dir)=nsamp(i_bin,j_bin,i_dir);
            counter=counter+1;
        end
    end
end


%% Calculate distance between trajectories with same duration but different directions
Delta_dist_all_same=[];
Ncomp_bins=Nbins-1;
counter_points_above=nan(Ncomp_bins*Ndir*Nbins*2,1);

Delta_distances=nan(Ncomp_bins*Ndir*Nbins*2,1);
All_dist_dur=nan(Ncomp_bins*Ndir*Nbins*2,1);
All_dist_dir=nan(Ncomp_bins*Ndir*Nbins*2,1);

max_transparency=max(nsamp,[],'all','omitnan');
icomp=1;
for i_bin=1:Nbins
    
    for i_dir=1:Ndir
        idx1=(idx_dir==i_dir & idx_duration==i_bin);
        for j_dir=1:Ndir
            idx2=(idx_dir==j_dir & idx_duration==i_bin);
            recurrence=pdist2(score(idx1,1:ndim),score(idx2,1:ndim));
            Hdist_dir(i_dir,j_dir,i_bin)=max([min(recurrence),min(recurrence,[],2)']);
        end
        
        tmp=circshift(Hdist_dir(i_dir,:,i_bin),round(Ndir/2)-i_dir);
        
        dist_bin=Hdist_dur(i_bin,:,i_dir);
        dist_bin=dist_bin(~isnan(dist_bin));
        
        counter_points_above((icomp-1)*Ncomp_bins*2+1:icomp*Ncomp_bins*2)=[tmp(Ndir/2-1)>dist_bin,tmp(Ndir/2+1)>dist_bin]';
        
        Delta_distances((icomp-1)*Ncomp_bins*2+1:icomp*Ncomp_bins*2)=[tmp(Ndir/2-1)-dist_bin,tmp(Ndir/2+1)-dist_bin]';
        
        All_dist_dur((icomp-1)*Ncomp_bins*2+1:icomp*Ncomp_bins*2)=[dist_bin,dist_bin]';
        All_dist_dir((icomp-1)*Ncomp_bins*2+1:icomp*Ncomp_bins*2)=[ones(Ncomp_bins,1)*tmp(Ndir/2-1);ones(Ncomp_bins,1)*tmp(Ndir/2+1)];
        
        icomp=icomp+1;
        
        if do_plot
            if ~isnan(Hdist_same(i_bin,i_dir))
            subplot(2,3,2)
            hold on
            scatter(dist_bin,Hdist_same(i_bin,i_dir)*ones(1,size(dist_bin,2)),8,ColourArea,'filled')
            Delta_dist_all_same=[Delta_dist_all_same,dist_bin-Hdist_same(i_bin,i_dir)*ones(1,size(dist_bin,2))];
            end
            
            subplot(2,3,1)
            hold on
      
            transparency=nsamp(i_bin,:,i_dir)./max_transparency;
            transparency=transparency(~isnan(transparency));
            
            for iplot=1:size(dist_bin,2)
                s=scatter(dist_bin(iplot),tmp(Ndir/2-1),8,ColourArea,'filled');
                alpha(s,transparency(iplot))
                s=scatter(dist_bin(iplot),tmp(Ndir/2+1),8,ColourArea,'filled');
                alpha(s,transparency(iplot))
                
            end
        end
    end
    
end

total_fraction=sum(counter_points_above)/sum(~isnan(counter_points_above(:)));

if do_plot
    subplot(2,3,1)
    plot([0 0.02],[0 0.02],'r')
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dir),'.k');
    errorbar(mean(All_dist_dur),mean(All_dist_dir),std(All_dist_dur),'horizontal','.k');
    
    box off
    xlabel('Distance durations')
    ylabel('Distance direction')
    
end

end