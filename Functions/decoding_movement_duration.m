function [median_total_error,std_error,estimated_duration]=decoding_movement_duration(score,idx_dir,idx_duration,ref_bin,Ndim,Nbins,do_plot)
%% decoding_movement_duration predicts the duration of trajectories
%% corresponding to the same movement direction than the sample trajectory
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
% ref_bin: number of the duration bin to be used as reference
%
% Ndim: number of dimensions of the trajectories
%
% Nbins: number of durations bins 
%
% do_plot: 1- plots the average (across directions of the) time indices of the test trajectories vs the time
%             indices of the reference trajectories
%          0- omits the plot above
%
% OUTPUTS
% 
% median_total_error: Median error in the duration prediction across
% all predicted durations
%
% std_error: standard deviation of the error in the duration prediction across
% all predicted durations
%
% estimated_duration: predicted duration for all test trajectories [ms]
%
% 28/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
colour_dir=hsv(Ndir);
colour_dur=plasma(Nbins);


test_bins=1:Nbins;
test_bins(ref_bin)=[];
Ntest_bins=numel(test_bins);
traj_length=nan(1,Ntest_bins);
fraction=zeros(Ndir,Ntest_bins);
estimated_duration=zeros(Ndir,Ntest_bins);
idx_tmp=cell(3,1);
error_duration=zeros(Ndir,Ntest_bins);



ref_duration=sum((idx_dir==1)& (idx_duration==ref_bin));

for i_bin=1:Ntest_bins
    bin_i=test_bins(i_bin);
    traj_length(i_bin)=sum((idx_dir==1)& (idx_duration==bin_i));
end


for i_dir=1:Ndir
    idx=((idx_dir==i_dir)& (idx_duration==ref_bin));
    %reference
    ref_traj=score(idx,1:Ndim);
    clear idx
    
    for i_bin=1:Ntest_bins
        bin_i=test_bins(i_bin);
        idx=((idx_dir==i_dir) & (idx_duration==bin_i));
        test_traj=score(idx,1:Ndim);
        %% test with random data
        %test_traj=rand(size(test_traj,1),Ndim);
        [fraction(i_dir,i_bin),idx_tmp{i_bin}(:,i_dir)]=compute_fraction_speed(ref_traj,test_traj);
        estimated_duration(i_dir,i_bin)=ref_duration/fraction(i_dir,i_bin);
        clear idx
    end
    error_duration(i_dir,:)=abs(estimated_duration(i_dir,:)-traj_length);
    
    
    subplot(2,3,4)
    plot(traj_length,estimated_duration(i_dir,:),'.','Color',colour_dir(i_dir,:))
    hold on
    
    subplot(2,3,5)
    plot(traj_length,100*abs(estimated_duration(i_dir,:)-traj_length)./traj_length,'.','Color',colour_dir(i_dir,:))
    hold on
    
end


std_error=std(error_duration(:));
median_total_error=median(error_duration(:));

if do_plot
    
    subplot(2,3,3)
    for i_bin=1:Ntest_bins
        y = mean(idx_tmp{i_bin},2)'; % your mean vector
        x = 1:numel(y);
        std_dev =std(idx_tmp{i_bin},[],2)';% std across directions
        curve1 = y + std_dev;
        curve2 = y - std_dev;
        test_traj = [x, fliplr(x)];
        inBetween = [curve1, fliplr(curve2)];
        fill(test_traj, inBetween, colour_dur(test_bins(i_bin),:))
        hold on
        plot(x, y, 'LineWidth', 2,'Color',colour_dur(test_bins(i_bin),:));
        plot([0 traj_length(i_bin)],[0 ref_duration],'k--')
    end
    plot([0 ref_duration],[0 ref_duration],'k--')
    xlim([0 traj_length(end)+100])
    box off
    xlabel('Time test trajectory')
    ylabel('Time ref trajectory')
end
end


