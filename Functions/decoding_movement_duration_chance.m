function error_by_dir_shift=decoding_movement_duration_chance(score,idx_dir,idx_duration,ref_bin,Ndim)
%% decoding_movement_duration_chance predicts the duration of trajectories in score using as reference the trajectory of duration bin ref_bin
% The output of this function is a matrix of the prediction error when the
% reference and sample trajectories are different. 
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
%
% OUTPUTS
%
% error_by_dir_shift: Matrix containing the prediction error for all test
% trajectories [ms]. Rows are directions shifted so that the middle row
% corresponds to ~0 angle difference and the first and end rows are the most
% distant direction bins (~-+180^o). Columns are the number of comparisons
% made (Ntest_bins*Ndir).
%
%
% 28/05/2023
% Andrea Colins Rodriguez

Ndir=max(idx_dir);
test_bins=1:max(idx_duration);
test_bins(ref_bin)=[];
ref_duration=sum((idx_dir==1)& (idx_duration==ref_bin));
Ntest_bins=numel(test_bins);
traj_length=nan(1,Ntest_bins);

for i_bin=1:Ntest_bins
    bin_i=test_bins(i_bin);
    traj_length(i_bin)=sum((idx_dir==1)& (idx_duration==bin_i));
end

fraction=zeros(Ndir,Ntest_bins);
estimated_duration=zeros(Ndir,Ntest_bins);
error_by_dir_shift=nan(Ndir,Ntest_bins*Ndir);
counter=1;
for i_dir=1:Ndir
    idx=(idx_dir==i_dir)& (idx_duration==ref_bin);
    %reference
    ref_traj=score(idx,1:Ndim);
    clear idx
    
    for i_bin=1:Ntest_bins
        bin_i=test_bins(i_bin);
        
        for j_dir=1:Ndir
            
            idx=(idx_dir==j_dir) & (idx_duration==bin_i);
            
            sample_traj=score(idx,1:Ndim);
            fraction(i_dir,i_bin,j_dir)=compute_fraction_speed(ref_traj,sample_traj);
            estimated_duration(i_dir,i_bin,j_dir)=ref_duration/fraction(i_dir,i_bin,j_dir);
            clear idx
        end

        error_by_dir=abs(squeeze(estimated_duration(i_dir,i_bin,:))-traj_length(i_bin));
        error_by_dir_shift(:,counter)=circshift(error_by_dir,round(Ndir/2)-i_dir);
        counter=counter+1;
        
    end
    
    
end

end