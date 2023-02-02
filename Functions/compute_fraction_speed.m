function [slope,idx_test]=compute_fraction_speed(ref_traj,sample_traj)
%% compute_fraction_speed computes the ratio between the speeds of the reference and the sample trajectory 
%
% INPUTS 
%
% ref_traj: reference trajectory. Rows are time bins, columns are
% dimensions.
% 
% sample_traj: sample trajectory. Rows are time bins, columns are
% dimensions.
% 
% OUTPUTS
%
% slope: ratio between speed of the sample and reference trajectories.
% slope>1 => sample is faster that the reference trajectory. 
%
% idx_test: idx_test(i) is the bin in the reference trajectory
% that is closest to the point i in the test trajectory. 
%
% 28/01/2023
% Andrea Colins Rodriguez

S_ref=size(ref_traj,1);
S_sample=size(sample_traj,1);
idx_test=nan(S_sample,1);

% Both trajectories recur at the some point. We need to limit the search of
% the indeces for the start and end points of the sample trajectory.
for i=1:S_sample
   
    
    if i<=200
        % just look the in the first half of the reference
        from=1;
        to=round(S_ref/2);
    elseif i>S_ref-200
        % just look the in the last half of the reference
        from=round(S_ref/2);
        to=S_ref;
    else
        % use the whole reference trajectory
        from=1;
        to=S_ref;
    end

    dist_i=pdist2(ref_traj(from:to,:),sample_traj(i,:));
    
       [~,idx_test(i)]=min(dist_i);
       idx_test(i)=idx_test(i)+from-1;
end


%% Linear regression
ft = fittype('a*x+b');
f1 = fit((1:S_sample)',idx_test,ft,'StartPoint',[0 0]);

slope=f1.a;

if slope<0
    % slope should not be negative!! (it means theres no correlation between trajectories)
    % if this happens, show the user
    disp('Ratio between speeds is negative. Sample and reference trajectories seem to have different geometry')
end

end