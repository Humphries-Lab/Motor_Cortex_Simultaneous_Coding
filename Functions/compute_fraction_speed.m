function [slope,idx_test]=compute_fraction_speed(xref,xsample)
%% compute_fraction_speed computes the ratio between the speeds of the reference and the sample trajectory 
%
% INPUTS 
%
% xref: reference trajectory. Rows are time bins, columns are dimensions
% 
% xsample: sample trajectory. Rows are time bins, columns are dimensions
% 
% OUTPUTS
%
% slope: ratio between speed of the sample and reference trajectories.
% slope>1=> sample is faster
%
% idx_test: idx_test(i) is the time at which the reference trajectory
% is closest to the point i in the test trajectory. 
%
% 28/01/2023
% Andrea Colins Rodriguez

S_ref=size(xref,1);
S_sample=size(xsample,1);
idx_test=nan(S_sample,1);

% Both trajectories recur at the some point. We need to limit the search of
% the indeces for the start and end points of the sample trajectory.
for i=1:S_sample
   
    xref_tmp=nan(size(xref));
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
    
    xref_tmp(from:to,:)=xref(from:to,:);
    dist_i=pdist2(xref_tmp,xsample(i,:));
    
    [~,idx_test(i)]=min(dist_i);

end


%% Linear regression
ft = fittype('a*x+b');
f1 = fit((1:S_sample)',idx_test,ft,'StartPoint',[0 0]);

slope=f1.a;

if slope<0
    % slope should not be negative!! (it means theres no correlation between trajectories)
    % if this happens, show the user
    disp('Ratio between speeds is negative')
end

end