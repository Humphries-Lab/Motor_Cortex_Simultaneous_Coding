function sp_matrix=spikest2vector(units,startt,endt)
%% spikes2vector transforms the time of the spikes in the population units to a matrix format
%INPUTS
%
% units: structure containing the time of the spikes of all the units in
% the population recording. units has dimensions 1 x number of units and
% the time of spikes (in seconds) is in the field ts.
%
% startt = selected start time (in seconds) of the neural activity
% 
% endt = selected end time (in seconds) of the neural activity
%
%
% OUTPUT
%
% sp_matrix = matrix containing the spikes of each neuron of the population
% (in ms). sp_matrix has dimensions of [n units, time (in ms)]. 
%
% 23/01/2023
% Andrea Colins Rodriguez


% rows contain units and columns represent t between startt and endt
sp_matrix=zeros(size(units,2),round((endt-startt)*1000));


for unit=1:size(units,2)

    idxsp=(units(unit).ts>=startt) & (units(unit).ts<endt);
    
    % convert time of the spike to ms from startt
    tsp=round((units(unit).ts(idxsp)-startt)*1000)+1;
    sp_matrix(unit,tsp)=1;
    
    % in the odd case that there are two spikes in the same ms time bin,
    % let the user know
    if numel(unique(tsp))<numel(tsp)
        disp('double spike')
        N = histcounts(tsp,0:max(tsp));
        double_spike=find(N>1);
    end
end

end