function ISI=compute_ISI(units,startt,endt)
%% compute_ISI calculates the Inter-Spike-Inverval of the population (units) between times startt and endt
%
%INPUTS
%
% units: structure containing the time of the spikes of all the units in
% the population recording. units has dimensions 1 x number of units and
% the time of spikes (in seconds) is in the field ts.
%
% startt= selected start time (in seconds) of the neural activity
% 
% endt= selected end time (in seconds) of the neural activity
%
%
% OUTPUT
%
% ISI= array containing the Inter-Spike-Intervals of the population (in ms)
%
% 23/01/2023
% Andrea Colins Rodriguez


ISI=[];

for iunit=1:size(units,2)
    idxsp=find(((units(iunit).ts>=startt) & (units(iunit).ts<endt)));
    ISI_iunit=diff(units(iunit).ts(idxsp));
    ISI=[ISI;ISI_iunit(:)];
end

ISI=ISI*1000; % return output in ms

end