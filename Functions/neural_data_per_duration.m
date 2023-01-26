function [Neural_info,Mov_params,Reach_info]=neural_data_per_duration(Session,Area,sigma_filter,t_from,t_upto,event,duration_range)
%% neural_data_per_duration selects the neural activity and the movement parameters for all movements within the duration range
%
% INPUTS
%
% Session= Name of the session to be analysed.
% e.g 'MC_S1_raw.mat'
%
% Area= Name of the area to be analysed in the Session.
% e.g 'M1'
%
% sigma_filter= size of the gaussian filter to convolve the spike trains with.
%
% t_from= start time of the neural activity relative to movement onset [S]
% e.g t_from=-0.5
%
% t_upto= end time of the neural activity relative to the movement onset [S]
% e.g t_from=0.6
%
% event = Indicates the event to which the neural activity will be aligned
%         1- Target onset
%         2- Movement onset
%
% duration_range= range of the durations of the movements [S] that will be
% selected. e.g [0.2 0.3]
%
% OUTPUTS
%
% Neural_info = structure containing the following fields
%                .FR = [nunits ntimebins n movements] each matrix contain the filtered spike trains for each selected movement.
%                .spikes = [nunits ntimebins n movements] each matrix contain the spike trains for each selected movement.
%
% Mov_params = structure containing the following fields
%                .distance= distance of all selected movement [cm].
%                .duration = duration of all selected movement [ms].
%                .max_distance = maximum speed of all selected movement [cm/s].
%                .position= x and y coordinates of the origin target [cm].
%
% Reach_info = structure containing the following fields                
%                .trial_number= trial number of each selected movement.
%                .reach_number= number of the movement within the trial (1-4).
%
% 25/01/2023
% Andrea Colins Rodriguez

load(Session,Area,'cont','trial_table2')
if strcmp(Area,'PMd')
    neural_data=PMd.units;
elseif strcmp(Area,'M1')
    neural_data=M1.units;
end

event_t=event-2;
Ntrial=1:size(trial_table2,1);

ntimebins=round((t_upto-t_from)*1000);

Neural_info.FR=zeros(size(neural_data,2),ntimebins);
Neural_info.spikes_matrix=zeros(size(neural_data,2),ntimebins);

% Define filter w for the neural activity
x = -5*sigma_filter:1:5*sigma_filter';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
shiftbase=floor(numel(x)/2);
w = (1/(sqrt(2*pi*sigma_filter^2)))*exp(-((x.^2*(1/(2*sigma_filter^2))))); % y-axis values of the Gaussian
w = w ./sum(w);

counter=1;
for i=1:numel(Ntrial)
    
    %reach
    mov_length=trial_table2(Ntrial(i),[3 8 13 18]+1)-trial_table2(Ntrial(i),[3 8 13 18]);
    ntarget_time=trial_table2(Ntrial(i),[3 8 13 18]+event_t)+t_from-shiftbase/1000;
    nonset_time=trial_table2(Ntrial(i),[3 8 13 18]+event_t)+t_upto+shiftbase/1000;
    
    % For the first movement, define the position of target 0 as the
    % position of the cursor at time of target onset 
    
    %Note: hand position time start at t==1s, that is why I need to
    %substract 1 when computing hand position
    if isnan(ntarget_time(1))
        xfirst=nan;
        yfirst=nan;
    else
        xfirst= cont.pos(round((trial_table2(Ntrial(i),2)-1)*1000),1);
        yfirst= cont.pos(round((trial_table2(Ntrial(i),2)-1)*1000),2);
    end
    
    ntargetx=[xfirst trial_table2(Ntrial(i),[5 10 15 20])];
    ntargety=[yfirst trial_table2(Ntrial(i),[6 11 16 21])];
    
    for target=1:4
        startt=ntarget_time(target);
        endt=nonset_time(target);
        
        % Select only movements which timing is defined and their duration
        % are within the selected range
        if ~isnan(endt) && ~isnan(startt) && mov_length(target)>duration_range(1) && mov_length(target)<=duration_range(2)
            Mov_params.direction(counter)=find_upcoming_direction([ntargetx(target) ntargety(target)],[ntargetx(target+1) ntargety(target+1)]);
            
            if ~isnan(Mov_params.direction(counter))
                matrix=spikest2vector(neural_data,startt,endt);
                
                
                % filter the neural activity
                filt= filter(w,1,matrix,[],2);
                matrix2=filt(:,2*shiftbase+1:end);
                spikes=matrix(:,shiftbase+1:size(matrix2,2)+shiftbase);
                
                % spikes and matrix2 should have ntimebins columns
                if size(matrix2,2)<ntimebins
                    matrix2=[matrix2,matrix2(:,end)];
                    spikes=[spikes,spikes(:,end)];
                elseif size(matrix2,2)>ntimebins
                    matrix2=matrix2(:,1:ntimebins);
                    spikes=spikes(:,1:ntimebins);
                end
                
                % trial and reach info
                Reach_info.trial_number(counter)=Ntrial(i);
                Reach_info.reach_number(counter)=target;
                
                % movement parameters
                Mov_params.duration(counter)=mov_length(target);
                Mov_params.position(counter,:)= [ntargetx(target) ntargety(target)];
                Mov_params.distance(counter)=sqrt((ntargetx(target)-ntargetx(target+1)).^2+ (ntargety(target)-ntargety(target+1)).^2);
               
                start2=round((startt-t_from+shiftbase/1000-1)*1000);  
                speed=sqrt(sum(cont.vel(start2:start2+round(mov_length(target)*1000),:).^2,2));
                Mov_params.max_speed(counter)=max(speed);

                % Neural activity 
                Neural_info.FR(:,:,counter)=matrix2;
                Neural_info.spikes_matrix(:,:,counter)=spikes;
               
                counter=counter+1;
            else
                Mov_params.direction(counter)=[];
            end
        end
    end
end

end