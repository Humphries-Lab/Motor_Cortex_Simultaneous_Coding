function [Neural_info,Mov_params,Reach_info]=neural_data_per_duration_normalised(session,Area,sigma_filter,t_1,t_2,event,duration_range)
%% neural_data_per_duration selects the neural activity and the movement parameters for all movements within the duration range
% INPUTS
%
% OUTPUTS
%

load(session,Area,'cont','trial_table2')
if strcmp(Area,'PMd')
    neural_data=PMd.units;
elseif strcmp(Area,'M1')
    neural_data=M1.units;
end


event_t=event-2;
Ntrial=1:size(trial_table2,1);


Reach_info.trial_idx=[];
ntimebins=600; %all outputs are scaled to have n timebins
normalised_t=linspace(0,1,ntimebins);
Neural_info.FR=zeros(size(neural_data,2),ntimebins);

% Define filter w for the neural activity
x = -5*sigma_filter:1:5*sigma_filter';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
shiftbase=floor(numel(x)/2);
w = (1/(sqrt(2*pi*sigma_filter^2)))*exp(-((x.^2*(1/(2*sigma_filter^2))))); % y-axis values of the Gaussian
w = w ./sum(w);

counter=1;
for i=1:numel(Ntrial)
    
    %reach
    mov_length=trial_table2(Ntrial(i),[3 8 13 18]+1)-trial_table2(Ntrial(i),[3 8 13 18]);
    ntarget_time=trial_table2(Ntrial(i),[3 8 13 18]+event_t)+t_1-shiftbase/1000;
    nonset_time=trial_table2(Ntrial(i),[3 8 13 18]+event_t)+t_2+shiftbase/1000;
    
    movement_init_t=round((trial_table2(Ntrial(i),[3 8 13 18])-1)*1000);
    movement_end_t=round((trial_table2(Ntrial(i),[3 8 13 18]+1)-1)*1000);
    
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
                FR=filt(:,2*shiftbase+1:end);
                % scale neural activity to have n timebins
                xtime2=linspace(0,1,size(FR,2));
                Neural_info.FR(:,:,counter)=interp1(xtime2,FR',normalised_t)';
                
                % trial and reach info
                Reach_info.trial_idx=[Reach_info.trial_idx,ones(1,size(FR,2))*counter];
                Reach_info.trial_number(counter)=Ntrial(i);
                Reach_info.reach_number(counter)=target;
                
                % movement parameters
                Mov_params.duration(counter)=mov_length(target);
                Mov_params.position(counter,:)= [ntargetx(target) ntargety(target)];
                Mov_params.distance(counter)=sqrt((ntargetx(target)-ntargetx(target+1)).^2+ (ntargety(target)-ntargety(target+1)).^2);
                
                speed=sqrt(sum(cont.vel(movement_init_t(target):movement_end_t(target),:).^2,2));
                Mov_params.max_speed(counter)=max(speed);
                
                counter=counter+1;
            else
                Mov_params.direction(counter)=[];
            end
        end
    end
end

end