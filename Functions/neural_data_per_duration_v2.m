function [condition_matrix,direction,total_matrix,trial_idx,trial_number,reach_number,mov_dist,mov_duration,max_speed,position]=neural_data_per_duration_v2(session,Area,sigma_filter,t_1,t_2,event,duration_range)
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


total_matrix=[];
trial_idx=[];
ntimebins=round((t_2-t_1)*1000);

condition_matrix=zeros(size(neural_data,2),ntimebins);


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
            direction(counter)=find_upcoming_direction([ntargetx(target) ntargety(target)],[ntargetx(target+1) ntargety(target+1)]);
            
            if ~isnan(direction(counter))
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
                trial_idx=[trial_idx,ones(1,size(matrix2,2))*counter];
                trial_number(counter)=Ntrial(i);
                reach_number(counter)=target;
                
                % movement parameters
                mov_duration(counter)=mov_length(target);
                position(counter,:)= [ntargetx(target) ntargety(target)];
                mov_dist(counter)=sqrt((ntargetx(target)-ntargetx(target+1)).^2+ (ntargety(target)-ntargety(target+1)).^2);
               
                start2=round((startt-t_1+shiftbase/1000-1)*1000);  
                speed=sqrt(sum(cont.vel(start2:start2+round(mov_length(target)*1000),:).^2,2));
                max_speed(counter)=max(speed);

                % Neural activity 
                condition_matrix(:,:,counter)=matrix2;
                total_matrix(:,:,counter)=spikes;
               
                counter=counter+1;
            else
                direction(counter)=[];
            end
        end
    end
end

end