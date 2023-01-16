function [condition_matrix,direction,total_matrix,trial_idx,trial_number,reach_number,dist_mov_dir,mov_duration,max_speed,prep_duration,position]=neural_data_per_duration_normalised(session,Area,ms,t_1,t_2,event,from_to)
load(session)
if strcmp(Area,'PMd')
    neural_data=PMd.units;
elseif strcmp(Area,'M1')
    neural_data=M1.units;
end
final_length=600;
normalised_t=linspace(0,1,final_length);
event_t=event-2;
Ntrial=1:size(trial_table,1);
Ntrial(isnan(sum(trial_table(:,[5 6 10 11 15 16 20 21]),2)))=[];
trial_table=trial_table2;

total_matrix=[];
trial_idx=[];
counter=1;
% w = gausswin(ms);
% w=w./max(cumsum(w));
x = -5*ms:1:5*ms';  % x-axis values of the discretely sampled Gaussian, out to 5xSD
shiftbase=floor(numel(x)/2);
w = (1/(sqrt(2*pi*ms^2)))*exp(-((x.^2*(1/(2*ms^2))))); % y-axis values of the Gaussian
w = w ./sum(w); 

for i=1:numel(Ntrial)
    
    %reach
    mov_length=trial_table(Ntrial(i),[3 8 13 18]+1)-trial_table(Ntrial(i),[3 8 13 18]);
    prep_length=trial_table(Ntrial(i),[3 8 13 18])-trial_table(Ntrial(i),[3 8 13 18]-1);
    

    ntarget_time=trial_table(Ntrial(i),[3 8 13 18]+event_t)+t_1-shiftbase/1000;
    nonset_time=trial_table(Ntrial(i),[3 8 13 18]+event_t+1)+t_2+shiftbase/1000;
    
    
    movement_init_t=round((trial_table(Ntrial(i),[3 8 13 18])-1)*1000);
    movement_end_t=round((trial_table(Ntrial(i),[3 8 13 18]+1)-1)*1000);
    %Note: hand position time start at t==1s, that is why I need to
    %substract 1 when computing hand position
    if isnan(ntarget_time(1))
        xfirst=nan;
        yfirst=nan;
    else
        xfirst= cont.pos(round((trial_table(Ntrial(i),2)-1)*1000),1);
        yfirst= cont.pos(round((trial_table(Ntrial(i),2)-1)*1000),2);
    end
    
    ntargetx=[xfirst trial_table(Ntrial(i),[5 10 15 20])];
    ntargety=[yfirst trial_table(Ntrial(i),[6 11 16 21])];
      %speed_target=trial_table(Ntrial(i),[4 9 14 19]);
    for target=1:4
        startt=ntarget_time(target);
        endt=nonset_time(target);
        
        if ~isnan(endt) && ~isnan(startt) && mov_length(target)>from_to(1) && mov_length(target)<=from_to(2)% && ((endt-startt)<=1)
            direction_tmp=find_upcoming_direction([ntargetx(target) ntargety(target)],[ntargetx(target+1) ntargety(target+1)]);
            if ~isnan(direction_tmp)
                direction(counter)=direction_tmp;
                matrix=spikest2vector(neural_data,startt,endt);
                %matrix_reach=spikest2vector(PMd.units,nonset_time2b(target)-2*ms/1000,endt2);
                %shuffle_test
                %                   p=randperm(size(matrix,2));
                %                   matrix=matrix(:,p);
                
                filt= filter(w,1,matrix,[],2);
                 matrix2=filt(:,2*shiftbase+1:end);
%                 [size(matrix2,2) size(matrix,2) (mov_length(target)+0.3+0.2)*1000 (endt-startt)*1000]
%                 %debugging
%                   plot(matrix(1,:))
%                   hold on
%                    plot(matrix2(1,:))
%                    pause

                    
                % compute max speed
                speed_target= sqrt(sum(cont.vel(movement_init_t(target):movement_end_t(target),:).^2,2));
                
                %save everything
                total_matrix=[total_matrix,matrix2];
                %total_matrix_reach=[total_matrix_reach,matrix_reach2];
                trial_idx=[trial_idx,ones(1,size(matrix2,2))*counter];
                position(counter,:)=cont.pos(movement_init_t(target),:);
                trial_number(counter)=Ntrial(i);
                reach_number(counter)=target;
                mov_duration(counter)=mov_length(target);
                prep_duration(counter)=prep_length(target);
                max_speed(counter)=max(speed_target);
                dist_mov_dir(counter)=sqrt((ntargetx(target)-ntargetx(target+1)).^2+ (ntargety(target)-ntargety(target+1)).^2);
              
               %%%% normalising time
               
               %%original
               %condition_matrix(:,:,counter)=matrix2;
               
               %normalised
               xtime2=linspace(0,1,size(matrix2,2));
               condition_matrix(:,:,counter)= interp1(xtime2,matrix2',normalised_t)';
               
               

                counter=counter+1;
            end
        end
    end
end

end