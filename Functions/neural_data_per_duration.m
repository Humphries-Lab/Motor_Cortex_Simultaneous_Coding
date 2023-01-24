function [condition_matrix,direction,total_matrix,trial_idx,trial_number,reach_number,dist_mov_dir,mov_duration,max_speed,position]=neural_data_per_duration(session,Area,ms,t_1,t_2,event,from_to,doplot)
load(session)
if strcmp(Area,'PMd')
    neural_data=PMd.units;
elseif strcmp(Area,'M1')
    neural_data=M1.units;
end

event_t=event-2;
Ntrial=1:size(trial_table,1);
trial_table=trial_table2;
%trial_table(:,column_start)=trial_table(:,column_start)+0.096;
%trial_table(:,column_start(2:4))=trial_table(:,column_start(2:4))-0.1;

total_matrix=[];
trial_idx=[];
segment1=round((t_2-t_1)*1000);

condition_matrix=zeros(size(neural_data,2),segment1);
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
    ntarget_time=trial_table(Ntrial(i),[3 8 13 18]+event_t)+t_1-shiftbase/1000;
    nonset_time=trial_table(Ntrial(i),[3 8 13 18]+event_t)+t_2+shiftbase/1000;
    
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
    
    for target=1:4
        startt=ntarget_time(target);
        endt=nonset_time(target);
        
        if ~isnan(endt) && ~isnan(startt) && mov_length(target)>from_to(1) && mov_length(target)<=from_to(2)% && ((endt-startt)<=1)
            direction(counter)=find_upcoming_direction([ntargetx(target) ntargety(target)],[ntargetx(target+1) ntargety(target+1)]);
            if ~isnan(direction(counter))
                matrix=spikest2vector(neural_data,startt,endt);
                %matrix_reach=spikest2vector(PMd.units,nonset_time2b(target)-2*ms/1000,endt2);
                %shuffle_test
                %                   p=randperm(size(matrix,2));
                %                   matrix=matrix(:,p);
                
                
                
                filt= filter(w,1,matrix,[],2);
                matrix2=filt(:,2*shiftbase+1:end);
                spikes=matrix(:,shiftbase+1:size(matrix2,2)+shiftbase);
                %debugging
                %                  size(matrix2)
                %                    plot(matrix(1,:))
                %                    hold on
                %                     plot(matrix2(1,:))
                %                     pause
                
                %                matrix2= filter(w,1,matrix')';
                %                 spikes=matrix(:,1:size(matrix,2)-ms-1);
                %                 spikes(:,1:ms)=[];
                %                 matrix2(:,size(matrix,2)-ms:end)=[];
                %                 matrix2(:,1:ms)=[];
                
                if size(matrix2,2)<segment1
                    matrix2=[matrix2,matrix2(:,end)];
                    spikes=[spikes,spikes(:,end)];
                elseif size(matrix2,2)>segment1
                    matrix2=matrix2(:,1:segment1);
                    spikes=spikes(:,1:segment1);
                end
                
                
                %total_matrix_reach=[total_matrix_reach,matrix_reach2];
                trial_idx=[trial_idx,ones(1,size(matrix2,2))*counter];
                trial_number(counter)=Ntrial(i);
                reach_number(counter)=target;
                mov_duration(counter)=mov_length(target);
                start2=round((startt-t_1+shiftbase/1000-1)*1000);
  
                speed=sqrt(sum(cont.vel(start2:start2+round(mov_length(target)*1000),:).^2,2));
                max_speed(counter)=max(speed);
                dist_mov_dir(counter)=sqrt((ntargetx(target)-ntargetx(target+1)).^2+ (ntargety(target)-ntargety(target+1)).^2);
                condition_matrix(:,:,counter)=matrix2;
                total_matrix(:,:,counter)=spikes;
                position(counter,:)= [ntargetx(target) ntargety(target)];
                %condition_matrix_reach(:,:,counter)=matrix_reach2;
                counter=counter+1;
            else
                direction(counter)=[];
            end
        end
    end
end
%average for condition
% direction2=ceil(10*(direction+pi)/(2*pi));
% average_cond=[];
% for i=1:10
%     average_cond=[average_cond,mean(condition_matrix(:,:,direction2==i),3)];
%
% end
%soft normalisation
% total_matrix=total_matrix'./repmat(range(average_cond')+5,size(total_matrix,2),1);
% tmp_range=repmat(range(average_cond')+5,size(average_cond,2),1);
% average_cond=average_cond'./repmat(range(average_cond')+5,size(average_cond,2),1);
end