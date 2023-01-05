function distance_position_all_sessions(session,Area,threshold,Ndir,Nbins,do_plot,t_1,t_2,from)

fig1=figure;
embedding_dim=zeros(size(session));
for i=1:size(session,2)
    i
    variance=embedding_dimensions2(session{i}, Area{i},do_plot,Ndir,Nbins,t_1(i),t_2(i,:),from);
    embedding_dim(i)=find(cumsum(variance)>threshold,1,'First');
    pause
end
figure(fig1)
subplot(2,1,1)
plot([0 100],[threshold threshold])
figure(fig1)
subplot(2,1,2)
hold on
plot(embedding_dim,'k')
box off
xlabel('N session')
ylabel('Embedding dimensions')
end


function [variance,t_1_new,t_2_new]=embedding_dimensions2(session, Area,do_plot,Ndir,Nbins,t_1,t_2,from)
%% Average across bins of segments durations (200,300,400,500 ms)
event=2;
if do_plot
    figure
    colour_dir=hsv(Ndir);
end

if strcmp(Area,'PMd')
    load(session,'PMd','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(PMd.units,startt,endt);
    %t_1=-0.5;
else
    load(session,'M1','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(M1.units,startt,endt);
    %t_1=-0.5;
end

ms=nanmedian(ISI); %round(ms*sqrt(12))
t=-pi:0.1:pi;
xcircle=sin(t);
ycircle=cos(t);
% t_2=from+0.1+0.5;%movement duration+ 200 ms
average_cond_1=[];
idx_dir=[];
idx_duration=[];
colour_plasma=plasma(Nbins);
counter=1;
nsamples_condition=zeros(Ndir,Nbins);
Neural_position=[];
idx_pos=[];
idx_dir2=[];
for i=1:Nbins
    from_to=[from(i) from(i)+0.1];
    [condition_matrix{i},direction,~,~,~,reach_number{i},dist_mov_dir{i},mov_duration{i},max_speed{i},position{i}]=neural_data_per_duration(session,Area,ms,t_1,t_2(i),event,from_to,0);
    xtime=round(t_1*1000):1:round(t_2(i)*1000-1);
    %take average by direction
    direction1=ceil(Ndir*(direction+pi)/(2*pi));
    
    for cond=1:Ndir
        nsamples_condition(cond,i)=sum(direction1==cond);
        subplot(2,2,3)
        if nsamples_condition(cond,i)>=2
            %plot(mean(position{i}(direction1==cond,1)),mean(position{i}(direction1==cond,2)),'o','Color',colour_plasma(i,:))
            plot(position{i}(direction1==cond,1),position{i}(direction1==cond,2),'o','Color',colour_dir(cond,:))
            hold on
            %Divide these samples into four groups and do PCA only there,
            %if there is a significant difference it has to be there
            %pause
            xsep=(min(position{i}(direction1==cond,1))+max(position{i}(direction1==cond,1)))/2;
            ysep=(min(position{i}(direction1==cond,2))+max(position{i}(direction1==cond,2)))/2;
            %group 1
            xy=position{i};
            cond1=direction1'==cond & xy(:,1)>xsep & xy(:,2)>ysep;
            cond2=direction1'==cond & xy(:,1)<xsep & xy(:,2)>ysep;
            cond3=direction1'==cond & xy(:,1)<xsep & xy(:,2)<ysep;
            cond4=direction1'==cond & xy(:,1)>xsep & xy(:,2)<ysep;
            
            g1=position{i}(cond1,1:2);
            g2=position{i}(cond2,1:2);
            g3=position{i}(cond3,1:2);
            g4=position{i}(cond4,1:2);
            
            %              plot(g1(:,1),g1(:,2),'ob')
            %              plot(g2(:,1),g2(:,2),'or')
            %              plot(g3(:,1),g3(:,2),'og')
            %              plot(g4(:,1),g4(:,2),'oc')
            xlim([-10 10])
            ylim([-10 10])
            pause
            hold off
            idx_dir2=[idx_dir2;zeros(round((t_2(i)-t_1)*1000)*4,1)+cond];
            
            idx_pos=[idx_pos;zeros(round((t_2(i)-t_1)*1000),1)+1;zeros(round((t_2(i)-t_1)*1000),1)+2;zeros(round((t_2(i)-t_1)*1000),1)+3;zeros(round((t_2(i)-t_1)*1000),1)+4];
            
            Neural_position=[Neural_position,mean(condition_matrix{i}(:,:,cond1),3),mean(condition_matrix{i}(:,:,cond2),3),mean(condition_matrix{i}(:,:,cond3),3),mean(condition_matrix{i}(:,:,cond4),3)];
            
            
            average_cond_1=[average_cond_1,mean(condition_matrix{i}(:,:,direction1==cond),3)];
            idx_dir=[idx_dir;zeros(round((t_2(i)-t_1)*1000),1)+cond];
            idx_duration=[idx_duration;zeros(round((t_2(i)-t_1)*1000),1)+i];
            counter=counter+1;
            
        else
            cond
            [session, Area]
            pause
        end
        
        
    end
    [~,score]=pca(Neural_position');
    subplot(2,2,2)
    hold on
    for cond=1:Ndir
        plot3(score(idx_dir2==cond,1),score(idx_dir2==cond,2),score(idx_dir2==cond,3),'Color',colour_dir(cond,:))
    end
    
    subplot(2,2,4)
    plot3(score(idx_pos==1,1),score(idx_pos==1,2),score(idx_pos==1,3),'b')
    hold on
    plot3(score(idx_pos==2,1),score(idx_pos==2,2),score(idx_pos==2,3),'r')
    plot3(score(idx_pos==3,1),score(idx_pos==3,2),score(idx_pos==3,3),'g')
    plot3(score(idx_pos==4,1),score(idx_pos==4,2),score(idx_pos==4,3),'c')
    pause
    hold off
    Neural_position=[];
    idx_pos=[];
    
    
    if do_plot
        
        
        %         subplot(4,5,14)
        %         hold on
        %         plot(dist_mov_dir{i},mov_duration{i},'.','Color',colour_plasma(i,:))
        %         xlabel('Distance [cm]')
        %         ylabel('Duration [ms]')
        %
        %         subplot(4,5,15)
        %         hold on
        %         plot(dist_mov_dir{i},max_speed{i},'.','Color',colour_plasma(i,:))
        %         xlabel('Distance [cm]')
        %         ylabel('Max speed [cm/s]')
        
        %         subplot(4,5,19)
        %         hold on
        %         plot(max_speed{i},mov_duration{i},'.','Color',colour_plasma(i,:))
        %         xlabel('Max speed [cm/s]')
        %         ylabel('Duration [ms]')
        
    end
end
R1=corr([dist_mov_dir{1:Nbins}]',[mov_duration{1:Nbins}]');
R2=corr([dist_mov_dir{1:Nbins}]',[max_speed{1:Nbins}]');
R1=corr([max_speed{1:Nbins}]',[mov_duration{1:Nbins}]');

%% Compute a common subspace
% soft normalization
%delete neurons that have minimal number of spikes
delete_units=sum(average_cond_1,2)<5;
average_cond_1(delete_units,:)=[];
normalisation=range(average_cond_1')+5;
average_cond_1=average_cond_1'./repmat(normalisation,size(average_cond_1,2),1);
%% Project same directions coloring different times
[coeffs,score,~,~,variance]=pca(average_cond_1);


end