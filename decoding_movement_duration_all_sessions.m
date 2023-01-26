function decoding_movement_duration_all_sessions(session,Area,Ndir,Nbins)
figure
ref_bin=2;
Nsessions=size(session,2);
ndim=zeros(Nsessions,1);
threshold=80;
kneigh=1;
%mean_total_error=zeros(Nsessions,Ndim);
estimated_dur_all=zeros(Ndir*Nsessions,Nbins-1);
error_duration=zeros(Ndir*Nsessions,Nbins-1);
M1=nan(Ndir*Nsessions,1);
mean_total_error=zeros(Nsessions,1);
std_error=zeros(Nsessions,1);
for j=1:Nsessions
    load(['../Output_files/PCA_' session{j}(1:end-4) '_' Area{j} '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto')
   
    %load(['scores_LDS_diff_duration_newfilter_' session{j} '_' Area{j} '.mat'],'score','idx_dir','idx_duration','variance','t_1','t_2','nsamples_condition')
    ndim(j)=find(cumsum(variance)>threshold,1,'First');
    if strcmp(Area{j},'M1')
        traj_length_M1=(abs(t_from)+t_upto)*1000;
        traj_length_M1(ref_bin)=[];
    else
        traj_length_PMd=(abs(t_from)+t_upto)*1000;
        traj_length_PMd(ref_bin)=[];
    end
    
    if j==1
        do_plot=1;
    else
        do_plot=0;
    end
    subplot(2,3,4)
    [mean_total_error(j),error_duration_tmp,std_error(j),estimated_duration]=decoding_movement_duration(score,idx_dir,idx_duration,ref_bin,ndim(j),Nbins,kneigh,do_plot);
    
    error_duration(Ndir*(j-1)+1:Ndir*j,:)=error_duration_tmp;
    estimated_dur_all(Ndir*(j-1)+1:Ndir*j,:)=estimated_duration;
    %error_per_sample=100*abs(estimated_duration-repmat(traj_length_M1,Ndir,1))./repmat(traj_length_M1,Ndir,1);
    

    if strcmp(Area{j},'M1')
        M1(Ndir*(j-1)+1:Ndir*j)=1;
    end
    
end

subplot(2,3,4)
errorbar(traj_length_M1,mean(estimated_dur_all(M1==1,:)),std(estimated_dur_all(M1==1,:)),'Color',[85 30 116]./256)
errorbar(traj_length_PMd,mean(estimated_dur_all(M1~=1,:)),std(estimated_dur_all(M1~=1,:)),'Color',[89 156 153]./256)
box off
ylabel('Predicted trajectory duration [ms]')
xlabel('Actual trajectory duration [ms]')
plot(traj_length_M1,traj_length_M1,'k')
xlim([traj_length_M1(1)-100 traj_length_PMd(end)+100])

subplot(2,3,5)
error_prop_M1=100*abs(estimated_dur_all(M1==1,:)-repmat(traj_length_M1,sum(M1==1),1))./repmat(traj_length_M1,sum(M1==1),1);
error_prop_PMd=100*abs(estimated_dur_all(M1~=1,:)-repmat(traj_length_PMd,sum(M1~=1),1))./repmat(traj_length_PMd,sum(M1~=1),1);
relative_error_M1=abs(error_prop_M1);
relative_error_PMd=abs(error_prop_PMd);
errorbar(traj_length_M1,mean(abs(error_prop_M1)),std(abs(error_prop_M1)),'Color',[85 30 116]./256)
errorbar(traj_length_PMd,mean(abs(error_prop_PMd)),std(abs(error_prop_PMd)),'Color',[89 156 153]./256)
%[mean(relative_error_M1(:)) mean(relative_error_PMd(:))]
text(700,100,['Mean M1 = ' num2str(mean(relative_error_M1(:)),2)])
text(1000,100,['Mean PMd = ' num2str(mean(relative_error_PMd(:)),2)])
box off
ylabel('Relative error [%]')
xlabel('Actual trajectory duration [ms]')
xlim([traj_length_M1(1)-100 traj_length_PMd(end)+100])

% subplot(3,3,5)
% errorbar(1:Nsessions,mean_total_error,std_error)
% hold on
% box off
% ylabel('Median error [ms]')
% xlim([0.8 Nsessions+0.2])

%% shuffle test: Shuffle direction label to test how breaking the hypothesis of geometry changes the results
%aim: produce a plot that shows diff angle vs error
clear error_duration_tmp
angle_diff=linspace(0,360,Ndir+1)-180;
angle_diff(end)=[];
angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); %center at zero
error_dir=zeros(Ndir*Nsessions*3,Ndir);
for j=1:Nsessions
    load(['../Output_files/PCA_' session{j}(1:end-4) '_' Area{j} '.mat'],'score','idx_dir','idx_duration','variance')
   
    %load(['scores_LDS_diff_duration_newfilter_' session{j} '_' Area{j} '.mat'],'score','idx_dir','idx_duration','variance')
    ndim(j)=find(cumsum(variance)>threshold,1,'First');

    [error_duration_tmp(:,j),error_by_dir_shift]=decoding_movement_duration_chance(score,idx_dir,idx_duration,ref_bin,ndim(j),kneigh);
    error_dir(Ndir*(j-1)*3+1:Ndir*j*3,:)=error_by_dir_shift';
    
    if strcmp(Area{j},'M1')
        colour_Area=[85 30 116]./256;
    else
        colour_Area=[89 156 153]./256;
    end
    
     subplot(2,3,6)
        plot(angle_diff,median(error_by_dir_shift,2),'Color',colour_Area)
        hold on
        box off
        ylabel('Median error')
        xlabel('Angle difference')
    
    
end

end