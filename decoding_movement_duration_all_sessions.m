function decoding_movement_duration_all_sessions(Sessions,Areas,threshold,Ndir,Nbins)
%% decoding_movement_duration_all_sessions predicts the duration of the neural trajectories of each recordings based on their relative speed to a reference trajectory
%
% INPUTS
%
% Sessions: cell array containing the names of the sessions to be analysed.
% e.g {'MC_S1_raw.mat','MC_S2_raw.mat'}
%
% Areas: cell array containing the names of the areas to be analysed in
% each session. Areas and Sessions must have the same number of elements.
% e.g {'M1','M1'}
%
% threshold: percentage of the variance to be explained by the first nPCs
%
% Ndir: number of directions to bin the movements
%
% Nbins: number of durations to bin the movements
%
% 28/05/2023
% Andrea Colins Rodriguez

figure
ref_bin=2;
Nsessions=size(Sessions,2);
ndim=zeros(Nsessions,1);


estimated_dur_all=zeros(Ndir*Nsessions,Nbins-1);
M1=nan(Ndir*Nsessions,1);
mean_total_error=zeros(Nsessions,1);
std_error=zeros(Nsessions,1);
colour_dir=hsv(Ndir);

test_bins=1:Nbins;
test_bins(ref_bin)=[];
Ntest_bins=numel(test_bins);



for isession=1:Nsessions
    load(['../Output_files/PCA_' Sessions{isession}(1:end-4) '_' Areas{isession} '.mat'],'score','idx_dir','idx_duration','variance','t_from','t_upto')

    ndim(isession)=find(cumsum(variance)>threshold,1,'First');
    for i_bin=1:Ntest_bins
        bin_i=test_bins(i_bin);
        traj_length(i_bin)=sum((idx_dir==1)& (idx_duration==bin_i));
    end
    if strcmp(Areas{isession},'M1')
        traj_length_M1=(abs(t_from)+t_upto)*1000;
        traj_length_M1(ref_bin)=[];
    else
        traj_length_PMd=(abs(t_from)+t_upto)*1000;
        traj_length_PMd(ref_bin)=[];
    end

    if isession==1
        do_plot=1;
    else
        do_plot=0;
    end

    [mean_total_error(isession),std_error(isession),estimated_duration]=decoding_movement_duration(score,idx_dir,idx_duration,ref_bin,ndim(isession),Nbins,do_plot);
    subplot(2,3,4)
    hold on
    for i_dir=1:Ndir
        if strcmp(Areas{isession},'M1')
            plot(traj_length-450,estimated_duration(i_dir,:)-450,'.','Color',colour_dir(i_dir,:))
        else
            plot(traj_length-500,estimated_duration(i_dir,:)-500,'.','Color',colour_dir(i_dir,:))
        end
    end

    estimated_dur_all(Ndir*(isession-1)+1:Ndir*isession,:)=estimated_duration;


    if strcmp(Areas{isession},'M1')
        M1(Ndir*(isession-1)+1:Ndir*isession)=1;
    end

end

subplot(2,3,4)
errorbar(traj_length_M1-450,mean(estimated_dur_all(M1==1,:))-450,std(estimated_dur_all(M1==1,:)),'Color',[85 30 116]./256)
errorbar(traj_length_PMd-500,mean(estimated_dur_all(M1~=1,:))-500,std(estimated_dur_all(M1~=1,:)),'Color',[89 156 153]./256)
box off
ylabel('Predicted trajectory duration [ms]')
xlabel('Actual trajectory duration [ms]')
plot(traj_length_M1-450,traj_length_M1-450,'k')
xlim([traj_length_M1(1)-600 traj_length_PMd(end)-400])

subplot(2,3,5)
error_prop_M1=100*abs(estimated_dur_all(M1==1,:)-repmat(traj_length_M1,sum(M1==1),1))./repmat(traj_length_M1,sum(M1==1),1);
error_prop_PMd=100*abs(estimated_dur_all(M1~=1,:)-repmat(traj_length_PMd,sum(M1~=1),1))./repmat(traj_length_PMd,sum(M1~=1),1);
relative_error_M1=abs(error_prop_M1);
relative_error_PMd=abs(error_prop_PMd);
errorbar(traj_length_M1,mean(abs(error_prop_M1)),std(abs(error_prop_M1)),'Color',[85 30 116]./256)
errorbar(traj_length_PMd,mean(abs(error_prop_PMd)),std(abs(error_prop_PMd)),'Color',[89 156 153]./256)

text(700,100,['Mean M1 = ' num2str(mean(relative_error_M1(:)),2)])
text(1000,100,['Mean PMd = ' num2str(mean(relative_error_PMd(:)),2)])
box off
ylabel('Relative error [%]')
xlabel('Actual trajectory duration [ms]')
xlim([traj_length_M1(1)-100 traj_length_PMd(end)+100])


%% shuffle test: Shuffle direction label to test how breaking the hypothesis of geometry changes the results
% aim: to produce a plot that shows the difference in movement angle vs error of
% decoding
clear error_duration_tmp
angle_diff=linspace(0,360,Ndir+1)-180;
angle_diff(end)=[];
angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); %center at zero

for isession=1:Nsessions
    load(['../Output_files/PCA_' Sessions{isession}(1:end-4) '_' Areas{isession} '.mat'],'score','idx_dir','idx_duration','variance')

    ndim(isession)=find(cumsum(variance)>threshold,1,'First');

    error_by_dir_shift=decoding_movement_duration_chance(score,idx_dir,idx_duration,ref_bin,ndim(isession));


    if strcmp(Areas{isession},'M1')
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