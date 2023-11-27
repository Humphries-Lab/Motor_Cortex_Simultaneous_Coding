function [variance,mov_distance,mov_duration,max_speed,Dur_var_exp]=embedding_dimensions_dpca(Session,Area,Ndir,Nbins,t_from,t_upto,edges_dur_bin,do_plot)
%% embedding_dimensions performs PCA on the neural data from Session and
%% Area binning the movements into Ndir directions and Nbins durations
% This function also saves the PCA results on a file for later analyses
%
% INPUTS
%
% Session: Name of the session to be analysed.
% e.g 'MC_S1_raw.mat'
%
% Area: Name of the area to be analysed in the Session.
% e.g 'M1'
%
% Ndir: number of directions to bin the movements
%
% Nbins: number of durations to bin the movements
%
% t_from: start time of the neural activity relative to movement onset [S]
% e.g t_from=-0.5
%
% t_upto: end time of the neural activity relative to movement end [S]
% e.g t_from=0.3
%
% edges_dur_bin= array containing the edges of each duration bin [S]
%
% do_plot= 1- plots the average FR of the population for each duration bin
%             and the PCA trajectories for the first duration bin
%
%          0- omit the above plot
% OUTPUTS
%
% variance: array containing the variance explained by each PC
%
% mov_distance= cell array containing the movement distance of all
% movements selected for PCA. Each cell contains an array with the distance of each movement within a duration bin
%
% mov_duration= cell array containing the movement duration of all
% movements selected for PCA. Each cell contains an array with the duration of each movement within a duration bin
%
% max_speed= cell array containing the maximum speed reached of all
% movements selected for PCA. Each cell contains an array with the maximum speed of each movement within a duration bin
%
% 24/05/2023
% Andrea Colins Rodriguez

ms=1000;

if do_plot
    figure
    colour_dir=hsv(Ndir);
    colour_plasma=plasma(Nbins);
end
load(Session,Area,'trial_table2','cont')
startt=trial_table2(1,1);
endt=trial_table2(size(trial_table2,1),22);

if strcmp(Area,'PMd')
    neural_data=PMd.units;
else
    neural_data=M1.units;
end

ISI=compute_ISI(neural_data,startt,endt);
sigma_filter=round(median(ISI));


nsamples_condition=zeros(Ndir,Nbins);
mov_distance=cell(Nbins,1);
mov_duration=cell(Nbins,1);
max_speed=cell(Nbins,1);
dur_binsize=(edges_dur_bin(2)-edges_dur_bin(1));
t_upto=edges_dur_bin(1:Nbins)+dur_binsize+t_upto;
total_time_bins=sum(round((t_upto(1:Nbins)-t_from)*ms))*Ndir;

average_cond=nan(size(neural_data,2),total_time_bins);
idx_dir=nan(total_time_bins,1);
idx_duration=nan(total_time_bins,1);
counter=0;
Neural_all=cell(Nbins,1);
% Extract movements and related neural activity for eac duration bin
for i_dur=1:Nbins

    current_dur_bin=[edges_dur_bin(i_dur) edges_dur_bin(i_dur+1)]; % select movements in this range of duration only

    [Neural_info,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);
    %[Neural_info,Mov_params]=neural_data_per_duration_normalised(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);

    %Save Mov params for later analyses
    mov_distance{i_dur}=Mov_params.distance;
    mov_duration{i_dur}=Mov_params.duration;
    max_speed{i_dur}=Mov_params.max_speed;

    Neural_all{i_dur}=Neural_info;
    % bin by direction
    direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));


    for i_dir=1:Ndir
        nsamples_condition(i_dir,i_dur)=sum(direction1==i_dir);
        %close_dir{i_dur,i_dir}=direction1==i_dir;
        if nsamples_condition(i_dir,i_dur)>=2


            timebins=round((t_upto(i_dur)-t_from)*ms);
            average_cond(:,counter+1:counter+timebins)=mean(Neural_info.FR(:,:,direction1==i_dir),3);
            idx_dir(counter+1:counter+timebins)=zeros(timebins,1)+i_dir;
            idx_duration(counter+1:counter+timebins)=zeros(timebins,1)+i_dur;

            %             if i_dur==1
            %             FR_dPCA(:,i_dir,:)=mean(Neural_info.FR(:,:,direction1==i_dir),3);
            %             end

            counter=counter+timebins;
        else
            disp([Session ' ' Area ' direction = ' num2str(i_dir) ' has less than 2 samples'])
            keyboard
        end

    end

end

%% dpca bit
%- 1 direction
%- 2 duration
%- 3 time
% normalise FR

delete_units=sum(average_cond,2)<5;
average_cond(delete_units,:)=[];
normalisation=repmat(range(average_cond,2)+5,1,size(average_cond,2));
average_cond=average_cond./normalisation;

for i_dur=1:Nbins
    for i_dir=1:Ndir
        % force to take the same number of points for all conditions
        idx=idx_dir==i_dir & idx_duration==i_dur;
        % normalising the trajectories by their length
        FR_dPCA(:,i_dir,i_dur,:)=interp1(linspace(0,1,sum(idx==1))',average_cond(:,idx)',linspace(0,1,600)')';
    end
end

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}};
margNames = {'Direction', 'Duration','Condition-independent','D/D Interaction'};
margColours =  [0.5 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1];

[W,V,whichMarg] = dpca(FR_dPCA, 15, 'combinedParams', combinedParams);
explVar = dpca_explainedVariance(FR_dPCA, W, V, ...
    'combinedParams', combinedParams);
variance=explVar.cumulativeDPCA;

NdimdPCA=15;%find(variance>=80,1,'First') 
if isempty(NdimdPCA)
    NdimdPCA=15;
end
% Variance explained by direction component in %
W=W(:,1:NdimdPCA);
V=V(:,1:NdimdPCA);
whichMarg=whichMarg(:,1:NdimdPCA);

Dur_var_exp=[0 0 0];
if sum(whichMarg==1)>0
    Dur_var_exp(1)=sum(explVar.componentVar(whichMarg==1));
end
if sum(whichMarg==2)>0
    Dur_var_exp(2)=sum(explVar.componentVar(whichMarg==2));
end
if sum(whichMarg==3)>0
    Dur_var_exp(3)=sum(explVar.componentVar(whichMarg==3));
end
disp(['Variance explained by duration Component = ',num2str(Dur_var_exp)])

dpca_plot(FR_dPCA, W, V, @dpca_plot_default, ...
    'explainedVar', explVar,...
    'time', linspace(0,1,size(FR_dPCA,4)),                        ...
    'timeEvents', 0,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg);

end
