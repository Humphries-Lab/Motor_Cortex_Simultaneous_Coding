function r=embedding_dimensions_upper_bound(Session,Area,Ndir,Nbins,threshold,t_from,t_upto,edges_dur_bin)
%% embedding_dimensions_upper_bound computes the upper bound of the coefficient of determination
% between the trajectories within the same bin direction and duration from
% Session and Area
% This function also saves the results (r) on a file for later analyses
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
% r= coefficient of determination for each pair of subsampled trajectories
%
%
% 24/05/2023
% Andrea Colins Rodriguez



%colour_dir=hsv(Ndir);

load(Session,Area,'trial_table2','cont')
load(['../Output_files/PCA_' Session(1:end-4) '_' Area '.mat'],'coeffs','delete_units','normalisation','mean_norm','variance')
ndim=find(cumsum(variance)>threshold,1,'First');

startt=trial_table2(1,1);
endt=trial_table2(size(trial_table2,1),22);

if strcmp(Area,'PMd')
    neural_data=PMd.units;
else
    neural_data=M1.units;
end

ISI=compute_ISI(neural_data,startt,endt);
sigma_filter=round(median(ISI));


dur_binsize=(edges_dur_bin(2)-edges_dur_bin(1));
t_upto=edges_dur_bin(1:Nbins)+dur_binsize+t_upto;

% mean_nsample=mean(nsamples_condition);
r=[];
Hdist_same=nan(Nbins,Ndir);

% Extract movements and related neural activity for each duration bin
for i_dur=1:Nbins
    
    current_dur_bin=[edges_dur_bin(i_dur) edges_dur_bin(i_dur+1)]; % select movements in this range of duration only
    
    [Neural_info,Mov_params]=neural_data_per_duration(cont,trial_table2,neural_data,sigma_filter,t_from,t_upto(i_dur),current_dur_bin);
    
    
    
    % bin by direction
    direction1=ceil(Ndir*(Mov_params.direction+pi)/(2*pi));
    
    
    for i_dir=1:Ndir
        
        
        %% calculate upper bound of similarity between single trajectories
        
        %1) select trajectories that are very similar
         close_dir=direction1==i_dir;

       [r_tmp,Hdist_same(i_dur,i_dir)]=upper_bound(close_dir,Neural_info,delete_units,normalisation,mean_norm,coeffs,ndim);
        r=[r,r_tmp];
           
    end
    
end
save(['../Output_files/PCA_' Session(1:end-4) '_' Area '.mat'],'Hdist_same','r','-append')
end

function [r,Hdist_same]=upper_bound(close_dir,Neural_info,delete_units,normalisation,mean_norm,coeffs,ndim)

% within bin
idx_sim=find(close_dir);

total_samples=numel(idx_sim);
nsamp=floor(total_samples/2);

tmp=nan(1,nsamp);
r=nan(1,nsamp);

for i_rand=1:nsamp
    idx_sim=idx_sim(randperm(total_samples));
    
    FR_sim1=mean(Neural_info.FR(:,:,idx_sim(1:nsamp)),3);
    FR_sim2=mean(Neural_info.FR(:,:,idx_sim(end-nsamp+1:end)),3);
    
    % 2) project into subspace
    FR_sim1(delete_units,:)=[];
    FR_sim1=FR_sim1'./repmat(normalisation,size(FR_sim1,2),1);
    FR_sim1=FR_sim1-mean_norm;
    proj_FR1=FR_sim1*coeffs;
    
    FR_sim2(delete_units,:)=[];
    FR_sim2=FR_sim2'./repmat(normalisation,size(FR_sim2,2),1);
    FR_sim2=FR_sim2-mean_norm;
    proj_FR2=FR_sim2*coeffs;
    
    %3) calculate R2
    r(i_rand)=R2(proj_FR1(:,1:ndim),proj_FR2(:,1:ndim));
    
    % 4) calculate the haussdorf distance
    recurrence=pdist2(proj_FR1(:,1:ndim),proj_FR2(:,1:ndim));
    tmp(i_rand)=max([min(recurrence),min(recurrence,[],2)']);
end



Hdist_same=mean(tmp);
end

function r=R2(X,Y)
%% R2 Computes the coefficient of determination between n-dimensional
%% trajectories X and Y
X2=reshape(X,numel(X),1);
Y2=reshape(Y,numel(Y),1);
ymean=mean(Y2);
Sres=sum((Y2-X2).^2);
Stot=sum((Y2-ymean).^2);
r=1-Sres/Stot;
end