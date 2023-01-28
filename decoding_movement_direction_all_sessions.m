function decoding_movement_direction_all_sessions(Sessions,Areas,threshold,Ndir,k_fold,nrep,shuffle)
%% decoding_movement_direction_all_sessions predicts the movement direction bin and angle from the neural activity
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
% Ndir: number of direction to bin the movements
%
% k_fold: number of folds for Cross-Validation (CV)
%
% nrep: Number of repetitions of CV
%
% shuffle: 1- shuffle the labels/direction bins
%          0- don't shuffle
%
%% 24/05/2023
% Andrea Colins Rodriguez

%% Overview
%   step 1: extract the information from all trials
%   step 2: each movement has an ID. Divide the total number of
%   trials into Ndir bins. Select a random sample of the movements so that all bins have the same number
%   of samples.
%   step 3: for those selected movements, split the data to do CV
%   step 4: train the naive bayes model inside CV
% step 5: test using 1 class prediction accuracy and the estimation of the
% angle 

colour_dir=hsv(Ndir);
median_error=zeros(size(Sessions,2),1);
median_error_angle_bin=zeros(size(Sessions,2),1);
std_error_bin=zeros(size(Sessions,2),1);
std_error=zeros(size(Sessions,2),1);
Acc=zeros(size(Sessions,2),1);
acc_error_mean=zeros(size(Sessions,2),1);
acc_error_std=zeros(size(Sessions,2),1);
Colour_M1=[85 30 116]./256;
Colour_PMd=[89 156 153]./256;

if shuffle==0
    figure
end

for isession=1:size(Sessions,2)
    disp(['Starting decoding direction for recording ' num2str(isession)])
    %% Step 1
    [angle_dir,idx_dir,neural_mov]=extract_movement_for_decoding(Sessions{isession}, Areas{isession},Ndir);
    
    %% Step 2
    Nsample=zeros(Ndir,1);
    
    for i_dir=1:Ndir
        Nsample(i_dir)=sum(idx_dir==i_dir);
    end
    
    New_N=min(Nsample);
    selected_angle=nan(New_N*Ndir,1);
    selected_class=nan(New_N*Ndir,1);
    selected_neural=nan(New_N*Ndir,size(neural_mov,1));
    
    % transform FR to a matrix [n mov, units]
    neural_mov=squeeze(mean(neural_mov,2))';
    xangle=zeros(1,Ndir);

    for i_dir=1:Ndir
        idx_i=find(idx_dir==i_dir);
        %randon permutation to break temporal patterns throughout the session
        p_tmp = randperm(numel(idx_i));
        selected_idx=idx_i(p_tmp(1:New_N));
        
        selected_angle((i_dir-1)*New_N+1:New_N*i_dir)=angle_dir(selected_idx)';
        selected_class((i_dir-1)*New_N+1:New_N*i_dir)=idx_dir(selected_idx)';
        selected_neural((i_dir-1)*New_N+1:New_N*i_dir,:)=neural_mov(selected_idx,:);
        
        xangle(i_dir)=mean(angle_dir(selected_idx));
        
    end
    
    %% Step 3
    %%%% repeat CV Nrep times to make a robust statistic analysis

    counter=1;
    Acc_per_cv=zeros(k_fold,nrep);
    
    %% Shuffle test: permute selected_class to estimate change levels
    if shuffle
        selected_class=selected_class(randperm(numel(selected_class)));
    end
    
            predictedLabels=[];
        trueLabels=[];
        predictedThetas=[];
        trueThetas=[];
        predictedThetas_bin=[];
        
    for i_rep=1:nrep
        
        % Partition data (stratify per class)
        c = cvpartition(selected_class,'KFold',k_fold);
        %preprocess all data
        average_cond_1=selected_neural./repmat(range(selected_neural)+5,size(selected_neural,1),1);
        

        
        Acc_per_fold=zeros(k_fold,1);
        mean_error_iter=zeros(k_fold,1);
        std_error_iter=zeros(k_fold,1);
        Ndim=zeros(k_fold,1);
        
        for k=1:k_fold
            %% Step 4
            idx_training = training(c,k);
            idx_testing = test(c,k);
            traning_neural=average_cond_1(idx_training,:);
            testing_neural=average_cond_1(idx_testing,:);
            training_class=selected_class(idx_training);
            testing_class=selected_class(idx_testing);
            training_angles=selected_angle(idx_training);
            testing_angles=selected_angle(idx_testing);
            xangle_training=zeros(1,Ndir);
            
            for i_dir=1:Ndir
                xangle_training(i_dir)=mean(training_angles(training_class==i_dir));
            end

            [coeff,score,~,~,variance]=pca(traning_neural);
            Ndim(k)=find(cumsum(variance)>threshold,1,'First');
            scores_test=(testing_neural-repmat(mean(traning_neural),size(testing_neural,1),1))*coeff(:,1:Ndim);
           
            if isession==11 && ~shuffle
                if k==1 && i_rep==1
                    subplot(2,3,1)
                     hold on
                    for i_dir=1:Ndir
                        plot(score(training_class==i_dir,1),score(training_class==i_dir,2),'.','Color',colour_dir(i_dir,:),'MarkerSize',12)
                        plot(scores_test(testing_class==i_dir,1),scores_test(testing_class==i_dir,2),'.','MarkerSize',18,'Color',colour_dir(i_dir,:)/2)
                    end
                end
            end
            %train Naive Bayes
            Mdl = fitcnb(score(:,1:Ndim),training_class);
            
            
            %% Step 5
            % predict classes and angles on the test set
            
            [predictedLabels_iter,thetas_iter,angle_bin]=predict_from_bayes(scores_test,Mdl,xangle_training);
            predictedLabels=[predictedLabels;predictedLabels_iter];
            trueLabels=[trueLabels;testing_class];
            predictedThetas=[predictedThetas;thetas_iter];
            predictedThetas_bin=[predictedThetas_bin;angle_bin];
            trueThetas=[trueThetas;testing_angles];
            
            C = confusionmat(testing_class,predictedLabels_iter);
            Acc_per_fold(k)=sum(diag(C))/numel(testing_class);
            error_angle_iter=abs(angdiff(thetas_iter,testing_angles)*180/pi);
            mean_error_iter(k)=median(error_angle_iter);
            std_error_iter(k)=std(error_angle_iter);
            counter=counter+1;
        end
        Acc_per_cv(:,i_rep)=Acc_per_fold;

    end
    
    acc_error_mean(isession)=mean(Acc_per_cv(:));
    acc_error_std(isession)=std(Acc_per_cv(:));
    error_angle=abs(angdiff(predictedThetas,trueThetas)*180/pi);
    error_angle_bin=abs(angdiff(predictedThetas_bin,trueThetas)*180/pi);
    
    median_error_angle_bin(isession)=median(error_angle_bin);
    std_error_bin(isession)=std(error_angle_bin);
    median_error(isession)=median(error_angle);
    std_error(isession)=std(error_angle);
    
    C_total = confusionmat(trueLabels,predictedLabels);
    Acc(isession)=sum(diag(C_total))/numel(trueLabels);
    fig1=figure;
    
    cm = confusionchart(trueLabels,predictedLabels,'Normalization','row-normalized');
    NormalizedValues=cm.NormalizedValues;
    
    close(fig1)
    

    
    if  isession==11 && ~shuffle
        subplot(2,3,2)
        imagesc(NormalizedValues)
        colormap gray
        title(['Acc = ' num2str(Acc(isession))])
        colorbar
        caxis([0 1])
        box off
        xlabel('Predicted direction')
        ylabel('True direction')
        box off
    end
    
end


%% Summary for all sessions
if ~shuffle
    M1=strcmp(Areas,'M1');
    subplot(2,3,6)
    hold on
    errorbar(1:sum(M1),median_error(M1),std_error(M1),'.','Color',Colour_M1)
    errorbar((1:sum(~M1))+sum(M1),median_error(~M1),std_error(~M1),'.','Color',Colour_PMd)
    box off
    xlabel('Session')
    ylabel('Absolute angle error [^o]')
    
    
    disp('---------------------------------------')
    disp(['Median error predicted bin = ' num2str(mean(median_error_angle_bin),2)])
    disp(['Median error predicted angle = ' num2str(mean(median_error),2)])
    disp('---------------------------------------')

    
    subplot(2,3,5)
    hold on
    errorbar(1:sum(M1),acc_error_mean(M1),acc_error_std(M1),'.','Color',Colour_M1) 

    errorbar((1:sum(~M1))+sum(M1),acc_error_mean(~M1),acc_error_std(~M1),'.','Color',Colour_PMd)
    text(1,0.4,['Mean M1 =' num2str(mean(acc_error_mean(M1)),3)],'FontSize',8) 
    text(10,0.4,['Mean PMd =' num2str(mean(acc_error_mean(~M1)),3)],'FontSize',8)
    box off
    hold on
    xlabel('Recording')
    ylabel('Average Accuracy per fold')
    ylim([0 1])
    
else
    colour_shuffle=[0.5 0.5 0.5];
    subplot(2,3,6)
    hold on
    errorbar(1:size(Sessions,2),median_error,std_error,'.','Color',colour_shuffle)
    
    
    subplot(2,3,5)
    hold on
    errorbar(1:size(Sessions,2),acc_error_mean,acc_error_std,'Color',colour_shuffle)
    plot([1 size(Sessions,2)],[1/8 1/8],'k')
end

end

function [angle_dir,class_dir,neural_mov]=extract_movement_for_decoding(session, Area,Ndir)
%% extract_movement_for_decoding extracts all movements and their respective
%% *preparatory* neural activity.
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
% OUTPUTS
%
% angle_dir: angle of each movment [rad]
% 
% class_dir: number of class/ direction bin of each movement
%
% neural_mov: Preparatory activity (upto 200 ms before mov onset) of each
% movement
%
% 27/05/2023
% Andrea Colins Rodriguez

event=2;
if strcmp(Area,'PMd')
    load(session,'PMd','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(PMd.units,startt,endt);
else
    load(session,'M1','trial_table2')
    startt=trial_table2(1,1);
    endt=trial_table2(size(trial_table2,1),22);
    ISI=compute_ISI(M1.units,startt,endt);
end

sigma_filter=round(median(ISI));

t_from=-0.2;
t_upto=0;
from_to=[0.2 1.2];

[Neural_info,Mov_params]=neural_data_per_duration(session,Area,sigma_filter,t_from,t_upto,event,from_to);
neural_mov=Neural_info.FR;
angle_dir=Mov_params.direction;
class_dir=ceil(Ndir*(angle_dir+pi)/(2*pi));

end


function [pred_bin,pred_theta,mean_angle_bin]=predict_from_bayes(score_test,Mdl,xangle)
%% predict_from_bayes predicts the direction bin corresponding to the neural activity in sample
%
% INPUTS
%
% score_test: projection of the test data into the subspace. Rows are
% samples, columns are neurons
%
% Mdl: Bayes Model trained with the training dataset
%
% xangle: Mean angle of each direction bin used during training
%
% OUTPUTS
% 
% pred_bin: predicted direction bin for each score_test 
%
% pred_theta: predicted angle movement for each score_test
%
% mean_angle_bin: mean angle direction of pred_bin
%
% 27/05/2023
% Andrea Colins Rodriguez


[pred_bin,Posterior] = predict(Mdl,score_test);
[x,y] = pol2cart(xangle,Posterior);
mean_angle_bin=xangle(pred_bin)';
[pred_theta,~] = cart2pol(sum(x,2),sum(y,2));
end