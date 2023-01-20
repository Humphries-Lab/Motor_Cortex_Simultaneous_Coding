function decoding_movement_direction_all_sessions(session,Area,threshold,Ndir,k_fold,nrep,shuffle,do_plot)
%step 1: extract the information from all trials
%step 2: each movement is going to have an ID. Divide the total number of
%trials into Ndir classes. Make sure that all classes have the same number
%of samples
%step 3: for those selected movements, split the data to do CV
% step 4: train the naive bayes model inside CV
% step 5: test using 1 class prediction accuracy and the estimation of the
% angle (Check with Mark that the last one makes sense)

colour_dir=hsv(Ndir);
median_error=zeros(size(session,2),1);
std_error=zeros(size(session,2),1);
Acc=zeros(size(session,2),1);
acc_error_mean=zeros(size(session,2),1);
acc_error_std=zeros(size(session,2),1);
for isession=1:size(session,2)
    disp(['Starting session ' num2str(isession)])
    %% Step 1
    [angle_dir,i_dir,neural_mov]=extract_movement_for_decoding(session{isession}, Area{isession},Ndir);
    
    % %% Step 2
    Nsample=zeros(Ndir,1);
    
    for i=1:Ndir
        Nsample(i)=sum(i_dir==i);
    end
    
    New_N=min(Nsample);
    selected_angle=[];
    selected_class=[];
    selected_neural=[];
    
    neural_mov=squeeze(mean(neural_mov,2))';
    %selected_class=i_dir;
    xangle=zeros(1,Ndir);

    for i=1:Ndir
        idx_i=find(i_dir==i);
        %randon permutation to break temporal patterns throughout the session
        p_tmp = randperm(numel(idx_i));
        selected_idx=idx_i(p_tmp(1:New_N));
        selected_angle=[selected_angle;angle_dir(selected_idx)'];
        selected_class=[selected_class;i_dir(selected_idx)'];
        xangle(i)=mean(angle_dir(selected_idx));
        
        selected_neural=[selected_neural;neural_mov(selected_idx,:)];
        
        
    end
    
    %% Step 3
    %%%% repeat CV 100 times to make a robust statistic analysis
    %always fun
    
    %% Shuffle test: permute selected_class to estimate change levels
    if shuffle
        selected_class=selected_class(randperm(numel(selected_class)));
    end
    counter=1;
    Acc_per_cv=zeros(k_fold,nrep);
    for r=1:nrep
        
        % Partition data (stratify per class)
        c = cvpartition(selected_class,'KFold',k_fold);
        %preprocess all data
        average_cond_1=selected_neural./repmat(range(selected_neural)+5,size(selected_neural,1),1);
        
        predictedLabels=[];
        trueLabels=[];
        predictedThetas=[];
        trueThetas=[];
        predictedThetas_bin=[];
        
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
            for i=1:Ndir
                xangle_training(i)=mean(training_angles(training_class==i));
                %         Nsample_training(i)=sum(training_class==i);
                %         Nsample_testing(i)=sum(testing_class==i);
            end
            %         subplot(2,2,3)
            %         plot(Nsample_training)
            %         hold on
            %         plot(Nsample_testing)
            
            [coeff,score,~,~,variance]=pca(traning_neural);
            Ndim(k)=find(cumsum(variance)>threshold,1,'First');
            scores_test=(testing_neural-repmat(mean(traning_neural),size(testing_neural,1),1))*coeff(:,1:Ndim);
           
            if do_plot && isession==11 && ~shuffle
                if k==1 && r==1
                    subplot(2,3,1)
                     hold on
                    for i=1:Ndir
                        plot(score(training_class==i,1),score(training_class==i,2),'.','Color',colour_dir(i,:),'MarkerSize',12)
                        plot(scores_test(testing_class==i,1),scores_test(testing_class==i,2),'.','MarkerSize',18,'Color',colour_dir(i,:)/2)
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
        Acc_per_cv(:,r)=Acc_per_fold;

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
    
    if do_plot && isession==11 && ~shuffle
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
    M1=strcmp(Area,'M1');
    subplot(2,3,6)
    hold on
    errorbar(1:sum(M1),median_error(M1),std_error(M1),'.m')
    errorbar((1:sum(~M1))+sum(M1),median_error(~M1),std_error(~M1),'.b')
    
    [mean(median_error(M1)) mean(median_error(~M1)) mean(median_error_angle_bin(M1)) mean(median_error_angle_bin(~M1))]
    box off
    xlabel('Session')
    ylabel('Absolute angle error [degrees]')
    
    
    subplot(2,3,5)
    hold on
    errorbar(1:sum(M1),acc_error_mean(M1),acc_error_std(M1),'.m')
    [mean(acc_error_mean(M1)) mean(acc_error_mean(~M1))]
    errorbar((1:sum(~M1))+sum(M1),acc_error_mean(~M1),acc_error_std(~M1),'.b')
    box off
    hold on
    xlabel('Session')
    ylabel('Average Accuracy per fold')
    ylim([0 1])
    
else
    colour_shuffle=[0.5 0.5 0.5];
    subplot(2,3,6)
    hold on
    errorbar(1:size(session,2),median_error,std_error,'.','Color',colour_shuffle)
    
    
    subplot(2,3,5)
    hold on
    errorbar(1:size(session,2),acc_error_mean,acc_error_std,'Color',colour_shuffle)
    plot([1 size(session,2)],[1/8 1/8],'k')
end

end

function [angle_dir,class_dir,neural_mov]=extract_movement_for_decoding(session, Area,Ndir)
%colour_dir=hsv(Ndir);
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

ms=nanmedian(ISI); %round(ms*sqrt(12))
t_1=-0.2;
t_2=0;
from_to=[0.2 1.2];
[neural_mov,angle_dir]=neural_data_per_duration(session,Area,ms,t_1,t_2,event,from_to,0);
class_dir=ceil(Ndir*(angle_dir+pi)/(2*pi));

end


function [label,theta,mean_angle_bin]=predict_from_bayes(sample,Mdl,xangle)
debugging=0;
[label,Posterior] = predict(Mdl,sample);
[x,y] = pol2cart(xangle,Posterior);
mean_angle_bin=xangle(label)';
[theta,~] = cart2pol(sum(x,2),sum(y,2));
if debugging
    figure
    for i=1:numel(theta)
    mean_angle_bin(i)*180/pi
    theta(i)*180/pi
    polarplot(xangle,ones(numel(xangle),1))
    hold on
    polarplot(xangle,Posterior(1,:))
    polarplot([theta(i) theta(i)],[0 1],'k')
    polarplot([mean_angle_bin(i) mean_angle_bin(i)],[0 0.8],'r')
    pause
    cla
    end
end

end