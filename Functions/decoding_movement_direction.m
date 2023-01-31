function [Acc_per_cv,error_angle,error_angle_bin]=decoding_movement_direction(selected_class,selected_angle,selected_neural,threshold,colour_dir,k_fold,Nrep,shuffle,do_plot)
Ndir=max(selected_class);
Acc_per_cv=zeros(k_fold,Nrep);
%% Shuffle test: permute selected_class to estimate change levels
if shuffle
    selected_class=selected_class(randperm(numel(selected_class)));
end

%predictedLabels=[];
% trueLabels=[];
% predictedThetas=[];
% trueThetas=[];
% predictedThetas_bin=[];
trueLabels=nan(50000,1);
predictedThetas=nan(50000,1);
trueThetas=nan(50000,1);
predictedThetas_bin=nan(50000,1);
predictedLabels=nan(50000,1);

counter=1;
counterend=0;
for i_rep=1:Nrep
    
    % Partition data (stratify per class)
    c = cvpartition(selected_class,'KFold',k_fold);
    % preprocess all data
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
        
        if do_plot
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
        
        counterend=size(predictedLabels_iter,1)+counterend;
        
        predictedLabels(counter:counterend)=predictedLabels_iter;
        trueLabels(counter:counterend)=testing_class;
        predictedThetas(counter:counterend)=thetas_iter;
        predictedThetas_bin(counter:counterend)=angle_bin;
        trueThetas(counter:counterend)=testing_angles;

        %predictedLabels=[predictedLabels;predictedLabels_iter];
%         trueLabels=[trueLabels;testing_class];
%         predictedThetas=[predictedThetas;thetas_iter];
%         predictedThetas_bin=[predictedThetas_bin;angle_bin];
%         trueThetas=[trueThetas;testing_angles];
        
        counter=counterend+1;
        
        
        C = confusionmat(testing_class,predictedLabels_iter);
        Acc_per_fold(k)=sum(diag(C))/numel(testing_class);
        error_angle_iter=abs(angdiff(thetas_iter,testing_angles)*180/pi);
        mean_error_iter(k)=median(error_angle_iter);
        std_error_iter(k)=std(error_angle_iter);
        
    end
    Acc_per_cv(:,i_rep)=Acc_per_fold;
    
end
predictedLabels(counter:end)=[];
trueLabels(counter:end)=[];
predictedThetas(counter:end)=[];
trueThetas(counter:end)=[];
predictedThetas_bin(counter:end)=[];

error_angle=abs(angdiff(predictedThetas,trueThetas)*180/pi);
error_angle_bin=abs(angdiff(predictedThetas_bin,trueThetas)*180/pi);

if  do_plot
    
    C_total = confusionmat(trueLabels,predictedLabels);
    Acc=sum(diag(C_total))/numel(trueLabels);
    
    fig1=figure;
    
    cm = confusionchart(trueLabels,predictedLabels,'Normalization','row-normalized');
    NormalizedValues=cm.NormalizedValues;
    
    close(fig1)
    
    subplot(2,3,2)
    imagesc(NormalizedValues)
    colormap gray
    title(['Acc = ' num2str(Acc)])
    colorbar
    caxis([0 1])
    box off
    xlabel('Predicted direction')
    ylabel('True direction')
    box off
end
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