function [error_all,error_by_dir_shift]=decoding_movement_duration_chance(score,idx_dir,idx_duration,ref_bin,Ndim,kneigh)
%% approach to decode duration
%% 1) select half of the trials (we need examples of all durations and directions)
%% 2) train the decoder
%% 3) test in the other half


Ndir=max(idx_dir);
error_all=[];
test_bins=1:4;
test_bins(ref_bin)=[];
traj_ref=sum((idx_dir==1)& (idx_duration==ref_bin));
Ntest_bins=numel(test_bins);
for i_bin=1:Ntest_bins
     bin_i=test_bins(i_bin);
traj_length(i_bin)=sum((idx_dir==1)& (idx_duration==bin_i));
end

fraction=zeros(Ndir,Ntest_bins);
estimated_duration=zeros(Ndir,Ntest_bins);
start=1;%round(t_1*-1000);
counter=1;
for i_dir=1:Ndir
    idx=find((idx_dir==i_dir)& (idx_duration==ref_bin));
    %reference/ longest/slowest
    x1=score(idx(1:end),1:Ndim);
    clear idx
    
    for i_bin=1:Ntest_bins
        bin_i=test_bins(i_bin);
        for j_dir=1:Ndir
        idx=find((idx_dir==j_dir) & (idx_duration==bin_i));
        x2=score(idx(start:end),1:Ndim);
        %% shuffle
        %x2=x2(randperm(size(x2,1)),:);
        [fraction(i_dir,i_bin,j_dir),idx_tmp{i_bin}(:,i_dir)]=compute_fraction_speed_general_v3(x1,x2,kneigh);
        estimated_duration(i_dir,i_bin,j_dir)=traj_ref/fraction(i_dir,i_bin,j_dir);
        clear idx
        end
        if i_dir~=j_dir
        error_all=[error_all;abs(squeeze(estimated_duration(i_dir,i_bin,:))-traj_length(i_bin))];
        end
        error_by_dir=abs(squeeze(estimated_duration(i_dir,i_bin,:))-traj_length(i_bin));
        error_by_dir_shift(:,counter)=circshift(error_by_dir,round(Ndir/2)-i_dir);
        counter=counter+1;

    end
    
    
end

end

function fraction=compute_fraction_speed_general(xref,xsample)
S_ref=size(xref,1);
S_sample=size(xsample,1);
max_diff=abs(S_ref-S_sample);
idx=nan(S_sample,1);
for i=1:S_sample
    %find only in a time bubble of size (at least) max_diff
    if  i<=max_diff
        from=i;
        to=round(S_ref/2);
        if to>S_ref
            to=S_ref;
        end
    elseif i>=S_ref-max_diff
           from=S_ref-round(S_ref/2);
           to=S_ref; 
    else
        from=i-max_diff;
        to=i+max_diff;
    end

%     x1, the longest, is the reference
    xref_tmp=nan(size(xref));
    xref_tmp(from:to,:)=xref(from:to,:);
    dist_i=pdist2(xref_tmp,xsample(i,:));

    [~,idx(i)]=min(dist_i);
%     [~,idx_tmp]=mink(dist_i,10);
%     idx(i)=mean(idx_tmp);

%     if i>1 && abs(idx(i)-idx(i-1))>max_diff && i>max_diff
%         [i from to abs(idx(i)-idx(i-1)) max_diff]
%         plot(idx)
%       pause
%          
%     end
%    [i from to idx(i)]    
%     subplot(2,1,1)
%     plot(xref_tmp(from:to,1),xref_tmp(from:to,2),'.-r')
%     hold on
%     plot(xsample(1:i,1),xsample(1:i,2),'Color','b')
%     hold off
%     [from to]
%     pause(0.01)

end

% subplot(2,1,2)
% plot(idx)
% hold on 
%  plot([0 S_sample],[0 S_sample])
%  pause
%% linear regression
%ft = fittype('a*x+b');
%f1 = fit((1:i)',idx(1:i),ft);
% 
       plot(idx)
       hold on 
%        plot(1:i,(1:i)*f1.a+f1.b)

%fraction=f1.a;
 pause

%% just the mean
fraction=mean(diff(idx));
%pause
end

function fraction=compute_fraction_speed_general_v2(xref,xsample)
S_ref=size(xref,1);
S_sample=size(xsample,1);
max_diff=abs(S_ref-S_sample);
idx=nan(S_sample,1);
for i=1:S_sample
    %find only in a time bubble of size (at least) max_diff
    from=1;
    to=size(xref,1);

%     x1, the longest, is the reference
    xref_tmp=nan(size(xref));
    xref_tmp(from:to,:)=xref(from:to,:);
    dist_i=pdist2(xref_tmp,xsample(i,:));

    [~,idx(i)]=min(dist_i);
   
%     if i>1 && abs(idx(i)-idx(i-1))>max_diff && i>max_diff
%         [i from to abs(idx(i)-idx(i-1)) max_diff]
%         plot(idx)
%       pause
%          
%     end
%    [i from to idx(i)]    
%     subplot(2,1,1)
%     plot(xref_tmp(from:to,1),xref_tmp(from:to,2),'.-r')
%     hold on
%     plot(xsample(1:i,1),xsample(1:i,2),'Color','b')
%     hold off
%     [from to]
%     pause(0.1)

end

% subplot(2,1,2)
plot(idx)
hold on 
 plot([0 S_sample],[0 S_sample])
 pause
%% linear regression
%ft = fittype('a*x+b');
%f1 = fit((1:i)',idx(1:i),ft);
% 
%        plot(idx)
%        hold on 
%        plot(1:i,(1:i)*f1.a+f1.b)

%fraction=f1.a;
% pause

%% just the mean
diff_idx=diff(idx);
diff_idx(abs(diff_idx)>S_ref/2)=[];
plot(diff_idx)
pause
fraction=mean(diff_idx);
%fraction=mean(diff(idx))
%pause
end

function [fraction,idx]=compute_fraction_speed_general_v3(xref,xsample,k)
S_ref=size(xref,1);
S_sample=size(xsample,1);
idx=nan(S_sample,1);
% from=1;
% to=S_ref;
for i=1:S_sample
    %find only in a time bubble of size (at least) max_diff
   
%     x1, the longest, is the reference
    xref_tmp=nan(size(xref));
    if i<200
        from=1;
        to=round(S_ref/2);
    elseif i>S_ref-200
        from=round(S_ref/2);
        to=S_ref;
    else
        from=1;
        to=S_ref;
    end
    xref_tmp(from:to,:)=xref(from:to,:);
    dist_i=pdist2(xref_tmp,xsample(i,:));

    %[~,idx(i)]=min(dist_i);
    [~,idx_tmp]=mink(dist_i,k);
    idx(i)=median(idx_tmp);
%     if i>1 && abs(idx(i)-idx(i-1))>max_diff && i>max_diff
%         [i from to abs(idx(i)-idx(i-1)) max_diff]
%         plot(idx)
%       pause
%          
%     end
%    [i from to idx(i)]    
%     subplot(2,1,1)
%     plot(xref_tmp(from:to,1),xref_tmp(from:to,2),'.-r')
%     hold on
%     plot(xsample(1:i,1),xsample(1:i,2),'Color','b')
%     hold off
%     [from to]
%    pause(0.01)

end

% subplot(2,1,2)
% plot(idx)
% hold on 
%  plot([0 S_sample],[0 S_sample])
%  pause
%% linear regression
ft = fittype('a*x+b');
f1 = fit((1:i)',idx(1:i),ft,'StartPoint',[0 0]);
% 
%        plot(idx)
%        hold on 
%        plot(1:i,(1:i)*f1.a+f1.b)

fraction=f1.a;
% pause

%% just the mean
%fraction=mean(diff(idx));
%pause
end
