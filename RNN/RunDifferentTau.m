function RunDifferentTau
%% Script to test MCtx RNN SOC model response to variation in time constant of the RNN
% Current version: Changing the tau of the RNN temporally scale the
% trajectories preserving their shape.

% Mark Humphries & Andrea Colins Rodriguez 28/07/22

clearvars; close all

run MCtx_RNN_SOC_common_parameters


%% parameters unqiue to this set of simulations
pars.rampTimeMin = 1000;
pars.rampTimeMax = 2000;
pars.rampTimeStep = 500;

pars.Tmax = pars.duration(1) + pars.rampTimeMax + max(pars.ts)/pars.dt*10; % set same max duration, to cope with longest ramp

pars.ts2=pars.ts*0.8;
pars.ts3=pars.ts*1.2;

%% Simulate network to actually see the trajectories
ts=[pars.ts2 pars.ts pars.ts3];

tSteps = round(pars.Tmax ./ pars.dt);   % how many time-steps

% set up parameter range to explore
rampTimeRange = pars.rampTimeMin:pars.rampTimeStep:pars.rampTimeMax;

% storage
nRampTimes = numel(rampTimeRange);
Data.temporalPCs = zeros(nRampTimes,tSteps,pars.nPCs);
rStore = zeros(nRampTimes * tSteps,Net.N);

%find subspace for all directions first
fig2=figure;
[coeffs,mean_data,~,initialProjection,threshold]=find_dir_subspace(pars,optimW,Net);

Ndir=10;
colour_dir=hsv(Ndir);
colour_dur=plasma(nRampTimes);
pars.angle = linspace(0,2*pi,Ndir+1); % want 11 of these... % different rotations of that vector: matching different reach angles in discrete analysis
pars.angle(end) = [];               % to get 10 reach angles, equally separated

for iAngle=1:Ndir
    rotatedProjection = rotate_n_dimensional_vector(initialProjection,pars.angle(iAngle)); % rotate current projection vector of ramp
    output_dur=zeros(nRampTimes,1);
    output_speed=zeros(nRampTimes,1);
    
    for iRampTime = 1:nRampTimes
        
        %% change the value of the time constant
        pars.ts=ts(iRampTime);
        
        % storage
        a = zeros(Net.N,tSteps);
        r = zeros(Net.N,tSteps);
        
        % create ramped input for this duration
        duration = [pars.duration(1) pars.duration(1) + rampTimeRange(2)];
        tInput = [round(duration(1) ./pars.dt), round((duration(2) + pars.fall_duration) ./ pars.dt)];    % how long is the input in timesteps?
        
        % generate ramped input - to maximum value
        % from Hennequin et al 2014: exponential rise to max, then exponential
        % decay
        
        [I_ramp] = makeRampInput(pars.dt,duration(end)-duration(1),pars.fall_duration,pars.ts_rise,pars.ts_fall,pars.ramp_max);
        
        
        
        % run simulation
        for iT = 2:tSteps
            % background input - set up here to add noise later
            I = zeros(Net.N,1) + pars.I_background;
            
            % add ramped input
            if iT > tInput(1) && iT <= tInput(2)
                I = I + I_ramp(iT - tInput(1)) .* rotatedProjection;
                
            end
            % update activity
            a(:,iT) = a(:,iT-1) + pars.dt * (-a(:,iT-1) + I + optimW * r(:,iT-1)) ./ pars.ts;
            
            % update output
            r(:,iT) = neuron_output(a(:,iT),pars.output_type,pars.output_arg1,pars.output_arg2);
            
        end
        
        rActive = r(:,1:tSteps);        % indexing here to later allow for turning off early
        
        
        % store outputs, concatenated
        rStore(1+(iRampTime-1)*tSteps:iRampTime*tSteps,:) = rActive';  % columns are neurons
        
    end
    
    %% process all simulation results
    % compute PCA on all concatenated outputs
    %[Vectors,temporalPCs,EigVals] = pca(rStore);   % do PCA on all neuron outputs
    
    rStore=rStore-repmat(mean_data,size(rStore,1),1);
    temporalPCs=rStore*coeffs;
    
    figure(fig2)
    
    
    upto=2700;
    % store top N PCs, for each Duration
    for iRampTime= 1:nRampTimes
        Data.temporalPCs(iRampTime,:,:) = temporalPCs(1+(iRampTime-1)*tSteps:iRampTime*tSteps,1:pars.nPCs);
        if iAngle==1
            subplot(2,3,3)
            hold on
        plot3(Data.temporalPCs(iRampTime,1000:upto,1),Data.temporalPCs(iRampTime,1000:upto,2),Data.temporalPCs(iRampTime,1000:upto,3),'Color',colour_dur(iRampTime,:),'LineWidth',2)
        end
        %% store data for angle only for one tau
    if iRampTime==2
           
        Data.temporalPCs_angle(iAngle,:,:) = temporalPCs(1+(iRampTime-1)*tSteps:iRampTime*tSteps,1:pars.nPCs);
    end
    
    end
    
    
    %% Estimate output duration using recurrence
    for iRampTime= 1:nRampTimes
        
        % using length of the vector
        norms=vecnorm(squeeze(Data.temporalPCs(iRampTime,:,1:pars.nPCs))');
        
        idx_start=find(norms>=threshold,1,'First');
        idx_end=find(norms>=threshold,1,'Last');
        
        
        output_dur(iRampTime)=(idx_end-idx_start).*pars.dt;
        output_speed(iRampTime)=mean(sqrt(sum(diff(squeeze(Data.temporalPCs(iRampTime,idx_start:idx_end,1:pars.nPCs))).^2,2)));
        
        
        
    end
    
    subplot(2,3,6),
    hold on
    plot(ts,output_dur,'Color',colour_dir(iAngle,:))
    
end

%% finish the plots
subplot(2,3,1)
box off
ylabel('PC 1')

subplot(2,3,4)
box off
ylabel('PC 2')
xlabel('Time [ms]')

subplot(2,3,2)
box off
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

subplot(2,3,3)
box off
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
title ('Network time constant [ms]')
legend(num2str(ts(1)),num2str(ts(2)),num2str(ts(3)))
    

subplot(2,3,6)
box off
xlabel('Network time constant [ms]')
ylabel('Estimated duration [ms]')



%% Calculate the distance between trajectories
% over all angles, compute Euclidean distance between trajectories of
% different angles
Data.maxdistance = zeros(numel(pars.angle));
for iAngle = 1:numel(pars.angle)
    for jAngle = iAngle+1:numel(pars.angle)
        
        % compute Euclidean distance
        for iPC = 1:pars.nPCs
            % for each PC. compute the difference between the
            % time-points at each angle
            diffsPCs(iPC,:) = squeeze(bsxfun(@minus,Data.temporalPCs_angle(iAngle,:,iPC),Data.temporalPCs_angle(jAngle,:,iPC)));
        end
        % compute Euclidean distance from those differences
        distancePCs(iAngle,jAngle,:) = sqrt(sum(diffsPCs.^2));
        
        % save max distance for analysis
        Data.maxdistance(iAngle,jAngle) = max(distancePCs(iAngle,jAngle,:));
        Data.maxdistance(jAngle,iAngle) = Data.maxdistance(iAngle,jAngle);
    end
end


% plot max distance between trajectories as a function of the angle between
% the inputs; one line per input
for iAngle = 1:numel(pars.angle)
    for jAngle = 1:numel(pars.angle)
        angle_difference_matrix(iAngle,jAngle) = pars.angle(iAngle) - pars.angle(jAngle);
    end
end

angle_diff=linspace(0,360,Ndir+1)-180;
angle_diff(end)=[];
angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); %center at zero
for iAngle = 1:numel(pars.angle)
    tmp(iAngle,:)=circshift(Data.maxdistance(iAngle,:),Ndir/2-iAngle);
    %plot(angle_diff,circshift(Data.maxdistance(iAngle,:),Ndir/2-iAngle),'.-'); hold on
end
subplot(2,3,5)
cla
errorbar(angle_diff,mean(tmp),std(tmp),'k','LineWidth',2)
box off
xlabel('Angle between input vectors (radians)')
ylabel('Mean distance between projections')
xlim([-180 180])
end