function RunDifferentTau(InputRange)
%% Function to test MCtx RNN SOC model response to variation in time constant of the RNN
% Current version: Changing the tau of the RNN temporally scale the
% trajectories preserving their shape.
%
% INPUT
% 
% InputRange: range of parameters of the vector input
%
% Note: to modify any other parameter, edit MCtx_RNN_SOC_common_parameters
%
% 31/01/2023
% Mark Humphries
% Andrea Colins Rodriguez

clearvars -except InputRange

run MCtx_RNN_SOC_common_parameters

% parameters of the vector input
pars.range = InputRange;

%% parameters unqiue to this set of simulations
pars.Tmax = pars.duration(1) + pars.rampTimeMax + max(pars.ts)/pars.dt*10; % set same max duration, to cope with longest ramp

ts=[0.8 1 1.2]*pars.ts; % select values of tau to evaluate

plot_traj_upto=2350;

tSteps = round(pars.Tmax ./ pars.dt);   % how many time-steps

Ndir=10;

% storage
nTau = numel(ts);
Data.temporalPCs = zeros(nTau,tSteps,pars.nPCs);

% create ramped input:
% generate ramped input - to maximum value
% from Hennequin et al 2014: exponential rise to max, then exponential
% decay

% time duration of input too
tInput = round([pars.duration(1) pars.duration(2) + pars.fall_duration]./pars.dt);

I_ramp = makeRampInput(pars.dt,pars.duration(end)-pars.duration(1),pars.fall_duration,pars.ts_rise,pars.ts_fall,pars.ramp_max);

fig2=figure;

colour_dir=My_hsv(Ndir);
colour_dur=plasma(nTau);
pars.angle = linspace(0,2*pi,Ndir+1); % want 11 of these... % different rotations of that vector: matching different reach angles in discrete analysis
pars.angle(end) = [];               % to get 10 reach angles, equally separated

%find subspace for all directions first
[coeffs,mean_data,initialProjection,threshold]=find_dir_subspace(pars,optimW,I_ramp,colour_dir);


for iAngle=1:Ndir
    rotatedProjection = rotate_n_dimensional_vector(initialProjection,pars.angle(iAngle)); % rotate current projection vector of ramp
    output_dur=zeros(nTau,1);
    
    for iTau = 1:nTau
        
        % storage
        a = zeros(Net.N,tSteps);
        r = zeros(Net.N,tSteps);
        
        
        
        % run simulation
        for iT = 2:tSteps
            % background input - set up here to add noise later
            I = zeros(Net.N,1) + pars.I_background;
            
            % add ramped input
            if iT > tInput(1) && iT <= tInput(2)
                I = I + I_ramp(iT - tInput(1)) .* rotatedProjection;
                
            end
            % update activity
            a(:,iT) = a(:,iT-1) + pars.dt * (-a(:,iT-1) + I + optimW * r(:,iT-1)) ./ ts(iTau);
            
            % update output
            r(:,iT) = neuron_output(a(:,iT),pars.output_type,pars.output_arg1,pars.output_arg2);
            
        end
        
        rActive = r(:,1:tSteps);        % indexing here to later allow for turning off early
        
        %% project onto common subspace
        
        rActive_center=rActive'-repmat(mean_data,tSteps,1);
        temporalPCs=rActive_center*coeffs;
        
        figure(fig2)
        
        % store top N PCs, for each Duration
        Data.temporalPCs(iTau,:,:) = temporalPCs(:,1:pars.nPCs);
        
        %% Estimate output duration using recurrence
        
        norms=vecnorm(temporalPCs(:,1:pars.nPCs)');
        idx_start=find(norms>=threshold,1,'First');
        idx_end=find(norms>=threshold,1,'Last');
        
        output_dur(iTau)=(idx_end-idx_start).*pars.dt;
        
        % plots examples
        if iAngle==1
            subplot(2,3,3)
            hold on
            plot3(Data.temporalPCs(iTau,tInput(1):plot_traj_upto,1),Data.temporalPCs(iTau,tInput(1):plot_traj_upto,2),Data.temporalPCs(iTau,tInput(1):plot_traj_upto,3),'Color',colour_dur(iTau,:),'LineWidth',2)
        end
        %% store data for angle only for one tau
        if iTau==2
            Data.temporalPCs_angle(iAngle,:,:) = temporalPCs(:,1:pars.nPCs);
        end
    end
    
    subplot(2,3,6),
    hold on
    plot(ts,output_dur,'Color',colour_dir(iAngle,:),'LineWidth',2)
    
    
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
view(-42,  44.3)

subplot(2,3,3)
box off
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
title ('Network time constant [ms]')
legend(num2str(ts(1)),num2str(ts(2)),num2str(ts(3)))
view(184, -12)

subplot(2,3,6)
box off
xlabel('Network time constant [ms]')
ylabel('Estimated duration [ms]')



%% Calculate the distance between trajectories
% over all angles, compute Euclidean distance between trajectories of
% different angles
Data.maxdistance = zeros(Ndir);
distancePCs=nan(Ndir,Ndir,tSteps);

for iAngle = 1:Ndir
    for jAngle = iAngle+1:Ndir
        diffsPCs=nan(pars.nPCs,tSteps);
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


angle_diff=pars.angle-pi;
angle_diff=angle_diff+(angle_diff(2)-angle_diff(1)); % center at zero
distance_shifted=zeros(Ndir,Ndir);
for iAngle = 1:numel(pars.angle)
    distance_shifted(iAngle,:)=circshift(Data.maxdistance(iAngle,:),Ndir/2-iAngle);
end

subplot(2,3,5)
errorbar(angle_diff*180/pi,mean(distance_shifted),std(distance_shifted),'k','LineWidth',2)
box off
xlabel('Angle between input vectors [^o]')
ylabel('Mean distance between projections')
xlim([-180 180])
end