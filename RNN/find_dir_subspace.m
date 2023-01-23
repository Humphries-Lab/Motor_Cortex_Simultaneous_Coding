function [coeffs,mean_data,std_data,initialProjection,threshold]=find_dir_subspace(pars,optimW,Net)
Ndir=10;
colour_dir=hsv(Ndir);
pars.angle = linspace(0,2*pi,Ndir+1); % want 11 of these... % different rotations of that vector: matching different reach angles in discrete analysis
pars.angle(end) = [];               % to get 10 reach angles, equally separated

%% run network
tSteps = round(pars.Tmax ./ pars.dt);   % how many time-steps

% time duration of input too
tInput = [round(pars.duration(1) ./pars.dt), round((pars.duration(2) + pars.fall_duration) ./ pars.dt)];

% create ramped input
[I_ramp] = makeRampInput(pars.dt,pars.duration(end)-pars.duration(1),pars.fall_duration,pars.ts_rise,pars.ts_fall,pars.ramp_max);
% storage
nAngles = numel(pars.angle);
Data.temporalPCs = zeros(pars.n_input_batches,nAngles,tSteps,pars.nPCs);
initialProjection = pars.range(1) + (pars.range(2) - pars.range(1)) * rand(Net.N,1);

rStore = zeros(numel(pars.angle) * tSteps,Net.N);
for iAngle = 1:numel(pars.angle)
    % transform input vector
    rotatedProjection = rotate_n_dimensional_vector(initialProjection,pars.angle(iAngle)); % rotate current projection vector of ramp
    
    % run simulation
    a = zeros(Net.N,tSteps);
    r = zeros(Net.N,tSteps);
    
    for iT = 2:tSteps
        
        % background input
        I = zeros(Net.N,1) + pars.I_background;
        % add ramped input
        if iT > tInput(1) && iT <= tInput(2)
            I = I + I_ramp(iT - tInput(1)) .* rotatedProjection;
        end
        
        % update activity
        a(:,iT) = a(:,iT-1) + pars.dt * (-a(:,iT-1) + I + optimW*r(:,iT-1)) ./ pars.ts;
        
        % update output
        r(:,iT) = neuron_output(a(:,iT),pars.output_type,pars.output_arg1,pars.output_arg2);
        
    end
    
    % store outputs, concatenated
    rStore(1+(iAngle-1)*tSteps:iAngle*tSteps,:) = r(:,1:tSteps)';  % columns are neurons
    
end

mean_data=mean(rStore);
std_data=std(rStore);
% compute PCA on all concatenated outputs
[coeffs,temporalPCs] = pca(rStore);   % do PCA on all neuron outputs
temporalPCs2=[];
% store top N PCs, for each angle

for iAngle = 1:numel(pars.angle)
    
    temporalPCs2=[temporalPCs2;temporalPCs(1001+(iAngle-1)*tSteps:4000+(iAngle-1)*tSteps,1:pars.nPCs)];
    subplot(2,3,1)
    plot(temporalPCs(1+(iAngle-1)*tSteps:iAngle*tSteps,1),'Color',colour_dir(iAngle,:))
    hold on
    subplot(2,3,4)
    plot(temporalPCs(1+(iAngle-1)*tSteps:iAngle*tSteps,2),'Color',colour_dir(iAngle,:))
    hold on
    
    subplot(2,3,2)
    plot3(temporalPCs(1+(iAngle-1)*tSteps:iAngle*tSteps,1),temporalPCs(1+(iAngle-1)*tSteps:iAngle*tSteps,2),temporalPCs(1+(iAngle-1)*tSteps:iAngle*tSteps,3),'Color',colour_dir(iAngle,:))
    hold on

end
dist_all_points=pdist(temporalPCs2(:,1:pars.nPCs));
threshold=prctile(dist_all_points,20);
end