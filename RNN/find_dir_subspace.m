function [coeffs,mean_data,initialProjection,threshold]=find_dir_subspace(pars,optimW,I_ramp,colour_dir)
%% find_dir_subspace defines a common subspace for trajectories across different directions
% INPUTS
%
% pars: parameters of the simulation set on MCtx_RNN_SOC_common_parameters
% 
%optimW: Weights of the network
%
%I_ramp: stimulus to the network
%
%colour_dir: colour map describing the colour of each direction
%
%OUTPUTS
%
% coeffs: Vectors defining the common subspace (PCs)
%
% mean_data: Mean FR across directions
% 
% initialProjection: Initial projection of the input into the network
%
% threshold: distance equivalent to the theta percentil of all distances in
% the data
%
% 31/07/2023
% Andrea Colins
% Mark Humphries

Nneurons=size(optimW,1);
Ndir=size(colour_dir,1);

%% run network
tSteps = round(pars.Tmax ./ pars.dt);   % how many time-steps

% % time duration of input too
 tInput = [round(pars.duration(1) ./pars.dt), round((pars.duration(2) + pars.fall_duration) ./ pars.dt)];

initialProjection = pars.range(1) + (pars.range(2) - pars.range(1)) * rand(Nneurons,1);

rStore = zeros(numel(pars.angle) * tSteps,Nneurons);

for iAngle = 1:Ndir
    % transform input vector
    rotatedProjection = rotate_n_dimensional_vector(initialProjection,pars.angle(iAngle)); % rotate current projection vector of ramp
    
    % run simulation
    a = zeros(Nneurons,tSteps);
    r = zeros(Nneurons,tSteps);
    
    for iT = 2:tSteps
        
        % background input
        I = zeros(Nneurons,1) + pars.I_background;
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
% compute PCA on all concatenated outputs
[coeffs,temporalPCs] = pca(rStore);   % do PCA on all neuron outputs
Npoints_per_dir=(tSteps-tInput(1));
temporalPCs2=nan(Ndir*Npoints_per_dir,pars.nPCs);
% store top N PCs, for each angle

for iAngle = 1:numel(pars.angle)
    
    temporalPCs2(1+(iAngle-1)*Npoints_per_dir:iAngle*Npoints_per_dir,:)=temporalPCs(tInput(1)+1+(iAngle-1)*tSteps:tSteps+(iAngle-1)*tSteps,1:pars.nPCs);
    
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
threshold=prctile(dist_all_points,pars.threshold_rec);
end