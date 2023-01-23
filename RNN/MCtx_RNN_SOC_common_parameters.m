% script to set common parameters for motor cortex RNN model

% where to export figure panels
% if ispc
%     exportpath = 'C:\Users\lpzmdh\Dropbox\Projects\Transitions in neural activity\RNN_model_of_MCtx_dynamics\SOC model networks';
% else
%     exportpath = '/Users/mqbssmhg/Dropbox/Projects/Transitions in neural activity/RNN_model_of_MCtx_dynamics/SOC model networks';
% end
rng('default')
%% load SOC wiring
if ispc
    network_path = '';
else
    network_path = '';
end

networkID = 'SOC_Network_20211006T143906';%'SOC_Network_20210927T141434';
load([network_path networkID]);


%% initialise inputs

% stepped inputs driving movement
pars.n_input_batches = 1; % 10;         % repeats of sampling initial input vector

% parameters of the ramp input - pars.duration sets length of ramp rise
pars.duration = [1000,2000];               % time-window of input (ms)
pars.fall_duration = 200; % duration of ramp fall - milliseconds
pars.ts_rise =200;         % time constants of ramp
pars.ts_fall = 0.01;
pars.ramp_max = 10;         % maximum value of the ramp
pars.input_projection = rand(Net.N,1); % projection of input to RNN units - overwritten when rotating
pars.N=200;
% parameters of the vector input
pars.range = [0 1];

% spontaneous input?
pars.I_background = 0.1;

%% initialise simulation parameters
% neuron
pars.ts = 200;     % slow dynamics, 200 ms in Hennequin et al 2014

% solution
pars.dt = 1;     % time-step, milliseconds
pars.Tmax = pars.duration(end) + pars.ts/pars.dt*3;      % length of simulation, milliseconds

% initialise output function
pars.output_type = 'linear';
pars.output_arg1 = 30;            % baseline output
pars.output_arg2 = 100;           % maximum output
% storage
pars.nPCs = 5;      % how many PCs of output to compute distances on

%% common analysis parameters
pars.epsilon = 1e-6;  % threshold for "no change"

