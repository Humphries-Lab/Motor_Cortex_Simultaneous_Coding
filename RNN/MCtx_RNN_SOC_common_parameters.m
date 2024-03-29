% script to set common parameters for motor cortex RNN model 

rng('default')

networkID = 'SOC_Network_20211006T143906';
load(networkID);


%% initialise inputs

% stepped inputs driving movement

% parameters of the ramp input - pars.duration sets length of ramp rise
pars.duration = [1000,2000];               % time-window of input (ms)
pars.fall_duration = 200; % duration of ramp fall - milliseconds
pars.ts_rise =200;         % time constants of ramp
pars.ts_fall = 0.1;
pars.ramp_max = 10;         % maximum value of the ramp
pars.input_projection = rand(Net.N,1); % projection of input to RNN units - overwritten when rotating
pars.N=200;
pars.rampTimeMax = 1000; % time of maximum value


% spontaneous input?
pars.I_background = 0.1;

%% initialise simulation parameters
% neuron
pars.ts = 200;     % slow dynamics, 200 ms in Hennequin et al 2014

% solution
pars.dt = 1;     % time-step, milliseconds
pars.Tmax = pars.duration(end) + pars.ts/pars.dt*3;      % length of simulation, milliseconds

% initialise output function
pars.output_type = 'tanh_baseline';
pars.output_arg1 = 5;            % baseline output
pars.output_arg2 = 70;           % maximum output
% storage
pars.nPCs = 5;      % how many PCs of output to compute distances on

%% common analysis parameters
pars.epsilon = 1e-6;  % threshold for "no change"

% threshold (percentil) to define recurrence of trajectories
pars.threshold_rec=20;