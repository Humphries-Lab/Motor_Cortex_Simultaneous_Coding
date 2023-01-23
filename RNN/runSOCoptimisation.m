%% script to run SOC optimisation and save results
% Mark Humphries 
clearvars; close all

%% initialise
% parameterise network
Net.N = 200;
Net.p = 0.1;
Net.fractionExcite = 0.5;
Net.ratioI = 3;
Net.R = 10;

% parameterise algorithm
Optim.C = 1.5;
Optim.B = 0.2;
Optim.learningRate = 10;
Optim.convergenceThreshold = 0.01;
Optim.fractionI = 0.5;

% create network
[W,Icells] = initialise_SOC_Weight_Matrix(Net.N,Net.p,Net.fractionExcite,Net.ratioI,Net.R);

% get eigenvalues of that network
initial_eigs = eig(W,'vector');
f = figure;
plot(real(initial_eigs),imag(initial_eigs),'.','Color',[0.7 0.7 0.7]); hold on

%% optimise W

[optimW,max_real_eig] = minimiseMaximumRealEigenvalue(W,Icells,Optim);


% save resulting network
fileID = datestr(now,30);

save(['.\SOC_Network_' fileID],'optimW','Net','Optim');

%% visualise results
% plot eigenvalues
SOCeigs = eig(optimW,'vector'); 
plot(real(SOCeigs),imag(SOCeigs),'.','Color',[0.7 0.2 0.4]); hold on

% add lines through origin to plot
line([0 0],f.CurrentAxes.YLim,'Color','k');
line(f.CurrentAxes.XLim,[0 0],'Color','k');
xlabel('Real(\lambda)'); ylabel('Im(\lambda)')

% plot minimisation curve
figure
plot(max_real_eig)
xlabel('Iteration'); ylabel('Maximum real eigenvalue of W')
