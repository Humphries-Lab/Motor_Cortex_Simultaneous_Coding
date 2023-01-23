function [W,indexInhibit] = initialise_SOC_Weight_Matrix(nNeurons,pConnection,fractionE,ratioI,spectralRadius)

% initialise_SOC_Weight_Matrix create Dale's law weights ready for SOC optimisation
% [W,indexI] = initialise_SOC_Weight_Matrix(N,p,f,gamma,R) creates the NxN weight
% matrix W, defined by:
%   p: the connection density (0,1];
%   f: the fraction of excitatory neurons [0,1]
%   gamma: the ratio of inhibitory:excitatory weight strength
%   R: the desired spectral radius
% 
% Also returns a vector indexI of the indices of the columns containing inhibitory neurons
%
% 24/9/2021 Initial version
% Mark Humphries

% create connections
W = double(rand(nNeurons) < pConnection);

% select excitatory cells
indexExcite = randi(nNeurons,round(fractionE*nNeurons),1);
indexInhibit = setdiff(1:nNeurons,indexExcite);

% create weights
w_zero = spectralRadius / sqrt(pConnection * (1-pConnection) * (1+ratioI^2)/2);

% keyboard

W(:,indexExcite) = W(:,indexExcite) * w_zero / sqrt(nNeurons); % use adjaceny matrix as mask to create weights
W(:,indexInhibit) = W(:,indexInhibit) * -(ratioI * w_zero) / sqrt(nNeurons); 

