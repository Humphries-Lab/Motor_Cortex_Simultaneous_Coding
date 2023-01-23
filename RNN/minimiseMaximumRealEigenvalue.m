function [newW,varargout] = minimiseMaximumRealEigenvalue(W,ICells,Options)

% minimiseMaximumRealEigenvalue optimises weight matrix to be stable
% newW = minimiseMaximumRealEigenvalue(W,IndexI,Options) iteratively alters
% the original NxN weight matrix W to minimise its maximum real eigenvalue
% and so enforce stability of RNN dynamics using W. 
%   IndexI : vector of indices of columns of W that are inhibitory
%   Optimisation defined by fields of Options struct:
%       .C : the ratio of the smoothed maximum eigenvalue to the actual
%       maximum (C > 1); set to 1.5 in Hennequin et al 2014
%       .B : the minimum value of the smoothed maximum eigenvalue; set to 
%       0.2 in Hennequin et al 2014
%       rate: the learning rate; set to 10 in Hennequin et al 2014
%       Threshold: threshold for convergence of maximum eigenvalue (change
%       on consecutive iterations of algorithm)
%       IndexI : index of all inhibitory neurons
%       fI : fraction of potential inhibitory connections that can be
%       modified
% 
% Notes:
%    (1) Source is Hennequin et al (2014) Neuron
%    (2) Initial version omits two constraint steps: normalisation of
%    weights; and limiting of size of potential inhibitory connections to
%    change
%
% 24/9/21 Initial version
% Mark Humphries

C = Options.C;
B = Options.B;
fractionI = Options.fractionI;
learningRate = Options.learningRate;
convergenceThreshold = Options.convergenceThreshold;

% choose set of inhibitory connections to modify
nNeurons = size(W,1);
AllowedI = logical(zeros(nNeurons)); 

% first add all connections that exist!
% AllowedI(:,ICells) = W(:,ICells) < 0; 

% Compute how many I cells we need to assign
possibleIconnections = nNeurons * numel(ICells);        % total number of I connections    
nIConnectionsToChoose = possibleIconnections * fractionI;   % how many of these to assign to weight changes
ixIconnections = find(W(:,ICells) < 0);         % indices of currently existing I connections, to will be changed
nRemainingIConnectionsToChoose = nIConnectionsToChoose - numel(ixIconnections);  % how many connections are left to assign to weight change

% find the indices of those that could be assigned
ixUnconnectedIconnections = find(W(:,ICells) == 0);        % indices of all possible I connectionsixNotConnected = setdiff(ixPossibleConnections,ixIconnections);  % indices of remaining unassigned I connections to choose from

% randomly choose from that potential group
ixChosenIconnections = ixUnconnectedIconnections(randi(numel(ixUnconnectedIconnections),nRemainingIConnectionsToChoose,1));

%keyboard
% write temporary mask the same size as Nx(number of I cells), to write to
% NxN mask
tmp = zeros(nNeurons,numel(ICells)); 
tmp(ixChosenIconnections) = 1;  % write in chosen unconnected
tmp(ixIconnections) = 1;        % write in already connected
AllowedI(:,ICells) = tmp;       % assign to NxN mask

% sanity check
% ixAllowedI = find(AllowedI);
% ixExistingI = find(W < 0);
% c = setdiff(ixExistingI,ixAllowedI); % should be empty

% run algorithm
blnConverged = 0;
iteration = 1;
newW = W;

while ~blnConverged
    % compute current smoothed spectral abcissa = estimate of max real eigenvalue
    % max_eig(iteration) = real(eigs(newW,1,'largestreal'));  % make sure we use only the real part...
    eigs = eig(newW,'vector');
    max_eig(iteration) = max(real(eigs));
    alpha_smooth(iteration) = max(C*max_eig(iteration),max_eig(iteration)+B);  % current estimate
    
    % plot(alpha_smooth); drawnow
    % check convergence
    if iteration > 1
        if abs(alpha_smooth(iteration) - alpha_smooth(iteration-1)) < convergenceThreshold
            blnConverged = 1;
        end
    end
    
    if ~blnConverged
        % do gradient descent
        Identity = eye(nNeurons);
        shiftedWeights = newW - alpha_smooth(iteration) .* Identity;
        % lyap(A,Q) solves: A*X + X*A' + Q = 0
        % but equation for Q is in form A'Q + QA = -2I; so we convert
        % to AQ + QA' = -2I by transposing A
        Q = lyap(shiftedWeights',2*Identity);  
        P = lyap(shiftedWeights,2*Identity);

        % find gradient
        QP = Q*P;
        gradient = QP / trace(QP);       
        
        % update weights
        newWeights = newW(AllowedI) - learningRate .* gradient(AllowedI);

        % clip any positive ones
        newWeights(newWeights > 0) = 0;

        % update matrix
        newW(AllowedI)  = newWeights; % does this work?
        
        iteration = iteration + 1;
    end
    
end
varargout{1} = max_eig;


