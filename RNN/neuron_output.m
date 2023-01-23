function output = neuron_output(a,output_type,varargin)

% NEURON_OUTPUT nonlinear output functions for rate-coding models
% O = NEURON_OUTPUT(A,TYPE), given a vector A of activation values,
% transforms them using the nonlinear output function of TYPE:
%       'ReLu' : rectified linear output - output in [0,inf)
%                   NEURON_OUTPUT(...,SLOPE): optional slope parameter SLOPE
%       'tanh' : classic tanh output - output in [-1,1]
%       'tanh_baseline' : tanh output as deviation from baseline:
%                  NEURON_OUTPUT(...,R0,RMAX): baseline R0, max rate RMAX
%                  NEURON_OUTPUT(...,R0,RMAX,GAIN): optional gain parameter
%       'sigmoid' : classic sigmoid with slope:
%                  NEURON_OUTPUT(...,SLOPE): slope parameter SLOPE
%       'power' : rectified power-law scaling output (supralinear) - output in [0,inf] 
%                  NEURON_OUTPUT(...,K,N): scaled by K, raised by power N
%       'linear' : no output function at all!
%
% Note that any optional parameter can be a scalar (applied to all A) or a
% vector of length(A), each value applied to that unit's activation in A
%
% Returns:
% O : the vector of resulting outputs
%
% Notes:
%   ReLu: slope parameter is m*a (default is m=1)
%   sigmoid: from Murray & Escola 2017: 1/(1+exp-SLOPE*x)
%   tanh_baseline: see that function; default gain = 1
%   Power-law output from Rubin et al (2015) Neuron
%
% 27/9/2021: moved tanh_baseline function into this one; and added option
% for gain in tanh_baseline
% 17/9/2021: added option for linear output
% 14/9/2021: added sigmoid, and slope to ReLu
% 4/7/2021 : initial version
% Mark Humphries

nNeurons = numel(a);

switch output_type
    case 'linear'
        output = a;
    case 'ReLu'
        m = ones(nNeurons,1);
        if nargin > 2
            m = m .* varargin{1};
        end
        output = m .* a .* (a > 0);
    case 'tanh'
        output = tanh(a);
    case 'tanh_baseline'
        gain = ones(nNeurons,1);
        if nargin > 4
            gain = gain.* varargin{3};
        end
        output = tanh_baseline(a,varargin{1},varargin{2},gain);
    case 'sigmoid'
        output = 1./ (1 + exp(-varargin{1} .* a));
    case 'power'
        output = varargin{1} * (a .* (a > 0)) .^varargin{2};
    otherwise
        error('Unknown neuron output function');
end

function r = tanh_baseline(x,r_zero,r_max,gain)

% TANH_BASELINE sigmoidal neuron output function with baseline
% R = TANH_BASELINE(X,R0,RMAX,G) given vector of inputs X and desired
% baseline (R0) and maximum (RMAX) firing rates (in spikes per second), and
% output gain vector G:
% returns the output R, after transforming through a tanh function with added baseline
%
% Reference:
% Hennequin et al (2014) Neuron, 82, 1394-1406
% Stroud et al (2018) Nature Neuroscience 21, 1774 (added Gain parameter)
%
% Mark H, 22/7/21

r = zeros(size(x));

r(x < 0) = r_zero * tanh(gain(x<0).*x(x<0) / r_zero);
r(x >= 0) =  (r_max - r_zero) * tanh(gain(x>=0).*x(x>=0) / (r_max - r_zero));