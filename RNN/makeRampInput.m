function [I_ramp,tRamp] = makeRampInput(timeStep,tRise,tFall,tsRise,tsFall,maxRamp)
% MAKERAMPINPUT make a scalar exponential rising ramp and fall
% [I,T] = MAKERAMPINPUT(timeSTEP,tRISE,tFALL,tsRISE,tsFALL,MAX)
%  computes an exponentially rising and falling ramp I, on time-base T,
%  defined by:
%   timeSTEP: duration of simulation time-step - defines units of all other
%   time parameters, so if in ms, then below in ms
%   tRISE: duration of rise (likely ms)
%   tFALL: duration of fall
%   tsRISE: time constant of rise
%   tsFALL: time constant of fall
%   MAX: maximum value of I, at peak of ramp (arbitrary units)
%
% 30/9/21 Initial version
%
% Mark Humphries 


t_ramp_rise = 0:timeStep:tRise; % times of rise component - centred on zero, so the curve starts near zero
t_ramp_fall = 0:timeStep:tFall; % times of fall component
A = maxRamp ./ exp(t_ramp_rise(end)/tsRise);  % exponential scaling parameter, to get desired maximum
I_ramp = A.*exp(t_ramp_rise/tsRise);    % rising part of input
I_ramp = [I_ramp maxRamp.*exp(-t_ramp_fall/tsFall)];    % concatenate falling part, starting at same maximum value

tRamp = [t_ramp_rise t_ramp_fall + t_ramp_rise(end) + timeStep];

end