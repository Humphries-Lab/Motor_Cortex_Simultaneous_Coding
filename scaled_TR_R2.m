function R2_scaled=scaled_TR_R2(X,Y,normalised_t)
%% scaled_TR_R2 scales trajectories X and Y to a normalised_t lenght and
%% computes their coefficient of determination
% 
% INPUTS
%
% X: n-dimensional trajectory. Rows are times, columns are dimensions 
%
% Y: n-dimensional trajectory. Rows are times, columns are dimensions 
%
% normalised_t: vector of timebins to normalise trajectories 
%
% OUTPUTS
%
% SI: Coefficient of determination between the scaled trajectories.
%
% 24/05/2023
% Andrea Colins Rodriguez

%normalise time
Xtime=linspace(0,1,size(X,1));
Ytime=linspace(0,1,size(Y,1));
% interpolate to match all vector lengths
Xnorm= interp1(Xtime,X,normalised_t);
Ynorm= interp1(Ytime,Y,normalised_t);

R2_scaled=R2(Xnorm,Ynorm);
end
function r=R2(X,Y)
%% R2 Computes the coefficient of determination between n-dimensional
%% trajectories X and Y
X2=reshape(X,numel(X),1);
Y2=reshape(Y,numel(Y),1);
ymean=mean(Y2);
Sres=sum((Y2-X2).^2);
Stot=sum((Y2-ymean).^2);
r=1-Sres/Stot;
end