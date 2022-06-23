function [d,s] = computeCohen_d_rh(x1, x2)
% 
% call: d = computeCohen_d(x1, x2, varargin)
% 
% EFFECT SIZE of the difference between the two 
% means of two samples, x1 and x2 (that are vectors), 
% computed as "Cohen's d". 
% 
% If x1 and x2 can be either two independent or paired 
% samples, and should be treated accordingly:
%  
%   d = computeCohen_d(x1, x2, 'independent');  [default]
%   d = computeCohen_d(x1, x2, 'paired');
% 
% Note: according to Cohen and Sawilowsky:
%
%      d = 0.01  --> very small effect size
%      d = 0.20  --> small effect size
%      d = 0.50  --> medium effect size
%      d = 0.80  --> large effect size
%      d = 1.20  --> very large effect size
%      d = 2.00  --> huge effect size
%
%
% Ruggero G. Bettinardi (RGB)
% Cellular & System Neurobiology, CRG
% -------------------------------------------------------------------------------------------
%
% Code History:
%
% 25 Jan 2017, RGB: Function is created
  

% basic quantities:
n1       = size(x1,2);
n2       = size(x2,2);
mean_x1  = nanmean(x1,2);
mean_x2  = nanmean(x2,2);
var_x1   = nanvar(x1,0,2);
var_x2   = nanvar(x2,0,2);
meanDiff = (mean_x1 - mean_x2);



sv1      = ((n1-1)*var_x1);
sv2      = ((n2-1)*var_x2);
numer    =  sv1 + sv2;
denom    = (n1 + n2 - 2);
pooledSD =  sqrt(numer / denom); % pooled Standard Deviation
s        = pooledSD;             % re-name
d        =  meanDiff ./ s;        % Cohen's d (for independent samples)
