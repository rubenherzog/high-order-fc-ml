function [tc,dtc] = tc_dtc_from_covmat_gpu(covmat,N,all_min_1,bc1,bcN_min_1,bcN,ent_fun)
% function [tc,dtc] = tc_dtc_from_covmat_gpu(covmat,N,bc1,bcN_min_1,bcN,ent_fun)
% Computes the total correlation of gaussian data given their covariance
% matrix 'covmat'.
%
% INPUTS
% covmat = N x N covariance matrix
% all_min_1 = index of all minus one system
% bc1 = bias corrector for N=1
% bcN_min_1 = bias corrector for N-1 system
% bcN = bias corrector for N system
%
% OUTPUT
% tc = total correlation of the system with covariance matrix covmat.
% dtc
% entropy of multivaraite gaussian distribution, x is dimensionsionality
% and y is the variables variance of the covariance matrix determinant.

% ent_fun = @(x,y) 0.5.*log((2*pi*exp(1)).^(x).*y);
% all_min_1=(arrayfun(@(x) setdiff(1:N,x),1:N,'uni',0)');
detmv = det(covmat); % determinant
detmv_min_1=(cellfun(@(x) det(covmat(x,x)),all_min_1));
single_vars = diag(covmat); % variance of single variables

var_ents= ent_fun(1,single_vars) - bc1;
sys_ent = ent_fun(N,detmv) - bcN;
ent_min_one = ent_fun(N-1,detmv_min_1) - bcN_min_1;

tc = sum(var_ents) - sys_ent;
dtc = sum(ent_min_one) - (N-1).*sys_ent;
% o_info = tc - dtc;
% s_info = tc + dtc;