function [multi_order_hoi] = multi_order_hoi_meas_lin_v2(gc_covmat,sub_t,nplets_ids,par_op)
% nplets_ids = norders x 1 cell array with ids of nplets
% all orders have the same numbe of nplets
% nplets_ids is a cell array with all the nplets
ent_fun = @(x,y) 0.5.*log((2*pi*exp(1)).^(x).*y);
nnplets = size(nplets_ids,1);
% nnplets = length(nplets_ids);
nsubs = length(sub_t);
tc = zeros(nnplets,nsubs);
dtc = tc;

if par_op    
    parfor i=1:nnplets
        sel_ids = nplets_ids(i,:);
%         sel_ids = nplets_ids{i};
        thisn = size(sel_ids,2);
        all_min_1=(arrayfun(@(x) setdiff(1:thisn,x),1:thisn,'uni',0)');
        for s=1:nsubs
            % bias corrector
            bc_thisn = gaussian_ent_biascorr(thisn,sub_t(s));
            bc_nmin1 = gaussian_ent_biascorr(thisn-1,sub_t(s));
            bcn1 = gaussian_ent_biascorr(1,sub_t(s));
            
            % Computing measure for each n-plet
            sel_covmat = gc_covmat(sel_ids,sel_ids,s); % selection of covmat
            
            [tc(i,s),dtc(i,s)] = tc_dtc_from_covmat_gpu(sel_covmat,...
                thisn,all_min_1,bcn1,bc_nmin1,bc_thisn,ent_fun);
            
        end
    end
else
    
    for i=1:nnplets
        sel_ids = nplets_ids(i,:);
%         sel_ids = nplets_ids{i};
        thisn = size(sel_ids,2);
        all_min_1=(arrayfun(@(x) setdiff(1:thisn,x),1:thisn,'uni',0)');
        for s=1:nsubs
            % bias corrector
            bc_thisn = gaussian_ent_biascorr(thisn,sub_t(s));
            bc_nmin1 = gaussian_ent_biascorr(thisn-1,sub_t(s));
            bcn1 = gaussian_ent_biascorr(1,sub_t(s));
            
            % Computing measure for each n-plet
            sel_covmat = gc_covmat(sel_ids,sel_ids,s); % selection of covmat
            
            [tc(i,s),dtc(i,s)] = tc_dtc_from_covmat_gpu(sel_covmat,...
                thisn,all_min_1,bcn1,bc_nmin1,bc_thisn,ent_fun);
            
        end
    end
end
o_inf = tc-dtc;
s_inf = tc + dtc;

multi_order_hoi = cat(3,tc,dtc,o_inf,s_inf);