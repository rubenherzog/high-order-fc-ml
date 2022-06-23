function [best_d,best_id_perk] = hoi_greedy_cohen_d_v2(maxk,gc_covmat,subt,meas,minmax,sub_id)
%%
% sel_subs = strcmp(fmri.sub_cond,condnames{c+1}) | strcmp(fmri.sub_cond,condnames{1});
% c=1;
% sel_subs = strcmp(eeg.sub_cond,condnames{c+1}) | strcmp(eeg.sub_cond,condnames{1});
% gc_covmat = covmats{2}(:,:,sel_subs);
% meas = 3; %
% minmax = 'min';
% maxk = 10;
% subt = sub_t{2}(sel_subs);
% sub_id = 0 for control, 1 for condition

%
% ent_fun = @(x,y) 0.5.*log((2*pi*exp(1)).^(x).*y);
% nsubs = size(subt,1);
N = size(gc_covmat,1);
int_k = 2:maxk;
norders = length(int_k);

% Creating pairwise interactions
pw_ints=nchoosek(1:N,2);
npw = length(pw_ints);

% Computing measures for pairwise interactions
[multi_order_hoi] = multi_order_hoi_meas_lin_v2(gc_covmat,subt,pw_ints,1);
tc = multi_order_hoi(:,:,1); % always starts with tc
% Computing Effect Sizes
x0 = tc(:,sub_id==0); % CN
x1 = tc(:,sub_id==1); % Condition
[this_d,~] = computeCohen_d_rh(x1, x0);



best_id_perk = cell(npw,norders);
best_id_perk(:,1) = num2cell(pw_ints,2);

best_d = zeros(npw,norders);
best_d(:,1) = this_d;

% Computing max or minima
k2 = norders;
switch minmax
    case 'min'
%         tic
%         parforange = 1:floor(npw/10);
        parforange = 1:npw;
        
        parfor p=parforange  % looping over pairs
            % Looping over orders
            
            this_best_regs = pw_ints(p,:);
            for k=2:k2
                % new ids
                new_ids = cat(2,repmat(this_best_regs,[N,1]),(1:N)');
                new_ids(this_best_regs,:) =[]; % removing repeated ids
                % Computing hoi measures
                [multi_order_hoi] = multi_order_hoi_meas_lin_v2(gc_covmat,subt,new_ids,0);
                this_meas = multi_order_hoi(:,:,meas);
                % Computing cohen's d
                x0 = this_meas(:,sub_id==0); % CN
                x1 = this_meas(:,sub_id==1); % Condition
                [this_d,~] = computeCohen_d_rh(x1, x0);                
                % Compute min
                [this_best_d,this_best_id] = min(this_d);
                this_best_regs = new_ids(this_best_id,:);
                % storing
                best_id_perk{p,k} = this_best_regs;
                best_d(p,k) = this_best_d;
                
                
            end
            
        end
%         toc
        
    case 'max'
%         tic
%         parforange = 1:floor(npw/10);
        parforange = 1:npw;
        
        parfor p=parforange  % looping over pairs
            % Looping over orders
            
            this_best_regs = pw_ints(p,:);
            for k=2:k2
                % new ids
                new_ids = cat(2,repmat(this_best_regs,[N,1]),(1:N)');
                new_ids(this_best_regs,:) =[]; % removing repeated ids
                % Computing hoi measures
                [multi_order_hoi] = multi_order_hoi_meas_lin_v2(gc_covmat,subt,new_ids,0);
                this_meas = multi_order_hoi(:,:,meas);
                % Computing cohen's d
                x0 = this_meas(:,sub_id==0); % CN
                x1 = this_meas(:,sub_id==1); % Condition
                [this_d,~] = computeCohen_d_rh(x1, x0);                
                % Compute max
                [this_best_d,this_best_id] = max(this_d);
                this_best_regs = new_ids(this_best_id,:);
                % storing
                best_id_perk{p,k} = this_best_regs;
                best_d(p,k) = this_best_d;
                
                
            end
            
        end
%         toc
    otherwise
        error('Wrong minmax parameter. Should be min or max')
end
% 
% figure,plot(int_k,best_d,'o-');hold on
% plot(int_k,max(best_d),'r*-');hold on
% plot(int_k,min(best_d),'b*-');hold on