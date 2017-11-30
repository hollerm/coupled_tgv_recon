function [u,tgv,datafit,gap,par_out,operator_pars] = tgv_recon_joint(filename,par_in)

%%Set Parameter##################################################################################################################
%All parameters defined here can be set via the test script by adding the parameter and the value to the "par_in" cell array

%Steplength
sig = sqrt(1/12);%12
tau = sqrt(1/12);%12
s_t_ratio = 1;

%Steplength
dx = 1;
dy = 1;
dz = 1;

%Stopping rule   
stop_rule = 'iteration';
stop_par = 200;

%Evaluate optimality (set 0 for speed)
eval_opt = 0;

%TGV weights
alph0 = sqrt(2);
alph1 = 1;

%Regularization parameter
ld = 1;

%Result file for intermediate saving
save_inter = 0;
result_file = 'testresult';

%Operator rescaling ( [gaussian_rescale,poisson_rescale] )
op_rescale = [3,10]; %If scalar one, no rescaling is performed.

%Voctor norm coupling for gradient (and symmetrized gradient)
%'nuc'... nuclear norm,
%'frob' .... frobenius norm,
%'sep' ... separate channels
vec_norm = 'frob';

%Name of global parfun: If set, uses this parfun instead of the one stored in the file
global_parfun = 0; %If = 0, parfun is taken from file

%Data subset
data_subset = 0; %Select only a subset of components, e.g. data_subset = [1 3 5] selects only components 1,3,5



%%Read parameter#################################################################################################################

%Input: par_in
%Generate list of parameters
vars = whos;
for l=1:size(vars,1)
    par_list{l,1} = vars(l).name;
end

%Set parameters according to list
for l=1:size(par_in,1);
    valid = false;
    for j=1:size(par_list,1); if strcmp(par_in{l,1},par_list{j,1})
            valid = true;
            eval([par_in{l,1},'=','par_in{l,2}',';']);
    end; end
    if valid == false; warning(['Unexpected parameter at ',num2str(l)]); end
end

%Stepsize ratio update
if abs(sig/tau - s_t_ratio^2) > 10^(-8)
    if abs(sig/tau - 1) < 10^(-8)
        sig = sig*s_t_ratio;
        tau = tau/s_t_ratio;
    else
        warning('Incompatible stepsize initialization and choice:')
        display(['sig/tau: ',num2str(sig/tau)]);
        display(['s_t_ratio^2: ',num2str(s_t_ratio^2)]);
    end
end


%%Load Data######################################################################################################################

load([filename,'.mat']);

%Set data subset
if data_subset == 0
    noc = length(data);
    data_subset = [1:noc];
else
    noc = length(data_subset);
end 

%Set operator rescaling
mr_scale = op_rescale(1);
pt_scale = op_rescale(2);
op_rescale = ones(1,noc);


for idx=1:noc
    
    %Set parfun
    if global_parfun == 0
        parfun = str2func(data{data_subset(idx)}.parfun);
    else
        parfun = str2func(global_parfun);
    end
    
    %Sensitivities in separate file
    if isfield(data{data_subset(idx)},'sens_source')
        
        if ~exist('imgSens')
            load(data{data_subset(idx)}.sens_source);
        end
        
        %Expecting operator_pars.ntype and operator_pars.norm_fac the operator to have norm \approx 1
        fout = parfun(data{data_subset(idx)},imgSens);
            
    %Sensitivities in data file
    else
        %Expecting operator_pars.ntype and operator_pars.norm_fac the operator to have norm \approx 1
        fout = parfun(data{data_subset(idx)});
        
    end
    
    %Set required function output    
    K{idx} = fout.K;
    Kt{idx} = fout.Kt;
    f0{idx} = fout.f0;
    %Set optional Optional output
    if ~isfield(fout,'c0')
        fout.c0 = 0;
    end
    c0{idx} = fout.c0;
    
    %Optional output
    if ~isfield(fout,'operator_pars')
        fout.operator_pars = struct;
    end
    if ~isfield(fout.operator_pars,'op_norm_fac')
        fout.operator_pars.op_norm_fac = 1.0;
    end
    if ~isfield(fout.operator_pars,'ntype')
        fout.operator_pars.ntype = 'gaussian';
    end
    operator_pars{idx} = fout.operator_pars;

        
    %Rescale operators
    if strcmp(operator_pars{idx}.ntype,'gaussian')
        op_rescale(idx) = mr_scale;
    elseif strcmp(operator_pars{idx}.ntype,'poisson')
        op_rescale(idx) = pt_scale;
    end
    K{idx} = @(x) op_rescale(idx)*K{idx}(x);
    Kt{idx} = @(x) op_rescale(idx)*Kt{idx}(x);
        
    
    %Rescale data such that the max is \approx 100
    mtol = 0.8;
    u0_tmp = abs( Kt{idx}(f0{idx} - c0{idx} ) );
    u0_tmp = u0_tmp( u0_tmp>mtol*max(u0_tmp(:)) );
    mx = sum(u0_tmp)/max(size(u0_tmp));
    data_rescale(idx) = (100/mx);

    f0{idx} = f0{idx}*data_rescale(idx);
    c0{idx} = c0{idx}*data_rescale(idx);
 
    
end    
    
%%Multiply by hard-coded factor to equalize cost for all norm-types. Factors were obtained by comparing the norms on randomly generated data
norm_fac_frob = [1.0000    0.8248    0.7534    0.7263    0.7126    0.7041    0.6986    0.6945    0.6915    0.6894    0.6876   0.6862    0.6846    0.6836    0.6828    0.6819    0.6813    0.6807    0.6802    0.6796    0.6794    0.6788   0.6785    0.6782    0.6778    0.6776    0.6773    0.6771    0.6769    0.6767    0.6765    0.6763    0.6761   0.6759    0.6759    0.6757    0.6756    0.6754    0.6754    0.6753 ];

norm_fac_sep = [1.0000    1.1413    1.2685    1.4076    1.5412    1.6664    1.7843    1.8949    2.0003    2.1012    2.1973 2.2899    2.3772    2.4628    2.5459    2.6257    2.7038    2.7792    2.8530    2.9246    2.9952    3.0632 3.1304    3.1960    3.2601    3.3233    3.3845    3.4457    3.5056    3.5643    3.6221    3.6786    3.7350 3.7900    3.8444    3.8979    3.9509    4.0031    4.0551    4.1059 ];



%Update regularization parameters
if length(ld) ~= noc
    if length(ld) == 1
        display('Using automatic parameter reweighting');
        ld = ld*ones(1,noc);
        %Rescaling with 10^6 to get human readable range
        ld = 10^(6)*ld./data_rescale;
    else
        warning('Inconsistend regularization parameter length and number of components, setting all regularization parameters to 1');
        ld = ones(1,noc);
    end
end

if noc <=40
    if strcmp(vec_norm,'frob') 
        norm_fac = norm_fac_frob(noc);
        ld = ld*norm_fac;
        
    elseif strcmp(vec_norm,'sep')
        norm_fac = norm_fac_sep(noc);
        ld = ld*norm_fac;

    end
else
    warning('Cost equalization for norms only implemented for up to 40 components. No cost equalization performed.');
end



%TGV-reconstruction

%%Initialization#################################################################################################################
%display('Initializing data...')

% Get data sizes from first channel
%noc = number_of_components
tmp = Kt{1}(f0{1} - c0{1});
[n,m,k] = size(tmp);



% Primal iterate, extragradient and dual iterate
x = zeros(n,m,k,4,noc); 
ext = zeros(n,m,k,4,noc); 

z = zeros(n,m,k,9,noc); 

r = cell(1,noc);
for ii=1:noc
    r{ii} = zeros(size(f0{ii}));
    %Sanity check
    if strcmp(operator_pars{ii}.ntype,'poisson') && min(f0{ii}(:)) < 0
        error('Poisson data must not be negativ');
    end
    if strcmp(operator_pars{ii}.ntype,'poisson')  && min(f0{ii}(:)) <= 0
        display('setting positive data')
        f0{ii}(f0{ii}<0.00001) = 0.00001;
    end


end

%Give warning when using the nuclear norm with single component
if strcmp(vec_norm,'nuc') && noc == 1
    warning('Incompatible parameter choice: Nuclear norm coupling with single-channel image. Using Frobenius-norm coupling.');
    vec_norm = 'frob';
end


% Initialization 
x(:,:,:,1,1) = tmp;
for ii=2:noc
    x(:,:,:,1,ii) = Kt{ii}(f0{ii} - c0{ii});
end


% Initialize extragradient
ext = x;

% Initialize tgv and datafit
tgv=zeros(stop_par+1,1);
datafit = zeros(stop_par+1,noc);
gap = zeros(stop_par+1,1); %Relative primal-dual gap


if eval_opt
    tgv(1) = get_tgv_joint(x,alph1,alph0,vec_norm,dx,dy,dz);
    datafit(1,:) = get_datafit_rec(x,ld,f0,c0,K,operator_pars);
    gap(1) = ( tgv(1)+ sum(datafit(1,:)) + gap_data_dual(r,ld,f0,c0,operator_pars) + gap_reg_dual(z,r,Kt,operator_pars,dx,dy,dz) )/(n*m*k*noc);
end

%%Main loop######################################################################################################################
kk=0;
cont =true;
display('Starting iterations...')
fprintf('Completed iteration ');
fprintf('%-10d',0)
while cont
    
    
    % Algorithmic issues
    % ------------------
    % Dual updates
    % ------------------
	%display('Dual step...')
    % Update TGV dual variables
    z = z + sig*cat(4, fgrad_3(ext(:,:,:,1,:)) - ext(:,:,:,2:4,:) , sym_bgrad_3(ext(:,:,:,2:4,:)) );
    
    % TGV Dual Prox
    z = prox_dual(z,alph1,alph0,vec_norm);


    
    
    %Update dual data variables
    for ii=1:noc
        r{ii} = r{ii} + sig*K{ii}(ext(:,:,:,1,ii));
        r{ii} = data_prox_dual(r{ii} + sig*c0{ii},f0{ii},ld(ii),sig,operator_pars{ii}.ntype);
    end


    % ------------------
    % Primal updates
    % ------------------
	%display('Primal step...')

	 ext=x-tau*cat(4,-bdiv_3(z(:,:,:,1:3,:),dx,dy,dz),-z(:,:,:,1:3,:)-fdiv_3(z(:,:,:,4:9,:),dx,dy,dz));



	 
	 for ii=1:noc
	 
	    ext(:,:,:,1,ii) = ext(:,:,:,1,ii) - tau*Kt{ii}(r{ii});
     
        % Additional prox für poisson data
        if strcmp(operator_pars{ii}.ntype,'poisson')
            
            %Bound to real data (only necessary for nuc)
            if strcmp(vec_norm,'nuc')
                ext(:,:,:,1,ii) = real(ext(:,:,:,1,ii));
            end

            %Bound to positive values
            tmp = ext(:,:,:,1,ii);
            tmp(tmp<0) = 0;
            ext(:,:,:,1,ii) = tmp;
        end
    end
    

    % Extragradient
	%display('Extragradient...')
    x=2*ext - x;


    [x,ext] = deal(ext,x); %Swap extragradient and primal variable

    % Adapt stepsize: Automatic stepsize adaption. For computational reasons stepsize is not always checked
    %   ->Adapt if convergence problems occure
    if (kk<50) || (rem(kk,50) == 0)
        [sig,tau] = steps_tgv_recon_joint(ext-x,sig,tau,s_t_ratio,K,dx,dy,dz);
        %display(['Sig: ',num2str(sig)])
    end
    
    kk = kk+1; %Update counter

    % Save intermediate result
    if save_inter ~= 0
        if (rem(kk,save_inter) == 0)
        
            for idx = 1:noc
                u{idx} = squeeze(x(:,:,:,1,idx))*op_rescale(idx)*operator_pars{idx}.op_norm_fac/data_rescale(idx);
                
                %Optional: Transpose Result    
                if isfield(operator_pars{idx},'transpose_result')
                    for ss = 1:size(u{idx},3)
                        u{idx}(:,:,ss) = u{idx}(:,:,ss).'; 
                    end
                end
            end
             
            save([result_file,'_inter_',num2str(kk),'.mat'],'u','tgv','datafit','gap','par_in','operator_pars');
            clear u
        end
    end
    
    %Stopping rule evaluation
    if eval_opt
        tgv(kk+1) = get_tgv_joint(x,alph1,alph0,vec_norm,dx,dy,dz);
        datafit(kk+1,:) = get_datafit_rec(x,ld,f0,c0,K,operator_pars);
        gap(kk+1) = ( tgv(kk+1)+ sum(datafit(kk+1,:)) + gap_data_dual(r,ld,f0,c0,operator_pars) + gap_reg_dual(z,r,Kt,operator_pars,dx,dy,dz) )...
                /(n*m*k*noc);
    end
    
    %Show iteration
    
    if rem(kk,100) == 0;
        fprintf('\b\b\b\b\b\b\b\b\b\b');
        fprintf('%-10d',kk)    
    end

    %Iteration:
    if strcmp(stop_rule,'iteration')
        if kk>=stop_par
            cont = false;
        end
    end
end


%%Postprocessing##################################################################################################################
display(' ')
display(['Total number of iterations:       ',num2str(kk)])


%Set output
for idx = 1:noc

    u{idx} = squeeze(x(:,:,:,1,idx))*op_rescale(idx)*operator_pars{idx}.op_norm_fac/data_rescale(idx);

    %Optional: Transpose Result    
    if isfield(operator_pars{idx},'transpose_result')
        for ss = 1:size(u{idx},3)
            u{idx}(:,:,ss) = u{idx}(:,:,ss).'; 
        end
    end
end
            
%Write parameter-------------------------------
%Input: kk (iteration number)-------------------
psz = size(par_list,1);
for l=1:psz
    par_out{l,1} = par_list{l,1};
    eval(['par_out{l,2} = ',par_list{l,1},';'])
end
par_out{psz+1,1} = 'iteration_nr'; par_out{psz+1,2}=kk;
par_out{psz+2,1} = mfilename;
%Output: par-----------------------------------
%----------------------------------------------	

end %End of function
