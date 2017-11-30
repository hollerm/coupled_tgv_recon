%% Demo script for MRI-PET reconstruction using TGV regularization with frob/sep/nuc norm coupling.
% -------------------------
% florian.knoll@nyumc.org
% martin.holler@uni-graz.at
% kristian.bredies@uni-graz.at
% thomas.koesters@nyumc.org
% 21.8.2017
% -------------------------
% If you consider this code to be useful for your research, please cite [1,2].
% 
% [1] F Knoll, M Holler, T Koesters, R Otazo, K Bredies, DK Sodickson:
% Joint MR-PET reconstruction using a multi-channel image regularizer.
% IEEE Transactions on Medical Imaging 36: 1-16 (2017)
% (http://ieeexplore.ieee.org/document/7466848/)
%
% [2] M. Holler, R. Huber, F. Knoll: Coupled regularization with multiple data discrepancies.
% Submitted, 2017.

clear all; close all; clc;
iptsetpref('ImshowBorder','loose');

%% Add to path
addpath(genpath('./utils/'));

%% Reconstruction parameters
% Please see readme.txt and publication [1] for an explanation and the
% recommended settings for the reconstruction hyperparameters. Experiment
% specific parameters are set directly at the test case switch
recon_pars = {  'vec_norm','nuc';...    % Vector norm: Options sep (no coupling), frob (frobenius norm coupling), nuc (nuclear norm coupling) are available
                'stop_par',1000;...     % Number of iterations: 1000 for guaranteed convergence, in practice no real change happens after 500
                'eval_opt',1;...        % Evaluate approximate primal-dual gap (costly)
                'save_inter',0;...      % Save intermediate reconstructions, = 0 means don't save, >0 gives rate
                's_t_ratio',1;...       % Primal-dual step ratio: Can be left unchanged
};

%% Select test case
testcase = 'digital_brain_phantom_mr_denoising'; % Use ld=1e-4 for nuclear norm coupling and (no) operator scaling of [1 1]
%testcase = 'digital_brain_phantom_mr_pet_2d';    % Use ld=1e-3 for nuclear norm coupling and (no) operator scaling of [1 1]
%testcase = 'mprage_2013_11_08';                  % Use ld=[10,90] for nuclear norm coupling and the default operator scaling of [3 10]


% List of test cases
switch testcase
    % simulation 2d denoising
    case 'digital_brain_phantom_mr_denoising'
        script_pars.datapath = './data/digital_brain_phantom_mr_denoising/';
        script_pars.source_data = 'digital_brain_phantom_mrsnr100';
        script_pars.result_path = 'results/digital_brain_phantom_mr_denoising/';
        recon_pars(end+1,:) = {'ld',[1e-4]};        % Operator is identity, use only a single parameter here, data terms are weighted automatically according to their range  
        recon_pars(end+1,:) = {'op_rescale',[1 1]}; % no rescaling
    % simulation 2d PET-MRI reconstruction
    case 'digital_brain_phantom_mr_pet_2d'
        script_pars.datapath = './data/digital_brain_phantom_mr_pet_2d/';
        script_pars.source_data = 'digital_brain_phantom_mr_pet_R4_mrsnr200_min10';
        script_pars.result_path = 'results/digital_brain_phantom_mr_pet_2d/';
        recon_pars(end+1,:) = {'ld',[1e-3]};        % Use only a single parameter here, data terms are weighted automatically according to their range  
        recon_pars(end+1,:) = {'op_rescale',[1 1]}; % no rescaling
    % in-vivo PET-MRI reconstruction
    case 'mprage_2013_11_08' 
        script_pars.datapath = './data/mprage_2013_11_08/';
        script_pars.source_data = 'rawdata_mprage_fdg_2013';
        if ~exist([script_pars.datapath,script_pars.source_data,'.mat'])
            error('3D patient data not in repository. Please download the data at http://cai2r.net/resources/software')
        end
        script_pars.result_path = 'results/results_mprage_2013_11_08/';
        recon_pars(end+1,:) = {'ld',[1,90]};        % Regularization parameters for MR and PET operators
        recon_pars(end+1,:) = {'op_rescale',[3 10]}; % Rescaling for [Gaussian, Poisson] operator. Default for real PET-MRI measurement operators [3 10]
end

%Generate filename
script_pars.outfilename = [script_pars.source_data,'_',recon_pars{1,2},'_ld',num2str(recon_pars{6,2}),'_iter',num2str(recon_pars{2,2})];
script_pars.outfilename = strrep(script_pars.outfilename,'.','_');

%Add filename to input for intermediate saving
recon_pars(end+1,:) = {'result_file',[script_pars.result_path,script_pars.outfilename]};


%% TGV Reconstruction
display('Performing TGV reconstruction...');

% Add datapath
addpath(script_pars.datapath);

% Set result paths
if ~exist(script_pars.result_path,'dir'); mkdir(script_pars.result_path); end

% Do reconstruction
t1 = clock;
[u,tgv,datafit,gap,par_out,operator_pars] = tgv_recon_joint(script_pars.source_data,recon_pars);
compTime = etime(clock,t1);
disp(['Computation time: ', num2str(compTime/60),' min']);

%% Display
if strcmp(testcase,'digital_brain_phantom_mr_denoising')
    load([script_pars.datapath,script_pars.source_data]);
    data_disp = cat(2,normalize01(data{1}.data),normalize01(data{2}.data),normalize01(data{3}.data));
    figure,imshow(data_disp,[]);
    title('noisy input');
    u_disp = cat(2,normalize01(u{1}),normalize01(u{2}),normalize01(u{3}));
    figure,imshow(u_disp,[]);
    title(['tgv ', recon_pars{1,2}, ' denoised']);
elseif strcmp(testcase,'digital_brain_phantom_mr_pet_2d')
    u_disp = cat(2,normalize01(u{1}),normalize01(u{2}),normalize01(u{3}),normalize01(u{4}));
    figure,imshow(u_disp,[]); title(['mr-pet tgv ', recon_pars{1,2}]);
else
    figure,imshow(abs(u{1}(:,:,110)),[]); title(['mr tgv ', recon_pars{1,2}]);
    figure,imshow(abs(u{2}(:,:,110)),[]); title(['pet tgv ', recon_pars{1,2}]);
end

%% Save results
save([script_pars.result_path,script_pars.outfilename],'u','tgv','datafit','gap','script_pars','par_out','operator_pars')

% Remove datapath
rmpath(script_pars.datapath);


