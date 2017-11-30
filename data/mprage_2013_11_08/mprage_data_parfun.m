function [fout] = mprage_data_parfun(data)

%% Set dimensions------------------------------------------------

% Activate remapping of MR
parm.mr_remap = 1;

% Set scanner specific parameter
parm.SCANNERTYPE = 8;
parm.FWHM = 4.5;
parm.EMRECON_VERSION = EMrecon_Version();
parm.CONVOLVE = 1;
parm.MAX_RING_DIFF        = 60; % All PET data
parm.VERBOSE = 0;
parm.dispSlice =  120;

parm_MRFOV = parm;

%% Spatial image parameters of the two data sets
disp('Defining spatial image parameters of the two data sets');
% setup pet image
parm.pet.size = [344 344 127];              % fixed, should not be touched
parm.pet.reso = [2.08626 2.08626 2.03125];  % fixed, should not be touched
parm.pet.offset = [2.0865 0.72428 1.77867]; % AMIDE Dicom: [x,y,z], reconqueue recon
parm.pet.flip = [1 -1 -1];

% setup mr image from Rawdata
parm.mr.size = [190 190 208];
parm.mr.reso = [1.1458 1.1458 1.2000];
parm.mr.offset = [0.1091 21.8456 -2.4662]; % sPosition: [dSag,-dCor,dTraCenter]
parm.mr.flip = [1 -1 -1];
parm.mr.oversampledkspace = 0;

% offset adjustment
parm.mr.offset(3) = - parm.mr.offset(3);    % required since PET is flipped wrt to z

% Now adjust for gantry offset
parm.pet.gantry_offset_x = 0.3365;
parm.pet.gantry_offset_y = -0.7243;
parm.pet.gantry_offset_z = -2.7943;
parm.pet.horizontal_bed = 0;
parm.pet.vertical_bed = 0;
       
parm.pet.offset(1) = parm.pet.offset(1) + parm.pet.gantry_offset_x;
parm.pet.offset(2) = parm.pet.offset(2) + parm.pet.gantry_offset_y;
parm.pet.offset(3) = parm.pet.offset(3) + parm.pet.gantry_offset_z;

[parm] = calculate_common_grid(parm);

%% Build a PET Operator for MR FOV
parm_MRFOV.SIZE_X = parm.mr.size(1);
parm_MRFOV.SIZE_Y = parm.mr.size(2);
parm_MRFOV.SIZE_Z = parm.mr.size(3);
parm_MRFOV.FOV_X =  parm.mr.size(1).*parm.mr.reso(1);
parm_MRFOV.FOV_Y = parm.mr.size(2).*parm.mr.reso(2);
parm_MRFOV.FOV_Z = parm.mr.size(3).*parm.mr.reso(3);
parm_MRFOV.OFFSET_X = -parm_MRFOV.FOV_X/2.0;
parm_MRFOV.OFFSET_Y = parm_MRFOV.OFFSET_X;
parm_MRFOV.OFFSET_Z = -parm_MRFOV.FOV_Z/2.0;

% Rescale factor for PET operator
parm.pet_rescale = 1/130.4179; % Norm was estimated before

%Rescale factor for MRI operator
parm.mri_rescale = 1;

% Transpose the final result or not
parm.transposeResult = 1;

%% Create operators
% MR
if strcmp(data.type,'mri')

    [nRO,nPE,nSl,nCh]=size(data.data);
    parm.FT = p2DFT(data.mask,[nRO,nPE],1,2); 

    display('Using remap strategy for MR Operator')
    K = @(x) parm.mri_rescale*forward_operator_MR(x,data.imgSens,parm);
    Kt = @(y) parm.mri_rescale*adjoint_operator_MR(y,data.imgSens,parm);
    
    % Set operator parameter
    % Required
    operator_pars.ntype = 'gaussian';
    operator_pars.op_norm_fac = parm.mri_rescale;
    % Optional
    operator_pars.fwd_op = 'forward_operator_MR';
    operator_pars.adj_op = 'adjoint_operator_MR';
    operator_pars.FT = 'p2DFT';
    operator_pars.transpose_result = 1;
    operator_pars.parm = parm;    
    
    % Set data
    f0 = data.data;
    c0 = 0;
    
% PET
elseif strcmp(data.type,'pet')

    % Set operator parameter
    % Required
    operator_pars.ntype = 'poisson';
    operator_pars.op_norm_fac = parm.pet_rescale;
    % Optional
    operator_pars.parm = parm;
    operator_pars.transpose_result = 1;
    
    expmKu = exp(-data.K_mu);

    K = @(x) parm.pet_rescale*EMrecon_ForwardProj_ND(parm_MRFOV,flipdim(x,3)).*expmKu;
    Kt = @(y) parm.pet_rescale*flipdim(EMrecon_BackProj_ND(parm_MRFOV,y.*expmKu),3);

    % Set data
    f0 = data.data;
    c0 = data.c0;
end


%Set required output
fout.K = K;
fout.Kt = Kt;
fout.f0 = f0;

%Set optional output
fout.operator_pars = operator_pars;
fout.c0 = c0;



