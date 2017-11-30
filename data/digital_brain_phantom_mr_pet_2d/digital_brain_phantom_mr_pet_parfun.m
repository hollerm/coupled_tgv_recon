function [fout] = digital_brain_phantom_mr_pet_parfun(data,imgSens)

%Set scanner specific parameter
parm.SCANNERTYPE = 1;
parm.FWHM = 4.5;
parm.EMRECON_VERSION = EMrecon_Version();
parm.CONVOLVE = 1;
    
parm.PROJECTIONS = 176;
parm.ANGLES = parm.PROJECTIONS * 3;
parm.BIN_WIDTH = 1.0;
parm.RADIUS = 200.0;
% setup reconstruction volume
parm.FOV_X = 176.0;
parm.FOV_Y = parm.FOV_X;
parm.FOV_Z = 1.0;
parm.OFFSET_X = -parm.FOV_X/2.0;
parm.OFFSET_Y = -parm.FOV_Y/2.0;
parm.OFFSET_Z = -parm.FOV_Z/2.0;
parm.SIZE_X = 176;
parm.SIZE_Y = parm.SIZE_X;
parm.SIZE_Z = 1;
parm.VERBOSE = 0; % 0: No output. 1: Output
    
% setup mr image from Rawdata
parm.mr.size = [176,176,1];
parm.mr.reso = [1 1 1];

% We need this just to be consistent with in-vivo data
parm_MRFOV = parm;

%Rescale factor for PET operator
parm.pet_rescale = 1/298.1508; % Norm was estimated before

%Rescale factor for MRI operator
parm.mri_rescale = 1/2.8276; % Norm was estimated before

% Transpose the final result or not
% parm.transposeResult = 0;

%% Create operators
% MR
if strcmp(data.type,'mri')

    [nRO,nPE,nCh]=size(data.data);
    parm.FT = p2DFT(data.mask,[nRO,nPE],1,2); 

    display('Use 2D non-remap MR Operator')
    K = @(x) parm.mri_rescale*forward_operator_MR_2d(x,imgSens,parm.FT);
    Kt = @(y) parm.mri_rescale*adjoint_operator_MR_2d(y,imgSens,parm.FT);

    % Set operator parameter
    % Required
    operator_pars.ntype = 'gaussian';
    operator_pars.op_norm_fac = parm.mri_rescale;
    % Optional
    operator_pars.fwd_op = 'forward_operator_MR';
    operator_pars.adj_op = 'adjoint_operator_MR';
    operator_pars.FT = 'p2DFT';
    operator_pars.parm = parm;    

    % Set data
    f0 = data.data;

% PET
elseif strcmp(data.type,'pet')

    % Set operator parameter
    % Required
    operator_pars.ntype = 'poisson';
    operator_pars.op_norm_fac = parm.pet_rescale;
    %Optional
    operator_pars.parm = parm;
    
    K = @(x) parm.pet_rescale*EMrecon_ForwardProj_ND(parm_MRFOV,x);
    Kt = @(y) parm.pet_rescale*(EMrecon_BackProj_ND(parm_MRFOV,y));

    %Set data
    f0 = data.data;
    
end

%Set required output
fout.K = K;
fout.Kt = Kt;
fout.f0 = f0;

%Set optional output
fout.operator_pars = operator_pars;






