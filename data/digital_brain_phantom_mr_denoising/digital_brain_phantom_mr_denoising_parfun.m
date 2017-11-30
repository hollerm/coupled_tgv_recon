function [fout] = digital_brain_phantom_mr_denoising_parfun(data)

%Size
% img_size = [data{1}.nR,data{1}.nC]; %Must be even

if strcmp(data.type,'mri')
    %Set operator parameter
    operator_pars.ntype = 'gaussian';
    
    %Set operator
    K = @(x) x;
    Kt = @(x) x;
    
    %Set data
    f0 = data.data;

        
elseif strcmp(data.type,'pet')

    %Set operator parameter
    operator_pars.ntype = 'poisson';

    K = @(x) x;
    Kt = @(x) x;
    
    %Set data
    f0 = data.data;

end


%Set required output
fout.K = K;
fout.Kt = Kt;
fout.f0 = f0;

%Set optional output
fout.operator_pars = operator_pars;
