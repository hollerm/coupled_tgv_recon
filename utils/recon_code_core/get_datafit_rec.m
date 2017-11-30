

function [z] = get_datafit_rec(u,ld,f0,c0,K,op_pars);


%Evaluate data fidelity
noc = size(u,5);

z = zeros(1,noc);

for ii=1:noc

    Ku = K{ii}(u(:,:,:,1,ii)) + c0{ii};
    u0 = f0{ii};

    if strcmp(op_pars{ii}.ntype,'gaussian')
		
	    z(ii) = (ld(ii)/2)* sum(abs(Ku(:) - u0(:) ).^2);
	
    elseif strcmp(op_pars{ii}.ntype,'poisson')



    %Get u0*log(u0/Ku) in a stable way
    
    %Mask of zero data
    zdata = (u0==0);
    
    %Problems occure only if Ku is zero where data is non-zero
	if min(Ku(~zdata))<=0
	    warning(['B(u)<=0 for u the pet image at positions with nonzero data: Min-value: ',num2str(min(Ku(:))),', n-of min: ',num2str(nnz(Ku(:) == min(Ku(:)) ))])
	    display('-> Cropping to 1e-10 for datafit evaluation...')
        mask = (Ku <=0) & (~zdata);
	    Ku(mask) = 10^(-10);
	end
	
    %Get second part of data fidelity only where data is non-zero
	tmp = zeros(size(u0));
    tmp(~zdata) = u0(~zdata).*log( u0(~zdata)./Ku(~zdata) );
    
    %Set data fidelity
	z(ii) = ld(ii)*(sum(    Ku(:) + tmp(:) ));
	
    else
	    error('Unknown noise type');
    end

end
