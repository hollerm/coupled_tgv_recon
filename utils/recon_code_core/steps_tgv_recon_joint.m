%Function to get the adaptive stepsize for tgv mri-pet reconstruction

function [sig,tau] = steps_tgv_recon_joint(x,sig,tau,s_t_ratio,K,dx,dy,dz)

	[n,m,k,~,noc] = size(x);
	
	
	%Get Kx
	%Kx = zeros(n,m,k,5,2);    
	Kx = cat(4, fgrad_3(x(:,:,:,1,:),dx,dy,dz) - x(:,:,:,2:4,:), sym_bgrad_3(x(:,:,:,2:4,:),dx,dy,dz) );
	
	Kx2 = 0;
	for ii=1:noc
	    tmp = abs(K{ii}(x(:,:,:,1,ii))).^2;
	    Kx2 = Kx2 + sum(tmp(:));
	end
    

    	%Get |Kx|
    	nKx = sqrt(	sum(sum(sum( ...
    	            sum( sum( abs(Kx(:,:,:,1:6,:)).^2 , 4 ) + 2*sum( abs(Kx(:,:,:,7:9,:)).^2 , 4 ) , 5) ...
    	            ))) + Kx2 );
    	
    	%Get |x|
    	nx = sqrt(sum( abs(x(:)).^2 ));

    %Set |x| / |Kx|
    tmp = (nx/nKx);
    theta = 0.95;
    
    %Check convergence condition
    if sig*tau > tmp^2 %If stepsize is too large
        if theta^(2)*sig*tau < tmp^2 %Check if minimal decrease satisfies condition
            sig = theta*sig;
            tau = theta*tau;
        else                        %If not, decrease further
            sig = tmp*s_t_ratio;
            tau = tmp/s_t_ratio;
        end
    end
    
    %sig = tmp;
    %tau = tmp;
    
end
