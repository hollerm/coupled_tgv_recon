
function [z] = prox_dual(z,alph1,alph0,norm_type)

    if strcmp(norm_type,'frob')
    
    
    
        max_norm = max(1, frob_norm_2(z(:,:,:,1:3,:))./alph1 );        
        z(:,:,:,1:3,:) = bsxfun(@rdivide,z(:,:,:,1:3,:),max_norm);
		

		max_norm = max(1,frob_norm_3(z(:,:,:,4:9,:))./alph0);
		z(:,:,:,4:9,:) = bsxfun(@rdivide,z(:,:,:,4:9,:),max_norm);
                        
                        
    elseif strcmp(norm_type,'nuc')
        
        % 3D data branch
        if size(z,3) > 1
        
            if size(z,5) > 2 %n-components
	            z(:,:,:,1:3,:) = project_spectral_ball_3x3(z(:,:,:,1:3,:), alph1);
	        else %2 components
	            [z(:,:,:,1:3,1), z(:,:,:,1:3,2)] = project_spectral_ball_2x2(z(:,:,:,1:3,1), z(:,:,:,1:3,2), alph1,4);
	        end
	    % 2D data branch
	    else
            
            [z(:,:,:,1,:), z(:,:,:,2,:)] = project_spectral_ball_2x2(z(:,:,:,1,:), z(:,:,:,2,:), alph1,5);
    
	    end
	    
	    
	    %Frobenius Prox of symetrized gradient
        max_norm = max(1,frob_norm_3(z(:,:,:,4:9,:))./alph0);
		z(:,:,:,4:9,:) = bsxfun(@rdivide,z(:,:,:,4:9,:),max_norm);
		
    elseif strcmp(norm_type,'sep')

		for jj=1:size(z,5)

			max_norm = max(1, frob_norm_2(z(:,:,:,1:3,jj))./alph1 );        
            z(:,:,:,1:3,jj) = bsxfun(@rdivide,z(:,:,:,1:3,jj),max_norm);
		
		    max_norm = max(1,frob_norm_3(z(:,:,:,4:9,jj))./alph0);
		    z(:,:,:,4:9,jj) = bsxfun(@rdivide,z(:,:,:,4:9,jj),max_norm);
		    
		end

        
    else
	    error('Wrong norm type for dual step');
    end

end



