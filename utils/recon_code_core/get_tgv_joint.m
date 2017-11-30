

function [tgv] = get_tgv_joint(x,alph1,alph0,norm_type,dx,dy,dz)


        %Calculate TGV norm of current iterate
        x0 =  cat(4,fgrad_3(x(:,:,:,1,:),dx,dy,dz)-x(:,:,:,2:4,:),sym_bgrad_3(x(:,:,:,2:4,:),dx,dy,dz));

	
        if strcmp(norm_type,'frob')                
        	tgv = alph1*(sum(sum(sum( frob_norm_2( x0(:,:,:,1:3,:) ) )))) + alph0*(sum(sum(sum( frob_norm_3( x0(:,:,:,4:9,:) ) ))));
        	
        elseif strcmp(norm_type,'nuc')

             % 3D data branch
            if size(x,3) > 1
        
                if size(x,5) > 2 %n-components
	                [~,nuc_norm] = project_spectral_ball_3x3(x0(:,:,:,1:3,:), alph1,1);
	            else %2 components
	                [~,~,nuc_norm] = project_spectral_ball_2x2(x0(:,:,:,1:3,1), x0(:,:,:,1:3,2), alph1,4,1);
	            end
	        % 2D data branch
	        else
            
            %Squeeze singelton dimensions and project
	        [~,~,nuc_norm] = project_spectral_ball_2x2(x0(:,:,:,1,:), x0(:,:,:,2,:), alph1,5,1);

            	        
	    end
            


        	tgv = alph1*sum( nuc_norm(:) ) + alph0*(sum(sum(sum( frob_norm_3( x0(:,:,:,4:9,:) ) ))));
        	
        	
       elseif strcmp(norm_type,'sep')
       
            tgv = 0;
            for ii=1:size(x,5);
                tgv = tgv + alph1*(sum(sum(sum( frob_norm_2( x0(:,:,:,1:3,ii) ) )))) + alph0*(sum(sum(sum( frob_norm_3( x0(:,:,:,4:9,ii) ) ))));
            end
            
        else
        	error('Wrong dual norm type');
        end
        

	
end
