

function [gap] = gap_reg_dual(y,r,Kt,op_pars,dx,dy,dz)



    [n,m,k,~,noc] = size(y);
    Ksy = zeros(n,m,k,4,noc);
    
    for ii=1:noc
        Ksy(:,:,:,:,ii) = cat(4, -bdiv_3(y(:,:,:,1:3,ii),dx,dy,dz) + Kt{ii}(r{ii}) , -y(:,:,:,1:3,ii) - fdiv_3(y(:,:,:,4:9,ii),dx,dy,dz) );

        %Compensate for restriction of pet image to positive reals
        if strcmp(op_pars{ii}.ntype,'poisson')
            tmp = Ksy(:,:,:,1,ii);
            tmp(tmp<=0) = 0; %Also correct when using the nuclear norm
            Ksy(:,:,:,1,ii) = tmp;
        end            
        
    end

gap = sum(abs(Ksy(:)));



