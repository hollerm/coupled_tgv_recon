

function [gap] = gap_data_dual(r,ld,f0,c0,op_pars)


gap = 0;
noc = length(r);

for ii=1:noc
    
    u0 = f0{ii};
    u = r{ii};
    
    if strcmp(op_pars{ii}.ntype,'gaussian')
        
        tmp = sum( real(u0(:)).*real(u(:)) )  + sum( imag(u0(:)).*imag(u(:)) );

        gap = gap + tmp  + sum(abs(u(:)).^2)/(2*ld(ii));
        
        
        
    elseif strcmp(op_pars{ii}.ntype,'poisson')
        
        tmp = (u-ld(ii));
        

        %Mask of zero pet data
        zdata = (u0==0);



        if max(tmp(~zdata)) >= 0
            warning(['Pet variable above ld by ',num2str( max(u(:)) - ld(ii) ),' at ',num2str( nnz(tmp(~zdata)>=0) ) ,' positions with nonzero data'])
            
            display('-> Reducing pet variable')
            tmp( ~zdata & (tmp>-10^(-10)) ) = -10^(-10);    
        end
        
        
        tmp(~zdata) = -ld(ii)./tmp(~zdata);




        gap = gap + sum( ld(ii)*u0(~zdata).*( log(tmp(~zdata)) - 1) ) - sum(u(:).*c0{ii}(:));

        %display(['Gap data dual: ',num2str(gap)]);    
    end


end



    
    
