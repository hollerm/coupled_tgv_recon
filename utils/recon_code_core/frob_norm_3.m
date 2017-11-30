%Function to calculate Frobenius-norm of a 3-vector of nxmx3 matrices
%Return value nxm matrix containing norm of elements

function [norm] = frob_norm_3(x)

x = abs(x);

norm = sqrt( sum( sum( x(:,:,:,1:3,:).^2 ,4) + 2*sum( x(:,:,:,4:6,:).^2 ,4), 5) );
              
end
