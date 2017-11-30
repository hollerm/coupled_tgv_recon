%Function to calculate Frobenius-norm of a 3-vector of nxmx2 matrices
%Return value is nxm matrix with norm of elements

function [norm] = frob_norm_2(x)

    x = abs(x);
    norm = sqrt( sum(sum( x.^2 ,5),4) );


end
