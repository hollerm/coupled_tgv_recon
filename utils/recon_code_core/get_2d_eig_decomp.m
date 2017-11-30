

function [sigma1,sigma2,v1,v2] = get_2d_eig_decomp(a,b,c);
zthresh = 10*eps;

    % compute eigenvalues
    d = a + c;
    e = a - c;
    f = sqrt(4*abs(b).^2 + e.^2);
    
    sigma1 = (d+f)/2;
    sigma2 = (d-f)/2;

    %Compute eigenvectors
    g = (e + f)/2;
    h = (e - f)/2;
    
    k = sqrt(g.^2 + abs(b).^2);
%    k(k < zthresh) = 1;
    l = sqrt(h.^2 + abs(b).^2);
 %   l(l < zthresh) = 1;
    
    
   
    %Get eigenvectors
    v1 = bsxfun(@rdivide,cat(2,g,conj(b)),k);
    v2 = bsxfun(@rdivide,cat(2,h,conj(b)),l);
    
    zmask = abs(b)<zthresh;
    sigma1(zmask) = a(zmask);
    sigma2(zmask) = c(zmask);
    v1([zmask,false(length(zmask),1)]) = 1;
    v1([false(length(zmask),1),zmask,]) = 0;
    v2([zmask,false(length(zmask),1)]) = 0;
    v2([false(length(zmask),1),zmask,]) = 1;

    

