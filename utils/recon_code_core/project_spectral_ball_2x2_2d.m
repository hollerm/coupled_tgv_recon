function [qx, qy,nrm,sigma1,sigma2] = project_spectral_ball_2x2_2d(px, py, lambda,norm_only)
%Note: The function has been adapted to deal with complex valued input, see handwritten notes and testsvd.m for a test


%Set default to projection
if ~exist('norm_only') || isempty(norm_only)
    norm_only = 0;
end


    P = size(px,3);

    % build A'*A
    a = sum(abs(px).^2, 3);
    b = sum(conj(px).*py, 3);
    c = sum(abs(py).^2, 3);
    

    % compute eigenvalues
    d = a + c;
    e = a - c;
    f = sqrt(4*abs(b).^2 + e.^2);
    
    
    sigma1 = sqrt((d+f)/2);
    sigma2 = sqrt((d-f)/2);
    
    if norm_only
        qx = 0;
        qy = 0;
        nrm = abs(sigma1) + abs(sigma2); %Due to numerical errors, the svds are sometimes complex with angle ~ e-18
    else
    
        nrm = 0;
    
        g = (e + f)/2;
        h = (e - f)/2;
    
        k = g.^2 + abs(b).^2;
        k(k == 0) = 1;
        l = h.^2 + abs(b).^2;
        l(l == 0) = 1;
    
        % modify singular values
        fac1 = 1./max(1, sigma1/lambda);
        fac2 = 1./max(1, sigma2/lambda);
        
        v1 = repmat(g, [1 1 P]).*px + repmat(conj(b), [1 1 P]).*py;
        v2 = repmat(h, [1 1 P]).*px + repmat(conj(b), [1 1 P]).*py;
        
        
        qx = repmat(fac1.*g, [1 1 P]).*v1./repmat(k, [1 1 P]) ...
            + repmat(fac2.*h, [1 1 P]).*v2./repmat(l, [1 1 P]);
        qy = repmat(fac1.*b, [1 1 P]).*v1./repmat(k, [1 1 P]) ...
            + repmat(fac2.*b, [1 1 P]).*v2./repmat(l, [1 1 P]);
            
    end
    
end
