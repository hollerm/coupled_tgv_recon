%%Authors:
%%Kristian Bredies and Martin Holler

function [qx, qy,nrm,sigma1,sigma2] = project_spectral_ball_2x2(px, py, lambda,vdim,norm_only)
%Note: The function has been adapted to deal with complex valued input, see handwritten notes and testsvd.m for a test



%Set default to projection
if ~exist('norm_only') || isempty(norm_only)
    norm_only = 0;
end
if ~exist('vdim') || isempty(vdim)
    vdim=4;
    warning('project_spectral_ball_2x2: Dimension not given, assuming 3 spatial dimensions')
end


%The function assumes vdim to be the vector-dimension of px, py, e.g. vdim=3 says px,py have two spatial dimensions

    P = size(px,vdim);
    repmatdims = ones(1,vdim);
    repmatdims(end) = P;

    
    % build A'*A
    a = sum(abs(px).^2, vdim);
    b = sum(conj(px).*py, vdim);
    c = sum(abs(py).^2, vdim);
    
    

    % compute eigenvalues
    d = a + c;
    e = a - c;
    f = sqrt(4*abs(b).^2 + e.^2);
    
    %Set singular values    
    sigma1 = sqrt((d+f)/2);
    sigma2 = sqrt((d-f)/2);


    if norm_only
        qx = 0;
        qy = 0;
        nrm = abs(sigma1) + abs(sigma2); %Due to numerical errors, the svds are sometimes complex with angle ~ e-18
    else
    
        nrm = 0;
        
        %Modify singular values
        fac1 = 1./max(1, sigma1/lambda);
        fac2 = 1./max(1, sigma2/lambda);

        
        %Compute entries of V such that A = UDV*
        g = (e + f)/2;
        h = (e - f)/2;
        
        k = g.^2 + abs(b).^2;
        %k(k == 0) = 1;
        l = h.^2 + abs(b).^2;
        %l(l == 0) = 1;
        
      
        %Get AV (without normalization)
        v1 = repmat(g, repmatdims).*px + repmat(conj(b), repmatdims).*py;
        v2 = repmat(h, repmatdims).*px + repmat(conj(b), repmatdims).*py;
        
        %Get (AV) (SV*) (where S= diag(fac1,fac2) is the prox on the diagonal singular value matrix)
        qx = repmat(fac1.*g, repmatdims).*v1./repmat(k, repmatdims) ...
            + repmat(fac2.*h, repmatdims).*v2./repmat(l, repmatdims);
        qy = repmat(fac1.*b, repmatdims).*v1./repmat(k, repmatdims) ...
            + repmat(fac2.*b, repmatdims).*v2./repmat(l, repmatdims);


        %Correction for special case of orthogonal columns in A (see handwritten notes)
            idb = abs(b) == 0; %Case of orthogonal columns in original matrix
            idbr =  repmat( idb,repmatdims); %Vector extension of above
            idpx = sum(abs(px),vdim)==0; %Subcase of first vector being zero
     
            %Only if px=0, fac need to be modified        
            fac2(idpx) = fac1(idpx);
            fac1(idpx) = 1;
            
            
            %Prox in case of orthogonal columns of A      
            if size(px,1)*size(px,2)>1
                qx(idbr) = px(idbr).*repmat(fac1(idb),[P 1]);
                qy(idbr) = py(idbr).*repmat(fac2(idb),[P 1]);
            elseif ~isempty(px(idbr)) %special case for 1x1 array
                qx(idbr) = px(idbr).*repmat(fac1(idb),repmatdims);
                qy(idbr) = py(idbr).*repmat(fac2(idb),repmatdims);
            end
            
    end
    
end
