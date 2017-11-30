%%Author: Martin Holler
%%If you use this code (also outside the provided reconstruction utility) please cite:
%M. Holler, R. Huber, F. Knoll: Coupled regularization with multiple data discrepancies. Submitted, 2017.


%Input: z = (n,m,k,3,noc) where [n,m,k] are the spatial dimensions and noc is any number of components, scalar alpha
%Ouput: Input with singular values projected to [0,alpha]
function [z,nrm] = project_spectral_ball_3x3(z,alpha,norm_only);

%Set default to projection
if ~exist('norm_only') || isempty(norm_only)
    norm_only = 0;
end
nrm = 0;

verb = 0;

zthresh = 10*eps;

n = size(z,1);
m = size(z,2);
k = size(z,3);

%Matrix to be multiplied with z for projection
A11 = zeros(n,m,k);
A12 = zeros(n,m,k);
A13 = zeros(n,m,k);
A22 = zeros(n,m,k);
A23 = zeros(n,m,k);
A33 = zeros(n,m,k);


%Build matrix A = z*z'
%A =   [ a  r  t
%        r' b  d
%        t' d' c ]


time = cputime;
    a = sum(abs(z(:,:,:,1,:)).^2,5);
    b = sum(abs(z(:,:,:,2,:)).^2,5);
    c = sum(abs(z(:,:,:,3,:)).^2,5);

    r = sum(z(:,:,:,1,:).*conj(z(:,:,:,2,:)),5);
    t = sum(z(:,:,:,1,:).*conj(z(:,:,:,3,:)),5);
    d = sum(z(:,:,:,2,:).*conj(z(:,:,:,3,:)),5);

%disp(['-----Building A*A: ',num2str(cputime-time)]);    



%Preprocessing to set t to zero by H'*A*H with H = [c s 0 ; conj(s) -conj(c) 0 ; 0 0 1];
if verb
tic
end

    
    ndt = abs(d).^2 + abs(t).^2;
    maskndt = ndt > zthresh;
    
        %########BRANCH########
        time = cputime;
        %Get result directly for ~maskndt
        ld1 = c(~maskndt);
        if verb
        display(['Branch |d|+|t|=0: ',num2str(length(ld1))]);
        end
        U = zeros(length(ld1),3,3);
        U(:,1,1) = 0; U(:,2,1) = 0; U(:,3,1) = 1;

        [ld2,ld3,U(:,1:2,2),U(:,1:2,3)] = get_2d_eig_decomp(a(~maskndt),r(~maskndt),b(~maskndt));

        ld1 = sqrt(ld1);ld2 = sqrt(ld2);ld3 = sqrt(ld3);
        if norm_only
            nrm = nrm + sum( abs(ld1(:))+abs(ld2(:))+abs(ld3(:)));
        else
            %Insert 
            [A11,A12,A13,A22,A23,A33] = insert_UFU(ld1,ld2,ld3,U,A11,A12,A13,A22,A23,A33,~maskndt,alpha,verb);
            if verb
            disp(['-----Branch1: ',num2str(cputime-time)]);    
            end
        end
        %######################
        
    time = cputime;
    %Precomputation
    t1 = b.*t - d.*r;
    t2 = a.*d - t.*conj(r);
    
    %Update elements
    b = ( d.*(b.*conj(d) + r.*conj(t)) + t.*(a.*conj(t) + conj(d).*conj(r)) )./ndt;
    a = ( conj(t).*(t1) + conj(d).*(t2) )./ndt;
    r = ( d.*(t1) - t.*(t2) )./ndt;
    
    %Set elements of H
    ndt(~maskndt) = 1;
    ndt = sqrt(ndt);
    cs = conj(d)./ndt;
    ss = -t./ndt;
    
    d = ( - abs(d).^2 - abs(t).^2 )./ndt; %Never zero if t and d are not zero
    

    %Result:
    %A =   [ a  r  0
    %        r' b  d
    %        0  d' c ]
    if verb
    disp(['-----Elimination of t: ',num2str(cputime-time)]);    
    end
       
    %Case r = 0
    maskr = (abs(r)<zthresh) & maskndt;
    

        %########BRANCH########
        time = cputime;
        %Get result directly for maskr
        ld1 = a(maskr);
        if verb
        display(['Branch r=0: ',num2str(length(ld1))]);
        end
        U = zeros(length(ld1),3,3);
        U(:,1,1) = 1; U(:,2,1) = 0; U(:,3,1) = 0;
        [ld2,ld3,U(:,2:3,2),U(:,2:3,3)] = get_2d_eig_decomp(b(maskr),d(maskr),c(maskr));
        
        %Invert preprocessing 
        U(:,1:2,:) = cat(2,bsxfun(@times,U(:,1,:),cs(maskr)) + bsxfun(@times,U(:,2,:),ss(maskr)),bsxfun(@times,U(:,1,:),conj(ss(maskr))) - bsxfun(@times,U(:,2,:),conj(cs(maskr))));
        
        %Insert 
        ld1 = sqrt(ld1);ld2 = sqrt(ld2);ld3 = sqrt(ld3);
        if norm_only
            nrm = nrm + sum( abs(ld1(:))+abs(ld2(:))+abs(ld3(:)));
        else
            [A11,A12,A13,A22,A23,A33] = insert_UFU(ld1,ld2,ld3,U,A11,A12,A13,A22,A23,A33,maskr,alpha,verb);
            if verb
            disp(['-----Branch 2: ',num2str(cputime-time)]);    
            end
        end
        %######################
    
    time =cputime;
    %Reduce elements
    a = a(maskndt & ~maskr);
    b = b(maskndt & ~maskr);
    c = c(maskndt & ~maskr);
    r = r(maskndt & ~maskr);
    d = d(maskndt & ~maskr);
    
    
    cs = cs(maskndt & ~maskr);
    ss = ss(maskndt & ~maskr);
    if verb
    disp(['-----Element reduction: ',num2str(cputime-time)]);    
    end

   %Result:
    %A =   [ a  r  0
    %        r' b  d
    %        0  d' c ]
    % with d ~=0, r~=0

time = cputime;
%Get eigenvalues of A
    c2 = -a-b-c;
    c1 = a.*b + a.*c + b.*c - abs(r).^2 - abs(d).^2;
    c0 = a.*(abs(d).^2) +  c.*(abs(r).^2) - a.*b.*c ;

    p = c2.^2 - 3.*c1;
    q = -(27/2)*c0 - c2.^3 + (9/2)*c2.*c1;

    phi = sqrt(27*( (1/4)*c1.^2.*(p-c1) + c0.*(q+(27/4).*c0) ));

    mask = q<0;
    phi(~mask) = atan(phi(~mask)./q(~mask))./3;
    phi(mask) = (atan(phi(mask)./q(mask)) + pi)./3;
        
    sc1 = cos(phi+(2*pi./3));
    sc2 = cos(phi-(2*pi./3));
    
    x1 = -2*(sc1+sc2); % 2*cos(phi)*cos(2pi/3) = sc1 + sc2, cos(2pi/3) = -0.5
    x2 = 2*sc1;
    x3 = 2*sc2;
    %x1 = 2*cos(phi);
    %x2 = 2*cos(phi+(2*pi./3));
    %x3 = 2*cos(phi-(2*pi./3));

    
    %Eigenvalues
    p=sqrt(p)./3;
    c2 = c2./3;
    ld1 = p.*x1 - c2;
    ld2 = p.*x2 - c2;
    ld3 = p.*x3 - c2;
%    ld1 = (sqrt(p).*x1./3) - c2./3;
%    ld2 = (sqrt(p).*x2./3) - c2./3;
%    ld3 = (sqrt(p).*x3./3) - c2./3;




if verb
disp(['-----Eigenvalues: ',num2str(cputime-time)]);    
end


        if norm_only
            nrm = nrm + sum( abs(sqrt(ld1(:)))+abs(sqrt(ld2(:)))+abs(sqrt(ld3(:))));
        else    


        time = cputime;
        %Get eigenvectors
            U = zeros(length(ld1),3,3);
            %1st eigenvector. Note: Since we have a tridiagonal matrix with r and d nonzero, the first two columns must be
            %lin. indep., hence the eigenvector as below is well defined

                %U(:,:,1) = cross( [a-ld1,conj(r),zeros(length(ld1),1)],[r,b-ld1,conj(d)] ); %Slower than costum implementation
                U(:,:,1) = conj(cross_prod( a-ld1,conj(r),0,r,b-ld1,conj(d),3));

                %Normalize
                U(:,:,1) = bsxfun(@rdivide,U(:,:,1),sqrt(sum(abs(U(:,:,1)).^2,2)));

                
            %2nd eigenvector. Note: The third component of u1 cannot be zero, hence the two vectors below ar lin.indep. again
                U(:,:,2) = conj(cross_prod( a-ld2,conj(r),0,U(:,1,1),U(:,2,1),U(:,3,1),3));
                %Normalize
                U(:,:,2) = bsxfun(@rdivide,U(:,:,2),sqrt(sum(abs(U(:,:,2)).^2,2)));



            %3rd eigenvector. u1 and u2 are linearly independent and orthogonal:
                U(:,:,3) = conj( cross_prod( U(:,1,1),U(:,2,1),U(:,3,1),U(:,1,2),U(:,2,2),U(:,3,2) ,3) );

        if verb
        disp(['-----Eigenvectors: ',num2str(cputime-time)]);    
        end

        %Singular values
            ld1 = sqrt(ld1);
            ld2 = sqrt(ld2);
            ld3 = sqrt(ld3);


        time =cputime;
        %Invert preprocessing
            %H*U
            U(:,1:2,:) = cat(2,bsxfun(@times,U(:,1,:),cs) + bsxfun(@times,U(:,2,:),ss),bsxfun(@times,U(:,1,:),conj(ss)) - bsxfun(@times,U(:,2,:),conj(cs)));
        if verb
        disp(['-----Inverse preproc.: ',num2str(cputime-time)]);    
        end

        if verb
        toc
        end

        time = cputime;
        %Insert 
        [A11,A12,A13,A22,A23,A33] = insert_UFU(ld1,ld2,ld3,U,A11,A12,A13,A22,A23,A33,maskndt & ~maskr,alpha,verb);
        if verb
        disp(['-----Insertion: ',num2str(cputime-time)]);    
        end


        time = cputime;
        %(U*FU)*A
            z = cat(4,sum( bsxfun(@times,z,cat(4,A11,A12,A13)) , 4),sum( bsxfun(@times,z,cat(4,conj(A12),A22,A23)) , 4),sum( bsxfun(@times,z,cat(4,conj(A13),conj(A23),A33 )) , 4) );
        if verb
        disp(['-----Projection: ',num2str(cputime-time)]);    
        end
        
    end

end



function [A11,A12,A13,A22,A23,A33] = insert_UFU(ld1,ld2,ld3,U,A11,A12,A13,A22,A23,A33,mask,alpha,verb);

%Project
    %fac*U'
    time = cputime;
    fac = cat(2,min(1,alpha./ld1),min(1,alpha./ld2),min(1,alpha./ld3));
    FU = bsxfun(@times,permute(conj(U),[1 3 2]),fac);
    if verb
    disp(['-----Insertion1: ',num2str(cputime-time)]);    
    end
    %U*FU
        time = cputime;
    if size(U,1)>1 | size(U,1)==0
        a11 = sum(squeeze(U(:,1,:)).*FU(:,:,1),2);
        a12 = sum(squeeze(U(:,1,:)).*FU(:,:,2),2);
        a13 = sum(squeeze(U(:,1,:)).*FU(:,:,3),2);
        a22 = sum(squeeze(U(:,2,:)).*FU(:,:,2),2);
        a23 = sum(squeeze(U(:,2,:)).*FU(:,:,3),2);
        a33 = sum(squeeze(U(:,3,:)).*FU(:,:,3),2);
    else
        a11 = squeeze(FU(:,:,1))*squeeze(U(:,1,:));
        a12 = squeeze(FU(:,:,1))*squeeze(U(:,2,:));
        a13 = squeeze(FU(:,:,1))*squeeze(U(:,3,:));
        a22 = squeeze(FU(:,:,2))*squeeze(U(:,2,:));
        a23 = squeeze(FU(:,:,2))*squeeze(U(:,3,:));
        a33 = squeeze(FU(:,:,3))*squeeze(U(:,3,:));
    end

    if verb
    disp(['-----Insertion2: ',num2str(cputime-time)]);    
    end
    
    
    
    
%Insert
    time = cputime;
    A11(mask) = a11;
    A12(mask) = a12;
    A13(mask) = a13;
    A22(mask) = a22;
    A23(mask) = a23;
    A33(mask) = a33;
    if verb
    disp(['-----Insertion3: ',num2str(cputime-time)]);    
    end

end
