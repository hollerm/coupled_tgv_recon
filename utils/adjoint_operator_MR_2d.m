function img = adjoint_operator_MR_2d(data,sens,FT)
% 
% Adjoint operator non-Cartesian imaging for multiple coils
% Last Change: 06/11/2013
% Florian Knoll (florian.knoll@nyumc.org)
% 

[nR,nC,nCh] = size(sens);
img = zeros(nR,nC);
for jj = 1:nCh
    img = img + conj(sens(:,:,jj)).*(FT'*(data(:,:,jj)));
end
% img = max(0,real(img));

%Slightly faster:
%img = sum( bsxfun(@times,conj(sens),FT'*bsxfun(@times,data,sqrt(w) ) ) , 3);

