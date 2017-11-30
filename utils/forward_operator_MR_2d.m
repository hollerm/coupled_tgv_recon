function data = forward_operator_MR_2d(img,sens,FT)
% 
% Forward operator non-Cartesian imaging for multiple coils
% 
% Last Change: 06/11/2013
% Florian Knoll (florian.knoll@nyumc.org)
% 


[~,~,nCh]=size(sens);
%data = zeros(numSamplesOnSpoke,numSpokes,nCh);
data = zeros(size(sens));

for ii = 1:nCh
    data(:,:,ii) = FT*(img.*sens(:,:,ii));
end

%Slightly faster:
%size(FT*bsxfun(@times,img,sens))

%data = bsxfun( @times,FT,bsxfun(@times,img,sens) );
