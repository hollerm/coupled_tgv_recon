function img = adjoint_operator_MR(data,sens,parm)
% Adjoint operator for multiple coils and multiple slices
% Last Change: 23/04/2015
% Florian Knoll (florian.knoll@nyumc.org)

[nRO,nPE,nSl,nCh] = size(sens);
for ss = 1:nSl
    img_temp = zeros(nRO,nPE);
    for cc = 1:nCh
        img_temp = img_temp + conj(sens(:,:,ss,cc)).*(parm.FT'*(data(:,:,ss,cc)));
    end
    
    % Remove Oversampling in FE direction
    if parm.mr.oversampledkspace
        img(:,:,ss) = img_temp(nRO/4+1:nRO*3/4,:,:);
    else
        img(:,:,ss) = img_temp;
    end
    
    % Adjust rotation to make them consistent with PET images
    img(:,:,ss) = flipud(fliplr(img(:,:,ss).'));
    % img(:,:,ss) = img(:,:,ss).';
end

% Align with PET
img = zpad(img,parm.SIZE_X,parm.SIZE_Y,parm.SIZE_Z,1);
img = resample_image(parm.mr_common,parm.pet_common,img,'nearest');
img = crop(img,parm.mr.size(1),parm.mr.size(2),parm.mr.size(3));
