function data = forward_operator_MR(img,sens,parm)
% Forward operator for multiple coils and multiple slices
% Last Change: 23/04/2015
% Florian Knoll (florian.knoll@nyumc.org)

[nRO,nPE,nSl,nCh]=size(sens);

% Undo alignment with PET before transforming back to k-space
img = zpad(img,parm.SIZE_X,parm.SIZE_Y,parm.SIZE_Z,1);
img = resample_image(parm.pet_common,parm.mr_common,img,'nearest');    
img = crop(img,parm.mr.size(1),parm.mr.size(2),parm.mr.size(3));

data = zeros(nRO,nPE,nSl,nCh);
for ss = 1:nSl
    % Re-Adjust rotation to make them consistent with PET images
    img(:,:,ss) = fliplr(flipud(img(:,:,ss))).';
        
    % Add zeros for oversampling in FE direction
    if parm.mr.oversampledkspace
        img_temp = zeros(nRO,nPE);
        padx = (nRO-size(img,1))/2;
        pad = zeros(padx,nPE);
        img_temp = cat(1,pad,img(:,:,ss),pad);
    else
       img_temp = img(:,:,ss);
    end
       
    for cc = 1:nCh
        data(:,:,ss,cc) = parm.FT*(img_temp.*sens(:,:,ss,cc));
    end
end
