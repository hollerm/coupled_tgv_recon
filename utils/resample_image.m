function [resampled_image] = resample_image(parm_IN,parm_OUT,image,interpolation)

if (nargin == 3)
    interpolation = 'linear';
end

dx_in = parm_IN.FOV_X / parm_IN.SIZE_X;
dy_in = parm_IN.FOV_Y / parm_IN.SIZE_Y;
dx_out = parm_OUT.FOV_X / parm_OUT.SIZE_X;
dy_out = parm_OUT.FOV_Y / parm_OUT.SIZE_Y;

if (ndims(image) == 3)
    dz_in = parm_IN.FOV_Z / parm_IN.SIZE_Z;
    dz_out = parm_OUT.FOV_Z / parm_OUT.SIZE_Z;
end

if (parm_IN.FLIP_X == -1)
    image = flipdim(image,1);
end

if (parm_IN.FLIP_Y == -1)
    image = flipdim(image,2);
end

if (ndims(image) == 3)
    if (parm_IN.FLIP_Z == -1)
        image = flipdim(image,3);
    end
end

if (mod(parm_IN.SIZE_X,2) == 0)
    xgrid_in = parm_IN.INTP_OFFSET_X-(parm_IN.SIZE_X/2.0-0.5)*dx_in:dx_in:parm_IN.INTP_OFFSET_X+(parm_IN.SIZE_X/2.0-0.5)*dx_in;
else
    xgrid_in = parm_IN.INTP_OFFSET_X-((parm_IN.SIZE_X-1)/2.0)*dx_in:dx_in:parm_IN.INTP_OFFSET_X+((parm_IN.SIZE_X-1)/2.0)*dx_in;
end

if (mod(parm_IN.SIZE_Y,2) == 0)
    ygrid_in = parm_IN.INTP_OFFSET_Y-(parm_IN.SIZE_Y/2.0-0.5)*dy_in:dy_in:parm_IN.INTP_OFFSET_Y+(parm_IN.SIZE_Y/2.0-0.5)*dy_in;
else
    ygrid_in = parm_IN.INTP_OFFSET_Y-((parm_IN.SIZE_Y-1)/2.0)*dy_in:dy_in:parm_IN.INTP_OFFSET_Y+((parm_IN.SIZE_Y-1)/2.0)*dy_in;
end

if (ndims(image) == 3)
    if (mod(parm_IN.SIZE_Z,2) == 0)
        zgrid_in = parm_IN.INTP_OFFSET_Z-(parm_IN.SIZE_Z/2.0-0.5)*dz_in:dz_in:parm_IN.INTP_OFFSET_Z+(parm_IN.SIZE_Z/2.0-0.5)*dz_in;
    else
        zgrid_in = parm_IN.INTP_OFFSET_Z-((parm_IN.SIZE_Z-1)/2.0)*dz_in:dz_in:parm_IN.INTP_OFFSET_Z+((parm_IN.SIZE_Z-1)/2.0)*dz_in;
    end
end

if (mod(parm_OUT.SIZE_X,2) == 0)
    xgrid_out = parm_OUT.INTP_OFFSET_X-(parm_OUT.SIZE_X/2.0-0.5)*dx_out:dx_out:parm_OUT.INTP_OFFSET_X+(parm_OUT.SIZE_X/2.0-0.5)*dx_out;
else
    xgrid_out = parm_OUT.INTP_OFFSET_X-((parm_OUT.SIZE_X-1)/2.0)*dx_out:dx_out:parm_OUT.INTP_OFFSET_X+((parm_OUT.SIZE_X-1)/2.0)*dx_out;
end

if (mod(parm_OUT.SIZE_Y,2) == 0)
    ygrid_out = parm_OUT.INTP_OFFSET_Y-(parm_OUT.SIZE_Y/2.0-0.5)*dy_out:dy_out:parm_OUT.INTP_OFFSET_Y+(parm_OUT.SIZE_Y/2.0-0.5)*dy_out;
else
    ygrid_out = parm_OUT.INTP_OFFSET_Y-((parm_OUT.SIZE_Y-1)/2.0)*dy_out:dy_out:parm_OUT.INTP_OFFSET_Y+((parm_OUT.SIZE_Y-1)/2.0)*dy_out;
end

if (ndims(image) == 3)
    if (mod(parm_OUT.SIZE_Z,2) == 0)
        zgrid_out = parm_OUT.INTP_OFFSET_Z-(parm_OUT.SIZE_Z/2.0-0.5)*dz_out:dz_out:parm_OUT.INTP_OFFSET_Z+(parm_OUT.SIZE_Z/2.0-0.5)*dz_out;
    else
        zgrid_out = parm_OUT.INTP_OFFSET_Z-((parm_OUT.SIZE_Z-1)/2.0)*dz_out:dz_out:parm_OUT.INTP_OFFSET_Z+((parm_OUT.SIZE_Z-1)/2.0)*dz_out;
    end
end

if (ndims(image) == 3)
    [X_in,Y_in,Z_in] = meshgrid(ygrid_in,xgrid_in,zgrid_in);
    [X_out,Y_out,Z_out] = meshgrid(ygrid_out,xgrid_out,zgrid_out);
    
    resampled_image = interp3(X_in,Y_in,Z_in,image,X_out,Y_out,Z_out,interpolation);
else
    [X_in,Y_in] = meshgrid(ygrid_in,xgrid_in);
    [X_out,Y_out] = meshgrid(ygrid_out,xgrid_out);
    
    resampled_image = interp2(X_in,Y_in,image,X_out,Y_out,interpolation);
end

if (parm_OUT.FLIP_X == -1)
    resampled_image = flipdim(resampled_image,1);
end

if (parm_OUT.FLIP_Y == -1)
    resampled_image = flipdim(resampled_image,2);
end

if (ndims(image) == 3)
    if (parm_OUT.FLIP_Z == -1)
        resampled_image = flipdim(resampled_image,3);
    end
end

end