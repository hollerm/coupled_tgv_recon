

function [x] = cross_prod(u1,u2,u3,v1,v2,v3,catdim);



x = cat(catdim,u2.*v3 - u3.*v2,u3.*v1 - u1.*v3,u1.*v2 - u2.*v1);


