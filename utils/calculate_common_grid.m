function [parm] = calculate_common_grid(parm)

%% setup PET image
parm_IN.SIZE_X = parm.pet.size(1);
parm_IN.SIZE_Y = parm.pet.size(2);
parm_IN.SIZE_Z = parm.pet.size(3);

parm_IN.FOV_X = parm_IN.SIZE_X * parm.pet.reso(1);
parm_IN.FOV_Y = parm_IN.SIZE_Y * parm.pet.reso(2);
parm_IN.FOV_Z = parm_IN.SIZE_Z * parm.pet.reso(3);

parm_IN.FLIP_X = parm.pet.flip(1);
parm_IN.FLIP_Y = parm.pet.flip(2);
parm_IN.FLIP_Z = parm.pet.flip(3);

parm_IN.INTP_OFFSET_X = parm.pet.offset(1);
parm_IN.INTP_OFFSET_Y = parm.pet.offset(2);
parm_IN.INTP_OFFSET_Z = parm.pet.offset(3);

%% setup MR image
parm_OUT.SIZE_X = floor(parm_IN.FOV_X / parm.mr.reso(1));
parm_OUT.SIZE_Y = floor(parm_IN.FOV_Y / parm.mr.reso(2));
parm_OUT.SIZE_Z = floor(parm_IN.FOV_Z / parm.mr.reso(3));

% make sure that corresponding image dims are both even / odd
if (mod(parm.mr.size(1),2) ~= mod(parm_OUT.SIZE_X,2))
    parm_OUT.SIZE_X = parm_OUT.SIZE_X + 1;
end

if (mod(parm.mr.size(2),2) ~= mod(parm_OUT.SIZE_Y,2))
    parm_OUT.SIZE_Y = parm_OUT.SIZE_Y + 1;
end

if (mod(parm.mr.size(3),2) ~= mod(parm_OUT.SIZE_Z,2))
    parm_OUT.SIZE_Z = parm_OUT.SIZE_Z + 1;
end

parm_OUT.FOV_X = parm_OUT.SIZE_X * parm.mr.reso(1);
parm_OUT.FOV_Y = parm_OUT.SIZE_Y * parm.mr.reso(2);
parm_OUT.FOV_Z = parm_OUT.SIZE_Z * parm.mr.reso(3);

parm_OUT.FLIP_X = parm_IN.FLIP_X;
parm_OUT.FLIP_Y = parm_IN.FLIP_Y;
parm_OUT.FLIP_Z = parm_IN.FLIP_Z;

parm_OUT.INTP_OFFSET_X = parm_IN.INTP_OFFSET_X;
parm_OUT.INTP_OFFSET_Y = parm_IN.INTP_OFFSET_Y;
parm_OUT.INTP_OFFSET_Z = parm_IN.INTP_OFFSET_Z;

%% new PET parameter
parm.SIZE_X = parm_OUT.SIZE_X;
parm.SIZE_Y = parm_OUT.SIZE_Y;
parm.SIZE_Z = parm_OUT.SIZE_Z;

parm.FOV_X = parm_OUT.FOV_X;
parm.FOV_Y = parm_OUT.FOV_Y;
parm.FOV_Z = parm_OUT.FOV_Z;

parm.OFFSET_X = - parm.FOV_X / 2.0;
parm.OFFSET_Y = - parm.FOV_Y / 2.0;
parm.OFFSET_Z = - parm.FOV_Z / 2.0;

parm.OFFSET_X = parm.OFFSET_X - parm.pet.flip(1);
parm.OFFSET_Y = parm.OFFSET_Y - parm.pet.flip(2);
% there is a "-" due to the flip wrt z
parm.OFFSET_Z = parm.OFFSET_Z - parm.pet.flip(3);

%% For update of parameter files
parm.pet_common = parm_OUT;

parm.mr_common = parm_OUT;
parm.mr_common.FLIP_X = parm.mr.flip(1);
parm.mr_common.FLIP_Y = parm.mr.flip(2);
parm.mr_common.FLIP_Z = parm.mr.flip(3);
parm.mr_common.INTP_OFFSET_X = parm.mr.offset(1);
parm.mr_common.INTP_OFFSET_Y = parm.mr.offset(2);
parm.mr_common.INTP_OFFSET_Z = parm.mr.offset(3);

end