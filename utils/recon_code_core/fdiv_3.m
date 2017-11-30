%Function to calculate the divergence of a function mapping from R^3 to [6,nch], where nch is an arbitrary number of channels
%Input:             Function [n,m,o,6,nch]
%Output:            Div      [n,m,o,3,nch]
%Differences:       Forward Differences [ u(i) = u(i+1) - u(i) ]
%Bdry extension:    Zero boundary extension [ u(n+1) = 0 ]

%NOTE: This function is the negative adjoint to sym_bgrad_3.m


function [div]=d3_fdiv_2(v,dx,dy,dz)


%Set default steps if not given
if ~exist('dx') || isempty(dx)
    dx = 1;
end
if ~exist('dy') || isempty(dy)
    dy = 1;
end
if ~exist('dz') || isempty(dz)
    dz = 1;
end

%Allocate
div = zeros([size(v,1),size(v,2),size(v,3),3,size(v,5)]);


%dx(v1) + dy(v4) + dz(v5)
    %Dx of v1
    div(2:end-1,:,:,1,:) = ( v(3:end,:,:,1,:) - v(2:end-1,:,:,1,:) )./dx;
    div(1,      :,:,1,:) = ( v(2,    :,:,1,:)                      )./dx;
    div(end,    :,:,1,:) = (                  - v(end,    :,:,1,:) )./dx;


    %+Dy of v4
    div(:,2:end-1,:,1,:) = ( v(:,3:end,:,4,:) - v(:,2:end-1,:,4,:) )./dy + ...
    div(:,2:end-1,:,1,:);
    div(:,1,      :,1,:) = ( v(:,2,    :,4,:)                      )./dy + ...
    div(:,1,      :,1,:);
    div(:,end,    :,1,:) = (                  - v(:,end,    :,4,:) )./dy + ...
    div(:,end,    :,1,:);

    %+Dz of v5
    if size(v,3)>1
        div(:,:,2:end-1,1,:) = ( v(:,:,3:end,5,:) - v(:,:,2:end-1,5,:) )./dz + ...
        div(:,:,2:end-1,1,:);
        div(:,:,1,      1,:) = ( v(:,:,2,    5,:)                      )./dz + ...
        div(:,:,1,      1,:);
        div(:,:,end,    1,:) = (                  - v(:,:,end,    5,:) )./dz + ...
        div(:,:,end,    1,:);
    end

%dx(v4) + dy(v2) + dz(v6)
    %Dx of v4
    div(2:end-1,:,:,2,:) = ( v(3:end,:,:,4,:) - v(2:end-1,:,:,4,:) )./dx;
    div(1,      :,:,2,:) = ( v(2,    :,:,4,:)                      )./dx;
    div(end,    :,:,2,:) = (                  - v(end,    :,:,4,:) )./dx;

    %+Dy of v2
    div(:,2:end-1,:,2,:) = ( v(:,3:end,:,2,:) - v(:,2:end-1,:,2,:) )./dy + ...
    div(:,2:end-1,:,2,:);
    div(:,1,      :,2,:) = ( v(:,2,    :,2,:)                      )./dy + ...
    div(:,1,      :,2,:);
    div(:,end,    :,2,:) = (                  - v(:,end,    :,2,:) )./dy + ...
    div(:,end,    :,2,:);

    %+Dz of v6
    if size(v,3)>1
        div(:,:,2:end-1,2,:) = ( v(:,:,3:end,6,:) - v(:,:,2:end-1,6,:) )./dz + ...
        div(:,:,2:end-1,2,:);
        div(:,:,1,      2,:) = ( v(:,:,2,    6,:)                      )./dz + ...
        div(:,:,1,      2,:);
        div(:,:,end,    2,:) = (                  - v(:,:,end,    6,:) )./dz + ...
        div(:,:,end,    2,:);
    end

%dx(v5) + dy(v6) + dz(v3)
    if size(v,3)>1
        %Dx of v5
        div(2:end-1,:,:,3,:) = ( v(3:end,:,:,5,:) - v(2:end-1,:,:,5,:) )./dx;
        div(1,      :,:,3,:) = ( v(2,    :,:,5,:)                      )./dx;
        div(end,    :,:,3,:) = (                  - v(end,    :,:,5,:) )./dx;

        %+Dy of v6
        div(:,2:end-1,:,3,:) = ( v(:,3:end,:,6,:) - v(:,2:end-1,:,6,:) )./dy + ...
        div(:,2:end-1,:,3,:);
        div(:,1,      :,3,:) = ( v(:,2,    :,6,:)                      )./dy + ...
        div(:,1,      :,3,:);
        div(:,end,    :,3,:) = (                  - v(:,end,    :,6,:) )./dy + ...
        div(:,end,    :,3,:);

        %+Dz of v3
        div(:,:,2:end-1,3,:) = ( v(:,:,3:end,3,:) - v(:,:,2:end-1,3,:) )./dz + ...
        div(:,:,2:end-1,3,:);
        div(:,:,1,      3,:) = ( v(:,:,2,    3,:)                      )./dz + ...
        div(:,:,1,      3,:);
        div(:,:,end,    3,:) = (                  - v(:,:,end,    3,:) )./dz + ...
        div(:,:,end,    3,:);
    end




