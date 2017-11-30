%Function to calculate the Jacobian of a function mapping from R^3 to any number of channels nch
%Input:             Function [n,m,o,nch] (or possibly [n,m,o,1,nch]
%Output:            Jacobian [n,m,o,3,nch]
%Differences:       Forward Differences [ u(i) = u(i+1) - u(i) ]
%Bdry extension:    Normal boundary extension [ u(n+1) = u(n) ]

%Remark: Even though the order [n,m,o,nch,3] seems more appealing as output order at first glance, we use [n,m,o,3,nch] since with this order the
%code can be applied both to matrices of size [n,m,o,nch] as well as [n,m,o]

function [grad] = fgrad_3(u,dx,dy,dz)

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
[n,m,o,nch] = size(u);%Deals with singelton fourth dimension
grad = zeros([n,m,o,3,nch]);

%Set (1,0,0) derivatives
grad(1:end-1,:,:,1,:) =	( u(2:end,:,:,:) - u(1:end-1,:,:,:) )./dx;

%Set (0,1,0) derivatives
grad(:,1:end-1,:,2,:) =	( u(:,2:end,:,:) - u(:,1:end-1,:,:) )./dy;

%Set (0,0,1) derivatives
grad(:,:,1:end-1,3,:) = (	u(:,:,2:end,:) - u(:,:,1:end-1,:) )./dz;


