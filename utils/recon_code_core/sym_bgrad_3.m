%Function to calculate the Symmetrized gradient of a function mapping from R^3 to [3,nch], where nch is an arbitrary number of channels
%Input:             Function [n,m,o,3,nch]
%Output:            Symgrad  [n,m,o,6,nch] [Order of derivatives: (1,1),(2,2),(3,3),(1,2),(1,3),(2,3)]
%Differences:       Backward Differences [ u(i) = u(i) - u(i-1) ]
%Bdry extension:    Normal boundary extension [ u(n+1) = u(n) ]


function [grad]=sym_bgrad_3(v,dx,dy,dz)

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
grad = zeros([size(v,1),size(v,2),size(v,3),6,size(v,5)]);


%%Diagonal entries--------------------------------------------------------

	grad(2:end,:,:,1,:) = ( v(2:end,:,:,1,:) - v(1:end-1,:,:,1,:) )./dx; %Get dx of v1


	
	grad(:,2:end,:,2,:) = ( v(:,2:end,:,2,:) - v(:,1:end-1,:,2,:) )./dy; %Get dy of v2
	

	grad(:,:,2:end,3,:) = ( v(:,:,2:end,3,:) - v(:,:,1:end-1,3,:) )./dz; %Get dz of v3


%%Off-diagonal entries----------------------------------------------------

    %(1,2) entry
    grad(:,2:end,:,4,:) = ( v(:,2:end,:,1,:) - v(:,1:end-1,:,1,:) )./dy; %Get dy of v1
    
    grad(2:end,:,:,4,:) = ( v(2:end,:,:,2,:) - v(1:end-1,:,:,2,:) )./dx + ... %Add dx of v2
    grad(2:end,:,:,4,:);
    
    grad(:,:,:,4,:) = grad(:,:,:,4,:)./2;
    
    %(1,3) entry
    grad(:,:,2:end,5,:) = ( v(:,:,2:end,1,:) - v(:,:,1:end-1,1,:) )./dz; %Get dz of v1
    
    grad(2:end,:,:,5,:) = ( v(2:end,:,:,3,:) - v(1:end-1,:,:,3,:) )./dx + ... %Add dx of v3
    grad(2:end,:,:,5,:);
    
    grad(:,:,:,5,:) = grad(:,:,:,5,:)./2;


    %(2,3) entry
    grad(:,:,2:end,6,:) = ( v(:,:,2:end,2,:) - v(:,:,1:end-1,2,:) )./dz; %Get dz of v2
    
    grad(:,2:end,:,6,:) = ( v(:,2:end,:,3,:) - v(:,1:end-1,:,3,:) )./dy + ... %Add dy of v3
    grad(:,2:end,:,6,:);
    
    grad(:,:,:,6,:) = grad(:,:,:,6,:)./2;












