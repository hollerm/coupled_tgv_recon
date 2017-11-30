%Function to calculate the divergence of a function mapping from R^3 to [3,nch], where nch is an arbitrary number of channels
%Input:             Function [n,m,o,3,nch]
%Output:            Div      [n,m,o,1,nch]
%Differences:       Backward Differences [ u(i) = u(i) - u(i-1) ]
%Bdry extension:    Zero boundary extension [ u(n+1) = 0 ]

%NOTE: This function is the negative adjoint to fgrad_3.m

function [div] = bdiv_3(v,dx,dy,dz)

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
div = zeros([size(v,1),size(v,2),size(v,3),1,size(v,5)]); %Artificial singelton dimension for consistency with v


	%dx of v1
	div(1      ,:,:,:,:) =  ( v(1      ,:,:,1,:)                      )./dx;
	div(end    ,:,:,:,:) =  (                    - v(end-1  ,:,:,1,:) )./dx;
	div(2:end-1,:,:,:,:) =  ( v(2:end-1,:,:,1,:) - v(1:end-2,:,:,1,:) )./dx;

	%dy of v2
	div(:,1    ,:,:,:,:) =  ( v(:,1      ,:,2,:)                      )./dy + ...
	div(:,1    ,:,:,:,:);
	div(:,end  ,:,:,:,:) =  (                    - v(:,end-1  ,:,2,:) )./dy + ...
	div(:,end  ,:,:,:,:);
	div(:,2:end-1,:,:,:) =  ( v(:,2:end-1,:,2,:) - v(:,1:end-2,:,2,:) )./dy + ...
	div(:,2:end-1,:,:,:);

	%dz of v3
	if size(v,3)>1
	    div(:,:,1    ,:,:,:) =  ( v(:,:,1      ,3,:)                      )./dz + ...
	    div(:,:,1    ,:,:,:);
	    div(:,:,end  ,:,:,:) =  (                    - v(:,:,end-1  ,3,:) )./dz + ...
	    div(:,:,end  ,:,:,:);
	    div(:,:,2:end-1,:,:) =  ( v(:,:,2:end-1,3,:) - v(:,:,1:end-2,3,:) )./dz + ...
	    div(:,:,2:end-1,:,:);
    end    
    
    %Remove artificial singelton dimension
    %div = squeeze(div);
	
	
	
