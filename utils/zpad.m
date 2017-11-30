function res = zpad(u,nx,ny,nz,nt)
    %  res = zpad(u,nx,ny)
    %  Zero pads a 2D matrix around its center.
    %
    %  res = zpad(u,nx,ny,nz,nt)
    %  Zero pads a 4D matrix around its center

    if nargin < 2
        error('must have a target size')
    end

    if nargin == 2
        s = nx;
    end

    if nargin == 3
        s = [nx,ny];
    end

    if nargin == 4
        s = [nx,ny,nz];
    end

    if nargin == 5
        s = [nx,ny,nz,nt];
    end

    m = size(u);
    if length(m) < length(s)
	    m = [m, ones(1,length(s)-length(m))];
    end
	
    if sum(m==s)==length(m)
	res = u;
	return;
    end

    res = zeros(s);
    
    for n=1:length(s)
	    idx{n} = floor(s(n)/2)+1+ceil(-m(n)/2) : floor(s(n)/2)+ceil(m(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s)=u;',cmd);
    eval(cmd);
end