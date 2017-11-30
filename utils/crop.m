function res = crop(u,nx,ny,nz,nt)
    %  res = crop(u,nx,ny)
    %  crops a 2D matrix around its center.
    %
    %  res = crop(u,nx,ny,nz,nt)
    %  crops a 4D matrix around its center

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
    if length(s) < length(m)
	    s = [s, ones(1,length(m)-length(s))];
    end
	
    if sum(m==s)==length(m)
        res = u;
        return;
    end
 
    for n=1:length(s)
	    idx{n} = floor(m(n)/2)+1+ceil(-s(n)/2) : floor(m(n)/2)+ceil(s(n)/2);
    end

    % this is a dirty ugly trick
    cmd = 'res = u(idx{1}';
    for n=2:length(s)
    	cmd = sprintf('%s,idx{%d}',cmd,n);
    end
    cmd = sprintf('%s);',cmd);
    eval(cmd);
end