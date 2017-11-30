

function [nrm] = power_it(K,Kt,dims,n)

x = rand(dims);



for i=1:n
    x = Kt(K(x));
    nrm = sqrt(sum(abs(x(:)).^2));
    x = x./nrm;
    if rem(i,10)==0
        display(['At it: ',num2str(i)]);
    end
end

%Nrm converges to |K|^2, hence:

nrm = sqrt(nrm);

