

function [test_par_list] = flatten_pars(test_pars);


%Generate Grid for pars
for ii = 1:size(test_pars,1)
    a{ii} = [1:length(test_pars{ii,2})];
end


N = numel(a);
v = cell(N,1);
[v{:}] = ndgrid(a{:});
res = reshape(cat(N+1,v{:}),[],N);

%Generate flattened list
test_par_list = cell(size(res,1),1);
for ii=1:size(res,1)
    
    val = test_pars{1,2}(res(ii,1));
        %Hack: Conversion in case of string parameter
        if iscell(val)
            val = val{1};
        end
   
    test_par_list{ii} = {test_pars{1,1}, val};
    
    for jj=2:size(res,2)
        
        val = test_pars{jj,2}(res(ii,jj));
        %Hack: Conversion in case of string parameter
        if iscell(val)
            val = val{1};
        end
        
        test_par_list{ii} = [ test_par_list{ii} ; { test_pars{jj,1} , val } ];
        
    end
    
end


