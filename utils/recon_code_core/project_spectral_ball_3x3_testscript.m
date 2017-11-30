

display('################Sanity Test################');


A = cell(1,1);
A{1} = randi(512,3,5)+ i*randi(512,3,5) -255;
A{2} = [0 0 0; 5 10 0.5; -2 -4 1];
A{3} = [0 0 0; 0 0 0; -2 -4 1];
A{4} = [5000 0 0; 0 2000 0; 0 0 1];
A{5} = [0 5000 0; 0 2000 0; 1 1 1];
A{6} = [100 0 0; 0 0 0; 0 0 0];



for jj=1:length(A);

    z = zeros(2,3,2,size(A{jj},1),size(A{jj},2));
    for ii=1:2;for ll=1:3;for kk=1:2
        z(ii,ll,kk,:,:) = A{jj};
    end;end;end
    
    disp('################')
    disp('Testing A*A^t=: ')
    disp(num2str(A{jj}*A{jj}'))
    disp(' ')

    Asv = svd(A{jj});
    alpha = 0.9*max(Asv);
    [z] = project_spectral_ball_3x3(z,alpha);
    Bsv = svd(squeeze(z(1,1,1,:,:)));
    disp(' ')
    disp(['Orig:  ',num2str(Asv')]);
    disp(['Tresh: ',num2str(alpha)]);
    disp(['Prox:  ',num2str(Bsv')]);
    Asv(Asv>alpha) = alpha;
    dif = sum(abs(Asv-Bsv));
    if dif < 10^(-10)
        disp(['-----Passed with error ',num2str(dif),'-----']);
    else
        warning(['Incorrect result: Error: ',num2str(dif)]);
    end
    disp(' ');
    disp(' ');
    disp(' ');
    
end



display('################Speed Test#################');


sn=200;
sm=200;
sz=200;
noc=5;

disp(' ');
disp('-----Costum implementation-----')
z = rand(sn,sm,sz,3,noc) + i*rand(sn,sm,sz,3,noc);
alpha = 0.8;
tic
z = project_spectral_ball_3x3(z,alpha);
toc

if 1
disp(' ');
disp('-----2xn prox implementation-----')
z1 = rand(sn,sm,sz,noc) + i*rand(sn,sm,sz,noc);
z2 = rand(sn,sm,sz,noc) + i*rand(sn,sm,sz,noc);
alpha = 0.8;
tic
[z1,z2] = project_spectral_ball_2x2(z1,z2,alpha);
toc
end

if 0
disp(' ');
disp('-----Naive loop for Matlab SVD-----')
z = rand(sn,sm,sz,3,noc) + i*rand(sn,sm,sz,3,noc);
tic
for ii=1:sn;for ll=1:sm;for kk=1:sz

    [U,D,V] = svd( squeeze(z(ii,ll,kk,:,:)),'econ');
    S = diag([1 1 1]);
    S(D>alpha) = alpha./D(D>alpha);
    z(ii,ll,kk,:,:) = U*S*U'*squeeze(z(ii,ll,kk,:,:));
    
end;end;end
toc
end


