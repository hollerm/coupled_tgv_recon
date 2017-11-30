

%3-dimensional test----------------------------------------
n=10;
m=150;
k=12;
nch = 5;


dx = 0.1;
dy = 2;
dz = 1.3;


%Vector operators
x = rand(n,m,k,nch);
y = rand(n,m,k,3,nch);

gx = fgrad_3(x,dx,dy,dz);
divy = bdiv_3(y,dx,dy,dz);

s1 = gx(:)'*y(:);
s2 = -x(:)'*divy(:);

display('############## Vector operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);




%Matrix operators
x = rand(n,m,k,3,nch);
y = rand(n,m,k,6,nch);



gx = sym_bgrad_3(x,dx,dy,dz);
divy = fdiv_3(y,dx,dy,dz);



tmp = gx(:,:,:,1:3,:).*y(:,:,:,1:3,:) + 2*gx(:,:,:,4:6,:).*y(:,:,:,4:6,:);

s1 = sum(tmp(:));
s2 = -x(:)'*divy(:);

display('############## Matrix operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);



%2-dimensional test----------------------------------------
n=13;
m=120;
k=1;
nch = 3;


dx = 0.1;
dy = 2;
dz = 1.3;


%Vector operators
x = rand(n,m,k,nch);
y = rand(n,m,k,3,nch);

gx = fgrad_3(x,dx,dy,dz);
divy = bdiv_3(y,dx,dy,dz);

gx = squeeze(gx(:,:,1,1:2,:));
sy = squeeze(y(:,:,1,1:2,:));
s1 = gx(:)'*sy(:);
s2 = -x(:)'*divy(:);

display('############## Vector operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);




%Matrix operators
x = rand(n,m,k,3,nch);
y = rand(n,m,k,6,nch);



gx = sym_bgrad_3(x,dx,dy,dz);
divy = fdiv_3(y,dx,dy,dz);



tmp = sum(gx(:,:,:,1:2,:).*y(:,:,:,1:2,:),4) + 2*gx(:,:,:,4,:).*y(:,:,:,4,:);

s1 = sum(tmp(:));
s2 = -x(:)'*divy(:);

display('############## Matrix operators ##############')
display(['S1:      ',num2str(s1)]);
display(['S2:      ',num2str(s2)]);
display(['Dif:     ',num2str(abs(s1-s2))]);
display(['Rel-Dif: ',num2str(abs(s1-s2)/(n*m*k*nch))]);

