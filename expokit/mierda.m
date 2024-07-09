dim=3;


clear all 
A=randn(6);
v=ones(6,1);
[w,err] = expv( 0.5, A, v,1.0e-7,dim);
[expm(0.5*A)*v w]



% clear all 
% A=randn(6);
% v=ones(6,1);
% [w,err,mierda,VM,bt,expH,Tstep] = expvmio( 0.5, A, v,1.0e-7,dim);
% [w VM*bt*expH(:,1)]
% 
% [w1,err,mierda,VM1,bt1,expH1,Tstep1] = expvmio( 1, A, v,1.0e-7,dim);
% expH2=expH*expH;
% [w1 VM*bt*expH2(:,1)]
% 
% [expm(2*Tstep*A)*v VM*bt*expH2(:,1)]


% clear all 
% A=randn(6);
% v=ones(6,1);
% [w,err,mierda,VM,bt,expH,Tstep] = expvmio( 0.25, A, v,1.0e-7,dim);
% [w VM*bt*expH(:,1)]
% 
% [w1,err,mierda,VM1,bt1,expH1,Tstep1] = expvmio( 0.5, A, v,1.0e-7,dim);
% expH2=expH*expH;
% [w1 VM*bt*expH2(:,1)]

% clear all 
% A=randn(6);
% v=ones(6,1);
% [w1,err,mierda,VM,bt,expH] = expvmio( 0.25, A, v,1.0e-7,dim);
% [w1 VM*bt*expH(:,1)]
% [w2,err] = expv( 0.5, A, v,1.0e-7,dim);
% expH2=expH*expH;
% [w2 VM*bt*expH2(:,1)]