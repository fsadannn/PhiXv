close all
clear all

n = 500;
file='Brusselator.txt';
debug = 0;
fid = fopen(file, 'w');
f=@f_brusselator;
J=@J_brusselator;

N=n;
x0=[1+sin((2*pi/(N+1))*(1:N)),repmat(3,1,N)];
IT=[0 1];
tt=f(0,x0(:));
tt=J(0,x0(:));
clear tt;
exp_4 = exp4set('Complex','off',...
    'AbsTol',1.0e-12,'RelTol',1.0e-9,'MatrixFunctions','arnoldi',...
    'CheckDimensions','off',...
    'ClearInternalData','on','OutputFcn',@outf);
alphac = 0.02;
aa=strsplit(file,'.');
fprintf(fid,'Problem: %s\n',aa{1});
fprintf('Problem: %s\n',aa{1});
fprintf(fid,'Param: Ne@%d\n',2*N);
fprintf(fid,'Param: IT@%.16f@%.16f\n',IT(1),IT(2));
fprintf(fid,'Param: \\alpha@%.2f\n',alphac);
fprintf(fid,'Precision: refine\n');
fprintf(fid,'Integrator: Exp4\n');
tic
[t, y] = exp4(f, J,IT, x0,exp_4);
toct=toc;
fprintf(fid,'t = %f\n',toct);
save('br.mat')
inf = outf(0,0,'get');
fprintf(fid,'Steps = %d\n',inf.stats.integratorSteps);
fprintf(fid,'Failed = %d\n',inf.stats.rejectedSteps);
fprintf(fid,'f-Evals = %d\n',inf.stats.RHSEv);
fprintf(fid,'J-Evals = %d\n',inf.stats.JacEv);
fprintf(fid,'MD = %d\n',inf.stats.matFun.NofDiag);
fprintf(fid,'EL = %d\n',inf.stats.matFun.NofDiag);
fprintf(fid,'K-Dim-Total = %d\n',inf.stats.matFun.NofKrySteps);
maxdim = 1;
mindim = Inf;
for field = fieldnames(inf.stats.matFun)'
    fieldName = field{1};
    if isfield(inf.stats.matFun.(fieldName), 'data')
        dimm = inf.stats.matFun.(fieldName).data(2, :);
        index=find(dimm>0);
        if ~isempty(index)
            maxdim = max(maxdim,max(dimm(index)));
            mindim = min(mindim,min(dimm(index)));
        end
    end
end
fprintf(fid,'K-Dim-Min = %d\n',mindim);
fprintf(fid,'K-Dim-Max = %d\n',maxdim);
outf(0,0,'clean');
[T,Y] = exactsol(f,J,t,x0);
Y = real(Y);
rel=RelError(Y,y);
fprintf(fid,'RelativeError = %.16f\n',rel);


