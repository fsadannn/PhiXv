function [phi, err, dim, nexpo] = k_phi_fj_adapt2(func,A,h,t,y,m,rtol,atol,kdmax,kdmin)
%PHI1LLDP_hJ_f Summary of this function goes here
%   Detailed explanation goes here
fy = func(t,y);
n=sum(size(A));
m=min(m,n);
gap = length(fy);
gap = n-gap;
b = [A(:,end);zeros(gap,1)];
b(end-1)=1;

% the coeficient is in rfacs( m )
% but is expensive call a function
% and coeficient are hardcode
rfacmin=1;
fac=1/log(2);
gamma = 0.001;
% gamma = 0.005;

nexpo=0;
% 
% tolr=[1.0e-6 1.0e-9 1.0e-12 1.0e-15];
% Table=[2 3 4 4;3 4 5 6];

btol = 2*eps;
mb = m;
k1 = 3;
ww = atol+rtol.*abs(b);

normy = norm(y);

% hsubspace = h;
delta = sqrt((1+normy)*eps);
idelta = 1/(2*delta);

normb = norm(b);
V = zeros(n,m+1);
H = zeros(m+3,m+3);

% A=hsubspace.*A;
rndoff= eps;



% build base of hA and b

% begin Arnoldi
V(:,1)=(1/normb).*b;
% test if norm(b) is 0
if any(isnan(V(:,1)))
    phi = sparse(n,1);
    err = 0;
    warning('llint:phi1LLDP_hJ_f_na',...
       'First vector of krylov subspace seems to be 0.')
    return;
end

for j = 1:m
%     p = A*V(:,j);
%     p = [idelta.*(func(t,y+delta.*V(1:end-gap,j))-fy)+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
    p = [idelta.*(func(t,y+delta.*V(1:end-gap,j))-func(t,y-delta.*V(1:end-gap,j)))+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
    s = V(:,1:j);
    H(1:j,j) = s.'*p;
    p = p - s*H(1:j,j);
    s = norm(p);
    if s < btol
        k1 = 0;
        mb = j;
        break
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s).*p;
end
% end Arnoldi
% using scaling invariance property of Arnoldi
% rescaling H to have H for the subspace of A and b
% H=(1/hsubspace).*H;
% build \hat{H}
if k1 == 0
    if mb>1
        mb=mb-1;
        warning('llint:phi1LLDP_hJ_f_na',...
       'Breakdown at dimension 1.')
    end
    m=mb;
    hk = H(m+1,m);
    H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
    % V=V(:,1:m+1);
else
    hk = H(m+1,m);
    H(m+1,m) = 0;
end
H(1,m+1) = 1;
H(m+1,m+2) = 1; H(m+2,m+3) = 1;
% avm1dot=A*V(:,m+1);
% avm1dot = [idelta.*(func(t,y+delta.*V(1:end-gap,m+1))-fy)+A*V(end-gap+1:end,m+1);V(end-gap+2:end,m+1);0];
avm1dot = [idelta.*(func(t,y+delta.*V(1:end-gap,m+1))-func(t,y-delta.*V(1:end-gap,m+1)))+A*V(end-gap+1:end,m+1);V(end-gap+2:end,m+1);0];

work=1;

while work
    % select p-p of Pade
%     nhC = h*norm(H,'inf');
%     col = find(tolr>=rtol,1,'last');
%     fil = (nhC>=1)+1;
%     pd = Table(fil,col);
% 
%     % scaling calculation
%     [~,e] = log2(nhC);
%     s = max(0,e+1);
%     pd=6;
    
    % exponential calculation
%     M1 = expm64v4(h.*H,pd,s);
    M1 = expm(h.*H);
    nexpo = nexpo + 1;
    
    % calating \hat{E}
    beta = normb;
    %error relative
    % the divsion by h is because Av_{m+1} is in reality
    % hAv{m+1}(avm1dot) and need rescaling
    err_rel=sqrt((1/n)*sum((((hk*M1(m,m+3)*beta).*avm1dot)./ww).^2));
    if err_rel/gamma>=1 && m<kdmax && k1~=0
        % the coeficient is in rfacs( m )
        % but is expensive call a function
        % and coeficient are hardcode
        %  [ rfacmin,rfacmax ] = rfacs( m );
        rfacmax=max(1,m/3);
        knew =  log(err_rel/gamma)*fac;
        knew = ceil(m + min(rfacmax,max(knew,rfacmin)));
        knew = max(kdmin, min(kdmax,knew));

        H = [H(1:m,1:m),zeros(m,knew-m+3);zeros(knew-m+3,knew+3)];
        H(m+1,m) = hk;
%         if size(V,2)<knew+1
            V = [V,zeros(n,knew-size(V,2)+1)];
%         end
        mtemp=m+1;
        m=knew;
        mb=m;
        k1 = 3;
        
        j=mtemp;
        s = V(:,1:j);
        H(1:j,j) = s.'*avm1dot;
        avm1dot = avm1dot - s*H(1:j,j);
        s = norm(avm1dot);
        if s < btol
            k1 = 0;
            mb = j;
            mtemp=m;
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s)*avm1dot;

        for j = mtemp+1:m
%             p = A*V(:,j);
%             p =  [idelta.*(func(t,y+delta.*V(1:end-gap,j))-fy)+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
            p =  [idelta.*(func(t,y+delta.*V(1:end-gap,j))-func(t,y-delta.*V(1:end-gap,j)))+A*V(end-gap+1:end,j);V(end-gap+2:end,j);0];
            s = V(:,1:j);
            H(1:j,j) = s.'*p;
            p = p - s*H(1:j,j);
            s = norm(p);
            if s < btol
                k1 = 0;
                mb = j;
                break;
            end
            H(j+1,j) = s;
            V(:,j+1) = (1/s)*p;
        end
%         H=(1/hsubspace).*H;
        % build \hat{H}
        if k1 == 0
            if mb>1
                mb=mb-1;
            end
            m=mb;
            hk = H(m+1,m);
            H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
        else
            hk = H(m+1,m);
            H(m+1,m) = 0;
        end
        H(1,m+1) = 1;
        H(m+1,m+2) = 1; H(m+2,m+3) = 1;
%         avm1dot=A*V(:,m+1);
%         avm1dot = [idelta.*(func(t,y+delta.*V(1:end-gap,m+1))-fy)+A*V(end-gap+1:end,m+1);V(end-gap+2:end,m+1);0];
        avm1dot = [idelta.*(func(t,y+delta.*V(1:end-gap,m+1))-func(t,y-delta.*V(1:end-gap,m+1)))+A*V(end-gap+1:end,m+1);V(end-gap+2:end,m+1);0];
    else
        break;
    end
end

M1(m+1,m+1) = hk*M1(m,m+2);
M1(m+2,m+1) = hk*M1(m,m+3);

mx = m + 1;

avm1dot = V(:,1:mx);
% Matix projection
phi = avm1dot*(beta*M1(1:mx,mx));

err = max(err_rel,rndoff);
dim = m;
end

