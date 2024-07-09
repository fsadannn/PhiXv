%  phi1_single computes an approximation of h*phi(h*A)*u,
%  for a general matrix A using Krylov subspace projection techniques.
%  Here, phi(z) = (exp(z)-I)/z and this phis are use in LLDP scheme.
%  It does not compute the matrix functions in isolation but instead,
%  it computes directly the action of these functions on the
%  operand vectors. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%  It only compute one matrix exponential and use the relation betwin
%  the coefficients and the exponential properties to compute the
%  others phi.
function phi = phi1_single(A, u, h, m)

n=length(u);

btol  = 2*eps;
% mb    = m;
beta = norm(u);
% k1 = 3;
% begin Arnoldi
V = zeros(n,m+1);
H = zeros(m+2,m+2);
V(:,1) = (1/beta).*u;

for j = 1:m
    p = h.*(A*V(:,j));
    s = V(:,1:j);
    H(1:j,j) = s.'*p;
    p = p - s*H(1:j,j);
    s = norm(p);
    if s < btol
%         k1 = 0;
%         mb = j;
        break;
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s).*p;
end


hk = H(m+1,m);
H(m+1,m) = 0;

H(1,m+1) = h;
H(m+1,m+2) = h;

% nhC=norm(H,'inf');
% scaling calculation
% [~,e] = log2(nhC);
% s = max(0,e+1);
% exponential calculation
% M1 = expm64v4(H,6,s);
% M1 = expmc_v4(H,6,s);
M1 = expm(H);

% calating \hat{E}
M1(m+1,m+1) = hk*M1(m,m+2);

mx = m + 1;
% Matix projection

phi = V*(beta*M1(1:mx,mx));

end