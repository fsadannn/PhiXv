%  phi1LLDP computes an approximation of phi(A)*u, phi(1/5*A)*u
%  phi(3/10*A)*u, phi(4/5*A)*u, phi(8/9*A)*u
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
function phi = k_exp(A, u, h, m)

n=length(u);

btol  = 2*eps;
mb    = m;
beta = norm(u);

% begin Arnoldi
V = zeros(n,m+1);
H = zeros(m+1,m+1);
V(:,1) = (1/beta).*u;

if any(isnan(V(:,1)))
   phi = sparse(length(u),1);
   return;
end

for j = 1:m
    p = h.*(A*V(:,j));
    s = V(:,1:j);
    H(1:j,j) = s.'*p;
    p = p - s*H(1:j,j);
    s = norm(p);
    if s < btol
        mb = j;
        break;
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s).*p;
end

M1 = expm(H(1:mb,1:mb));

phi = V(:,1:mb)*(beta.*M1(1:mb,1));

end