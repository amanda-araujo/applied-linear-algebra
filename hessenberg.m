function [Q, H] = hessenberg(A)
% Find the Hessenberg decomposition of A=Q*H*Q'
% [Q H]=hessenberg(A)
%  
% Input:
% A: a n by n matrix
% 
% Output:
% Q: a n by n orthogonal matrix
% H: a n by n upper Hessenberg matrix

n = size(A);
n = n(1);

% Initialization
Q = eye(n); % identity
H = A;      % matrix A

% Loop over columns of h from 1 to n - 2
for k = 1 : n-2

  % Extraction of the subvector x: (n - k) x 1 from H for the current column k
  x = H(k+1 : n, k);

  % Initializing vector u
  u = zeros(n-k, 1);

  % Computation of the first element of u (norm)
  u(1) = -(sign(x(1))*(x(1)~=0) + (x(1)==0)) * (x'*x)^.5; 
  
  % Computation of the Householder vector v
  v = (x - u);
  v = v/(v'*v)^.5; % normalization

  % Construction of the Householder matrix W
  W=[eye(k) zeros(k,n-k);zeros(n-k,k) eye(n-k)-2*v*v'];

  % Application of the Householder transformation
  H = W*H*W';
  H

  % Update of the orthogonal matrix Q
  Q = Q*W';

end 

% If only one argument is given, return H
if nargout()<=1
  Q = H;
end

end