function [Q, R] = hess_qr_givens(A)
% Hessenberg matrix QR decomposition via Givens rotations

n = size(A, 1);

% Initializing
Q = eye(n);
R = A;

for j = 1 : n - 1

    % Compute Givens rotation parameters
    [c, s] = givens(R(j, j), R(j + 1, j));

    % Givens rotation matrix of two rows of A
    G = [c s; -s c];

    % Apply Givens rotation to eliminate R(i, k)
    R(j : j+1, j : n) = G' * R(j : j+1, j : n);

    % Atualization of orthogonal matrix Q = G1 * G2 ... * Gt
    % Apply Givens rotation to accumulate Q
    Q(:, j : j+1) = Q(:, j : j+1) * G;
    
    % Display
    Q
    R

end

end