% Define os autovalores desejados
lambda1 = 3 + 4i;
lambda2 = 3 - 4i;
lambda3 = 5 + 6i;
lambda4 = 5 - 6i;

% Cria matrizes de rotação complexas para os autovalores
A1 = [real(lambda1), -imag(lambda1); imag(lambda1), real(lambda1)];
A2 = [real(lambda3), -imag(lambda3); imag(lambda3), real(lambda3)];

% Combina as matrizes 2x2 em uma matriz 4x4
A = blkdiag(A1, A2);

% Ajuste o quadrante inferior esquerdo para ter números iguais
equal_number = 1;
A(3:4, 1:2) = equal_number;

% Exibe a matriz resultante
disp('Matriz A com dois pares de autovalores complexos:');
disp(A);

% Verifica os autovalores para confirmar
eigenvalues = eig(A);
disp('Autovalores de A:');
disp(eigenvalues);
