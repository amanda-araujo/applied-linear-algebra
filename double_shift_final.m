function [B, erros] = double_shift_final(H, nmax, tol)
% Algoritmo QR com deslocamento duplo (double shift) para matriz Hessenberg
% 
% Input:
% H: matriz n x n Hessenberg
% nmax: nº de iterações
% tol: tolerância
%
% Output:
% B: matriz n x n (aproximadamente) triangular superior >>> autovalores na 
%    diagonal
% erros: a norma da parte subdiagonal da matriz, erro de cada iteração 

    n = size(H, 1); % tamanho da matriz H

    if nargin < 2
        nmax = 1000; % Define o número máximo de iterações, se não especificado
    end
    if nargin < 3
        tol = 1e-10; % Define a tolerância para convergência, se não especificada
    end
    
    % Inicialização 
    iter = 0;
    erros = []; % vetor de erros

    % Loop para percorrer a matriz de baixo para cima
    for k = n:-1:2
        
        I = eye(k); % matriz identidade de tamanho k

        % Loop enquanto o elemento subdiagonal for maior que a tolerância
        while abs(H(k,k-1)) >= tol * (abs(H(k,k)) + abs(H(k-1,k-1)))
            % Incrementa o contador de iterações
            iter = iter + 1;
            
            % Verifica se o número máximo de iterações foi atingido
            if iter > nmax
                warning('Número máximo de iterações atingido');
                break;
            end

            % Parte do deslocamento duplo
            if k > 2 && abs(H(k-1,k-2)) <= tol * (abs(H(k-1,k-1)) + abs(H(k-2,k-2)))
                
                % Calcula os autovalores da submatriz 2x2
                %lambda = eig(H(k-1:k,k-1:k)); >> raízes pol grau 2
                a = 1;
                b = -(H(k-1,k-1) + H(k,k));
                c = H(k-1,k-1)*H(k,k) - H(k-1,k)*H(k,k-1);
                delta = b^2 - 4*a*c;
                lambda1 = (-b + sqrt(delta)) / (2*a);
                lambda2 = (-b - sqrt(delta)) / (2*a);
                lambda = [lambda1; lambda2];
                
                % Aplica o deslocamento duplo para cada autovalor
                for i = 1:2
                    mu = lambda(i);
                    [Q, R] = qr(H(1:k,1:k) - mu * I);
                    H(1:k,1:k) = R * Q + mu * I;
                end

            % Parte do deslocamento simples
            else
                % Usa o elemento diagonal como deslocamento
                mu = H(k,k);
                [Q, R] = qr(H(1:k,1:k) - mu * I);
                H(1:k,1:k) = R * Q + mu * I;
            end

            % Calcula e armazena o erro
            error = norm(tril(H,-1));
            erros = [erros; error];
        end

        % Define o elemento subdiagonal como zero
        H(k,k-1) = 0;
        % Atualiza a matriz B com a matriz H modificada
        B = H;
    end
end
