function [c, s] = givens(a, b)
% Compute Givens rotation parameters c and s

if b == 0
    c = 1;
    s = 0;

else
    if abs(b) > abs(a)
        tau = -a / b;
        s = 1 / sqrt(1 + tau^2);
        c = s * tau;
    else
        tau = -b / a;
        c = 1 / sqrt(1 + tau^2);
        s = c * tau;
    end

end

end