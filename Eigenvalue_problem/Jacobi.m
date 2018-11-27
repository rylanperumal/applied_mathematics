function [out, theta] = Jacobi(A)
%
% INPUT -> Takes in a symmetric matrix A
% OUTPUT -> Reduced Matrix using the Jacobi method
% Can add tolerence instead of for loop to stop iterations
% Tolerence sum of off diagonal elements < 10^(-3)
%
val = sum(sum(A - diag(diag(A))));
while (val > (10 ^ -4) || val < - (10 ^ -4))
    B = A - diag(diag(A));
    maxVal = max(abs(B(:)));
    P = eye(length(A));
    if maxVal == 0
        out = A;
        return;
    end
    [row, col] = find(abs(B) == maxVal);
    r = row(end);
    c = col(end);
    
    if(A(c, c) == A(r, r))
        if(A(r, c) > 0)
            theta = (1 / 2) * (pi / 2);
        else
           theta = - (1 / 2) * (pi / 2); 
        end
    else
        theta = (1 / 2) * (atan((2 * A(r, c)) / (A(c, c) - A(r,r))));
    end
%     theta = 1 / 2 * (atan((2 * A(r, c)) / (A(c, c) - A(r, r))));
    
    P(r, r) = cos(theta);
    P(r, c) = sin(theta);
    P(c, r) = -sin(theta);
    P(c, c) = cos(theta);
    
    A = transpose(P) * A * P;
    val = sum(sum(A - diag(diag(A))));
end
out = A;
end
