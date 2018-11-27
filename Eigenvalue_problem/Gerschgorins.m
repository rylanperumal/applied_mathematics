function [center, radius] = Gerschgorins(A)
%
% This function computes the center and radius using Gerschgorin's Theorem
% OUTPUT -> Returns center vector and radius vector
% INPUT -> Matrix A which is nxn
%
rad = zeros(length(A), 1);
cen = zeros(length(A), 1);
for i = 1:length(A)
 for j = 1:length(A)
    if(i == j)
        cen(i) = A(i, j);
    else
        rad(i) = rad(i) + abs(A(i, j));
    end
 end
end
center = cen;
radius = rad;
end