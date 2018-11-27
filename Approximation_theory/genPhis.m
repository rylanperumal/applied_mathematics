function phin = genPhis(n, x)
    % This fucntion takes 'n' and a variable x. Returns a function phi_n(x)
    % makes use of function genTns(n, x), makes the set orthonormal
    if n == 0
        phin = 1/sqrt(pi) * ones(1, length(x));
        return
    end
    phin = sqrt(2 / pi) * genTns(n, x);
end