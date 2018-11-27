function Tn = genTns(n, x)
    % This codes computes the nth order ChebyShev ploynomial denoted in the
    % notes by T_n(x), always work in vectors
    if n == 0
        % we must return a vector
        Tn = ones(1, length(x)); % row vector
        return
    end
    if n == 1
        % we want to return x
        Tn = x;
        return
    end
    Tn = ((2 * x.* genTns(n - 1, x)) - genTns(n-2, x)); % recursive call
end