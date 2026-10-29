function C = kmeans_simple(X, K)
% k-means sur les colonnes de X : retourne les K centres (colonnes de C)
rng(0);
C = X(:, randperm(size(X,2), K));
for it = 1:100
    D = sum(X.^2, 1)' + sum(C.^2, 1) - 2*(X'*C);   % distances au carre
    [~, lab] = min(D, [], 2);
    Cn = C;
    for k = 1:K
        if any(lab == k), Cn(:,k) = mean(X(:, lab == k), 2); end
    end
    if isequal(Cn, C), break; end
    C = Cn;
end
end
