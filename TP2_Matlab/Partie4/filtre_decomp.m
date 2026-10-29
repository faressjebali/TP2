function [w, rd1] = filtre_decomp(rx, V, U, modele)
% Filtre de Wiener avec r_d1 = a v_i, ou r_x ~ a v_i + b u_j (a, b >= 0)
% R_x mesure (modele = false) ou R_x = toeplitz(a v_i + b u_j) (modele = true)
best = inf;  rd1 = zeros(size(rx));  rxm = rd1;
for i = 1:size(V,2)
    for j = 1:size(U,2)
        A = [V(:,i) U(:,j)];
        c = A \ rx;
        if any(c < 0)                % contrainte a, b >= 0 : un seul representant
            a = max(0, V(:,i)'*rx / (V(:,i)'*V(:,i)));
            b = max(0, U(:,j)'*rx / (U(:,j)'*U(:,j)));
            if norm(rx - a*V(:,i)) < norm(rx - b*U(:,j)), c = [a; 0]; else, c = [0; b]; end
        end
        res = norm(rx - A*c);
        if res < best
            best = res;  rd1 = c(1) * V(:,i);  rxm = A * c;
        end
    end
end
if modele, w = toeplitz(rxm) \ rd1; else, w = toeplitz(rx) \ rd1; end
if any(~isfinite(w)), w = zeros(size(rx)); end   % trame vide
end
