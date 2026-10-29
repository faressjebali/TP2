function dhat = wiener_ola(x, P, L, filtre)
% Filtre de Wiener par trames (recouvrement 50 %, overlap-add)
% filtre : fonction (rx, idx) -> w, filtre a appliquer sur la trame d'indices idx
x = x(:);
N   = length(x);
H   = L/2;
win = hann(L, 'periodic');
dhat = zeros(N,1);

for s = 1:H:N-L+1
    idx = s : s+L-1;
    rx = autocorr_biased(x(idx), P);
    w  = filtre(rx, idx);            % filtre de la trame (r_d connu ou estime)

    deb = max(1, s-P+1);
    y = filter(w, 1, x(deb:idx(end)));
    dhat(idx) = dhat(idx) + win .* y(end-L+1:end);
end
end
