function R = vecteurs_rd(a, P, L)
% Vecteurs r_d (retards 0..P-1) des trames actives de a (pas L/2)
seuil = 0.01 * mean(a.^2);
R = [];
for s = 1:L/2:length(a)-L+1
    t = a(s:s+L-1);
    if mean(t.^2) > seuil            % trames de silence ignorees
        R(:, end+1) = autocorr_biased(t, P); %#ok<AGROW>
    end
end
end
