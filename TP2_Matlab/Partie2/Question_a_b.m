clear; clc; close all;
N     = 1e4;          
n     = (1:N)';
w0    = 0.05*pi;
phi   = pi/4;         
P     = 10;           
sigv2 = 0.1;          

d = sin(n*w0 + phi);
v = sqrt(sigv2)*randn(N,1);
x = d + v;

rx = xcorr(x, P-1, 'biased');  rx = rx(P:end);
rd = xcorr(d, P-1, 'biased');  rd = rd(P:end);
Rx = toeplitz(rx);            
wopt = Rx \ rd;              
dhat = filter(wopt, 1, x);

%% Question (a)
idx = 101:200;
figure;
plot(idx, d(idx), 'b', idx, x(idx), 'Color', [0.7 0.7 0.7]); hold on;
plot(idx, dhat(idx), 'r', 'LineWidth', 1.5);
grid on; xlabel('n'); ylabel('Amplitude');
legend('d(n) : signal propre', 'x(n) : signal bruité', '$\hat{d}(n)$ : signal débruité', 'Interpreter', 'latex');
title('Débruitage par filtre de Wiener (P = 10, \sigma_v^2 = 0.1)');

%% Question (b)
Rd0 = rd(1); 
EQM_min = Rd0 - wopt' * rd;
fprintf('Valeur de l''EQM_min théorique (pour P = %d, sigma_v^2 = %.2f) : %.4f\n', P, sigv2, EQM_min);



