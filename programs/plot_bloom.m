function [tt, N, P] = plot_bloom(I, q, N_0, P_0, step)
    days = 365;
    timesteps = days/step;
    tt = linspace(0, days, timesteps);
    f_N = @(n, p) I - n .* p + q;
    f_P = @(n, p) n .* p - p;
    
    N = N_0; 
    P = P_0;
    for t = 1 : timesteps
        N(t+1) = max(N(t) + f_N(N(t), P(t)), 0);
        P(t+1) = max(P(t) + f_P(N(t), P(t)),0);
    end

    hold on;
    plot(tt, N, 'b', tt, P, 'r');
    legend('Nutrient (N)', 'Phytoplankton (P)');
    xlabel('Time');
    ylabel('Population');
    hold off;
end

% Huppert Lotka-Volterra model
N_0 = 5e-4;
P_0 = 5e-2;
q = 0;
I = 7.5e-2;
% Heggerud light intensity model
z = 7;
H = 120;
K_bg = 0.3;
k = 4e-4;
% seasonality effect on light intensity
S_I = 

mu = growth(N, P, z, H, S_I, K_bg, k, 'All');

f_N = @(N) I - mu - N * q;
f_P = @(N, P) mu - P;