function lotka_volterra()
    X0 = [5e-2 5e-4]';
    Xeq = [0.1 / 1 7.5e-3 * 1 / (1 * 0.1)]
    
    function dx = rhs(~, x)
        N = x(1);
        P = x(2);
        I = 7.5e-4 * 1 / 0.1;
        dN = I - 1 * N * P;
        dP = 1 * N * P - 1 * P;
        dx = [dN; dP];
    end

    [t, X] = ode45(@(t,x) rhs(t, x), 0:2000, Xeq.');
    n = X(:,1);
    p = X(:,2);
    hold on;
    plot(t, n / 5, 'DisplayName', 'Nutrient');
    plot(t, p, 'DisplayName', 'Phytoplankton');
    legend();
    hold off;
    % ylim([min(min(N),min(P)) max(max(N),max(P))]);
end