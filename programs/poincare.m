function poincare()
    % Parameters
    a = 7.5e-3; b = 1; c = 1; d = 0.1; %#ok<NASGU>
    T = 1.0;  % 1 year

    % Seasonal forcing: ensure it stays > 0
    Sbar = 0.6;  A = 0.4;   % so S(t) in [0.2, 1.0]
    S = @(t) Sbar + A*sin(2*pi*t);

    % Vector field
    f = @(t,x) [ a - S(t)*x(1)*x(2);
                 S(t)*x(1)*x(2) - d*x(2) ];

    % Poincare map: integrate one period
    Pmap = @(x0) one_period_map(f, x0, T);

    % --- Find a fixed point of the Poincare map (a periodic orbit) ---
    xguess = [0.2; 0.05];    % initial guess, adjust if needed
    opts_fsolve = optimoptions('fsolve','Display','iter','FunctionTolerance',1e-12);

    xstar = fsolve(@(x) Pmap(x) - x, xguess, opts_fsolve);
    fprintf('Fixed point of Poincare map x* = [%.6g, %.6g]^T\n', xstar(1), xstar(2));

    % --- Estimate Jacobian of Poincare map at xstar (finite differences) ---
    J = poincare_jacobian(Pmap, xstar);
    eigJ = eig(J);
    fprintf('Floquet multipliers (approx) = %.6g%+.6gi, %.6g%+.6gi\n', ...
        real(eigJ(1)), imag(eigJ(1)), real(eigJ(2)), imag(eigJ(2)));
    fprintf('Magnitudes = %.6g, %.6g\n', abs(eigJ(1)), abs(eigJ(2)));

    % --- Show convergence by iterating the map ---
    x = [0.35; 0.08];
    K = 365;
    X = zeros(2,K+1); X(:,1)=x;
    for k=1:K
        x = Pmap(x);
        X(:,k+1)=x;
    end

    figure;
    plot(0:K, X(1,:), '-o'); hold on;
    plot(0:K, X(2,:), '-o');
    xlabel('Year (Poincare iterate k)');
    ylabel('State at t = kT');
    legend('N(kT)','P(kT)','Location','best');
    grid on;
end

function xT = one_period_map(f, x0, T)
    % Integrate from t=0 to t=T; return x(T)
    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);
    sol = ode45(f, [0 T], x0, opts);
    xT = deval(sol, T);
end

function J = poincare_jacobian(Pmap, x)
    % Finite-difference Jacobian of the Poincare map at x
    n = numel(x);
    J = zeros(n);
    fx = Pmap(x);
    h = 1e-7;  % step size

    for j=1:n
        ej = zeros(n,1); ej(j)=1;
        fplus  = Pmap(x + h*ej);
        fminus = Pmap(x - h*ej);
        J(:,j) = (fplus - fminus) / (2*h);
    end
end