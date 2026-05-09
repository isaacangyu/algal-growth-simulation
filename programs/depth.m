function depth()
    z       = 8.0;
    Iin     = 120.0;
    HL = 120.0;
    Kbg     = 0.35;
    kP      = 8.0;
    nQuad   = 150; % discrete quadrature
    s = linspace(0, z, nQuad);

    P = 5e-4;
    L = @(s,P) Iin .* exp((-Kbg - kP .* P) .* s);
    figure; 
    plot(s, L(s,P), 'DisplayName', 'L');
    integrand = @(s,P) L(s,P) ./ (L(s,P) + HL);
    g = @(P) trapz(s, integrand(s,P)) / z;
    fprintf('g with fixed P %d\n', g(P));

    numPs = 10;
    gPs = zeros(numPs);
    Ps = linspace(1e-4, 1e-3, 10);
    for i = 1:numPs
        gPs(i) = trapz(s, integrand(s,Ps(i))) / z;
    end
    figure;
    plot(Ps, gPs, 'DisplayName', 'g');

    
    % hold on;
    % title('Light intensity vs. Areal/Volumetric growth');
    % legend();
    % hold off;
end