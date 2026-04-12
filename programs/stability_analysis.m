function out = stability_analysis(params, opts, doPrint)
% STABILITY_ANALYSIS  Compute equilibrium and linearized stability data.
%
% out = stability_analysis(params, opts, doPrint)
%
% Uses the effective baseline growth factor implied by SPECIFIC_GROWTH.
% For seasonal terms, the control/baseline value is used so that the
% equilibrium corresponds to the theoretical reference state.

if nargin < 3
    doPrint = true;
end

muEff = specific_growth(0, 0, 0, opts, 'baseline');

if muEff <= 0
    error('Effective growth multiplier must be positive for equilibrium analysis.');
end
if isfield(opts, 'useTemperature') && opts.useTemperature
    Nstar = params.d / (params.c * muEff);
    Pstar = params.a * params.c / (params.b * params.d);
% For e = 0, this matches the equilibrium used in the manuscript.
% if params.e == 0
%     Nstar = params.d / (params.c * muEff);
%     Pstar = params.a * params.c / (params.b * params.d);
% else
%     % General equilibrium when e ~= 0.
%     Nstar = params.d / (params.c * muEff);
%     Pstar = (params.a - params.e * Nstar) / (params.b * muEff * Nstar);
% end

J = [ -params.b * muEff * Pstar - params.e,  -params.b * muEff * Nstar; ...
       params.c * muEff * Pstar,            params.c * muEff * Nstar - params.d ];

trJ = trace(J);
detJ = det(J);
eigsJ = eig(J);

if detJ < 0
    className = 'saddle';
elseif trJ < 0 && detJ > 0
    disc = trJ^2 - 4*detJ;
    if disc < 0
        className = 'stable spiral sink';
    elseif disc > 0
        className = 'stable nodal sink';
    else
        className = 'degenerate stable node';
    end
elseif trJ > 0 && detJ > 0
    className = 'unstable';
else
    className = 'non-hyperbolic / borderline';
end

out = struct();
out.muEff = muEff;
out.Nstar = Nstar;
out.Pstar = Pstar;
out.J = J;
out.trace = trJ;
out.det = detJ;
out.eigs = eigsJ;
out.classification = className;
out.stabilityLine = struct('x', Nstar, 'y', Pstar);

if doPrint
    fprintf('\n--- Stability analysis ---\n');
    fprintf('mu_eff = %.6f\n', muEff);
    fprintf('Equilibrium: N* = %.6f, P* = %.6f\n', Nstar, Pstar);
    fprintf('Trace(J*)   = %.6f\n', trJ);
    fprintf('Det(J*)     = %.6f\n', detJ);
    fprintf('Eigenvalues = [%.6f%+.6fi, %.6f%+.6fi]\n', ...
        real(eigsJ(1)), imag(eigsJ(1)), real(eigsJ(2)), imag(eigsJ(2)));
    fprintf('Class       = %s\n', className);
end
end
