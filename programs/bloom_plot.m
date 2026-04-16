function bloom_plot(factors)
% BLOOM_PLOT  Simulate and plot a nutrient-phytoplankton Lotka-Volterra model.
%
% This reproduces the two figures in the prompt:
%   1) phase portrait in (N,P) with the theoretical equilibrium lines and
%      dots placed at whole-year intervals,
%   2) time series for nutrient and phytoplankton.
%
% The code is modular: growth terms are delegated to SPECIFIC_GROWTH and
% linearized stability information is handled by STABILITY_ANALYSIS.
%
% MATLAB R2024a compatible.

%% ----------------------- model parameters ------------------------------
% NOTE:
% The text in the manuscript writes a = 7.5e3, but the stated equilibrium
% P* = 7.5e-3 and the plotted scales are consistent with a = 7.5e-4.
% We therefore use the scaled value below to match the figures.
params.a = 7.5e-4;   % constant nutrient input
params.b = 1.0;      % exchange-rate modifier in N equation
params.c = 1.0;      % exchange-rate modifier in P equation
params.d = 0.1;      % proportional phytoplankton death rate
params.e = 0.0;      % proportional nutrient loss rate

%% ----------------------- growth configuration --------------------------
function opts = opts_from_factors(factors)
    factors = upper(string(factors));

    % ---------- master switches ----------
    hasT = contains(factors, "T");
    hasL = contains(factors, "L");
    hasG = contains(factors, "G");
    hasS = contains(factors, "S");
    hasD = contains(factors, "D");

    if factors == "BASE"
        hasT = false;
        hasL = false;
        hasG = false;
        hasS = false;
        hasD = false;
    end

    % ---------- main mode flags ----------
    opts.useTemperature         = hasT;
    opts.useLight               = hasL;

    % growth-specific flags apply only to enabled drivers
    opts.useTemperatureSpecific = hasT && hasG;
    opts.useLightSpecific       = hasL && hasG;

    % seasonality applies only to enabled drivers
    opts.useTemperatureSeasonal = hasT && hasS;
    opts.useLightSeasonal       = hasL && hasS;

    % depth only makes sense for light
    opts.useLightDepth          = hasL && hasD;

    % if both light and temperature are on, default to multiplicative coupling
    opts.combineAdditively      = false;
end
opts = opts_from_factors(factors);
% Temperature parameters
opts.temp.muConst = 0.5;    % half-saturating constant temperature growth
opts.temp.Tmin    = 0.0;
opts.temp.Topt    = 11.0;
opts.temp.Tmax    = 22.0;
opts.temp.Tamp    = 7.0;    % seasonal amplitude if enabled
opts.temp.period  = 365.0;  % days
opts.temp.phase   = 0.0;    % use 0 to start near spring crossing

% Light parameters
opts.light.muConst = 0.5;
opts.light.HL      = 120.0;  % half-saturation value
opts.light.I0      = 120.0;  % baseline intensity
opts.light.Iamp    = 60.0;
opts.light.period  = 365.0;
opts.light.phase   = 0.0;

% Depth / Beer-Lambert parameters
opts.depth.z       = 8.0;    % epilimnion depth
opts.depth.Iin     = 120.0;
opts.depth.Kbg     = 0.35;
opts.depth.kP      = 8.0;
opts.depth.nQuad   = 150;

%% ----------------------- time span / initial data ----------------------
baseX0 = [7.5e-2, 7.5e-3]; % 5e-2, 5e-4
scales = [1, 10, 100, -2];
X0 = scales(1) .* baseX0;

tFinal = 2500;                     % days

eq = stability_analysis(params, opts, false);

%% ----------------------- labeling ----------------------
function label = build_model_label(opts)
    if opts.useTemperature && opts.useLight
        if opts.combineAdditively
            base = 'Light+Temperature';
        else
            base = 'Light*Temperature';
        end
    elseif opts.useTemperature
        base = 'Temperature';
    elseif opts.useLight
        base = 'Light';
    else
        base = 'Base';
    end

    tags = {};
    if opts.useTemperatureSpecific || opts.useLightSpecific
        tags{end+1} = 'SG';
    end
    if opts.useLightDepth
        tags{end+1} = 'D';
    end
    if opts.useTemperatureSeasonal || opts.useLightSeasonal
        tags{end+1} = 'S';
    end

    if isempty(tags)
        label = base;
    else
        label = sprintf('%s (%s)', base, strjoin(tags, ', '));
    end
end

%% ----------------------- plotting -------------------------------------
fig = figure(1);
clf(fig);   % clear old contents
set(fig, 'Color', 'w', 'Position', [100 100 1100 430]);

% ---------------- phase portrait ----------------
subplot(1,2,1); hold on; box on;
[t, X] = ode45(@(t,x) lv_rhs(t, x, params, opts), 0:tFinal, X0.');
fprintf('%d\n',  size(X));
N = X(:,1);
P = X(:,2);
N = N / 5;
plot(N, P, 'LineWidth', 0.9);

% place dots at yearly marks
plotYears = 1:365:tFinal;
Nt = N(plotYears);
Pt = P(plotYears);
scatter(Nt, Pt, 'o', 'MarkerFaceColor', 'b');

% theoretical equilibrium lines
fprintf('%d\n', eq.Nstar / 5);
xline(eq.Nstar / 5, '-k', 'LineWidth', 1.0);
yline(eq.Pstar, '-k', 'LineWidth', 1.0);
plot(eq.Nstar, eq.Pstar, 'rp', 'MarkerSize', 8, 'MarkerFaceColor', 'r');

xlabel('Nutrient (N)');
ylabel('Phytoplankton (P)');
xlim([0 max(N)]);
ylim([0 max(P)]);
title(build_model_label(opts));

%-----------------analysis---------------------
nPeriod = estimate_period(t, N);
pPeriod = estimate_period(t, P);
fprintf('nPeriod %d pPeriod %d\n', nPeriod, pPeriod);

% ---------------- time series ----------------
subplot(1,2,2); hold on; box on;
yline(eq.Nstar / 5, '-b', 'LineWidth', 1.0, 'HandleVisibility', 'off'); 
yline(eq.Pstar, '-r', 'LineWidth', 1.0, 'HandleVisibility', 'off');
plot(t, N / 5, 'LineWidth', 1.0, 'DisplayName', 'Nutrient (N)');
plot(t, P, 'LineWidth', 1.0, 'DisplayName', 'Phytoplankton (P)');
for yr = 365:365:tFinal
    xline(yr, '-', 'Color', 0.82*[1 1 1], 'HandleVisibility','off');
end
xlabel('Time (Days)');
ylabel('Population');
% xlim([0 tFinal]);
ylim([0 max(max(N),max(P))]);
legend('Location','northeast');
title(build_model_label(opts));

hold off;
end

function T = estimate_period(t, y)
% Estimate period from successive local maxima.
% Returns mean spacing between peaks.

    y = y(:);
    t = t(:);

    peakIdx = [];
    for i = 2:length(y)-1
        if y(i) > y(i-1) && y(i) >= y(i+1)
            peakIdx(end+1) = i; %#ok<AGROW>
        end
    end

    % ignore first transient peaks if many exist
    if numel(peakIdx) >= 4
        peakIdx = peakIdx(2:end);
    end

    if numel(peakIdx) < 2
        T = NaN;
        return;
    end

    peakTimes = t(peakIdx);
    T = mean(diff(peakTimes));
end

function dx = lv_rhs(t, x, params, opts)
% Right-hand side of the nutrient-phytoplankton model.
N = x(1);
P = x(2);

% Optional guard against numerically induced negative states, remove w/ -2x
% Nsafe = max(N, 0);
% Psafe = max(P, 0);

mu = specific_growth(t, N, P, opts);
mu
dN = params.a - params.b * mu * N * P - params.e * N;
dP = params.c * mu * N * P - params.d * P;

dx = [dN; dP];
end
