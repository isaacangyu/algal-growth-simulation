function bloom_plot(factors, delta, tFinal, phasePlot)
params.a = 7.5e-4;   % constant nutrient input
params.b = 1.0;      % exchange-rate modifier in N equation
params.c = 1.0;      % exchange-rate modifier in P equation
params.d = 0.1;      % proportional phytoplankton death rate
params.I = params.a * params.c / params.d^2;

%% ----------------------- growth configuration --------------------------
function opts = opts_from_factors(factors)
    factors = upper(string(factors));

    opts.flags.hasT = contains(factors, "T");
    opts.flags.hasL = contains(factors, "L");
    opts.flags.hasG = contains(factors, "G");
    opts.flags.hasS = contains(factors, "S");
    opts.flags.hasD = contains(factors, "D");

    opts.useTemperatureSpecific = opts.flags.hasT && opts.flags.hasG;
    opts.useLightSpecific       = opts.flags.hasL && opts.flags.hasG;
    opts.useTemperatureSeasonal = opts.flags.hasT && opts.flags.hasS;
    opts.useLightSeasonal       = opts.flags.hasL && opts.flags.hasS;
    opts.useLightDepth          = opts.flags.hasL && opts.flags.hasD;

    % if both light and temperature are on, default to multiplicative coupling
    % opts.combineAdditively      = false;
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
opts.depth.z       = 7.0; 
opts.depth.Iin     = 120.0;
opts.depth.Kbg     = 0.35;
opts.depth.kP      = 8.0;
opts.depth.nQuad   = 150;

%% ----------------------- setup ---------------------
eq = stability_analysis(params, opts, false);

baseX0 = [0.025 0.03]; 
% baseX0 = [eq.Nstar eq.Pstar];
% scales = [1, 10, 100, -2];
X0 = delta .* baseX0;

function dx = lv_rhs(t, x, params, opts)
    N = x(1);
    P = x(2);
    
    % Optional guard against numerically induced negative states
    % N = max(N, 0);
    % P = max(P, 0);
    
    mu = specific_growth(t, N, P, opts)
    % dN = params.a - params.b * mu * N * P;
    % dP = params.c * mu * N * P - params.d * P;
    I = params.a * params.c / params.d^2;
    dN = I - mu * N * P;
    dP = mu * N * P - P;
    
    f_scale = 1/10;
    dN = dN .* f_scale;
    dP = dP .* f_scale;
    
    dx = [dN; dP];
end

%% ----------------------- simulation -------------------------------------
fig = figure(1);
clf(fig);
set(fig, 'Color', 'w', 'Position', [100 100 1100 430]);

% ---------------- phase portrait ----------------
if phasePlot
    subplot(1,2,1); hold on; box on;
    h = 0.1;
    [t, X] = ode45(@(t,x) lv_rhs(t, x, params, opts), 0:h:tFinal, X0.');
    steps = length(X)
    N = X(:,1);
    P = X(:,2);
    N_scale = 1/10;
    N = N * N_scale;
    plot(N, P);
    
    % place dots at yearly marks
    period = opts.temp.period;
    plotYears = 1:period/h:steps;
    Nt = N(plotYears);
    Pt = P(plotYears);
    scatter(Nt, Pt, 'o', 'MarkerFaceColor', 'b');
    
    % nullclines (for non-autononmous mu)
    minN = min(N);
    minP = min(P);
    maxN = max(N);
    maxP = max(P);
    stepN = (maxN-minN)/100;
    n = minN:stepN:maxN;
    random_t = 300;
    mu = specific_growth(random_t, N, P, opts);
    n_nullcline = @(n) N_scale .* params.I ./ (mu .* n);
    p_nullcline = N_scale;
    plot(n, n_nullcline(n), '--b');
    xline(p_nullcline, '--r');
    
    % theoretical equilibrium lines (should be intersection)
    eq.Nstar = eq.Nstar / 10;
    fprintf('Nstar %d\n', eq.Nstar);
    fprintf('Pstar %d\n', eq.Pstar);
    % xline(eq.Nstar, '-k');
    % yline(eq.Pstar, '-k');
    plot(eq.Nstar, eq.Pstar, 'rp', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    xlabel('Nutrient (N)');
    ylabel('Phytoplankton (P)');
    xlim([minN maxN]);
    ylim([minP maxP]);
    title(build_model_label(opts));
end

% ---------------- time series ----------------
subplot(1,2,1+phasePlot); hold on; box on;
yline(eq.Nstar, '--b', 'HandleVisibility', 'off'); 
yline(eq.Pstar, '--r', 'HandleVisibility', 'off');
plot(t, N, 'DisplayName', 'Nutrient (N)');
plot(t, P, 'DisplayName', 'Phytoplankton (P)');
for yr = 365:365:tFinal
    xline(yr, '-', 'Color', 0.82*[1 1 1], 'HandleVisibility','off');
end
xlabel('Time (Days)');
ylabel('Population');
xlim([0 tFinal]);
ylim([0 max(maxN, maxP)]);
legend('Location','northeast');
title(build_model_label(opts));

%-----------------period---------------------
% subplot(1,3,2); hold on; box on;
nPeriod = estimate_period(t, N);
pPeriod = estimate_period(t, P);
fprintf('nPeriod %d pPeriod %d\n', round(nPeriod), round(pPeriod));
% plot(1:length(nPeriod), nPeriod);
% plot(1:length(pPeriod), pPeriod);

hold off;

% Save to ../plots (folder must exist)
filePath = filepath(opts, delta);
exportgraphics(fig, filePath, 'Resolution', 300);
fprintf('Plot saved to: %s\n', filePath);
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
    % if numel(peakIdx) >= 4
    %     peakIdx = peakIdx(2:end);
    % end
    % 
    % if numel(peakIdx) < 2
    %     T = NaN;
    %     return;
    % end

    peakTimes = t(peakIdx);
    T = diff(peakTimes);
    % T = mean(diff(peakTimes));
end

function label = build_model_label(opts)
    % if opts.flags.hasT && opts.flags.hasL
    %     if opts.combineAdditively
    %         base = 'Light+Temperature';
    %     else
    %         base = 'Light*Temperature';
    %     end
    if opts.flags.hasT
        base = 'Temperature';
    end
    if opts.flags.hasL
        base = 'Light';
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

function filePath = filepath(opts, delta)
    % First part: assume hasT or hasL is true
    % if opts.flags.hasT && opts.flags.hasL
    %     prefix = "temp_light";
    if opts.flags.hasT
        prefix = "temp";
    else
        prefix = "light";
    end
    
    % Assemble name parts
    nameParts = [prefix];
    if opts.flags.hasG, nameParts(end+1) = "growth"; end
    if opts.flags.hasS, nameParts(end+1) = "seasonality"; end
    if opts.flags.hasD, nameParts(end+1) = "depth"; end
    nameParts(end+1) = sprintf("%gx", delta);
    
    fileName = strjoin(nameParts, "_") + ".png";

    filePath = fullfile("../plots", fileName);
end