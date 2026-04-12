function mu = specific_growth(t, N, P, opts, mode)
% SPECIFIC_GROWTH  Effective growth multiplier for temperature/light suites.
%
% mu = specific_growth(t, N, P, opts)
% mu = specific_growth(t, N, P, opts, mode)
%
% Inputs
%   t    : time in days
%   N,P  : state variables (P used in depth-dependent light attenuation)
%   opts : configuration structure
%   mode : 'dynamic'  -> evaluate full model at time t (default)
%          'baseline' -> evaluate seasonality-controlled average reference
%
% Supported switches
%   Temperature: constant / specific-growth / seasonality
%   Light:       constant / specific-growth / depth / seasonality
%
% Conventions requested in the prompt:
%   - For seasonality, control against the baseline.
%   - For temperature, the non-seasonal specific-growth baseline is
%     integral(mu_T)/length over the survivable interval.
%   - For light depth, the baseline is the average over the epilimnion.

if nargin < 5
    mode = 'dynamic';
end

% ---------- temperature contribution ----------
muT = 1.0;
if isfield(opts, 'useTemperature') && opts.useTemperature
    if isfield(opts, 'useTemperatureSpecific') && opts.useTemperatureSpecific
        if strcmpi(mode, 'baseline') || ~opts.useTemperatureSeasonal
            muT = temp_specific_baseline(opts.temp);
        else
            T = opts.temp.Topt + opts.temp.Tamp .* ...
                sin(2*pi*t/opts.temp.period + opts.temp.phase);
            muT = temp_specific_curve(T, opts.temp);
        end
    else
        if strcmpi(mode, 'baseline') || ~opts.useTemperatureSeasonal
            muT = opts.temp.muConst;
        else
            % Seasonal control around the same baseline mean.
            muT = opts.temp.muConst + opts.temp.muConst .* ...
                sin(2*pi*t/opts.temp.period + opts.temp.phase);
            muT = max(muT, 0);
            muT = min(muT, 1);
        end
    end
end

% ---------- light contribution ----------
muL = 1.0;
if isfield(opts, 'useLight') && opts.useLight
    if isfield(opts, 'useLightDepth') && opts.useLightDepth
        % Beer-Lambert depth average or Monod(depth-light) average.
        if strcmpi(mode, 'baseline')
            Puse = 0;
        else
            Puse = max(P, 0);
        end
        muL = depth_light_average(Puse, opts);
    elseif isfield(opts, 'useLightSpecific') && opts.useLightSpecific
        if strcmpi(mode, 'baseline') || ~opts.useLightSeasonal
            muL = light_specific_curve(opts.light.I0, opts.light.HL);
        else
            I = opts.light.I0 + opts.light.Iamp .* ...
                sin(2*pi*t/opts.light.period + opts.light.phase);
            I = max(I, 0);
            muL = light_specific_curve(I, opts.light.HL);
        end
    else
        if strcmpi(mode, 'baseline') || ~opts.useLightSeasonal
            muL = opts.light.muConst;
        else
            muL = opts.light.muConst + opts.light.muConst .* ...
                sin(2*pi*t/opts.light.period + opts.light.phase);
            muL = max(muL, 0);
            muL = min(muL, 1);
        end
    end
end

% ---------- combine terms ----------
if isfield(opts, 'useTemperature') && opts.useTemperature && ...
   isfield(opts, 'useLight') && opts.useLight
    if isfield(opts, 'combineAdditively') && opts.combineAdditively
        mu = 0.5 * (muT + muL);
    else
        mu = muT .* muL;
    end
elseif isfield(opts, 'useTemperature') && opts.useTemperature
    mu = muT;
elseif isfield(opts, 'useLight') && opts.useLight
    mu = muL;
else
    mu = 1.0;
end

mu = max(0, min(1, mu));
end

function mu = temp_specific_curve(T, temp)
% Left-skewed beta-like curve on [Tmin, Tmax], peaking at Topt.
Tmin = temp.Tmin; Topt = temp.Topt; Tmax = temp.Tmax;
mu = zeros(size(T));
mask = (T > Tmin) & (T < Tmax);

% Exponents chosen so the left tail is broader than the right tail.
alpha = 1.6;
beta  = 0.9;

x = (T(mask) - Tmin) ./ (Tmax - Tmin);
xopt = (Topt - Tmin) / (Tmax - Tmin);

raw = (x ./ xopt) .^ alpha .* ((1 - x) ./ (1 - xopt)) .^ beta;
mu(mask) = raw;
mu = max(mu, 0);
mu = mu ./ max(mu(:) + eps);
end

function muBar = temp_specific_baseline(temp)
% Average specific growth over the survivable interval.
Tf = @(T) temp_specific_curve(T, temp);
muBar = integral(Tf, temp.Tmin, temp.Tmax, 'ArrayValued', true) / (temp.Tmax - temp.Tmin);
end

function mu = light_specific_curve(I, HL)
% Monod kinetics.
mu = I ./ (I + HL);
end

function mu = depth_light_average(P, opts)
% Average light or Monod-transformed light over depth.
z   = opts.depth.z;
Iin = opts.depth.Iin;
Kbg = opts.depth.Kbg;
kP  = opts.depth.kP;
HL  = opts.light.HL;
ns  = opts.depth.nQuad;

s = linspace(0, z, ns);
L = Iin .* exp((-Kbg - kP .* P) .* s);

if isfield(opts, 'useLightSpecific') && opts.useLightSpecific
    integrand = L ./ (L + HL);
else
    integrand = L ./ max(Iin, eps);
end

mu = trapz(s, integrand) / z;
mu = max(0, min(1, mu));
end
