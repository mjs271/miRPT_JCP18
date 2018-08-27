%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm as
% a function of the parameter, kappa.
% This script, will generate Figures A.13 and A.14 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-3;
dt = 1e-1;
maxtime = 1e0;
omega = 1;

% number of mobile and immobile particles
factor = 10;
Nm = 1e3;
Ni = ceil(Nm / factor);

% number of kappa values to test
num = 9;

% vector holding values of kappa to be tested (print to screen)
kappavec = linspace(0.9 / (num), 0.9, num)

% variance for Gaussian IC
sigma = 1e0 * sqrt(4 * D);

nsteps = floor(maxtime / dt);

stab_condM = (1 ./ Nm).^2 ./ (kappavec .* D * dt)
stab_condI = (1 ./ Ni).^2 ./ ((1 - kappavec) .* D * dt)
stab_condAvg = (stab_condM + stab_condI) / 2
stab_condGeo = sqrt(stab_condM .* stab_condI)

%% loop over values of kappa

% initialize arrays to store the errors as we go
%     row 1 is RMSE, row 2 is infinity norm
WIWMerr = zeros(2, num);
DMerr = zeros(2, num);

% immobile particle locations
imp = linspace(0, omega, Ni);

% mobile particle locations (randomly- or evenly-spaced)
mop = sort(omega * rand(1, Nm));
% mop = linspace(0, omega, Nm);

for i = 1 : num
    
%     assign current value of kappa
    kappa = kappavec(i);
    
%     LH heaviside IC
%     massmob = zeros(Nm, 1);
%     massmob(mop < 0.5) = 1;
%     analytic = 1 - normcdf(mop, 0.5, sqrt(2 * D * maxtime))';

%     gaussian IC
    massmob = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 * omega - mop').^2 / (2 * sigma^2)));
    analytic = (1 / sqrt(2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((0.5 * omega - mop').^2 / (2 * (sigma^2 + 2 * D * maxtime))));
    
%     Pairwise distance matrices
%     this has dimension Ni x Nm--shape corresponds to Wmmat
    dist = abs(bsxfun(@minus, mop, imp'));
%     this has dimension Nm x Nm--shape corresponds to DMmat
    mobdist = abs(bsxfun(@minus, mop, mop'));

%     diffusion operator matrix
    DMmat = (1 / sqrt(4 * pi * D * dt)) * exp(-((mobdist).^2 / (4 * D * dt)));
    DMmat = DMmat * diag(1./(sum(DMmat)));

    %     miRPT mass transfer matrices (W_I and W_M) calculatyed below
%         note that they are simulating different amounts of diffusion that
%         is scaled by kappa

%     encounter probability matrix
    Pmat = (1 / sqrt(kappa * 4 * pi * D * dt)) * exp(-((dist).^2 / (kappa * 4 * D * dt)));
    Pmat = Pmat * diag(1./(sum(Pmat)));

    WMmat = Pmat * diag(1./(sum(Pmat)));
    
%     encounter probability matrix
    Pmat = (1 / sqrt((1 - kappa) * 4 * pi * D * dt)) * exp(-((dist).^2 / ((1 - kappa) * 4 * D * dt)));
    Pmat = Pmat * diag(1./(sum(Pmat)));
    
    WImat = Pmat' * diag(1./(sum(Pmat, 2)));
    
%     make the requisite number of mass transfers by applying the weight
%     matrices/diffusion operator the approproate number of times
    WIWMmass = (WImat * WMmat)^nsteps * massmob;
    DMmass = DMmat^nsteps * massmob;
    
%     track the errors for each refinement
    WIWMerr(:, i) = ([sqrt(mean((WIWMmass - analytic).^2)); norm(WIWMmass - analytic, inf)]);
    DMerr(:, i) = ([sqrt(mean((DMmass - analytic).^2)); norm(DMmass - analytic, inf)]);
    
end

%% Plots

% spatial mass plot for final time step
figure(1)
clf
plot(mop, massmob, '--', 'LineWidth', 1.5)
hold on
plot(mop, WIWMmass, 'gd')
plot(mop, DMmass, 'bo')
plot(mop, analytic, 'LineWidth', 1.5)
xlabel('\textbf{Particle Position} $(x)$','Interpreter','latex', 'FontSize', 18)
ylabel('\textbf{Mass}','Interpreter','latex', 'FontSize', 18)
legend({'\textbf{IC}','\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{Analytic Solution}'},'Interpreter','latex', 'FontSize', 16,'Location','northeast')

% error/stability condition plot
figure(62)
clf
[hAx,hLine1,hLine2] = plotyy([kappavec', kappavec'], [WIWMerr(1, :)', DMerr(1, :)'], [linspace(0, 1, 100)', linspace(0, 1, 100)'], [1.0 * ones(1, 100)', 0.25 * ones(1, 100)'], @semilogy);
hLine2(1).Color = 'k';
hLine2(1).LineWidth = 1.5;
hLine2(2).Color = 'k';
hLine2(2).LineWidth = 1.5;
hLine2(2).LineStyle = '--';
hLine1(1).Color = 'b';
hLine1(1).Marker = 'o';
hLine1(1).LineWidth = 1.5;
hLine1(2).Color = 'r';
hLine1(2).Marker = '^';
hLine1(2).LineWidth = 1.5;
hold(hAx(1),'on')
hold(hAx(2),'on')
scatter(hAx(2), kappavec, stab_condAvg, 70, [0 0.7 0.4], 'filled')
scatter(hAx(2), kappavec, stab_condGeo, 70, [0.4 0.1 0.6], 'filled', 's')
legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{1.0}', '\textbf{0.25}', '\textbf{Arithmetic SC} $(\bar \eta)$', '\textbf{Geometric SC} $(\hat \eta)$'},'Interpreter','latex', 'FontSize', 20,'Location','northwest')
hAx(2).YColor = [0.1500    0.1500    0.1500];
xlabel('$\kappa^{(M)}$','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(1), '\textbf{RMSE}','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(2), '\textbf{Stability Condition}','Interpreter','latex', 'FontSize', 18)
hAx(1).YLim = [1e-2, 7e0];
hAx(2).YLim = [0.0, 5.3e0];
hAx(1).XLim = [1e-1, 1e0];
hAx(2).XLim = [1e-1, 1e0];
hAx(2).YTick = [0.25 1.0 5.0];
hAx(2).YTickLabel = {'0.25', '1.0', '5.0'};

