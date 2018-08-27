%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm in
% terms of variance growth over time.
% This script, will generate (something similar to) Figures 9 - 11 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-3;
dt = 1e-1;
maxtime = 1e0;
kappa = 0.5;

% Ni = Nm / factor;
factor = 1;

% number of refinements to make
num = 1;

% variance for Gaussian IC
sigma = 1e0 * sqrt(4 * D);

% number of mobile particles (print to screen)
Nm = 1e2 * 2^4

Ni = Nm / factor;

% calculate and print stability condition to screen
stab_cond = (1 ./ ((min(Nm, Ni))).^2) ./ (D * dt)

% number of time steps
nsteps = floor(maxtime / dt);

% boolean for running test case with analytic solution (Gaussian of Heaviside)
%     or non-analytic "noisy box"
analytic_sol = true;

%% loop over values of Nm

imp = linspace(0, 1, Ni);
mop = sort(rand(1, Nm));

if analytic_sol

% %     LH heaviside IC
%     massmob = zeros(Nm, 1);
%     massmob(mop < 0.5) = 1;
%     analytic = 1 - normcdf(mop, 0.5, sqrt(2 * D * maxtime))';

%     gaussian IC
    massmob = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 - mop').^2 / (2 * sigma^2)));
    analytic = (1 / sqrt(2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((0.5 - mop').^2 / (2 * (sigma^2 + 2 * D * maxtime))));

else

%     nosy box IC
    width = 200;
    massmob = zeros(Nm, 1);
    massmob(Nm/2 - width : Nm/2 + width) = rand(size(massmob(Nm/2 - width : Nm/2 + width)));

end

% calculate spatial mean and variance using trapezoidal quadrature
ICsum = trapz(mop, massmob);
norm_IC = massmob / ICsum;
m1_IC = trapz(mop, norm_IC .* mop');
initvar = trapz(mop, norm_IC .* (mop'- m1_IC).^2);

% Pairwise distance matrices
% this has dimension Ni x Nm--shape corresponds to Wmmat
dist = abs(bsxfun(@minus, mop, imp'));
% this has dimension Nm x Nm--shape corresponds to DMmat
mobdist = abs(bsxfun(@minus, mop, mop'));

% diffusion operator matrix
DMmat = (1 / sqrt(4 * pi * D * dt)) * exp(-((mobdist).^2 / (4 * D * dt)));
DMmat = DMmat * diag(1./(sum(DMmat)));

% encounter probability matrix
Pmat = (1 / sqrt(kappa * 4 * pi * D * dt)) * exp(-((dist).^2 / (kappa * 4 * D * dt)));

% mass transfer matrices (W_I and W_M)
WMmat = Pmat * diag(1./(sum(Pmat)));
WImat = Pmat' * diag(1./(sum(Pmat, 2)));

% vectors to store mobile particle mass as we step through time for
% mobile/immobile and diffusion operator approaches
WIWMmass = massmob;
DMmass = massmob;

% vectors for storing the spatial moments
m1_WIWM = zeros(nsteps, 1);
m2_WIWM = zeros(nsteps, 1);
m1_DM = zeros(nsteps, 1);
m2_DM = zeros(nsteps, 1);

% time-stepping loop
for i = 1 : nsteps

%     make the mass transfers
    WIWMmass = (WImat * WMmat) * WIWMmass;
    DMmass = DMmat * DMmass;
    
%     calculate spatial mean and variance using trapezoidal quadrature
    WIWMsum = trapz(mop, WIWMmass);
    norm_WIWM = WIWMmass / WIWMsum;
    m1_WIWM(i) = trapz(mop, norm_WIWM .* mop');
    m2_WIWM(i) = trapz(mop, norm_WIWM .* (mop'- m1_WIWM(i)).^2);
    DMsum = trapz(mop, DMmass);
    norm_DM = WIWMmass / DMsum;
    m1_DM(i) = trapz(mop, norm_DM .* mop');
    m2_DM(i) = trapz(mop, norm_DM .* (mop'- m1_DM(i)).^2);

end

%%  Plots

% time vector for plotting
tvec = linspace(dt, maxtime, nsteps);

% spatial mass plots for final time step
if analytic_sol
    
    figure(1)
    clf
    plot(mop, massmob, '--', 'LineWidth', 1.5)
    hold on
    % plot(mop, WMmass)
    plot(mop, WIWMmass, 'gd')
    plot(mop, DMmass, 'bo')
    plot(mop, analytic, 'LineWidth', 1.5)
    xlabel('\textbf{Particle Position} $(x)$','Interpreter','latex', 'FontSize', 18)
    ylabel('\textbf{Mass}','Interpreter','latex', 'FontSize', 18)
    legend({'\textbf{IC}','\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{Analytic Solution}'},'Interpreter','latex', 'FontSize', 16,'Location','northeast')
    
else
    
    figure(1)
    clf
    plot(mop, massmob, '*', 'LineWidth', 1.5)
    hold on
    % plot(mop, WMmass)
    plot(mop, WIWMmass, 'gd')
    plot(mop, DMmass, 'bo')
    % plot(mop, analytic, 'LineWidth', 1.5)
    xlabel('\textbf{Particle Position} $(x)$','Interpreter','latex', 'FontSize', 18)
    ylabel('\textbf{Mass}','Interpreter','latex', 'FontSize', 18)
    legend({'\textbf{IC}','\textbf{miRPT}', '\textbf{Diffusion Operator}'},'Interpreter','latex', 'FontSize', 16,'Location','northeast')
    
end

% variance growth plot
figure(2)
clf
loglog(tvec, m2_WIWM - initvar, 'b-o', 'LineWidth', 1.5)
hold on
loglog(tvec, m2_DM - initvar, 'r-^', 'LineWidth', 1.5)
loglog(tvec, 2 * D * tvec, 'Color', 'g', 'LineWidth', 2.5)
ylabel('\textbf{Spatial Variance Growth} ','Interpreter','latex', 'FontSize', 18)
xlabel('\textbf{Time} $(t)$','Interpreter','latex', 'FontSize', 18)
legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', '$2 D t$'},'Interpreter','latex', 'FontSize', 16,'Location','best')
