%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm as
% a function of the number of mobile particles (Nm).
% This script, will generate Figures 4 and 5 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-3;
dt = 1e-1;
maxtime = 1e0;
kappa = 0.5;

% Ni = Nm / factor
factor = 2;

% number of refinements to make
num = 6;

% variance for Gaussian IC
sigma = 1e0 * sqrt(4 * D);

% vector holding number of mobile particles (print to screen)
Nvec = 1e2 * 2.^linspace(0, num - 1, num)

% calculate and print stability condition to screen
stab_cond = (1 ./ ((min(Nvec, Nvec ./ factor))).^2) ./ (D * dt)

% number of time steps
nsteps = floor(maxtime / dt);

%% loop over values of Nm

% initialize arrays to store the errors as we go
%     row 1 is RMSE, row 2 is infinity norm
WIWMerr = zeros(2, num);
DMerr = zeros(2, num);

for i = 1 : num
    
%     assign current number of mobile and immobile particles
    Nm = Nvec(i);
    Ni = Nm / factor;
    
%     immobile particle locations
    imp = linspace(0, 1, Ni);
%     mobile particle locations (randomly- or evenly-spaced)
%     mop = sort(rand(1, Nm));
    mop = linspace(0, 1, Nm);
    
% %     Heaviside IC
%     massmob = zeros(Nm, 1);
%     massmob(mop < 0.5) = 1;
%     analytic = 1 - normcdf(mop, 0.5, sqrt(2 * D * maxtime))';

%     gaussian IC
    massmob = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 - mop').^2 / (2 * sigma^2)));
    analytic = (1 / sqrt(2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((0.5 - mop').^2 / (2 * (sigma^2 + 2 * D * maxtime))));
    
%     Pairwise distance matrices
%     this has dimension Ni x Nm--shape corresponds to Wmmat
    dist = abs(bsxfun(@minus, mop, imp'));
%     this has dimension Nm x Nm--shape corresponds to DMmat
    mobdist = abs(bsxfun(@minus, mop, mop'));

%     diffusion operator matrix
    DMmat = (1 / sqrt(4 * pi * D * dt)) * exp(-((mobdist).^2 / (4 * D * dt)));
    DMmat = DMmat * diag(1./(sum(DMmat)));

%     encounter probability matrix
    Pmat = (1 / sqrt(kappa * 4 * pi * D * dt)) * exp(-((dist).^2 / (kappa * 4 * D * dt)));
    Pmat = Pmat * diag(1./(sum(Pmat)));

%     miRPT mass transfer matrices (W_I and W_M)
    WMmat = Pmat * diag(1./(sum(Pmat)));
    WImat = Pmat' * diag(1./(sum(Pmat, 2)));
    
%     make the requisite number of mass transfers by applying the weight
%     matrices/diffusion operator the approproate number of times
    WIWMmass = (WImat * WMmat)^nsteps * massmob;
    DMmass = DMmat^nsteps * massmob;
    
%     track the errors for each refinement
    WIWMerr(:, i) = ([sqrt(mean((WIWMmass - analytic).^2)); norm(WIWMmass - analytic, inf)]);
    DMerr(:, i) = ([sqrt(mean((DMmass - analytic).^2)); norm(DMmass - analytic, inf)]);
    
end

%%  Plots

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
figure(2)
clf
[hAx,hLine1,hLine2] = plotyy([Nvec', Nvec'], [WIWMerr(1, :)', DMerr(1, :)'], logspace(2, 4, 100)', 1.0 * ones(1, 100)', @loglog);
hLine2.Color = 'k';
hLine2.LineWidth = 1.5;
hLine1(1).Color = 'b';
hLine1(1).Marker = 'o';
hLine1(1).LineWidth = 1.5;
hLine1(2).Color = 'r';
hLine1(2).Marker = '^';
hLine1(2).LineWidth = 1.5;
hold(hAx(1),'on')
hold(hAx(2),'on')
scatter(hAx(2), Nvec, stab_cond, 70, [0 0.7 0.4], 'filled')
legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{1.0}', '\textbf{Stability Condition} $(\eta)$'},'Interpreter','latex', 'FontSize', 20,'Location','east')
hAx(2).YColor = [0.1500    0.1500    0.1500];
hAx(1).Box = 'on';
xlabel('$N_M$','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(1), '\textbf{RMSE}','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(2), '\textbf{Stability Condition}','Interpreter','latex', 'FontSize', 18)
