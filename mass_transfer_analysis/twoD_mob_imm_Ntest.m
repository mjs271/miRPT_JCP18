%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm as
% a function of the number of mobile particles (Nm) for a 2D domain.
% This script forms the distance matrices naively and results in dense
% matrices that will overflow memory for large particle numbers.
% This script, will generate Figure 7 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-2;
dt = 0.1;
maxtime = 1e0;
kappa = 0.5;

% filename to save workspace variables
filename = '2D_num6_D1e-3_dt1e-1_factor9e-1.mat';

% Ni = Nm / factor
factor = sqrt(0.5);

% number of refinements to make
num = 4;

% variance for Gaussian IC
sigma = 1e-1 * sqrt(4 * D);

% vector holding (square root of) number of mobile particles
Nvec = 4e1 * 2.^linspace(0, num - 1, num);
% print vector of total number of mobile particles to screen
Nvec2 = Nvec.^2

% calculate and print stability condition to screen
stab_cond = (1 ./ ((min(Nvec, ceil(Nvec * factor)))).^2) ./ (D * dt)

% number of time steps
nsteps = floor(maxtime / dt);

%% loop over values of Nm

% initialize arrays to store the errors as we go
%     row 1 is RMSE, row 2 is infinity norm
WIWMerr = zeros(2, num);
DMerr = zeros(2, num);

for i = 1 : num
    
%     assign current (square root of) number of mobile and immobile particles
    Nm = Nvec(i);
    Ni = ceil(Nm * factor);
    
%     linearly-spaced vectors used to generate the immobile and mobile
%     meshgrids
    imlin = linspace(0, 1, Ni);
    molin = linspace(0, 1, Nm);
    
%     spatial meshgrids for immobile and mobile particle locations
    [impX, impY] = meshgrid(imlin, imlin);
    [mopX, mopY] = meshgrid(molin, molin);

%     gaussian IC
    massmob = (1 / (2 * pi * sigma^2)) * exp(-((sqrt((0.5 - mopX).^2 + (0.5 - mopY).^2)).^2 / (2 * sigma^2)));
    analytic = (1 / (2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((sqrt((0.5 - mopX).^2 + (0.5 - mopY).^2)).^2 / (2 * (sigma^2 + 2 * D * maxtime))));
    
%     turn the meshgrids into long vectors
    vecmopX = reshape(mopX, Nm^2, 1);
    vecmopY = reshape(mopY, Nm^2, 1);
    vecimpX = reshape(impX, Ni^2, 1);
    vecimpY = reshape(impY, Ni^2, 1);
    massmobvec = reshape(massmob, Nm^2, 1);
    analyticvec = reshape(analytic, Nm^2, 1);
    
%     Pairwise distance matrices
%     this has dimension Ni x Nm--shape corresponds to Wmmat
    distX = abs(bsxfun(@minus, vecmopX', vecimpX));
    distY = abs(bsxfun(@minus, vecmopY', vecimpY));
    dist = sqrt(distX.^2 + distY.^2);
%     this has dimension Nm x Nm--shape corresponds to DMmat
    mobdistX = abs(bsxfun(@minus, vecmopX', vecmopX));
    mobdistY = abs(bsxfun(@minus, vecmopY', vecmopY));
    mobdist = sqrt(mobdistX.^2 + mobdistY.^2);

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
    WIWMmass = (WImat * WMmat)^nsteps * massmobvec;
    DMmass = DMmat^nsteps * massmobvec;
    
%     track the errors for each refinement
    WIWMerr(:, i) = ([sqrt(mean((WIWMmass - analyticvec).^2)); norm(WIWMmass - analyticvec, inf)]);
    DMerr(:, i) = ([sqrt(mean((DMmass - analyticvec).^2)); norm(DMmass - analyticvec, inf)]);
    
end

% reshape the vectors into matrices for plotting
WIWMmass_mat = reshape(WIWMmass, Nm, Nm);
DMmass_mat = reshape(DMmass, Nm, Nm);

%% save workspace and exit here for long/remote runs

% save(filename)
% exit

%% Plots

% spatial mass plot for final time step
figure(1)
clf
subplot(2, 2, 1)
surf(mopX, mopY, massmob)
shading interp
title('IC','Interpreter','latex', 'FontSize', 20)
subplot(2, 2, 2)
surf(mopX, mopY, analytic)
shading interp
title('Analytic','Interpreter','latex', 'FontSize', 20)
subplot(2, 2, 3)
surf(mopX, mopY, WIWMmass_mat)
shading interp
title('miRPT','Interpreter','latex', 'FontSize', 20)
subplot(2, 2, 4)
surf(mopX, mopY, DMmass_mat)
shading interp
title('Diffusion Operator','Interpreter','latex', 'FontSize', 20)

% error/stability condition plot
figure(2)
clf
[hAx,hLine1,hLine2] = plotyy([Nvec.^2', Nvec.^2'], [WIWMerr(1, :)', DMerr(1, :)'], [logspace(1, 4, 100)'], [1.0 * ones(1, 100)'], @loglog);
hLine2(1).Color = 'k';
hLine2(1).LineWidth = 1.5;
hLine1(1).Color = 'b';
hLine1(1).Marker = 'o';
hLine1(1).LineWidth = 1.5;
hLine1(2).Color = 'r';
hLine1(2).Marker = '^';
hLine1(2).LineWidth = 1.5;
hold(hAx(1),'on')
hold(hAx(2),'on')
plot(hAx(1), logspace(1, 4, 100), D * dt * ones(1, 100), '--', 'LineWidth', 1.5, 'Color', [0.6, 0.2, 0.8])
scatter(hAx(2), Nvec.^2, stab_cond, 70, [0 0.7 0.4], 'filled')
legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', ['$D \Delta t$ = ' num2str(D * dt)], '\textbf{1.0}', '$\eta := \frac{(L/\min(\sqrt{N_I}, \sqrt{N_M}))^2}{D \Delta t}$'},'Interpreter','latex', 'FontSize', 20,'Location','northeast')
hAx(2).YColor = [0.1500    0.1500    0.1500];
hAx(1).Box = 'on';
xlabel('$N_M$','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(1), '\textbf{RMSE}','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(2), '\textbf{Stability Condition}','Interpreter','latex', 'FontSize', 18)
hAx(2).YLim = [0.05, 1e2];

