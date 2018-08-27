%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm as
% a function of the number of mobile particles (Nm) for randomly-scattered
% particles.
% This script, will generate (something similar to) Figure 8 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-3;
dt = 1e-3;
maxtime = 1e0;
kappa = 0.5;

% filename to save workspace variables
filename = 'gaussian_num8_ens100_scattered.mat';

% Ni = Nm / factor;
factor = 1;

% number of refinements to make
num = 8;

% number of instances in ensemble run
ens = 100;

% variance for Gaussian IC
sigma = sqrt(4 * D);

% vector holding number of mobile particles (print to screen)
Nvec = 1e2 * 2.^linspace(0, num - 1, num);

% number of time steps
nsteps = floor(maxtime / dt);

% calculate and print stability condition to screen
stab_cond = (1 ./ ((min(Nvec, Nvec ./ factor))).^2) ./ (D * dt)

%% loop overensemble and values of Nm

% cell arrays holding final masses/locations for every instance, so they
% can rebinned post hoc and then calculate the errors based on the binning
% structure
WIWMmass = cell(num, ens);
DMmass = cell(num, ens);
mop = cell(num, ens);

% ensemble loop
for j = 1 : ens

%     mobile particle number (Nm) loop
    for i = 1 : num

%         assign current number of mobile and immobile particles
        Nm = Nvec(i);
        Ni = Nm / factor;

        imp = linspace(0, 1, Ni);
        mop{i, j} = sort(rand(1, Nm));

% %         LH heaviside IC
%         massmob = zeros(Nm, 1);
%         massmob(mop < 0.5) = 1;
%         analytic = 1 - normcdf(mop, 0.5, sqrt(2 * D * maxtime))';

%         gaussian IC
        massmob = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 - mop{i, j}').^2 / (2 * sigma^2)));

%         Pairwise distance matrices
%         this has dimension Ni x Nm--shape corresponds to Wmmat
        dist = abs(bsxfun(@minus, mop{i, j}, imp'));
        mobdist = abs(bsxfun(@minus, mop{i, j}, mop{i, j}'));

%             diffusion operator matrix
        DMmat = (1 / sqrt(4 * pi * D * dt)) * exp(-((mobdist).^2 / (4 * D * dt)));
        DMmat = DMmat * diag(1./(sum(DMmat)));

%         encounter probability matrix
        Pmat = (1 / sqrt(kappa * 4 * pi * D * dt)) * exp(-((dist).^2 / (kappa * 4 * D * dt)));
        Pmat = Pmat * diag(1./(sum(Pmat)));

%         Pairwise distance matrices
%         this has dimension Ni x Nm--shape corresponds to Wmmat
        WMmat = Pmat * diag(1./(sum(Pmat)));
        WImat = Pmat' * diag(1./(sum(Pmat, 2)));

%         make the requisite number of mass transfers by applying the weight
%         matrices/diffusion operator the approproate number of times
%                 save the masses in cell arrays: i = num, j = ens
        WIWMmass{i, j} = (WImat * WMmat)^nsteps * massmob;
        DMmass{i, j} = DMmat^nsteps * massmob;

    end
    
end

%% Check for NaN's

nannys = cell2mat(cellfun(@(x)any(isnan(x)),WIWMmass,'UniformOutput',false));
if (sum(sum(nannys)) > 0)
    warning('NaN values in input')
%     [Nindex, Ensdex] = ind2sub(size(WIWMmass), find(nannys > 0))
    figure(100)
    clf
    imagesc(nannys)
    colorbar
end

%% Bin the masses, since positions were random

mobx = cell(1, num);
initmass = cell(1, num);
analytic = cell(1, num);

% number of spatial bins
bins = 2.5e2;

for i = 1 : num
    
    mobx{i} = linspace(0, 1, bins);
    initmass{i} = (1 / sqrt(2 * pi * sigma^2)) * exp(-((0.5 - mobx{i}).^2 / (2 * sigma^2)));
    analytic{i} = (1 / sqrt(2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((0.5 - mobx{i}).^2 / (2 * (sigma^2 + 2 * D * maxtime))));
    
    for j = 1 : ens

        bin = floor(mop{i, j}./mobx{i}(2)) + 1;
        
        binmass = zeros(1, bins);
        for k = 1 : Nvec(i)
            binmass(bin(k)) = binmass(bin(k)) + WIWMmass{i, j}(k);
        end
        WIWMmass{i, j} = binmass / (Nvec(i) / bins);
        
        binmass = zeros(1, bins);
        for k = 1 : Nvec(i)
            binmass(bin(k)) = binmass(bin(k)) + DMmass{i, j}(k);
        end
        DMmass{i, j} = binmass / (Nvec(i) / bins);

    end
end

%% Average over the ensemble

WIWMmass_avg = cell(1, num);
WIWMmass_std = cell(1, num);
DMmass_avg = cell(1, num);
DMmass_std = cell(1, num);

dim = bins;

for i = 1 : num
    
    WIWMtemp = reshape([WIWMmass{i, :}], dim, ens);
    DMtemp = reshape([DMmass{i, :}], dim, ens);

    WIWMmass_avg{i} = mean(WIWMtemp, 2);
    WIWMmass_std{i} = std(WIWMtemp, 0, 2);
    DMmass_avg{i} = mean(DMtemp, 2);
    DMmass_std{i} = std(DMtemp, 0, 2);

end

%% Compute Error for Averages

WIWMerr = zeros(2, num);
DMerr = zeros(2, num);

for i = 1 : num
    
    WIWMerr(:, i) = ([sqrt(mean((WIWMmass_avg{i} - analytic{i}').^2)); norm(WIWMmass_avg{i} - analytic{i}', inf)]);
    DMerr(:, i) = ([sqrt(mean((DMmass_avg{i} - analytic{i}').^2)); norm(DMmass_avg{i} - analytic{i}', inf)]);
    
end

%% save workspace and exit here for long/remote runs

% save(filename)
% exit

%%  Plots

% spatial mass plot for final time step
figure(1)
clf
plot(mobx{end}, initmass{end}, '--', 'LineWidth', 1.5)
hold on
plot(mobx{end}, WIWMmass_avg{end}, 'gd')
plot(mobx{end}, DMmass_avg{end}, 'bo')
plot(mobx{end}, analytic{end}, 'r', 'LineWidth', 1.5)
xlabel('\textbf{Particle Position} $(x)$','Interpreter','latex', 'FontSize', 18)
ylabel('\textbf{Mass}','Interpreter','latex', 'FontSize', 18)
legend({'\textbf{IC}','\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{Analytic Solution}'},'Interpreter','latex', 'FontSize', 16,'Location','northeast')

% error/stability condition plot
figure(2)
clf
[hAx,hLine1,hLine2] = plotyy([Nvec', Nvec'], [WIWMerr(1, :)', DMerr(1, :)'], [logspace(2, 4.35, 100)'], [1.0 * ones(1, 100)'], @loglog);
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
plot(hAx(2), logspace(2, 4.35, 100), 0.25 * ones(1, 100), 'k--', 'LineWidth', 1.5)
scatter(hAx(2), Nvec(2 : end), stab_cond(2 : end), 70, [0 0.7 0.4], 'filled')
lgd = legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{1.0}', '\textbf{0.25}', '\textbf{Stability Condition} $(\eta)$'},'Interpreter','latex', 'FontSize', 20,'Location','northeast');
lgd.Position = lgd.Position + [1e-2 1e-2 0 0];
hAx(2).YColor = [0.1500    0.1500    0.1500];
hAx(1).Box = 'on';
xlabel('$N_M$','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(1), '\textbf{RMSE}','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(2), '\textbf{Stability Condition}','Interpreter','latex', 'FontSize', 18)
