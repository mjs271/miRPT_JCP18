%==========================================================================
% This script conducts an analysis of the miRPT mass-transfer algorithm as
% a function of the number of mobile particles (Nm) for a 2D domain.
% This script uses a kD tree to form the sparse distance matrices and
% should be quite memory efficient.
% This script, will generate Figure 7 in:
%     "A Lagrangian Method for Reactive Transport with Solid/Aqueous
%     Chemical Phase Interaction," JCP 2018.
%==========================================================================

% set constants
D = 1e-2;
dt = 1e-1;
maxtime = 1e-1;
kappa = 0.5;

% filename to save workspace variables
filename = '2D_kD.mat';

% Ni = Nm / factor
factor = sqrt(0.5);

% number of refinements to make
num = 4;

% variance for Gaussian IC
sigma = 1e0 * sqrt(4 * D);

% vector holding (square root of) number of mobile particles
Nvec = 1e1 * 2.^linspace(2, num + 1, num);
% print vector of total number of mobile particles to screen
Nvec2 = Nvec.^2

% calculate and print stability condition to screen
stab_cond = (1 ./ ((min(Nvec, Nvec * factor))).^2) ./ (D * dt)

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

%     make sure to clear variable as we go to manage memory
    clear imlin molin

%     gaussian IC
    massmob = (1 / (2 * pi * sigma^2)) * exp(-((sqrt((0.5 - mopX).^2 + (0.5 - mopY).^2)).^2 / (2 * sigma^2)));
    analytic = (1 / (2 * pi * (sigma^2 + 2 * D * maxtime))) * exp(-((sqrt((0.5 - mopX).^2 + (0.5 - mopY).^2)).^2 / (2 * (sigma^2 + 2 * D * maxtime))));

%     turn the meshgrids into long vectors
    vecmopX = reshape(mopX, Nm^2, 1);
    vecmopY = reshape(mopY, Nm^2, 1);
    vecimpX = reshape(impX, Ni^2, 1);
    vecimpY = reshape(impY, Ni^2, 1);
    
    clear impX impY
    
    massmobvec = reshape(massmob, Nm^2, 1);
    analyticvec = reshape(analytic, Nm^2, 1);
    
%     search radius for kD tree range search
    cutdist = 4 * sqrt(2 * D * dt);
%     max bucket size for kD tree range search (1% of particle number is ad
%     hoc value that seems to work well)
    BS = max(ceil(1e-1 * max(Nm + Ni)), 1);
%     conduct the range search for mobile-immobile distances
    [MIidx, MIr] = rangesearch([vecmopX, vecmopY], [vecimpX, vecimpY], cutdist, 'BucketSize', BS);
    
    clear vecimpX vecimpY
    
%     preallocate the vectors used to build the sparse matrix
    Nclose = sum(cellfun('length', MIidx));
    MIrow = zeros(1, Nclose);
    MIcol = zeros(1, Nclose);
    MIval = zeros(1, Nclose);

%     fill the vectors used to build the sparse matrix
    count = 1;
    for ii = 1 : Ni^2
       for jj = 1 : length(MIidx{ii})
           
           MIrow(count) = ii;
           MIcol(count) = MIidx{ii}(jj);
           MIval(count) = MIr{ii}(jj);
           count = count + 1;
           
       end
    end
    
    clear MIidx MIr
    
%     sparse distance matrix for mobile-immobile distances
    Sdist = sparse(MIrow, MIcol, MIval);
    
    clear MIrow MIcol MIval
    
%     conduct the range search for mobile-mobile distances
    [MMidx, MMr]=rangesearch([vecmopX, vecmopY], [vecmopX, vecmopY], cutdist, 'BucketSize', BS);
    
    clear vecmopX vecmopY
    
%     preallocate the vectors used to build the sparse matrix
    Nclose = sum(cellfun('length', MMidx));
    MMrow = zeros(1, Nclose);
    MMcol = zeros(1, Nclose);
    MMval = zeros(1, Nclose);
    
%     fill the vectors used to build the sparse matrix
    count = 1;
    for ii = 1 : Nm^2
       for jj = 1 : length(MMidx{ii})
           
           MMrow(count) = ii;
           MMcol(count) = MMidx{ii}(jj);
           MMval(count) = MMr{ii}(jj);
           count = count + 1;
           
       end
    end
    
    clear MMidx MMr
    
%     sparse distance matrix for mobile-immobile distances
    Smobdist = sparse(MMrow, MMcol, MMval);
    
    clear MMrow MMcol MMval

%     diffusion operator matrix (this way of doing things is slow but
%     should be relatively memory efficient)
    DMmat = Smobdist;
    enzies = DMmat ~= 0;
    DMmat(enzies) = (1 / sqrt(4 * pi * D * dt)) * exp(-((Smobdist(enzies)).^2 / (4 * D * dt)));
    clear Smobdist
    temp = sum(DMmat);
    DMmat = DMmat * spdiags(1./(temp(:)), 0, Nm^2, Nm^2);

    % encounter probability matrix
    Pmat = Sdist;
    enzies2 = Pmat ~= 0;
    Pmat(enzies2) = (1 / sqrt(kappa * 4 * pi * D * dt)) * exp(-((Sdist(enzies2)).^2 / (kappa * 4 * D * dt)));
    clear Sdist
    temp = sum(Pmat);
    temp2 = sum(Pmat, 2);
    WMmat = Pmat * spdiags(1./(temp(:)), 0, Nm^2, Nm^2);
    WImat = Pmat' * spdiags(1./(temp2(:)), 0, Ni^2, Ni^2);
    
    clear Pmat temp temp2
    
%     make the requisite number of mass transfers by applying the weight
%     matrices/diffusion operator the approproate number of times
    WIWMmass = (WImat * WMmat)^nsteps * massmobvec;
    DMmass = DMmat^nsteps * massmobvec;
    
    clear massmobvec WImat WMmat DMmat
    
%   make the mass vectors sparse so we can make sparse matrices for
%   plotting
    WIWMmass = sparse(WIWMmass);
    DMmass = sparse(DMmass);
    
%     track the errors for each refinement
    WIWMerr(:, i) = ([sqrt(mean((WIWMmass - analyticvec).^2)); norm(WIWMmass - analyticvec, inf)]);
    DMerr(:, i) = ([sqrt(mean((DMmass - analyticvec).^2)); norm(DMmass - analyticvec, inf)]);
    
%     save as we go, in case of memory overflow crash
    save(filename)
    fprintf('saved \n')
   
end

% reshape the vectors into matrices for plotting
WIWMmass_mat = reshape(WIWMmass, Nm, Nm);
DMmass_mat = reshape(DMmass, Nm, Nm);

clear WIWMmass DMmass

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
figure(62)
clf
[hAx,hLine1,hLine2] = plotyy([Nvec.^2', Nvec.^2'], [WIWMerr(1, :)', DMerr(1, :)'], [logspace(3, 5, 100)'], [1.0 * ones(1, 100)'], @loglog);
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
scatter(hAx(2), Nvec.^2, stab_cond, 70, [0 0.7 0.4], 'filled')
legend({'\textbf{miRPT}', '\textbf{Diffusion Operator}', '\textbf{1.0}', '$\eta := \frac{(L/\min(N_I, N_M))^2}{D \Delta t}$'},'Interpreter','latex', 'FontSize', 20,'Location','southwest')
hAx(2).YColor = [0.1500    0.1500    0.1500];
hAx(1).Box = 'on';
xlabel('$N_M$','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(1), '\textbf{RMSE}','Interpreter','latex', 'FontSize', 18)
ylabel(hAx(2), '\textbf{Stability Condition}','Interpreter','latex', 'FontSize', 18)

