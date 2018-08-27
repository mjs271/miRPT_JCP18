program dolo_ADRE
use ADRE_mod
use PhreeqcRM
implicit none

! ==============================================================================
!                           SIMULATION PARAMETERS
! ==============================================================================
integer, parameter                 :: nthreads = 0
    ! number of OpenMP threads for reaction module. <= 0 implies the number of
    ! threads equals the number of processors of the computer
double precision, parameter        :: maxtime = 60e3_dp ! 1000 MIN
double precision, parameter        :: Omega = 0.5d0 ! length of domain [m]
double precision, parameter        :: dx = 1d-3 ! spatial discretization parameter
integer, parameter                 :: ncell = nint(Omega/dx) - 1
    ! number of chemistry cells
    ! subtract 1 because won't be calculating chemistry for boundary cell
integer, parameter                 :: ntrans = ncell + 1
    ! number of cells for transport
double precision, parameter        :: dt = 1e1_dp ! time step
integer, parameter                 :: nsteps = nint(maxtime/dt) ! number of time steps
double precision, parameter        :: save_dt = dt * 1e3_dp
integer, parameter                 :: save_steps = nint(maxtime/(save_dt)) + 1
    ! time step for saving concentrations for plotting and number of writes to be done

! ==============================================================================
!                           PHYSICAL PARAMETERS
! ==============================================================================
double precision, parameter        :: darvel = 1.2e-5_dp ! Darcy flux [m/s]
double precision, parameter        :: porosity = 0.5d0
double precision, parameter        :: init_v = darvel/porosity
    ! porosity advective velocity
double precision, parameter        :: alpha_l = 0.005d0 ! longitudinal dispersivity
double precision, parameter        :: init_D = alpha_l*init_v
double precision                   :: v(ntrans - 1) = init_v, D(ntrans - 1) = init_D
    ! vectors for velocity and diffusion coefficients

! ==============================================================================
!                       GENERAL/PHREEQCRM VARIABLES
! ==============================================================================
double precision                   :: init_calcite, na_inflow, mg_inflow,&
                                      ca_inflow, cl_inflow, co2_inflow
double precision, dimension(ncell) :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0
double precision                   :: cur_time
integer                            :: i, j, m, id, status, ngrd
integer                            :: ncomp, nchem, so_col
integer                            :: ic1(ncell, 7)
integer, dimension(ncell)          :: mask
double precision, allocatable      :: bc_conc(:, :), comp_conc(:, :),&
                                      sout(:, :), concs(:, :),&
                                      plot_concs(: , :), plot_times(:)
integer                            :: bc1(1) ! this is for one boundary condition

cur_time = 0.0d0

! print some summary info, to begin with
print *, '======================================='
print *, 'grid_Pe =', (dx * init_v) / init_D
print *, 'CFL = ', (init_v*dt)/dx
print *, '1.0d0 / dx = ', 1.0d0 / dx
print *, 'dt = ', dt
print *, 'ncell (chemistry) = ', ncell
print *, '======================================='

! initial concentration in domain [mol/L]
init_calcite = .270865d0 ! = 270.865 mol/m^3

! inflow concentration of CO2
! this is just a big number so it is saturated
co2_inflow = 10.0d0 ! this is from Leal's input file

! inflow concentrations [mol/L]
na_inflow = 0.884882d0
mg_inflow = 0.0491601d0
ca_inflow = 0.00983202d0
cl_inflow = 1.00287d0

! write the PhreeqcRM input file
call phreeqc_input(init_calcite, na_inflow, mg_inflow, ca_inflow, cl_inflow,&
                   co2_inflow)

! ==============================================================================
!                   DEFINE PHYSICAL CONDITIONS (for PHREEQCRM)
! ==============================================================================
den_0 = 1.0d0 ! Density
prs_0 = 98.6923d0 ! Pressure 98.6923 atm = 100 bar
tmp_0 = 60.0d0 ! Temperature
sat_0 = 1.0d0 ! Saturation
por_0 = porosity ! Porosity
vol_0 = (Omega * 1e3_dp)/dble(ncell) ! Representative Volume of cell (L), default = 1

! print chemistry mask to print detailed output (headings) for only first element
mask = 0
mask(1) = 1

id = RM_Create(ncell, nthreads)
! local database location
status = RM_LoadDatabase(id, '/usr/local/share/doc/phreeqcrm/database/phreeqc.dat')
if (status < 0) then
    print *, 'Database Load Error'
    call exit(status)
endif

! initialize all the settings for the reaction module
status = RM_OpenFiles(id)
status = RM_SetRepresentativeVolume(id, vol_0)
status = RM_SetSaturation(id, sat_0)
status = RM_SetPorosity(id, por_0)
status = RM_SetTime(id, cur_time)
status = RM_SetTimeStep(id, dt)
status = RM_SetDensity(id, den_0)
status = RM_SetTemperature(id, tmp_0)
status = RM_SetPressure(id, prs_0)
status = RM_SetSelectedOutputOn(id, 1) ! turn on selected output to get dolomite/calcite values
status = RM_SetUnitsSolution(id, 2) ! 2 = mol/L
status = RM_SetUnitsPPassemblage(id, 0) ! 0 = mol/L of representative volume
status = RM_SetPrintChemistryMask(id, mask)
status = RM_SetPrintChemistryOn(id, 0, 0, 0)
status = RM_SetScreenOn(id, 0) ! turn off useless messages about rebalancing
status = RM_SetComponentH2O(id, 1)
    ! 0 implies that H2O is not given as a separate component--just H and O
    ! 1 is default with H2O as separate component
status = RM_RunFile(id, 1, 1, 1, 'dolomite_chem.in')

! the following isn't strictly necessary if you already know what your problem
! looks like, but it's quick and only done once.

ncomp = RM_FindComponents(id)
! check if the hard-coded nspec is correct
if (ncomp /= nspec) then
    print *, '***** hard-coded nspec is different from ncomp ******'
    call exit(status)
endif

! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= ncell) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)
if (nchem /= ncell) then
    print *, '****** Chem count differs from particle number ******'
    call exit(status)
endif

! ==============================================================================
!                       INITIAL/BOUNDARY CONDITIONS
! ==============================================================================
! ===  (1) SOLUTIONS, (2) EQUILIBRIUM_PHASES, (3) EXCHANGE,
! ===  (4) SURFACE, (5) GAS_PHASE, (6) SOLID_SOLUTIONS, and (7) KINETICS

ic1 = -1
ic1(:, 1) = 1
ic1(:, 2) = 1

allocate(bc_conc(1, ncomp))
    ! must be 2D array for module--requires singleton dimension in this case
bc1 = 0 ! corresponds to solution zero in input file

status = RM_InitialPhreeqc2Module(id, ic1)
! ==============================================================================

! the following commented blocks will give you the names of the components and
! headings for selected output.
! They are left here purely for reference, and will not run unless you
! instantiate the corresponding variables.

! ! this makes a list of the names of components
! allocate(comp_list(ncomp))
! ! print *, 'ncomp = ', ncomp
! do i = 1, ncomp
!     status = RM_GetComponent(id, i, tempname)
!     allocate(character(len_trim(tempname)) :: comp_list(i)%comp)
!     comp_list(i)%comp = trim(tempname)
!     print *, comp_list(i)%comp
! enddo

status = RM_RunCells(id)

! so_col = RM_GetSelectedOutputcolumncount(id)
! so_row = RM_GetSelectedOutputrowcount(id)
! allocate(sout(so_row, so_col))
! status = RM_GetSelectedOutput(id, sout)

! this makes a list of the names of selected output elements
! allocate(head_list(so_col))
! do i = 1, so_col
!     status = RM_GetSelectedOutputheading(id, i, tempname)
!     allocate(character(len_trim(tempname)) :: head_list(i)%head)
!     head_list(i)%head = trim(tempname)
!     print *, head_list(i)%head
! enddo

! get ion concentrations from the reaction module
allocate(comp_conc(ncell, ncomp))
status = RM_GetConcentrations(id, comp_conc)
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

! ==============================================================================
!               SELECTED OUTPUT (based on the PHREEQC input file)
! ==============================================================================
! == (1) pH, (2) calcite, (3) change calcite, (4) dolomite, (5), change dolomite
! ==============================================================================

! ==============================================================================
!            PLOT_CONCS (the quantities of interest for the problem)
! ==============================================================================
! == (1) pH, (2) calcite, (3) dolomite
! ==============================================================================

! ==============================================================================
!        CONCS (what is being transported and given to the reaction module)
! ==============================================================================
! dimension of concs is ntrans (# chemistry cells + 1 for boundary) x ncomp

! For this specific problem:

! If (RM_SetComponentH2O(id, 0)):
! (1) H, (2) O, (3) charge, (4) C, (5) Ca, (6) Cl, (7) Mg, (8) Na, (9) Si

! If (RM_SetComponentH2O(id, 1)):
! (1) H2O, (2) H, (3) O, (4) charge, (5) C, (6) Ca, (7) Cl, (8) Mg, (9) Na, (10) Si
! ==============================================================================

so_col = RM_GetSelectedOutputcolumncount(id)
allocate(sout(ncell, so_col))
status = RM_GetSelectedOutput(id, sout)

allocate(concs(ntrans, ncomp), plot_concs(ntrans, 3), plot_times(save_steps))
concs(1, :) = bc_conc(1, :)
concs(2 : ntrans, :) = comp_conc

plot_concs(1, :) = 0.0d0
plot_concs(2 : ntrans, 1) = sout(:, 1)
plot_concs(2 : ntrans, 2) = sout(:, 2) / vol_0
plot_concs(2 : ntrans, 3) = sout(:, 4) / vol_0

plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)

! this is the file that stores what comes from plot_concs
! the header gives the shape of the full array and the number of time steps
open (unit=12, file='time_concs.txt', action='write')
write (12, *) shape(plot_concs), save_steps
write (12, *) plot_concs

! counter for plot_times
j = 2

! time stepping loop
do m = 1, nsteps

    print *, 'step = ', m, ' of ', nsteps

    call advect(concs(:, :), v, dx, dt, ntrans)
    call diffuse(concs(:, :), D, dx, dt, ntrans)

    cur_time = cur_time + dt
    status = RM_SetTime(id, cur_time)
    comp_conc = concs(2 : ntrans, :)
    status = RM_SetConcentrations(id, comp_conc)
    status = RM_RunCells(id)
    status = RM_GetConcentrations(id, comp_conc)
    concs(2 : ntrans, :) = comp_conc

    if (cur_time >= plot_times(j)) then
        status = RM_GetSelectedOutput(id, sout)
        plot_concs(2 : ntrans, 1) = sout(:, 1)
        plot_concs(2 : ntrans, 2) = sout(:, 2) / vol_0
        plot_concs(2 : ntrans, 3) = sout(:, 4) / vol_0
        write (12, *) plot_concs
        j = j + 1
    endif

enddo

close (unit=12, status='keep')

status = RM_Destroy(id)
deallocate(concs, plot_concs, plot_times, sout, comp_conc, bc_conc)

end program dolo_ADRE
