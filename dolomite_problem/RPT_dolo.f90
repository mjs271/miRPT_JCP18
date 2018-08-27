program RPT_dolo
use RPT_mod
use PhreeqcRM
implicit none

! ==============================================================================
!                           SIMULATION PARAMETERS
! ==============================================================================
integer, parameter          :: nthreads = 0
    ! number of OpenMP threads for reaction module. <= 0 implies the number of
    ! threads equals the number of processors of the computer
double precision, parameter :: maxtime = 60d3 ! 1000 MIN
double precision, parameter :: lowx = 0.0d0, upx = 0.5d0, omega = upx - lowx
    ! domain length, lower, upper bounds
double precision, parameter :: dt = 1e0_dp ! time step
integer, parameter          :: nsteps = nint(maxtime/dt)! number of time steps
integer, parameter          :: nipart = 3000, ncellv = 3000
    ! number of immobile particles and velocity field cells
double precision, parameter :: dxv = omega/dble(ncellv), dx = omega/dble(nipart)
    ! dx for velocity grid, and for calculating boundary injection
integer, parameter          :: nptot = 6000
    ! number of mobile particles
double precision, parameter :: save_dt = dt * 1e2_dp
integer, parameter          :: save_steps = nint(maxtime/(save_dt)) + 1
    ! time step for saving concentrations for plotting and number of writes to be done

! ==============================================================================
!                           PHYSICAL PARAMETERS
! ==============================================================================
double precision, parameter :: darvel = 1.2e-5_dp ! Darcy flux [m/s]
double precision, parameter :: init_porosity = 0.5d0
double precision, parameter :: init_v = darvel/init_porosity
double precision, parameter :: alpha_l = 0.005d0 ! longitudinal dispersivity
double precision, parameter :: init_D = alpha_l * init_v
double precision            :: v(ncellv) = init_v, D(ncellv) = init_D, repvol
    ! vectors for velocity and diffusion coefficients

! ==============================================================================
!                               GENERAL VARIABLES
! ==============================================================================
type(mparticle)               :: mparts(nptot) ! mobile particles
type(iparticle)               :: iparts(nipart) ! immobile particles
integer                       :: indices(nptot) ! array for easy indexing
integer, allocatable          :: alive(:), nalive(:)
    ! arrays for indexing to alive particles
integer                       :: nactive, nexit
    ! number of active mobile particles and how many exited in time step
double precision              :: Pe_grid, CFL
double precision, allocatable :: Distmat(:, :)
    ! pairwise distance matrix for mass transfers

! ==============================================================================
!                          PHREEQCRM VARIABLES
! ==============================================================================
double precision                      :: init_calcite, na_inflow, mg_inflow,&
                                         ca_inflow, cl_inflow, co2_inflow
double precision, dimension(nipart)   :: den_0, prs_0, tmp_0, sat_0, por_0, vol_0
double precision                      :: cur_time
integer                               :: i, j, m, id, status, ngrd, nsolids
integer                               :: ncomp, nchem, so_col
integer                               :: ic1(nipart, 7)
integer                               :: mask(nipart)
double precision, allocatable         :: bc_conc(:, :), comp_conc(:, :),&
                                         sout(:, :), plot_concs(:, :), plot_times(:)
integer                               :: bc1(1) ! this is for one boundary condition

Pe_grid = (dx * init_v) / init_D
CFL = (init_v * dt) / dx

! print some summary info, to begin with
print *, '====================================================================='
print *, 'Pe_grid =', Pe_grid
print *, 'CFL = ', CFL
write (*, "(' 2 * CFL < Pe_grid < 2 / CFL = ', 3(es9.2, 1x))") 2 * CFL, Pe_grid, 2 / CFL
print *, 'kappa = ', kappa
print *, 'init_v = ', init_v
print *, 'init_D = ', init_D
print *, 'MT stability condition = ', (omega / dble(min(nipart, nptot)))**2 /&
                                        (4.0d0 * kappa * init_D * dt)
print *, 'repvol = ', (omega * 1e3_dp)/ dble(nipart)
print *, '====================================================================='

cur_time = 0.0d0
nsolids = 2 ! calcite and dolomite

! calculate representative volume for immobile particles (1e3 = L/m^3)
repvol = (omega * 1e3_dp)/ dble(nipart)

! make all initial particles active
mparts%active = .false.
nactive = nptot
mparts(1 : nptot)%active = .true.
indices = (/(i, i = 1, nptot)/) ! this is for easy array-indexing of mparts

! scatter the mobile particles randomly throughout domain
call init_random_seed()
call random_number(mparts(1 : nptot)%loc)
mparts(1 : nptot)%loc = (upx - lowx) * mparts(1 : nptot)%loc + lowx
! initialize all mobile particle's bin to -999 for error catching
mparts%bin = -999

! distribute immobile particles evenly
iparts%loc = (/(dble(i) * (omega / (dble(nipart) - 1.0d0)), i = 0, nipart - 1)/)
do i = 1, nipart
    iparts(i)%concs = 0.0d0
enddo

! initial concentration of calcite [mol]
init_calcite = 0.270865d0 ! = 270.865 mol/m^3

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
por_0 = init_porosity ! Porosity
vol_0 = repvol

! print chemistry mask to print detailed output (headings) for only first element
mask = 0
mask(1) = 1

id = RM_Create(nipart, nthreads)
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
status = RM_SetUnitsSSassemblage(id, 1) ! 1 = mol/L of representative volume
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
! get number of grid cells in model (this can change, so it needs to be checked)
ngrd = RM_GetGridCellCount(id)
if (ngrd /= nipart) then
    print *, '****** Cell count differs from particle number ******'
    call exit(status)
endif

nchem = RM_GetChemistryCellCount(id)
if (nchem /= nipart) then
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
! =======================================================================

! the following commented blocks will give you the names of the components and
! headings for selected output.
! They are left here purely for reference, and will not run unless you
! instantiate the corresponding variables.

! this makes a list of the names of components
! allocate(comp_list(ncomp))
! ! print *, 'ncomp = ', ncomp
! do i = 1, ncomp
!     status = RM_GetComponent(id, i, tempname)
!     allocate(character(len_trim(tempname)) :: comp_list(i)%comp)
!     comp_list(i)%comp = trim(tempname)
!     ! print *, comp_list(i)%comp
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
allocate(comp_conc(nipart, ncomp))
status = RM_GetConcentrations(id, comp_conc)

! transfer the concentration array to the immobile particles
call concs_to_iparts(comp_conc, iparts, nipart)

! get the boundary condition
status = RM_InitialPhreeqc2Concentrations(id, bc_conc, 1, bc1)

! this is the injection that corresponds to the FD simulation with the same discretization
bc_conc(1, :) = ((v(1) * dt)/dx) * bc_conc(1, :)

! add the boundary condition to the boundary immobile particle
iparts(1)%concs(:) = iparts(1)%concs(:) + bc_conc(1, :)

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
allocate(sout(nipart, so_col))
status = RM_GetSelectedOutput(id, sout)
allocate(plot_concs(nipart, 3), plot_times(save_steps))

! find out which bins the initial particles are in
allocate (alive(nactive))
alive = (/(i, i = 1, nactive)/)
mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1 ! bin based on velocity grid

do i = 1, nptot
    mparts(i)%concs = 0.0d0
enddo

! convert concentrations on immobile particles to masses
do i = 1, ncomp
    iparts(:)%concs(i) = iparts(:)%concs(i) * repvol
enddo

! get the pairwise distance matrix
allocate(Distmat(nipart, nactive))
call build_Distmat(iparts, mparts, nactive, alive, nipart, D, dt, Distmat)

! do the immobile to mobile transfer
call immobile2mobile(iparts, mparts, nactive, alive, nipart, D, dt, ncomp, Distmat)
deallocate(Distmat, alive)

plot_concs(: , 1) = sout(:, 1)
plot_concs(: , 2) = sout(:, 2) / repvol
plot_concs(: , 3) = sout(:, 4) / repvol
plot_times = (/((i - 1) * save_dt, i = 1, save_steps)/)

! this is the file that stores what comes from plot_concs
! the header gives the shape of the full array and the number of time steps
open (unit=12, file='time_concs.txt', action='write')
write (12, *) shape(plot_concs), save_steps
write (12, *) plot_concs
close (unit=12, status='keep')

! counter for plot_times
j = 2

! time stepping loop
do m = 1, nsteps
    print *, 'step = ', m, ' of ', nsteps

    ! get indices of alive particles
    nactive = count(mparts%active)
    allocate (alive(nactive))
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1 ! bin based on velocity grid

    call advect(mparts, v, dt, alive)
    call diffuse(mparts, nactive, D, dt, alive)

    ! impose reflecting/absorbing boundary conditions
    call reflectlow(mparts, lowx, alive)
    call absorbhigh(mparts, upx, alive, nactive, nexit)

    ! determine what particles are alive after transport
    deallocate (alive)
    nactive = count(mparts%active)
    allocate (alive(nactive))
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1

    ! convert masses on mobile particles to concentrations
    do i = 1, ncomp
        mparts(alive)%concs(i) = mparts(alive)%concs(i) / repvol
    enddo

    allocate(Distmat(nipart, nactive))
    call build_Distmat(iparts, mparts, nactive, alive, nipart, D, dt, Distmat)
    call mobile2immobile(iparts, mparts, nactive, alive, nipart, D, dt, ncomp, Distmat)
    deallocate(Distmat)

    ! transfer the concentrations in the particles array to the PHREEQCRM concs array
    call iparts_to_concs(comp_conc, iparts, nipart)

    cur_time = cur_time + dt
    status = RM_SetTime(id, cur_time)
    status = RM_SetConcentrations(id, comp_conc)
    status = RM_RunCells(id)
    status = RM_GetConcentrations(id, comp_conc)

    call concs_to_iparts(comp_conc, iparts, nipart)

    ! apply "wraparound" boundary condition
    ! move the particles that exited the upper boundary to the lower boundary
    ! and assign them zero mass
    if (m < nsteps) then
        if (nexit > 0) then
            if (nactive + nexit > nptot) then
                print *, '****ERROR: NOT ENOUGH PARTICLES ALLOCATED****'
                call exit(status)
            endif
            nactive = count(.not. mparts%active)
            allocate (nalive(nactive))
            nalive = pack(indices, .not. mparts%active)
            mparts(nalive)%active = .true.
            do i = 1, nactive
                ! this random perturbation keeps kD tree from getting
                    ! tripped up when many new injection particles have the same location
                mparts(nalive(i))%loc = 0.0d0 + rand() * 1.0d-7
                mparts(nalive(i))%concs = 0.0d0
            enddo
            deallocate (nalive)
        endif
    endif

    ! add boundary injection to first immobile particle
    iparts(1)%concs(:) = iparts(1)%concs(:) + bc_conc(1, :)

    deallocate (alive)
    nactive = count(mparts%active)
    allocate (alive(nactive)) ! maybe preallocate to avoid repeatedly doing this
    alive = pack(indices, mparts%active)
    mparts(alive)%bin = floor(mparts(alive)%loc/dxv) + 1

    ! convert concentrations on immobile particles to masses
    do i = 1, ncomp
        iparts(:)%concs(i) = iparts(:)%concs(i) * repvol
    enddo

    allocate(Distmat(nipart, nactive))
    call build_Distmat(iparts, mparts, nactive, alive, nipart, D, dt, Distmat)
    call immobile2mobile(iparts, mparts, nactive, alive, nipart, D, dt, ncomp, Distmat)
    deallocate(Distmat)

    if (cur_time >= plot_times(j)) then
        status = RM_GetSelectedOutput(id, sout)
        plot_concs(: , 1) = sout(:, 1)
        plot_concs(: , 2) = sout(:, 2) / repvol
        plot_concs(: , 3) = sout(:, 4) / repvol
        open (unit=12, file='time_concs.txt', action='write', status='old', access='append')
        write (12, *) plot_concs
        close (unit=12, status='keep')
        j = j + 1
    endif

    deallocate (alive)

enddo

status = RM_Destroy(id)
deallocate (bc_conc, comp_conc, sout, plot_concs, plot_times)

end program RPT_dolo
