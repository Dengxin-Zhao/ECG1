MODULE globvars
!===================================================
!This module contains global variables used in the whole program
!====================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!Numerical constants for accuracy
!===============================

IMPLICIT NONE
INTEGER,PARAMETER ::dp=8           !double precision
REAL(dp),PARAMETER::ZERO=0.E0_dp
REAL(dp),PARAMETER::ONE=1.E0_dp   
REAL(dp),PARAMETER::TWO=2.E0_dp
REAL(dp),PARAMETER::THREE=3.E0_dp
REAL(dp),PARAMETER::FOUR=4.E0_dp
REAL(dp),PARAMETER::FIVE=5.E0_dp
REAL(dp),PARAMETER::SIX=6.E0_dp 
REAL(dp),PARAMETER::SEVEN=7.E0_dp
REAL(dp),PARAMETER::EIGHT=8.E0_dp
REAL(dp),PARAMETER::NINE=9.E0_dp
REAL(dp),PARAMETER::TEN=10.0_dp
REAL(dp),PARAMETER::PI=3.14159265358979323846264338327950E0_dp  
REAL(dp),PARAMETER::SQRTPI=1.77245385090551602729816748334115E0_dp
REAL(dp),PARAMETER::EPS=EPSILON(ONE) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!system information
!======================

!particle system to be calculated
!number of letters must less than 20
CHARACTER(20)::Glob_particle_system

!coordinate case(1:Jacobi or 2:heavy-center)
INTEGER::Glob_coordinate_case

! the total angular momentum
INTEGER::Glob_angular_momentum(2)

!the energy level to be optimized
INTEGER::Glob_energy_level

!integer number for choosing which basis to use
INTEGER::Glob_basis_form

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!basic particle attributes
!=====================

!number of particles
INTEGER::Glob_Nparticle

!number of psudoparticles
INTEGER::Glob_Np

!mass of particles
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_mass

!charge of particles
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_charge

!the type of all particles(random integer numbers)
!the same number corresbonds to the same particle type
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_ptype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!some transformation matrix 
!=====================

!the tranformation matrix U and its inverse 
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_U,Glob_inv_U

!the symmetric matrix Lambda in kinetic term of Hamiltonian
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Lambda

!ri-x_N=wi*x 
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_wi 

!ri-rj=wij*x
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)::Glob_wij

!wi*wi' and  wij*wij'
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)::Glob_Jij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!permutation, spin and parity matrix
!=====================

!numbers of terms of permutation operator
INTEGER::Glob_Nperm

!numbers of spin configurations
INTEGER::Glob_Nspin

!the permutation matrix containing all Glob_Nperm 
!permutaions of particles index in colunms
INTEGER,ALLOCATABLE,DIMENSION(:,:)::Glob_P

!the permutation acting on Jaccobi or heavy particle coordinate
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)::Glob_Tp

!the parity of each permutation in Glob_P
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_permut_parity

!spin configuration(two different integer numbers)
!each number represents one spin direction
INTEGER,ALLOCATABLE,DIMENSION(:,:)::Glob_spin_config

!spin parity for each spin configuration
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_spin_parity

!the symmetry factor corresponding to each permutation
!spin_symmetry* permutation_symmetry
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_symmetry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!nonlinear parameters for all basis
!only storing the lower(uuper)tragular matrix Lk
!=====================

!numbers of nolinear parameters in Lk of each basis
!and matrix storing it
!Glob_NLk=Glob_Np(Glob_Np+1)/2
INTEGER::Glob_NLk
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Lk

!number of mk parameters and matrix storing it
INTEGER::Glob_Nmk
INTEGER,ALLOCATABLE,DIMENSION(:,:)::Glob_mk

!number of initial basis, number of final basis 
INTEGER::Glob_Nbasis_start,Glob_Nbasis_final

!basis increment number per time
INTEGER::Glob_increase_num

!case select global variables for initialize optimization parameters
INTEGER::Glob_init_paras

!overlap threshold
REAL(dp)::Glob_overlap_threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!parameters martix
!Hamilton H and overlap S matrix
!=======================

REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Hkl
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Skl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
!task switch vector
!Glob_task_onoff(i)=0: off
!Glob_task_onoff(i)=1: on
!=========================

!integer vector for switch the main task of the program
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_task_onoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!optimization Global parameters
!========================

!optimization strategy
INTEGER::Glob_opt_strategy 

!limit of total rounds of optimization
INTEGER::Glob_opt_rounds_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!gradient optimization
!========================

!basis number of one batch
INTEGER::Glob_gvm_batch

!the gradient optimization iteration max times
INTEGER::Glob_gvm_ITmax

!the gradient optimization rounds
INTEGER::Glob_gvm_rounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!stochastic optimization
!========================

!optimization times for one parameter 
INTEGER::Glob_svm_IT1

!optimization times for one basis
INTEGER::Glob_svm_IT2

!stochastic optimization rounds
INTEGER::Glob_svm_rounds

!Lk mode and mode parameter 
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_Lk_mode
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Lk_mode_para

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!Global variables required when performing physical output
!=======================

!Rmin-dR-Rmax when calculating correlation function
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_gij_R

!vector to justify  wheather to calculte i-j correlation function
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_gij_onoff

!vector to justify wheather to calculte i-j correlation function
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_rij_onoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!the Global variables used in parallel code (MPI)
!=======================

!set process 0 to be root process
INTEGER,PARAMETER::Glob_root=0

!number of all available process
INTEGER::Glob_numprocs

!MPI error
INTEGER::Glob_MPIerr

!MPI ststus
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_status

!running start time
REAL(dp)::Glob_start_time

!current number of basis number when chage basis
INTEGER::Glob_Nbasis_reach

!the optimized basis reached
INTEGER::Glob_opt_reach

!the energy reached
REAL(dp)::Glob_E_reach

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS


END MODULE globvars
