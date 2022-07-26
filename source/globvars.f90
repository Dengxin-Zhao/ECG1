MODULE globvars
!===================================================
!This module contains global variables used in the whole program
!====================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!Numerical constants for accuracy
!===============================

IMPLICIT NONE
INTEGER,PARAMETER ::dp=8                !precision
REAL(dp),PARAMETER::EPS=EPSILON(1.0_dp) !min value in precision 
REAL(dp),PARAMETER::HUG=HUGE(1.0_dp)    !max value
REAL(dp),PARAMETER::TNY=TINY(1.0_dp)    !min value
REAL(dp),PARAMETER::PI=3.14159265358979323846264338327950_dp  
REAL(dp),PARAMETER::SQRTPI=1.77245385090551602729816748334115_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!system information
!========================

!particle system to be calculated (number of letters must less than 50)
CHARACTER(50)::Glob_particle_system

!coordinate case(1:Jacobi or 2:heavy-center)
INTEGER::Glob_coordinate_case

!angular momentum and magnetic quantum number
INTEGER::Glob_LM(2)

!the energy level to be optimized
INTEGER::Glob_energy_level

!integer number for choosing which basis to use
INTEGER::Glob_basis_form

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!basic particle attributes
!========================

!number of particles
INTEGER::Glob_Nparticle

!number of psudoparticles
INTEGER::Glob_Np

!mass of particles
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_mass

!charge of particles
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_charge

!the type of all particles
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_ptype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!some transformation matrix 
!========================

!the tranformation matrix U and its inverse 
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_U,Glob_inv_U

!the symmetric matrix Lambda in kinetic term of Hamiltonian
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Lambda

!ri-x_N=wi*x and ri-rj=wij*x
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)::Glob_wij

!wi*wi' and  wij*wij'
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:,:)::Glob_Jij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!permutation, spin and parity matrix
!========================

!numbers of terms of permutation operator
INTEGER::Glob_Nperm

!numbers of spin configurations
INTEGER::Glob_Nspin

!the permutation matrix containing all Nperm 
!permutaions of particles index in colunms
INTEGER,ALLOCATABLE,DIMENSION(:,:)::Glob_P

!permutation acting on Jaccobi or heavy particle coordinate
REAL(dp),ALLOCATABLE,DIMENSION(:,:,:)::Glob_Tp

!parity of each permutation in Glob_P
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_permut_parity

!spin configuration
INTEGER,ALLOCATABLE,DIMENSION(:,:)::Glob_spin_config

!spin parity for each spin configuration
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_spin_parity

!symmetry factor corresponding to each permutation
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_symmetry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!basis parameters
!========================

!nolinear parameters Lk 
INTEGER::Glob_NLk
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Lk

!parameters mk
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_mk

!number of initial basis and final basis 
INTEGER::Glob_Nbasis_start,Glob_Nbasis_final

!basis increment number per time
INTEGER::Glob_increase_num

!case select global variables for initialize optimization parameters
INTEGER::Glob_init_paras

!overlap threshold
REAL(dp)::Glob_overlap_threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!Hamilton H and overlap S matrix
!========================

REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Hkl
REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Glob_Skl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!task switch vector
!Glob_task_onoff(i)=0: off
!Glob_task_onoff(i)=1: on
!========================

!integer vector for switch the task of the program
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_task_onoff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!optimization parameters
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!Global variables required when performing physical output
!========================

!Rmin-dR-Rmax when calculating correlation function
REAL(dp),ALLOCATABLE,DIMENSION(:)::Glob_gij_R

!vector to justify  wheather to calculte i-j correlation function
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_gij_onoff

!vector to justify wheather to calculte i-j correlation function
INTEGER,ALLOCATABLE,DIMENSION(:)::Glob_rij_onoff

!power of interparticle distance
REAL(dp)::Glob_rij_power

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!the Global variables used in parallel code (MPI)
!========================

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE globvars
