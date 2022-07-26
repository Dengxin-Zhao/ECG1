MODULE IO
!=================================================
!this module contains subroutine for input and output
!and allocation of basic viriables and matrices
!=================================================
!input.txt:
!  read the input parameters
!log.txt: 
!  record the optimization process
!parafile1.txt and parafile2.txt: 
!  save the optimized parameters in turn
!=================================================
USE MPI
USE globvars
USE symmetry
USE auxfun
USE svmopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!private variavles and matrices used in module IO.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
INTEGER,PRIVATE,PARAMETER::p_mode_paranum_max=5
!!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER,PRIVATE::p_Lk_mode_num
CHARACTER(100),PRIVATE::p_Lk_mode_Char(30)
CHARACTER(100),PRIVATE::p_gij_mode_Char
CHARACTER(100),PRIVATE::p_rij_mode_Char
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Lk_mode
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_Lk_mode_para
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to read from input.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE readfile()
!===============================
!read and allocate Global variables and matrix
!from input.txt
!===============================
IMPLICIT NONE
INTEGER::i,j,k
CHARACTER(100)::preChar(10)

OPEN(1,file='input.txt',action='read')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!particle system 
!======================

!line seperator
READ(1,*)preChar(1)

!name of particle system
!the number of letters must less than 20 
READ(1,*)preChar(1),GLob_particle_system

!system coordinate 
!=======================
!Jacobi and heavy-center
!1:Jacobi
!2:Heavy-center
!=======================
READ(1,*)preChar(1),preChar(2)
IF(preChar(2)(1:6)=='Jacobi')THEN
  Glob_coordinate_case=1
ELSEIF(preChar(2)(1:12)=='Heavy-center')THEN
  Glob_coordinate_case=2
ENDIF

!total angular momentum of the system
READ(1,*)preChar(1),Glob_LM(1:2)

!switch basis form used in the program
!=========================
!Glob_basis_form=0: L=0,M=0 
!Glob_basis_form=1: L=1,M=0
!=========================
IF((Glob_LM(1)==0).AND.(Glob_LM(2)==0))THEN
  Glob_basis_form=0
ELSEIF((Glob_LM(1)==1).AND.(Glob_LM(2)==0))THEN
   Glob_basis_form=1
ENDIF

!energy level to be optimized
READ(1,*)preChar(1),Glob_energy_level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!basic particle attribute
!=====================

!line separator
READ(1,*)preChar(1)

!read the number of particle
READ(1,*)preChar(1),Glob_Nparticle

!define psudoparticles
Glob_Np=Glob_Nparticle-1

!read the masses
ALLOCATE(Glob_mass(Glob_Nparticle))
READ(1,*)preChar(1),Glob_mass(1:Glob_Nparticle)

!read the charges
ALLOCATE(Glob_charge(Glob_Nparticle))
READ(1,*)preChar(1),Glob_charge(1:Glob_Nparticle)

!read the types of  particles
ALLOCATE(Glob_ptype(Glob_Nparticle))
READ(1,*)preChar(1),Glob_ptype(1:Glob_Nparticle)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!allocate transformation matrix associated with masses
!========================

ALLOCATE(Glob_U(Glob_Nparticle,Glob_Nparticle))
ALLOCATE(Glob_inv_U(Glob_Nparticle,Glob_Nparticle))
ALLOCATE(Glob_Lambda(Glob_Nparticle,Glob_Nparticle))
ALLOCATE(Glob_wij(Glob_Np,Glob_Nparticle,Glob_Nparticle))
ALLOCATE(Glob_Jij(Glob_Np,Glob_Np,Glob_Nparticle,Glob_Nparticle))

!form U, inv_U, Lambda matrix
CALL U_Lambda_fun()

!form wi and wij matrix
CALL wij_fun()            

!Jij matrix
CALL Jij_fun()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!spin, permutation and parity
!=======================

!read the numbers of spin configuration
READ(1,*)preChar(1),Glob_Nspin

!read the spin configuration
ALLOCATE(Glob_spin_parity(Glob_Nspin))
ALLOCATE(Glob_spin_config(Glob_Nparticle,Glob_Nspin))
DO i=1,Glob_Nspin
  READ(1,*)Glob_spin_parity(i),preChar(1)(1:1)&
  &,Glob_spin_config(1:Glob_Nparticle,i),preChar(1)(2:2)
ENDDO

!number of all permutations
Glob_Nperm=fac(Glob_Nparticle)  

ALLOCATE(Glob_P(Glob_Nparticle,Glob_Nperm))
ALLOCATE(Glob_Tp(Glob_Np,Glob_Np,Glob_Nperm))
ALLOCATE(Glob_permut_parity(Glob_Nperm))
ALLOCATE(Glob_symmetry(Glob_Nperm))

!generate permutation matrix Glob_P,Glob_Tp
CALL permut_mat_fun()

!form symmetry factor
CALL symmetry_fun()   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!general parameters for basis 
!=======================

!line separator
READ(1,*)preChar(1)

!number of optimization parameters
Glob_NLk=Glob_Np*(Glob_Np+1)/2

!read the initial numbers of basis
READ(1,*)preChar(1),Glob_Nbasis_start

!read the final numbers of basis
READ(1,*)preChar(1),Glob_Nbasis_final

!read the increse number of basis each time
READ(1,*)preChar(1),Glob_increase_num

!read the nonlinear parameter matrix
READ(1,*)preChar(1),Glob_init_paras

!read the overlap threshold
READ(1,*)preChar(1),Glob_overlap_threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
!allocate the private Global variables used in matelem and MPImat
!allocate L H, Smatrix
!=========================

ALLOCATE(Glob_mk(Glob_Nbasis_final))
ALLOCATE(Glob_Lk(Glob_NLk,Glob_Nbasis_final))
ALLOCATE(Glob_Hkl(Glob_Nbasis_final,Glob_Nbasis_final))
ALLOCATE(Glob_Skl(Glob_Nbasis_final,Glob_Nbasis_final))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!allocate task switch matrix
!0: the task is off
!1: the task is on
!==========================

ALLOCATE(Glob_task_onoff(2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!Global optimization parameters
!========================

!line separator
READ(1,*)preChar(1)

!read task on_off
READ(1,*)preChar(1),preChar(2)
IF(preChar(2)(1:2)=='on')THEN
  Glob_task_onoff(1)=1
ELSEIF(preChar(2)(1:3)=='off')THEN
  Glob_task_onoff(1)=0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(1)==0)THEN

DO i=1,10
  READ(1,*)
ENDDO
READ(1,*)preChar(1),k
DO i=1,k
  READ(1,*)
  READ(1,*)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSEIF(Glob_task_onoff(1)==1)THEN

!read the optimization strategy
READ(1,*)preChar(1),Glob_opt_strategy

!read the rounds limit of optimization
READ(1,*)preChar(2),Glob_opt_rounds_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!gvm optimization Global parameters
!========================

!line separator
READ(1,*)preChar(1)

!read the number of basis in one batch
READ(1,*)preChar(1),Glob_gvm_batch

!read the gvm-optimization iteration times 
READ(1,*)preChar(1),Glob_gvm_ITmax

!read the gvm-optimization rounds 
READ(1,*)preChar(1),Glob_gvm_rounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!svm optimization Global parameters
!========================

!line separator
READ(1,*)preChar(1)

!read svm iterations for one parameter
READ(1,*)preChar(1),Glob_svm_IT1

!read svm iterations for one basis 
READ(1,*)preChar(1),Glob_svm_IT2

!read svm-optimization rounds
READ(1,*)preChar(1),Glob_svm_rounds

!=======================
!read Lk mode and mode parameters
!=======================

READ(1,*)preChar(1),p_Lk_mode_num
CALL Lk_mode_read()

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!read physical output: 
!1: Correaltion function
!2: interparticle distance
!========================

!line seperator
READ(1,*)preChar(1)

!read wheather the task of physical output is on
READ(1,*)preChar(1),preChar(2)
IF(preChar(2)(1:2)=='on')THEN
  Glob_task_onoff(2)=1
ELSEIF(preChar(2)(1:3)=='off')THEN
  Glob_task_onoff(2)=0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(2)==0)THEN

DO i=1,6
  READ(1,*)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSEIF(Glob_task_onoff(2)==1)THEN

!=======================
!Correlation function 
!=======================

!line seperator
READ(1,*)preChar(1)

!read Correlation 
CALL gij_mode_read()

ALLOCATE(Glob_gij_R(3))
READ(1,*)preChar(1),Glob_gij_R(1),preChar(2)(1:1),Glob_gij_R(2)&
&,preChar(2)(2:2),Glob_gij_R(3)

!======================
!interparticle distance
!======================

!line seperator
READ(1,*)preChar(1)

!read expectation value rij
CALL rij_mode_read()

!read the power of rij
READ(1,*)preChar(1),Glob_rij_power

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_matelem()
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CLOSE(1)
RETURN
END SUBROUTINE readfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to write information to log.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE writefile()
!===============================
!write the alculation parameters into log.txt
!by Glob_root process
!===============================
IMPLICIT NONE
INTEGER::i,j,k,myid
CHARACTER(100)::st

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

IF(myid==Glob_root)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!CALL write_symmetry_fun()

OPEN(unit=2,file='log.txt')
WRITE(2,*)'=================================================='
WRITE(2,*)'optimization program is running'
WRITE(2,*)'the number of running process is:',Glob_numprocs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'===================system========================='

WRITE(2,*)'particle_system:',' ',GLob_particle_system(1:20)

IF(Glob_coordinate_case==1)THEN
  WRITE(2,*)'system_coordinate:',' ','Jacobi'
ELSEIF(Glob_coordinate_case==2)THEN
  WRITE(2,*)'system_coordinate:',' ','Heavy-center'
ENDIF

WRITE(2,*)'opt_angular_momentum:',GLob_LM(1:2)

WRITE(2,*)'opt_energy_level:',Glob_energy_level

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'================attributes========================='

WRITE(2,*)'particles:',Glob_Nparticle

st="(????????f21.10)"
WRITE(st(2:9),fmt="(TL1,I8.8)")Glob_Nparticle

WRITE(2,*)'masses:'
WRITE(2,st)Glob_mass(1:Glob_Nparticle)

WRITE(2,*)'charges:'
WRITE(2,st)Glob_charge(1:Glob_Nparticle)

st="(????????I8)"
WRITE(st(2:9),fmt="(TL1,I8.8)")Glob_Nparticle

WRITE(2,*)'particle_type:'
WRITE(2,st)Glob_ptype(1:Glob_Nparticle)

st="(f10.5,A5,????????I8,A5)"
WRITE(st(11:18),fmt="(TL1,I8.8)")Glob_Nparticle

WRITE(2,*)'spin_config_num:',Glob_Nspin
DO i=1,Glob_Nspin
WRITE(2,st)Glob_spin_parity(i),'(',Glob_spin_config(1:Glob_Nparticle,i),')'
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'======================basis========================='

WRITE(2,*)'start_basis_size:',Glob_Nbasis_start

WRITE(2,*)'final_basis_size:',Glob_Nbasis_final

WRITE(2,*)'increase_per_time:',Glob_increase_num

WRITE(2,*)'init_nonlinear_para:',Glob_init_paras

WRITE(2,*)'overlap_threshold:',Glob_overlap_threshold

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(1)==1)THEN

WRITE(2,*)'==================opt_para=========================='

WRITE(2,*)'task_on_off:',' ','on'

WRITE(2,*)'opt_strategy:',Glob_opt_strategy

WRITE(2,*)'rounds_limit:',Glob_opt_rounds_limit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'=============================='

WRITE(2,*)'gvm_batch_num:',Glob_gvm_batch

WRITE(2,*)'gvm_opt_ITmax:',Glob_gvm_ITmax

WRITE(2,*)'gvm_opt_rounds:',Glob_gvm_rounds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'==============================='

WRITE(2,*)'svm_opt_IT1:',Glob_svm_IT1

WRITE(2,*)'svm_opt_IT2:',Glob_svm_IT2

WRITE(2,*)'svm_opt_rounds:',Glob_svm_rounds

WRITE(2,*)'Lk_mode_num:',p_Lk_mode_num

st="(????????f10.5)"
WRITE(st(2:9),fmt="(TL1,I8.8)")p_mode_paranum_max

DO i=1,p_Lk_mode_num
WRITE(2,'(I4,A2,A100)')p_Lk_mode(i),'  ',p_Lk_mode_Char(i)  
WRITE(2,st)(p_Lk_mode_para(j,i),j=1,p_mode_paranum_max)
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(2)==1)THEN

WRITE(2,*)'================structure_info===================='

WRITE(2,*)'task_on_off:',' ','on'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'==============================='

WRITE(2,*)'correlation_Gij:'

WRITE(2,'(A17,A2,A100)')'correlation_Gij:','  ',p_gij_mode_Char

WRITE(2,'(A10,f10.5,A5,f10.5,A5,f10.5)')'meshgrid:',Glob_gij_R(1)&
&,':',Glob_gij_R(2),':',Glob_gij_R(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'==============================='

WRITE(2,*)'particle_distance:'

WRITE(2,'(A19,A2,A100)')'particle_distance:','  ',p_rij_mode_Char

WRITE(2,*)'power_of_rij:',Glob_rij_power

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

WRITE(2,*)'=================================================='
WRITE(2,*)
WRITE(2,*)'CALCULATION STARTS:'
WRITE(2,*)

CLOSE(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE writefile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!coordinate transformation matrix 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE  U_Lambda_fun()
!=============================================
!generate the transformation matrix Glob_U and its inverse Glob_inv_U 
!then claculate the Lambda matrix
!=============================================
!input:
!  Glob_mass
!output:
!  Glob_U,Glob_inv_U
!  Glob_Lambda
!=============================================
IMPLICIT NONE
INTEGER::i,j,k
REAL(dp)::mass_sum(Glob_Nparticle)
INTEGER::case_num

!m1,m1+m2,m1+m2+m3,m1+m2+m3+m4,...
mass_sum(1)=Glob_mass(1)
DO i=2,Glob_Nparticle
  mass_sum(i)=mass_sum(i-1)+Glob_mass(i)
ENDDO

SELECT CASE(Glob_coordinate_case)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!select coordinate system
!case(1):Jacobi
!case(2):Heavy-center
!=====================
CASE(1)!jacobi coordinate
    
!U
Glob_U(:,:)=0.0_dp
Glob_U(1,1)=1.0_dp
DO i=1,Glob_Nparticle-1
  Glob_U(i,i+1)=-1.0_dp
ENDDO

DO i=2,Glob_Nparticle
  DO j=1,i
    Glob_U(i,j)=Glob_mass(j)/mass_sum(i)
  ENDDO
ENDDO

!inv_U
Glob_inv_U(:,:)=0.0_dp
DO i=1,Glob_Nparticle
  Glob_inv_U(i,Glob_Nparticle)=1.0_dp
ENDDO

DO i=2,Glob_Nparticle
  Glob_inv_U(i,i-1)=-mass_sum(i-1)/mass_sum(i)
ENDDO

DO i=1,Glob_Nparticle-1
  DO j=i,Glob_Nparticle-1
    Glob_inv_U(i,j)=Glob_mass(j+1)/mass_sum(j+1)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)!heavy particle coordinate
  
!U
Glob_U(:,:)=0.0_dp
DO i=1,Glob_Nparticle-1
  Glob_U(i,i)=1.0_dp
  Glob_U(i,Glob_Nparticle)=-1.0_dp
ENDDO

DO i=1,Glob_Nparticle
  Glob_U(Glob_Nparticle,i)=Glob_mass(i)/mass_sum(Glob_Nparticle)
ENDDO

!inv_U
Glob_inv_U(:,:)=0.0_dp
DO i=1,Glob_Nparticle
  Glob_inv_U(i,i)=1.0_dp
  Glob_inv_U(i,Glob_Nparticle)=1.0_dp
ENDDO

DO i=1,Glob_Nparticle
  DO j=1,Glob_Nparticle-1
    Glob_inv_U(i,j)=Glob_inv_U(i,j)-Glob_mass(j)/mass_sum(Glob_Nparticle)
  ENDDO
ENDDO
    
END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Lambda
DO i=1,Glob_Np
 DO j=1,Glob_Np
   Glob_Lambda(i,j)=0.0_dp
   DO k=1,Glob_Nparticle
   Glob_Lambda(i,j)=Glob_Lambda(i,j)+Glob_U(i,k)*Glob_U(j,k)/Glob_mass(k)
   ENDDO
 ENDDO
ENDDO

!put 1/2 into Lambda
Glob_Lambda(:,:)=Glob_Lambda(:,:)*0.5_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE U_Lambda_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate wi and wij matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE wij_fun()
!=============================================
!calculate the matrix used to from ri-x_(N+1)
!and reletive coordinate ri-rj
!=============================================
!input:
!  Glob_U,Glob_inv_U
!output:
!  wi:  ri-x_(N+1)=wi*x
! wij: ri-rj=wij*x
!=============================================
IMPLICIT NONE
INTEGER::i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Glob_wij(:,:,:)=0.0_dp

!wi
DO i=1,Glob_Nparticle
  DO k=1,Glob_Np
    Glob_wij(k,i,i)=Glob_inv_U(i,k)
  ENDDO
ENDDO

!wij
DO i=1,Glob_Nparticle-1
  DO j=i+1,Glob_Nparticle
    DO k=1,Glob_Np
      Glob_wij(k,i,j)=Glob_inv_U(i,k)-Glob_inv_U(j,k)
    ENDDO
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE wij_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate Jji matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Jij_fun()
!==========================================
!calculate the Jij matrix:
! Jii=wi*wi'; Jij=wij*wij'
!==========================================
!input:
!  Glob_wij
!output:
!  Glob_Jij
!==========================================
IMPLICIT NONE
INTEGER::ii,jj,i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Glob_Jij(:,:,:,:)=0.0_dp

!Jii(i=j)
DO ii=1,Glob_Nparticle
    
  DO i=1,Glob_Np
    DO j=1,Glob_Np
    Glob_Jij(i,j,ii,ii)=Glob_wij(i,ii,ii)*Glob_wij(j,ii,ii)
    ENDDO
  ENDDO
  
ENDDO

!Jij(i<j)
DO ii=1,Glob_Nparticle-1
  DO jj=ii+1,Glob_Nparticle
      
    DO i=1,Glob_Np
      DO j=1,Glob_Np
        Glob_Jij(i,j,ii,jj)=Glob_wij(i,ii,jj)*Glob_wij(j,ii,jj)
      ENDDO
    ENDDO
    
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Jij_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function used for read parameters form input.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_mode_read()
!===============================================
!read the svm mode and corresbonding mode parameters
!of Lk with format <ij-kl~mn>:
!ij-kl represents the the matrix element at (i,j) and
!(k,l) position. kl~mn represents the matrix element 
!at (k,l) and (m,n) position as well as the position 
!between them. (the read format only suits where particle 
!number less than 10)
!===============================================
!output:
!  Glob_Lk_mode and Glob_Lk_mode_para
!===============================================
IMPLICIT NONE
INTEGER::i,j,k,im
INTEGER::index1,index2,index3,index4

ALLOCATE(p_Lk_mode(p_Lk_mode_num))
ALLOCATE(p_Lk_mode_para(p_mode_paranum_max,p_Lk_mode_num))

ALLOCATE(Glob_Lk_mode(Glob_NLk))
ALLOCATE(Glob_Lk_mode_para(p_mode_paranum_max,Glob_NLk))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

loop1:DO im=1,p_Lk_mode_num

READ(1,*)p_Lk_mode(im),p_Lk_mode_Char(im)
READ(1,*)(p_Lk_mode_para(i,im),i=1,p_mode_paranum_max)

loop2:DO i=2,100,3
    
IF(p_Lk_mode_Char(im)(i+2:i+2)=='~')THEN
    
  READ(p_Lk_mode_Char(im)(i:i),*)index1
  READ(p_Lk_mode_Char(im)(i+1:i+1),*)index2 
  IF(index1>=index2)PAUSE'Lk mode allocation error'

  index3=index2-index1
  IF(index1>1)THEN
  DO j=1,index1-1
  index3=index3+(Glob_Nparticle-j)
  ENDDO
  ENDIF
  
  READ(p_Lk_mode_Char(im)(i+3:i+3),*)index1
  READ(p_Lk_mode_Char(im)(i+4:i+4),*)index2
  IF(index1>=index2)PAUSE'Lk mode allocation error'

  index4=index2-index1
  IF(index1>1)THEN
  DO j=1,index1-1
  index4=index4+(Glob_Nparticle-j)
  ENDDO
  ENDIF

  DO j=index3,index4-1
    DO k=1,p_mode_paranum_max
      Glob_Lk_mode_para(k,j)=p_Lk_mode_para(k,im)
    ENDDO
    Glob_Lk_mode(j)=p_Lk_mode(im)
  ENDDO
  
ELSE
  
  READ(p_Lk_mode_Char(im)(i:i),*)index1
  READ(p_Lk_mode_Char(im)(i+1:i+1),*)index2
  IF(index1>=index2)PAUSE'Mode allocation error'

  index3=index2-index1
  IF(index1>1)THEN
  DO j=1,index1-1
  index3=index3+(Glob_Nparticle-j)
  ENDDO
  ENDIF

  DO k=1,p_mode_paranum_max
    Glob_Lk_mode_para(k,index3)=p_Lk_mode_para(k,im)
  ENDDO
  Glob_Lk_mode(index3)=p_Lk_mode(im)
  
  IF(p_Lk_mode_Char(im)(i+2:i+2)=='>')THEN
    EXIT loop2
  ELSEIF(p_Lk_mode_Char(im)(i+2:i+2)=='-')THEN
    CYCLE loop2
  ELSE
    PAUSE'Lk mode allocation error'
  ENDIF
    
ENDIF

ENDDO loop2

ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Lk_mode_read 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gij_mode_read()
!==============================================
!read the calculation parameters of correlation function gij
!==============================================
IMPLICIT NONE
INTEGER::i,j,k,L
INTEGER::index1,index2,index3,index4
CHARACTER(100)::preChar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

READ(1,*)preChar,p_gij_mode_Char

ALLOCATE(Glob_gij_onoff(Glob_Nparticle*(Glob_Nparticle+1)/2))
DO i=1,Glob_Nparticle*(Glob_Nparticle+1)/2
  Glob_gij_onoff(i)=0
ENDDO

loop:DO i=2,100,3

IF(p_gij_mode_Char(i+2:i+2)=='~')THEN
    
READ(p_gij_mode_Char(i:i),*)index1
READ(p_gij_mode_Char(i+1:i+1),*)index2 
IF(index1>index2)PAUSE'gij mode allocation error'
  
index3=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
    index3=index3+(Glob_Nparticle-j+1)
  ENDDO
ENDIF
    
READ(p_gij_mode_Char(i+3:i+3),*)index1
READ(p_gij_mode_Char(i+4:i+4),*)index2
IF(index1>index2)PAUSE'gij mode allocation error'
    
index4=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
    index4=index4+(Glob_Nparticle-j+1)
  ENDDO
ENDIF

DO j=index3,index4-1
  Glob_gij_onoff(j)=1
ENDDO
    
ELSE
    
READ(p_gij_mode_Char(i:i),*)index1
READ(p_gij_mode_Char(i+1:i+1),*)index2
IF(index1>index2)PAUSE'gij mode allocation error'

index3=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
  index3=index3+(Glob_Nparticle-j+1)
  ENDDO
ENDIF
  
Glob_gij_onoff(index3)=1
    
IF(p_gij_mode_Char(i+2:i+2)=='>')THEN
  EXIT loop
ELSEIF(p_gij_mode_Char(i+2:i+2)=='-')THEN
  CYCLE loop
ELSE
  PAUSE'gij mode allocation error'
ENDIF

ENDIF

ENDDO loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gij_mode_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_mode_read()
!=========================================
!read the calculation parameters of expectation value of rij
!=========================================
IMPLICIT NONE
INTEGER::i,j,k,L
INTEGER::index1,index2,index3,index4
CHARACTER(100)::preChar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

READ(1,*)preChar,p_rij_mode_Char

ALLOCATE(Glob_rij_onoff(Glob_Nparticle*(Glob_Nparticle+1)/2))
DO i=1,Glob_Nparticle*(Glob_Nparticle+1)/2
  Glob_rij_onoff(i)=0
ENDDO

loop:DO i=2,100,3

IF(p_rij_mode_Char(i+2:i+2)=='~')THEN
    
READ(p_rij_mode_Char(i:i),*)index1
READ(p_rij_mode_Char(i+1:i+1),*)index2 
IF(index1>index2)PAUSE'rij mode allocation error'
  
index3=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
    index3=index3+(Glob_Nparticle-j+1)
  ENDDO
ENDIF
    
READ(p_rij_mode_Char(i+3:i+3),*)index1
READ(p_rij_mode_Char(i+4:i+4),*)index2
IF(index1>index2)PAUSE'rij mode allocation error'
    
index4=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
    index4=index4+(Glob_Nparticle-j+1)
  ENDDO
ENDIF

DO j=index3,index4-1
  Glob_rij_onoff(j)=1
ENDDO
    
ELSE
    
READ(p_rij_mode_Char(i:i),*)index1
READ(p_rij_mode_Char(i+1:i+1),*)index2
IF(index1>index2)PAUSE'rij mode allocation error'

index3=index2-index1+1
IF(index1>1)THEN
  DO j=1,index1-1
  index3=index3+(Glob_Nparticle-j+1)
  ENDDO
ENDIF
  
Glob_rij_onoff(index3)=1
    
IF(p_rij_mode_Char(i+2:i+2)=='>')THEN
  EXIT loop
ELSEIF(p_rij_mode_Char(i+2:i+2)=='-')THEN
  CYCLE loop
ELSE
  PAUSE'rij mode allocation error'
ENDIF

ENDIF

ENDDO loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN  
END SUBROUTINE rij_mode_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!initial optimization parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_paras()
!=============================================
!generate or read the nolinear parameters 
!Glob_init_paras=0: generate new basis by basis_increase
!Glob_init_paras=1: read from parafile1.txt 
!Glob_init_paras=2: read from parafile2.txt
!=============================================
IMPLICIT NONE
INTEGER::i,j,k,myid
INTEGER::sendnum

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

SELECT CASE(Glob_init_paras)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)!generate random parameters 

Glob_Nbasis_reach=0
Glob_opt_reach=0
CALL basis_increase(0,Glob_Nbasis_start,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)!read from parafile1.txt

IF(myid==Glob_root)THEN     
 CALL read_parafile1()
ENDIF    
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_BCAST(Glob_Nbasis_reach,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(Glob_opt_reach,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

IF(Glob_basis_form==1)THEN
  sendnum=Glob_Nbasis_final
  CALL MPI_BCAST(Glob_mk,sendnum,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
ENDIF

sendnum=Glob_NLk*Glob_Nbasis_final
CALL MPI_BCAST(Glob_Lk,sendnum,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)!read from parafile2.txt

IF(myid==Glob_root)THEN     
  CALL read_parafile2()
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_BCAST(Glob_Nbasis_reach,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(Glob_opt_reach,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

IF(Glob_basis_form==1)THEN
  sendnum=Glob_Nbasis_final
  CALL MPI_BCAST(Glob_mk,sendnum,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
ENDIF

sendnum=Glob_NLk*Glob_Nbasis_final
CALL MPI_BCAST(Glob_Lk,sendnum,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE init_paras

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!read the parameters from parafile1.txt and parafile2.txt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_parafile1()
!================================================
!read parameters from parafile1.txt
!================================================
IMPLICIT NONE
INTEGER::i,j,k

CHARACTER(100)::preChar(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(unit=3,file='parafile1.txt')

READ(3,*)preChar(1)

READ(3,*)preChar(1),preChar(2)

READ(3,*)preChar(1),preChar(2)

READ(3,*)preChar(1),preChar(2)

READ(3,*)preChar(1),Glob_Nbasis_reach
IF(Glob_Nbasis_reach/=Glob_Nbasis_start)THEN
  PAUSE'basis number in parafile1.txt differs from that in input.txt'
ENDIF

READ(3,*)preChar(1),Glob_opt_reach

READ(3,*)preChar(1),Glob_E_reach

READ(3,*)preChar(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_basis_form==0)THEN

  DO i=1,Glob_Nbasis_reach
    READ(3,*)Glob_Lk(1:Glob_NLk,i)
  ENDDO

ELSEIF(Glob_basis_form==1)THEN

  DO i=1,Glob_Nbasis_reach
    READ(3,*)Glob_mk(i),Glob_Lk(1:Glob_NLk,i)
  ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CLOSE(3)
RETURN
END SUBROUTINE read_parafile1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_parafile2()
!================================================
!read parameters from parafile.txt
!================================================
INTEGER::i,j,k

CHARACTER(100)::preChar(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(unit=4,file='parafile2.txt')

READ(4,*)prechar(1)

READ(4,*)preChar(1),preChar(2)

READ(4,*)preChar(1),preChar(2)

READ(4,*)preChar(1),preChar(2)

READ(4,*)preChar(1),Glob_Nbasis_reach
IF(Glob_Nbasis_reach/=Glob_Nbasis_start)THEN
  PAUSE'basis number in parafile2.txt differs from that in input.txt'
ENDIF

READ(4,*)preChar(1),Glob_opt_reach

READ(4,*)preChar(1),Glob_E_reach

READ(4,*)prechar(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_basis_form==0)THEN

  DO i=1,Glob_Nbasis_reach
    READ(4,*)Glob_Lk(1:Glob_NLk,i)
  ENDDO

ELSEIF(Glob_basis_form==1)THEN

  DO i=1,Glob_Nbasis_reach
    READ(4,*)Glob_mk(i),Glob_Lk(1:Glob_NLk,i)
  ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CLOSE(4)
RETURN
END SUBROUTINE read_parafile2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE IO
