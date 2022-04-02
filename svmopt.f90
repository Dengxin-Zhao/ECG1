MODULE svmopt
!=================================
!this module contains subroutines for 
!stochastic optimization 
!=================================
USE MPI
USE globvars
USE auxfun
USE eigen
USE MPImat
USE gvmopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!privete parameters used in svmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!seeds number for generating random number between [0,1]
!=========================
INTEGER,PRIVATE::Lkseed1=-60067
INTEGER,PRIVATE::Lkseed2=-41077
INTEGER,PRIVATE::mkseed=-50031
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!generation loop limits
!==========================
INTEGER,PRIVATE,PARAMETER::generate_limit=10
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to increase the basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE basis_increase(Nb_pre,Nb_after,gvm_IT1,gvm_IT2,svm_IT1,svm_IT2)
!===============================
!this subroutine increase the basis size by choose
!the new stochastic parameters lowering the energy
!====================
!strategy: 
!generate sets of parameters randomly one by one chose the one that
!lower the energy, then optimize the new basis gvm_IT1 times with BFGs
!meanwhilde optimize the new basis with stochastic approach  
!generate a new basis accordingly, optimize all new basis
!from Nb_pre->Nb_now ITmax times with batch BFGs
!====================
!input:
!  Nb_pre: the basis number before increase
!  Nb_after: the basis number after increase
! gvm_IT1: optimization times for one increased basis with GVM 
! gvm_IT2: optimization times of all increaed dN basis with GVM 
!          after increase the basis from Nb_pre to Nb_now
! svm_IT1: optimization times of each parameter for one increased basis with SVM 
! svm_IT2: optimization times of each basis for one increased basis with SVM
!===============================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nb_pre,Nb_after
INTEGER,INTENT(IN)::gvm_IT1,gvm_IT2
INTEGER,INTENT(IN)::svm_IT1,svm_IT2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER,PARAMETER::Nb_init=60   !initial basis size
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
REAL(dp)::Eval,Ci(Nb_after)

INTEGER::i,j,k,L,myid,ERR
INTEGER::nb,Nb_start,count,write_index
INTEGER::mk,bk_cal(2)

REAL(dp)::randnum,alpha(Glob_NLk),vechL(Glob_NLk)

CHARACTER(14)::time_st

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

IF(Nb_pre>=Nb_init)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!IF(Nb_pre>=Nb_init),diagonalize the Nb_pre basis once to 
!obtain the energy to start with
!===============================

Nb_start=Nb_pre
    
ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!IF(Nb_pre<Nb_init),increase the basis size to the amount
!able to be diagonalized: Nb_init
!===============================
    
Nb_start=Nb_init

ERR=1
DO WHILE(ERR/=0)
ERR=0

init_loop:DO nb=Nb_pre+1,Nb_init

CALL pvars_MPImat(nb,1) 

!==========================
!generate basis
!==========================
ERR=1
generate_loop:DO WHILE(ERR/=0)
ERR=0

IF(myid==Glob_root)THEN
    
!generate mk
IF(Glob_basis_form==1)THEN
  randnum=randnum_fun(mkseed)
  mk=1+floor(randnum*Glob_Np)
ENDIF

!generate alpha
DO i=1,Glob_NLk
  randnum=ZERO
  CALL Lk_mode_randnum(i,Glob_Lk_mode(i),randnum)
  alpha(i)=randnum
ENDDO

!transform alpha to L
CALL alpha_to_L(alpha,vechL,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop

IF(Glob_basis_form==1)THEN
  CALL MPI_BCAST(mk,Glob_Nmk,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
  Glob_mk(1,nb)=mk
ENDIF

CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,nb)=vechL(:)

!check its overlap with other basis 
CALL overlap_check_MPIfun(nb,nb,ERR)

ENDDO generate_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(nb,0) 

ENDDO init_loop

CALL pvars_MPImat(Nb_start,1) 
bk_cal(1)=1
bk_cal(2)=Nb_start
CALL E_rows_MPIfun(Nb_start,bk_cal,Nb_start,Eval,Ci(1:Nb_start),ERR,1)
CALL pvars_MPImat(Nb_start,0)

ENDDO 
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nb_start,1) 
CALL E_rows_MPIfun(Nb_start,bk_cal,Nb_start,Eval,Ci(1:Nb_start),ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before generating new basis'
CALL pvars_MPImat(Nb_start,0) 

IF(myid==Glob_root)THEN 
  IF(Nb_pre==0)THEN
    OPEN(unit=2,file='log.txt',position='append')
    WRITE(2,*)'=================================================' 
    WRITE(2,*)'initalize basis to:',Glob_Nbasis_start
    WRITE(2,*)'================================================='
  ELSE
    OPEN(unit=2,file='log.txt',position='append')
    WRITE(2,*)'=================================================' 
    WRITE(2,'(A20,I8,A6,I8)')'increase basis from',Nb_pre,'to',Nb_after
    WRITE(2,*)'================================================='
  ENDIF
ENDIF

IF(myid==Glob_root)THEN 
  WRITE(2,*)'start energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,A10,A14,A11)')'time','count','basis_num','energy'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!if Nb_pre<Nb_init then optimize the initial basis 
!========================
count=0
write_index=0

IF(Nb_pre<Nb_init)THEN
    
CALL pvars_MPImat(Nb_start,1) 
CALL pvars_eigen(Nb_start,1)
CALL pvars_gvmopt(1,Glob_NLk,1)

basis_loop1:DO nb=Nb_pre+1,Nb_start
count=count+1
bk_cal(1)=nb
bk_cal(2)=nb

CALL one_diag_init_MPIfun(nb,Nb_start)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(nb,Nb_start,Eval,Ci(1:Nb_start))
ENDIF

CALL Lk_gvm_MPIfun(1,bk_cal,Nb_start,gvm_IT1,Eval,Ci(1:Nb_start)) 

CALL Lk_svm_MPIfun(nb,Nb_start,svm_IT1,svm_IT2,Eval,Ci(1:Nb_start))  

Glob_Nbasis_reach=Nb_start
Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)

!write present energy in log file
WRITE(2,'(A14,I10,I14,f28.16)')time_st,count,nb,Eval

!save in parafile1.txt and parafile2.txt
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) !write in parafile1.txt 
ELSE
  CALL write_parafile(2) !write in parafile2.txt 
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO basis_loop1

CALL pvars_gvmopt(1,Glob_NLk,0)
CALL pvars_gvmopt(Nb_start,Glob_NLk,1)

bk_cal(1)=1
bk_cal(2)=Nb_start
CALL Lk_gvm_MPIfun(Nb_start,bk_cal,Nb_start,gvm_IT2,Eval,Ci(1:Nb_start))

CALL pvars_gvmopt(Nb_start,Glob_NLk,0)
CALL pvars_MPImat(Nb_start,0) 
CALL pvars_eigen(Nb_start,0)

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!====================
!increse basis and optimize it
!====================

IF(Nb_after>Nb_init)THEN
    
CALL pvars_gvmopt(1,Glob_NLk,1)

basis_loop2:DO nb=Nb_start+1,Nb_after 
count=count+1

CALL pvars_MPImat(nb,1)  
CALL pvars_eigen(nb,1)
CALL one_diag_init_MPIfun(nb,nb)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL basis_generate_MPIfun(nb,nb,Eval)

bk_cal(1)=nb
bk_cal(2)=nb
CALL E_rows_MPIfun(1,bk_cal,nb,Eval,Ci(1:nb),ERR,1)

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(nb,nb,Eval,Ci(1:nb))
ENDIF

CALL Lk_gvm_MPIfun(1,bk_cal,nb,gvm_IT1,Eval,Ci(1:nb))

CALL Lk_svm_MPIfun(nb,nb,svm_IT1,svm_IT2,Eval,Ci(1:nb))


Glob_Nbasis_reach=nb
Glob_opt_reach=nb
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)

!write present energy in log file
WRITE(2,'(A14,I10,I14,f28.16)')time_st,count,nb,Eval

!save in parafile1.txt and parafile2.txt
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) !write in parafile1.txt 
ELSE
  CALL write_parafile(2) !write in parafile2.txt 
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(nb,0) 
CALL pvars_eigen(nb,0)

ENDDO basis_loop2

CALL pvars_gvmopt(1,Glob_NLk,0)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!optimzize all generated basis ITmax2 times
IF(gvm_IT2/=0)THEN
  IF(myid==Glob_root)THEN
    WRITE(2,*)'================================================='
    WRITE(2,*)'start optimizing all new basis with BFGs, iterations:',gvm_IT2
  ENDIF
  
  bk_cal(1)=Nb_pre+1
  bk_cal(2)=Nb_after
  
  CALL pvars_MPImat(Nb_after,1) 
  CALL pvars_gvmopt(Nb_after-Nb_pre,Glob_NLk,1)
  CALL Lk_gvm_MPIfun(Nb_after-Nb_pre,bk_cal,Nb_after,gvm_IT2,Eval,Ci(1:Nb_after))
  CALL pvars_MPImat(Nb_after,0) 
  CALL pvars_gvmopt(Nb_after-Nb_pre,Glob_NLk,0)
  
  IF(myid==Glob_root)THEN
    WRITE(2,*)'energy reached:',Eval
  ENDIF
ENDIF

Glob_E_reach=Eval

IF(myid==Glob_root)THEN 
IF(Nb_pre==0)THEN
  WRITE(2,*)'================================================='
  WRITE(2,*)'basis initialization finished, basis reached:',Nb_after
  WRITE(2,*)'=================================================' 
ELSE
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'basis increment finished, basis reached:',Nb_after
  WRITE(2,*)'=================================================' 
ENDIF
  WRITE(2,*)
ENDIF
CLOSE(2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE basis_increase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to increase basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE basis_generate_MPIfun(bk,Nbasis,Eval)
!================================
!generate one new basis 
!stategy:
!generate sets of parameters randomly one
!by one chose the one that lower the energy
!input:
!    bk: basis index to generate
!Nbasis: basis number after generating
!inout:
!  Eval: eigen value
!output:
!  new parameters stored in Glob_mk and Glob_Lk
!================================
!this subroutine should
!CALL pvars_eigen(Nbasis,1)
!before calling it 
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,Nbasis
REAL(dp),INTENT(INOUT)::Eval

INTEGER::i,j,k,L,myid,ERR
INTEGER::mk,bk_cal(2)

REAL(dp)::randnum
REAL(dp)::E_temp,C_temp(Nbasis)
REAL(dp)::alpha(Glob_NLk),vechL(Glob_NLk)

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

bk_cal(1)=bk
bk_cal(2)=bk

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ERR=1
generate_loop:DO WHILE(ERR/=0)
ERR=0

IF(myid==Glob_root)THEN

!generate mk
IF(Glob_basis_form==1)THEN
  randnum=randnum_fun(mkseed)
  mk=1+floor(randnum*Glob_Np)
ENDIF

!generate alpha
DO i=1,Glob_NLk
  randnum=ZERO
  CALL Lk_mode_randnum(i,Glob_Lk_mode(i),randnum)
  alpha(i)=randnum
ENDDO

!transform alpha to L
CALL alpha_to_L(alpha,vechL,ERR) 

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop

IF(Glob_basis_form==1)THEN
  CALL MPI_BCAST(mk,Glob_Nmk,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
  Glob_mk(1,bk)=mk
ENDIF

CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechL(:)

!overlap check
CALL overlap_check_MPIfun(bk,Nbasis,ERR)
IF(ERR/=0)CYCLE generate_loop

!check whether energy is lower
CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_temp,C_temp,ERR,0)
IF(E_temp>=Eval)THEN
  ERR=1
  CYCLE generate_loop
ENDIF

ENDDO generate_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Eval=E_temp

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE basis_generate_MPIfun 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!stockastic optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svm_opt(Nbasis,IT1,IT2,Nrounds)
!================================
!this subroutine optimizes the nonlinear parameters
!with stochastic approach
!====================
!strategy: 
!first generate random numbers of given modes as an new
!nolinear parameter of the basis to form a new A matrix, second
!check whether it's positive definite, then decompose it into 
!Lmatrix and check whether its overlap with others is less than
!a given threshold.If both requirement is fufilled, calculate
!new eigenvalues. If the new eigenvalue is better than the 
!best one ever found,the new parameter is accepted.
!==================== 
!input:
!Nbasis: number of basis
!   IT1: optimization times for one parameter in the optimized basis
!   IT2: optimization times for one basis
!Nrounds: perform svm optimization Nround for all basis
!=================================
INTEGER,INTENT(IN)::Nbasis
INTEGER,INTENT(IN)::IT1,IT2
INTEGER,INTENT(IN)::Nrounds

INTEGER::i,j,k,myid,ERR,index
INTEGER::nr,nb,bk,nb_start
INTEGER::bk_cal(2)
INTEGER::count,write_index,generate_index

REAL(dp)::Eval,Ci(Nbasis)

CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,1)
CALL pvars_eigen(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================
!full diagonalization once first
!and initialize Glob_S and Glob_H matrix
!================================

bk_cal(1)=1
bk_cal(2)=Nbasis
CALL E_rows_MPIfun(Nbasis,bk_cal,Nbasis,Eval,Ci,ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before all_svm'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'stochastic optimization start:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'================================================='
  WRITE(2,*)'start energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,2A10,A14,A11)')'time','count','round','basis','energy'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!stochastic optimization
!=============================

count=0
write_index=0
optimize_loop:DO nr=1,Nrounds

nb_start=1
IF(nr==1)THEN
  nb_start=Glob_opt_reach
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

basis_loop:DO nb=nb_start,Nbasis
count=count+1

CALL one_diag_init_MPIfun(nb,Nbasis)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
  
!==========================
!optimization of mk
!==========================
IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(nb,Nbasis,Eval,Ci)
ENDIF

!===========================
!optimization of Lk with BFGs algrithm
!===========================
CALL Lk_svm_MPIfun(nb,Nbasis,IT1,IT2,Eval,Ci)

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)

!write present energy in log file
WRITE(2,'(A14,2I10,I14,f28.16)')time_st,count,nr,nb,Eval

!save in parafile1.txt and parafile2.txt
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) !write in parafile1.txt 
ELSE
  CALL write_parafile(2) !write in parafile2.txt 
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO basis_loop
ENDDO optimize_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN  !start this optimization
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'stochastic optimization finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,0)
CALL pvars_eigen(Nbasis,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE svm_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gvm_svm_opt(Nbs,Nbasis,gvm_ITmax,svm_IT1,svm_IT2,Nrounds)
!===================================
!optimize each batch of basis with SVM and GVM at same time
!input:
!    Nbasis: the present total basis number
!Nbasis: present total basis number
! gvm_ITmax: optimization times for one increased basis with GVM 
! svm_IT1: optimization times of each parameter for one increased basis with SVM 
! svm_IT2: optimization times of each basis for one increased basis with SVM
! Nrounds: total Nrounds of optimization
!===================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,Nbasis
INTEGER,INTENT(IN)::gvm_ITmax,svm_IT1,svm_IT2
INTEGER,INTENT(IN)::Nrounds

INTEGER::i,j,k,myid,ERR
INTEGER::nr,nb,bk,count,nb_start
INTEGER::Nbatches,bk_cal(2)
INTEGER,ALLOCATABLE,DIMENSION(:,:)::basis_index

REAL(dp)::Eval,Ci(Nbasis)

INTEGER::write_index
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================
!initialization
!=================================

CALL pvars_MPImat(Nbasis,1)
CALL pvars_eigen(Nbasis,1)
CALL pvars_gvmopt(Nbs,Glob_NLk,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

!total batches
k=mod(Nbasis,Nbs)
Nbatches=1+(Nbasis-k)/Nbs
IF(k==0)THEN
  Nbatches=Nbatches-1
ENDIF

!basis index of each batch
ALLOCATE(basis_index(2,Nbatches))
basis_index(1,Nbatches)=1+(Nbasis-Nbs)
basis_index(2,Nbatches)=Nbasis
IF(Nbatches>1)THEN
DO i=1,Nbatches-1
  j=1+(i-1)*Nbs
  k=j+Nbs-1
  basis_index(1,i)=j
  basis_index(2,i)=k
ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!====================================
!full diagonalization once first
!and initialize Glob_Skl and Glob_Hkl matrix
!====================================
bk_cal(1)=1
bk_cal(2)=Nbasis
CALL E_rows_MPIfun(Nbasis,bk_cal,Nbasis,Eval,Ci,ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before all_gvm_svm'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'gvm+svm optimization start:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'================================================='
  WRITE(2,*)'start energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,2A10,A14,A11)')'time','count','round','basis','energy'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!gvm-svm optimizations
!==============================

count=0
write_index=0
optimize_loop:DO nr=1,Nrounds

nb_start=1
IF(nr==1)THEN
  k=mod(Glob_opt_reach,Nbs)
  nb_start=1+(Glob_opt_reach-k)/Nbs
  IF(k==0)THEN
    nb_start=nb_start-1
  ENDIF
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!optimize each basis batch 
!===============================

batch_loop:DO nb=nb_start,Nbatches 
count=count+1
bk_cal(:)=basis_index(:,nb)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!gvm optimization
!===============================

IF(Glob_basis_form==1)THEN
  CALL mk_gvm_MPIfun(Nbs,bk_cal,Nbasis,Eval,Ci)
ENDIF

CALL Lk_gvm_MPIfun(Nbs,bk_cal,Nbasis,gvm_ITmax,Eval,Ci)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!svm optimization
!===============================

DO bk=bk_cal(1),bk_cal(2) 

CALL one_diag_init_MPIfun(bk,Nbasis)

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(bk,Nbasis,Eval,Ci)
ENDIF

CALL Lk_svm_MPIfun(bk,Nbasis,svm_IT1,svm_IT2,Eval,Ci)

ENDDO 

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)

!write present energy in log file
WRITE(2,'(A14,2I10,I14,f28.16)')time_st,count,nr,nb,Eval

!save in savefile1.txt and parafile.txt
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) !write in parafile1.txt 
ELSE
  CALL write_parafile(2) !write in parafile2.txt 
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO batch_loop
ENDDO optimize_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN  !start this optimization
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'stochastic optimization finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,0)
CALL pvars_eigen(Nbasis,0)
CALL pvars_gvmopt(Nbs,Glob_NLk,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE gvm_svm_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_svm_MPIfun(bk,Nbasis,IT1,IT2,Eval,Ci)  
!===================================
!optimize one basis with SVM approach 
!input:
!    bk: the basis to be optimized
!Nbasis: present total basis number
!   IT1: optimization times for one parameter in the optimized basis
!   IT2: optimization times for one basis
!output: 
!  Eval: eigen value
!    Ci: eigen vector
!======================
!this subroutine should
!CALL pvars_MPImat(Nbasis,1)
!before using it
!===================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,Nbasis
INTEGER,INTENT(IN)::IT1,IT2
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)

INTEGER::i,j,k,myid,ERR
INTEGER::t1,t2,ip,generate_index,bk_cal(2)

REAL(dp)::randnum
REAL(dp)::E_best,C_best(Nbasis)
REAL(dp)::E_temp,C_temp(Nbasis)
REAL(dp)::alpha_best(Glob_NLk),L_best(Glob_NLk)
REAL(dp)::L_old(Glob_NLk)
REAL(dp)::alpha(Glob_NLk),vechL(Glob_NLk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(IT2==0)RETURN
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

bk_cal(1)=bk
bk_cal(2)=bk

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==================================
!stochastic optimization 
!==================================

E_best=Eval
L_best(:)=Glob_Lk(:,bk)
L_old(:)=L_best(:)
CALL L_to_alpha(L_best,alpha_best)

!optimize nb_th basis IT2 times
IT2_loop:DO t2=1,IT2

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!IT1==0 : do not optimize each parameter seperately 
!=======================
IF(IT1==0)THEN
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!generate nolinear parameter Lk
!===========================

generate_index=0
generate_loop1:DO WHILE(generate_index<generate_limit)
generate_index=generate_index+1
ERR=0

IF(myid==Glob_root)THEN
    
!generate alpha
DO i=1,Glob_NLk
  randnum=alpha_best(i)
  CALL Lk_mode_randnum(i,Glob_Lk_mode(i),randnum)
  alpha(i)=randnum
ENDDO

!transform alpha to L
CALL alpha_to_L(alpha,vechL,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop1

!check overlap 
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechL(:)
CALL overlap_check_MPIfun(bk,Nbasis,ERR)
IF(ERR==0)EXIT generate_loop1

ENDDO generate_loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(generate_index<generate_limit)THEN

CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_temp,C_temp,ERR,0)
IF(E_temp<E_best.AND.ERR==0)THEN
  E_best=E_temp
  alpha_best(:)=alpha(:)
  L_best(:)=vechL(:)
ELSE
  alpha(:)=alpha_best(:)
ENDIF

ELSE
    
  alpha(:)=alpha_best(:)
  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!IT1/=0 : optimize each parameter seperately IT1 timies 
!=======================
ELSE

alpha(:)=alpha_best(:)

!generate loop
ip_loop:DO ip=1,Glob_NLk  

! iterations for optimizing each parameter
IT1_loop:DO t1=1,IT1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
! generate nonlinear parameter
!==============================

generate_index=0
generate_loop2:DO WHILE(generate_index<generate_limit)
generate_index=generate_index+1
ERR=0

IF(myid==Glob_root)THEN
    
!generate new parameter vechA(ip)
randnum=alpha_best(ip)
CALL Lk_mode_randnum(ip,Glob_Lk_mode(ip),randnum)
alpha(ip)=randnum

!transfor alpha to L
CALL alpha_to_L(alpha,vechL,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop2

!check overlap 
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr) 
CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechL(:)

CALL overlap_check_MPIfun(bk,Nbasis,ERR)
IF(ERR==0)EXIT generate_loop2

ENDDO generate_loop2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(generate_index<generate_limit)THEN

CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_temp,C_temp,ERR,0)
IF(E_temp<E_best.AND.ERR==0)THEN
  E_best=E_temp
  alpha_best(ip)=alpha(ip)
  L_best(:)=vechL(:)
ELSE
  alpha(ip)=alpha_best(ip)
ENDIF

ELSE
    
  alpha(ip)=alpha_best(ip)
  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO IT1_loop
ENDDO ip_loop
ENDIF
ENDDO IT2_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Glob_Lk(:,bk)=L_best(:)
CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_best,C_best,ERR,1)
IF(E_best<Eval.AND.ERR==0)THEN
  Eval=E_best
  Ci(:)=C_best(:)
ELSE
  Glob_Lk(:,bk)=L_old(:)
  CALL E_rows_MPIfun(1,bk_cal,Nbasis,Eval,Ci,ERR,1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_svm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!optimization of mk in svmopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mk_svm_MPIfun(bk,Nbasis,Eval,Ci)
!================================
!Optimization of parameter mk
!strategy:
!calculate the eigen values for each mk(1->Glob_Nparticle)
!and find the one making the value lowest 
!input:
!    bk: basis index
!Nbasis：present basis number
!  Eval: eigen value
!    Ci: eigen vector
!output:
!  Eval: eigen value
!    Ci: eigen vector
!  best mk stored in Glob_mk(bk)
!===============================
!this subroutine should 
!CALL pvars_MPImat(Nbasis,1)
!CALL pvars_eigen(Nbasis,1)
!before using it
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,Nbasis
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)
  
INTEGER::i,j,k,L,myid,ERR,ip
INTEGER::mk_best,mk_old,bk_cal(2)
REAL(dp)::E_best,C_best(Nbasis)
REAL(dp)::E_temp,C_temp(Nbasis)
  
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
  
bk_cal(1)=bk
bk_cal(2)=bk
!CALL E_rows_MPIfun(1,bk_cal,Nbasis,Eval,Ci,ERR,0)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!optimization
!============================
  
E_best=Eval
mk_best=Glob_mk(1,bk)
mk_old=mk_best
  
ip_loop:DO ip=1,Glob_Np
  
Glob_mk(1,bk)=ip
CALL overlap_check_MPIfun(bk,Nbasis,ERR)
  
IF(ERR==0)THEN
  
  CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_temp,C_temp,ERR,0)
  IF(E_temp<E_best.AND.ERR==0)THEN
    E_best=E_temp
    mk_best=ip
  ELSE
    Glob_mk(1,bk)=mk_best
  ENDIF
  
ELSE
    
  Glob_mk(1,bk)=mk_best
  CYCLE ip_loop
  
ENDIF
  
ENDDO ip_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
Glob_mk(1,bk)=mk_best
CALL E_rows_MPIfun(1,bk_cal,Nbasis,E_temp,C_temp,ERR,1)
  
IF(E_temp<Eval.AND.ERR==0)THEN
  Eval=E_temp
  Ci(:)=C_temp(:)
ELSE
  Glob_mk(1,bk)=mk_old
  CALL E_rows_MPIfun(1,bk_cal,Nbasis,Eval,Ci,ERR,1)
ENDIF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE mk_svm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!random number generating function for alpah->Ak->Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_mode_randnum(ith,mode,num)
!===============================
!generate random number of given mode of alpha
!input:
!  ith: the ith parameter
! mode: random generating mode
!inout:
!  num: if previous nolinear parameter is not require, 
!       num is noly out; otherwise num is inout
!===============================
IMPLICIT NONE
INTEGER,INTENT(IN)::ith
INTEGER,INTENT(IN)::mode
REAL(dp),INTENT(INOUT)::num

INTEGER::i,j,k
REAL(dp)::x1,x2,y,z
REAL(dp)::a(5)

DO i=1,5
  a(i)=Glob_Lk_mode_para(i,ith)
ENDDO

SELECT CASE(mode)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
!======================
!num=1.d0/(a1+a2*x1)**2
!======================
x1=randnum_fun(Lkseed1)
num=1.d0/(a(1)+a(2)*x1)**2
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
!======================
!num=a1*(a2/a1)**x1  
!======================
x1=randnum_fun(Lkseed1)
num=a(1)*(a(2)/a(1))**x1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(3)
!=======================
!num=a1*tan(x1*PI/2)**a2   
!=======================
x1=randnum_fun(Lkseed1) 
num=a(1)*tan(x1*PI/2.d0)**a(2)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(4)   
!=======================
!y=tan((x1-0.5)*PI)
!num=a1*y*abs(y)**(a2-1)
!=======================
x1=randnum_fun(Lkseed1)
y=tan((x1-0.5d0)*PI)
num=a(1)*y*dabs(y)**(a(2)-1.d0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(5)   
!========================
!Local optimization choice : Can only be used for
!improving existing basism, not for enlarging basis
!y=(1-a1)
!z=(1+a2)
!num=pre_para*y*(z/y)**x1
!========================
x1=randnum_fun(Lkseed1)
x2=num
y=(1.d0-a(1))
z=(1.d0+a(2))
IF(num==ZERO)THEN
  x2=randnum_fun(Lkseed2)
  x2=-1.d0+2.d0*x2
ENDIF
num=x2*y*(z/y)**x1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(6)
!========================
!num=a1*a2*(a3/a2)**x1
!where a1 is the sign number
!a1=+1: for attractive interactions 
!a1=-1: for attractive interactions 
!========================   
x1=randnum_fun(Lkseed1)
num=a(1)*a(2)*(a(3)/a(2))**x1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(7)
!========================
!0<x1<1: y=(a5/a1)**2
!1<x1<2: y=(a5/a2)**2
!2<x1<3: y=(a5/a3)**2
!num=y*tan(x2*PI/2)**a4	
!========================
x1=randnum_fun(Lkseed1)
x1=x1*3
x2=randnum_fun(Lkseed2)

IF((0.d0<x1).AND.(x1<1.d0))THEN
  y=(a(5)/a(1))**2
ELSEIF((1.d0<x1).AND.(x1<2.d0))THEN
  y=(a(5)/a(2))**2
ELSEIF((2.d0<x1).AND.(x1<3.d0))THEN
  y=(a(5)/a(3))**2
ENDIF

num=y*tan(x2*PI/2.d0)**a(4)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

RETURN
END SUBROUTINE Lk_mode_randnum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!exchange function of alpha and L
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE alpha_to_L(alpha,vechL,ERR)
!============================
!first transform alpha ij (pair correlation) to A matrix (nonlinear parameters 
!in Jacobi or heavy particle coordinite), then check wheather it's
!positive definite, if it's not, return ERR=1 otherwise ERR=0,
!and decompose A into lower triangular matrix L (vechL)
!input:
!   alpha: pair correlation parameters 
!output:
!  vechL: lower triangular matrix L  corresbonding to A matrix
!  ERR:
!     0:A is positive definite
!     1:A is not positive definite
!============================
IMPLICIT NONE
REAL(dp),INTENT(IN)::alpha(Glob_NLk)
REAL(dp),INTENT(OUT)::vechL(Glob_NLk)
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L,index1,index2
REAL(dp)::temp
REAL(dp)::vechA(Glob_NLk)

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!transform alpha to A
index1=0
DO k=1,Glob_Np
  DO L=k,Glob_Np
    index1=index1+1
    temp=ZERO
    index2=0
    DO i=1,Glob_Nparticle-1
      DO j=i+1,Glob_Nparticle
        index2=index2+1
        temp=temp+alpha(index2)*Glob_wij(k,i,j)*Glob_wij(L,i,j)
      ENDDO
    ENDDO
    vechA(index1)=temp
  ENDDO
ENDDO

!check wheather A matrix is positive definite
CALL PD_check_fun(Glob_Np,vechA,ERR)
IF(ERR/=0)THEN
  ERR=1
  RETURN
ENDIF

!decompose A matrix into L(vechL)
CALL Ckolesky_decompose(Glob_Np,vechA,vechL)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE alpha_to_L

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE L_to_alpha(vechL,alpha)
!============================
!first transform L to A matrix : A=L*L^T, then transform A matrix
!to alpha
!input:
!  vechL: lower triangular matrix L  corresbonding to A matrix
!output:
!   alpha: pair correlation parameters 
!============================
IMPLICIT NONE
REAL(dp),INTENT(IN)::vechL(Glob_NLk)
REAL(dp),INTENT(OUT)::alpha(Glob_NLk)

INTEGER::i,j,k,L,index
REAL(dp)::temp
REAL(dp)::Amatrix(Glob_Np,Glob_Np),Lmatrix(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

index=0
DO i=1,Glob_Np
  DO j=i,Glob_Np
    index=index+1
    Lmatrix(i,j)=ZERO
    Lmatrix(j,i)=vechL(index)     !lower
  ENDDO
ENDDO

!form A matrix
DO i=1,Glob_Np
  DO j=i,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+Lmatrix(i,k)*Lmatrix(j,k)
    ENDDO
    Amatrix(i,j)=temp
    Amatrix(j,i)=temp
  ENDDO
ENDDO

!transform A to alpha
index=0
DO i=1,Glob_Nparticle-1
  DO j=i+1,Glob_Nparticle
    index=index+1
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_U(k,i)*Amatrix(k,L)*Glob_U(L,j)
      ENDDO
    ENDDO
    alpha(index)=-temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE L_to_alpha

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE svmopt
