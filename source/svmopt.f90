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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to increase the basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE basis_increase(Nb_pre,Nb_after,gvm_IT1,gvm_IT2,svm_IT1,svm_IT2)
!===================================================
!this subroutine increase the basis size by choosing
!the new stochastic parameters lowering the energy
!===================================================
!strategy: 
!generate sets of parameters randomly one by one chose the one that
!lower the energy, then optimize the new basis gvm_IT1 times with gvm
!meanwhilde optimize the new basis with stochastic approach  
!generate a new basis accordingly, finally optimize all new basis
!from Nb_pre->Nb_after ITmax times with gvm 
!===================================================
!input:
!  Nb_pre: basis number before increase
!Nb_after: basis number after increase
! gvm_IT1: optimization times for new basis with gvm 
! gvm_IT2: optimization times of all dN new basis with gvm 
!          after increase the basis from Nb_pre to Nb_now
! svm_IT1: optimization times of each parameter of new basis with svm
! svm_IT2: optimization times of each basis with svm
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nb_pre,Nb_after
INTEGER,INTENT(IN)::gvm_IT1,gvm_IT2,svm_IT1,svm_IT2
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER,PARAMETER::Nb_init=60   !initial basis size
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER::i,j,k,L,myid,ERR,count,write_index
INTEGER::nb,Nb_start,mk,bk_cal(2)
REAL(dp)::Eval,Ci(Nb_after)
REAL(dp)::randnum,alpha(Glob_NLk),vechL(Glob_NLk)
CHARACTER(14)::time_st

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Nb_pre>=Nb_init)THEN

Nb_start=Nb_pre
    
ELSE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!IF(Nb_pre<Nb_init),increase the basis size 
!to the amount able to be diagonalized: Nb_init
!===============================
    
Nb_start=Nb_init
ERR=1
DO WHILE(ERR/=0)
ERR=0

init_loop:DO nb=Nb_pre+1,Nb_init
CALL pvars_MPImat(nb,1) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!generate Lk
DO i=1,Glob_NLk
  randnum=0.0_dp
  CALL Lk_mode_randnum_fun(i,randnum)
  alpha(i)=randnum
ENDDO
CALL alpha_to_L_fun(alpha,vechL,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop

IF(Glob_basis_form==1)THEN
  CALL MPI_BCAST(mk,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
  Glob_mk(nb)=mk
ENDIF

CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,nb)=vechL(:)

!check overlap
CALL Skl_check_MPIfun(nb,nb,ERR)

ENDDO generate_loop
CALL pvars_MPImat(nb,0) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO init_loop

CALL pvars_MPImat(Nb_start,1) 
bk_cal(1)=1
bk_cal(2)=Nb_start
CALL Eval_MPIfun(Nb_start,Nb_start,bk_cal,Eval,Ci(1:Nb_start),ERR,1)
CALL pvars_MPImat(Nb_start,0)

ENDDO 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF

CALL pvars_MPImat(Nb_start,1) 
bk_cal(1)=1
bk_cal(2)=Nb_start
CALL Eval_MPIfun(Nb_start,Nb_start,bk_cal,Eval,Ci(1:Nb_start),ERR,1)
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
!==========================
!if Nb_pre<Nb_init then optimize the initial basis 
!==========================
count=0
write_index=0
IF(Nb_pre<Nb_init)THEN
    
CALL pvars_MPImat(Nb_start,1) 
CALL pvars_eigen(Nb_start,1)
CALL pvars_gvmopt(1,1)

basis_loop1:DO nb=Nb_pre+1,Nb_start
count=count+1
bk_cal(1)=nb
bk_cal(2)=nb

CALL one_diag_init_MPIfun(Nb_start,nb)
!CALL Eval_MPIfun(Nb_start,1,bk_cal,Eval,Ci(1:Nb_start),ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(Nb_start,nb,Eval,Ci(1:Nb_start))
ENDIF

!CALL Lk_gvm_MPIfun(Nb_start,1,bk_cal,gvm_IT1,Eval,Ci(1:Nb_start))  
CALL Lk_svm_MPIfun(Nb_start,nb,svm_IT1,svm_IT2,Eval,Ci(1:Nb_start))

Glob_Nbasis_reach=Nb_start
Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)
WRITE(2,'(A14,I10,I14,f28.16)')time_st,count,nb,Eval
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1)
ELSE
  CALL write_parafile(2)
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO basis_loop1

CALL pvars_gvmopt(1,0)
CALL pvars_gvmopt(Nb_start,1)
bk_cal(1)=1
bk_cal(2)=Nb_start
CALL Eval_MPIfun(Nb_start,Nb_start,bk_cal,Eval,Ci(1:Nb_start),ERR,1)
CALL Lk_gvm_MPIfun(Nb_start,Nb_start,bk_cal,gvm_IT2,Eval,Ci(1:Nb_start))
CALL pvars_gvmopt(Nb_start,0)

CALL pvars_MPImat(Nb_start,0) 
CALL pvars_eigen(Nb_start,0)

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!increse basis and optimize it
!============================

IF(Nb_after>Nb_init)THEN
CALL pvars_gvmopt(1,1)

basis_loop2:DO nb=Nb_start+1,Nb_after 
count=count+1
bk_cal(1)=nb
bk_cal(2)=nb

CALL pvars_MPImat(nb,1)  
CALL pvars_eigen(nb,1)
CALL one_diag_init_MPIfun(nb,nb)
CALL basis_generate_MPIfun(nb,nb,Eval)
!CALL Eval_MPIfun(nb,1,bk_cal,Eval,Ci(1:nb),ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(nb,nb,Eval,Ci(1:nb))
ENDIF

!CALL Lk_gvm_MPIfun(nb,1,bk_cal,gvm_IT1,Eval,Ci(1:nb))
CALL Lk_svm_MPIfun(nb,nb,svm_IT1,svm_IT2,Eval,Ci(1:nb))

Glob_Nbasis_reach=nb
Glob_opt_reach=nb
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

CALL time_cal_fun(time_st)
WRITE(2,'(A14,I10,I14,f28.16)')time_st,count,nb,Eval
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) 
ELSE
  CALL write_parafile(2) 
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CALL pvars_MPImat(nb,0) 
CALL pvars_eigen(nb,0)
ENDDO basis_loop2

CALL pvars_gvmopt(1,0)
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
  CALL pvars_gvmopt(Nb_after-Nb_pre,1)
  CALL Lk_gvm_MPIfun(Nb_after,Nb_after-Nb_pre,bk_cal,gvm_IT2,Eval,Ci(1:Nb_after))
  CALL pvars_MPImat(Nb_after,0) 
  CALL pvars_gvmopt(Nb_after-Nb_pre,0)
  
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE basis_increase

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to increase basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE basis_generate_MPIfun(Nbasis,bk,Eval)
!===================================================
!generate one new basis with stategy:
!generate sets of parameters randomly one
!by one chose the one that lower the energy
!===================================================
!input:
!Nbasis: basis number after generating
!    bk: basis index to generate
!inout:
!  Eval: eigen value
!output:
!  new parameters stored in Glob_mk and Glob_Lk
!===================================================
!this subroutine should CALL pvars_eigen(Nbasis,1)
!and CALL one_diag_init_MPIfun(Nbasis,bk)   
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,bk
REAL(dp),INTENT(INOUT)::Eval

INTEGER::i,j,k,L,myid,ERR
INTEGER::mk,bk_cal(2)
REAL(dp)::randnum,Etemp,Ctemp(Nbasis)
REAL(dp)::alpha(Glob_NLk),vechL(Glob_NLk)

bk_cal(1)=bk
bk_cal(2)=bk
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ERR=1
loop:DO WHILE(ERR/=0)

ERR=1
generate_loop:DO WHILE(ERR/=0)
ERR=0

IF(myid==Glob_root)THEN

!generate mk
IF(Glob_basis_form==1)THEN
  randnum=randnum_fun(mkseed)
  mk=1+floor(randnum*Glob_Np)
ENDIF

!generate Lk
DO i=1,Glob_NLk
  randnum=0.0_dp
  CALL Lk_mode_randnum_fun(i,randnum)
  alpha(i)=randnum
ENDDO
CALL alpha_to_L_fun(alpha,vechL,ERR) 

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop

IF(Glob_basis_form==1)THEN
  CALL MPI_BCAST(mk,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
  Glob_mk(bk)=mk
ENDIF

CALL MPI_BCAST(vechL,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechL(:)

!overlap check
CALL Skl_check_MPIfun(bk,Nbasis,ERR)
IF(ERR/=0)CYCLE generate_loop

!check whether energy is lower
CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,0)
IF(Etemp>Eval)THEN
  ERR=1
  CYCLE generate_loop
ENDIF

ENDDO generate_loop

CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,1)
IF(Etemp>Eval)ERR=1

ENDDO loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Eval=Etemp

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE basis_generate_MPIfun 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gvm_svm_opt(Nbasis,Nbs,gvm_ITmax,svm_IT1,svm_IT2,Nrounds)
!===================================================
!optimize each batch of basis with SVM and GVM at same time
!===================================================
!input:
!   Nbasis: present basis number
!      Nbs: batch basis number
!gvm_ITmax: optimization times for new basis with gvm
!  svm_IT1: optimization times of each parameter of one basis with svm 
!  svm_IT2: optimization times of each basis with svm
!  Nrounds: total rounds of optimization
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs
INTEGER,INTENT(IN)::gvm_ITmax,svm_IT1,svm_IT2,Nrounds

INTEGER::i,j,k,myid,ERR
INTEGER::Nbatches,bk_cal(2),nr,nb,bk,count,nb_start,write_index
INTEGER,ALLOCATABLE,DIMENSION(:,:)::basis_index
REAL(dp)::Eval,Ci(Nbasis)
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================
!initialization
!=================================

CALL pvars_MPImat(Nbasis,1)
CALL pvars_eigen(Nbasis,1)
CALL pvars_gvmopt(Nbs,1)
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

bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before all_gvm_svm'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'gvm+svm optimization starts:'
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!optimize each basis batch 
!===============================

batch_loop:DO nb=nb_start,Nbatches 
count=count+1
bk_cal(:)=basis_index(:,nb)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!gvm optimization
!===============================

IF(Glob_basis_form==1)THEN
  CALL mk_gvm_MPIfun(Nbasis,Nbs,bk_cal,Eval,Ci)
ENDIF

CALL Lk_gvm_MPIfun(Nbasis,Nbs,bk_cal,gvm_ITmax,Eval,Ci)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!svm optimization
!===============================

DO bk=bk_cal(1),bk_cal(2) 

CALL one_diag_init_MPIfun(Nbasis,bk)

IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(Nbasis,bk,Eval,Ci)
ENDIF

CALL Lk_svm_MPIfun(Nbasis,bk,svm_IT1,svm_IT2,Eval,Ci)

ENDDO 

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

CALL time_cal_fun(time_st)
WRITE(2,'(A14,2I10,I14,f28.16)')time_st,count,nr,nb,Eval
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1)
ELSE
  CALL write_parafile(2)
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO batch_loop
ENDDO optimize_loop

IF(myid==Glob_root)THEN 
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'gvm+svm optimization finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,0)
CALL pvars_eigen(Nbasis,0)
CALL pvars_gvmopt(Nbs,0)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE gvm_svm_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!stockastic optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svm_opt(Nbasis,IT1,IT2,Nrounds)
!===================================================
!this subroutine optimizes the nonlinear parameters 
!with stochastic approach
!===================================================
!input:
! Nbasis: present basis number
!    IT1: optimization times for one parameter of one basis
!    IT2: optimization times for one basis
!Nrounds: total rounds of svm optimization
!===================================================
INTEGER,INTENT(IN)::Nbasis,IT1,IT2,Nrounds

INTEGER::i,j,k,myid,ERR,index
INTEGER::nr,nb,bk,nb_start,bk_cal(2)
INTEGER::count,write_index,generate_index
REAL(dp)::Eval,Ci(Nbasis)
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,1)
CALL pvars_eigen(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before all_svm'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'svm optimization starts:'
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

CALL one_diag_init_MPIfun(Nbasis,nb)
  
!==========================
!optimization of mk
!==========================
IF(Glob_basis_form==1)THEN
  CALL mk_svm_MPIfun(Nbasis,nb,Eval,Ci)
ENDIF

!===========================
!optimization of Lk with BFGs algrithm
!===========================
CALL Lk_svm_MPIfun(Nbasis,nb,IT1,IT2,Eval,Ci)

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 
  
CALL time_cal_fun(time_st)
WRITE(2,'(A14,2I10,I14,f28.16)')time_st,count,nr,nb,Eval
write_index=write_index+1
IF(write_index<=10)THEN
  CALL write_parafile(1) 
ELSE
  CALL write_parafile(2)
  IF(write_index==20)THEN
      write_index=0
  ENDIF
ENDIF

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO basis_loop
ENDDO optimize_loop

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'svm optimization finished'
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!optimization of mk in svmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mk_svm_MPIfun(Nbasis,bk,Eval,Ci)
!===================================================
!optimization of parameter mk with strategy:
!calculate the eigen values for each mk(1->Glob_Np)
!and accept the one making eigen value lowest as its new value 
!===================================================
!input:
!Nbasis: resent basis number
!    bk: basis index
!inout:
!  Eval: eigen value
!    Ci: eigen vector
!===================================================
!this subroutine should CALL pvars_MPImat(Nbasis,1)
!and CALL pvars_eigen(Nbasis,1)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,Nbasis
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)
  
INTEGER::i,j,k,L,myid,ERR,ip
INTEGER::mk_best,mk_old,bk_cal(2)
REAL(dp)::Ebest,Cbest(Nbasis),Etemp,Ctemp(Nbasis)
  
bk_cal(1)=bk
bk_cal(2)=bk
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!CALL Eval_MPIfun(Nbasis,1,bk_cal,Eval,Ci,ERR,0)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!optimization
!============================
  
Ebest=Eval
mk_best=Glob_mk(bk)
mk_old=mk_best
  
ip_loop:DO ip=1,Glob_Np
IF(ip/=mk_old)THEN

Glob_mk(bk)=ip
CALL Skl_check_MPIfun(Nbasis,bk,ERR)
  
IF(ERR==0)THEN
  
  CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,0)
  IF(Etemp<Ebest.AND.ERR==0)THEN
    Ebest=Etemp
    mk_best=ip
  ELSE
    Glob_mk(bk)=mk_best
  ENDIF
  
ELSE
    
  Glob_mk(bk)=mk_best
  CYCLE ip_loop
  
ENDIF
  
ENDIF
ENDDO ip_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
Glob_mk(bk)=mk_best
CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,1)
  
IF(Etemp<Eval.AND.ERR==0)THEN
  Eval=Etemp
  Ci(:)=Ctemp(:)
ELSE
  Glob_mk(bk)=mk_old
  CALL Eval_MPIfun(Nbasis,1,bk_cal,Eval,Ci,ERR,1)
  IF(ERR/=0)PAUSE'error in mk_svm_MPIfun'
ENDIF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE mk_svm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_svm_MPIfun(Nbasis,bk,IT1,IT2,Eval,Ci)  
!===================================================
!optimize parameter Lk of one basis with SVM approach 
!===================================================
!input:
!Nbasis: present basis number
!    bk: the basis to be optimized
!   IT1: optimization times for one parameter in bk basis
!   IT2: optimization times for bk basis
!inout: 
!  Eval: eigen value
!    Ci: eigen vector
!===================================================
!this subroutine should CALL pvars_MPImat(Nbasis,1)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,bk,IT1,IT2
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)

INTEGER::i,j,k,myid,ERR
INTEGER::t1,t2,ip,generate_index,bk_cal(2)
REAL(dp)::randnum,Ebest,Cbest(Nbasis),Etemp,Ctemp(Nbasis)
REAL(dp)::alpha_best(Glob_NLk),Lk_best(Glob_NLk)
REAL(dp)::Lk_old(Glob_NLk),alpha(Glob_NLk),vechLk(Glob_NLk)

IF(IT2==0)RETURN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bk_cal(1)=bk
bk_cal(2)=bk
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==================================
!stochastic optimization 
!==================================

Ebest=Eval
Lk_best(:)=Glob_Lk(:,bk)
Lk_old(:)=Lk_best(:)
CALL L_to_alpha_fun(Lk_best,alpha_best)

!optimize nb_th basis IT2 times
IT2_loop:DO t2=1,IT2

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!IT1==0 : do not optimize each parameter seperately 
!===========================
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
    
!generate Lk
DO i=1,Glob_NLk
  randnum=alpha_best(i)
  CALL Lk_mode_randnum_fun(i,randnum)
  alpha(i)=randnum
ENDDO
CALL alpha_to_L_fun(alpha,vechLk,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop1

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(vechLk,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechLk(:)

!check overlap 
CALL Skl_check_MPIfun(Nbasis,bk,ERR)
IF(ERR==0)EXIT generate_loop1

ENDDO generate_loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(generate_index<generate_limit)THEN

CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,0)
IF(Etemp<Ebest.AND.ERR==0)THEN
  Ebest=Etemp
  alpha_best(:)=alpha(:)
  Lk_best(:)=vechLk(:)
ENDIF
 
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!IT1/=0 : optimize each parameter seperately IT1 timies 
!=======================
ELSE

alpha(:)=alpha_best(:)

ip_loop:DO ip=1,Glob_NLk  
IT1_loop:DO t1=1,IT1
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
! generate nonlinear parameter
!==============================

generate_index=0
generate_loop2:DO WHILE(generate_index<generate_limit)
generate_index=generate_index+1
ERR=0

IF(myid==Glob_root)THEN
    
!generate Lk
randnum=alpha_best(ip)
CALL Lk_mode_randnum_fun(ip,randnum)
alpha(ip)=randnum
CALL alpha_to_L_fun(alpha,vechLk,ERR)

ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(ERR,1,MPI_INTEGER,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
IF(ERR/=0)CYCLE generate_loop2

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr) 
CALL MPI_BCAST(vechLk,Glob_NLk,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
Glob_Lk(:,bk)=vechLk(:)

!check overlap 
CALL Skl_check_MPIfun(Nbasis,bk,ERR)
IF(ERR==0)EXIT generate_loop2

ENDDO generate_loop2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(generate_index<generate_limit)THEN

CALL Eval_MPIfun(Nbasis,1,bk_cal,Etemp,Ctemp,ERR,0)

IF(Etemp<Ebest.AND.ERR==0)THEN
  Ebest=Etemp
  alpha_best(ip)=alpha(ip)
  Lk_best(:)=vechLk(:)
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

Glob_Lk(:,bk)=Lk_best(:)
CALL Eval_MPIfun(Nbasis,1,bk_cal,Ebest,Cbest,ERR,1)

IF(Ebest<Eval.AND.ERR==0)THEN
  Eval=Ebest
  Ci(:)=Cbest(:)
ELSE
  Glob_Lk(:,bk)=Lk_old(:)
  CALL Eval_MPIfun(Nbasis,1,bk_cal,Eval,Ci,ERR,1)
  IF(ERR/=0)PAUSE'error in Lk_svm_MPIfun'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_svm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!random number generating function for alpah->Ak->Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_mode_randnum_fun(ith,num)
!===================================================
!generate random number of given mode of alpha
!===================================================
!input:
!  ith: the ith parameter
!inout:
!  num: random number of alpha(ith)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ith
REAL(dp),INTENT(INOUT)::num

INTEGER::i,j,k
REAL(dp)::x1,x2,y,z
REAL(dp)::a(5)

DO i=1,5
  a(i)=Glob_Lk_mode_para(i,ith)
ENDDO

SELECT CASE(Glob_Lk_mode(ith))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
IF(num==0.0_dp)THEN
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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

RETURN
END SUBROUTINE Lk_mode_randnum_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!transform pair correlation alpha to L_matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE alpha_to_L_fun(alpha,vechL,ERR)
!===================================================
!transform pair correlation alpha to L matrix:
!first transform alpha_ij to A_matrix, then check 
!wheather it's positive definite, if it's not, 
!return ERR=1, otherwise ERR=0, finally decompose 
!A_martix into lower triangular matrix L
!===================================================
!input:
! alpha: pair correlation parameters 
!output:
! vechL: lower triangular matrix L
!   ERR:
!       0: A is positive definite
!       1: A is not positive definite
!===================================================
IMPLICIT NONE
REAL(dp),INTENT(IN)::alpha(Glob_NLk)
REAL(dp),INTENT(OUT)::vechL(Glob_NLk)
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L,index1,index2
REAL(dp)::temp
REAL(dp)::vechA(Glob_NLk)

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!transform alpha to A
index1=0
DO k=1,Glob_Np
  DO L=k,Glob_Np
    index1=index1+1
    temp=0.0_dp
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE alpha_to_L_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!transform L_matrix to pair correlation alpha
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE L_to_alpha_fun(vechL,alpha)
!===================================================
!transform L matrix into pair correlation alpha:
!first transform L to A matrix: A=L*L^T, then 
!transform A matrix to alpha
!===================================================
!input:
!  vechL: lower triangular matrix L  corresbonding to A matrix
!output:
!   alpha: pair correlation parameters 
!===================================================
IMPLICIT NONE
REAL(dp),INTENT(IN)::vechL(Glob_NLk)
REAL(dp),INTENT(OUT)::alpha(Glob_NLk)

INTEGER::i,j,k,L,index
REAL(dp)::temp
REAL(dp)::Amatrix(Glob_Np,Glob_Np),Lmatrix(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

index=0
DO i=1,Glob_Np
  DO j=i,Glob_Np
    index=index+1
    Lmatrix(i,j)=0.0_dp
    Lmatrix(j,i)=vechL(index)     !lower
  ENDDO
ENDDO

!form A matrix
DO i=1,Glob_Np
  DO j=i,Glob_Np
    temp=0.0_dp
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
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_U(k,i)*Amatrix(k,L)*Glob_U(L,j)
      ENDDO
    ENDDO
    alpha(index)=-temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE L_to_alpha_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE svmopt
