PROGRAM MAIN
!====================================================
!structure calculation program based on ECG basis
!====================================================
!ECG basis form:
!L=0,M=0: psi=exp(-rAr)
!L=1,M=0: psi=Z_mk*exp(-rAr)
!====================================================
!written by @ Dengxin Zhao
!====================================================
USE MPI
USE globvars
USE IO
USE MPImat
USE gvmopt
USE svmopt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!variables used in main.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
INTEGER::i,j,k,L,myid
INTEGER::it,bi,bj,dN
INTEGER::ktime,Nb_last

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!prgram initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
!MPI initialize
!=========================

CALL MPI_INIT(Glob_MPIerr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,Glob_numprocs,Glob_MPIerr)
ALLOCATE(Glob_status(MPI_STATUS_SIZE))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
!initialization
!=========================

!read the input parameters from input.txt
CALL readfile()

!write neccesary input parameters into log.txt
CALL writefile()

!initialize basis optimization parameters
CALL init_paras()

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL CPU_TIME(Glob_start_time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!optimization strategy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(1)==1)THEN
SELECT CASE(Glob_opt_strategy)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!strategy 1:
!   If don't increase basis, then optimize with GVM only
!   If increase basis, optimize with GVM after each increasing
!==============================
CASE(1)

IF((Glob_increase_num==0).OR.(Glob_Nbasis_start==Glob_Nbasis_final))THEN

DO i=1,Glob_opt_rounds_limit
  IF(i>1)Glob_opt_reach=1
  CALL gvm_opt(Glob_Nbasis_start,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 
ENDDO

ELSE

!increase time
dN=Glob_Nbasis_final-Glob_Nbasis_start
Nb_last=mod(dN,Glob_increase_num)
ktime=(dN-Nb_last)/Glob_increase_num
    
increase_loop1:DO it=1,ktime

!bi->bj basis
bi=Glob_Nbasis_start+(it-1)*Glob_increase_num
bj=bi+Glob_increase_num
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!GVM optimization 
Glob_opt_reach=1
CALL gvm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 

ENDDO increase_loop1

IF(Nb_last/=0)THEN
  
!bi->bj basis
bi=Glob_Nbasis_final-Nb_last
bj=Glob_Nbasis_final
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)
  
!GVM optimization
Glob_opt_reach=1
CALL gvm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 

ENDIF

DO i=1,Glob_opt_rounds_limit
  Glob_opt_reach=1
  CALL gvm_opt(Glob_Nbasis_final,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 
ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!!strategy 2:
!   If don't increase basis, then optimize with SVM only
!   If increase basis, optimize with SVM after each increasing
!===============================
CASE(2)
 
IF((Glob_increase_num==0).OR.(Glob_Nbasis_start==Glob_Nbasis_final))THEN
 
DO i=1,GLob_opt_rounds_limit
  IF(i>1)Glob_opt_reach=1
  CALL svm_opt(Glob_Nbasis_start,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds)
ENDDO

ELSE

!increase time
dN=Glob_Nbasis_final-Glob_Nbasis_start
Nb_last=mod(dN,Glob_increase_num)
ktime=(dN-Nb_last)/Glob_increase_num       

increase_loop2:DO it=1,ktime
 
!bi->bj basis
bi=Glob_Nbasis_start+(it-1)*Glob_increase_num
bj=bi+Glob_increase_num
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!SVM optimization
Glob_opt_reach=1
CALL svm_opt(bj,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds) 

ENDDO increase_loop2

IF(Nb_last/=0)THEN
    
!bi->bj basis 
bi=Glob_Nbasis_final-Nb_last
bj=Glob_Nbasis_final
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!SVM optimization
Glob_opt_reach=1
CALL svm_opt(bj,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds)
  
ENDIF

DO i=1,Glob_opt_rounds_limit
  Glob_opt_reach=1
  CALL svm_opt(Glob_Nbasis_final,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds) 
ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!straregy 3:
! optimize the basis with both SVM and GVM
! perform both optimization one by one: GVM-SVM-GVM-SVM
!============================
CASE(3)
        
IF((Glob_increase_num==0).OR.(Glob_Nbasis_start==Glob_Nbasis_final))THEN
 
DO i=1,Glob_opt_rounds_limit
  IF(i>1)Glob_opt_reach=1
  CALL gvm_opt(Glob_Nbasis_start,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds)
  Glob_opt_reach=1
  CALL svm_opt(Glob_Nbasis_start,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds)
ENDDO

ELSE

!increase time
dN=Glob_Nbasis_final-Glob_Nbasis_start
Nb_last=mod(dN,Glob_increase_num)
ktime=(dN-Nb_last)/Glob_increase_num         

increase_loop3:DO it=1,ktime
 
!bi->bj basis
bi=Glob_Nbasis_start+(it-1)*Glob_increase_num
bj=bi+Glob_increase_num
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!GVM optimization
Glob_opt_reach=1
CALL gvm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds)

!SVM optimization
Glob_opt_reach=1
CALL svm_opt(bj,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds) 

ENDDO increase_loop3


IF(Nb_last/=0)THEN
    
!bi->bj basis 
bi=Glob_Nbasis_final-Nb_last
bj=Glob_Nbasis_final
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!GVM optimization
Glob_opt_reach=1
CALL gvm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 

!SVM optimization
Glob_opt_reach=1
CALL svm_opt(bj,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds)
  
ENDIF

DO i=1,Glob_opt_rounds_limit
  Glob_opt_reach=1   
  CALL gvm_opt(Glob_Nbasis_final,Glob_gvm_batch,Glob_gvm_ITmax,Glob_gvm_rounds) 
  Glob_opt_reach=1
  CALL svm_opt(Glob_Nbasis_final,Glob_svm_IT1,Glob_svm_IT2,Glob_svm_rounds)
ENDDO

ENDIF   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!straregy 4:
! optimize each batch of basis with SVM-GVM
!============================
CASE(4)

IF((Glob_increase_num==0).OR.(Glob_Nbasis_start==Glob_Nbasis_final))THEN
 
DO i=1,Glob_opt_rounds_limit
  IF(i>1)Glob_opt_reach=1
  CALL gvm_svm_opt(Glob_Nbasis_start,Glob_gvm_batch,Glob_gvm_ITmax,&
  &Glob_svm_IT1,Glob_svm_IT2,max(Glob_gvm_rounds,Glob_svm_rounds))  
ENDDO 

ELSE

!increase time
dN=Glob_Nbasis_final-Glob_Nbasis_start
Nb_last=mod(dN,Glob_increase_num)
ktime=(dN-Nb_last)/Glob_increase_num    

increase_loop4:DO it=1,ktime
 
!bi->bj basis
bi=Glob_Nbasis_start+(it-1)*Glob_increase_num
bj=bi+Glob_increase_num
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!GVM+SVM optimization
Glob_opt_reach=1
CALL gvm_svm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,&
&Glob_svm_IT1,Glob_svm_IT2,max(Glob_gvm_rounds,Glob_svm_rounds))   

ENDDO increase_loop4

IF(Nb_last/=0)THEN
    
!bi->bj basis 
bi=Glob_Nbasis_final-Nb_last
bj=Glob_Nbasis_final
CALL basis_increase(bi,bj,Glob_gvm_ITmax,Glob_gvm_ITmax,Glob_svm_IT1,Glob_svm_IT2)

!GVM+SVM optimization
Glob_opt_reach=1
CALL gvm_svm_opt(bj,Glob_gvm_batch,Glob_gvm_ITmax,&
&Glob_svm_IT1,Glob_svm_IT2,max(Glob_gvm_rounds,Glob_svm_rounds)) 
  
ENDIF
  
DO i=1,Glob_opt_rounds_limit
  Glob_opt_reach=1
  CALL gvm_svm_opt(Glob_Nbasis_final,Glob_gvm_batch,Glob_gvm_ITmax,&
    &Glob_svm_IT1,Glob_svm_IT2,max(Glob_gvm_rounds,Glob_svm_rounds))  
ENDDO

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!structure information output
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(Glob_task_onoff(2)==1)THEN
!========================
!Correlation function gij 
!========================

IF(Glob_basis_form==0)THEN
  CALL gij_0_MPIfun(Glob_Nbasis_final)
ELSEIF(Glob_basis_form==1)THEN
  CALL gij_1_MPIfun(Glob_Nbasis_final)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!interparticle distance <rij>
!========================

IF(Glob_basis_form==0)THEN
  CALL rij_0_MPIfun(Glob_Nbasis_final,Glob_rij_power)
ELSEIF(Glob_basis_form==1)THEN
  CALL rij_1_MPIfun(Glob_Nbasis_final,Glob_rij_power)
ENDIF


ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_FINALIZE(Glob_MPIerr)
END PROGRAM MAIN