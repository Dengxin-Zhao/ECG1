MODULE gvmopt
!===========================================================
!this module contians subroutines for BFGs optimization 
!with parallel approach
!===========================================================
USE MPI
USE auxfun
USE globvars
USE MPImat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!private variables and matrix used in gvmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================================
!these virables should be recalculated 
!when the batch basis number is changed
!========================================
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!matrix related to process distribution
!==========================
INTEGER,PRIVATE::p_Nmax
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Nproc,p_displs
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:,:,:)::p_paracal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!variables used in BFGs optimization
!==========================
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_x1_each,p_x2_each
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dx_each,p_p_each
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_H_ddf_each
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_uterm_each
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_hessin_each !hessian matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dx            !x_(i+1)-x_i
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_df1           !df_(i+1)
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_df2           !df_i
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_ddf           !df_(i+1)-df_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_H_ddf         !Hessin*(df_i+1-df_i)
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_uterm         !additional term of BFGs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate and deallocate variables used in gvmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE pvars_gvmopt(Nbs,case_num)
!================================
!allocate and deallocate private varibles used in gvmopt.f90
!input:
!      Nbs: number of basis in each batch
!   Nbasis: present basis number
! case_num:
!           1: allocate the variables
!           0: deallocate the variables
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,case_num

INTEGER::i,j,k,L,index1,index2
INTEGER::paranum,N_last,N_each

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

SELECT CASE(case_num)   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)!allocate

paranum=Glob_NLk*Nbs               !numbers of parameters in each batch
N_last=mod(paranum,Glob_numprocs)
N_each=(paranum-N_last)/Glob_numprocs
p_Nmax=N_each+1
IF(N_last==0)THEN
p_Nmax=p_Nmax-1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(p_x1_each(p_Nmax),p_x2_each(p_Nmax))
ALLOCATE(p_dx_each(p_Nmax),p_p_each(p_Nmax))
ALLOCATE(p_H_ddf_each(p_Nmax),p_uterm_each(p_Nmax))
ALLOCATE(p_hessin_each(paranum,p_Nmax))
!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(p_dx(paranum))
ALLOCATE(p_df1(paranum))
ALLOCATE(p_df2(paranum))
ALLOCATE(p_ddf(paranum))
ALLOCATE(p_H_ddf(paranum))
ALLOCATE(p_uterm(paranum))
!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(p_Nproc(Glob_numprocs))
ALLOCATE(p_displs(Glob_numprocs))
ALLOCATE(p_paracal(2,p_Nmax,Glob_numprocs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!numbers of parameters manipulated by each process
DO i=1,Glob_numprocs   
  IF(i<=N_last)THEN
    p_Nproc(i)=N_each+1   
  ELSE
    p_Nproc(i)=N_each 
  ENDIF
ENDDO

!displacement when gatherv from every process
p_displs(1)=0
IF(Glob_numprocs>1)THEN
DO i=2,Glob_numprocs
  p_displs(i)=p_displs(i-1)+p_Nproc(i-1)
ENDDO
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!location of parameters calculated by each process
p_paracal=0
index1=1 !process sign
index2=0 !parameter number sign
DO i=1,Nbs
  DO j=1,Glob_NLk
    index2=index2+1
    p_paracal(1,index2,index1)=i
    p_paracal(2,index2,index1)=j
    IF((index2==(N_each+1)).AND.(index1<=N_last))THEN
      index1=index1+1
      index2=0
    ELSEIF((index2==N_each).AND.(index1>N_last))THEN
      index1=index1+1
      index2=0
    ENDIF     
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)!deallocate
    

DEALLOCATE(p_x1_each,p_x2_each)
DEALLOCATE(p_dx_each,p_p_each)
DEALLOCATE(p_H_ddf_each,p_uterm_each)
DEALLOCATE(p_hessin_each)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DEALLOCATE(p_dx)
DEALLOCATE(p_df1)
DEALLOCATE(p_df2)
DEALLOCATE(p_ddf)
DEALLOCATE(p_H_ddf)
DEALLOCATE(p_uterm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DEALLOCATE(p_Nproc)
DEALLOCATE(p_displs)
DEALLOCATE(p_paracal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE pvars_gvmopt
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!main subroutine performing gradient optimization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gvm_opt(Nbasis,Nbs,ITmax,Nrounds)
!===================================================
!this subroutine optimize the nonlinear parameters
!of each batch one by one with BFGs algrithm
!===================================================
!input:
! Nbasis: present basis number
!    Nbs: batch number
!  ITmax: max iterations of gvm
!Nrounds: optimization rounds of gvm
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs,ITmax,Nrounds

INTEGER::i,j,k,L,myid,ERR
INTEGER::Nbatches,bk_cal(2),nr,nb,it,nb_start,count,write_index
INTEGER,ALLOCATABLE,DIMENSION(:,:)::basis_index
REAL(dp)::Eval,Ci(Nbasis)
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!initialization
!===============================

CALL pvars_MPImat(Nbasis,1)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)
IF(ERR/=0)PAUSE'eigen solution error before gvm_opt'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'==================================================' 
  WRITE(2,*)'gvm optimization with BFGs start:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'batch basis number:',Nbs
  WRITE(2,*)'=================================================='
  WRITE(2,*)'start energy:',Eval
  WRITE(2,*)'=================================================='
  WRITE(2,'(A14,2A10,A14,A11)')'time','count','round','batch_num','energy'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================
!BFGs optimization
!================================

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

batch_loop:DO nb=nb_start,Nbatches 
bk_cal(:)=basis_index(:,nb)
count=count+1

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
  
!==========================
!optimization of mk
!==========================
IF(Glob_basis_form==1)THEN
  CALL mk_gvm_MPIfun(Nbasis,Nbs,bk_cal,Eval,Ci)
ENDIF

!===========================
!optimization of Lk with BFGs algrithm
!===========================
CALL Lk_gvm_MPIfun(Nbasis,Nbs,bk_cal,ITmax,Eval,Ci)

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'==================================================' 
  WRITE(2,*)'gvm optimization with BFGs finished'
  WRITE(2,*)'=================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,0)
CALL pvars_gvmopt(Nbs,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE gvm_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine for optimization of mk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mk_gvm_MPIfun(Nbasis,Nbs,bk_cal,Eval,Ci)
!===================================================
!optimization of parameter mk with strategy:
!calculate the eigen values for each mk(1->Glob_Nparticle)
!and accept the one making the eigen value lowest as new value 
!===================================================
!input:
!Nbasis: present basis number
!   Nbs: number of batch
!bk_cal: basis index: bk_cal(1)->bk_cal(2)
!inout:
!  Eval: eigen value
!    Ci: eigen vector
!===================================================
!this subroutine should CALL pvars_MPImat(Nbasis,1)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs,bk_cal(2)
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)

INTEGER::i,j,k,L,myid,ERR
INTEGER::bi,bk,ip,basis_cal(2),mk_best,mk_old
REAL(dp)::Ebest,Cbest(Nbasis),Etemp,Ctemp(Nbasis)

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!CALL Eval_MPIfun(Nbs,bk_cal,Nbasis,Eval,Ci,ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!optimization
!============================

Ebest=Eval

bk_loop:DO bk=bk_cal(1),bk_cal(2)

basis_cal(1)=bk
basis_cal(2)=bk
mk_best=Glob_mk(bk)
mk_old=mk_best

ip_loop:DO ip=1,Glob_Np
IF(ip/=mk_old)THEN

Glob_mk(bk)=ip
CALL Skl_check_MPIfun(Nbasis,bk,ERR)

IF(ERR==0)THEN

  CALL Eval_MPIfun(Nbasis,1,basis_cal,Etemp,Ctemp,ERR,1)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Glob_mk(bk)=mk_best
CALL Eval_MPIfun(Nbasis,1,basis_cal,Etemp,Ctemp,ERR,1)

IF(Etemp<Eval.AND.ERR==0)THEN
  Eval=Etemp
  Ci(:)=Ctemp(:)
ELSE
  Glob_mk(bk)=mk_old
  CALL Eval_MPIfun(Nbasis,1,basis_cal,Eval,Ci,ERR,1)
  IF(ERR/=0)PAUSE'error in mk_gvm_MPIfun'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO bk_loop

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE mk_gvm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient optimization of Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_gvm_MPIfun(Nbasis,Nbs,bk_cal,ITmax,f1,C1)
!===================================================
!this subroutine optimize the nonlinear parameters Lk
!with BFGs algrithm
!===================================================
!input:
!Nbasis: present basis number
!   Nbs: batch number
!bk_cal: basis index: bk_cal(1)->bk_cal(2) 
! ITmax: max iterations of gvm
!inout:
!    f1: eigen value
!    C1: eigen vector
!===================================================
!this subroutine should CALL pvars_MPImat(Nbasis,1)
!and CALL pvars_gvmopt(Nbs,Glob_NLk,1)
!=================================================== 
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs,bk_cal(2),ITmax
REAL(dp),INTENT(INOUT)::f1,C1(Nbasis)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
REAL(dp),PARAMETER::TOL_dx=3.0_dp*EPS !convergence criterion on dx
REAL(dp),PARAMETER::TOL_df=3.0_dp*EPS !convergence criterion on df
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
REAL(dp)::stpmax                      !max step allowed
REAL(dp)::ddf_dx                      !(df_i+1-df_i)(x_(i+1)-x_i)
REAL(dp)::ddf_H_ddf                   !(df_i+1-df_i)Hessin(df_i+1-df_i)
REAL(dp)::sum_ddf,sum_dx              !(df_i+1-df_i)^2,(x_(i+1)-x_i)^2
REAL(dp)::temp,test,f2,C2(Nbasis) 

INTEGER::i,j,k,L,myid,ERR
INTEGER::it,bi,bj,Nproc,paranum
INTEGER,ALLOCATABLE,DIMENSION(:,:)::paracal              

IF(ITmax==0)RETURN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

!distribution of process
Nproc=p_Nproc(myid+1)      
ALLOCATE(paracal(2,p_Nmax))         
paracal(:,:)=p_paracal(:,:,myid+1)

!bi->bj basis
bi=bk_cal(1)
bj=bk_cal(2)
paranum=Nbs*Glob_NLk

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!initialization of BFGs
!==============================

!f1 and df1
!CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f1,C1,ERR,1)
CALL dE_dLk_MPIfun(Nbasis,Nbs,bk_cal,f1,C1,p_df1)

temp=0.0_dp
p_hessin_each(:,:)=0.0_dp
IF(Nproc/=0)THEN    
DO i=1,Nproc
    
  j=bi+paracal(1,i)-1
  k=paracal(2,i)
  p_x1_each(i)=Glob_Lk(k,j)          !x1
  temp=temp+p_x1_each(i)**2
  
  j=paracal(1,i)
  k=paracal(2,i)+(j-1)*Glob_NLk
  p_hessin_each(k,i)=1.0_dp          !hessin=I
  
ENDDO
ENDIF

test=temp
CALL MPI_ALLREDUCE(test,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

!max step allowed (avoid unreasonable steps)
stpmax=100.E0_dp*max(dsqrt(temp),dble(paranum))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================
!optimization loop with BFGs
!=================================
BFGs_loop:DO it=1,ITmax
    
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!================================
!new direction: p_p=-Hessin*p_df
!================================
IF(Nproc/=0)THEN
DO i=1,Nproc

  temp=0.0_dp
  DO k=1,paranum
    temp=temp-p_hessin_each(k,i)*p_df1(k)
  ENDDO
  p_p_each(i)=temp

ENDDO
ENDIF

!================================
!scale direction p if it is too big
!================================
temp=0.0_dp
IF(Nproc/=0)THEN 
DO i=1,Nproc
  temp=temp+p_p_each(i)**2
ENDDO
ENDIF

test=temp
CALL MPI_ALLREDUCE(test,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

IF(Nproc/=0)THEN
    
temp=dsqrt(temp)
IF(temp>stpmax)THEN
  DO i=1,Nproc
    p_p_each(i)=p_p_each(i)*stpmax/temp
  ENDDO
ENDIF

ENDIF

!==============================
!step search: x1->x2 (x2=x1+alpha*p)
!==============================
CALL Lk_stepsrch_MPIfun(Nbasis,Nbs,bk_cal,f1,C1,f2,C2,TOL_dx)

!dx
p_dx_each=p_x2_each-p_x1_each 

!=============================
!dx check: max(dx)<TOL_dx -> exit BFGs_loop 
!=============================
test=0.0_dp
IF(Nproc/=0)THEN
DO i=1,Nproc

  temp=dabs(p_dx_each(i))/max(dabs(p_x2_each(i)),1.0_dp)
  IF(temp>test)THEN
    test=temp
  ENDIF

ENDDO
ENDIF

temp=test
CALL MPI_ALLREDUCE(temp,test,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,Glob_MPIerr)
IF(test<TOL_dx)EXIT BFGs_loop 

!=============================
!calculate new gradient df2
!=============================

CALL dE_dLk_MPIfun(Nbasis,Nbs,bk_cal,f2,C2,p_df2)

!=============================
!df check: max(df2)<TOL_df -> exit BFGs_loop
!=============================
test=0.0_dp
IF(Nproc/=0)THEN 
DO i=1,Nproc

  j=paracal(1,i)
  k=paracal(2,i)+Glob_NLk*(j-1)
  temp=dabs(p_df2(k))*max(dabs(p_x2_each(i)),1.0_dp)/max(f2,1.0_dp)
  IF(temp>test)THEN 
    test=temp
  ENDIF

ENDDO
ENDIF

temp=test
CALL MPI_ALLREDUCE(temp,test,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,Glob_MPIerr)
IF(test<TOL_df)EXIT BFGs_loop 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================
!renew hessin matrix procedure
!================================

p_ddf(:)=p_df2(:)-p_df1(:)   

ddf_dx=0.0_dp
ddf_H_ddf=0.0_dp  
sum_ddf=0.0_dp
sum_dx=0.0_dp

IF(Nproc/=0)THEN
DO i=1,Nproc
    
  !Hessin*(df_i+1-df_i)
  temp=0.0_dp  
  DO k=1,paranum
    temp=temp+p_hessin_each(k,i)*p_ddf(k)
  ENDDO
  p_H_ddf_each(i)=temp

  j=paracal(1,i)
  k=paracal(2,i)+Glob_NLk*(j-1)
  
  !(df_i+1-df_i)(x_(i+1)-x_i)
  ddf_dx=ddf_dx+p_ddf(k)*p_dx_each(i)  
  
  !(df_i+1-df_i)Hessin(df_i+1-df_i)
  ddf_H_ddf=ddf_H_ddf+p_ddf(k)*p_H_ddf_each(i)

  sum_ddf=sum_ddf+p_ddf(k)**2

  sum_dx=sum_dx+p_dx_each(i)**2
  
ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

temp=ddf_dx
CALL MPI_ALLREDUCE(temp,ddf_dx,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

temp=ddf_H_ddf
CALL MPI_ALLREDUCE(temp,ddf_H_ddf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

temp=sum_ddf
CALL MPI_ALLREDUCE(temp,sum_ddf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

temp=sum_dx
CALL MPI_ALLREDUCE(temp,sum_dx,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

!==============================
!skip update hessin if ddf_dx is not sufficiently positive
!==============================
IF(ddf_dx>dsqrt(EPS*sum_ddf*sum_dx))THEN
 
IF(Nproc/=0)THEN   
DO i=1,Nproc
  p_uterm_each(i)=p_dx_each(i)/ddf_dx-p_H_ddf_each(i)/ddf_H_ddf
ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(p_dx_each,Nproc,MPI_DOUBLE_PRECISION,p_dx,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(p_uterm_each,Nproc,MPI_DOUBLE_PRECISION,p_uterm,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(p_H_ddf_each,Nproc,MPI_DOUBLE_PRECISION,p_H_ddf,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)


IF(Nproc/=0)THEN
DO i=1,Nproc

  j=paracal(1,i)
  k=paracal(2,i)+(j-1)*Glob_NLk
  DO j=1,paranum
    p_hessin_each(j,i)=p_hessin_each(j,i)+p_dx(k)*p_dx(j)/ddf_dx-p_H_ddf(k)*p_H_ddf(j)/ddf_H_ddf&
    &+ddf_H_ddf*p_uterm(k)*p_uterm(j)
  ENDDO

ENDDO
ENDIF

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!update
!============================
p_x1_each=p_x2_each
f1=f2
C1=C2
p_df1=p_df2    

ENDDO BFGs_loop

CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f1,C1,ERR,1)
IF(ERR/=0)PAUSE'error in Lk_gvm_MPIfun'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_gvm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!step search function for parameter Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_stepsrch_MPIfun(Nbasis,Nbs,bk_cal,f1,C1,f2,C2,TOL_dx)
!===================================================
!This subroutine search the step of parameters with 
!backtracking algrithm
!===================================================
!input:
!Nbasis: present basis number
!   Nbs: batch number 
!bk_cal: basis index: bk_cal(1)->bk_cal(2)
!    f1: old eigen value
!    C1: old eigen vector
!TOL_dx: tolerance of dx 
!output:
!     f2: new objective value
!     C2: new eigen vector corresponding to f2 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs,bk_cal(2)
REAL(dp),INTENT(IN)::f1,C1(Nbasis),TOL_dx
REAL(dp),INTENT(OUT)::f2,C2(Nbasis)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
REAL(dp),PARAMETER::d1=1.d-4  !sufficient decrease coef
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER::i,j,k,L,myid,index,ERR,bk,Nproc
INTEGER,ALLOCATABLE,DIMENSION(:,:)::paracal

REAL(dp)::df_p,alpha_min,alpha1,alpha2,ftemp
REAL(dp)::g1,g2,a,b,test,temp,Lk_old(Glob_NLk*Nbs),x2(Glob_NLk*Nbs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

!distribution of process
Nproc=p_Nproc(myid+1)
ALLOCATE(paracal(2,p_Nmax))
paracal(:,:)=p_paracal(:,:,myid+1)

!store old Lk
index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Lk_old(index)=Glob_Lk(k,bk)
  ENDDO
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================
!df*p check: df*p<0
!================================

df_p=0.0_dp
IF(Nproc/=0)THEN
DO i=1,Nproc
  j=paracal(1,i)
  k=paracal(2,i)+(j-1)*Glob_NLk
  df_p=df_p+p_df1(k)*p_p_each(i)
ENDDO
ENDIF

temp=df_p
CALL MPI_ALLREDUCE(temp,df_p,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)
IF(df_p>=0.0_dp)PAUSE'roundoff problem in step search'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!calculate minimum step allowed
!===============================

test=0.0_dp
IF(Nproc/=0)THEN 
DO i=1,Nproc
  temp=dabs(p_p_each(i))/max(dabs(p_x1_each(i)),1.0_dp)
  IF(temp>test)THEN
    test=temp
  ENDIF
ENDDO
ENDIF

temp=test
CALL MPI_ALLREDUCE(temp,test,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,Glob_MPIerr)
alpha_min=TOL_dx/test  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!backtracking initialization
!===========================

alpha1=1.0_dp
IF(Nproc/=0)THEN  
DO i=1,Nproc
  p_x2_each(i)=p_x1_each(i)+alpha1*p_p_each(i)
ENDDO
ENDIF

CALL MPI_ALLGATHERV(p_x2_each,Nproc,MPI_DOUBLE_PRECISION,x2,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Glob_Lk(k,bk)=x2(index)
  ENDDO
ENDDO
CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f2,C2,ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!backtracking procedure
!===============================
step_loop:DO WHILE(f2>(f1+d1*alpha1*df_p)) !sufficient decrease condition

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
    
!=============================
!convergence on alpha: if alpha<alpha_min then return
!=============================
IF((alpha1<alpha_min).OR.(ERR/=0))THEN

IF(Nproc/=0)THEN
DO i=1,Nproc
  p_x2_each(i)=p_x1_each(i)
ENDDO
ENDIF

index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Glob_Lk(k,bk)=Lk_old(index)
  ENDDO
ENDDO
CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f2,C2,ERR,1)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!backtracking
!==========================

IF(alpha1==1.0_dp)THEN
  alpha2=-df_p/(2.0_dp*(f2-f1-df_p))
ELSE
  g1=f2-f1-alpha1*df_p
  g2=ftemp-f1-alpha2*df_p
  a=(g1/alpha1**2-g2/alpha2**2)/(alpha1-alpha2)
  b=(-alpha2*g1/alpha1**2+alpha1*g2/alpha2**2)/(alpha1-alpha2)
  
  IF(a==0.0_dp)THEN
    alpha2=-df_p/(2.0_dp*b)
  ELSE 
    temp=b*b-3.0_dp*a*df_p
    IF(temp<0.0_dp)THEN
      alpha2=alpha1*0.5_dp
    ELSEIF(b<=0.0_dp)THEN
      alpha2=(-b+dsqrt(temp))/(3.0_dp*a)
    ELSE
      alpha2=-df_p/(b+dsqrt(temp))
    ENDIF
  ENDIF
  
  !alpha<=0.5 alpha1
  IF(alpha2>(alpha1*0.5_dp))THEN
    alpha2=alpha1*0.5_dp
  ENDIF
  
ENDIF

!update
temp=alpha1
alpha1=max(alpha2,alpha1*0.1_dp) !alpha>=0.1 alpha1
alpha2=temp
ftemp=f2   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!recalculate new value
!==============================

IF(Nproc/=0)THEN    
DO i=1,Nproc
  p_x2_each(i)=p_x1_each(i)+alpha1*p_p_each(i)
ENDDO
ENDIF

CALL MPI_ALLGATHERV(p_x2_each,Nproc,MPI_DOUBLE_PRECISION,x2,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Glob_Lk(k,bk)=x2(index)
  ENDDO
ENDDO
CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f2,C2,ERR,1)

ENDDO step_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!check overlap
DO bk=bk_cal(1),bk_cal(2)
  CALL Skl_check_MPIfun(Nbasis,bk,ERR)
  IF(ERR/=0)EXIT
ENDDO

IF(ERR/=0)THEN

IF(Nproc/=0)THEN
DO i=1,Nproc
  p_x2_each(i)=p_x1_each(i)
ENDDO
ENDIF

index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Glob_Lk(k,bk)=Lk_old(index)
  ENDDO
ENDDO
CALL Eval_MPIfun(Nbasis,Nbs,bk_cal,f2,C2,ERR,1)
IF(ERR/=0)PAUSE'error in Lk_stepsrch_MPIfun'

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_stepsrch_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine for writing parameters into parafile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_parafile(write_case)
!===================================================
!save nonlinear parameters into parafile1.txt and parafile2.txt
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::write_case

INTEGER::i,j,k
CHARACTER(100)::Rst

SELECT CASE(write_case)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
    
OPEN(unit=3,file='parafile1.txt')
WRITE(3,*)'========================================================'
WRITE(3,*)'particle_system:  ',GLob_particle_system
WRITE(3,*)'opt_angular_momentum:',GLob_LM(1:2)
WRITE(3,*)'opt_energy_level:',Glob_energy_level
WRITE(3,*)'present_basis_num:',Glob_Nbasis_reach
WRITE(3,*)'opt_basis_reached:',Glob_opt_reach
WRITE(3,*)'opt_energy_reached:',Glob_E_reach
WRITE(3,*)'========================================================'

IF(Glob_basis_form==0)THEN
  Rst="(????????f26.16)"
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(3,Rst)Glob_Lk(1:Glob_NLk,i)
  ENDDO
ELSEIF(Glob_basis_form==1)THEN
  Rst="(????????I8,????????f26.16)"
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")1
  WRITE(Rst(13:20),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(3,Rst)Glob_mk(i),Glob_Lk(1:Glob_NLk,i)
  ENDDO
ENDIF

CLOSE(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)

OPEN(unit=4,file='parafile2.txt')
WRITE(4,*)'========================================================'
WRITE(4,*)'particle_system:  ',GLob_particle_system
WRITE(4,*)'opt_angular_momentum:',GLob_LM(1:2)
WRITE(4,*)'opt_energy_level:',Glob_energy_level
WRITE(4,*)'present_basis_num:',Glob_Nbasis_reach
WRITE(4,*)'opt_basis_reached:',Glob_opt_reach
WRITE(4,*)'opt_energy_reached:',Glob_E_reach
WRITE(4,*)'========================================================'

IF(Glob_basis_form==0)THEN
  Rst="(????????f26.16)"
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(4,Rst)Glob_Lk(1:Glob_NLk,i)
  ENDDO
ELSEIF(Glob_basis_form==1)THEN
  Rst="(????????I8,????????f26.16)"
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")1
  WRITE(Rst(13:20),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(4,Rst)Glob_mk(i),Glob_Lk(1:Glob_NLk,i)
  ENDDO
ENDIF

CLOSE(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SELECT

RETURN
END SUBROUTINE write_parafile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE gvmopt
