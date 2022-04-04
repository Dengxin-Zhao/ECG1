MODULE gvmopt
!===========================================================
!this module contians subroutines for BFGs optimization 
!with parallel approach
!===========================================================
USE MPI
USE auxfun
USE globvars
USE MPImat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!private variables and matrix used in gvmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Nproc
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_displs
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:,:,:)::p_para_cal
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
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dx       !x_(i+1)-x_i
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_df1      !df_(i+1)
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_df2      !df_i
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_ddf      !df_(i+1)-df_i
!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_H_ddf    !Hessin*(df_i+1-df_i)
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_uterm    !additional term of BFGs 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate and deallocate variables used in gvmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE pvars_gvmopt(Nbs,para_each,case_num)
!================================
!allocate and deallocate private varibles used in gvmopt.f90
!input:
!     Nbs: number of basis in each batch
!para_each: number of parameters in each basis
!  Nbasis: present basis number
!case_num: 1: allocate the variables
!          0: deallocate the variables
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs
INTEGER,INTENT(IN)::para_each
INTEGER,INTENT(IN)::case_num

INTEGER::i,j,k,L
INTEGER::para_num
INTEGER::N_last,N_each
INTEGER::index1,index2

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
SELECT CASE(case_num)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1) !allocate

para_num=para_each*Nbs               !numbers of parameters in each batch
N_last=mod(para_num,Glob_numprocs)
N_each=(para_num-N_last)/Glob_numprocs
p_Nmax=N_each+1
IF(N_last==0)THEN
p_Nmax=p_Nmax-1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(p_x1_each(p_Nmax),p_x2_each(p_Nmax))
ALLOCATE(p_dx_each(p_Nmax),p_p_each(p_Nmax))
ALLOCATE(p_H_ddf_each(p_Nmax),p_uterm_each(p_Nmax))
ALLOCATE(p_hessin_each(para_num,p_Nmax))
!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(p_dx(para_num))
ALLOCATE(p_df1(para_num))
ALLOCATE(p_df2(para_num))
ALLOCATE(p_ddf(para_num))
ALLOCATE(p_H_ddf(para_num))
ALLOCATE(p_uterm(para_num))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(p_Nproc(Glob_numprocs))
ALLOCATE(p_displs(Glob_numprocs))
ALLOCATE(p_para_cal(2,p_Nmax,Glob_numprocs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!location of parameters calculated by each process
p_para_cal=0
index1=1 !process sign
index2=0 !parameter number sign
DO i=1,Nbs
  DO j=1,para_each
    index2=index2+1
    p_para_cal(1,index2,index1)=i
    p_para_cal(2,index2,index1)=j
    IF((index2==(N_each+1)).AND.(index1<=N_last))THEN
      index1=index1+1
      index2=0
    ELSEIF((index2==N_each).AND.(index1>N_last))THEN
      index1=index1+1
      index2=0
    ENDIF     
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)!deallocate
    

DEALLOCATE(p_x1_each,p_x2_each)
DEALLOCATE(p_dx_each,p_p_each)
DEALLOCATE(p_H_ddf_each,p_uterm_each)
DEALLOCATE(p_hessin_each)
!!!!!!!!!!!!!!!!!!!!!!
DEALLOCATE(p_dx)
DEALLOCATE(p_df1)
DEALLOCATE(p_df2)
DEALLOCATE(p_ddf)
DEALLOCATE(p_H_ddf)
DEALLOCATE(p_uterm)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DEALLOCATE(p_Nproc)
DEALLOCATE(p_displs)
DEALLOCATE(p_para_cal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE pvars_gvmopt
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    
SUBROUTINE gvm_opt(Nbs,Nbasis,ITmax,Nrounds)
!========================================
!this subroutine optimize the nonlinear parameters
!of each batch one by one with BFGs algrithm 
!and write the optimization procedure in log.txt
!record every optimized result into parafile1.txt and parafile2.txt
!======================
!input:
!      Nbs: optimize Nbs basis every time
!   Nbasis: the basis number to be optimized now
!           (Nbasis basis counting from 1 in Glob_Lk)
!    ITmax: max iterations for every batch
!  Nrounds: the rounds of BFGs optimizations of present N basis
!output:
!  eignvalues are stored in log.txt and nonlinear parameters
!  are stored in savefile.txt 
!========================================   
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,Nbasis
INTEGER,INTENT(IN)::ITmax,Nrounds

INTEGER::i,j,k,L,myid,ERR
INTEGER::nr,nb,it,nb_start,count,write_index
INTEGER::Nbatches,para_num,bk_cal(2)
INTEGER,ALLOCATABLE,DIMENSION(:,:)::basis_index

REAL(dp)::Eval,Ci(Nbasis),temp

CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!initialization
!===============================

CALL pvars_MPImat(Nbasis,1)
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
IF(ERR/=0)PAUSE'eigen solution error before gvm_opt'

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'gvm optimization with BFGs start:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'batch basis number:',Nbs
  WRITE(2,*)'================================================='
  WRITE(2,*)'start energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,2A10,A14,A11)')'time','count','round','batch_num','energy'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

batch_loop:DO nb=nb_start,Nbatches 
bk_cal(:)=basis_index(:,nb)
count=count+1

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
  
!==========================
!optimization of mk
!==========================
IF(Glob_basis_form==1)THEN
  CALL mk_gvm_MPIfun(Nbs,bk_cal,Nbasis,Eval,Ci)
ENDIF

!===========================
!optimization of Lk with BFGs algrithm
!===========================
CALL Lk_gvm_MPIfun(Nbs,bk_cal,Nbasis,ITmax,Eval,Ci)

Glob_opt_reach=bk_cal(2)
Glob_E_reach=Eval

!!!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!write present energy in log.txt
CALL time_cal_fun(time_st)
WRITE(2,'(A14,2I10,I14,f28.16)')time_st,count,nr,nb,Eval

!write parameters into parafile.txt
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO batch_loop
ENDDO optimize_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'gvm optimization with BFGs finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,0)
CALL pvars_gvmopt(Nbs,Glob_NLk,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE gvm_opt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_gvm_MPIfun(Nbs,bk_cal,Nbasis,ITmax,f1,C1)
!========================================
!this subroutine optimize the nonlinear parameters
!======================
!input:
!      Nbs: optimize Nbs basis every time
!   bk_cal: basis index bk_cal(1)->bk_cal(2) 
!   Nbasis: the basis number to be optimized now
!    ITmax: max iterations for every batch
!       f1: the initial eigenvalue
!       C1: initial eigenvector
!inout:
!       f1: optimized eigenvalue
!       C1: optimized eigenvector
!========================
!this subroutine should 
!CALL pvars_MPImat(Nbasis,1)
!CALL pvars_gvmopt(Nbs,Glob_NLk,1)
!before using it.
!========================================   
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,Nbasis,ITmax
INTEGER,INTENT(IN)::bk_cal(2)
REAL(dp),INTENT(INOUT)::f1,C1(Nbasis)

REAL(dp),PARAMETER::TOL_dx=FIVE*EPS   !convergence criterion on dx
REAL(dp),PARAMETER::TOL_df=FIVE*EPS   !convergence criterion on df

REAL(dp)::stpmax                      !max step allowed
REAL(dp)::ddf_dx                      !(df_i+1-df_i)(x_(i+1)-x_i)
REAL(dp)::ddf_H_ddf                   !(df_i+1-df_i)Hessin(df_i+1-df_i)
REAL(dp)::sum_ddf,sum_dx              !(df_i+1-df_i)^2,(x_(i+1)-x_i)^2
REAL(dp)::temp,test,f2,C2(Nbasis) 

INTEGER::i,j,k,L,myid,ERR
INTEGER::it,bi,bj,Nproc,para_num
INTEGER,ALLOCATABLE,DIMENSION(:,:)::para_cal              

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(ITmax==0)RETURN
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

!distribution of process
Nproc=p_Nproc(myid+1)      
ALLOCATE(para_cal(2,p_Nmax))         
para_cal(:,:)=p_para_cal(:,:,myid+1)

!bi->bj basis
bi=bk_cal(1)
bj=bk_cal(2)
para_num=Nbs*Glob_NLk

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==============================
!initialization of BFGs
!==============================

!f1 and df1
!CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,f1,C1,ERR,1)
CALL dE_dLk_MPIfun(Nbs,bk_cal,Nbasis,f1,C1,p_df1)

temp=ZERO
p_hessin_each(:,:)=ZERO
IF(Nproc/=0)THEN    
DO i=1,Nproc
    
  j=bi+para_cal(1,i)-1
  k=para_cal(2,i)
  p_x1_each(i)=Glob_Lk(k,j)          !x1
  temp=temp+p_x1_each(i)**2
  
  j=para_cal(1,i)
  k=para_cal(2,i)+(j-1)*Glob_NLk
  p_hessin_each(k,i)=ONE             !hessin=I
  
ENDDO
ENDIF

test=temp
CALL MPI_ALLREDUCE(test,temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

!max step allowed (avoid unreasonable steps)
stpmax=100.E0_dp*max(dsqrt(temp),dble(para_num))

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

  temp=ZERO
  DO k=1,para_num
    temp=temp-p_hessin_each(k,i)*p_df1(k)
  ENDDO
  p_p_each(i)=temp

ENDDO
ENDIF

!================================
!scale direction p if it is too big
!================================
temp=ZERO
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
CALL Lk_stepsrch_MPIfun(Nbs,bk_cal,Nbasis,f1,C1,f2,C2,TOL_dx)

!dx
p_dx_each=p_x2_each-p_x1_each 

!=============================
!dx check: max(dx)<TOL_dx -> exit BFGs_loop 
!=============================
test=ZERO
IF(Nproc/=0)THEN
DO i=1,Nproc

  temp=dabs(p_dx_each(i))/max(dabs(p_x2_each(i)),ONE)
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
CALL dE_dLk_MPIfun(Nbs,bk_cal,Nbasis,f2,C2,p_df2)

!=============================
!df check: max(df2)<TOL_df -> exit BFGs_loop
!=============================
test=ZERO
IF(Nproc/=0)THEN 
DO i=1,Nproc

  j=para_cal(1,i)
  k=para_cal(2,i)+Glob_NLk*(j-1)
  temp=dabs(p_df2(k))*max(dabs(p_x2_each(i)),ONE)/max(f2,ONE)
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

ddf_dx=ZERO
ddf_H_ddf=ZERO  
sum_ddf=ZERO
sum_dx=ZERO

IF(Nproc/=0)THEN
DO i=1,Nproc
    
  !Hessin*(df_i+1-df_i)
  temp=ZERO  
  DO k=1,para_num
    temp=temp+p_hessin_each(k,i)*p_ddf(k)
  ENDDO
  p_H_ddf_each(i)=temp

  j=para_cal(1,i)
  k=para_cal(2,i)+Glob_NLk*(j-1)
  
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

  j=para_cal(1,i)
  k=para_cal(2,i)+(j-1)*Glob_NLk
  DO j=1,para_num
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

CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,f1,C1,ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_gvm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!step search function for parameter Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Lk_stepsrch_MPIfun(Nbs,bk_cal,Nbasis,f1,C1,f2,C2,TOL_dx)
!==============================
!This subroutine search the step of parameters. Be aware 
!of the change of Glob_Lk Glob_Skl and Glob_Hkl matrix
!when calling E_rows_MPIfun, if no new steps are found
!they should keep their old values
!==============================
!input:
!    Nbs: the number of basis of present batch 
! bk_cal: basis index: bk_cal(1)->bk_cal(2)
! Nbasis: present basis number
!     f1: old objective value
!     C1: old eigen vector corresponding to f1 
!output:
!     f2: new objective value
!     C2: new eigen vector corresponding to f2 
! TOL_dx: tolerance of dx 
!==============================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,Nbasis
INTEGER,INTENT(IN)::bk_cal(Nbs)
REAL(dp),INTENT(IN)::TOL_dx
REAL(dp),INTENT(IN)::f1,C1(Nbasis)
REAL(dp),INTENT(OUT)::f2,C2(Nbasis)
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
REAL(dp),PARAMETER::d1=1.d-4  !sufficient decrease coef
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

INTEGER::i,j,k,L,myid,ERR
INTEGER::bk,Nproc,index
INTEGER,ALLOCATABLE,DIMENSION(:,:)::para_cal

REAL(dp)::alpha_min             
REAL(dp)::df_p,alpha1,alpha2,f_temp
REAL(dp)::g1,g2,a,b
REAL(dp)::test,temp

REAL(dp)::Lk_old(Glob_NLk*Nbs)
REAL(dp)::x2(Glob_NLk*Nbs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

!distribution of process
Nproc=p_Nproc(myid+1)
ALLOCATE(para_cal(2,p_Nmax))
para_cal(:,:)=p_para_cal(:,:,myid+1)

!store old Lk
index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Lk_old(index)=Glob_Lk(k,bk)
  ENDDO
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!================================
!df*p check: df*p<0
!================================

df_p=ZERO
IF(Nproc/=0)THEN
DO i=1,Nproc

  j=para_cal(1,i)
  k=para_cal(2,i)+(j-1)*Glob_NLk
  df_p=df_p+p_df1(k)*p_p_each(i)

ENDDO
ENDIF

temp=df_p
CALL MPI_ALLREDUCE(temp,df_p,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

IF(df_p>=ZERO)PAUSE'roundoff problem in step search'

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!calculate minimum step allowed
!===============================

test=ZERO
IF(Nproc/=0)THEN 
DO i=1,Nproc
  
  temp=dabs(p_p_each(i))/max(dabs(p_x1_each(i)),ONE)
  IF(temp>test)THEN
    test=temp
  ENDIF

ENDDO
ENDIF

temp=test
CALL MPI_ALLREDUCE(temp,test,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,Glob_MPIerr)

alpha_min=TOL_dx/test  

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!backtracking initialization
!===========================

alpha1=ONE
IF(Nproc/=0)THEN  
DO i=1,Nproc
    
  p_x2_each(i)=p_x1_each(i)+alpha1*p_p_each(i)

ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(p_x2_each,Nproc,MPI_DOUBLE_PRECISION,x2,p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

index=0
DO bk=bk_cal(1),bk_cal(2)
  DO k=1,Glob_NLk
    index=index+1
    Glob_Lk(k,bk)=x2(index)
  ENDDO
ENDDO
  
!check overlap
DO bk=bk_cal(1),bk_cal(2)
  CALL overlap_check_MPIfun(bk,Nbasis,ERR)
  IF(ERR/=0)EXIT
ENDDO
IF(ERR/=0)THEN
  f2=1.d5
ELSE    
  CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,f2,C2,ERR,1)
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,f2,C2,ERR,1)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!backtracking
!==========================

IF(alpha1==ONE)THEN
  alpha2=-df_p/(TWO*(f2-f1-df_p))
ELSE
  g1=f2-f1-alpha1*df_p
  g2=f_temp-f1-alpha2*df_p
  a=(g1/alpha1**2-g2/alpha2**2)/(alpha1-alpha2)
  b=(-alpha2*g1/alpha1**2+alpha1*g2/alpha2**2)/(alpha1-alpha2)
  
  IF(a==ZERO)THEN
    alpha2=-df_p/(TWO*b)
  ELSE 
    temp=b*b-THREE*a*df_p
    IF(temp<0)THEN
      alpha2=alpha1/TWO
    ELSEIF(b<=0)THEN
      alpha2=(-b+dsqrt(temp))/(THREE*a)
    ELSE
      alpha2=-df_p/(b+dsqrt(temp))
    ENDIF
  ENDIF
  
  !alpha <= 0.5 alpha1
  IF(alpha2>(alpha1/TWO))THEN
    alpha2=alpha1/TWO
  ENDIF
  
ENDIF

!update
temp=alpha1
alpha1=max(alpha2,alpha1/TEN) !alpha>=0.1 alpha1
alpha2=temp
f_temp=f2   

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

!check overlap
DO bk=bk_cal(1),bk_cal(2)
  CALL overlap_check_MPIfun(bk,Nbasis,ERR)
  IF(ERR/=0)THEN
    f2=1.d5
    CYCLE step_loop
  ENDIF
ENDDO

CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,f2,C2,ERR,1)
IF(ERR/=0)f2=1.d5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO step_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Lk_stepsrch_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine for optimization of mk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mk_gvm_MPIfun(Nbs,bk_cal,Nbasis,Eval,Ci)
!================================
!optimization of parameter mk
!strategy:
!calculate the eigen values for each mk(1->Glob_Nparticle)
!and find the one making the value lowest 
!input:
!   Nbs: number of batch
!bk_cal: basis index: bk_cal(1)->bk_cal(2)
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
!before using it
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbs,Nbasis
INTEGER,INTENT(IN)::bk_cal(2)
REAL(dp),INTENT(INOUT)::Eval,Ci(Nbasis)

INTEGER::i,j,k,L,myid,ERR
INTEGER::bi,bk,ip,basis_cal(2)
INTEGER::mk_best,mk_old
REAL(dp)::E_best,C_best(Nbasis)
REAL(dp)::E_temp,C_temp(Nbasis)

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!CALL E_rows_MPIfun(Nbs,bk_cal,Nbasis,Eval,Ci,ERR,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!optimization
!============================

E_best=Eval

bk_loop:DO bk=bk_cal(1),bk_cal(2)

basis_cal(1)=bk
basis_cal(2)=bk
mk_best=Glob_mk(1,bk)
mk_old=mk_best

ip_loop:DO ip=1,Glob_Np

Glob_mk(1,bk)=ip
CALL overlap_check_MPIfun(bk,Nbasis,ERR)

IF(ERR==0)THEN

  CALL E_rows_MPIfun(1,basis_cal,Nbasis,E_temp,C_temp,ERR,1)
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
CALL E_rows_MPIfun(1,basis_cal,Nbasis,E_temp,C_temp,ERR,1)

IF(E_temp<Eval.AND.ERR==0)THEN
  Eval=E_temp
  Ci(:)=C_temp(:)
ELSE
  Glob_mk(1,bk)=mk_old
  CALL E_rows_MPIfun(1,basis_cal,Nbasis,Eval,Ci,ERR,1)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDDO bk_loop

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE mk_gvm_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine for writing parameters into parafile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_parafile(write_case)
!================================OK
!save nonlinear parameters into savefile1.txt and savefile2.txt
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::write_case

INTEGER::i,j,k
CHARACTER(100)::Rst

SELECT CASE(write_case)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
    
OPEN(unit=3,file='parafile1.txt')
WRITE(3,*)'========================================================'
WRITE(3,*)'particle_system:  ',GLob_particle_system
WRITE(3,*)'opt_angular_momentum:',GLob_angular_momentum(1:2)
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
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")Glob_Nmk
  WRITE(Rst(13:20),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(3,Rst)Glob_mk(1:Glob_Nmk,i),Glob_Lk(1:Glob_NLk,i)
  ENDDO
ENDIF

CLOSE(3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)

OPEN(unit=4,file='parafile2.txt')
WRITE(4,*)'========================================================'
WRITE(4,*)'particle_system:  ',GLob_particle_system
WRITE(4,*)'opt_angular_momentum:',GLob_angular_momentum(1:2)
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
  WRITE(Rst(2:9),fmt="(TL1,I8.8)")Glob_Nmk
  WRITE(Rst(13:20),fmt="(TL1,I8.8)")Glob_NLk
  DO i=1,Glob_Nbasis_reach
    WRITE(4,Rst)Glob_mk(1:Glob_Nmk,i),Glob_Lk(1:Glob_NLk,i)
  ENDDO
ENDIF

CLOSE(4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
END SELECT

RETURN
END SUBROUTINE write_parafile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE gvmopt
