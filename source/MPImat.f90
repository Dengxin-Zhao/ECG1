MODULE MPImat
!==============================================
!this module contains subroutines computing 
!H and S matrix, dH and dS matrix with MPI parallel approach
!==============================================
USE MPI
USE globvars
USE matelem
USE eigen
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!private variables and matrix used in MPImat.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
INTEGER,PRIVATE::p_Nmax
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Nproc
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_displs
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_bl_cal
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate and deallocate variables used in gvmopt.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pvars_MPImat(Nbasis,case_num)
!===================================================
!allocate the MPI allocation variables used in forming 
!H,S matrix and gradient matrix.
!===================================================
!input:
!  Nbasis: present basis number
!case_num: select case
!output:
!  p_Nmax: max column numbers minipulated by each basis
! p_Nproc: column numbers minipulated by each process
!p_displs: displacement when gatherv 
!p_bl_cal: columns calculated by each pocess
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis
INTEGER,INTENT(IN)::case_num

INTEGER::i,j,k,L
INTEGER::index1,index2
INTEGER::N_last,N_each

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
SELECT CASE(case_num)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1) 
!===========================
!allocate Nbasis to each each process
!===========================

N_last=mod(Nbasis,Glob_numprocs)
N_each=(Nbasis-N_last)/Glob_numprocs
p_Nmax=N_each+1
IF(N_last==0)THEN
p_Nmax=p_Nmax-1
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ALLOCATE(p_Nproc(Glob_numprocs))
ALLOCATE(p_displs(Glob_numprocs))
ALLOCATE(p_bl_cal(p_Nmax,Glob_numprocs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!numbers of parameters manipulated by each procs
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

!basis calculated by each process
p_bl_cal=0
index1=1 !process sign
index2=0 !basis sign
DO i=1,Nbasis
  index2=index2+1
  p_bl_cal(index2,index1)=i
  IF((index2==(N_each+1)).AND.(index1<=N_last))THEN
    index1=index1+1
    index2=0
  ELSEIF((index2==N_each).AND.(index1>N_last))THEN
    index1=index1+1
    index2=0
  ENDIF     
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)  !deallocate

DEALLOCATE(p_displs)
DEALLOCATE(p_Nproc)
DEALLOCATE(p_bl_cal)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE pvars_MPImat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculaiton of Hamilton and overlap matrix and  diagonalization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Eval_MPIfun(Nbasis,Nbs,bk_cal,Eval,Ci,ERR,case_num)
!===================================================
!recalculate specified rows (and columns) of H and S matrix 
!bk_cal(1)->bk_cal(2)
!===================================================
!input:
! Nbasis: present numbers of basis
!     Nb: numbers of basis to be recalculated
! bk_cal: basis index from bk_cal(1)->bk_cal(2) to recalculated
!output:
!   Eval: eigenvalue energy
!     Ci: the eigenvector corresbonding to Eval
!  ERR: 
!    0: within overlap threshold
!    1: exceed overlap threshold
!case_num: 
!    0: one_diag_MPIfun
!    1: full_diag_MPIfun
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbs,bk_cal(2),case_num
REAL(dp),INTENT(OUT)::Eval,Ci(Nbasis)
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,myid
INTEGER::bk,bl,Nproc,bl_cal(p_Nmax)
REAL(dp)::send_S(p_Nmax),send_H(p_Nmax) 

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

Nproc=p_Nproc(myid+1)
bl_cal(:)=p_bl_cal(:,myid+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!calculation of Hkl and Skl
!==========================

bk_loop:DO bk=bk_cal(1),bk_cal(2)

IF(Nproc/=0)THEN   
DO i=1,Nproc
    
  bl=bl_cal(i)
  IF(Glob_basis_form==0)THEN
    CALL HS_0_fun(bk,bl,send_S(i),send_H(i))
  ELSEIF(Glob_basis_form==1)THEN
    CALL HS_1_fun(bk,bl,send_S(i),send_H(i))
  ENDIF

ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(send_S,Nproc,MPI_DOUBLE_PRECISION,Glob_Skl(:,bk),p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLGATHERV(send_H,Nproc,MPI_DOUBLE_PRECISION,Glob_Hkl(:,bk),p_Nproc,p_displs,&
&MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,Glob_MPIerr)

DO i=1,Nbasis
  Glob_Skl(bk,i)=Glob_Skl(i,bk)
  Glob_Hkl(bk,i)=Glob_Hkl(i,bk) 
ENDDO

ENDDO bk_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!calculatioin of eigenvalue
!==========================

SELECT CASE(case_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)

CALL one_diag_MPIfun(Nbasis,bk_cal(1),Eval,ERR)
Ci(:)=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
  
CALL full_diag_MPIfun(Nbasis,Eval,Ci,ERR)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE Eval_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap check function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_check_MPIfun(Nbasis,bk,ERR)
!===================================================
!check wheather the overlap exceeds the threshold
!===================================================
!input:
!Nbasis: present numbers of basis
!    bk: the bk_th row of overlap S matrix
!output:
!  ERR
!      0: overlap within overlap threshold 
!      1: overlap exceed overlap threshold
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,bk
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,myid
INTEGER::bl,Nproc,bl_cal(p_Nmax)

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

ERR=0
Nproc=p_Nproc(myid+1)
bl_cal(:)=p_bl_cal(:,myid+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!check the overlap 
!============================

IF(Nproc/=0)THEN
DO i=1,Nproc
    
  bl=bl_cal(i)
  IF(Glob_basis_form==0)THEN
    CALL Skl_0_check_fun(bk,bl,ERR)
  ELSEIF(Glob_basis_form==1)THEN
     CALL Skl_1_check_fun(bk,bl,ERR)
  ENDIF

  IF(ERR/=0)EXIT

ENDDO
ENDIF

k=ERR
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_ALLREDUCE(k,ERR,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE Skl_check_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient of energy with respect to Lk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dE_dLk_MPIfun(Nbasis,Nbk,bk_cal,Eval,Ci,dE)
!===================================================
!calculate the energy gradient with respect to Lk
!===================================================
!input:
!Nbasis: present basis nunber
!    Nb: basis number
!bk_cal: basis index from bk_cal(1)->bk_cal(2)
!  Eval: present eigenvalue 
!    Ci: eigenvector corresbonding to Eval
!output:
!    dE: gradient of energy with respect to Lk 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,Nbk,bk_cal(2)
REAL(dp),INTENT(IN)::Eval,Ci(Nbasis)
REAL(dp),INTENT(OUT)::dE(Nbk*Glob_NLk)

INTEGER::i,j,k,myid,index
INTEGER::bk,bl,Nproc,bl_cal(p_Nmax)
REAL(dp)::dS_dLk(Glob_NLk),dH_dLk(Glob_NLk)
REAL(dp)::send_dE(Glob_NLk),recv_dE(Glob_NLk)

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

Nproc=p_Nproc(myid+1)
bl_cal(:)=p_bl_cal(:,myid+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=========================
!calculation of gradient of E 
!=========================

index=0
bk_loop:DO bk=bk_cal(1),bk_cal(2)
    
send_dE(:)=0.0_dp
recv_dE(:)=0.0_dp

IF(Nproc/=0)THEN   
DO i=1,Nproc
    
  bl=bl_cal(i)
  IF(Glob_basis_form==0)THEN
    CALL dHS_dLk_0_fun(bk,bl,dS_dLk,dH_dLk)
  ELSEIF(Glob_basis_form==1)THEN
    CALL dHS_dLk_1_fun(bk,bl,dS_dLk,dH_dLk)
  ENDIF
  
  IF(bk==bl)THEN
    send_dE(:)=send_dE(:)+Ci(bk)*Ci(bk)*(dH_dLk(:)-Eval*dS_dLk(:))
  ELSE
    send_dE(:)=send_dE(:)+2.0_dp*Ci(bk)*Ci(bl)*(dH_dLk(:)-Eval*dS_dLk(:))
  ENDIF
  
ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_ALLREDUCE(send_dE,recv_dE,Glob_NLk,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,Glob_MPIerr)

DO i=1,Glob_NLk
  index=index+1
  dE(index)=recv_dE(i)
ENDDO

ENDDO bk_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE dE_dLk_MPIfun   
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function gij for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
SUBROUTINE gij_0_MPIfun(Nbasis)
!===================================================
!thia subroutine calculates the Correlation function gij
!and store them into gij.txt
!===================================================
!input:
!  Nbasis: present basis number
!output:
! correlation function gij stored in gij.txt
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis

INTEGER::i,j,k,L,index,myid,ERR,count
INTEGER::ii,jj,bk,bl,Nproc,bk_cal(2)
INTEGER,ALLOCATABLE::bl_cal(:)

REAL(dp)::Eval,Ci(Nbasis)
REAL(dp)::R,temp
REAL(dp)::gij(GLob_Nparticle*(Glob_Nparticle+1)/2)

CHARACTER(14)::time_st
CHARACTER(100)::Rst
Rst="(????????f26.15)"
WRITE(Rst(2:9),fmt="(TL1,I8.8)")1+GLob_Nparticle*(Glob_Nparticle+1)/2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

Nproc=p_Nproc(myid+1)
ALLOCATE(bl_cal(p_Nmax))
bl_cal(:)=p_bl_cal(:,myid+1)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!calculate the eigen value and eigen vector
bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  OPEN(unit=6,file='gij.txt')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of correaltion function starts:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'present energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,2A10,A28)')'time','count','R','status'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

count=0
node_loop:DO R=Glob_gij_R(1),Glob_gij_R(3),Glob_gij_R(2)
count=count+1

index=0
ii_loop:DO ii=1,Glob_Nparticle
jj_loop:DO jj=ii,Glob_Nparticle
index=index+1
gij(index)=0.0_dp
IF(Glob_gij_onoff(index)/=0)THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bk_loop:DO bk=1,Nbasis

IF(Nproc/=0)THEN   
DO i=1,Nproc
    
  bl=bl_cal(i)

  CALL gij_kl_0_fun(bk,bl,ii,jj,R,temp)
  
  gij(index)=gij(index)+4.0_dp*PI*R*R*Ci(bk)*Ci(bl)*temp

ENDDO
ENDIF

ENDDO bk_loop

!sum all partial density and reduce to Glob_root
temp=gij(index)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_REDUCE(temp,gij(index),1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDIF
ENDDO jj_loop
ENDDO ii_loop

!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)
  
!write present status in log file
WRITE(2,'(A14,I10,f26.15,A11)')time_st,count,R,'OK'

WRITE(6,Rst)R,gij(1:index)

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO node_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of gij finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  CLOSE(2)
  CLOSE(6)
ENDIF
CALL pvars_MPImat(Nbasis,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,GLob_MPIerr)

RETURN
END SUBROUTINE gij_0_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function gij for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gij_1_MPIfun(Nbasis)
!===================================================
!this subroutine calculates the correlation function
!for basis with L=1,M=0
!===================================================
!input:
!  Nbasis: present basis number
!output:
! correlation function stored in gij.txt
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis

INTEGER::i,j,k,L,index,myid,ERR,count
INTEGER::ii,jj,bk,bl,Nproc,bk_cal(2)
INTEGER,ALLOCATABLE::bl_cal(:)
  
REAL(dp)::Eval,Ci(Nbasis)
REAL(dp)::Rxy,Rz,temp
REAL(dp)::gij(Glob_Nparticle*(Glob_Nparticle+1)/2)

CHARACTER(14)::time_st
CHARACTER(100)::Rst
Rst="(????????f26.15)"
WRITE(Rst(2:9),fmt="(TL1,I8.8)")2+GLob_Nparticle*(Glob_Nparticle+1)/2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL pvars_MPImat(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)

Nproc=p_Nproc(myid+1)
ALLOCATE(bl_cal(p_Nmax))
bl_cal(:)=p_bl_cal(:,myid+1)

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!calculate the eigen value and eigen vector
bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  OPEN(unit=6,file='gij.txt')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of correaltion function starts:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'present energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,2A10,A26,A30)')'time','count','Rxy','Rz','status'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

count=0
node_loop1:DO Rxy=Glob_gij_R(1),Glob_gij_R(3),Glob_gij_R(2)
node_loop2:DO Rz=-Glob_gij_R(3)/2.0_dp,Glob_gij_R(3)/2.0_dp,Glob_gij_R(2)
count=count+1

index=0
ii_loop:DO ii=1,Glob_Nparticle
jj_loop:DO jj=ii,Glob_Nparticle
index=index+1
gij(index)=0.0_dp

IF(Glob_gij_onoff(index)/=0)THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

bk_loop:DO bk=1,Nbasis

IF(Nproc/=0)THEN   
DO i=1,Nproc
    
  bl=bl_cal(i)

  CALL gij_kl_1_fun(bk,bl,ii,jj,Rxy,Rz,temp)
  
  gij(index)=gij(index)+2.0_dp*PI*Rxy*Ci(bk)*Ci(bl)*temp

ENDDO
ENDIF

ENDDO bk_loop

!sum all part and reduce to Glob_root
temp=gij(index)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_REDUCE(temp,gij(index),1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDIF
ENDDO jj_loop
ENDDO ii_loop

!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)
  
!write present status in log file
WRITE(2,'(A14,I10,2f26.15,A10)')time_st,count,Rxy,Rz,'OK'

WRITE(6,Rst)Rxy,Rz,gij(1:index)

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO node_loop2
ENDDO node_loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of gij finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  CLOSE(2)
  CLOSE(6)
ENDIF
CALL pvars_MPImat(Nbasis,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,GLob_MPIerr)

RETURN
END SUBROUTINE gij_1_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!interparticle distance <rij> for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_0_MPIfun(Nbasis,v)
!===================================================
!this subroutine calculates the expectation value
!of interparticle distance for basis with L=0,M=0
!===================================================
!input:
!  Nbasis: present basis number
!       v: power of interparticle distance
!output:
!  the expetation value of rij stored in log.txt
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis
REAL(dp)::v

INTEGER::i,j,k,L,ii,jj,index,myid,ERR
INTEGER::bk,bl,Nproc,bk_cal(2)
INTEGER,ALLOCATABLE::bl_cal(:)
    
REAL(dp)::Eval,Ci(Nbasis)
REAL(dp)::rij,temp
  
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
CALL pvars_MPImat(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
  
Nproc=p_Nproc(myid+1)
ALLOCATE(bl_cal(p_Nmax))
bl_cal(:)=p_bl_cal(:,myid+1)
  
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!calculate the eigen value and eigen vector
bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)
  
IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of <rij> starts:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'present energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,*)'power of rij:',v
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,3A10,A11)')'time','count','pk','pl','<rij>'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

index=0
ii_loop:DO ii=1,Glob_Nparticle
jj_loop:DO jj=ii,Glob_Nparticle
index=index+1
rij=0.0_dp
IF(Glob_rij_onoff(index)/=0)THEN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

bk_loop:DO bk=1,Nbasis
  
IF(Nproc/=0)THEN   
DO i=1,Nproc
      
  bl=bl_cal(i)
  
  CALL rij_kl_0_fun(bk,bl,ii,jj,v,temp)
    
  rij=rij+Ci(bk)*Ci(bl)*temp
  
ENDDO
ENDIF
  
ENDDO bk_loop
  
!sum all partial density and reduce to Glob_root
temp=rij
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_REDUCE(temp,rij,1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

ENDIF
!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!
IF(myid==Glob_root)THEN 

!time calculate  
CALL time_cal_fun(time_st)

!write present status in log file
WRITE(2,'(A14,3I10,f26.16)')time_st,index,ii,jj,rij
  
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
ENDDO jj_loop
ENDDO ii_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of <rij> finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


IF(myid==Glob_root)CLOSE(2)
CALL pvars_MPImat(Nbasis,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE rij_0_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!interparticle distance <rij> for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_1_MPIfun(Nbasis,v)
!===================================================
!this subroutine calculates the expectation value
!of interparticle distance for basis with L=1,M=0
!===================================================
!input:
!  Nbasis: present basis number
!output:
!  the expetation value of rij stored in log.txt
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis
REAL(dp),INTENT(IN)::v
  
INTEGER::i,j,k,L,index,myid,ERR
INTEGER::ii,jj,bk,bl,Nproc,bk_cal(2)
INTEGER,ALLOCATABLE::bl_cal(:)
      
REAL(dp)::Eval,Ci(Nbasis)
REAL(dp)::rij,temp
    
CHARACTER(14)::time_st

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
CALL pvars_MPImat(Nbasis,1)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
    
Nproc=p_Nproc(myid+1)
ALLOCATE(bl_cal(p_Nmax))
bl_cal(:)=p_bl_cal(:,myid+1)
        
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!calculate the eigen value and eigen vector
bk_cal(1)=1
bk_cal(2)=Nbasis
CALL Eval_MPIfun(Nbasis,Nbasis,bk_cal,Eval,Ci,ERR,1)

IF(myid==Glob_root)THEN 
  OPEN(unit=2,file='log.txt',position='append')
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of <rij> starts:'
  WRITE(2,*)'present total basis number:',Nbasis
  WRITE(2,*)'present energy:',Eval
  WRITE(2,*)'================================================='
  WRITE(2,*)'power of rij:',v
  WRITE(2,*)'================================================='
  WRITE(2,'(A14,3A10,A11)')'time','count','pk','pl','<rij>'
ENDIF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
index=0
ii_loop:DO ii=1,Glob_Nparticle
jj_loop:DO jj=ii,Glob_Nparticle
index=index+1
rij=0.0_dp

IF(Glob_rij_onoff(index)/=0)THEN
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

bk_loop:DO bk=1,Nbasis
    
IF(Nproc/=0)THEN   
DO i=1,Nproc
        
  bl=bl_cal(i)
    
  CALL rij_kl_1_fun(bk,bl,ii,jj,v,temp)
      
  rij=rij+Ci(bk)*Ci(bl)*temp
    
  ENDDO
  ENDIF
    
ENDDO bk_loop
    
!sum all partial density and reduce to Glob_root
temp=rij
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_REDUCE(temp,rij,1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
    
!!!!!!!!!!!!!!output!!!!!!!!!!!!!!!!!!!

ENDIF

IF(myid==Glob_root)THEN 
  
CALL time_cal_fun(time_st)
WRITE(2,'(A14,3I10,f26.16)')time_st,index,ii,jj,rij
    
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ENDDO jj_loop
ENDDO ii_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)THEN
  WRITE(2,*)'=================================================' 
  WRITE(2,*)'calculation of <rij> finished'
  WRITE(2,*)'================================================='
  WRITE(2,*)
  CLOSE(2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(myid==Glob_root)CLOSE(2)
CALL pvars_MPImat(Nbasis,0)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

RETURN
END SUBROUTINE rij_1_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE MPImat