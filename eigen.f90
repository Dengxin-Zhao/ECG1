MODULE eigen
!=========================================================
!this module contains subroutines for solving the eigen equation
!=========================================================
USE MPI
USE globvars
!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Ei,p_hk
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_Ci
!!!!!!!!!!!!!!!!!!!!!!!!!
INTEGER,PRIVATE::p_Nmax 
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:)::p_Nproc,p_displs
INTEGER,PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_bk_cal
!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!private variables used in eigen.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pvars_eigen(Nbasis,case_num)
!==========================================
!allocate and deallocate the variables used in eigen.f90
!==========================================
INTEGER,INTENT(IN)::Nbasis,case_num

INTEGER::i,j,k,index1,index2
INTEGER::N_each,N_last

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr) 

SELECT CASE(case_num)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
  
ALLOCATE(p_Ei(Nbasis),p_hk(Nbasis))
ALLOCATE(p_Ci(Nbasis,Nbasis))
ALLOCATE(p_Nproc(Glob_numprocs))
ALLOCATE(p_displs(Glob_numprocs))

!allocate Nbasis to each process
N_last=mod(Nbasis,Glob_numprocs)
N_each=(Nbasis-N_last)/Glob_numprocs
p_Nmax=N_each+1
IF(N_last==0)THEN
  p_Nmax=p_Nmax-1
ENDIF
    
ALLOCATE(p_bk_cal(p_Nmax,Glob_numprocs))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
p_bk_cal=0
index1=1 !process sign
index2=0 !basis sign
DO i=1,Nbasis
  index2=index2+1
  p_bk_cal(index2,index1)=i
  IF((index2==(N_each+1)).AND.(index1<=N_last))THEN
    index1=index1+1
    index2=0
  ELSEIF((index2==N_each).AND.(index1>N_last))THEN
    index1=index1+1
    index2=0
  ENDIF     
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)
  
DEALLOCATE(p_Ei)
DEALLOCATE(p_Ci)
DEALLOCATE(p_hk)
DEALLOCATE(p_Nproc)
DEALLOCATE(p_displs)
DEALLOCATE(p_bk_cal)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE  pvars_eigen

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!full diagonalization  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE full_diag_MPIfun(Nbasis,Ei,Ci,ERR)
!===================================================
!this subroutine perform full diagonalization 
!of eigenvalue problem: HC=ESC
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis
REAL(dp),INTENT(OUT)::Ei,Ci(Nbasis)
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k
INTEGER::myid,numprocs,Ntemp
INTEGER::E_target_num
INTEGER::MT,MB

INTEGER::ICTXT,Npx,Npy
INTEGER::NPCOL,NPROW,MYCOL,MYROW
INTEGER::MMYCOL,MMYROW
INTEGER::IXROW,IXLOC,IYROW,IYLOC  

INTEGER::LOCR,LOCC,LOCRC
INTEGER::LOCR0,LOCRC0
INTEGER,ALLOCATABLE,DIMENSION(:,:)::LOCDD 

REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Hmat,Smat,Cmat
REAL(dp),ALLOCATABLE,DIMENSION(:)::Emat
REAL(dp),ALLOCATABLE,DIMENSION(:)::Ctemp
INTEGER DESC_H(9),DESC_S(9),DESC_C(9)

REAL(dp)::ABSTOL
INTEGER::E_num,C_num

INTEGER::LWORK,LIWORK,INFO
INTEGER,ALLOCATABLE,DIMENSION(:)::IWORK,IFAIL,ICLUSTR
REAL(dp),ALLOCATABLE,DIMENSION(:)::WORK,GAP 

INTEGER,EXTERNAL::NUMROC
REAL(dp),EXTERNAL::PDLAMCH

ERR=0
E_target_num=1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!process block division
!======================

CALL BLACS_PINFO(myid,numprocs)
k=INT(SQRT(REAL(numprocs)))!	  

!rows of the process grid: Npx
!columns of the process grid: Npy
loop1:DO i=k,1,-1
loop2: DO j=k,numprocs
    IF((i*j).EQ.numprocs)THEN
      Npx=i
	  Npy=j	  
      EXIT loop1
    ENDIF
  ENDDO loop2
ENDDO loop1

!initialize block
!=======================
!Gets values that BLACS use for internal defaults.
!1: integer handle indicating the system context
!2: indicate what BLACS internal should be return
!3: value of BLACS internal
!=======================
CALL BLACS_GET(-1,0,ICTXT)  

!=======================
!Assigns available processes into BLACS process grid.
!input:
!1: integer handle indicating the system context
!output:
!1: integer handle to the created BLACS context
!2: indicates how to map processes to BLACS grid
!3: indicates how many process rows the process grid should contain.
!4: indicates how many process columns the process grid should contain.
!=======================
CALL BLACS_GRIDINIT(ICTXT,'Row-major',Npx, Npy)

!=======================
!Returns information on the current grid.
!input:
!1: integer handle that indicates the context
!output:
!2: number of process rows in the current process grid.
!3: number of process columns in the current process grid.
!4: row coordinate of the calling process in the process grid
!5: column coordinate of the calling process in the process grid.
!=======================
CALL BLACS_GRIDINFO(ICTXT,NPROW,NPCOL,MYROW,MYCOL)

ALLOCATE(ICLUSTR(2*Npx*Npy),GAP(Npx*Npy),LOCDD(2,numprocs-1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!allocate variables used in PDSYGVX
!======================

MT=Nbasis
LWORK=MT*(10+MT)
LIWORK=6*MT
ALLOCATE(IFAIL(MT),IWORK(LIWORK),WORK(LWORK))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!block size MB
IF((MT/Npy.GE.4).AND.(MT/Npy.LT.8))THEN
  MB=4
ELSEIF((MT/Npy.GE.8).AND.(MT/Npy.LT.16))THEN
  MB=8
ELSEIF((MT/Npy.GE.16).AND.(MT/Npy.LT.64))THEN
  MB=16
ELSEIF(MT/Npy.GE.64)THEN
  MB=64
ENDIF	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==========================
!Computes the number of rows or columns of a distributed matrix owned by the process.
!1: The number of rows/columns in distributed matrix
!2: Block size, size of the blocks the distributed matrix is split into.
!3: The coordinate of the process whose local array row or column is to be determined
!4: The coordinate of the process that possesses the first row or
!   column of the distributed matrix
!5: The total number processes over which the matrix is distributed.
!==========================
LOCR=NUMROC(MT,MB,MYROW,0,NPROW)
LOCC=NUMROC(MT,MB,MYCOL,0,NPCOL)	  
LOCRC=LOCR*LOCC

ALLOCATE(Emat(MT))

!local matrix of very process
ALLOCATE(Hmat(LOCR,LOCC),Smat(LOCR,LOCC),Cmat(LOCR,LOCC))	

!!!!!!!!!!!!!!!!!!!!!!!!!!
!store the block size of every process in process 0
!LOCDD matrix

IF(myid.NE.0)THEN
    
  CALL MPI_SEND(LOCR,1,MPI_INTEGER,0,21,MPI_COMM_WORLD,Glob_MPIerr) 
  CALL MPI_SEND(LOCRC,1,MPI_INTEGER,0,22,MPI_COMM_WORLD,Glob_MPIerr)
  
ELSEIF(myid.EQ.0)THEN
  Ntemp=0
  DO i=1,numprocs-1
    CALL MPI_RECV(j,1,MPI_INTEGER,i,21,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
    CALL MPI_RECV(k,1,MPI_INTEGER,i,22,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
    LOCDD(1,i)=j
    LOCDD(2,i)=k
    IF(Ntemp<k)THEN
      Ntemp=k
    ENDIF
  ENDDO

ALLOCATE(Ctemp(Ntemp))

ENDIF	    

!============================
!initializes the array descriptor for distributed matrix
!1: the array descriptor of a distributed matrix to be set
!2: the number of rows in the distributed matrix
!3: the number of columns in the distributed matrix
!4: the blocking factor used to distribute the rows of the matrix
!5: the blocking factor used to distribute the columns of the matrix
!6: the process row over which the first row of the matrix is distributed. 0 <= IRSRC < NPROW.
!7: the process column over which the first column of the matrix is distributed. 0 <= ICSRC < NPCOL
!8: integer handle that indicates the context
!9: the leading dimension of the local array storing the local
!   blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)). LOCr() denotes
!   the number of rows of a global dense matrix that the process in a grid
!   receives after data distributing.
!10: return information 
!============================
CALL DESCINIT(DESC_H,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)
CALL DESCINIT(DESC_S,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)
CALL DESCINIT(DESC_C,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!distribute H S matrix elements into every process block
DO j=1,MT	  
DO i=1,j
  
IXROW=MOD((i-1)/MB,NPROW)
IXLOC=((i-1)/(MB*NPROW))*MB+MOD(i-1,MB)+1
IYROW=MOD((j-1)/MB,NPCOL)
IYLOC=((j-1)/(MB*NPCOL))*MB+MOD(j-1,MB)+1
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IF((MYROW.EQ.IXROW).AND.(MYCOL.EQ.IYROW))THEN
  Hmat(IXLOC,IYLOC)=Glob_Hkl(i,j)
  Smat(IXLOC,IYLOC)=Glob_Skl(i,j)
ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ENDDO
ENDDO	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ABSTOL=2.d0*PDLAMCH(ICTXT,'U')
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL PDSYGVX(&
& 1,'V','I','U',MT,&    !case selection
& Hmat,1,1,DESC_H,&     !local pieces of H(local)
& Smat,1,1,DESC_S,&     !local pieces of S(local)
& 0.0,1.0,Glob_energy_level,Glob_energy_level,& !RANGE='I',energy level to be optimized
& ABSTOL,E_num,C_num,&  !precision and number of eigenvalues and eigenvectors found   
& Emat,0.0001,&         !eigen value(all)
& Cmat,1,1,DESC_C,&     !eigen vector(local)
& WORK,LWORK,IWORK,LIWORK,IFAIL,ICLUSTR,GAP,INFO) !other paras

IF(myid==0)WRITE(*,*)
IF(INFO/=0)THEN
  IF(myid==0)WRITE(*,*)INFO,Emat(E_target_num)
  ERR=1
  Ei=1.d5
  Ci(:)=0.d0
  RETURN
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!form target eigenvalue Ei and eigenvector Ci from process blocks
!========================	 

IF(myid==0)THEN
    
!target Energy Ei  
Ei=Emat(E_target_num)
CALL MPI_BCAST(Ei,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eigen vector from 0 prcess bolck

DO i=1,MT
        
IXROW=MOD((i-1)/MB,NPROW)
IXLOC=((i-1)/(MB*NPROW))*MB+MOD(i-1,MB)+1
IYROW=MOD((E_target_num-1)/MB,NPCOL)
IYLOC=((E_target_num-1)/(MB*NPCOL))*MB+MOD(E_target_num-1,MB)+1
      
IF((MYROW.EQ.IXROW).AND.(MYCOL.EQ.IYROW))THEN
  Ci(i)=Cmat(IXLOC,IYLOC)
ENDIF

ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eigen vector from other prcess bolck

DO i=1,numprocs-1

LOCR0=LOCDD(1,i)
LOCRC0=LOCDD(2,i)

CALL MPI_RECV(Ctemp,LOCRC0,MPI_DOUBLE_PRECISION,i,11,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
CALL MPI_RECV(MMYROW,1,MPI_INTEGER,i,12,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
CALL MPI_RECV(MMYCOL,1,MPI_INTEGER,i,13,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)	  
      
 DO j=1,MT
   IXROW=MOD((j-1)/MB,NPROW)
   IXLOC=((j-1)/(MB*NPROW))*MB+MOD(j-1,MB)+1
   IYROW=MOD((E_target_num-1)/MB,NPCOL)
   IYLOC=((E_target_num-1)/(MB*NPCOL))*MB+MOD(E_target_num-1,MB)+1

   IF((MMYROW.EQ.IXROW).AND.(MMYCOL.EQ.IYROW))THEN
    Ci(j)=Ctemp(IXLOC+(IYLOC-1)*LOCR0)
   ENDIF

  ENDDO

ENDDO

!bcast Ci to other process
CALL MPI_BCAST(Ci,Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
ELSE
    
!target energy
CALL MPI_BCAST(Ei,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!send local eigenvector to process 0
CALL MPI_SEND(Cmat,LOCRC,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_SEND(MYROW,1,MPI_INTEGER,0,12,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_SEND(MYCOL,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,Glob_MPIerr)

!target eigenvector
CALL MPI_BCAST(Ci,Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)
    
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF((dabs(Ei)>1.d4).OR.(isnan(Ei)))THEN
  ERR=1
  Ei=1.d5
  Ci(:)=0.d0
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BLACS_GRIDEXIT(ICTXT)	
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE full_diag_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!diagonalization when only one basis varies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE one_diag_MPIfun(Nbasis,nb,Ei,ERR)
!===================================================
!this subroutine perform diagnalization when only
!one basis varies diag_MPIfun should be called to 
!derive p_Ei and p_Ci before calling this subroutine
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,nb
REAL(dp),INTENT(OUT)::Ei
INTEGER,INTENT(OUT)::ERR
!@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER,PARAMETER::bracket_num1=10
INTEGER,PARAMETER::bracket_num2=10
!@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER::i,j,k,ib,myid,ERROR,Nproc,bk_cal(p_Nmax)
REAL(dp)::Sk(Nbasis),Hk(Nbasis)
REAL(dp)::temp,norm,root,x1,x2,dx,Echeck

REAL(dp)::Sk_send(p_Nmax),Hk_send(p_Nmax),hkk_send(p_Nmax),hkk_nb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_COMM_RANK(MPI_COMM_WORLD,myid,Glob_MPIerr)
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

ERR=0
IF(Glob_energy_level>=nb)THEN
 Echeck=p_Ei(Glob_energy_level+1)
ELSE
  Echeck=p_Ei(Glob_energy_level)
ENDIF
IF((dabs(Echeck)>1.d4).OR.(isnan(Echeck)))THEN
  ERR=1
 Ei=1.d5 
 RETURN
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

Nproc=p_Nproc(myid+1)
bk_cal(:)=p_bk_cal(:,myid+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

norm=0.0_dp
hkk_nb=0.0_dp

IF(Nproc/=0)THEN
DO ib=1,Nproc 
k=bk_cal(ib)
IF(k/=nb)THEN

!Sk=<psi_k|phi_nb>=sum_{i/=nb}C_i^k<phi_i|psi_nb>
temp=0.0_dp
DO i=1,Nbasis   
  temp=temp+p_Ci(i,k)*Glob_Skl(i,nb)
ENDDO
Sk_send(ib)=temp

!Hk=<psi_k|H|phi_nb>=sum_{i/=nb}C_i^k<phi_i|H|psi_nb>
temp=0.0_dp
DO i=1,Nbasis
  temp=temp+p_Ci(i,k)*Glob_Hkl(i,nb)
ENDDO
Hk_send(ib)=temp

!norm  
norm=norm-Sk_send(ib)*Sk_send(ib)
  
!hkk(k/=nb)=<psi_k|psi_nb>
hkk_send(ib)=(Hk_send(ib)-p_Ei(k)*Sk_send(ib))

!hkk(k=nb)=<psi_nb|psi_nb>
hkk_nb=hkk_nb-2.d0*Hk_send(ib)*Sk_send(ib)+p_Ei(k)*Sk_send(ib)*Sk_send(ib)

ELSE
    
  Sk_send(ib)=0.0_dp
  Hk_send(ib)=0.0_dp
  hkk_send(ib)=0.0_dp
  
ENDIF
ENDDO
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_GATHERV(Sk_send,Nproc,MPI_DOUBLE_PRECISION,Sk,p_Nproc,&
     &p_displs,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_GATHERV(Hk_send,Nproc,MPI_DOUBLE_PRECISION,Hk,p_Nproc,&
     &p_displs,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

CALL MPI_GATHERV(hkk_send,Nproc,MPI_DOUBLE_PRECISION,p_hk,p_Nproc,&
     &p_displs,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

temp=norm
CALL MPI_REDUCE(temp,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)

temp=hkk_nb
CALL MPI_REDUCE(temp,hkk_nb,1,MPI_DOUBLE_PRECISION,MPI_SUM,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)


IF(myid==Glob_root)THEN
    
norm=norm+Glob_Skl(nb,nb)

DO i=1,Nbasis
  p_hk(i)=p_hk(i)/dsqrt(norm)
ENDDO
p_hk(nb)=(hkk_nb+Glob_Hkl(nb,nb))/norm

ENDIF
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=======================
!search for the roots of characteristic function
!by Glob_root process
!=======================
IF(myid==Glob_root)THEN  

!initialize
IF(Glob_energy_level==1)THEN
  x2=p_Ei(Glob_energy_level)
  IF(nb==1)THEN
    x2=p_Ei(Glob_energy_level+1)
  ENDIF
  dx=0.005d0
ELSE
  IF(Glob_energy_level<nb)THEN
    x1=p_Ei(Glob_energy_level-1)
    x2=p_Ei(Glob_energy_level)
  ELSEIF((Glob_energy_level-1)>=nb)THEN
    x1=p_Ei(Glob_energy_level)
    x2=p_Ei(Glob_energy_level+1)    
  ELSE
    x1=p_Ei(Glob_energy_level-1)
    x2=p_Ei(Glob_energy_level+1)   
  ENDIF
  dx=(x2-x1)/dble(bracket_num2)
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(Glob_energy_level==1)THEN  

bracket_loop1:DO ib=1,bracket_num1
!@@@@@@@@@@@@@@@@@@@@
x1=x2-dble(2.d0**ib)*dx
!@@@@@@@@@@@@@@@@@@@@
CALL root_fun(nb,Nbasis,x1,x2,root,ERROR)
IF(ERROR==0)THEN
  Ei=root
  EXIT bracket_loop1
ELSE
  Ei=1.d5 
  CYCLE bracket_loop1
ENDIF
ENDDO bracket_loop1

ELSE

bracket_loop2:DO ib=1,bracket_num2
!@@@@@@@@@@@@@@@@@
x1=x2-dble(ib)*dx
!@@@@@@@@@@@@@@@@@
CALL root_fun(nb,Nbasis,x1,x2,root,ERROR)
IF(ERROR==0)THEN
  Ei=root
  EXIT bracket_loop2
ELSE
  Ei=1.d5 
  CYCLE bracket_loop2
ENDIF
ENDDO bracket_loop2

ENDIF

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_BCAST(Ei,1,MPI_DOUBLE_PRECISION,Glob_root,MPI_COMM_WORLD,Glob_MPIerr)
ERR=0
IF((dabs(Ei)>1.d4).OR.(isnan(Ei)))THEN
  ERR=1
  Ei=1.d5 
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE one_diag_MPIfun  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!diagonalization initialization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE one_diag_init_MPIfun(Nbasis,nb)
!===================================================
!this subroutine perform diagonalization of Nbasis except the nb_th
!to provide eigen value p_Ei and eigen vector p_Ci for one_diag_MPIfun
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::Nbasis,nb

INTEGER::i,j,k,ip,myid,numprocs,Ntemp
INTEGER::MT,MB

INTEGER::ICTXT,Npx,Npy
INTEGER::NPCOL,NPROW,MYCOL,MYROW
INTEGER::MMYCOL,MMYROW
INTEGER::IXROW,IXLOC,IYROW,IYLOC  

INTEGER::LOCR,LOCC,LOCRC
INTEGER::LOCR0,LOCRC0
INTEGER,ALLOCATABLE,DIMENSION(:,:)::LOCDD 

REAL(dp),ALLOCATABLE,DIMENSION(:,:)::Hmat,Smat,Cmat
REAL(dp),ALLOCATABLE,DIMENSION(:)::Emat
REAL(dp),ALLOCATABLE,DIMENSION(:)::Ctemp
INTEGER DESC_H(9),DESC_S(9),DESC_C(9)

REAL(dp)::ABSTOL
INTEGER::E_num,C_num

INTEGER::LWORK,LIWORK,INFO
INTEGER,ALLOCATABLE,DIMENSION(:)::IWORK,IFAIL,ICLUSTR
REAL(dp),ALLOCATABLE,DIMENSION(:)::WORK,GAP 

INTEGER,EXTERNAL::NUMROC
REAL(dp),EXTERNAL::PDLAMCH

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!process block division
!======================

CALL BLACS_PINFO(myid,numprocs)
k=INT(SQRT(REAL(numprocs)))	  

!rows and columns of the process grid: Npx,Npy
loop1:DO i=k,1,-1
loop2: DO j=k,numprocs
    IF((i*j).EQ.numprocs)THEN
      Npx=i
	    Npy=j	  
      EXIT loop1
    ENDIF
  ENDDO loop2
ENDDO loop1

!initialize block
!=======================
!Gets values that BLACS use for internal defaults.
!1: integer handle indicating the system context
!2: indicate what BLACS internal should be return
!3: value of BLACS internal
!=======================
CALL BLACS_GET(-1,0,ICTXT)  

!=======================
!Assigns available processes into BLACS process grid.
!input:
!1: integer handle indicating the system context
!output:
!1: integer handle to the created BLACS context
!2: indicates how to map processes to BLACS grid
!3: indicates how many process rows the process grid should contain.
!4: indicates how many process columns the process grid should contain.
!=======================
CALL BLACS_GRIDINIT(ICTXT,'Row-major',Npx, Npy)

!=======================
!Returns information on the current grid.
!input:
!1: integer handle that indicates the context
!output:
!2: number of process rows in the current process grid.
!3: number of process columns in the current process grid.
!4: row coordinate of the calling process in the process grid
!5: column coordinate of the calling process in the process grid.
!=======================
CALL BLACS_GRIDINFO(ICTXT,NPROW,NPCOL,MYROW,MYCOL)

ALLOCATE(ICLUSTR(2*Npx*Npy),GAP(Npx*Npy),LOCDD(2,numprocs-1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!allocate variables used in PDSYGVX
!======================

MT=Nbasis-1
LWORK=MT*(10+MT)
LIWORK=6*MT
ALLOCATE(IFAIL(MT),IWORK(LIWORK),WORK(LWORK))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!block size MB
IF((MT/Npy.GE.4).AND.(MT/Npy.LT.8))THEN
  MB=4
ELSEIF((MT/Npy.GE.8).AND.(MT/Npy.LT.16))THEN
  MB=8
ELSEIF((MT/Npy.GE.16).AND.(MT/Npy.LT.64))THEN
  MB=16
ELSEIF(MT/Npy.GE.64)THEN
  MB=64
ENDIF	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!==========================
!Computes the number of rows or columns of a distributed matrix owned by the process.
!1: The number of rows/columns in distributed matrix
!2: Block size, size of the blocks the distributed matrix is split into.
!3: The coordinate of the process whose local array row or column is to be determined
!4: The coordinate of the process that possesses the first row or
!   column of the distributed matrix
!5: The total number processes over which the matrix is distributed.
!==========================
LOCR=NUMROC(MT,MB,MYROW,0,NPROW)
LOCC=NUMROC(MT,MB,MYCOL,0,NPCOL)	  
LOCRC=LOCR*LOCC

ALLOCATE(Emat(MT))

!local matrix of very process
ALLOCATE(Hmat(LOCR,LOCC),Smat(LOCR,LOCC),Cmat(LOCR,LOCC))	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!store the block size of every process in process 0
!LOCDD matrix

IF(myid.NE.0)THEN
    
  CALL MPI_SEND(LOCR,1,MPI_INTEGER,0,21,MPI_COMM_WORLD,Glob_MPIerr) 
  CALL MPI_SEND(LOCRC,1,MPI_INTEGER,0,22,MPI_COMM_WORLD,Glob_MPIerr)
  
ELSEIF(myid.EQ.0)THEN
  Ntemp=0
  DO i=1,numprocs-1
    CALL MPI_RECV(j,1,MPI_INTEGER,i,21,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
    CALL MPI_RECV(k,1,MPI_INTEGER,i,22,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
    LOCDD(1,i)=j
    LOCDD(2,i)=k
    IF(Ntemp<k)THEN
      Ntemp=k
    ENDIF
  ENDDO

ALLOCATE(Ctemp(Ntemp))

ENDIF	    

!============================
!initializes the array descriptor for distributed matrix
!1: the array descriptor of a distributed matrix to be set
!2: the number of rows in the distributed matrix
!3: the number of columns in the distributed matrix
!4: the blocking factor used to distribute the rows of the matrix
!5: the blocking factor used to distribute the columns of the matrix
!6: the process row over which the first row of the matrix is distributed. 0 <= IRSRC < NPROW.
!7: the process column over which the first column of the matrix is distributed. 0 <= ICSRC < NPCOL
!8: integer handle that indicates the context
!9: the leading dimension of the local array storing the local
!   blocks of the distributed matrix. LLD >= MAX(1,LOCr(M)). LOCr() denotes
!   the number of rows of a global dense matrix that the process in a grid
!   receives after data distributing.
!10: return information 
!============================
CALL DESCINIT(DESC_H,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)
CALL DESCINIT(DESC_S,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)
CALL DESCINIT(DESC_C,MT,MT,MB,MB,0,0,ICTXT,LOCR,INFO)	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!distribute H S matrix elements into every process block
DO j=1,MT	  
DO i=1,j
    
IXROW=MOD((i-1)/MB,NPROW)
IXLOC=((i-1)/(MB*NPROW))*MB+MOD(i-1,MB)+1
IYROW=MOD((j-1)/MB,NPCOL)
IYLOC=((j-1)/(MB*NPCOL))*MB+MOD(j-1,MB)+1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IF((MYROW.EQ.IXROW).AND.(MYCOL.EQ.IYROW))THEN
  IF(j<nb)THEN
    Hmat(IXLOC,IYLOC)=Glob_Hkl(i,j)
    Smat(IXLOC,IYLOC)=Glob_Skl(i,j)
  ELSE
    IF(i<nb)THEN
      Hmat(IXLOC,IYLOC)=Glob_Hkl(i,j+1)
      Smat(IXLOC,IYLOC)=Glob_Skl(i,j+1)
    ELSE
      Hmat(IXLOC,IYLOC)=Glob_Hkl(i+1,j+1)
      Smat(IXLOC,IYLOC)=Glob_Skl(i+1,j+1)
    ENDIF
  ENDIF
ENDIF	
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
ENDDO
ENDDO	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ABSTOL=2.d0*PDLAMCH(ICTXT,'U')

CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
CALL PDSYGVX(&
& 1,'V','A','U',MT,&    !case selection
& Hmat,1,1,DESC_H,&     !local pieces of H(local)
& Smat,1,1,DESC_S,&     !local pieces of S(local)
& 0.0,1.0,1,1,&         !RANGE='A',they are not referenced
& ABSTOL,E_num,C_num,&  !precision and number of eigenvalues and eigenvectors found   
& Emat,0.0001,&         !eigen value(all)
& Cmat,1,1,DESC_C,&     !eigen vector(local)
& WORK,LWORK,IWORK,LIWORK,IFAIL,ICLUSTR,GAP,INFO) !other paras

IF(myid==0)WRITE(*,*)
IF(INFO/=0)THEN
WRITE(*,*)INFO,Emat(1)
PAUSE'one_diag_MPIfun error'
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!========================
!form target eigenvalue Ei and eigenvector Ci from process blocks
!========================	 

IF(myid==0)THEN
    
!target Energy Ei 
IF(nb>1)THEN
p_Ei(1:nb-1)=Emat(1:nb-1)
ENDIF
p_Ei(nb)=0.0_dp
IF(nb<Nbasis)THEN
  p_Ei(nb+1:Nbasis)=Emat(nb:Nbasis-1)
ENDIF

!bcast E to other process
CALL MPI_BCAST(p_Ei,Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eigen vector from 0 prcess bolck
DO j=1,MT
DO i=1,MT
        
IXROW=MOD((i-1)/MB,NPROW)
IXLOC=((i-1)/(MB*NPROW))*MB+MOD(i-1,MB)+1
IYROW=MOD((j-1)/MB,NPCOL)
IYLOC=((j-1)/(MB*NPCOL))*MB+MOD(j-1,MB)+1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   
IF((MYROW.EQ.IXROW).AND.(MYCOL.EQ.IYROW))THEN
  IF(j<nb)THEN
    IF(i<nb)THEN
      p_Ci(i,j)=Cmat(IXLOC,IYLOC)
    ELSE
      p_Ci(i+1,j)=Cmat(IXLOC,IYLOC)
    ENDIF
  ELSE
    IF(i<nb)THEN
      p_Ci(i,j+1)=Cmat(IXLOC,IYLOC)
    ELSE
      p_Ci(i+1,j+1)=Cmat(IXLOC,IYLOC)
    ENDIF
  ENDIF
ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!eigen vector from other prcess bolck

DO ip=1,numprocs-1

LOCR0=LOCDD(1,ip)
LOCRC0=LOCDD(2,ip)

CALL MPI_RECV(Ctemp,LOCRC0,MPI_DOUBLE_PRECISION,ip,11,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
CALL MPI_RECV(MMYROW,1,MPI_INTEGER,ip,12,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)
CALL MPI_RECV(MMYCOL,1,MPI_INTEGER,ip,13,MPI_COMM_WORLD,Glob_status,Glob_MPIerr)	  
      
DO j=1,MT
DO i=1,MT
       
IXROW=MOD((i-1)/MB,NPROW)
IXLOC=((i-1)/(MB*NPROW))*MB+MOD(i-1,MB)+1
IYROW=MOD((j-1)/MB,NPCOL)
IYLOC=((j-1)/(MB*NPCOL))*MB+MOD(j-1,MB)+1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
IF((MMYROW.EQ.IXROW).AND.(MMYCOL.EQ.IYROW))THEN
  IF(j<nb)THEN
    IF(i<nb)THEN
      p_Ci(i,j)=Ctemp(IXLOC+(IYLOC-1)*LOCR0)
    ELSE
      p_Ci(i+1,j)=Ctemp(IXLOC+(IYLOC-1)*LOCR0)
    ENDIF
  ELSE
    IF(i<nb)THEN
      p_Ci(i,j+1)=Ctemp(IXLOC+(IYLOC-1)*LOCR0)
    ELSE
      p_Ci(i+1,j+1)=Ctemp(IXLOC+(IYLOC-1)*LOCR0)
    ENDIF  
  ENDIF    
ENDIF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

ENDDO
ENDDO

ENDDO

p_Ci(1:Nbasis,nb)=0.0_dp
p_Ci(nb,1:Nbasis)=0.0_dp

!bcast Ei to other process
CALL MPI_BCAST(p_Ci,Nbasis*Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
ELSE  !other procs
    
!target energy
CALL MPI_BCAST(p_Ei,Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)

!send local eigenvector to process 0
CALL MPI_SEND(Cmat,LOCRC,MPI_DOUBLE_PRECISION,0,11,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_SEND(MYROW,1,MPI_INTEGER,0,12,MPI_COMM_WORLD,Glob_MPIerr)
CALL MPI_SEND(MYCOL,1,MPI_INTEGER,0,13,MPI_COMM_WORLD,Glob_MPIerr)

!target eigenvector
CALL MPI_BCAST(p_Ci,Nbasis*Nbasis,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,Glob_MPIerr)
    
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL BLACS_GRIDEXIT(ICTXT)	
CALL MPI_BARRIER(MPI_COMM_WORLD,Glob_MPIerr)
RETURN
END SUBROUTINE one_diag_init_MPIfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!characteristic function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION D_fun(nb,Nbasis,E)
!===================================================
!characteristic function of E
!===================================================
INTEGER,INTENT(IN)::nb,Nbasis
REAL(dp),INTENT(IN)::E
REAL(dp)::D_fun

INTEGER::k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

D_fun=p_hk(nb)-E
DO k=1,Nbasis
IF(k/=nb)THEN
  D_fun=D_fun-p_hk(k)*p_hk(k)/(p_Ei(k)-E)
ENDIF
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END FUNCTION D_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!search the root of characteristic function 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE root_fun(nb,Nbasis,x1,x2,root,ERR)
!===================================================
!search for the zeros for function D_fun
!===================================================
INTEGER,INTENT(IN)::nb,Nbasis
REAL(dp),INTENT(IN)::x1,x2
REAL(dp),INTENT(OUT)::root
INTEGER,INTENT(OUT)::ERR
!@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER,PARAMETER::ITmax=300 !the max search time 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@
INTEGER::i,j,k,it,case1,case2
REAL(dp)::z1,z2,z3,z12,z13,z23,w,zz
REAL(dp)::g1,g2,g3,f2,temp

ERR=0

z1=x1
z2=0.5d0*(x1+x2)
z3=x2

f2=D_fun(nb,Nbasis,z2)
IF(dabs(f2)<=EPS)THEN
  root=z2
  RETURN
ENDIF

g1=1.0_dp/D_fun(nb,Nbasis,z1)
g2=1.0_dp/f2
g3=1.0_dp/D_fun(nb,Nbasis,z3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!======================
!bisection loop
!=====================
loop:DO it=1,ITmax

case1=1
case2=1

IF((g1<0.0_dp.AND.g2<0.0_dp).OR.(g1>0.0_dp.AND.g2>0.0_dp))THEN
    case1=0
ENDIF
IF((g2<0.0_dp.AND.g3<0.0_dp).OR.(g2>0.0_dp.AND.g3>0.0_dp))THEN
    case2=0
ENDIF

IF(case1+case2==0)THEN
  root=z1
  IF(it==1)THEN
    ERR=1
  ENDIF
  RETURN
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

z12=z1-z2
z23=z2-z3
z13=z1-z3
w=g1*z23-g2*z13+g3*z12
IF(w/=0.0_dp)THEN
  zz=(g1*z1*z23-g2*z2*z13+g3*z3*z12)/w
ELSE
  zz=1.0d16
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!root is between z1 and z2
!==========================
IF(case1==1)THEN
    
z3=z2
g3=g2
IF(zz<z1.OR.zz>z2)THEN
  z2=0.5d0*(z1+z2) 
ELSE
  z2=zz
ENDIF

f2=D_fun(nb,Nbasis,z2)
IF(dabs(f2)<=EPS)THEN
  root=z2
  RETURN
ENDIF
g2=1.d0/f2
    
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!root is between z2 and z3
!==========================    
IF(case2==1)THEN    
   
z1=z2
g1=g2
IF(zz<z2.OR.zz>z3)THEN
  z2=0.5d0*(z2+z3) 
ELSE
  z2=zz
ENDIF

f2=D_fun(nb,Nbasis,z2)
IF(dabs(f2)<=EPS)THEN
  root=z2
  RETURN
ENDIF
g2=1.d0/f2

ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(dabs(zz-z3)<=EPS)THEN
  root=zz
  RETURN
ENDIF

ENDDO loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ERR=1
root=1.d5

RETURN
END SUBROUTINE root_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE  eigen
