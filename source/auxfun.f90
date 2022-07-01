MODULE auxfun
!==========================================================
!this module contains auxiliary functions and subroutines for matrix manipulation
!they are for general funcitons,no Global variables should be used here
!except the ones regarding the accuracy like dp
!==========================================================
USE globvars

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!determinant and inverse of real symmetric matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
SUBROUTINE det_inv_fun(N,A,det_A,inv_A)
!===================================================
!calculate the determinant and inverse of symmetric and positive
!definite Hermitian matrix A based on Cholesky decomposition
!===================================================
!input:
!  N: size of matrix A
!  A: the symmetric positive definite matrix
!output:
!  det_A: determinant of A
!  inv_A: inverse of A
!=================================================== 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dp),INTENT(IN)::A(N,N)
REAL(dp),INTENT(OUT)::det_A
REAL(dp),INTENT(OUT)::inv_A(N,N)

INTEGER::i,j,k
REAL(dp)::temp
REAL(dp)::L(N,N),inv_L(N,N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(N==1)THEN

det_A=A(1,1)
inv_A=1.0_dp/A(1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE

L=0.0_dp

!first column of L
L(1,1)=dsqrt(A(1,1))
DO i=2,N
  L(i,1)=A(i,1)/L(1,1)
ENDDO
det_A=A(1,1)

!second to N column
DO j=2,N   !j_th column
  DO i=j,N !i_th row
    temp=A(i,j)
    DO k=1,j-1
      temp=temp-L(i,k)*L(j,k)
    ENDDO
    IF(i==j)THEN
      L(i,j)=dsqrt(temp)
      det_A=det_A*temp    !det_A=prod_i(L_ii*L_ii)
    ELSE
      L(i,j)=temp/L(j,j)
    ENDIF
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!inv_L
inv_L=0.0_dp
inv_L(1,1)=1.0_dp/L(1,1)
DO i=2,N               !i_th row
  inv_L(i,i)=1.0_dp/L(i,i)!diagonal
  DO j=1,i-1          !j_th column
    temp=0.0_dp
    DO k=j,i-1
      temp=temp-L(i,k)*inv_L(k,j)
    ENDDO
    inv_L(i,j)=temp/L(i,i)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!inv_A=transpose(inv_L)*inv_L
inv_A=0.0_dp
DO i=1,N
  DO j=i,N
    temp=0.0_dp
    DO k=1,N
      temp=temp+inv_L(k,i)*inv_L(k,j)
    ENDDO
    inv_A(i,j)=temp
    inv_A(j,i)=temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF

RETURN
END SUBROUTINE det_inv_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ckolesky decomposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Ckolesky_decompose(N,vechA,vechL)
!====================================================
!Cholesky decomposition of symmetric matrix A
!====================================================
! input:
!  N: the size of matrix A
!  vechA: symmetric matrix to be decomposed
!output:
!  vechL: the lower trangular matrix(A=L*L^T)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dp),INTENT(IN)::vechA(N*(N+1)/2)
REAL(dp),INTENT(OUT)::vechL(N*(N+1)/2)

INTEGER::i,j,k,index
REAL(dp)::temp,A(N,N),L(N,N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(N==1)THEN

vechL(1)=dsqrt(vechA(1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE

!vechA->A
index=0
DO i=1,N
  DO j=i,N
    index=index+1
    A(i,j)=vechA(index)
    A(j,i)=vechA(index)
  ENDDO
ENDDO

!first column of L
L(1,1)=dsqrt(A(1,1))
DO i=2,N
  L(i,1)=A(i,1)/L(1,1)
ENDDO

!second to last column of L
DO j=2,N   !j_th column
  DO i=j,N !i_th row
    temp=A(i,j)
    DO k=1,j-1
      temp=temp-L(i,k)*L(j,k)
    ENDDO
    IF(i==j)THEN
      L(i,j)=dsqrt(temp)
    ELSE
      L(i,j)=temp/L(j,j)
    ENDIF
  ENDDO
ENDDO

!L->vechL
index=0
DO i=1,N
  DO j=i,N
    index=index+1
    vechL(index)=L(j,i)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF

RETURN
END SUBROUTINE Ckolesky_decompose

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!inverse of lower trangular matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE inv_L_lower(N,L,inv_L)
!===================================================
!calculate the inverse of lower trangular matrix L
!===================================================
!input:
!  N: matrix size 
!  L: lower trangular matrix L(N,N)
!output:
!inv_L: the inverse of Ll
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dp),INTENT(IN)::L(N,N)
REAL(dp),INTENT(OUT)::inv_L(N,N)

INTEGER::i,j,k
REAL(dp)::temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(N==1)THEN

inv_L(1,1)=1.0_dp/L(1,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE

inv_L(:,:)=0.0_dp
inv_L(1,1)=1.0_dp/L(1,1)
DO i=2,N               !i_th row
  inv_L(i,i)=1.0_dp/L(i,i)!diagonal
  DO j=1,i-1           !j_th column
    temp=0.0_dp
    DO k=j,i-1
      temp=temp-L(i,k)*inv_L(k,j)
    ENDDO
    inv_L(i,j)=temp/L(i,i)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ENDIF

RETURN
END SUBROUTINE inv_L_lower
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!check positive_definite
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PD_check_fun(N,vechA,ERR)
!===================================================
!check wheather symmetric matrix A is positive-definite
!by ckecking wheather the diagonal element is larger than 0
!when performing Colesky decomposition
!===================================================
!input:
!  N: matrix dimension of A
!  A: squre symmetric matrix to be checked
!output: 
!ERR: 
!    0: positive definite; 
!    1: non-positive definite matrix
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL(dp),INTENT(IN)::vechA(N*(N+1)/2)
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,index
REAL(dp)::sum,Atemp(N,N),p(N)

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

index=0
DO i=1,N
  DO j=i,N
    index=index+1
    Atemp(i,j)=vechA(index)
    Atemp(j,i)=vechA(index)
  ENDDO
ENDDO

DO i=1,N
  DO j=i,N
    sum=Atemp(i,j)
    DO k=i-1,1,-1
      sum=sum-Atemp(i,k)*Atemp(j,k)
    ENDDO
    IF(i==j)THEN
      IF(sum<=0.d0)THEN
        ERR=1
        RETURN
      ENDIF
      p(i)=dsqrt(sum)
    ELSE
      Atemp(j,i)=sum/p(i)
    ENDIF
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE PD_check_fun
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!generation of random number
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION randnum_fun(inum)
!===================================================
!generate Long period (>10^18) random number with 
!Bays-Durham shuffle and added safeguards.
!===================================================
!input:
!  inum: negative integer initial number(fixed in a series generation) 
!===================================================
IMPLICIT NONE
INTEGER::inum
REAL(8)::randnum_fun
  
INTEGER,PARAMETER::M1=2147483563,M2=2147483399,M3=M1-1
INTEGER,PARAMETER::A1=40014,A2=40692
INTEGER,PARAMETER::Q1=53668,Q2=52774
INTEGER,PARAMETER::R1=12211,R2=3791
INTEGER,PARAMETER::NTAB=32,NDIV=1+M3/NTAB
REAL(8),PARAMETER::AM=1.d0/M1
REAL(8),PARAMETER::EPS=1.2d-7,RNMX=1.d0-EPS
  
INTEGER::i,j,k
INTEGER,SAVE::num,IY,IV(NTAB) 
  
DATA num/123456789/,IV/NTAB*0/,IY/0/
  
IF(inum<=0)THEN
  inum=max(-inum,1)
  num=inum
  DO i=NTAB+8,1,-1
    j=inum/Q1
    inum=A1*(inum-j*Q1)-j*R1
    IF(inum<0) inum=inum+M1
    IF(i<=NTAB) IV(i)=inum
  ENDDO
  IY=IV(1)  
ENDIF
  
k=inum/Q1
inum=A1*(inum-k*Q1)-k*R1
IF(inum<0)THEN
  inum=inum+M1
ENDIF
  
k=num/Q2
num=A2*(num-k*Q2)-k*R2
IF(num<0)THEN
  num=num+M2
ENDIF
  
i=1+IY/NDIV
IY=IV(i)-num
IV(i)=inum
IF(IY<1)THEN
  IY=IY+M3
ENDIF
randnum_fun=min(AM*IY,RNMX)
  
RETURN
END FUNCTION randnum_fun 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine calculate the system time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE time_cal_fun(time_st)
!===================================================
!calculate and write the time into character time_st 
!must be a character of 14 and the Glob_start time 
!must be allocated
!===================================================
!output:
!  the time wirtten as a character
!===================================================
CHARACTER(14),INTENT(OUT)::time_st

REAL(dp)::time_now
INTEGER::time
INTEGER::temp
INTEGER::hour_now,minite_now,second_now

time_st="????????:??:??"

!time difference
CALL CPU_TIME(time_now)
time=nint(time_now-Glob_time_start)

temp=mod(time,3600)
hour_now=(time-temp)/3600
second_now=mod(temp,60)
minite_now=(temp-second_now)/60

!hour
WRITE(time_st(1:8),fmt="(TL1,I8)")hour_now

!minite
WRITE(time_st(10:11),fmt="(TL1,I2.2)")minite_now

!second
WRITE(time_st(13:14),fmt="(TL1,I2.2)")second_now

RETURN
END SUBROUTINE time_cal_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!factorial function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION fac(k)
!===================================================
!calculate the factorial of k:  k!
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::k
INTEGER::fac

INTEGER::i

fac=1
DO i=1,k
    fac=fac*i
ENDDO

RETURN
END FUNCTION fac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!double-factorial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION dfac(k)
!===================================================
!calculate the double-factorial of k:  k!!
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::k
INTEGER::dfac

INTEGER::i

i=k
dfac=k
DO WHILE(i>2)
i=i-2
dfac=dfac*i
ENDDO

RETURN
END FUNCTION dfac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!delta function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
FUNCTION delta(i,j)
!===================================================
!delta function
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::i,j
REAL(dp)::delta

delta=0.0_dp
IF(i==j)THEN 
delta=1.0_dp
ENDIF

RETURN
END FUNCTION delta   
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE auxfun