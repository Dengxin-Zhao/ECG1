MODULE symmetry
!====================================================
!this module contains subroutines for symmetry operator
!====================================================
USE globvars
USE auxfun

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate the symmetry factor due to spin and permutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE symmetry_fun()
!================================================
!generate the symmetry factor due to spin and permutation parity
!================================================
!output:
!  Glob_symmetry(Glob_Nperm): symmetry factor
!================================================
IMPLICIT NONE
INTEGER::i,j,ip,mi,mj
REAL(dp)::temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO ip=1,Glob_Nperm

Glob_symmetry(ip)=0.0_dp

DO i=1,Glob_Nspin
DO j=1,Glob_Nspin
    
temp=Glob_permut_parity(ip)*Glob_spin_parity(i)*Glob_spin_parity(j)
    
DO mi=1,Glob_Nparticle
  mj=Glob_P(mi,ip)
  temp=temp*delta(Glob_ptype(mi),Glob_ptype(mj))&
  &*delta(Glob_spin_config(mi,i),Glob_spin_config(mj,j))
ENDDO

Glob_symmetry(ip)=Glob_symmetry(ip)+temp

ENDDO
ENDDO

ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE symmetry_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to write symmetry factor into file 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_symmetry_fun()
!===================================================
!This program write symmetry factors into symetry.txt
!===================================================
IMPLICIT NONE
INTEGER::i,ip
CHARACTER(100)::IRAst1,IRAst2
IRAst1="(f10.5,A1,????????I8,A5)"
IRAst2="(????????f10.5)"

OPEN(unit=20,file='symmetry.txt')
WRITE(IRAst1(11:18),fmt="(TL1,I8.8)")Glob_Nparticle
WRITE(IRAst2(2:9),fmt="(TL1,I8.8)")Glob_Np

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO ip=1,Glob_Nperm
  IF(dabs(Glob_symmetry(ip))>0.0_dp)THEN
    WRITE(20,IRAst1)Glob_symmetry(ip),'(',Glob_P(1:Glob_Nparticle,ip),')'
    WRITE(20,*)'=================================='
    DO i=1,Glob_Np
      WRITE(20,IRAst2)Glob_Tp(i,1:Glob_Np,ip)
    ENDDO
     WRITE(20,*)'=================================='
  ENDIF
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CLOSE(20)
RETURN
END SUBROUTINE write_symmetry_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!calculate the permutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE permut_mat_fun()
!==================================================
!generate the permutation operator and its parity
!==================================================
!output:
!            Glob_P: all permutations
!           Glob_Tp: permutations matrix
!Glob_permut_parity: parity of each permutation 
!==================================================
IMPLICIT NONE
INTEGER::ip,i,j,k,l,pk,ipermut(Glob_Nparticle)
REAL(dp)::iparity,temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

loop:DO ip=1,Glob_Nperm
    
CALL permut_fun(ip,Glob_Nparticle,ipermut)
CALL permut_parity_fun(ipermut,Glob_Nparticle,iparity)

Glob_P(:,ip)=ipermut(:)        !ip permutation
Glob_permut_parity(ip)=iparity !ip permutation parity

!permutation matrix
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
	  DO k=1,Glob_Nparticle
      pk=ipermut(k)
	    temp=temp+Glob_U(i,k)*Glob_inv_U(pk,j)
	  ENDDO
    Glob_Tp(i,j,ip)=temp
  ENDDO
ENDDO
  
ENDDO loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE permut_mat_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!generate premutation of N particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE permut_fun(k,N,ipermut)
!==============================================
!generate all permutations of N particles
!(all rearangement of N numbers from 1 to N)
!==============================================
!input:
!  k: the k_th permutation
!  N: particle number
!output:
!  ipermut: column vector containing the k_th permutation(ipermut(i)->j)
!==============================================
IMPLICIT NONE
INTEGER,INTENT(IN)::k,N
INTEGER,INTENT(OUT)::ipermut(N)

INTEGER::i,j,is,io,iv(N+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO i=1,N
  iv(i)=i
ENDDO

io=k-1
DO i=N-1,1,-1

is=io/fac(i)+1
io=MOD(io,fac(i))
ipermut(N-i)=iv(is)
  
DO j=is,i
  iv(j)=iv(j+1)
ENDDO

ENDDO

ipermut(N)=iv(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE permut_fun
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!parity of permutation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  
SUBROUTINE permut_parity_fun(ipermut,N,iparity)
!=================================================
!the parity of permutation ia
!=================================================
!input:
!  ipermut:permutation
!   N:particle number
!output:
!  iparity:parity of ipermut permutation
!=================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ipermut(N)
REAL(dp),INTENT(OUT)::iparity

INTEGER::i,j,M,temp,itemp(N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO i=1,N
  itemp(i)=ipermut(i)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

M=0

DO j=1,N
DO i=j,N
    
IF(itemp(i)==j.AND.i/=j)THEN
  temp=itemp(j)
  itemp(j)=itemp(i)
  itemp(i)=temp
  M=M+1
ENDIF   
  
ENDDO
ENDDO  

iparity=(-1.0_dp)**M 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE permut_parity_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE symmetry