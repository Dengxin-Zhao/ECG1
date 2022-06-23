MODULE symmetry
!====================================================
!this module contains subroutines for symmetry operator
!====================================================
  USE globvars
    
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE symmetry_fun()
!===========================OK
!generate the symmetry factor due to spin and permutation parity
!output:
!  Glob_symmetry(Glob_Nperm): symmetry factor
!===========================
IMPLICIT NONE
INTEGER::i,j,k,mi,mj
REAL(dp)::temp

DO i=1,Glob_Nperm
  Glob_symmetry(i)=ZERO
  DO j=1,Glob_Nspin
    DO k=1,Glob_Nspin
      temp=Glob_permut_parity(i)*Glob_spin_parity(j)*Glob_spin_parity(k)
      DO mi=1,Glob_Nparticle
        mj=Glob_P(mi,i)
        temp=temp*delta(Glob_ptype(mi),Glob_ptype(mj))* &
        & delta(Glob_spin_config(mi,j),Glob_spin_config(mj,k))
      ENDDO
      Glob_symmetry(i)=Glob_symmetry(i)+temp
    ENDDO
  ENDDO
ENDDO

RETURN
END SUBROUTINE symmetry_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE write_symmetry_fun()
!==============================
!This program write symmetry factors into 
!symetry.txt file
!==============================
IMPLICIT NONE
INTEGER::ip
CHARACTER(100)::IRAst
IRAst="(f10.5,A1,????????I8,A5)"

OPEN(unit=20,file='symmetry.txt')
WRITE(IRAst(11:18),fmt="(TL1,I8.8)")Glob_Nparticle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO ip=1,Glob_Nperm
  IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    WRITE(20,IRAst)Glob_symmetry(ip),'(',Glob_P(1:Glob_Nparticle,ip),')'
  ENDIF
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CLOSE(20)
RETURN
END SUBROUTINE write_symmetry_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE permutation()
!===========================OK
!generate the permutation operator and its parity
!output:
!  Glob_P: all permutations
!  Glob_Tp: permutations acting on relative coordinate
!  Glob_permut_parity: parity of each permutation 
!===========================
IMPLICIT NONE
INTEGER::ip,i,j,k,l,pk
REAL(dp)::iparity
INTEGER::ipermut(Glob_Nparticle)

loop:DO ip=1,Glob_Nperm
    
  CALL permut_fun(ip,Glob_Nparticle,ipermut)
  CALL permut_parity_fun(ipermut,Glob_Nparticle,iparity)
  Glob_permut_parity(ip)=iparity
  Glob_P(:,ip)=ipermut(:)  ! ip_th permutation among Glob_Nparticle particles

!permutation Glob_Tp acting on Jaccobi coordinate x
  DO i=1,Glob_Np
    DO j=1,Glob_Np
      Glob_Tp(i,j,ip)=ZERO
	  DO k=1,Glob_Nparticle
        pk=ipermut(k)
	    Glob_Tp(i,j,ip)=Glob_Tp(i,j,ip)+Glob_U(i,k)*Glob_inv_U(pk,j)
	  ENDDO
    ENDDO
  ENDDO
  
ENDDO loop

END SUBROUTINE permutation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE permut_fun(k,N,ipermut)
!============================OK
!generate all permutations of N particles
!(all rearangement of N numbers from 1 to N)
!input:
!  k:the k_th permutation
!  N:particle number
!output:
!  ipermut:column vector containing the k_th permutation(ipermut(i)->j)
!============================
IMPLICIT NONE
INTEGER::k,N
INTEGER::ipermut(N)
INTEGER::i,j,is,io
INTEGER::iv(N+1)

ipermut=ZERO
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
ipermut(N)=IV(1)
RETURN

END SUBROUTINE permut_fun
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
	  
SUBROUTINE permut_parity_fun(ipermut,N,iparity)
!============================OK
!the parity of permutation ia
!input:
!  ipermut:permutation
!   N:particle number
!output:
!  iparity:parity of ipermut permutation
!============================
IMPLICIT NONE
INTEGER::N
INTEGER::ipermut(N)
INTEGER::itemp(N)
INTEGER::i,j,M,temp
REAL(dp)::iparity

DO i=1,N
  itemp(i)=ipermut(i)
ENDDO
!!!!!!!!!!!!!!!!!!!!
M=0
DO i=1,N
  DO j=i,N
    IF(itemp(j)==i .AND. i/=j)THEN
      temp=itemp(i)
      itemp(i)=itemp(j)
      itemp(j)=temp
      M=M+1
    ENDIF   
  ENDDO
ENDDO  
iparity=(-ONE)**M 
RETURN
END SUBROUTINE permut_parity_fun
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION fac(k)
!======================OK
!calculate the factorial of k:  k!
!======================
IMPLICIT NONE
INTEGER::k,fac
INTEGER::i
fac=1
DO i=1,k
    fac=fac*i
ENDDO
RETURN
END FUNCTION fac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION dfac(k)
!====================OK
!calculate the double-factorial of k:  k!!
!====================
IMPLICIT NONE
INTEGER::k,dfac
INTEGER::i
i=k
dfac=k
DO WHILE(i>2)
i=i-2
dfac=dfac*i
ENDDO
RETURN
END FUNCTION dfac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
FUNCTION delta(i,j)
!===============================OK
!delta function
!===============================
IMPLICIT NONE
INTEGER::i,j
REAL(dp)::delta
delta=ZERO
IF(i==j)THEN 
delta=ONE
ENDIF
RETURN
END FUNCTION delta   
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE symmetry