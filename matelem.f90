MODULE matelem
!======================================================
!this module contains subroutines computing the matrix
!element of Hamiltonian H and overlap S
!======================================================
USE globvars
USE auxfun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!define private matrix and variables used in module matelem.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================================
!private variable used when calculating matrix element(marked by p_ )
!these variables require allocations once at the begining
!of the whole program. 
!be aware of the order of calculating each variable.
!=================================
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables used for all basis form
!=============================
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_Lk,p_Ll
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_Ak,p_Al
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tAl,p_tAkl
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_inv_tAkl
REAL(dp),PRIVATE::p_det_tAkl
!!!!!!!!!!!!!!!!!!!!!
REAL(dp),PRIVATE::p_Skl,p_Tkl,p_Vkl_Coulomb
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dSkl_dLk,p_dHkl_dLk
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dTkl_dLk,p_dVkl_Coulomb_dLk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables only used for basis with L=0,M=0
!=============================
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tr_invtAkl_Jij 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables only used for basis with L=0,M=0
!=============================
INTEGER,PRIVATE::p_mk,p_ml
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_vk,p_vl,p_tvl
REAL(dp),PRIVATE::p_tau1,p_tau2,p_tau3
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_eta1,p_eta2
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tKkl
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_inv_Akk,p_inv_All
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tAl_Lambda_Ak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate private variables used in module matelem.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE pvars_matelem()
!==============================
!allocate matrix used in module matelem.f90
!only need allocate once at the begining of the whole program
!since they are used in the whole time
!==============================
IMPLICIT NONE

SELECT CASE(Glob_basis_form)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)

ALLOCATE(p_Lk(Glob_Np,Glob_Np),p_Ll(Glob_Np,Glob_Np))
ALLOCATE(p_Ak(Glob_Np,Glob_Np),p_Al(Glob_Np,Glob_Np))
ALLOCATE(p_tAl(Glob_Np,Glob_Np),p_tAkl(Glob_Np,Glob_Np))
ALLOCATE(p_inv_tAkl(Glob_Np,Glob_Np))
ALLOCATE(p_dSkl_dLk(Glob_NLk))
ALLOCATE(p_dTkl_dLk(Glob_NLk))
ALLOCATE(p_dVkl_Coulomb_dLk(Glob_NLK))
ALLOCATE(p_dHkl_dLk(Glob_NLk))
!!!!!!!!!!!!!!!!!!!!
ALLOCATE(p_tr_invtAkl_Jij(Glob_Nparticle,Glob_Nparticle))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)

ALLOCATE(p_Lk(Glob_Np,Glob_Np),p_Ll(Glob_Np,Glob_Np))
ALLOCATE(p_Ak(Glob_Np,Glob_Np),p_Al(Glob_Np,Glob_Np))
ALLOCATE(p_tAl(Glob_Np,Glob_Np),p_tAkl(Glob_Np,Glob_Np))
ALLOCATE(p_inv_tAkl(Glob_Np,Glob_Np))
ALLOCATE(p_dSkl_dLk(Glob_NLk))
ALLOCATE(p_dTkl_dLk(Glob_NLk))
ALLOCATE(p_dVkl_Coulomb_dLk(Glob_NLk))
ALLOCATE(p_dHkl_dLk(Glob_NLk))
!!!!!!!!!!!!!!!!!!!!
ALLOCATE(p_vk(Glob_Np),p_vl(Glob_Np),p_tvl(Glob_Np))
ALLOCATE(p_eta1(Glob_Nparticle,GLob_Nparticle))
ALLOCATE(p_eta2(Glob_Nparticle,GLob_Nparticle))
ALLOCATE(p_tKkl(Glob_Np,Glob_Np))
ALLOCATE(p_inv_Akk(Glob_Np,Glob_Np),p_inv_All(Glob_Np,Glob_Np))
ALLOCATE(p_tAl_Lambda_Ak(Glob_Np,Glob_Np))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT 

RETURN
END SUBROUTINE pvars_matelem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!matrix element for Hamilton and overlap for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HS_0_0_fun(bk,bl,Skl,Hkl)
!===============================OK
!calculate the matrix element of Hamitonian
!for basis with L=0,M=0
!input:
!  bk,bl: integer index
!output:
!  Skl: overlap of bk_th and bl_th basis
!  Hkl: ovarlap of Hamitonian between the bk_th and bl_th basis
!===============================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::Skl,Hkl

INTEGER::i,j,k,L,ip
REAL(dp)::temp

Skl=ZERO
Hkl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters Lk,Ll,Ak,Al from Glob_LK and Glob_mk
CALL paras_0_0_fun(bk,bl)  

!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
    
!=====================
!if symmetry factor=0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 

!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl,Tkl,Vkl
CALL Skl_0_0_fun()
CALL Tkl_0_0_fun()
CALL Vkl_Coulomb_0_0_fun()

!Skl,Hkl
Skl=Skl+p_Skl*Glob_symmetry(ip)
Hkl=Hkl+(p_Tkl+p_Vkl_Coulomb)*Glob_symmetry(ip)

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE HS_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient matrix element for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gradHS_0_0_fun(bk,bl,dSkl_dLk,dHkl_dLk)
!=============================
!calculate the gradient of overlap element
!for basis with L=0,M=0
!input:
!  bk,bl: basis index
!output:
!  dSkl_dLk: gradient of overlap Skl
!  dHkl_dLk: gradient of Hamilton
!=============================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::dSkl_dLk(Glob_NLk)
REAL(dp),INTENT(OUT)::dHkl_dLk(Glob_NLk)
  
INTEGER::i,j,k,L,ip
REAL(dp)::temp

dSkl_dLk=ZERO
dHkl_dLk=ZERO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_0_0_fun(bk,bl)  
  
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
  
!=====================
!if symmetry factor=0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
  
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
  
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
  
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl,Tkl,Vkl
CALL Skl_0_0_fun()
CALL Tkl_0_0_fun()
CALL Vkl_Coulomb_0_0_fun()

!dSkl_dLk,dTkl_dLk,dVkl_Coulomb_dLk
CALL dSkl_dLk_0_0_fun(ip,bk,bl)
CALL dTkl_dLk_0_0_fun(ip,bk,bl)
CALL dVkl_Coulomb_dLk_0_0_fun(ip,bk,bl)
    
dSkl_dLk=dSkl_dLk+p_dSkl_dLk*Glob_symmetry(ip)
dHkl_dLk=dHkl_dLk+(p_dTkl_dLk+p_dVkl_Coulomb_dLk)*Glob_symmetry(ip)
  
ENDIF
ENDDO perm_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE  gradHS_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function Gij for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE  gij_kl_0_0_fun(bk,bl,ii,jj,R,gij)
!===================================
!this subroutine calculate the pair correlation between 
!bk and bl basis gij(bk,bl)=<bk|delta(r-x)|bl>
!input: 
! bk,bl: basis index
! ii,jj: particle index
!     R: radius
!output:
!   gij: the pair corelation 
!===================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(IN)::ii,jj
REAL(dp),INTENT(IN)::R
REAL(dp),INTENT(OUT)::gij
  
INTEGER::i,j,k,L,ip
REAL(dp)::temp,temp1,temp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters Lk,Ll,Ak,Al from Glob_Lk and Glob_mk
CALL paras_0_0_fun(bk,bl)  
    
gij=ZERO
  
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
        
!=====================
!if symmetry factor < 0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
    
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!temp1=tr[in_tAkl*Jij]
temp1=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp1=temp1+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
temp1=dabs(temp1)    
  
!temp2=exp(-x^2/tr[inv_tAkl*Jij])
temp2=dexp(-R**2/temp1)
  
!Skl,Tkl,Vkl
CALL Skl_0_0_fun()
  
!Skl
gij=gij+p_Skl*Glob_symmetry(ip)*temp2/(PI*SQRTPI)/(temp1*dsqrt(temp1))

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gij_kl_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!expectation value of interparticle distance for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_kl_0_0_fun(bk,bl,ii,jj,rij)
!=================================
!this subroutine calculate the expectation value of interparticle distance
!for basis with L=0,M=0
  !input: 
! bk,bl: basis index
! ii,jj: particle index
!output:
!   rij: expectation value of interparticle diatance 
!=================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(IN)::ii,jj
REAL(dp),INTENT(OUT)::rij

INTEGER::i,j,k,L,ip
REAL(dp)::temp

rij=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_0_0_fun(bk,bl)  
    
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
        
!=====================
!if symmetry factor < 0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
    
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl
CALL Skl_0_0_fun()

!temp=tr[in_tAkl*Jij]
temp=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
temp=dabs(temp)   

!<rij>
rij=rij+p_Skl*(TWO/SQRTPI)/temp

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE rij_kl_0_0_fun


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap check function for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_0_0_check_fun(bk,bl,ERR)
!========================
!ckeck the overlap of bk-bl basis
!for basis with L=0,M=0 
!input:
!bk,bl: basis index
!output:
!  ERR: 
!    0: within overlap threshold
!    1: exceed overlap threshold
!=========================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L
INTEGER::ip
REAL(dp)::temp

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_0_0_fun(bk,bl)  

!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
    
!=====================
!if symmetry factor=0,then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 

!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl
CALL Skl_0_0_fun()

IF((dabs(p_Skl)>Glob_overlap_threshold).AND.(bk/=bl))THEN
  ERR=1
  EXIT perm_loop
ENDIF

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_0_0_check_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!parameter extraction function for L=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE paras_0_0_fun(bk,bl)
!================================OK
!extract related parameter matrix of the 
!k_th and l_th basis with L=0,M=0
!input:
!  bk,bl
!output:
!  p_Lk,p_Ll,p_Ak,p_Al
!=================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,index
REAL(dp)::vechLk(Glob_NLk),vechLl(Glob_NLk)
REAL(dp)::temp1,temp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extrat Lk and Ll from Glob_Lk
vechLk(:)=Glob_Lk(:,bk)
vechLl(:)=Glob_Lk(:,bl)

!allocate Lk and Ll
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    p_Lk(j,i)=ZERO
    p_Ll(j,i)=ZERO
    p_Lk(i,j)=vechLk(index) !lower trangle 
    p_Ll(i,j)=vechLl(index) !lower trangle
  ENDDO
ENDDO

!form Ak,Al matrix
DO i=1,Glob_Np
  DO j=i,Glob_Np
   temp1=ZERO
   temp2=ZERO
   DO k=1,Glob_Np
    temp1=temp1+p_Lk(i,k)*p_Lk(j,k)
    temp2=temp2+p_Ll(i,k)*p_Ll(j,k)
   ENDDO
   p_Ak(i,j)=temp1
   p_Ak(j,i)=temp1
   p_Al(i,j)=temp2
   p_Al(j,i)=temp2
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE paras_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap Skl and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_0_0_fun()
!==================================
!calculate the overlap
!output:
! p_Skl: overlap 
!==================================
IMPLICIT NONE
INTEGER::i,j
REAL(dp)::det_Lk,det_Ll
REAL(dp)::temp

p_Skl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!det_Lk,det_Ll
det_Lk=ONE
det_Ll=ONE
DO i=1,Glob_Np
  det_Lk=det_Lk*p_Lk(i,i)
  det_Ll=det_Ll*p_Ll(i,i)
ENDDO

!Skl
temp=dabs(det_Lk*det_Ll)/p_det_tAkl
IF(dabs(p_det_tAkl)<=ZERO)PAUSE'det_tAkl==0'

p_Skl=(TWO**(THREE*Glob_Np/TWO))*temp*dsqrt(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dSkl_dLk_0_0_fun(ip,bk,bl)
!===============================OK
!calculate gradient of overlap Skl
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
!  p_dSkl_dLk: gradient of overlap Skl
!===============================
INTEGER,INTENT(IN)::ip,bk,bl
  
INTEGER::i,j,k,L,index
REAL(dp)::temp
REAL(dp)::inv_Lk(Glob_Np,Glob_Np)
REAL(dp)::inv_tAkl_ttAkl(Glob_Np,Glob_Np)
  
p_dSkl_dLk=ZERO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!bk==bl ->inv_tAkl_ttAkl=inv_tAkl+inv_ttAkl
!bk!=bl ->inv_tAkl_ttAkl=inv_tAkl
!===========================
  
!inv_Lk
CALL inv_L_lower(Glob_Np,p_Lk,inv_Lk)
  
!(bk!=bl)
inv_tAkl_ttAkl=p_inv_tAkl
  
!(bk==bl)
IF(bk==bl)THEN
      
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*p_inv_tAkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
  inv_tAkl_ttAkl(i,j)=inv_tAkl_ttAkl(i,j)+temp
  ENDDO
ENDDO 
  
inv_Lk=TWO*inv_Lk
  
ENDIF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
!dSkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
      temp=temp+inv_tAkl_ttAkl(i,k)*p_Lk(k,j)
    ENDDO
    index=index+1
    p_dSkl_dLk(index)=(THREE/TWO)*p_Skl*(inv_Lk(j,i)-TWO*temp)
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dSkl_dLk_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kinetic T_kl and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tkl_0_0_fun()
!==============================
!calculate the overlap of kinetic energy 
!output:
!  p_Tkl: kinetic energy 
!==============================
IMPLICIT NONE
INTEGER::i,j,k
REAL(dp)::temp
REAL(dp)::temp1(Glob_Np,Glob_Np)
REAL(dp)::temp2(Glob_Np,Glob_Np)

p_Tkl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!inv_tAkl*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_Ak(k,j)
    ENDDO
    temp1(i,j)=temp
  ENDDO
ENDDO

!tAl*inv_tAkl*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+p_tAl(i,k)*temp1(k,j)
    ENDDO
    temp2(i,j)=temp
  ENDDO
ENDDO

!Tkl=6*S_kl*tr[Lambda*tAl*inv_tAkl*Ak]
temp=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+Glob_Lambda(i,k)*temp2(k,i)
  ENDDO
ENDDO
p_Tkl=SIX*p_Skl*temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Tkl_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dTkl_dLk_0_0_fun(ip,bk,bl)
!==================================
!calculate the gradient of kinetic energy
!input:
!    ip: permutation index
! bk,bl: basis index
!output:
!p_dTkl_dLl: gradient of kinetic energy 
!==================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl
  
INTEGER::i,j,k,L,index
REAL(dp)::temp
REAL(dp)::temp1(Glob_Np,Glob_Np)
REAL(dp)::temp2(Glob_Np,Glob_Np)
REAL(dp)::F_G(Glob_Np,Glob_Np)
  
p_dTkl_dLk=ZERO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!F_G: bk==bl -> F_G=F+G
!     bk!=bl -> F_G=F
!===============================
  
!inv_tAkl*tAl*Lambda
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,GLob_Np
       temp=temp+p_inv_tAkl(i,k)*p_tAl(k,L)*Glob_Lambda(L,j)
      ENDDO
    ENDDO
    temp1(i,j)=temp
  ENDDO
ENDDO

!F=inv_tAkl*tAl*Lambda*tAl*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,GLob_Np
        temp=temp+temp1(i,k)*p_tAl(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    F_G(i,j)=temp
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!(bk=bl)
IF(bk==bl)THEN

!inv_tAkl*Ak*Lambda
DO j=1,GLob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*p_Ak(k,L)*Glob_Lambda(L,j)
      ENDDO
    ENDDO
    temp1(i,j)=temp
  ENDDO
ENDDO

!inv_tAkl*Ak*Lambda*Ak*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+temp1(i,k)*p_Ak(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    temp2(i,j)=temp
  ENDDO
ENDDO

!G=Tp(inv_tAkl*Ak*Lambda*Ak*inv_tAkl)Tp'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*temp2(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    F_G(i,j)=F_G(i,j)+temp
  ENDDO
ENDDO
  
ENDIF
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dTkl_dLl
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
      temp=temp+F_G(i,k)*p_Lk(k,j)
    ENDDO
    index=index+1
    p_dTkl_dLk(index)=(p_Tkl/p_Skl)*p_dSkl_dLk(index)&
    &+THREE*FOUR*p_Skl*temp
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dTkl_dLk_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!coulumb potential and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Vkl_Coulomb_0_0_fun()
!================================
!calculate Coulomb potential energy
!output:
!  p_Vkl_Coulomb
!================================
IMPLICIT NONE
INTEGER::i,j,k
INTEGER::ii,jj
REAL(dp)::temp,Vkl_1

p_Vkl_Coulomb=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!interactive Coulomb potential between each particle: Vkl_1
!============================

Vkl_1=ZERO

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!Tr[inv_tAkl*Jij]
temp=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
temp=dabs(temp)            
p_tr_invtAkl_Jij(ii,jj)=temp   !store Tr[inv_tAkl*Jij] 

Vkl_1=Vkl_1+Glob_charge(ii)*Glob_charge(jj)/dsqrt(temp)
   
ENDDO loop2
ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==================================
!sum of all Coulomb potantial 
!==================================

p_Vkl_Coulomb=(TWO/SQRTPI)*p_Skl*Vkl_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Vkl_Coulomb_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dVkl_Coulomb_dLk_0_0_fun(ip,bk,bl)
!==============================OK
!calculate the gradient of Coulomb potential energy 
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
!  p_dVkl_Coulomb_dL1:gradient of Coulomb potential
!==============================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
INTEGER::ii,jj
REAL(dp)::temp
REAL(dp)::dVkl_1_dLk(Glob_NLk)
REAL(dp)::temp1(Glob_Np,Glob_Np)
REAL(dp)::tQ_ttQ(Glob_Np,Glob_Np)

p_dVkl_Coulomb_dLk=ZERO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!interactive Coulomb potential between each particle: dVkl_1_dLk
!=====================

dVkl_1_dLk=ZERO

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle
    
!tQ=inv_tAkl*Jij*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    tQ_ttQ(i,j)=temp
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!(bk=bl)
IF(bk==bl)THEN
    
!ttQ=P*tQ*P'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tQ_ttQ(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    temp1(i,j)=temp
  ENDDO
ENDDO
  
!tQ_ttQ=tQ+ttQ
tQ_ttQ=tQ_ttQ+temp1
  
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!dVkl_1dLK=sum{qi*qj*tr[inv_tAkl*Jij]^(-3/2)*vech(tQ_ttQ*Lk)}
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
     temp=temp+tQ_ttQ(i,k)*p_Lk(k,j)
    ENDDO
    index=index+1
    dVkl_1_dLk(index)=dVkl_1_dLk(index)+Glob_charge(ii)*Glob_charge(jj)&
    &*(p_tr_invtAkl_Jij(ii,jj)**(-THREE/TWO))*temp
  ENDDO
ENDDO

ENDDO loop2
ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=================
!sum of all Coulomb potential gradient
!=================

p_dVkl_Coulomb_dLk=(p_Vkl_Coulomb/p_Skl)*p_dSkl_dLk+&
&(TWO/SQRTPI)*p_Skl*dVkl_1_dLk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dVkl_Coulomb_dLk_0_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Hamilton and overlap for L=1,M=0 basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HS_1_0_fun(bk,bl,Skl,Hkl)
!===================================
!Hamiton and overlap matrix element for basis with L=1,M=0
!input:
!  bk,bl: basis index
!output:
!  Skl: overlap
!  Hkl: Hamilton 
!===================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::Skl,Hkl

INTEGER::i,j,k,L,ip
REAL(dp)::temp

Skl=ZERO
Hkl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract mk,ml and matrix Lk,Ll,Ak,Al from GLob_Lk
CALL paras_1_0_fun(bk,bl)

perm_loop:DO ip=1,GLob_Nperm

!=====================
!if symmetry factor=0, then the corresbonding terms need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 

!tvl
DO i=1,Glob_Np  
  p_tvl(i)=Glob_Tp(p_ml,i,ip)
ENDDO

!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl,Tkl,Vkl
CALL Skl_1_0_fun()
CALL Tkl_1_0_fun()
CALL Vkl_Coulomb_1_0_fun()

!Skl,Hkl
Skl=Skl+p_Skl*Glob_symmetry(ip)
Hkl=Hkl+(p_Tkl+p_Vkl_Coulomb)*Glob_symmetry(ip)

ENDIF
ENDDO perm_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE HS_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient of Hamiton and overlap for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gradHS_1_0_fun(bk,bl,dSkl_dLk,dHkl_dLk)
!==============================
!calculate the gradient of Hamilton and overlap 
!for basis with L=1,M=0
!input:
! bk,bl: basis index
!output:
! dSkl_dLk: gradient of overlap
! dHkl_dLk: gradient of Hamiton
!==============================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::dSkl_dLk(Glob_NLk)
REAL(dp),INTENT(OUT)::dHkl_dLk(Glob_NLk)
    
INTEGER::i,j,k,L,ip
REAL(dp)::temp
  
dSkl_dLk=ZERO
dHkl_dLk=ZERO
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_0_0_fun(bk,bl)  
    
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
    
!=====================
!if symmetry factor=0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
    
!tvl
DO i=1,Glob_Np  
  p_tvl(i)=Glob_Tp(p_ml,i,ip)
ENDDO

!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
  
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
    
!Skl,Tkl,Vkl
CALL Skl_1_0_fun()
CALL Tkl_1_0_fun()
CALL Vkl_Coulomb_1_0_fun()
  
!dSkl_dLk,dTkl_dLk,dVkl_Coulomb_dLk
CALL dSkl_dLk_1_0_fun(ip,bk,bl)
CALL dTkl_dLk_1_0_fun(ip,bk,bl)
CALL dVkl_Coulomb_dLk_1_0_fun(ip,bk,bl)

dSkl_dLk=dSkl_dLk+p_dSkl_dLk*Glob_symmetry(ip)
dHkl_dLk=dHkl_dLk+(p_dTkl_dLk+p_dVkl_Coulomb_dLk)*Glob_symmetry(ip)

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gradHS_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gij_kl_1_0_fun(bk,bl,ii,jj,Rxy,Rz,gij)
!==================================
!this subroutine calculates the correaltion function at x
!for basis with L=1,M=0
!input: 
! bk,bl: basis index
! ii,jj: particle index
!   Rxy: sqrt(x**2+y**2)
!    Rz: z component of R
!output:
!   gij: correlation function coresbonding to (R,z)
!==================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(IN)::ii,jj
REAL(dp),INTENT(IN)::Rxy,Rz
REAL(dp),INTENT(OUT)::gij

INTEGER::i,j,k,L,ip
REAL(dp)::eta1,eta2,temp,temp1,temp2

REAL(dp)::R
R=dsqrt(Rxy**2+Rz**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_1_0_fun(bk,bl)  
    
gij=ZERO
  
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
        
!=====================
!if symmetry factor=0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
  
!tvl
DO i=1,Glob_Np  
  p_tvl(i)=Glob_Tp(p_ml,i,ip)
ENDDO
    
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl and tau3
CALL Skl_1_0_fun()

!eta1
eta1=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta1=eta1+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!eta2
eta2=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta2=eta2+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
  
temp1=(ONE+(eta2/eta1/p_tau3)*(TWO*(Rz**2)/eta1-ONE))
temp2=dexp(-R**2/eta1)

gij=gij+Glob_symmetry(ip)*p_Skl*((PI*eta1)**(-THREE/TWO))*&
&temp1*temp2

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gij_kl_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!expectation value of interparticle distace for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_kl_1_0_fun(bk,bl,ii,jj,rij)
!==================================
!this subroutine calculate the expectation value of 
!interparticle diatance for basis with L=1,M=0
!input:
!  bk,bl: basis index
!  ii,jj: particle index
!output:
!    rij: expectation value of interparticle diatance
!==================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(IN)::ii,jj
REAL(dp),INTENT(OUT)::rij

INTEGER::i,j,k,L,ip
REAL(dp)::eta1,eta2,temp

rij=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters Lk,Ll,Ak,Al from Glob_Lk
CALL paras_1_0_fun(bk,bl)  
  
!sum of symmetry loop
perm_loop:DO ip=1,Glob_Nperm
        
!=====================
!if symmetry factor=0, then the corresbonding term need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 
  
!tvl
DO i=1,Glob_Np  
  p_tvl(i)=Glob_Tp(p_ml,i,ip)
ENDDO
    
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl and tau3
CALL Skl_1_0_fun()

!eta1
eta1=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta1=eta1+p_inv_tAKl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!eta2
eta2=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta2=eta2+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

temp=ONE+(eta2/eta1/p_tau3)/THREE

rij=rij+Glob_symmetry(ip)*p_Skl*(TWO/SQRTPI)*temp*dsqrt(eta1)

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE rij_kl_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap threshold check fun for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_1_0_check_fun(bk,bl,ERR)
!=================================
!check wheather the overlap have exceeded the threshold 
!input:
!  bk,bl: basis index
!=================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L,ip
REAL(dp)::Skl,temp

Skl=ZERO
ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract mk,ml and matrix Lk,Ll,Ak,Al from GLob_Lk
CALL paras_1_0_fun(bk,bl)

perm_loop:DO ip=1,Glob_Nperm

!=====================
!if symmetry factor=0, then the corresbonding terms need not be calculated
!=====================
IF(dabs(Glob_symmetry(ip))>ZERO)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(k,i,ip)*p_Al(k,L)*Glob_Tp(L,j,ip)
      ENDDO
    ENDDO
    p_tAl(i,j)=temp
    p_tAl(j,i)=temp
  ENDDO
ENDDO 

!tvl
DO i=1,Glob_Np  
  p_tvl(i)=Glob_Tp(p_ml,i,ip)
ENDDO
  
!tAkl=Ak+tAl
p_tAkl=p_Ak+p_tAl
  
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl
CALL Skl_1_0_fun()

IF((dabs(p_Skl)>Glob_overlap_threshold).AND.(bk/=bl))THEN
  ERR=1
  EXIT perm_loop
ENDIF

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_1_0_check_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!parameters extraction function for L=1,M=0 basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE paras_1_0_fun(bk,bl)
!===================================OK
!extract related parameter matrix of the 
!k_th and l_th basis with L=1,M=0
!input:
!  bk,bl: basis index
!output:
!  p_mk,p_ml,p_vk,p_vl,p_Lk,p_Ll,p_Ak,p_Al
!===================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,index
REAL(dp)::vechLk(Glob_NLk),vechLl(Glob_NLk)
REAL(dp)::temp1,temp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract mk and ml from GLob_mk
p_mk=Glob_mk(1,bk)
p_ml=Glob_mk(1,bl)

p_vk=ZERO
p_vl=ZERO
p_vk(p_mk)=ONE
p_vl(p_ml)=ONE


!extract Lk and Ll from Glob_Lk
vechLk(:)=Glob_Lk(1:GLob_NLk,bk)
vechLl(:)=Glob_Lk(1:GLob_NLk,bl)

!allocate Lk and Ll
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    p_Lk(j,i)=ZERO
    p_Ll(j,i)=ZERO
    p_Lk(i,j)=vechLk(index) !lower trangle 
    p_Ll(i,j)=vechLl(index) !lower trangle
  ENDDO
ENDDO

!form Ak,Al matrix
DO i=1,Glob_Np
  DO j=i,Glob_Np
   temp1=ZERO
   temp2=ZERO
   DO k=1,Glob_Np
    temp1=temp1+p_Lk(i,k)*p_Lk(j,k)
    temp2=temp2+p_Ll(i,k)*p_Ll(j,k)
   ENDDO
   p_Ak(i,j)=temp1
   p_Ak(j,i)=temp1
   p_Al(i,j)=temp2
   p_Al(j,i)=temp2
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE paras_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap S_kl and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_1_0_fun()
!=============================OK
!overlap element with L=1,M=0
!output:
!  p_Skl
!=============================
IMPLICIT NONE

INTEGER::i,j,k,L
REAL(dp)::temp1,temp2
REAL(dp)::det_Lk,det_Ll
REAL(dp)::inv_Lk(Glob_Np,Glob_Np),inv_Ll(Glob_Np,Glob_Np)

p_Skl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!det_Lk,det_Ll
det_Lk=ONE
det_Ll=ONE
DO i=1,Glob_Np
  det_Lk=det_Lk*p_Lk(i,i)
  det_Ll=det_Ll*p_Ll(i,i)
ENDDO

!inv_Lk,inv_Ll
CALL inv_L_lower(Glob_Np,p_Lk,inv_Lk)
CALL inv_L_lower(Glob_Np,p_Ll,inv_Ll)

!inv_Akk,inv_All
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp1=ZERO
    temp2=ZERO
    DO k=1,Glob_Np
      temp1=temp1+inv_Lk(k,i)*inv_Lk(k,j)
      temp2=temp2+inv_Ll(k,i)*inv_Ll(k,j)
    ENDDO
    temp1=temp1/TWO
    temp2=temp2/TWO
    p_inv_Akk(i,j)=temp1
    p_inv_Akk(j,i)=temp1
    p_inv_All(i,j)=temp2
    p_inv_All(j,i)=temp2
  ENDDO
ENDDO

!vk'*inv_Akk*vk
temp1=p_inv_Akk(p_mk,p_mk)

!vl'*inv_All*vl
temp2=p_inv_All(p_ml,p_ml)

!tau3=tr[tvl*vk'*inv_tAkl]=vk'*inv_tAkl*tvl
p_tau3=ZERO
DO i=1,Glob_Np
  p_tau3=p_tau3+p_inv_tAkl(p_mk,i)*p_tvl(i)
ENDDO

!Skl
p_Skl=(TWO**(THREE*Glob_Np/TWO))*&
&((dabs(det_Lk*det_Ll)/p_det_tAkl)**(THREE/TWO))*&
&(p_tau3/dsqrt(temp1*temp2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dSkl_dLk_1_0_fun(ip,bk,bl)
!=============================OK
!gradient of overlap Skl
!for basis with L=1,M=0
!input:
!     ip: permutation index
!  bk,bl: basis index
!output:
! p_dSkl: gradient of Skl
!=============================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
REAL(dp)::temp,temp1(Glob_Np,Glob_Np)
REAL(dp)::tFkl_tGkl(Glob_Np,Glob_Np)
REAL(dp)::Fkk_Fll(Glob_Np,Glob_Np)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Fkk
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=p_inv_Akk(i,p_mk)*p_inv_Akk(p_mk,j)/p_inv_Akk(p_mk,p_mk)&
    &+THREE*p_inv_Akk(i,j)/TWO
    Fkk_Fll(i,j)=temp
    Fkk_Fll(j,i)=temp
  ENDDO
ENDDO

!tFkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    tFkl_tGkl(i,j)=p_tKkl(i,j)/p_tau3+THREE*p_inv_tAkl(i,j)/TWO
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!！
IF(bk==bl)THEN

!Fkk+Fll
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=p_inv_All(i,p_ml)*p_inv_All(p_ml,j)/p_inv_All(p_ml,p_ml)&
    &+THREE*p_inv_All(i,j)/TWO
    Fkk_Fll(i,j)=Fkk_Fll(i,j)+temp
    Fkk_Fll(j,i)=Fkk_Fll(j,i)+temp
  ENDDO
ENDDO
  
!tGkl=Tp*tFkl*Tp'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tFkl_tGkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    temp1(i,j)=temp
  ENDDO
ENDDO

!tFkl+tGkl
DO j=1,Glob_Np
  DO i=1,Glob_Np    
    tFkl_tGkl(i,j)=tFkl_tGkl(i,j)+temp1(i,j)
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dSkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
      temp=temp+(TWO*Fkk_Fll(i,k)-tFkl_tGkl(i,k)-tFkl_tGkl(k,i))*p_Lk(k,j)
    ENDDO
    index=index+1
    p_dSkl_dLk(index)=p_Skl*temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dSkl_dLk_1_0_fun  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kinetic energy and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tkl_1_0_fun()
!=============================OK
!kinetic energy 
!for basis with L=1,M=0
!output:
! p_Tkl: kinetic energy
!=============================
IMPLICIT NONE
INTEGER::i,j,k,L
REAL(dp)::temp

p_Tkl=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!tAl*Lambda*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_tAl(i,k)*Glob_Lambda(k,L)*p_Ak(L,j)
      ENDDO
    ENDDO
    p_tAl_Lambda_Ak(i,j)=temp
  ENDDO
ENDDO

!tau1 and tau2
p_tau1=ZERO
p_tau2=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    p_tau1=p_tau1+p_inv_tAkl(i,k)*p_tAl_Lambda_Ak(k,i)
    p_tau2=p_tau2+p_tKkl(i,k)*p_tAl_Lambda_Ak(k,i)
  ENDDO
ENDDO

!Tkl
p_Tkl=p_Skl*(SIX*p_tau1+FOUR*p_tau2/p_tau3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Tkl_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dTkl_dLk_1_0_fun(ip,bk,bl)
!===============================OK
!gradient of kinetic energy
!for basis with L=1,M=0
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
! p_dTkl_dLk:gradient of kinetic energy
!==============================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
REAL(dp)::temp,temp1,temp2,temp3(Glob_Np,Glob_Np)
REAL(dp)::tAl_Lambda_tAl(Glob_Np,Glob_Np)
REAL(dp)::Ak_Lambda_Ak(Glob_Np,Glob_Np)
REAL(dp)::tUkl_tWkl(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!tAl*Lambda*tAl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_tAl(i,k)*Glob_Lambda(k,L)*p_tAl(L,j)
      ENDDO
    ENDDO
    tAl_Lambda_tAl(i,j)=temp
  ENDDO
ENDDO

!tUkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    temp1=ZERO
    temp2=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*tAl_Lambda_tAl(k,L)*p_inv_tAkl(L,j)
        temp1=temp1+p_tKkl(i,k)*tAl_Lambda_tAl(k,L)*p_inv_tAkl(L,j)
        temp2=temp2+p_inv_tAkl(i,k)*p_tAl_Lambda_Ak(k,L)*p_tKkl(L,j)
      ENDDO
    ENDDO
    tUkl_tWkl(i,j)=SIX*temp+(FOUR/p_tau3)*(temp1-temp2)&
    &+FOUR*p_tau2*p_tKkl(i,j)/(p_tau3**2)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN

!Ak*Lambda*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_Ak(i,k)*Glob_Lambda(k,L)*p_Ak(L,j)
      ENDDO
    ENDDO
    Ak_Lambda_Ak(i,j)=temp
  ENDDO
ENDDO

!tWkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    temp1=ZERO
    temp2=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Ak_Lambda_Ak(k,L)*p_inv_tAkl(L,j)
        temp1=temp1+p_inv_tAkl(i,k)*Ak_Lambda_Ak(k,L)*p_tKkl(L,j)
        temp2=temp2+p_tKkl(i,k)*p_tAl_Lambda_Ak(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    temp3(i,j)=SIX*temp+(FOUR/p_tau3)*(temp1-temp2)+&
    &FOUR*p_tau2*p_tKkl(i,j)/(p_tau3**2)
  ENDDO
ENDDO

!tUkl+tWkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*temp3(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    tUkl_tWkl(i,j)=tUkl_tWkl(i,j)+temp
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dTkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
      temp=temp+(tUkl_tWkl(i,k)+tUkl_tWkl(k,i))*p_Lk(k,j)
    ENDDO
    index=index+1
    p_dTkl_dLk(index)=p_Tkl*p_dSkl_dLk(index)/p_Skl+&
    & p_Skl*temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dTkl_dLk_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!coulumb potential energy and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Vkl_Coulomb_1_0_fun()
!============================OK
!Coulomb potential energy Vkl_Coulomb 
!for basis with L=1,M=0
!output:
!p_Vkl_Coulomb: Coulomb potential energy
!============================
IMPLICIT NONE
INTEGER::i,j,k,L
INTEGER::ii,jj
REAL(dp)::temp,Vkl_1

p_Vkl_Coulomb=ZERO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!interactive Coulomb potential between each particle: Vkl_1
!============================

Vkl_1=ZERO

loop1:DO ii=1,GLob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!eta1
temp=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_inv_tAKl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
p_eta1(ii,jj)=temp

!eta2
temp=ZERO
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
p_eta2(ii,jj)=temp

Vkl_1=Vkl_1+Glob_charge(ii)*Glob_charge(jj)*&
&(ONE-(p_eta2(ii,jj)/p_eta1(ii,jj)/p_tau3)/THREE)/dsqrt(p_eta1(ii,jj))

ENDDO loop2
ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!==========================
!sum of all Coulomb potential energy
!==========================

p_Vkl_Coulomb=(TWO/SQRTPI)*p_Skl*Vkl_1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Vkl_Coulomb_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dVkl_Coulomb_dLk_1_0_fun(ip,bk,bl)
!=================================OK
!the gradient of Coulomb potential energy
!for basis with L=1,M=0
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
!  p_dVkl_Coulomb_dLk: gradient of Coulomb potential energy
!================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
INTEGER::ii,jj
REAL(dp)::temp,temp1,temp2,temp3(Glob_Np,Glob_Np)
REAL(dp)::dVkl_1_dLk(Glob_NLk)
REAL(dp)::tQkl_tDkl(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!gradient of interactive Coulomb potential between each particle: dVkl_1_dLk
!============================

dVKL_1_dLk=ZERO

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!tQkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    temp1=ZERO
    temp2=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
        temp1=temp1+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_tKkl(L,j)
        temp2=temp2+p_tKkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    tQkl_tDkl(i,j)=((ONE-p_eta2(ii,jj)/p_eta1(ii,jj)/p_tau3)*temp/TWO&
    &+(temp1+temp2)/p_tau3/THREE&
    &-p_eta2(ii,jj)*p_tKkl(i,j)/THREE/(p_tau3**2))&
    &/(p_eta1(ii,jj)**(THREE/TWO))
  ENDDO
ENDDO

IF(bk==bl)THEN

!tDkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=ZERO
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tQkl_tDkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    temp3(i,j)=temp
  ENDDO
ENDDO

!tQkl+tDkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    tQkl_tDkl(i,j)=tQkl_tDkl(i,j)+temp3(i,j)
  ENDDO
ENDDO

ENDIF

!dVkl_1_dLK
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=ZERO
    DO k=j,Glob_Np
      temp=temp+(tQkl_tDkl(i,k)+tQkl_tDkl(k,i))*p_Lk(k,j)
    ENDDO
    index=index+1
    dVkl_1_dLK(index)=dVkl_1_dLk(index)&
    &+Glob_charge(ii)*Glob_charge(jj)*temp
  ENDDO
ENDDO

ENDDO loop2
ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!sum of all gradient of Coulomb potential
!============================

p_dVkl_Coulomb_dLk=(p_Vkl_Coulomb/p_Skl)*p_dSkl_dLk+&
(TWO/SQRTPI)*p_Skl*dVkl_1_dLk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dVkl_Coulomb_dLk_1_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE matelem
