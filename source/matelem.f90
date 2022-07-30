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
!================================================
!private variable used when calculating matrix element
!(marked by p_ ). These variables require allocations 
!once at the begining of the whole program. 
!================================================
IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables used for all basis form
!=============================
REAL(dp),PRIVATE::p_det_tAkl
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_Lk,p_Ll,p_Ak,p_Al
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tAl,p_tAkl,p_inv_tAkl
!!!!!!!!!!!!!!!!!!!!!
REAL(dp),PRIVATE::p_Skl,p_Tkl,p_Vkl_Coulomb
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dSkl_dLk,p_dHkl_dLk
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_dTkl_dLk,p_dVkl_coulomb_dLk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables only used for basis with L=0,M=0
!=============================
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tr_invtAkl_Jij 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!matrix and variables only used for basis with L=1,M=0
!=============================
INTEGER,PRIVATE::p_mk,p_ml
REAL(dp),PRIVATE::p_tau1,p_tau2,p_tau3
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:)::p_vk,p_vl,p_tvl
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_eta1,p_eta2
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_inv_Akk,p_inv_All
REAL(dp),PRIVATE,ALLOCATABLE,DIMENSION(:,:)::p_tAl_Lambda_Ak,p_tKkl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate private variables used in module matelem.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
SUBROUTINE pvars_matelem()
!====================================================
!allocate matrix used in module matelem.f90 only need 
!allocate once at the begining of the whole program 
!since they are used in the whole time
!====================================================
IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=============================
!allocate matrix and variables used for all basis form
!=============================

ALLOCATE(p_Lk(Glob_Np,Glob_Np),p_Ll(Glob_Np,Glob_Np))
ALLOCATE(p_Ak(Glob_Np,Glob_Np),p_Al(Glob_Np,Glob_Np))
ALLOCATE(p_tAl(Glob_Np,Glob_Np),p_tAkl(Glob_Np,Glob_Np),p_inv_tAkl(Glob_Np,Glob_Np))
ALLOCATE(p_dSkl_dLk(Glob_NLk),p_dHkl_dLk(Glob_NLk))
ALLOCATE(p_dTkl_dLk(Glob_NLk),p_dVkl_coulomb_dLk(Glob_NLK))

SELECT CASE(Glob_basis_form)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)

ALLOCATE(p_tr_invtAkl_Jij(Glob_Nparticle,Glob_Nparticle))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)

ALLOCATE(p_vk(Glob_Np),p_vl(Glob_Np),p_tvl(Glob_Np))
ALLOCATE(p_eta1(Glob_Nparticle,GLob_Nparticle))
ALLOCATE(p_eta2(Glob_Nparticle,GLob_Nparticle))
ALLOCATE(p_inv_Akk(Glob_Np,Glob_Np),p_inv_All(Glob_Np,Glob_Np))
ALLOCATE(p_tAl_Lambda_Ak(Glob_Np,Glob_Np),p_tKkl(Glob_Np,Glob_Np))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SELECT 

RETURN
END SUBROUTINE pvars_matelem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!matrix element for Hamilton and overlap for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HS_0_fun(bk,bl,Skl,Hkl)
!====================================================
!calculate the matrix element of Hamitonian for basis 
!with L=0,M=0
!====================================================
!input:
!  bk,bl: integer index
!output:
!  Skl: overlap of bk_th and bl_th basis
!  Hkl: Hamitonian of the bk_th and bl_th basis
!====================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::Skl,Hkl

INTEGER::i,j,k,L,ip
REAL(dp)::temp

Skl=0.0_dp
Hkl=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_0_fun(bk,bl)  

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl,Tkl,Vkl_Coulomb
CALL Skl_0_fun()
CALL Tkl_0_fun()
CALL Vkl_Coulomb_0_fun()

!Skl,Hkl
Skl=Skl+Glob_symmetry(ip)*p_Skl
Hkl=Hkl+Glob_symmetry(ip)*(p_Tkl+p_Vkl_Coulomb)

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE HS_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap check function for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_0_check_fun(bk,bl,ERR)
!====================================================
!check the overlap of bk-bl basis for basis form L=0,M=0 
!====================================================
!input:
!  bk,bl: basis index
!output:
!  ERR: 
!      0: within overlap threshold
!      1: exceed overlap threshold
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L,ip
REAL(dp)::temp

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_0_fun(bk,bl)  

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl
CALL Skl_0_fun()
IF(dabs(p_Skl)>Glob_overlap_threshold)THEN
  IF(bk/=bl)THEN
    ERR=1
    RETURN
  ENDIF
ENDIF

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_0_check_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient matrix element for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dHS_dLk_0_fun(bk,bl,dSkl_dLk,dHkl_dLk)
!==================================================
!calculate the gradient of H and S with respect to 
!Lk for basis with L=0,M=0
!==================================================
!input:
!  bk,bl: basis index
!output:
!  dSkl_dLk: gradient of overlap Skl
!  dHkl_dLk: gradient of Hamilton
!==================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::dSkl_dLk(Glob_NLk),dHkl_dLk(Glob_NLk)
  
INTEGER::i,j,k,L,ip
REAL(dp)::temp

dSkl_dLk(:)=0.0_dp
dHkl_dLk(:)=0.0_dp
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!extract parameters
CALL paras_0_fun(bk,bl)  
  
!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
  
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl,Tkl,Vkl_Coulomb
CALL Skl_0_fun()
CALL Tkl_0_fun()
CALL Vkl_Coulomb_0_fun()

!dSkl_dLk,dTkl_dLk,dVkl_Coulomb_dLk
CALL dSkl_dLk_0_fun(ip,bk,bl)
CALL dTkl_dLk_0_fun(ip,bk,bl)
CALL dVkl_Coulomb_dLk_0_fun(ip,bk,bl)
    
dSkl_dLk(:)=dSkl_dLk(:)+Glob_symmetry(ip)*p_dSkl_dLk(:)
dHkl_dLk(:)=dHkl_dLk(:)+Glob_symmetry(ip)*(p_dTkl_dLk(:)+p_dVkl_Coulomb_dLk(:))
  
ENDIF
ENDDO perm_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE  dHS_dLk_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!parameter extraction function for L=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE paras_0_fun(bk,bl)
!==================================================
!extract related parameters of the bk and bl basis 
!==================================================
!input:
!  bk,bl: basis index
!output:
!  p_Lk,p_Ll,p_Ak,p_Al
!==================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,index
REAL(dp)::temp1,temp2
REAL(dp)::vechLk(Glob_NLk),vechLl(Glob_NLk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extrat Lk and Ll from Glob_Lk
vechLk(:)=Glob_Lk(:,bk)
vechLl(:)=Glob_Lk(:,bl)

!allocate Lk and Ll
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    p_Lk(j,i)=0.0_dp
    p_Ll(j,i)=0.0_dp
    p_Lk(i,j)=vechLk(index) !lower trangle 
    p_Ll(i,j)=vechLl(index) !lower trangle
  ENDDO
ENDDO

!form Ak,Al matrix
DO j=1,Glob_Np
  DO i=j,Glob_Np
   temp1=0.0_dp
   temp2=0.0_dp
   DO k=1,j
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
END SUBROUTINE paras_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap Skl and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_0_fun()
!===================================================
!calculate the overlap matrix element for basis with L=0,M=0
!===================================================
!output: p_Skl
!===================================================
IMPLICIT NONE
INTEGER::i,j
REAL(dp)::det_Lk,det_Ll,temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(p_det_tAkl==0.0_dp)PAUSE'det_tAkl==0'

!det_Lk,det_Ll
det_Lk=1.0_dp
det_Ll=1.0_dp
DO i=1,Glob_Np
  det_Lk=det_Lk*p_Lk(i,i)
  det_Ll=det_Ll*p_Ll(i,i)
ENDDO

!Skl
temp=ABS(det_Lk*det_Ll)/p_det_tAkl
p_Skl=2.0_dp**(1.5_dp*Glob_Np)*temp*SQRT(temp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dSkl_dLk_0_fun(ip,bk,bl)
!===================================================
!calculate gradient of overlap Skl with respect to Lk
!===================================================
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
!  p_dSkl_dLk
!===================================================
INTEGER,INTENT(IN)::ip,bk,bl
  
INTEGER::i,j,k,L,index
REAL(dp)::temp,inv_Lk(Glob_Np,Glob_Np),inv_tAkl_ttAkl(Glob_Np,Glob_Np)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===========================
!bk!=bl: inv_tAkl_ttAkl=inv_tAkl
!bk==bl: inv_tAkl_ttAkl=inv_tAkl+inv_ttAkl
!===========================
  
CALL inv_L_lower(Glob_Np,p_Lk,inv_Lk)
  
DO j=1,Glob_Np
  DO i=1,Glob_Np
    inv_tAkl_ttAkl(i,j)=p_inv_tAkl(i,j)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN
      
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*p_inv_tAkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    inv_tAkl_ttAkl(i,j)=inv_tAkl_ttAkl(i,j)+temp
    inv_tAkl_ttAkl(j,i)=inv_tAkl_ttAkl(i,j)
  ENDDO
ENDDO 

DO j=1,Glob_Np
  DO i=1,Glob_Np  
    inv_Lk(i,j)=2.0_dp*inv_Lk(i,j)
  ENDDO
ENDDO
  
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!dSkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
      temp=temp+inv_tAkl_ttAkl(i,k)*p_Lk(k,j)
    ENDDO

    p_dSkl_dLk(index)=1.5_dp*p_Skl*(inv_Lk(j,i)-2.0_dp*temp)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dSkl_dLk_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kinetic energy and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tkl_0_fun()
!===================================================
!calculate the overlap of kinetic energy for basis with L=0,M=0
!===================================================
!output: p_Tkl
!===================================================
IMPLICIT NONE
INTEGER::i,j,k,L
REAL(dp)::temp,Atemp(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!tAl*inv_tAkl*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_tAl(i,k)*p_inv_tAkl(k,L)*p_Ak(L,j)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO

!Tkl=6*S_kl*tr[Lambda*tAl*inv_tAkl*Ak]
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+Glob_Lambda(i,k)*Atemp(k,i)
  ENDDO
ENDDO
p_Tkl=6.0_dp*p_Skl*temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Tkl_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dTkl_dLk_0_fun(ip,bk,bl)
!===================================================
!calculate the gradient of kinetic energy
!===================================================
!input:
!    ip: permutation index
! bk,bl: basis index
!output:
!  p_dTkl_dLk 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip,bk,bl
  
INTEGER::i,j,k,L,index
REAL(dp)::temp,Atemp(Glob_Np,Glob_Np),Btemp(Glob_Np,Glob_Np),F_G(Glob_Np,Glob_Np)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!===============================
!bk!=bl: F_G=F
!bk==bl: F_G=F+G
!===============================
  
!inv_tAkl*tAl*Lambda
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,GLob_Np
       temp=temp+p_inv_tAkl(i,k)*p_tAl(k,L)*Glob_Lambda(L,j)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO

!F=inv_tAkl*tAl*Lambda*tAl*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,GLob_Np
        temp=temp+Atemp(i,k)*p_tAl(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    F_G(i,j)=temp
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN

!inv_tAkl*Ak*Lambda
DO j=1,GLob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*p_Ak(k,L)*Glob_Lambda(L,j)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO

!inv_tAkl*Ak*Lambda*Ak*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Atemp(i,k)*p_Ak(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    Btemp(i,j)=temp
  ENDDO
ENDDO

!G=Tp(inv_tAkl*Ak*Lambda*Ak*inv_tAkl)Tp'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*Btemp(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    F_G(i,j)=F_G(i,j)+temp
  ENDDO
ENDDO
  
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dTkl_dLl
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
      temp=temp+F_G(i,k)*p_Lk(k,j)
    ENDDO

    p_dTkl_dLk(index)=p_Tkl*p_dSkl_dLk(index)/p_Skl+12.0_dp*p_Skl*temp
  ENDDO
ENDDO
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dTkl_dLk_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!coulumb potential and its gradient for basis with L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Vkl_coulomb_0_fun()
!===================================================
!calculate Coulomb potential energy for basis with L=0,M=0
!===================================================
!output:
!  p_tr_inv_tCkl_Jij,p_Vkl_Coulomb
!===================================================
IMPLICIT NONE
INTEGER::i,j,k,ii,jj
REAL(dp)::temp,Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!interactive Coulomb potential between each particle: Vkl
!============================

Vkl=0.0_dp

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!Tr[inv_tAkl*Jij]
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO           
p_tr_invtAkl_Jij(ii,jj)=temp   !store Tr[inv_tAkl*Jij] 

Vkl=Vkl+Glob_charge(ii)*Glob_charge(jj)/dsqrt(temp)
   
ENDDO loop2
ENDDO loop1

p_Vkl_Coulomb=(2.0_dp/SQRTPI)*p_Skl*Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Vkl_Coulomb_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dVkl_Coulomb_dLk_0_fun(ip,bk,bl)
!===================================================
!calculate the gradient of Coulomb potential energy
!with respect to Lk
!===================================================
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
! p_dVkl_Coulomb_dLk
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,ii,jj,index
REAL(dp)::temp,Atemp(Glob_Np,Glob_Np)
REAL(dp)::tQ_ttQ(Glob_Np,Glob_Np),dVkl_dLk(Glob_NLk)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================
!interactive Coulomb potential between each particle: dVkl_dLk
!=====================

dVkl_dLk(:)=0.0_dp

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle
    
!tQ=inv_tAkl*Jij*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    tQ_ttQ(i,j)=temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN
    
!ttQ=P*tQ*P'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tQ_ttQ(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO
  
!tQ_ttQ=tQ+ttQ
DO j=1,Glob_Np
  DO i=1,Glob_Np
    tQ_ttQ(i,j)=tQ_ttQ(i,j)+Atemp(i,j)
  ENDDO
ENDDO
  
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!dVkl_dLK=sum{qi*qj*tr[inv_tAkl*Jij]^(-3/2)*vech(tQ_ttQ*Lk)}
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
     temp=temp+tQ_ttQ(i,k)*p_Lk(k,j)
    ENDDO

    dVkl_dLk(index)=dVkl_dLk(index)+Glob_charge(ii)*Glob_charge(jj)&
    &*temp*p_tr_invtAkl_Jij(ii,jj)**(-1.5_dp)
  ENDDO
ENDDO

ENDDO loop2
ENDDO loop1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

DO i=1,GLob_NLk
  p_dVkl_Coulomb_dLk(i)=p_Vkl_Coulomb*p_dSkl_dLk(i)/p_Skl&
  &+(2.0_dp/SQRTPI)*p_Skl*dVkl_dLk(i)
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dVkl_Coulomb_dLk_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function Gij for L=0,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE  gij_kl_0_fun(bk,bl,ii,jj,R,gij)
!===================================================
!this subroutine calculate the pair correlation between 
!bk and bl basis gij(bk,bl)=<bk|delta(r-x)|bl>
!===================================================
!input: 
! bk,bl: basis index
! ii,jj: particle index
!     R: radius
!output:
!   gij: the pair corelation 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl,ii,jj
REAL(dp),INTENT(IN)::R
REAL(dp),INTENT(OUT)::gij
  
INTEGER::i,j,k,L,ip
REAL(dp)::temp,temp1,temp2

gij=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters
CALL paras_0_fun(bk,bl)  
      
!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl,Tkl,Vkl
CALL Skl_0_fun()
  
!tr[in_tAkl*Jij]
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO  
  
!tr[in_tAkl*Jij]^(-3/2)
temp1=temp**(-1.5_dp)

!exp(-x^2/tr[inv_tAkl*Jij])
temp2=EXP(-R**2/temp)
    
!Skl
gij=gij+Glob_symmetry(ip)*p_Skl*temp1*temp2/(PI*SQRTPI)

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gij_kl_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!expectation value of interparticle distance for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_kl_0_fun(bk,bl,ii,jj,v,rij)
!===================================================
!this subroutine calculate the expectation value of 
!interparticle distance for basis with L=0,M=0
!<rij^v>=<bk|rij^v|bl>
!===================================================
!input: 
! bk,bl: basis index
! ii,jj: particle index
!     v: power of rij
!output:
!   rij: expectation value of interparticle diatance 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl,ii,jj
REAL(dp),INTENT(IN)::v
REAL(dp),INTENT(OUT)::rij

INTEGER::i,j,k,L,ip
REAL(dp)::temp

rij=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters
CALL paras_0_fun(bk,bl)  
    
!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)

!Skl
CALL Skl_0_fun()

![tr[in_tAkl*Jij]]^(v/2)
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
temp=temp**(v*0.5_dp)   

!<rij>
rij=rij+Glob_symmetry(ip)*p_Skl*(2.0_dp/SQRTPI)*GAMMA(0.5_dp*v+1.5_dp)*temp

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE rij_kl_0_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Hamilton and overlap for L=1,M=0 basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HS_1_fun(bk,bl,Skl,Hkl)
!===================================================
!Hamiton and overlap matrix element for basis with L=1,M=0
!===================================================
!input:
!  bk,bl: basis index
!output:
!  Skl: overlap
!  Hkl: Hamilton 
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::Skl,Hkl

INTEGER::i,j,k,L,ip
REAL(dp)::temp

Skl=0.0_dp
Hkl=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_1_fun(bk,bl)

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO

!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl,Tkl,Vkl
CALL Skl_1_fun()
CALL Tkl_1_fun()
CALL Vkl_coulomb_1_fun()

!Skl,Hkl
Skl=Skl+Glob_symmetry(ip)*p_Skl
Hkl=Hkl+Glob_symmetry(ip)*(p_Tkl+p_Vkl_coulomb)

ENDIF
ENDDO perm_loop
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE HS_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap threshold check fun for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_1_check_fun(bk,bl,ERR)
!===================================================
!ckeck the overlap of bk-bl basis for basis with L=1,M=0 
!===================================================
!input:
!bk,bl: basis index
!output:
!  ERR: 
!      0: within overlap threshold
!      1: exceed overlap threshold
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
INTEGER,INTENT(OUT)::ERR

INTEGER::i,j,k,L,ip
REAL(dp)::Skl,temp

ERR=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_1_fun(bk,bl)

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
  
!det_tAkl,inv_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl
CALL Skl_1_fun()

IF(dabs(p_Skl)>Glob_overlap_threshold)THEN
  IF(bk/=bl)THEN
    ERR=1
    RETURN
  ENDIF
ENDIF

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_1_check_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gradient of Hamiton and overlap for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dHS_dLk_1_fun(bk,bl,dSkl_dLk,dHkl_dLk)
!===================================================
!calculate the gradient of Hamilton and overlap 
!for basis with L=1,M=0
!===================================================
!input:
! bk,bl: basis index
!output:
! dSkl_dLk: gradient of overlap
! dHkl_dLk: gradient of Hamiton
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl
REAL(dp),INTENT(OUT)::dSkl_dLk(Glob_NLk),dHkl_dLk(Glob_NLk)
    
INTEGER::i,j,k,L,ip
REAL(dp)::temp
  
dSkl_dLk(:)=0.0_dp
dHkl_dLk(:)=0.0_dp
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
!extract parameters
CALL paras_1_fun(bk,bl)

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
    
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
  
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
    
!Skl,Tkl,Vkl
CALL Skl_1_fun()
CALL Tkl_1_fun()
CALL Vkl_coulomb_1_fun()
  
!dSkl_dLk,dTkl_dLk,dVkl_Coulomb_dLk
CALL dSkl_dLk_1_fun(ip,bk,bl)
CALL dTkl_dLk_1_fun(ip,bk,bl)
CALL dVkl_coulomb_dLk_1_fun(ip,bk,bl)

dSkl_dLk(:)=dSkl_dLk(:)+Glob_symmetry(ip)*p_dSkl_dLk(:)
dHkl_dLk(:)=dHkl_dLk(:)+Glob_symmetry(ip)*(p_dTkl_dLk(:)+p_dVkl_coulomb_dLk(:))

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dHS_dLk_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!parameters extraction function for L=1,M=0 basis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE paras_1_fun(bk,bl)
!===================================================
!extract related parameter matrix of the 
!k_th and l_th basis with L=1,M=0
!===================================================
!input:
!  bk,bl: basis index
!output:
!  p_mk,p_ml,p_vk,p_vl,p_Lk,p_Ll,p_Ak,p_Al
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,index
REAL(dp)::temp1,temp2,vechLk(Glob_NLk),vechLl(Glob_NLk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract mk and ml from GLob_mk
p_mk=Glob_mk(bk)
p_ml=Glob_mk(bl)

p_vk=0.0_dp
p_vl=0.0_dp
p_vk(p_mk)=1.0_dp
p_vl(p_ml)=1.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract Lk and Ll from Glob_Lk
vechLk(:)=Glob_Lk(1:GLob_NLk,bk)
vechLl(:)=Glob_Lk(1:GLob_NLk,bl)

!allocate Lk and Ll
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    p_Lk(j,i)=0.0_dp
    p_Ll(j,i)=0.0_dp
    p_Lk(i,j)=vechLk(index) !lower trangle 
    p_Ll(i,j)=vechLl(index) !lower trangle
  ENDDO
ENDDO

!form Ak,Al matrix
DO j=1,Glob_Np
  DO i=j,Glob_Np
   temp1=0.0_dp
   temp2=0.0_dp
   DO k=1,j
    temp1=temp1+p_Lk(i,k)*p_Lk(j,k)
    temp2=temp2+p_Ll(i,k)*p_Ll(j,k)
   ENDDO
   p_Ak(i,j)=temp1
   p_Ak(j,i)=temp1
   p_Al(i,j)=temp2
   p_Al(j,i)=temp2
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE paras_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!overlap S_kl and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Skl_1_fun()
!===================================================
!overlap element for basis with L=1,M=0
!===================================================
!output:
! p_inv_Akk,p_inv_All,p_tau3,p_Skl
!===================================================
IMPLICIT NONE

INTEGER::i,j,k,L
REAL(dp)::temp,temp1,temp2,det_Lk,det_Ll
REAL(dp)::inv_Lk(Glob_Np,Glob_Np),inv_Ll(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!det_Lk,det_Ll
det_Lk=1.0_dp
det_Ll=1.0_dp
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
    temp1=0.0_dp
    temp2=0.0_dp
    DO k=1,Glob_Np
      temp1=temp1+inv_Lk(k,i)*inv_Lk(k,j)
      temp2=temp2+inv_Ll(k,i)*inv_Ll(k,j)
    ENDDO
    temp1=0.5_dp*temp1
    temp2=0.5_dp*temp2
    p_inv_Akk(i,j)=temp1
    p_inv_Akk(j,i)=temp1
    p_inv_All(i,j)=temp2
    p_inv_All(j,i)=temp2
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!vk'*inv_Akk*vk
temp1=p_inv_Akk(p_mk,p_mk)

!vl'*inv_All*vl
temp2=p_inv_All(p_ml,p_ml)

!tau3=tr[tvl*vk'*inv_tAkl]=vk'*inv_tAkl*tvl
p_tau3=0.0_dp
DO i=1,Glob_Np
  p_tau3=p_tau3+p_inv_tAkl(p_mk,i)*p_tvl(i)
ENDDO

temp=(ABS(det_Lk*det_Ll)/p_det_tAkl)**(1.5_dp)
temp=temp*2.0_dp**(1.5_dp*Glob_Np)

!Skl
p_Skl=temp*(p_tau3/SQRT(temp1*temp2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Skl_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dSkl_dLk_1_fun(ip,bk,bl)
!===================================================
!gradient of overlap Skl for basis with L=1,M=0
!===================================================
!input:
!     ip: permutation index
!  bk,bl: basis index
!output:
! p_dSkl_dBk
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
REAL(dp)::temp,Atemp(Glob_Np,Glob_Np)
REAL(dp)::Fkk_Fll(Glob_Np,Glob_Np),tFkl_tGkl(Glob_Np,Glob_Np)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Fkk
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=p_inv_Akk(i,p_mk)*p_inv_Akk(p_mk,j)/p_inv_Akk(p_mk,p_mk)&
    &+1.5_dp*p_inv_Akk(i,j)
    Fkk_Fll(i,j)=temp
    Fkk_Fll(j,i)=temp
  ENDDO
ENDDO

!tFkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    tFkl_tGkl(i,j)=p_tKkl(i,j)/p_tau3+1.5_dp*p_inv_tAkl(i,j)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN

!Fkk+Fll
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=p_inv_All(i,p_ml)*p_inv_All(p_ml,j)/p_inv_All(p_ml,p_ml)&
    &+1.5_dp*p_inv_All(i,j)
    Fkk_Fll(i,j)=Fkk_Fll(i,j)+temp
    Fkk_Fll(j,i)=Fkk_Fll(i,j)
  ENDDO
ENDDO
  
!tGkl=Tp*tFkl*Tp'
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tFkl_tGkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO

!tFkl+tGkl
DO j=1,Glob_Np
  DO i=1,Glob_Np    
    tFkl_tGkl(i,j)=tFkl_tGkl(i,j)+Atemp(i,j)
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dSkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
      temp=temp+(2.0_dp*Fkk_Fll(i,k)-tFkl_tGkl(i,k)-tFkl_tGkl(k,i))*p_Lk(k,j)
    ENDDO

    p_dSkl_dLk(index)=p_Skl*temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
RETURN
END SUBROUTINE dSkl_dLk_1_fun  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!kinetic energy and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Tkl_1_fun()
!===================================================
!matrix elment of kinetic energy for basis with L=1,M=0
!===================================================
!output:
!  p_tKkl,p_tAl_Lambda_Ak,p_tau1,p_tau2,p_Tkl
!===================================================
IMPLICIT NONE
INTEGER::i,j,k,L
REAL(dp)::temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!tAl*Lambda*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_tAl(i,k)*Glob_Lambda(k,L)*p_Ak(L,j)
      ENDDO
    ENDDO
    p_tAl_Lambda_Ak(i,j)=temp
  ENDDO
ENDDO

!tau1 and tau2
p_tau1=0.0_dp
p_tau2=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    p_tau1=p_tau1+p_inv_tAkl(i,k)*p_tAl_Lambda_Ak(k,i)
    p_tau2=p_tau2+p_tKkl(i,k)*p_tAl_Lambda_Ak(k,i)
  ENDDO
ENDDO

!Tkl
p_Tkl=p_Skl*(6.0_dp*p_tau1+4.0_dp*p_tau2/p_tau3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Tkl_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dTkl_dLk_1_fun(ip,bk,bl)
!===================================================
!gradient of kinetic energy for basis with L=1,M=0
!===================================================
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
! p_Tkl_dLk
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip
INTEGER,INTENT(IN)::bk,bl

INTEGER::i,j,k,L,index
REAL(dp)::temp,temp1,temp2,Atemp(Glob_Np,Glob_Np)
REAL(dp)::tAl_Lambda_tAl(Glob_Np,Glob_Np),Ak_Lambda_Ak(Glob_Np,Glob_Np)
REAL(dp)::tUkl_tWkl(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!tAl*Lambda*tAl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
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
    temp=0.0_dp
    temp1=0.0_dp
    temp2=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*tAl_Lambda_tAl(k,L)*p_inv_tAkl(L,j)
        temp1=temp1+p_tKkl(i,k)*tAl_Lambda_tAl(k,L)*p_inv_tAkl(L,j)
        temp2=temp2+p_inv_tAkl(i,k)*p_tAl_Lambda_Ak(k,L)*p_tKkl(L,j)
      ENDDO
    ENDDO
    tUkl_tWkl(i,j)=6.0_dp*temp+(4.0_dp/p_tau3)*(temp1-temp2)&
    &+4.0_dp*p_tau2*p_tKkl(i,j)/(p_tau3**2)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN

!Ak*Lambda*Ak
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
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
    temp=0.0_dp
    temp1=0.0_dp
    temp2=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Ak_Lambda_Ak(k,L)*p_inv_tAkl(L,j)
        temp1=temp1+p_inv_tAkl(i,k)*Ak_Lambda_Ak(k,L)*p_tKkl(L,j)
        temp2=temp2+p_tKkl(i,k)*p_tAl_Lambda_Ak(k,L)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    Atemp(i,j)=6.0_dp*temp+(4.0_dp/p_tau3)*(temp1-temp2)&
    &+4.0_dp*p_tau2*p_tKkl(i,j)/(p_tau3**2)
  ENDDO
ENDDO

!tUkl+tWkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*Atemp(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    tUkl_tWkl(i,j)=tUkl_tWkl(i,j)+temp
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dTkl_dLk
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
      temp=temp+(tUkl_tWkl(i,k)+tUkl_tWkl(k,i))*p_Lk(k,j)
    ENDDO

    p_dTkl_dLk(index)=p_Tkl*p_dSkl_dLk(index)/p_Skl+p_Skl*temp
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dTkl_dLk_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!coulumb potential energy and its gradient for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE Vkl_coulomb_1_fun()
!===================================================
!Coulomb potential energy Vkl_Coulomb for basis with L=1,M=0
!===================================================
!output:
!  p_eta1,p_eta2,p_Vkl_Coulomb
!===================================================
IMPLICIT NONE
INTEGER::i,j,k,L,ii,jj
REAL(dp)::temp,Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!interactive Coulomb potential between each particle: Vkl
!============================

Vkl=0.0_dp

loop1:DO ii=1,GLob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!eta1
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
p_eta1(ii,jj)=temp

!eta2
temp=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    temp=temp+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
p_eta2(ii,jj)=temp

temp=p_eta1(ii,jj)**(-0.5_dp)
temp=temp*(1.0_dp-p_eta2(ii,jj)/p_eta1(ii,jj)/(p_tau3*3.0_dp))

Vkl=Vkl+Glob_charge(ii)*Glob_charge(jj)*temp

ENDDO loop2
ENDDO loop1

p_Vkl_coulomb=(2.0_dp/SQRTPI)*p_Skl*Vkl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE Vkl_coulomb_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE dVkl_coulomb_dLk_1_fun(ip,bk,bl)
!===================================================
!the gradient of Coulomb potential energy for basis with L=1,M=0
!===================================================
!input:
!   ip: permutation index
!bk,bl: basis index
!output:
!  p_dVkl_Coulomb_dLk
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::ip,bk,bl

INTEGER::i,j,k,L,ii,jj,index
REAL(dp)::temp,temp1,temp2,Atemp(Glob_Np,Glob_Np)
REAL(dp)::dVkl_dLk(Glob_NLk),tQkl_tDkl(Glob_Np,Glob_Np)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!============================
!gradient of interactive Coulomb potential between each particle: dVkl_dLk
!============================

dVkl_dLk(:)=0.0_dp

loop1:DO ii=1,Glob_Nparticle-1
loop2:DO jj=ii+1,Glob_Nparticle

!tQkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    temp1=0.0_dp
    temp2=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
        temp1=temp1+p_inv_tAkl(i,k)*Glob_Jij(k,L,ii,jj)*p_tKkl(L,j)
        temp2=temp2+p_tKkl(i,k)*Glob_Jij(k,L,ii,jj)*p_inv_tAkl(L,j)
      ENDDO
    ENDDO
    temp=0.5_dp*(1.0_dp-p_eta2(ii,jj)/p_eta1(ii,jj)/p_tau3)*temp
    temp=temp+(temp1+temp2)/p_tau3/3.0_dp
    temp=temp-p_eta2(ii,jj)*p_tKkl(i,j)/(p_tau3**2)/3.0_dp
    tQkl_tDkl(i,j)=temp*p_eta1(ii,jj)**(-1.5_dp)
  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(bk==bl)THEN

!tDkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      DO L=1,Glob_Np
        temp=temp+Glob_Tp(i,k,ip)*tQkl_tDkl(k,L)*Glob_Tp(j,L,ip)
      ENDDO
    ENDDO
    Atemp(i,j)=temp
  ENDDO
ENDDO

!tQkl+tDkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    tQkl_tDkl(i,j)=tQkl_tDkl(i,j)+Atemp(i,j)
  ENDDO
ENDDO

ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!dVkl_1_dLK
index=0
DO j=1,Glob_Np
  DO i=j,Glob_Np
    index=index+1
    
    temp=0.0_dp
    DO k=j,Glob_Np
      temp=temp+(tQkl_tDkl(i,k)+tQkl_tDkl(k,i))*p_Lk(k,j)
    ENDDO

    dVkl_dLk(index)=dVkl_dLk(index)+Glob_charge(ii)*Glob_charge(jj)*temp
  ENDDO
ENDDO

ENDDO loop2
ENDDO loop1

DO i=1,Glob_NLk
  p_dVkl_coulomb_dLk(i)=p_Vkl_coulomb*p_dSkl_dLk(i)/p_Skl&
  &+(2.0_dp/SQRTPI)*p_Skl*dVkl_dLk(i)
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE dVkl_coulomb_dLk_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!correlation function for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gij_kl_1_fun(bk,bl,ii,jj,Rxy,Rz,gij)
!===================================================
!this subroutine calculates the correaltion function
!for basis with L=1,M=0
!===================================================
!input: 
! bk,bl: basis index
! ii,jj: particle index
!   Rxy: sqrt(x**2+y**2)
!    Rz: z component of R
!output:
!   gij: correlation function coresbonding to (R,z)
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl,ii,jj
REAL(dp),INTENT(IN)::Rxy,Rz
REAL(dp),INTENT(OUT)::gij

INTEGER::i,j,k,L,ip
REAL(dp)::eta1,eta2,temp,temp1,temp2,R

gij=0.0_dp
R=SQRT(Rxy**2+Rz**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_1_fun(bk,bl)

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl and tau3
CALL Skl_1_fun()

!eta1
eta1=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta1=eta1+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!eta2
eta2=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta2=eta2+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO
  
temp1=1.0_dp+(eta2/eta1/p_tau3)*(2.0_dp*(Rz**2)/eta1-1.0_dp)
temp2=EXP(-R**2/eta1)*(PI*eta1)**(-1.5_dp)

gij=gij+Glob_symmetry(ip)*p_Skl*temp1*temp2

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE gij_kl_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!expectation value of interparticle distace for basis with L=1,M=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rij_kl_1_fun(bk,bl,ii,jj,v,rij)
!===================================================
!this subroutine calculate the expectation value of 
!interparticle diatance for basis with L=1,M=0
!===================================================
!input:
!  bk,bl: basis index
!  ii,jj: particle index
!      v: power of distance
!output:
!    rij: expectation value of interparticle diatance
!===================================================
IMPLICIT NONE
INTEGER,INTENT(IN)::bk,bl,ii,jj
REAL(dp),INTENT(IN)::v
REAL(dp),INTENT(OUT)::rij

INTEGER::i,j,k,L,ip
REAL(dp)::eta1,eta2,temp

rij=0.0_dp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!extract parameters
CALL paras_1_fun(bk,bl)

!permutation loop
perm_loop:DO ip=1,Glob_Nperm
IF(Glob_symmetry(ip)/=0.0_dp)THEN
        
!tAl=Tp'*Al*Tp
DO j=1,Glob_Np
  DO i=j,Glob_Np
    temp=0.0_dp
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
DO j=1,Glob_Np
  DO i=1,Glob_Np
    p_tAkl(i,j)=p_Ak(i,j)+p_tAl(i,j)
  ENDDO
ENDDO
    
!inv_tAkl,det_tAkl
CALL det_inv_fun(Glob_Np,p_tAkl,p_det_tAkl,p_inv_tAkl)
  
!Skl and tau3
CALL Skl_1_fun()

!eta1
eta1=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta1=eta1+p_inv_tAkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

!tKkl=inv_tAkl*tvl*vk'*inv_tAkl
DO j=1,Glob_Np
  DO i=1,Glob_Np
    temp=0.0_dp
    DO k=1,Glob_Np
      temp=temp+p_inv_tAkl(i,k)*p_tvl(k)*p_inv_tAkl(p_mk,j)
    ENDDO
    p_tKkl(i,j)=temp
  ENDDO
ENDDO

!eta2
eta2=0.0_dp
DO i=1,Glob_Np
  DO k=1,Glob_Np
    eta2=eta2+p_tKkl(i,k)*Glob_Jij(k,i,ii,jj)
  ENDDO
ENDDO

temp=1.0_dp+(v/3.0_dp)*(eta2/eta1/p_tau3)
temp=temp*GAMMA(0.5_dp*v+1.5_dp)*eta1**(0.5_dp*v)

rij=rij+Glob_symmetry(ip)*(2.0_dp/SQRTPI)*p_Skl*temp

ENDIF
ENDDO perm_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

RETURN
END SUBROUTINE rij_kl_1_fun

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE matelem
