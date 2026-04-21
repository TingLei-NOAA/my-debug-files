!#
!                                            *********************************
!                                            *      module ybfil             *
!                                            *      R. J. Purser             *
!                                            *      NOAA/NCEP/EMC            *
!                                            *      April 2026               *
!                                            *    jim.purser@noaa.gov        *
!                                            *********************************
!
! Codes for the beta radial filters.
! The filters invoke the aspect tensor information encoded by the 
! Cholesky lower-triangular factors, el, of the INVERSE aspect tensors.
! The routines, "cholaspect", convert (in place) the field of given
! aspect tensors A to the equivalent cholesky factors of A^(-1).
! The routines, "getlinesum" precompute the normalization coefficients for
! each line (row) of the implied matrix form of the beta filter so that the
! normalized line sum associated with each point of application becomes
! unity. This makes the application of each filter significantly faster
! than having to work out the normalization on the fly.
!
! COMPILE AFTER: {pmat}
!
!=============================================================================
module ybfil
!=============================================================================
use pkind, only: spi,dp
use pietc, only: u0,u1
implicit none
private
public:: &
     t22_to_3,t2_to_3,t3_to_22,t33_to_6,t3_to_6,t6_to_33,&
     t44_to_10,t4_to_10,t10_to_44,p,rsg,i2pair,i3pair,i4pair, &
     finp,inip,&
     rcalib,rbeta,rbetat,flip,flipt
integer(spi),parameter       :: nsres=50
real(dp),    parameter       :: deltabi=10! inverse of sres tables' delta-b
real(dp),parameter           :: eps=1.e-12_dp,emin=.3e0_dp,&
     om0=1.170561,om1=-4.345418,om2=2.583651 ! (om0,om1,om2 all computed in wasserstein2)
real(dp),dimension(0:nsres)  :: sres2,sres3
real(dp)                     :: rsg,rpp3o2,rpp4o2,rpp5o2,rpp6o2
integer(spi)                 :: p
integer(spi),dimension(2,0:2):: i2pair
integer(spi),dimension(2,6)  :: i3pair
integer(spi),dimension(2,10) :: i4pair

! p=2 s-residual for each tabulated scale b, s being the half-span:
data sres2/0.9999999E+00,& ! 1-epsilon (avoids null line filtering)
0.8390E+00,0.6805E+00,0.5256E+00,0.3757E+00,0.2328E+00,&
0.1004E+00,-.1556E-01,-.1031E+00,-.1345E+00,-.2307E-01,&
0.9016E-01,0.3778E-01,-.1302E-01,-.4396E-01,-.3795E-01,&
0.2998E-01,0.3412E-01,0.6966E-02,-.1677E-01,-.2450E-01,&
-.5024E-02,0.2284E-01,0.1201E-01,-.4674E-02,-.1482E-01,&
-.1056E-01,0.1156E-01,0.1199E-01,0.1285E-02,-.8219E-02,&
-.9839E-02,0.1731E-02,0.9906E-02,0.4143E-02,-.3733E-02,&
-.7635E-02,-.3235E-02,0.6907E-02,0.5222E-02,-.7235E-03,&
-.5267E-02,-.4615E-02,0.3483E-02,0.5193E-02,0.1212E-02,&
-.3152E-02,-.4450E-02,-.1756E-04,0.4429E-02,0.2340E-02/
! p=3 s-residual for each tabulated scale b, s being the half-span:
data sres3/0.9999999E+00,& ! 1-epsilon (avoids null line filtering)
0.8637E+00,0.7058E+00,0.5476E+00,0.3934E+00,0.2467E+00,&
0.1123E+00,-.1744E-02,-.8045E-01,-.8895E-01,0.3107E-01,&
0.3997E-01,0.8750E-02,-.1796E-01,-.2029E-01,0.7223E-02,&
0.1323E-01,0.2605E-02,-.7944E-02,-.6694E-02,0.4073E-02,&
0.5567E-02,0.2230E-03,-.4496E-02,-.2143E-02,0.2710E-02,&
0.2604E-02,-.6152E-03,-.2734E-02,-.4276E-03,0.1852E-02,&
0.1211E-02,-.8620E-03,-.1619E-02,0.2304E-03,0.1253E-02,&
0.4811E-03,-.8617E-03,-.8482E-03,0.4500E-03,0.8209E-03,&
0.8178E-04,-.7499E-03,-.3567E-03,0.4787E-03,0.5061E-03,&
-.1339E-03,-.5853E-03,-.6824E-04,0.4252E-03,0.2775E-03/

data p/2/ ! default
data i2pair/1,1, 2,2, 2,1/
data i3pair/1,1, 2,2, 3,3, 3,2, 3,1, 2,1/
!data i4pair/1,1, 2,2, 3,3, 4,4, 2,1, 3,1, 4,1, 4,3, 4,2, 3,2/
data i4pair/1,1, 2,2, 3,3, 4,4, 2,1, 3,1, 4,1, 3,2, 4,2, 4,3/
interface t22_to_3;    module procedure i22_to_3, r22_to_3;   end interface
interface t2_to_3;     module procedure i2_to_3,  r2_to_3;    end interface
interface t3_to_22;    module procedure i3_to_22, r3_to_22;   end interface
interface t33_to_6;    module procedure i33_to_6, r33_to_6;   end interface
interface t3_to_6;     module procedure i3_to_6,  r3_to_6;    end interface
interface t6_to_33;    module procedure i6_to_33, r6_to_33;   end interface
interface t44_to_10;   module procedure i44_to_10,r44_to_10;  end interface
interface t4_to_10;    module procedure i4_to_10, r4_to_10;   end interface
interface t10_to_44;   module procedure i10_to_44,r10_to_44;  end interface
!---
interface finp;        module procedure finp;                 end interface
interface inip;        module procedure inip;                 end interface
interface rcalib
   module procedure rcalib1,rcalib2,rcalib3,rcalib4
end interface rcalib
interface rbetat
   module procedure rbeta1t,vrbeta1t, rbeta2t,vrbeta2t, rbeta3t,vrbeta3t, rbeta4t,vrbeta4t
end interface
interface flipt
   module procedure flip1t,vflip1t,   flip2t,vflip2t,   flip3t,vflip3t,   flip4t,vflip4t
end interface
interface flip
   module procedure flip1,vflip1,     flip2,vflip2,     flip3,vflip3,     flip4,vflip4
end interface
interface rbeta
   module procedure rbeta1,vrbeta1,   rbeta2,vrbeta2,   rbeta3,vrbeta3,   rbeta4,vrbeta4
end interface
contains

!==============================================================================
subroutine i22_to_3(i22,i3)!                                         [t22_to_3]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(2,2),intent(in ):: i22
integer(spi),dimension(0:2),intent(out):: i3
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2; i3(L)=i22(i2pair(1,L),i2pair(2,L)); enddo
end subroutine i22_to_3
!==============================================================================
subroutine r22_to_3(r22,r3)!                                         [t22_to_3]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(2,2),intent(in ):: r22
real(dp),dimension(0:2),intent(out):: r3
!-----------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2; r3(L)=r22(i2pair(1,L),i2pair(2,L)); enddo
end subroutine r22_to_3

!==============================================================================
subroutine i2_to_3(i2,i3)!                                            [t2_to_3]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(2),intent(in ):: i2
integer(spi),dimension(3),intent(out):: i3
!------------------------------------------------------------------------------
call t22_to_3(outer_product(i2,i2),i3)
end subroutine i2_to_3
!==============================================================================
subroutine r2_to_3(r2,r3)!                                            [t2_to_3]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(2),intent(in ):: r2
real(dp),dimension(3),intent(out):: r3
!------------------------------------------------------------------------------
call t22_to_3(outer_product(r2,r2),r3)
end subroutine r2_to_3

!==============================================================================
subroutine i3_to_22(i3,i22)!                                         [t3_to_22]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(0:2),intent(in ):: i3
integer(spi),dimension(2,2),intent(out):: i22
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2
   i22(i2pair(1,L),i2pair(2,L))=i3(L)
   i22(i2pair(2,L),i2pair(1,L))=i3(L)
enddo
end subroutine i3_to_22
!==============================================================================
subroutine r3_to_22(r3,r22)!                                         [t3_to_22]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(0:2),intent(in ):: r3
real(dp),dimension(2,2),intent(out):: r22
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=0,2
   r22(i2pair(1,L),i2pair(2,L))=r3(L)
   r22(i2pair(2,L),i2pair(1,L))=r3(L)
enddo
end subroutine r3_to_22

!==============================================================================
subroutine i33_to_6(i33,i6)!                                         [t33_to_6]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(3,3),intent(in ):: i33
integer(spi),dimension(6)  ,intent(out):: i6
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6; i6(L)=i33(i3pair(1,L),i3pair(2,L)); enddo
end subroutine i33_to_6
!==============================================================================
subroutine r33_to_6(r33,r6)!                                         [t33_to_6]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(3,3),intent(in ):: r33
real(dp),dimension(6)  ,intent(out):: r6
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6; r6(L)=r33(i3pair(1,L),i3pair(2,L)); enddo
end subroutine r33_to_6

!==============================================================================
subroutine i3_to_6(i3,i6)!                                            [t3_to_6]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(3),intent(in ):: i3
integer(spi),dimension(6),intent(out):: i6
!------------------------------------------------------------------------------
call t33_to_6(outer_product(i3,i3),i6)
end subroutine i3_to_6
!==============================================================================
subroutine r3_to_6(r3,r6)!                                            [t3_to_6]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(3),intent(in ):: r3
real(dp),dimension(6),intent(out):: r6
!------------------------------------------------------------------------------
call t33_to_6(outer_product(r3,r3),r6)
end subroutine r3_to_6

!==============================================================================
subroutine i6_to_33(i6,i33)!                                         [t6_to_33]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(6),  intent(in ):: i6
integer(spi),dimension(3,3),intent(out):: i33
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6
   i33(i3pair(1,L),i3pair(2,L))=i6(L)
   i33(i3pair(2,L),i3pair(1,L))=i6(L)
enddo
end subroutine i6_to_33
!==============================================================================
subroutine r6_to_33(r6,r33)!                                         [t6_to_33]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(6),  intent(in ):: r6
real(dp),dimension(3,3),intent(out):: r33
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,6
   r33(i3pair(1,L),i3pair(2,L))=r6(L)
   r33(i3pair(2,L),i3pair(1,L))=r6(L)
enddo
end subroutine r6_to_33

!==============================================================================
subroutine i44_to_10(i44,i10)!                                      [t44_to_10]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(4,4),intent(in ):: i44
integer(spi),dimension(10) ,intent(out):: i10
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10; i10(L)=i44(i4pair(1,L),i4pair(2,L)); enddo
end subroutine i44_to_10
!==============================================================================
subroutine r44_to_10(r44,r10)!                                      [t44_to_10]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(4,4),intent(in ):: r44
real(dp),dimension(10) ,intent(out):: r10
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10; r10(L)=r44(i4pair(1,L),i4pair(2,L)); enddo
end subroutine r44_to_10

!==============================================================================
subroutine i4_to_10(i4,i10)!                                         [t4_to_10]
!==============================================================================
use pkind, only: spi
use pmat4, only: outer_product
implicit none
integer(spi),dimension(4), intent(in ):: i4
integer(spi),dimension(10),intent(out):: i10
!------------------------------------------------------------------------------
call t44_to_10(outer_product(i4,i4),i10)
end subroutine i4_to_10
!==============================================================================
subroutine r4_to_10(r4,r10)!                                         [t4_to_10]
!==============================================================================
use pkind, only: dp
use pmat4, only: outer_product
implicit none
real(dp),dimension(4), intent(in ):: r4
real(dp),dimension(10),intent(out):: r10
!------------------------------------------------------------------------------
call t44_to_10(outer_product(r4,r4),r10)
end subroutine r4_to_10

!==============================================================================
subroutine i10_to_44(i10,i44)!                                      [t10_to_44]
!==============================================================================
use pkind, only: spi
implicit none
integer(spi),dimension(10), intent(in ):: i10
integer(spi),dimension(4,4),intent(out):: i44
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10
   i44(i4pair(1,L),i4pair(2,L))=i10(L)
   i44(i4pair(2,L),i4pair(1,L))=i10(L)
enddo
end subroutine i10_to_44
!==============================================================================
subroutine r10_to_44(r10,r44)!                                      [t10_to_44]
!==============================================================================
use pkind, only: spi,dp
implicit none
real(dp),dimension(10), intent(in ):: r10
real(dp),dimension(4,4),intent(out):: r44
!------------------------------------------------------------------------------
integer(spi):: L
!==============================================================================
do L=1,10
   r44(i4pair(1,L),i4pair(2,L))=r10(L)
   r44(i4pair(2,L),i4pair(1,L))=r10(L)
enddo
end subroutine r10_to_44

!--------------------------------------------------------

!======================================================================= [finp]
subroutine finp
!==============================================================================
implicit none
p=2 ! revert to default value
end subroutine finp
!======================================================================= [inip]
subroutine inip(p_prescribe,ff)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u2,u3,u5,o2
implicit none
integer(spi),intent(in ):: p_prescribe
logical,     intent(out):: ff
!------------------------------------------------------------------------------
integer(spi),parameter    :: np=6
real(dp),parameter        :: u16=16,u3o2=u3*o2,u5o2=u5*o2
real(dp),dimension(np*4+1):: ffac
real(dp),dimension(np)    :: fac
integer(spi)              :: i
!==============================================================================
ff=(p_prescribe<1 .or. p_prescribe>np)
if(ff)then
   print'(" In inip; prescribed exponent p out of bounds")'
   return
endif
p=p_prescribe
rpp3o2=sqrt(p+u3o2)
rpp4o2=sqrt(p+u2)
rpp5o2=sqrt(p+u5o2)
rpp6o2=sqrt(p+u3)
ffac(1)=1
do i=3,p*4+1,2
   ffac(i)=ffac(i-2)*i
enddo
fac(1)=1
do i=2,p
   fac(i)=fac(i-1)*i
enddo
rsg=ffac(p*4+1)*u2**(p-1)/(ffac(p*2-1)*fac(p)*u16**p)
end subroutine inip

!# Calibration routines to initialize filter scale constants at each treated boundary,
!  the amplitude normalization coefficients, el(0,..), the Cholesky coefficients for
! each given overall aspect tensor, and the halo widths needed by the filters in
! each particular subdomain.
!======================================================================[rcalib]
subroutine rcalib1(Lx,mx,Lbx,mbx,asLbx,asmbx,as,xLb,xmb,el,hxm)!                            
!==============================================================================
! Convert the given field, as, of aspect tensors first into the corresponding
! scale b. Find the implied partial filter half-span s by the asymptotic
! formula, and apply the residual correction from the look-up table sres
! (if p=2, or 3) if b < 5 grid units, so that the filters' final normalized
! second moment will adequately match the intended aspect value as.
! into the equivalent field el(1,..)
! of Cholesky lower-triangular factors of the inverses of the aspect tensors
! in 1D, and put the exact normalizing coefficent into el(0,..)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u1,u2,o2
use pmat,  only: inv, l1lm
implicit none
integer(spi),                 intent(in   ):: lx,mx
logical,                      intent(in   ):: Lbx,mbx
real(dp),                     intent(in   ):: asLbx,asmbx
real(dp),    dimension(lx:mx),intent(in   ):: as
real(dp),                     intent(inout):: xLb,xmb
real(dp),dimension(0:1,lx:mx),intent(  out):: el
integer(spi),dimension(lx:mx),intent(  out):: hxm
!------------------------------------------------------------------------------
real(dp),dimension(Lx-mx:mx-Lx):: fs
real(dp)                       :: b,exx,f,r,rc,rrc,s,x
integer(spi)                   :: ib,ix,ixp,ixm,gx,gxm,gxn,Lxmix,mxmix
!==============================================================================
if(lbx)then; xLb=0; if(asLbx>u0)xLb=u1/sqrt(asLbx); endif! Lx-boundary, inverse filter spread
if(mbx)then; xmb=0; if(asmbx>u0)xmb=u1/sqrt(asmbx); endif! Mx-boundary, inverse filter spread 
do ix=Lx,mx
   b=sqrt(as(ix))
   s=rpp3o2*b! Asymptotic approximation to the single filter half-span, valid at large b
   r=b*deltabi
   ib=r
   if(ib<nsres)then ! Finesse the half-span for short lines by adding from the sres table:
      r=r-ib; rc=u1-r
      if    (p==2)then; s=s+rc*sres2(ib)+r*sres2(ib+1)
      elseif(p==3)then; s=s+rc*sres3(ib)+r*sres3(ib+1)
      endif
   endif
   exx=u1/s
   el(1,ix)=exx
   gxm=floor(u1/exx); hxm(ix)=gxm
   fs(-gxm:gxm)=0
   fs(0)=u1
   do gx=-gxm,-1
      rrc=u1-(gx*exx)**2
      f=rrc**p; fs(-gx)=f; fs(gx)=f
   enddo
   if(Lbx)then
      Lxmix=Lx-ix
      gxn=gxm+Lxmix
      do gx=1,gxn
         x=gx*xLb; if(x>u0)x=om0+x*(om1+x*om2)
         ixm=lxmix-gx
         ixp=lxmix+gx-1
         fs(ixp)=fs(ixp)+x*fs(ixm); fs(ixm)=0
      enddo
   endif
   if(mbx)then
      mxmix=mx-ix
      gxn=gxm-mxmix
      do gx=1,gxn
         x=gx*xmb; if(x>u0)x=om0+x*(om1+x*om2)
         ixp=mxmix+gx
         ixm=mxmix-gx+1
         fs(ixm)=fs(ixm)+x*fs(ixp); fs(ixp)=0
      enddo
   endif
   el(0,ix)=u1/sqrt(sum(fs(-gxm:gxm)**2))
enddo
end subroutine rcalib1
!===================================================================== [rcalib]
subroutine rcalib2(lx,mx, ly,my, lbx,mbx,lby,mby,asLbx,asmbx,asLby,asmby,as,&
     xLb,xmb,yLb,ymb,el,hxm,hym)
!==============================================================================
! Convert the given field, as, of aspect tensors into the equivalent field, el(1:,..)
! of Cholesky lower-triangular factors of the inverses of the aspect tensors
! in 2D, and put the exact normalizing coefficent into el(0,..)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u1,u2,o2
use pmat, only: inv, l1lm
implicit none
integer(spi),                       intent(in   ):: lx,mx, ly,my
logical,                            intent(in   ):: lbx,mbx,lby,mby
real(dp),                           intent(in   ):: asLbx,asmbx,asLby,asmby
real(dp),dimension(3  ,lx:mx,ly:my),intent(in   ):: as
real(dp),                           intent(inout):: xLb,xmb,yLb,ymb
real(dp),dimension(0:3,lx:mx,ly:my),intent(  out):: el
integer(spi),dimension(Lx:Mx,Ly:My),intent(  out):: hxm,hym
!------------------------------------------------------------------------------
real(dp),dimension(-mx+Lx:mx-Lx,-my+Ly:my-Ly):: fs
real(dp),dimension(2,2):: tas,tel22
real(dp),dimension(3)  :: tel
real(dp)               :: aa0,cx,exx,eyy,eyx,f,r,rrc,rrxc
integer(spi)           :: gx,gxmm,gxL,gxm,gy,gym,gxn,gyn,ix,ixm,ixp,iy,iym,iyp,&
     Lxmix,mxmix,Lymiy,mymiy,p2p4
!==============================================================================
p2p4=p*2+4
if(lbx)then; xLb=0; if(asLbx>u0)xLb=u1/sqrt(asLbx); endif! Lx-boundary inverse filter spread
if(mbx)then; xmb=0; if(asmbx>u0)xmb=u1/sqrt(asmbx); endif! Mx-boundary inverse filter spread 
if(lby)then; yLb=0; if(asLby>u0)xLb=u1/sqrt(asLby); endif! Ly-boundary inverse filter spread
if(mby)then; ymb=0; if(asmby>u0)xmb=u1/sqrt(asmby); endif! My-boundary inverse filter spread
fs=0
do iy=ly,my; do ix=lx,mx
   call t3_to_22(as(:,ix,iy),tas)
   tas=p2p4*tas/2; call inv(tas); call l1lm(tas,tel22)
   call t22_to_3(tel22,tel)
   el(1:3,ix,iy)=tel
   exx=tel(1); eyy=tel(2); eyx=tel(3)
   fs(0,0)=u1
   gym=floor(u1/eyy); gxmm=0
   lgy: do gy=-gym,0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      gxL=ceiling((cx-r)/exx); gxm=floor((cx+r)/exx); gxmm=max(-gxL,max(gxm,gxmm))
      do gx=gxL,gxm; ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         f=rrc**p
         fs(gx,gy)=f; fs(-gx,-gy)=f
      enddo! gx
   enddo lgy
   gxm=gxmm
   hxm(ix,iy)=gxm
   hym(ix,iy)=gym
! If applicable, treat domain boundary conditions   
   if(Lbx)then
      Lxmix=Lx-ix
      gxn=gxm+Lxmix
      do gx=1,gxn
         r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
         ixm=lxmix-gx
         ixp=lxmix+gx-1
         fs(ixp,-gym:gym)=fs(ixp,-gym:gym)+r*fs(ixm,-gym:gym); fs(ixm,-gym:gym)=0
      enddo
   endif
   if(mbx)then
      mxmix=mx-ix
      gxn=gxm-mxmix
      do gx=1,gxn
         r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
         ixp=mxmix+gx
         ixm=mxmix-gx+1
         fs(ixm,-gym:gym)=fs(ixm,-gym:gym)+r*fs(ixp,-gym:gym); fs(ixp,-gym:gym)=0
      enddo
   endif
   
   if(Lby)then
      Lymiy=Ly-iy
      gyn=gym+Lymiy
      do gy=1,gyn
         r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
         iym=lymiy-gy
         iyp=lymiy+gy-1
         fs(-gxm:gxm,iyp)=fs(-gxm:gxm,iyp)+r*fs(-gxm:gxm,iym); fs(-gxm:gxm,iym)=0
      enddo
   endif
   if(mby)then
      mymiy=my-iy
      gyn=gym-mymiy
      do gy=1,gyn
         r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
         iyp=mymiy+gy
         iym=mymiy-gy+1
         fs(-gxm:gxm,iym)=fs(-gxm:gxm,iym)+r*fs(-gxm:gxm,iyp); fs(-gxm:gxm,iyp)=0
      enddo
   endif
   el(0,ix,iy)=u1/sqrt(sum(fs(-gxm:gxm,-gym:gym)**2)); fs(-gxm:gxm,-gym:gym)=0
enddo;       enddo
end subroutine rcalib2
!===================================================================== [rcalib]
subroutine rcalib3(Lx,mx, Ly,my, Lz,mz, Lbx,mbx,Lby,mby,Lbz,mbz,&
     asLbx,asmbx,asLby,asmby,asLbz,asmbz,as, xLb,xmb,yLb,ymb,zLb,zmb,el,hxm,hym,hzm)
!==============================================================================
! Convert the given field, as, of aspect tensors into the equivalent field el(1:,..)
! of Cholesky lower-triangular factors of the inverses of the aspect tensors
! in 3D, and put the exact normalizing coefficient in el(0,..)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u2,o2
use pmat, only: inv, l1lm
implicit none
integer(spi),                             intent(in   ):: Lx,mx, Ly,my, Lz,mz
logical,                                  intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz
real(dp),                                 intent(in   ):: asLbx,asmbx,asLby,asmby,asLbz,asmbz
real(dp),dimension(  6,Lx:mx,Ly:my,Lz:mz),intent(in   ):: as
real(dp),                                 intent(inout):: xLb,xmb,yLb,ymb,zLb,zmb
real(dp),dimension(0:6,Lx:mx,Ly:my,Lz:mz),intent(  out):: el
integer(spi),dimension(Lx:mx,Ly:my,Lz:mz),intent(  out):: hxm,hym,hzm
!------------------------------------------------------------------------------
real(dp),dimension(-mx+Lx:mx-Lx,-my+Ly:my-Ly,-mz+Lz:mz-Lz):: fs
real(dp),dimension(3,3):: tas,teL33
real(dp),dimension(3)  :: g
real(dp),dimension(6)  :: teL
real(dp)               :: aa0,cx,cy,czx,exx,eyy,ezz,eyx,ezx,ezy,&
                          f,r,rrc,rrxc,rryc
integer(spi)           :: gx,gxL,gxm,gxn,gy,gyL,gym,gyn,gz,gzm,gzn,gxmm,gymm,&
     ix,ixp,ixm,iy,iyp,iym,iz,izp,izm,&
     Lxmix,mxmix,Lymiy,mymiy,Lzmiz,mzmiz,p2p5
!==============================================================================
p2p5=p*2+5
if(lbx)then; xLb=0; if(asLbx>u0)xLb=u1/sqrt(asLbx); endif! Lx-boundary inverse filter spread
if(mbx)then; xmb=0; if(asmbx>u0)xmb=u1/sqrt(asmbx); endif! Mx-boundary inverse filter spread 
if(lby)then; yLb=0; if(asLby>u0)xLb=u1/sqrt(asLby); endif! Ly-boundary inverse filter spread
if(mby)then; ymb=0; if(asmby>u0)xmb=u1/sqrt(asmby); endif! My-boundary inverse filter spread 
if(lbz)then; zLb=0; if(asLbz>u0)zLb=u1/sqrt(asLbz); endif! Lz-boundary inverse filter spread
if(mbz)then; zmb=0; if(asmbz>u0)zmb=u1/sqrt(asmbz); endif! Mz-boundary inverse filter spread
fs=0
do iz=Lz,mz; do iy=Ly,my; do ix=Lx,mx
   call t6_to_33(as(:,ix,iy,iz),tas)
   tas=p2p5*tas/2; call inv(tas); call L1Lm(tas,tel33)
   call t33_to_6(teL33,teL)
   eL(1:6,ix,iy,iz)=teL
   exx=tel(1);eyy=tel(2);ezz=tel(3);ezy=tel(4);ezx=tel(5);eyx=tel(6)
   fs(0,0,0)=u1
   gzm=floor(u1/ezz); ; gxmm=0; gymm=0
   lgz: do gz=-gzm,0; izp=iz+gz; izm=iz-gz
      rryc=abs(u1-(gz*ezz)**2); r=sqrt(rryc); cy=-gz*ezy; czx=-gz*ezx
      gyL=ceiling((cy-r)/eyy); gym=floor((cy+r)/eyy); gymm=max(-gyL,max(gym,gymm))
      do gy=gyL,gym; iyp=iy+gy; iym=iy-gy
         rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
         gxL=ceiling((cx-r)/exx); gxm=floor((cx+r)/exx); gxmm=max(-gxL,max(gxm,gxmm))
         do gx=gxL,gxm; ixp=ix+gx; ixm=ix-gx
            if(gz==0.and.gy==0.and.gx==0)exit lgz
            rrc=rrxc-(gx*exx-cx)**2
            f=rrc**p
            fs(gx,gy,gz)=f; fs(-gx,-gy,-gz)=f
         enddo! gx
      enddo! gy
   enddo lgz
   gxm=gxmm
   gym=gymm
   hxm(ix,iy,iz)=gxm
   hym(ix,iy,iz)=gym
   hzm(ix,iy,iz)=gzm
! If applicable, treat domain boundary conditions
   if(Lbx)then
      Lxmix=Lx-ix
      gxn=gxm+Lxmix
      do gx=1,gxn
         r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
         ixm=Lxmix-gx
         ixp=Lxmix+gx-1
         fs(ixp,-gym:gym,-gzm:gzm)=fs(ixp,-gym:gym,-gzm:gzm)+r*fs(ixm,-gym:gym,-gzm:gzm)
         fs(ixm,-gym:gym,-gzm:gzm)=0
      enddo
   endif
   if(mbx)then
      mxmix=mx-ix
      gxn=gxm-mxmix
      do gx=1,gxn
         r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
         ixp=mxmix+gx
         ixm=mxmix-gx+1
         fs(ixm,-gym:gym,-gzm:gzm)=fs(ixm,-gym:gym,-gzm:gzm)+r*fs(ixp,-gym:gym,-gzm:gzm)
         fs(ixp,-gym:gym,-gzm:gzm)=0
      enddo
   endif

   if(Lby)then
      Lymiy=Ly-iy
      gyn=gym+Lymiy
      do gy=1,gyn
         r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
         iym=Lymiy-gy
         iyp=Lymiy+gy-1
         fs(-gxm:gxm,iyp,-gzm:gzm)=fs(-gxm:gxm,iyp,-gzm:gzm)+r*fs(-gxm:gxm,iym,-gzm:gzm)
         fs(-gxm:gxm,iym,-gzm:gzm)=0
      enddo
   endif
   if(mby)then
      mymiy=my-iy
      gyn=gym-mymiy
      do gy=1,gyn
         r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
         iyp=mymiy+gy
         iym=mymiy-gy+1
         fs(-gxm:gxm,iym,-gzm:gzm)=fs(-gxm:gxm,iym,-gzm:gzm)+r*fs(-gxm:gxm,iyp,-gzm:gzm)
         fs(-gxm:gxm,iyp,-gzm:gzm)=0
      enddo
   endif
   
   if(Lbz)then
      Lzmiz=Lz-iz
      gzn=gzm+Lzmiz
      do gz=1,gzn
         r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
         izm=Lzmiz-gz
         izp=Lzmiz+gz-1
         fs(-gxm:gxm,-gym:gym,izp)=fs(-gxm:gxm,-gym:gym,izp)+r*fs(-gxm:gxm,-gym:gym,izm)
         fs(-gxm:gxm,-gym:gym,izm)=0
      enddo
   endif
   if(mbz)then
      mzmiz=mz-iz
      gzn=gzm-mzmiz
      do gz=1,gzn
         r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
         izp=mzmiz+gz
         izm=mzmiz-gz+1
         fs(-gxm:gxm,-gym:gym,izm)=fs(-gxm:gxm,-gym:gym,izm)+r*fs(-gxm:gxm,-gym:gym,izp)
         fs(-gxm:gxm,-gym:gym,izp)=0
      enddo
   endif
   
   el(0,ix,iy,iz)=u1/sqrt(sum(fs(-gxm:gxm,-gym:gym,-gzm:gzm)**2))
   fs(-gxm:gxm,-gym:gym,-gzm:gzm)=0
enddo;       enddo;       enddo
end subroutine rcalib3
!======================================================================[rcalib]
subroutine rcalib4(Lx,mx,Ly,my,Lz,mz,Lw,mw, Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw,&
     asLbx,asmbx,asLby,asmby,asLbz,asmbz,asLbw,asmbw,as,&
     xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb,el,hxm,hym,hzm,hwm)
!==============================================================================
! Convert the given field, as, of aspect tensors into the equivalent field, el(1:,..)
! of Cholesky lower-triangular factors of the inverses of the aspect tensors
! in 4D, and put the exact normalizing coefficent into el(0,..)
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u0,u2,o2
use pmat,  only: inv, L1Lm
implicit none
integer(spi),intent(in   ):: lx,mx, ly,my, lz,mz, lw,mw
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw
real(dp),    intent(in   ):: asLbx,asmbx,asLby,asmby,asLbz,asmbz,asLbw,asmbw
real(dp),dimension(  10,Lx:mx,Ly:my,Lz:mz,Lw:mw),intent(in   ):: as
real(dp),    intent(inout):: xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb
real(dp),dimension(0:10,Lx:mx,Ly:my,Lz:mz,Lw:mw),intent(  out):: el
integer(spi), dimension(Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(  out):: hxm,hym,hzm,hwm
!------------------------------------------------------------------------------
real(dp),dimension(-mx+Lx:mx-Lx,-my+Ly:my-Ly,-mz+Lz:mz-Lz,-mw+Lw:mw-Lw):: fs
real(dp),dimension(4,4):: tas,tel44
real(dp),dimension(10) :: tel
real(dp)               :: cx,cy,cz,cwy,cwx,czx,&
                          exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,&
                          f,r,rrc,rrxc,rryc,rrzc
integer(spi)           :: gx,gxL,gxm,gxn,gy,gyL,gym,gyn,gz,gzL,gzm,gzn,gw,gwm,gwn,&
     gxmm,gymm,gzmm,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm,iw,iwp,iwm,&
          Lxmix,mxmix,Lymiy,mymiy,Lzmiz,mzmiz,Lwmiw,mwmiw,p2p6
!==============================================================================
p2p6=p*2+6
if(lbx)then; xLb=0; if(asLbx>u0)xLb=u1/sqrt(asLbx); endif! Lx-boundary inverse filter spread
if(mbx)then; xmb=0; if(asmbx>u0)xmb=u1/sqrt(asmbx); endif! Mx-boundary inverse filter spread 
if(lby)then; yLb=0; if(asLby>u0)xLb=u1/sqrt(asLby); endif! Ly-boundary inverse filter spread
if(mby)then; ymb=0; if(asmby>u0)xmb=u1/sqrt(asmby); endif! My-boundary inverse filter spread 
if(lbz)then; zLb=0; if(asLbz>u0)zLb=u1/sqrt(asLbz); endif! Lz-boundary inverse filter spread
if(mbz)then; zmb=0; if(asmbz>u0)zmb=u1/sqrt(asmbz); endif! Mz-boundary inverse filter spread 
if(lbw)then; wLb=0; if(asLbw>u0)wLb=u1/sqrt(asLbw); endif! Lw-boundary inverse filter spread
if(mbw)then; wmb=0; if(asmbw>u0)wmb=u1/sqrt(asmbw); endif! Mw-boundary inverse filter spread 
fs=0
do iw=lw,mw; do iz=lz,mz; do iy=ly,my; do ix=lx,mx
   call t10_to_44(as(:,ix,iy,iz,iw),tas)
   tas=p2p6*tas/2; call inv(tas); call L1Lm(tas,tel44)
   call t44_to_10(tel44,tel)
   el(1:10,ix,iy,iz,iw)=tel
   exx=tel(1);eyy=tel(2);ezz=tel(3);eww=tel(4);eyx=tel(5)
   ezx=tel(6);ewx=tel(7);ezy=tel(8);ewy=tel(9);ewz=tel(10)
   fs(0,0,0,0)=u1
   gwm=floor(u1/eww); gxmm=0; gymm=0; gzmm=0
   lgw: do gw=-gwm,0; iwp=iw+gw; iwm=iw-gw
      rrzc=abs(u1-(gw*eww)**2); r=sqrt(rrzc); cz=-gw*ewz; cwy=-gw*ewy; cwx=-gw*ewx
      gzL=ceiling((cz-r)/ezz); gzm=floor((cz+r)/ezz); gzmm=max(-gzL,max(gzm,gzmm))
      do gz=-gzL,gzm; izp=iz+gz; izm=iz-gz
         rryc=abs(rrzc-(gz*ezz-cz)**2); r=sqrt(rryc); cy=cwy-gz*ezy; czx=cwx-gz*ezx
         gyL=ceiling((cy-r)/eyy); gym=floor((cy+r)/eyy); gymm=max(-gyL,max(gym,gymm))
         do gy=gyL,gym; iyp=iy+gy; iym=iy-gy
            rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
            gxL=ceiling((cx-r)/exx); gxm=floor((cx+r)/exx); gxmm=max(-gxL,max(gxm,gxmm))
            do gx=gxL,gxm; ixp=ix+gx; ixm=ix-gx
               if(gw==0.and.gz==0.and.gy==0.and.gx==0)exit lgw
               rrc=rrxc-(gx*exx-cx)**2
               f=rrc**p
               fs(gx,gy,gz,gw)=f; fs(-gx,-gy,-gz,-gw)=f
            enddo! gx
         enddo! gy
      enddo! gz
   enddo lgw
   gxm=gxmm
   gym=gymm
   gzm=gzmm
   hxm(ix,iy,iz,iw)=gxm
   hym(ix,iy,iz,iw)=gym
   hzm(ix,iy,iz,iw)=gzm
   hwm(ix,iy,iz,iw)=gwm

! If applicable, treat domain boundary conditions
   if(Lbx)then
      Lxmix=Lx-ix
      gxn=gxm+Lxmix
      do gx=1,gxn
         r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
         ixm=lxmix-gx
         ixp=lxmix+gx-1
         fs(ixp,-gym:gym,-gzm:gzm,-gwm:gwm)=fs(ixp,-gym:gym,-gzm:gzm,-gwm:gwm)&
              +r*fs(ixm,-gym:gym,-gzm:gzm,-gwm:gwm)
         fs(ixm,-gym:gym,-gzm:gzm,-gwm:gwm)=0
      enddo
   endif
   if(mbx)then
      mxmix=mx-ix
      gxn=gxm-mxmix
      do gx=1,gxn
         r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
         ixp=mxmix+gx
         ixm=mxmix-gx+1
         fs(ixm,-gym:gym,-gzm:gzm,-gwm:gwm)=fs(ixm,-gym:gym,-gzm:gzm,-gwm:gwm)&
              +r*fs(ixp,-gym:gym,-gzm:gzm,-gwm:gwm)
         fs(ixp,-gym:gym,-gzm:gzm,-gwm:gwm)=0
      enddo
   endif

   if(Lby)then
      Lymiy=Ly-iy
      gyn=gym+Lymiy
      do gy=1,gyn
         r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
         iym=lymiy-gy
         iyp=lymiy+gy-1
         fs(-gxm:gxm,iyp,-gzm:gzm,-gwm:gwm)=fs(-gxm:gxm,iyp,-gzm:gzm,-gwm:gwm)&
              +r*fs(-gxm:gxm,iym,-gzm:gzm,-gwm:gwm)
         fs(-gxm:gxm,iym,-gzm:gzm,-gwm:gwm)=0
      enddo
   endif
   if(mby)then
      mymiy=my-iy
      gyn=gym-mymiy
      do gy=1,gyn
         r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
         iyp=mymiy+gy
         iym=mymiy-gy+1
         fs(-gxm:gxm,iym,-gzm:gzm,-gwm:gwm)=fs(-gxm:gxm,iym,-gzm:gzm,-gwm:gwm)&
              +r*fs(-gxm:gxm,iyp,-gzm:gzm,-gwm:gwm)
         fs(-gxm:gxm,iyp,-gzm:gzm,-gwm:gwm)=0
      enddo
   endif

   if(Lbz)then
      Lzmiz=Lz-iz
      gzn=gzm+Lzmiz
      do gz=1,gzn
         r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
         izm=lzmiz-gz
         izp=lzmiz+gz-1
         fs(-gxm:gxm,-gym:gym,izp,-gwm:gwm)=fs(-gxm:gxm,-gym:gym,izp,-gwm:gwm)&
              +r*fs(-gxm:gxm,-gym:gym,izm,-gwm:gwm)
         fs(-gxm:gxm,-gym:gym,izm,-gwm:gwm)=0
      enddo
   endif
   if(mbz)then
      mzmiz=mz-iz
      gzn=gzm-mzmiz
      do gz=1,gzn
         r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
         izp=mzmiz+gz
         izm=mzmiz-gz+1
         fs(-gxm:gxm,-gym:gym,izm,-gwm:gwm)=fs(-gxm:gxm,-gym:gym,izm,-gwm:gwm)&
              +r*fs(-gxm:gxm,-gym:gym,izp,-gwm:gwm)
         fs(-gxm:gxm,-gym:gym,izp,-gwm:gwm)=0
      enddo
   endif

   if(Lbw)then
      Lwmiw=Lw-iw
      gwn=gwm+Lwmiw
      do gw=1,gwn
         r=gw*wLb; if(r>u0)r=om0+r*(om1+r*om2)
         iwm=lwmiw-gw
         iwp=lwmiw+gw-1
         fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwp)=fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwp)&
              +r*fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwm)
         fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwm)=0
      enddo
   endif
   if(mbw)then
      mwmiw=mw-iw
      gwn=gwm-mwmiw
      do gw=1,gwn
         r=gw*wmb; if(r>u0)r=om0+r*(om1+r*om2)
         iwp=mwmiw+gw
         iwm=mwmiw-gw+1
         fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwm)=fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwm)&
              +r*fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwp)
         fs(-gxm:gxm,-gym:gym,-gwm:gwm,iwp)=0
      enddo
   endif

   el(0,ix,iy,iz,iw)=u1/sqrt(sum(fs(-gxm:gxm,-gym:gym,-gzm:gzm,-gwm:gwm)**2))
   fs(-gxm:gxm,-gym:gym,-gzm:gzm,-gwm:gwm)=0
enddo;       enddo;       enddo;       enddo
end subroutine rcalib4

!# Adjoint beta filters:
!======================================================================[rbetat]
subroutine rbeta1T(hx,Lx,mx, el, a,b)
!=============================================================================
! Perform an ADJOINT beta-function filter in 1D.
!=============================================================================
integer(spi),                    intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:mx   ),intent(in   ):: el
real(dp),dimension(    Lx:mx   ),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp)    :: tafrow,tas
real(dp)    :: exx,rrc
integer(spi):: ix,ixp,ixm,gx
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tas=a(ix)*el(0,ix)
   b(ix)=b(ix)+tas
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tafrow=tas*rrc**p
      b(ixp)=b(ixp)+tafrow
      b(ixm)=b(ixm)+tafrow
   enddo
enddo
end subroutine rbeta1t
!======================================================================[rbetat]
subroutine vrbeta1T(nv, hx,lx,mx, el, a,b)
!=============================================================================
! Vector version of rbeta1t filtering nv fields at once.
!=============================================================================
integer(spi),                      intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(0:1,  Lx:mx),   intent(in   ):: el
real(dp),dimension(nv,   Lx:mx),   intent(in   ):: a
real(dp),dimension(nv,Lx-hx:mx+hx),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(nv):: tafrow,tas
real(dp)              :: exx,rrc
integer(spi)          :: ix,ixp,ixm,gx
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tas=a(:,ix)*el(0,ix)
   b(:,ix)=b(:,ix)+tas
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tafrow=tas*rrc**p
      b(:,ixp)=b(:,ixp)+tafrow
      b(:,ixm)=b(:,ixm)+tafrow
   enddo
enddo
end subroutine vrbeta1t
!======================================================================[rbetat]
subroutine rbeta2T(hx,lx,mx, hy,ly,my, el, a,b)
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 2D.
!=============================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my
real(dp),dimension(0:3,Lx:Mx,       Ly:My   ),intent(in   ):: el
real(dp),dimension(    Lx:mx,       Ly:my   ),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3):: tel
real(dp)             :: tafrow,tas
real(dp)             :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)         :: gx,gy,ix,ixp,ixm,iy,iyp,iym
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tas=a(ix,iy)*el(0,ix,iy)
   b(ix,iy)=b(ix,iy)+tas
Lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit Lgy
         rrc=rrxc-(gx*exx-cx)**2
         tafrow=tas*rrc**p
         b(ixp,iyp)=b(ixp,iyp)+tafrow
         b(ixm,iym)=b(ixm,iym)+tafrow
      enddo! gx
   enddo Lgy
enddo;  enddo! ix, iy
end subroutine rbeta2t
!======================================================================[rbetat]
subroutine vrbeta2T(nv,hx,Lx,mx, hy,Ly,my, el, a,b)
!=============================================================================
! Vector version of rbeta2t filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,Ly,my
real(dp),dimension( 0:3,Lx:Mx,Ly:My),          intent(in   ):: el
real(dp),dimension(nv,  Lx:mx,Ly:my),          intent(in   ):: a
real(dp),dimension(nv,Lx-hx:mx+hx,Ly-hy:my+hy),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3) :: tel
real(dp),dimension(nv):: tafrow,tas
real(dp)              :: cx,rrc,rrxc,exx,eyy,eyx,r
integer(spi)          :: ix,ixp,ixm,gx,iy,iyp,iym,gy
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tas=a(:,ix,iy)*el(0,ix,iy)
   b(:,ix,iy)=b(:,ix,iy)+tas
Lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit Lgy
         rrc=rrxc-(gx*exx-cx)**2
         tafrow=tas*rrc**p
         b(:,ixp,iyp)=b(:,ixp,iyp)+tafrow
         b(:,ixm,iym)=b(:,ixm,iym)+tafrow
      enddo! gx
   enddo Lgy
enddo; enddo ! ix, iy
end subroutine vrbeta2t
!======================================================================[rbetat]
subroutine rbeta3T(hx,Lx,mx, hy,Ly,my, hz,Lz,mz, el, a,b)
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 3D.
!=============================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,Ly,my,hz,Lz,mz
real(dp),dimension(0:6,Lx:Mx,       Ly:My,       Lz:Mz   ),intent(in   ):: el
real(dp),dimension(    Lx:mx,       Ly:my,       Lz:mz   ),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(6):: tel
real(dp)             :: tafrow,tas
real(dp)             :: cx,cy,czx,exx,eyy,ezz,eyx,ezx,ezy,r,rrc,rrxc,rryc
integer(spi)         :: gx,gy,gz,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:6,ix,iy,iz)
   exx=tel(1);eyy=tel(2);ezz=tel(3);ezy=tel(4);ezx=tel(5);eyx=tel(6)
   tas=a(ix,iy,iz)*el(0,ix,iy,iz)
   b(ix,iy,iz)=b(ix,iy,iz)+tas
lgz: do gz=ceiling(-u1/ezz),0; izp=iz+gz; izm=iz-gz
        rryc=abs(u1-(gz*ezz)**2); r=sqrt(rryc); cy=-gz*ezy; czx=-gz*ezx
      do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
         rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
         do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
            if(gz==0.and.gy==0.and.gx==0)exit lgz
            rrc=rrxc-(gx*exx-cx)**2
            tafrow=tas*rrc**p
            b(ixp,iyp,izp)=b(ixp,iyp,izp)+tafrow
            b(ixm,iym,izm)=b(ixm,iym,izm)+tafrow
         enddo! gx
      enddo! gy
   enddo lgz
enddo;  enddo;  enddo ! ix, iy, iz
end subroutine rbeta3t
!=====================================================================[rbetat]
subroutine vrbeta3T(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, el, a,b)
!=============================================================================
! Vector version of rbeta3t filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my,hz,lz,mz
real(dp),dimension(0:6,Lx:mx,Ly:my,Lz:mz),intent(in   ):: el
real(dp),dimension(nv, Lx:mx,Ly:my,Lz:mz),intent(in   ):: a
real(dp),dimension(nv,Lx-hx:mx+hx,Ly-hy:my+hy,Lz-hz:mz+hz),&
                                          intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(6) :: tel
real(dp),dimension(nv):: tafrow,tas
real(dp)              :: cx,cy,czx,exx,eyy,ezz,eyx,ezx,ezy,r,rrc,rrxc,rryc
integer(spi)          :: gx,gy,gz,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:6,ix,iy,iz)
   exx=tel(1);eyy=tel(2);ezz=tel(3);ezy=tel(4);ezx=tel(5);eyx=tel(6)
   tas=a(:,ix,iy,iz)*el(0,ix,iy,iz)
   b(:,ix,iy,iz)=b(:,ix,iy,iz)+tas
Lgz: do gz=ceiling(-u1/ezz),0; izp=iz+gz; izm=iz-gz
      rryc=abs(u1-(gz*ezz)**2); r=sqrt(rryc); cy=-gz*ezy; czx=-gz*ezx
      do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
         rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
         do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
            if(gz==0.and.gy==0.and.gx==0)exit Lgz
            rrc=rrxc-(gx*exx-cx)**2
            tafrow=tas*rrc**p
            b(:,ixp,iyp,izp)=b(:,ixp,iyp,izp)+tafrow
            b(:,ixm,iym,izm)=b(:,ixm,iym,izm)+tafrow
         enddo! gx
      enddo! gy
   enddo Lgz
enddo; enddo; enddo! ix, iy, iz
end subroutine vrbeta3t
!=============================================================================[rbetat]
subroutine rbeta4T(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el, a,b)
!=============================================================================
! Perform an ADJOINT radial beta-function filter in 4D.
!=============================================================================
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw
real(dp),dimension(0:10,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(     Lx:mx,Ly:my,Lz:mz,Lw:mw),intent(in   ):: a
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz,-hw+Lw:mw+hw),&
                                                intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(10):: tel
real(dp)              :: tafrow,tas
real(dp):: cx,cy,cz,cwy,cwx,czx,exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,&
           r,rrc,rrxc,rryc,rrzc
integer(spi) :: gx,gy,gz,gw,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm,iw,iwp,iwm
!=============================================================================
b=0
do iw=Lw,Mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:10,ix,iy,iz,iw)
   exx=tel(1);eyy=tel(2);ezz=tel(3);eww=tel(4);eyx=tel(5)
   ezx=tel(6);ewx=tel(7);ezy=tel(8);ewy=tel(9);ewz=tel(10)
   tas=a(ix,iy,iz,iw)*el(0,ix,iy,iz,iw)
   b(ix,iy,iz,iw)=b(ix,iy,iz,iw)+tas
Lgw: do gw=ceiling(-u1/eww),0; iwp=iw+gw; iwm=iw-gw
      rrzc=abs(u1-(gw*eww)**2); r=sqrt(rrzc); cz=-gw*ewz; cwy=-gw*ewy; cwx=-gw*ewx
      do gz=ceiling((cz-r)/ezz),floor((cz+r)/ezz); izp=iz+gz; izm=iz-gz
         rryc=abs(rrzc-(gz*ezz-cz)**2); r=sqrt(rryc); cy=cwy-gz*ezy; czx=cwx-gz*ezx
         do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
            rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
            do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
               if(gw==0.and.gz==0.and.gy==0.and.gx==0)exit Lgw
               rrc=rrxc-(gx*exx-cx)**2
               tafrow=tas*rrc**p
               b(ixp,iyp,izp,iwp)=b(ixp,iyp,izp,iwp)+tafrow
               b(ixm,iym,izm,iwm)=b(ixm,iym,izm,iwm)+tafrow
            enddo! gx
         enddo! gy
      enddo! gz
   enddo Lgw
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
end subroutine rbeta4t
!======================================================================[rbetat]
subroutine vrbeta4T(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el, a,b)
!=============================================================================
! Vector version of rbeta4t filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
real(dp),dimension(0:10,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(nv,  lx:mx,ly:my,lz:mz,lw:mw),intent(in   ):: a
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),&
                                               intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(10):: tel
real(dp),dimension(nv):: tafrow,tas
real(dp)    :: cx,cy,cz,cwy,cwx,czx,exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,&
               r,rrc,rrxc,rryc,rrzc
integer(spi):: ix,ixp,ixm,gx,iy,iyp,iym,gy,iz,izp,izm,gz,iw,iwp,iwm,gw
!=============================================================================
b=0
do iw=Lw,Mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:10,ix,iy,iz,iw)
   exx=tel(1);eyy=tel(2);ezz=tel(3);eww=tel(4);eyx=tel(5)
   ezx=tel(6);ewx=tel(7);ezy=tel(8);ewy=tel(9);ewz=tel(10)
   tas=a(:,ix,iy,iz,iw)*el(0,ix,iy,iz,iw)
   b(:,ix,iy,iz,iw)=b(:,ix,iy,iz,iw)+tas
Lgw: do gw=ceiling(-u1/eww),0; iwp=iw+gw; iwm=iw-gw
      rrzc=abs(u1-(gw*eww)**2); r=sqrt(rrzc); cz=-gw*ewz; cwy=-gw*ewy; cwx=-gw*ewx
      do gz=ceiling((cz-r)/ezz),floor((cz+r)/ezz); izp=iz+gz; izm=iz-gz
         rryc=abs(rrzc-(gz*ezz-cz)**2); r=sqrt(rryc); cy=cwy-gz*ezy; czx=cwx-gz*ezx
         do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
            rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
            do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
               if(gw==0.and.gz==0.and.gy==0.and.gx==0)exit Lgw
               rrc=rrxc-(gx*exx-cx)**2
               tafrow=tas*rrc**p
               b(:,ixp,iyp,izp,iwp)=b(:,ixp,iyp,izp,iwp)+tafrow
               b(:,ixm,iym,izm,iwm)=b(:,ixm,iym,izm,iwm)+tafrow
            enddo! gx
         enddo! gy
      enddo! gz
   enddo Lgw
enddo; enddo; enddo; enddo! ix, iy, iz, iw
end subroutine vrbeta4t

!# Routines for the treatment of adjoint filter results at domain boundaries:
!========================================================================[flipt]
subroutine flip1t(hx,Lx,mx,Lb,mb, xLb,xmb, a)
!===============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                   intent(in   ):: hx,Lx,mx
logical,                        intent(in   ):: Lb,mb
real(dp),                       intent(in   ):: xLb,xmb
real(dp),dimension(Lx-hx:mx+hx),intent(inout):: a
!-------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,ix,ixp,ixm,Lxm,mxp
!===============================================================================
if(Lb)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add; set the halo to 0
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp)=a(ixp)+r*a(ixm); a(ixm)=0
   enddo
endif
if(mb)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm)=a(ixm)+r*a(ixp); a(ixp)=0
   enddo
endif
end subroutine flip1t
!========================================================================[flipt]
subroutine vflip1t(nv,hx,Lx,mx,Lb,mb, xLb,xmb, a)
!===============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                      intent(in   ):: nv
integer(spi),                      intent(in   ):: hx,Lx,mx
logical,                           intent(in   ):: Lb,mb
real(dp),                          intent(in   ):: xLb,xmb
real(dp),dimension(nv,Lx-hx:mx+hx),intent(inout):: a
!-------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,ix,ixp,ixm,Lxm,mxp
!===============================================================================
if(Lb)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add; set the halo to 0
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp)=a(:,ixp)+r*a(:,ixm); a(:,ixm)=0
   enddo
endif
if(mb)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm)=a(:,ixm)+r*a(:,ixp); a(:,ixp)=0
   enddo
endif
end subroutine vflip1t
!=======================================================================[flipt]
subroutine flip2t(hx,Lx,mx, hy,Ly,my, Lbx,mbx,Lby,mby, xLb,xmb,yLb,ymb, a)
!==============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),intent(in   ):: hx,Lx,mx,hy,Ly,my
logical,     intent(in   ):: Lbx,mbx,Lby,mby
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(inout):: a
!------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,ixm,ixp,iym,iyp,kx,ky,Lxm,mxp,Lym,myp,nx,ny
!==============================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp)=a(kx:nx,iyp)+r*a(kx:nx,iym); a(kx:nx,iym)=0
   enddo!       gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym)=a(kx:nx,iym)+r*a(kx:nx,iyp); a(kx:nx,iyp)=0
   enddo!       gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny)=a(ixp,ky:ny)+r*a(ixm,ky:ny); a(ixm,ky:ny)=0
   enddo!       gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny)=a(ixm,ky:ny)+r*a(ixp,ky:ny); a(ixp,ky:ny)=0
   enddo!       gx
endif
end subroutine flip2t
!=======================================================================[flipt]
subroutine vflip2t(nv, hx,Lx,mx, hy,Ly,my, Lbx,mbx,Lby,mby, xLb,xmb,yLb,ymb, a)
!==============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),intent(in   ):: nv
integer(spi),intent(in   ):: hx,Lx,mx,hy,Ly,my
logical,     intent(in   ):: Lbx,mbx,Lby,mby
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(inout):: a
!------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,ixm,ixp,iym,iyp,kx,ky,Lxm,mxp,Lym,myp,nx,ny
!==============================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp)=a(:,kx:nx,iyp)+r*a(:,kx:nx,iym); a(:,kx:nx,iym)=0
   enddo!       gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym)=a(:,kx:nx,iym)+r*a(:,kx:nx,iyp); a(:,kx:nx,iyp)=0
   enddo!       gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny)=a(:,ixp,ky:ny)+r*a(:,ixm,ky:ny); a(:,ixm,ky:ny)=0
   enddo!       gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny)=a(:,ixm,ky:ny)+r*a(:,ixp,ky:ny); a(:,ixp,ky:ny)=0
   enddo!       gx
endif
end subroutine vflip2t

!============================================================================[flipt]
subroutine flip3t(hx,Lx,mx, hy,Ly,my, hz,Lz,mz,&
     Lbx,mbx,Lby,mby,Lbz,mbz, xLb,xmb,yLb,ymb,zLb,zmb, a)
!==================================================================================
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(inout):: a
!----------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,ixp,ixm,iyp,iym,izp,izm,kx,ky,kz, &
               Lxm,mxp,Lym,myp,Lzm,mzp,nx,ny,nz
!==================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz

if(Lbz)then; Lzm=Lz-1! reflect exterior halo <Lz inwards and add:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izp)=a(kx:nx,ky:ny,izp)+r*a(kx:nx,ky:ny,izm); a(kx:nx,ky:ny,izm)=0
   enddo!            gz
endif
if(mbz)then; mzp=mz+1! reflect exterior halo >mz inwards and add:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izm)=a(kx:nx,ky:ny,izm)+r*a(kx:nx,ky:ny,izp); a(kx:nx,ky:ny,izp)=0
   enddo!           gz
endif

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp,kz:nz)=a(kx:nx,iyp,kz:nz)+r*a(kx:nx,iym,kz:nz); a(kx:nx,iym,kz:nz)=0
   enddo!        gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym,kz:nz)=a(kx:nx,iym,kz:nz)+r*a(kx:nx,iyp,kz:nz); a(kx:nx,iyp,kz:nz)=0
   enddo!        gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny,kz:nz)=a(ixp,ky:ny,kz:nz)+r*a(ixm,ky:ny,kz:nz); a(ixm,ky:ny,kz:nz)=0
   enddo!     gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny,kz:nz)=a(ixm,ky:ny,kz:nz)+r*a(ixp,ky:ny,kz:nz); a(ixp,ky:ny,kz:nz)=0
   enddo!     gx
endif
end subroutine flip3t
!============================================================================[flipt]
subroutine vflip3t(nv, hx,Lx,mx, hy,Ly,my, hz,Lz,mz,&
     Lbx,mbx,Lby,mby,Lbz,mbz, xLb,xmb,yLb,ymb,zLb,zmb, a)
!==================================================================================
integer(spi),intent(in   ):: nv
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(inout):: a
!----------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,ixp,ixm,iyp,iym,izp,izm,kx,ky,kz, &
               Lxm,mxp,Lym,myp,Lzm,mzp,nx,ny,nz
!==================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz

if(Lbz)then; Lzm=Lz-1! reflect exterior halo <Lz inwards and add:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izp)=a(:,kx:nx,ky:ny,izp)+r*a(:,kx:nx,ky:ny,izm); a(:,kx:nx,ky:ny,izm)=0
   enddo!            gz
endif
if(mbz)then; mzp=mz+1! reflect exterior halo >mz inwards and add:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izm)=a(:,kx:nx,ky:ny,izm)+r*a(:,kx:nx,ky:ny,izp); a(:,kx:nx,ky:ny,izp)=0
   enddo!           gz
endif

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp,kz:nz)=a(:,kx:nx,iyp,kz:nz)+r*a(:,kx:nx,iym,kz:nz); a(:,kx:nx,iym,kz:nz)=0
   enddo!        gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym,kz:nz)=a(:,kx:nx,iym,kz:nz)+r*a(:,kx:nx,iyp,kz:nz); a(:,kx:nx,iyp,kz:nz)=0
   enddo!        gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny,kz:nz)=a(:,ixp,ky:ny,kz:nz)+r*a(:,ixm,ky:ny,kz:nz); a(:,ixm,ky:ny,kz:nz)=0
   enddo!     gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny,kz:nz)=a(:,ixm,ky:ny,kz:nz)+r*a(:,ixp,ky:ny,kz:nz); a(:,ixp,ky:ny,kz:nz)=0
   enddo!     gx
endif
end subroutine vflip3t
!===========================================================================[flipt]
subroutine flip4t(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw,  xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb,a)
!=================================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz,-hw+Lw:mw+hw),&
             intent(inout):: a
!---------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,gw,ixp,ixm,iyp,iym,izp,izm,iwp,iwm,&
               kx,ky,kz,kw,Lxm,mxp,Lym,myp,Lzm,mzp,Lwm,mwp,nx,ny,nz,nw
!=================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz; kw=Lw-hw; nw=mw+hw

if(Lbw)then; Lwm=Lw-1! reflect exterior halo <Lw inwards and add:
   do gw=1,hw; iwm=Lw-gw; iwp=Lwm+gw; r=gw*wLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,kz:nz,iwp)=a(kx:nx,ky:ny,kz:nz,iwp)+r*a(kx:nx,ky:ny,kz:nz,iwm)
      a(kx:nx,ky:ny,kz:nz,iwm)=0
   enddo!                     gw
endif
if(mbw)then; mwp=mw+1! reflect exterior halo >mw inwards and add:
   do gw=1,hw; iwp=mw+gw; iwm=mwp-gw; r=gw*wmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,kz:nz,iwm)=a(kx:nx,ky:ny,kz:nz,iwm)+r*a(kx:nx,ky:ny,kz:nz,iwp)
      a(kx:nx,ky:ny,kz:nz,iwp)=0
   enddo!                     gw
endif

if(Lbz)then; Lzm=Lz-1! reflect exterior halo <Lz inwards and add:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izp,kw:nw)=a(kx:nx,ky:ny,izp,kw:nw)+r*a(kx:nx,ky:ny,izm,kw:nw)
      a(kx:nx,ky:ny,izm,kw:nw)=0
   enddo!                  gz
endif
if(mbz)then; mzp=mz+1! reflect exterior halo >mz inwards and add:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izm,kw:nw)=a(kx:nx,ky:ny,izm,kw:nw)+r*a(kx:nx,ky:ny,izp,kw:nw)
      a(kx:nx,ky:ny,izp,kw:nw)=0
   enddo!                  gz
endif

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp,kz:nz,kw:nw)=a(kx:nx,iyp,kz:nz,kw:nw)+r*a(kx:nx,iym,kz:nz,kw:nw)
      a(kx:nx,iym,kz:nz,kw:nw)=0
   enddo!               gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym,kz:nz,kw:nw)=a(kx:nx,iym,kz:nz,kw:nw)+r*a(kx:nx,iyp,kz:nz,kw:nw)
      a(kx:nx,iyp,kz:nz,kw:nw)=0
   enddo!               gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny,kz:nz,kw:nw)=a(ixp,ky:ny,kz:nz,kw:nw)+r*a(ixm,ky:ny,kz:nz,kw:nw)
      a(ixm,ky:ny,kz:nz,kw:nw)=0
   enddo!            gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny,kz:nz,kw:nw)=a(ixm,ky:ny,kz:nz,kw:nw)+r*a(ixp,ky:ny,kz:nz,kw:nw)
      a(ixp,ky:ny,kz:nz,kw:nw)=0
   enddo!            gx
endif
end subroutine flip4t
!===========================================================================[flipt]
subroutine vflip4t(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw,  xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb,a)
!=================================================================================
integer(spi),intent(in   ):: nv
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz,-hw+Lw:mw+hw),&
             intent(inout):: a
!---------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,gw,ixp,ixm,iyp,iym,izp,izm,iwp,iwm,&
               kx,ky,kz,kw,Lxm,mxp,Lym,myp,Lzm,mzp,Lwm,mwp,nx,ny,nz,nw
!=================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz; kw=Lw-hw; nw=mw+hw

if(Lbw)then; Lwm=Lw-1! reflect exterior halo <Lw inwards and add:
   do gw=1,hw; iwm=Lw-gw; iwp=Lwm+gw; r=gw*wLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,kz:nz,iwp)=a(:,kx:nx,ky:ny,kz:nz,iwp)+r*a(:,kx:nx,ky:ny,kz:nz,iwm)
      a(:,kx:nx,ky:ny,kz:nz,iwm)=0
   enddo!                     gw
endif
if(mbw)then; mwp=mw+1! reflect exterior halo >mw inwards and add:
   do gw=1,hw; iwp=mw+gw; iwm=mwp-gw; r=gw*wmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,kz:nz,iwm)=a(:,kx:nx,ky:ny,kz:nz,iwm)+r*a(:,kx:nx,ky:ny,kz:nz,iwp)
      a(:,kx:nx,ky:ny,kz:nz,iwp)=0
   enddo!                     gw
endif

if(Lbz)then; Lzm=Lz-1! reflect exterior halo <Lz inwards and add:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izp,kw:nw)=a(:,kx:nx,ky:ny,izp,kw:nw)+r*a(:,kx:nx,ky:ny,izm,kw:nw)
      a(:,kx:nx,ky:ny,izm,kw:nw)=0
   enddo!                  gz
endif
if(mbz)then; mzp=mz+1! reflect exterior halo >mz inwards and add:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izm,kw:nw)=a(:,kx:nx,ky:ny,izm,kw:nw)+r*a(:,kx:nx,ky:ny,izp,kw:nw)
      a(:,kx:nx,ky:ny,izp,kw:nw)=0
   enddo!                  gz
endif

if(Lby)then; Lym=Ly-1! reflect exterior halo <Ly inwards and add:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp,kz:nz,kw:nw)=a(:,kx:nx,iyp,kz:nz,kw:nw)+r*a(:,kx:nx,iym,kz:nz,kw:nw)
      a(:,kx:nx,iym,kz:nz,kw:nw)=0
   enddo!               gy
endif
if(mby)then; myp=my+1! reflect exterior halo >my inwards and add:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym,kz:nz,kw:nw)=a(:,kx:nx,iym,kz:nz,kw:nw)+r*a(:,kx:nx,iyp,kz:nz,kw:nw)
      a(:,kx:nx,iyp,kz:nz,kw:nw)=0
   enddo!               gy
endif

if(Lbx)then; Lxm=Lx-1! reflect exterior halo <Lx inwards and add:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny,kz:nz,kw:nw)=a(:,ixp,ky:ny,kz:nz,kw:nw)+r*a(:,ixm,ky:ny,kz:nz,kw:nw)
      a(:,ixm,ky:ny,kz:nz,kw:nw)=0
   enddo!            gx
endif
if(mbx)then; mxp=mx+1! reflect exterior halo >mx inwards and add:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny,kz:nz,kw:nw)=a(:,ixm,ky:ny,kz:nz,kw:nw)+r*a(:,ixp,ky:ny,kz:nz,kw:nw)
      a(:,ixp,ky:ny,kz:nz,kw:nw)=0
   enddo!            gx
endif
end subroutine vflip4t

!# Routine to prepare exterior boundary halos for the application of direct filters:
!==========================================================================[flip]
subroutine flip1(hx,Lx,mx,Lb,mb, xLb,xmb, a)
!===============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                   intent(in   ):: hx,Lx,mx
logical,                        intent(in   ):: Lb,mb
real(dp),                       intent(in   ):: xLb,xmb
real(dp),dimension(Lx-hx:mx+hx),intent(inout):: a
!-------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: ix,ixp,ixm,gx,Lxm,mxp
!===============================================================================
if(Lb)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm)=r*a(ixp)
   enddo
endif
if(mb)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp)=r*a(ixm)
   enddo
endif
end subroutine flip1
!==========================================================================[flip]
subroutine vflip1(nv,hx,Lx,mx,Lb,mb, xLb,xmb, a)
!===============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                      intent(in   ):: nv
integer(spi),                      intent(in   ):: hx,Lx,mx
logical,                           intent(in   ):: Lb,mb
real(dp),                          intent(in   ):: xLb,xmb
real(dp),dimension(nv,Lx-hx:mx+hx),intent(inout):: a
!-------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: ix,ixp,ixm,gx,Lxm,mxp
!===============================================================================
if(Lb)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm)=r*a(:,ixp)
   enddo
endif
if(mb)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp)=r*a(:,ixm)
   enddo
endif
end subroutine vflip1
!=========================================================================[flip]
subroutine flip2(hx,Lx,mx,hy,Ly,my, Lbx,mbx,Lby,mby, xLb,xmb,yLb,ymb, a)
!==============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                                 intent(in   ):: hx,Lx,mx,hy,Ly,my
logical,                                      intent(in   ):: Lbx,mbx,Lby,mby
real(dp),                                     intent(in   ):: xLb,xmb,yLb,ymb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(inout):: a
!------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,ixm,ixp,iym,iyp,kx,ky,Lxm,Lym,mxp,myp,nx,ny
!==============================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xlb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny)=r*a(ixp,ky:ny)
   enddo! gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny)=r*a(ixm,ky:ny)
   enddo! gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym)=r*a(kx:nx,iyp)
   enddo!   gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp)=r*a(kx:nx,iym)
   enddo!   gy
endif
end subroutine flip2

!=========================================================================[flip]
subroutine vflip2(nv, hx,Lx,mx,hy,Ly,my, Lbx,mbx,Lby,mby, xLb,xmb,yLb,ymb, a)
!==============================================================================
use pkind, only: dp,spi
implicit none
integer(spi),                                 intent(in   ):: nv
integer(spi),                                 intent(in   ):: hx,Lx,mx,hy,Ly,my
logical,                                      intent(in   ):: Lbx,mbx,Lby,mby
real(dp),                                     intent(in   ):: xLb,xmb,yLb,ymb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(inout):: a
!------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,ixm,ixp,iym,iyp,kx,ky,Lxm,Lym,mxp,myp,nx,ny
!==============================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xlb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny)=r*a(:,ixp,ky:ny)
   enddo! gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny)=r*a(:,ixm,ky:ny)
   enddo! gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym)=r*a(:,kx:nx,iyp)
   enddo!   gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp)=r*a(:,kx:nx,iym)
   enddo!   gy
endif
end subroutine vflip2
!============================================================================[flip]
subroutine flip3(hx,Lx,mx, hy,Ly,my, hz,Lz,mz,&
     Lbx,mbx,Lby,mby,Lbz,mbz, xLb,xmb,yLb,ymb,zLb,zmb, a)
!==================================================================================
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(inout):: a
!----------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,ixp,ixm,iyp,iym,izp,izm,&
            kx,ky,kz,Lxm,mxp,Lym,myp,Lzm,mzp,nx,ny,nz
!==================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny,kz:nz)=r*a(ixp,ky:ny,kz:nz)
   enddo!    gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny,kz:nz)=r*a(ixm,ky:ny,kz:nz)
   enddo!    gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym,kz:nz)=r*a(kx:nx,iyp,kz:nz)
   enddo!        gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp,kz:nz)=r*a(kx:nx,iym,kz:nz)
   enddo!        gy
endif

if(Lbz)then; Lzm=Lz-1! reflect interior values into the exterior halo <Lz:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izm)=r*a(kx:nx,ky:ny,izp)
   enddo!            gz
endif
if(mbz)then; mzp=mz+1! reflect interior values into the exterior halo >mz:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izp)=r*a(kx:nx,ky:ny,izm)
   enddo!           gz
endif
end subroutine flip3
!============================================================================[flip]
subroutine vflip3(nv, hx,Lx,mx, hy,Ly,my, hz,Lz,mz,&
     Lbx,mbx,Lby,mby,Lbz,mbz, xLb,xmb,yLb,ymb,zLb,zmb, a)
!==================================================================================

integer(spi),intent(in   ):: nv
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(inout):: a
!----------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,ixp,ixm,iyp,iym,izp,izm,&
            kx,ky,kz,Lxm,mxp,Lym,myp,Lzm,mzp,nx,ny,nz
!==================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny,kz:nz)=r*a(:,ixp,ky:ny,kz:nz)
   enddo!    gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny,kz:nz)=r*a(:,ixm,ky:ny,kz:nz)
   enddo!    gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym,kz:nz)=r*a(:,kx:nx,iyp,kz:nz)
   enddo!        gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp,kz:nz)=r*a(:,kx:nx,iym,kz:nz)
   enddo!        gy
endif

if(Lbz)then; Lzm=Lz-1! reflect interior values into the exterior halo <Lz:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izm)=r*a(:,kx:nx,ky:ny,izp)
   enddo!            gz
endif
if(mbz)then; mzp=mz+1! reflect interior values into the exterior halo >mz:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izp)=r*a(:,kx:nx,ky:ny,izm)
   enddo!           gz
endif
end subroutine vflip3
!===========================================================================[flip]
subroutine flip4(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw, xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb, a)
!=================================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz,-hw+Lw:mw+hw),&
             intent(inout):: a
!---------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,gw,ixp,ixm,iyp,iym,izp,izm,iwp,iwm,&
               kx,ky,kz,kw,Lxm,Lym,Lzm,Lwm,mxp,myp,mzp,mwp,nx,ny,nz,nw
!=================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz; kw=Lw-hw; nw=mw+hw

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixm,ky:ny,kz:nz,kw:nw)=r*a(ixp,ky:ny,kz:nz,kw:nw)
   enddo!            gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(ixp,ky:ny,kz:nz,kw:nw)=r*a(ixm,ky:ny,kz:nz,kw:nw)
   enddo!            gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iym,kz:nz,kw:nw)=r*a(kx:nx,iyp,kz:nz,kw:nw)
   enddo!               gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,iyp,kz:nz,kw:nw)=r*a(kx:nx,iym,kz:nz,kw:nw)
   enddo!               gy
endif

if(Lbz)then; Lzm=Lz-1! reflect values into the exterior halo <Lz:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izm,kw:nw)=r*a(kx:nx,ky:ny,izp,kw:nw)
   enddo!                  gz
endif
if(mbz)then; mzp=mz+1! reflect values into the exterior halo >mz:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,izp,kw:nw)=r*a(kx:nx,ky:ny,izm,kw:nw)
   enddo!                  gz
endif

if(Lbw)then; Lwm=Lw-1! reflect values into the exterior halo <Lw:
   do gw=1,hw; iwm=Lw-gw; iwp=Lwm+gw; r=gw*wLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,kz:nz,iwm)=r*a(kx:nx,ky:ny,kz:nz,iwp)
   enddo!                     gw
endif
if(mbw)then; mwp=mw+1! reflect values into the exterior halo >mw:
   do gw=1,hw; iwp=mw+gw; iwm=mwp-gw; r=gw*wmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(kx:nx,ky:ny,kz:nz,iwp)=r*a(kx:nx,ky:ny,kz:nz,iwm)
   enddo!                     gw
endif
end subroutine flip4
!===========================================================================[flip]
subroutine vflip4(nv, hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, &
     Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw, xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb, a)
!=================================================================================

integer(spi),intent(in   ):: nv  
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
logical,     intent(in   ):: Lbx,mbx,Lby,mby,Lbz,mbz,Lbw,mbw
real(dp),    intent(in   ):: xLb,xmb,yLb,ymb,zLb,zmb,wLb,wmb
real(dp),dimension(nv,-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz,-hw+Lw:mw+hw),&
             intent(inout):: a
!---------------------------------------------------------------------------------
real(dp)    :: r
integer(spi):: gx,gy,gz,gw,ixp,ixm,iyp,iym,izp,izm,iwp,iwm,&
               kx,ky,kz,kw,Lxm,Lym,Lzm,Lwm,mxp,myp,mzp,mwp,nx,ny,nz,nw
!=================================================================================
kx=Lx-hx; nx=mx+hx; ky=Ly-hy; ny=my+hy; kz=Lz-hz; nz=mz+hz; kw=Lw-hw; nw=mw+hw

if(Lbx)then; Lxm=Lx-1! reflect values into the exterior halo <Lx:
   do gx=1,hx; ixm=Lx-gx; ixp=Lxm+gx; r=gx*xLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixm,ky:ny,kz:nz,kw:nw)=r*a(:,ixp,ky:ny,kz:nz,kw:nw)
   enddo!            gx
endif
if(mbx)then; mxp=mx+1! reflect values into the exterior halo >mx:
   do gx=1,hx; ixp=mx+gx; ixm=mxp-gx; r=gx*xmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,ixp,ky:ny,kz:nz,kw:nw)=r*a(:,ixm,ky:ny,kz:nz,kw:nw)
   enddo!            gx
endif

if(Lby)then; Lym=Ly-1! reflect values into the exterior halo <Ly:
   do gy=1,hy; iym=Ly-gy; iyp=Lym+gy; r=gy*yLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iym,kz:nz,kw:nw)=r*a(:,kx:nx,iyp,kz:nz,kw:nw)
   enddo!               gy
endif
if(mby)then; myp=my+1! reflect values into the exterior halo >my:
   do gy=1,hy; iyp=my+gy; iym=myp-gy; r=gy*ymb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,iyp,kz:nz,kw:nw)=r*a(:,kx:nx,iym,kz:nz,kw:nw)
   enddo!               gy
endif

if(Lbz)then; Lzm=Lz-1! reflect values into the exterior halo <Lz:
   do gz=1,hz; izm=Lz-gz; izp=Lzm+gz; r=gz*zLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izm,kw:nw)=r*a(:,kx:nx,ky:ny,izp,kw:nw)
   enddo!                  gz
endif
if(mbz)then; mzp=mz+1! reflect values into the exterior halo >mz:
   do gz=1,hz; izp=mz+gz; izm=mzp-gz; r=gz*zmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,izp,kw:nw)=r*a(:,kx:nx,ky:ny,izm,kw:nw)
   enddo!                  gz
endif

if(Lbw)then; Lwm=Lw-1! reflect values into the exterior halo <Lw:
   do gw=1,hw; iwm=Lw-gw; iwp=Lwm+gw; r=gw*wLb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,kz:nz,iwm)=r*a(:,kx:nx,ky:ny,kz:nz,iwp)
   enddo!                     gw
endif
if(mbw)then; mwp=mw+1! reflect values into the exterior halo >mw:
   do gw=1,hw; iwp=mw+gw; iwm=mwp-gw; r=gw*wmb; if(r>u0)r=om0+r*(om1+r*om2)
      a(:,kx:nx,ky:ny,kz:nz,iwp)=r*a(:,kx:nx,ky:ny,kz:nz,iwm)
   enddo!                     gw
endif
end subroutine vflip4

!# Direct beta filters:
!=========================================================================[rbeta]
subroutine rbeta1(hx,Lx,mx, el, a,b)
!===============================================================================
! Perform a direct beta-function filter in 1D.
!
! The input data occupy the extended region:
! Lx-hx <= jx <= mx+hx.
! The output data occupy the central region
! Lx <= ix <= Mx.
!===============================================================================
use pkind, only: dp,spi
use pietc, only: u1
implicit none
integer(spi),                    intent(in   ):: hx,Lx,mx
real(dp),dimension(0:1,Lx:mx   ),intent(in   ):: el
real(dp),dimension(-hx+Lx:mx+hx),intent(in   ):: a
real(dp),dimension(    Lx:mx   ),intent(  out):: b
!-------------------------------------------------------------------------------
real(dp)     :: tb
real(dp)     :: exx,rrc
integer(spi) :: gx,ix,ixp,ixm
!===============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(ix)
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**p*(a(ixp)+a(ixm))
   enddo
   b(ix)=tb*el(0,ix)
enddo
end subroutine rbeta1
!===================================================================[rbeta]
subroutine vrbeta1(nv,hx,lx,mx, el, a,b)
!================================================================================
! Vector version of rbeta1 filtering nv fields at once.
!=============================================================================
integer(spi),                      intent(in   ):: nv,hx,Lx,mx
real(dp),dimension(0:1,Lx:mx),     intent(in   ):: el
real(dp),dimension(nv,lx-hx:mx+hx),intent(in   ):: a
real(dp),dimension(nv, Lx:mx),     intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(nv):: tb
real(dp)              :: exx,rrc
integer(spi)          :: gx,ix,ixp,ixm
!=============================================================================
b=0
do ix=Lx,Mx
   exx=el(1,ix)
   tb=a(:,ix)
   do gx=ceiling(-u1/exx),-1; ixp=ix+gx; ixm=ix-gx
      rrc=u1-(gx*exx)**2
      tb=tb+rrc**p*(a(:,ixp)+a(:,ixm))
   enddo
   b(:,ix)=tb*el(0,ix)
enddo
end subroutine vrbeta1
!========================================================================[rbeta]
subroutine rbeta2(hx,Lx,mx, hy,Ly,my, el, a,b)
!==============================================================================
! Perform a direct radial beta-function filter in 2D.
!==============================================================================
use pkind, only: dp,spi
use pietc, only: u1
implicit none
integer(spi),                                 intent(in   ):: hx,Lx,mx,hy,Ly,my
real(dp),dimension(0:3,Lx:mx,       Ly:my   ),intent(in   ):: el
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy),intent(in   ):: a
real(dp),dimension(    Lx:mx,       Ly:my   ),intent(  out):: b
!------------------------------------------------------------------------------
real(dp),dimension(3):: tel
real(dp)             :: tb
real(dp)             :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)         :: gx,gy,ix,ixm,ixp,iy,iym,iyp
!==============================================================================
b=0
do iy=Ly,my; do ix=Lx,mx
   tel=el(1:3,ix,iy) 
   exx=tel(1); eyy=tel(2); eyx=tel(3)
   tb=a(ix,iy)
lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx);ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         tb=tb+rrc**p*(a(ixp,iyp)+a(ixm,iym))
      enddo! gx
   enddo lgy
   b(ix,iy)=tb*el(0,ix,iy)
enddo; enddo! ix, iy
end subroutine rbeta2
!====================================================================== [rbeta]
subroutine vrbeta2(nv,hx,lx,mx, hy,ly,my, el, a,b)
!=============================================================================
! Vector version of rbeta2 filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my
real(dp),dimension(0:3,Lx:Mx,Ly:My),           intent(in   ):: el
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy),intent(in   ):: a
real(dp),dimension(nv, lx:mx,ly:my),           intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(3) :: tel
real(dp),dimension(nv):: tb
real(dp)              :: cx,exx,eyy,eyx,r,rrc,rrxc
integer(spi)          :: gx,gy,ix,ixp,ixm,iy,iyp,iym
!=============================================================================
b=0
do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:3,ix,iy)
   exx=tel(1);eyy=tel(2);eyx=tel(3)
   tb=a(:,ix,iy)
lgy: do gy=ceiling(-u1/eyy),0; iyp=iy+gy; iym=iy-gy
      rrxc=abs(u1-(gy*eyy)**2); r=sqrt(rrxc); cx=-gy*eyx
      do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
         if(gy==0.and.gx==0)exit lgy
         rrc=rrxc-(gx*exx-cx)**2
         tb=tb+rrc**p*(a(:,ixp,iyp)+a(:,ixm,iym))
      enddo! gx
   enddo lgy
   b(:,ix,iy)=tb*el(0,ix,iy)
enddo;   enddo! ix, iy
end subroutine vrbeta2
!===========================================================================[rbeta]
subroutine rbeta3(hx,Lx,mx, hy,Ly,my, hz,Lz,mz, el,a,b)
!==================================================================================
! Perform a direct radial beta-function filter in 3D.
!==================================================================================
integer(spi),intent(in   ):: hx,Lx,mx, hy,ly,my, hz,lz,mz
real(dp),dimension(0:6,Lx:mx,       Ly:my,       Lz:mz),   intent(in   ):: el
real(dp),dimension(-hx+Lx:mx+hx,-hy+Ly:my+hy,-hz+Lz:mz+hz),intent(in   ):: a
real(dp),dimension(    Lx:mx,       Ly:my,       Lz:mz),   intent(  out):: b
!----------------------------------------------------------------------------------
real(dp),dimension(6):: tel
real(dp)             :: tb
real(dp)             :: cx,cy,czx,exx,eyy,ezz,eyx,ezx,ezy,r,rrc,rrxc,rryc
integer(spi)         :: gx,gy,gz,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm
!==================================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:6,ix,iy,iz)
   exx=tel(1);eyy=tel(2);ezz=tel(3);ezy=tel(4);ezx=tel(5);eyx=tel(6)
   tb=a(ix,iy,iz)
   lgz: do gz=ceiling(-u1/ezz),0; izp=iz+gz; izm=iz-gz
      rryc=abs(u1-(gz*ezz)**2); r=sqrt(rryc); cy=-gz*ezy; czx=-gz*ezx
      do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
         rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
         do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
            if(gz==0.and.gy==0.and.gx==0)exit lgz
            rrc=rrxc-(gx*exx-cx)**2
            tb=tb+rrc**p*(a(ixp,iyp,izp)+a(ixm,iym,izm))
         enddo! gx
      enddo! gy
   enddo lgz
   b(ix,iy,iz)=tb*el(0,ix,iy,iz)
enddo;   enddo;    enddo! ix, iy, iz
end subroutine rbeta3
!=======================================================================[rbeta]
subroutine vrbeta3(nv, hx,lx,mx, hy,ly,my, hz,lz,mz, el,a,b)
!=============================================================================
! Vector version of rbeta3 filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my,hz,lz,mz
real(dp),dimension(0:6,Lx:mx,Ly:my,Lz:mz),intent(in   ):: el
real(dp),dimension(nv,Lx-hx:mx+hx,Ly-hy:my+hy,Lz-hz:mz+hz),&
                                          intent(in   ):: a
real(dp),dimension(nv, Lx:mx,Ly:my,Lz:mz),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(6) :: tel
real(dp),dimension(nv):: tb
real(dp)              :: cx,cy,czx,exx,eyy,ezz,eyx,ezx,ezy,r,rrc,rrxc,rryc
integer(spi)          :: gx,gy,gz,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm
!=============================================================================
b=0
do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:6,ix,iy,iz)
   exx=tel(1);eyy=tel(2);ezz=tel(3);ezy=tel(4);ezx=tel(5);eyx=tel(6)
   tb=a(:,ix,iy,iz)
lgz: do gz=ceiling(-u1/ezz),0; izp=iz+gz; izm=iz-gz
      rryc=abs(u1-(gz*ezz)**2); r=sqrt(rryc); cy=-gz*ezy; czx=-gz*ezx
      do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
         rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
         do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
            if(gz==0.and.gy==0.and.gx==0)exit lgz
            rrc=rrxc-(gx*exx-cx)**2
            tb=tb+rrc**p*(a(:,ixp,iyp,izp)+a(:,ixm,iym,izm))
         enddo! gx
      enddo! gy
   enddo lgz
   b(:,ix,iy,iz)=tb*el(0,ix,iy,iz)
enddo;   enddo;    enddo! ix, iy, iz
end subroutine vrbeta3
!===========================================================================[rbeta]
subroutine rbeta4(hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,a,b)
!=================================================================================
! Perform a direct radial beta-function filter in 4D.
!=================================================================================
integer(spi),intent(in   ):: hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
real(dp),dimension(0:10,Lx:Mx,Ly:My,Lz:Mz,Lw:Mw),intent(in   ):: el
real(dp),dimension(lx-hx:mx+hx,ly-hy:my+hy,lz-hz:mz+hz,lw-hw:mw+hw),&
                                                 intent(in   ):: a
real(dp),dimension(     lx:mx,ly:my,lz:mz,lw:mw),intent(  out):: b
!---------------------------------------------------------------------------------
real(dp),dimension(10):: tel
real(dp)              :: tb
real(dp)              :: cx,cy,cz,cwy,cwx,czx,exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,&
                         r,rrc,rrxc,rryc,rrzc
integer(spi)          :: gx,gy,gz,gw,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm,iw,iwp,iwm
!=================================================================================
b=0
do iw=lw,mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:10,ix,iy,iz,iw)
   exx=tel(1);eyy=tel(2);ezz=tel(3);eww=tel(4);eyx=tel(5)
   ezx=tel(6);ewx=tel(7);ezy=tel(8);ewy=tel(9);ewz=tel(10)
   tb=a(ix,iy,iz,iw)
   lgw: do gw=ceiling(-u1/eww),0; iwp=iw+gw; iwm=iw-gw
      rrzc=abs(u1-(gw*eww)**2); r=sqrt(rrzc); cz=-gw*ewz; cwy=-gw*ewy; cwx=-gw*ewx
      do gz=ceiling((cz-r)/ezz),floor((cz+r)/ezz); izp=iz+gz; izm=iz-gz
         rryc=abs(rrzc-(gz*ezz-cz)**2); r=sqrt(rryc); cy=cwy-gz*ezy; czx=cwx-gz*ezx
         do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
            rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
            do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
               if(gw==0.and.gz==0.and.gy==0.and.gx==0)exit lgw
               rrc=rrxc-(gx*exx-cx)**2
               tb=tb+rrc**p*(a(ixp,iyp,izp,iwp)+a(ixm,iym,izm,iwm))
            enddo! gx
         enddo! gy
      enddo! gz
   enddo lgw
   b(ix,iy,iz,iw)=tb*el(0,ix,iy,iz,iw)
enddo;   enddo;   enddo;   enddo! ix, iy, iz, iw
end subroutine rbeta4
!=======================================================================[rbeta]
subroutine vrbeta4(nv,hx,lx,mx, hy,ly,my, hz,lz,mz, hw,lw,mw, el,a,b)
!=============================================================================
! Vector version of rbeta4 filtering nv fields at once.
!=============================================================================
integer(spi),intent(in   ):: nv,hx,Lx,mx,hy,ly,my,hz,lz,mz,hw,lw,mw
real(dp),dimension(0:10,Lx:mx,Ly:my,Lz:mz,Lw:mw),intent(in   ):: el
real(dp),dimension(nv,lx-hx:mx+hx,ly-hy:my+hy,&
     lz-hz:mz+hz,lw-hw:mw+hw),                   intent(in   ):: a
real(dp),dimension(nv,  Lx:mx,Ly:my,Lz:mz,Lw:mw),intent(  out):: b
!-----------------------------------------------------------------------------
real(dp),dimension(10):: tel
real(dp),dimension(nv):: tb
real(dp)              :: cx,cy,cz,cwy,cwx,czx,exx,eyy,ezz,eww,eyx,ezx,ewx,ezy,ewy,ewz,&
                         r,rrc,rrxc,rryc,rrzc
integer(spi)          :: gx,gy,gz,gw,ix,ixp,ixm,iy,iyp,iym,iz,izp,izm,iw,iwp,iwm
!=============================================================================
b=0
do iw=lw,mw; do iz=Lz,Mz; do iy=Ly,My; do ix=Lx,Mx
   tel=el(1:10,ix,iy,iz,iw)
   exx=tel(1);eyy=tel(2);ezz=tel(3);eww=tel(4);eyx=tel(5)
   ezx=tel(6);ewx=tel(7);ezy=tel(8);ewy=tel(9);ewz=tel(10)
   tb=a(:,ix,iy,iz,iw)
lgw: do gw=ceiling(-u1/eww),0; iwp=iw+gw; iwm=iw-gw
      rrzc=abs(u1-(gw*eww)**2); r=sqrt(rrzc); cz=-gw*ewz; cwy=-gw*ewy; cwx=-gw*ewx
      do gz=ceiling((cz-r)/ezz),floor((cz+r)/ezz); izp=iz+gz; izm=iz-gz
         rryc=abs(rrzc-(gz*ezz-cz)**2); r=sqrt(rryc); cy=cwy-gz*ezy; czx=cwx-gz*ezx
         do gy=ceiling((cy-r)/eyy),floor((cy+r)/eyy); iyp=iy+gy; iym=iy-gy
            rrxc=abs(rryc-(gy*eyy-cy)**2); r=sqrt(rrxc); cx=czx-gy*eyx
            do gx=ceiling((cx-r)/exx),floor((cx+r)/exx); ixp=ix+gx; ixm=ix-gx
               if(gw==0.and.gz==0.and.gy==0.and.gx==0)exit lgw
               rrc=rrxc-(gx*exx-cx)**2
               tb=tb+rrc**p*(a(:,ixp,iyp,izp,iwp)+a(:,ixm,iym,izm,iwm))
            enddo! gx
         enddo! gy
      enddo! gz
   enddo lgw
   b(:,ix,iy,iz,iw)=tb*el(0,ix,iy,iz,iw)
enddo;  enddo;  enddo;  enddo! ix, iy, iz, iw
end subroutine vrbeta4

end module ybfil
!#
