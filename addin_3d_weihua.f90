!      nf: the # of real values for an expansion of fictitious value, 
!          real value may repeat
!      ns: the # of real values for distretization of poisson equation 
!          at one point. tentative value ns=350 for 4th order, ns=100 for 2nd order
!      ls: the # of points used for derivatives to find fictitious points.
!      ms: the # of points used for derivatves in governing equations.
!      iord: define the order to search auxiliary points on lines that 
!            parallel to main direction
!            iord(:,1) is for interface between i and i+1,
!            iord(:,-1) is for interface between i-1 and i. where i 
!                       could be ix,jy or kz when interface is passes 
!                       through i, only iord(1:ls,-1) is needed
!      iside(nx*ny*nz) : indicator of side of point 1 inside, 0 oninterface, -1 outside
!      irr(nx*ny*nz): indicator of irregular 1 or irregular point 0
!      idre(3,7): the (i,j,k)s of the seven(x+/-, y+/-, z+/-, center) points 
!      idms(2,3*(ms-1)): the (i,j,k)s of surrounding points for (ms-1)th order 
!                         difference scheme
!      icolu(ns,nx,ny,nz) : for each point A(ix,jy,kz), store the (is-1)*ny*nz
!                           +(js-1)*nz+ks of ns points(is,js,ks) envolved in 
!                           point A's diffenrence scheme
!      coef(ns,nx,ny,nz) : together with icolu(ns,nx,ny,nz) and provides the 
!                           coefficient of point(is,js,ks)
!      licolu(nx,ny,nz) : the real total number of (is,js,ks) points for
!                         point A's difference schem, it should less or equal ns
!      cij(nx,ny,nz)     : added up all terms with jumpconditions for point 
!                          A's difference scheme.
!      iord(ls+1,-1:1) : define the order to search auxiliary points, 
!                        iord(:,1) is for interface between i and i+1, 
!                        iord(:,-1) is for interface between i-1 and i. 
!                        where i could be ix,jy or kz
!      idft(nx,ny,nz) :  the index of an irregular point in the array of all 
!                        irregular point, the value of idft might be multi-
!                        defined, which means the defined value are overwritten
!                        if the irregular point intersects interface
!                        more than once.
!      numir : the total number of fictitious values needed, counted by paired 
!              irregular points > the number of irregular point
!      nmax  : the total number of nonzero entry in bicg matrix
!      icoluf(numir,:,nf) : for each irregular point A(ir,jr,kr), stores the 
!                           (if-1)*ny*nz+(jf-1)*nz+kf of points(if,jf,kf), which
!                           are used for one of the expansions of the fictitious  
!                           values of point A. some parts of the icoluf are not 
!                           used because the multi-defining of idft
!      icoluf(:,1:7,nf) : the expansion can be obtained from 6 different direction.
!                         the direction are defined in idre(3,7). direction 2 is void
!                         which means the center in idre(3,7).
!      akf(numir,:,nf) : together with icoluf(numir,7,nf) and provide the coefficients
!                        of each point(if,jf,kf), similarly, some parts of the akf
!                        are not used because the multi-defining of idft
!      ija(nmax) : the index of nonzero entry of bicg matrix
!      sa(nmax) :  the nonzero entry of bicg matrix

!********************************************************************************
!addin.f90
!**********************************************************************************
!idre(1:3,1)=(/-1,0,0/)
!idre(1:3,2)=(/0,0,0/)
!idre(1:3,3)=(/1,0,0/)
!idre(1:3,4)=(/0,-1,0/)
!idre(1:3,5)=(/0,1,0/)
!idre(1:3,6)=(/0,0,-1/)
!idre(1:3,7)=(/0,0,1/)

!
!-------------------------------------------------------------------------------
!
subroutine settlecoefrp(ix,jy,kz,inext)
!-----------------------------------------------
! four points 1,2,9,10 and geometry are shown below
!             \
!              \  
! --*---9---10--\--1---2---*---*------------
!               |
!     inside   /  outside
!-----------------------------------------------------
use pbedata
use comdata
implicit double precision(a-h,o-z)
dimension infopt(19,2),ak(22,2),iun(4)

!print *,ix,jy,kz,inext, 'patch needed'
idx=idre(1,inext)
idy=idre(2,inext)
idz=idre(3,inext)

! decide whether fictious points f_10, f_1 and f_9/_2 available
iun=0 
!if(sum(iftpts(indirr(ix+idx,jy+idy,kz+idz),:,1))==0)  iun(1)=1  ! u_1^+ is unavailable
if(sum(abs(ftpts(indirr(ix+idx,jy+idy,kz+idz),:,:)))==0)  iun(1)=1  ! u_1^+ is unavailable

if(io(ix+2*idx,jy+2*idy,kz+2*idz)==-1)then
!    if(sum(iftpts(indirr(ix+2*idx,jy+2*idy,kz+2*idz),:,1))==0) iun(2)=1 
	iun(2)=1
endif
! u_2^- is unavailable, if u_2 is on grid, then it is available, but correction is needed.
if(io(ix-idx,jy-idy,kz-idz)==1)then
!    if(sum(iftpts(indirr(ix-idx,jy-idy,kz-idz),:,1))==0)  iun(3)=1  ! u_9^+ is unavailable
	iun(3)=1
endif

!if(sum(iftpts(indirr(ix,jy,kz),:,1))==0)    iun(4)=1 !        u_10^- is unavailable
if(sum(abs(ftpts(indirr(ix,jy,kz),:,:)))==0)    iun(4)=1 !        u_10^- is unavailable

if(iun(1)==0.and.iun(4)==0) then 
	!print *,ix,jy,kz,inext, 'returned'
	return
endif

!if (iun(2)+iun(3)==1 .and. sum(iun(:))<=3) then
!    print *,ix,jy,kz,inext,'patch 3 needed'
!    call replace3(ix,jy,kz,inext,iun,lpp,lnn)
    
!else
if(sum(iun)>=3)then              ! if three or more unknown
	itnew=itnew+1
	itpsnew(itnew,1:4)=(/ix,jy,kz,inext/)
	itpsnew(itnew,5:8)=iun(:)
	!print *,ix,jy,kz,inext,'three unknows'
elseif(sum(iun)==1)then          ! elseif two are available, one is unknown
	!print *,ix,jy,kz,inext,' patch1 needed'
	call replace1(ix,jy,kz,inext,iun)
elseif(sum(iun)==2)then			! elseif one are available, two are unknown
	!print *,ix,jy,kz,inext,'patch2 needed'
	call replace2(ix,jy,kz,inext,iun)
endif

!print *,ix*10000+jy*100+kz,inext,iun

return
end
!
!------------------------------------------------------------------------
!
subroutine replace1(ix,jy,kz,inext,iun)
use comdata
use pbedata
!---------------------------------------------------------------
!iuk1 tells which fictious point unknown
!ak1(:,1) the coefficience for f^- 
!ak1(:,2) the coefficience for f^+ 
!ak1(1:13,:) the coefficience of (u1-u10,w,v,dw) respectively
!ak(:,1) the coefficience for f_1 as output
!ak(:,2) the coefficience for f_9 as output
!au(2,11) the coefficience A of u_i, i.ne.2,9 and w,v,dw
!bu(2,4) the coefficience B of f_1,f_2,f_3(f_9),f_4(f_10)
!buk(2,2) the coefficience of unknown fictious points among f_1,f_2,f_3,f_4
!bkn(2,2)  the coefficience of known fictious points among f_1,f_2,f_3,f_4
!bukin(2,2) the inverse of buk(2,2)
!-------------------------------------------------------------
implicit double precision(a-h,o-z)
real*8 bu(4),au(3),weitn(0:2,0:1),weitp(0:2,0:1),akwv(3),x3p(0:2),ak4(nf)
integer ifour(4),ijkn(3),ikn(3),ixkn(3),jykn(3),kzkn(3),iun(4),iex(4),icolu4(nf)

akwv=0.d0
ikn(1)=0 ; ikn(2)=0
x3p=(/-dx,0.d0,dx/)
do i=1,4
    if(iun(i)==0)then
	    if(ikn(1)==0)then
		    ikn(1)=i
		elseif(ikn(2)==0)then
		    ikn(2)=i
		else
		    ikn(3)=i
		endif
	else
	    iuk1=i
	endif
enddo

xx=x(ix)
yy=y(jy)
zz=z(kz)

idx=idre(1,inext)
idy=idre(2,inext)
idz=idre(3,inext)

in=ix+idx
jn=jy+idy
kn=kz+idz

ifour(1)=id3d(in,jn,kn)              ! f_1
ifour(2)=id3d(ix+2*idx,jy+2*idy,kz+2*idz)  ! f_2
ifour(3)=id3d(ix-idx,jy-idy,kz-idz)      ! f_9
ifour(4)=id3d(ix,jy,kz)              ! f_10

iex(1)=7-inext;	iex(2)=2  ! 2 here means from useless direction
iex(3)=2;		iex(4)=inext

! to find u_o on the interface
t=dist_angle(5,inext,indirr(ix,jy,kz))

xo=xx+t*abs(idx)  ; yo=yy+t*abs(idy) ; zo=zz+t*abs(idz)

if (ibd==0) then
    weta=-uj(xo,yo,zo)                 ! weta=ujump
else
    if (icg==2) then
        weta=phi_star(xo,yo,zo)        ! jumps in mibpb_newReg
    else
        weta=0.d0    ! ujump=0.d0 in PBE (mibpb2,mibpb3)
    endif
endif

! the weight of point 9,10 and f^-
call weights(abs(t),x3p,2,2,1,weitp)        ! the order of weight is from inside to outside
! the weight of point f^+ and 1,2
call weights(abs(t)-dx,x3p,2,2,1,weitn)  ! the order of weight is from inside to outside


bu(1:4)=(/weitp(2,0),-weitn(2,0),weitp(0,0),-weitn(0,0)/)
au(1:3)=(/weitn(1,0),-weitp(1,0),weta/)
bi=bu(iuk1)

ixyz4=ifour(iuk1)
call line2cube(ix4,jy4,kz4,ixyz4)

ak4(:)=0.d0
icolu4(:)=-1
ak4(1:2)=(/au(1)/bi,au(2)/bi/)
icolu4(1)=id3d(ix+idx,jy+idy,kz+idz)    ! u_1^-
icolu4(2)=id3d(ix,jy,kz)            ! u_10^+
ipt=2

do i=1,3
	ijkn(i)=ifour(ikn(i))
	call line2cube(ic,jc,kc,ijkn(i))
	if(io(ic,jc,kc)==(-1)**(ikn(i)).or.io(ic,jc,kc)==0)then  ! known point has real value,

                                              ! don't need fiticious value
		ipt=ipt+1
		ak4(ipt)=-bu(ikn(i))/bi             ! the coefficient of known is stored
		icolu4(ipt)=ifour(ikn(i))       ! the location of known is stored 
	    akwv(i)=0.d0

	else
		call line2cube(ixkn(i),jykn(i),kzkn(i),ifour(ikn(i)))
		if(iftpts(indirr(ixkn(i),jykn(i),kzkn(i)),iex(ikn(i)),1)>0)then   ! check fictitious point from which direction
			iexch=iex(ikn(i))
		else
			jj=1 
			do while(jj<=6.and.iftpts(indirr(ixkn(i),jykn(i),kzkn(i)),jj,1)==0)
				jj=jj+1
			enddo
			iexch=jj
		endif
	 	do j=1,16
			ipt=ipt+1
			ii=iftpts(indirr(ixkn(i),jykn(i),kzkn(i)),iexch,3*j-2)
			jj=iftpts(indirr(ixkn(i),jykn(i),kzkn(i)),iexch,3*j-1)
			kk=iftpts(indirr(ixkn(i),jykn(i),kzkn(i)),iexch,3*j)
			icolu4(ipt)=id3d(ii,jj,kk)
			ak4(ipt)=-bu(ikn(i))/bi*ftpts(indirr(ixkn(i),jykn(i),kzkn(i)),iexch,j)
		enddo
		akwv(i)=akwv(i)-bu(ikn(i))/bi*ftc(indirr(ixkn(i),jykn(i),kzkn(i)),iexch)
	endif
enddo
 
ipt=ipt+1
icolu4(ipt)=0
ak4(ipt)=au(3)/bi+sum(akwv(1:3))   ! the sum of w,v,dw terms for unknown

!---------------------------------------------------
ixn2=ix+2*idre(1,inext)
jyn2=jy+2*idre(2,inext)
kzn2=kz+2*idre(3,inext)
if(iun(2)==0.and.io(ixn2,jyn2,kzn2)==0)then
    xn2=x(ixn2)
	yn2=y(jyn2)
	zn2=z(kzn2)
    if(ikn(1)==2)then
	    ichange=1
	elseif(ikn(2)==2)then
	    ichange=2
	else
	    ichange=3
	endif
    if (ibd==0) then
        ak4(ipt)=ak4(ipt)+uj(xn2,yn2,zn2)! ujump=0.d0 in PBE  
    else 
        if (icg==2) then
            ak4(ipt)=ak4(ipt)+phi_star(xn2,yn2,zn2)! ujump=phi_star in mibpb_newReg  
        endif
    endif
endif
!--------------------------------------------

if(iex(iuk1)==7)then
    ak4(:)=0.d0
	icolu4(:)=-1
endif

call line2cube(ic,jc,kc,ixyz4)
if(irrpts(ic,jc,kc)==1)then
    iph=iph+1
	idpatch(iph)=id3d(ix4,jy4,kz4)
	ikk=iex(iuk1)
	ikk=7-ikk
	iftpts(indirr(ix4,jy4,kz4),ikk,1)=-1
	!print *,ix4,jy4,kz4,ikk,'fixed by patch1'
	icpatch(iph,ikk,:)=icolu4(:)
	akpatch(iph,ikk,1:ipt-1)=ak4(1:ipt-1)
	pchftc(iph,ikk)=ak4(ipt)
endif

end
!
!------------------------------------------------------------------------
!
subroutine replace2(ix,jy,kz,inext,iun)
use pbedata
use comdata
!---------------------------------------------------------------
!iuk1, iuk2 tells which two fictious points unknown
!ak1(:,1) the coefficience for f^- 
!ak1(:,2) the coefficience for f^+ 
!ak1(1:13,:) the coefficience of (u1-u10,w,v,dw) respectively
!ak(:,1) the coefficience for f_1 as output
!ak(:,2) the coefficience for f_9 as output
!au(2,11) the coefficience A of u_i, i.ne.2,9 and w,v,dw
!bu(2,4) the coefficience B of f_1,f_2,f_3(f_9),f_4(f_10)
!buk(2,2) the coefficience of unknown fictious points among f_1,f_2,f_3,f_4
!bkn(2,2)  the coefficience of known fictious points among f_1,f_2,f_3,f_4
!bukin(2,2) the inverse of buk(2,2)
!-------------------------------------------------------------
implicit double precision(a-h,o-z)
!dimension infopt1(10,2),ak1(13,2)
dimension   au(2,15),bu(2,4),ikn(2),ijkn(2),ixkn(2),jykn(2),kzkn(2)&       
            ,buk(2,2),bkn(2,2),bukin(2,2),ifour(4),iun(4),iex(4),akwv(2,2)&
			,icolu3(nf),ak3(nf),icolu4(nf),ak4(nf),ictemp(nf),aktemp(nf)&
			,aktempn(nf)
integer iy2g(4),iy2gg(4)


akwv=0.d0
ikn(1)=0 ; iuk1=0
do i=1,4
    if(iun(i)==0)then
	    if(ikn(1)==0)then
		    ikn(1)=i
		else
		    ikn(2)=i
		endif
	else
	    if(iuk1==0)then
		    iuk1=i
		else
		    iuk2=i
		endif
	endif
enddo

i_trouble=itrbl

iforce=1
call findfctin(ix,jy,kz,inext,iforce)


if(i_trouble.ne.itrbl .or. iforce==2)then 
    !print *,'sub_patch2',ix,jy,kz,inext
    !call substitution(ix,jy,kz,inext)
    !return
	itnew=itnew+1
	itpsnew(itnew,1:4)=(/ix,jy,kz,inext/)
	itpsnew(itnew,5:8)=iun(:)
return	
!   print *,'patch2 doesnt work',ix,jy,kz,inext
endif

in=ix+idre(1,inext)
jn=jy+idre(2,inext)
kn=kz+idre(3,inext)

iy2g=(/3,4,1,2/)
do ii=1,4
i0=iftpts(indirr(ix,jy,kz),7-inext,3*iy2g(ii)-2)
j0=iftpts(indirr(ix,jy,kz),7-inext,3*iy2g(ii)-1)
k0=iftpts(indirr(ix,jy,kz),7-inext,3*iy2g(ii))
ifour(ii)=id3d(i0,j0,k0)   

enddo

! the direction of fictitious values of point 1 and 10 is given
iex(1)=7-inext;		iex(2)=7  ! 2 here means from useless direction, 7 for michael's code
iex(3)=7;			iex(4)=inext


aktemp(1:16)=ftpts(indirr(ix,jy,kz),7-inext,1:16)
aktemp(17)=ftc(indirr(ix,jy,kz),7-inext)
aktempn(1:16)=ftpts(indirr(in,jn,kn),inext,1:16)
aktempn(17)=ftc(indirr(in,jn,kn),inext)


do iii=1,16
i0=iftpts(indirr(ix,jy,kz),7-inext,iii*3-2)
j0=iftpts(indirr(ix,jy,kz),7-inext,iii*3-1)
k0=iftpts(indirr(ix,jy,kz),7-inext,iii*3)
ictemp(iii)=id3d(i0,j0,k0)
enddo

!print *,iftpts(indirr(ix,jy,kz),7-inext,:)-iftpts(indirr(in,jn,kn),inext,:)

!test$$$$$$$$$$$$$$$$$$$$
	utest=0.d0
	ixt=in;	jyt=jn;	kzt=kn;	inet=inext
	!print *,ixt,jyt,kzt,inet,io(ixt,jyt,kzt)
	do i=1,16
	   i0=iftpts(indirr(ixt,jyt,kzt),inet,i*3-2)
	   j0=iftpts(indirr(ixt,jyt,kzt),inet,i*3-1)
	   k0=iftpts(indirr(ixt,jyt,kzt),inet,i*3)
	   !print *,i,i0,j0,k0,ftpts(indirr(ixt,jyt,kzt),inet,i)  
	   ijk=id3d(i0,j0,k0)
	   if (i0==29 .and. j0==18 .and. k0==47) then
	   	    utest=utest+un(x(i0),y(j0),z(k0))*ftpts(indirr(ixt,jyt,kzt),inet,i)   
	   else
			utest=utest+uexa(ijk)*ftpts(indirr(ixt,jyt,kzt),inet,i)   
	   endif
	enddo
	   utest=utest+ftc(indirr(ixt,jyt,kzt),inet)
	   
	   
	   !???????????????????????????????????????????????????
	   if (io(ixt,jyt,kzt)==0 .or. io(ixt,jyt,kzt)==-1) then
		  uexact=up(x(ixt),y(jyt),z(kzt))
	   else
	   	  uexact=un(x(ixt),y(jyt),z(kzt))
	   endif
	   !print *,'fic value of point',real(uexact),real(utest),abs(real(uexact-utest))
!--------------------------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$

ftpts(indirr(ix,jy,kz),7-inext,:)=0.d0
ftc(indirr(ix,jy,kz),7-inext)=0.d0
iftpts(indirr(ix,jy,kz),7-inext,:)=0

ftpts(indirr(in,jn,kn),inext,:)=0.d0
ftc(indirr(in,jn,kn),inext)=0.d0
iftpts(indirr(in,jn,kn),inext,:)=0

iy2gg=(/2,1,3,4/)
au(1,1)=aktempn(iy2gg(1)) 
au(2,1)=aktemp(iy2gg(1))
au(1,2)=aktempn(iy2gg(3))
au(2,2)=aktemp(iy2gg(3))
au(1,3:15)=aktempn(5:17) 
au(2,3:15)=aktemp(5:17)
bu(1,:)=(/1.d0,-aktempn(iy2gg(4)),-aktempn(iy2gg(2)),0.d0/)
bu(2,:)=(/0.d0,-aktemp(iy2gg(4)),-aktemp(iy2gg(2)),1.d0/)

buk(1,:)=(/bu(1,iuk1),bu(1,iuk2)/)
buk(2,:)=(/bu(2,iuk1),bu(2,iuk2)/)
bkn(1,:)=(/bu(1,ikn(1)),bu(1,ikn(2))/)
bkn(2,:)=(/bu(2,ikn(1)),bu(2,ikn(2))/)
dbuk=buk(1,1)*buk(2,2)-buk(1,2)*buk(2,1)
bukin(1,:)=(/buk(2,2)/dbuk,-buk(1,2)/dbuk/)
bukin(2,:)=(/-buk(2,1)/dbuk,buk(1,1)/dbuk/)
au=matmul(bukin,au)
bkn=-matmul(bukin,bkn)

ixyz3=ifour(iuk1)
call line2cube(ix3,jy3,kz3,ixyz3)
ixyz4=ifour(iuk2)
call line2cube(ix4,jy4,kz4,ixyz4)

icolu3(:)=0 ; icolu4(:)=0
ak3(:)=0.d0  ; ak4(:)=0.d0

icolu3(1)=ictemp(iy2gg(1))  
icolu3(2)=ictemp(iy2gg(3))
icolu3(3:14)=ictemp(5:16)
ak3(1:14)=au(1,1:14)
ak4(1:14)=au(2,1:14)

ipt=14

do j=1,2
	ijkn(j)=ifour(ikn(j))
	call line2cube(ic,jc,kc,ijkn(j))
    if(io(ic,jc,kc)==(-1)**(ikn(j)).or.io(ic,jc,kc)==0)then  ! known point #1 has real value,
	                                     ! don't need fiticious value
		ipt=ipt+1
		ak3(ipt)=bkn(1,j)                ! the coefficient of known #1 is stored
		ak4(ipt)=bkn(2,j)
    	icolu3(ipt)=ifour(ikn(j))       ! the location of known#1 is stored 
		akwv(1,j)=0.d0
		akwv(2,j)=0.d0
	else
		if(abs(bkn(1,j))>100.d0.or.abs(bkn(2,j))>10)then
			itnew=itnew+1
			itpsnew(itnew,1:4)=(/ix,jy,kz,inext/)
			itpsnew(itnew,5:8)=iun(:)
			!print *,ix,jy,kz,inext,' pch2 err is big, use patch3'
			return
		endif
		
		call line2cube(ixkn(j),jykn(j),kzkn(j),ifour(ikn(j)))
		if(iftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iex(ikn(j)),1)>0)then
			iexch=iex(ikn(j))
		else
			jj=1 
			do while(jj<=6.and.iftpts(indirr(ixkn(j),jykn(j),kzkn(j)),jj,1)==0)
				jj=jj+1
			enddo
			iexch=jj
		endif
		do i=1,16
			ipt=ipt+1
			i0=iftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iexch,3*i-2)
			j0=iftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iexch,3*i-1)
			k0=iftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iexch,3*i)
			icolu3(ipt)=id3d(i0,j0,k0)
			ak3(ipt)=bkn(1,j)*ftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iexch,i)
			ak4(ipt)=bkn(2,j)*ftpts(indirr(ixkn(j),jykn(j),kzkn(j)),iexch,i)
		enddo
			akwv(1,j)=akwv(1,j)+bkn(1,j)*ftc(indirr(ixkn(j),jykn(j),kzkn(j)),iexch)
			akwv(2,j)=akwv(2,j)+bkn(2,j)*ftc(indirr(ixkn(j),jykn(j),kzkn(j)),iexch)
	endif
enddo

ipt=ipt+1
icolu3(ipt)=0
ak3(ipt)=au(1,15)+akwv(1,1)+akwv(1,2)   ! the sum of w,v,dw terms for unknown#1
ak4(ipt)=au(2,15)+akwv(2,1)+akwv(2,2)   ! the sum of w,v,dw terms for unknown#2

icolu4(1:ipt)=icolu3(1:ipt)

!---------------------------------------------------
ixn2=ix+2*idre(1,inext)
jyn2=jy+2*idre(2,inext)
kzn2=kz+2*idre(3,inext)

if(iun(2)==0.and.io(ixn2,jyn2,kzn2)==0)then
    xn2=x(ixn2)
	yn2=y(jyn2)
	zn2=z(kzn2)
    if(ikn(1)==2)then
	    ichange=1
	else
	    ichange=2
	endif
    if (ibd==0) then
        ak3(ipt)=ak3(ipt)+bkn(1,ichange)*uj(xn2,yn2,zn2)
	    ak4(ipt)=ak4(ipt)+bkn(2,ichange)*uj(xn2,yn2,zn2)    ! ujump=0.d0 in PBE
    else
        if (icg==2) then !mibpb3_newReg
            ak3(ipt)=ak3(ipt)+bkn(1,ichange)*phi_star(xn2,yn2,zn2)
	        ak4(ipt)=ak4(ipt)+bkn(2,ichange)*phi_star(xn2,yn2,zn2)    
        endif
    endif
endif

call line2cube(ic,jc,kc,ixyz3)
if(irrpts(ic,jc,kc)==1 .and. iex(iuk1) .ne. 7) then
	do i=1,iph
		if (idpatch(i)==id3d(ix3,jy3,kz3)) then
			ikk=iex(iuk1)
			ikk=7-ikk
			iftpts(indirr(ix3,jy3,kz3),ikk,1)=-1
			icpatch(i,ikk,:)=icolu3(:)
			akpatch(i,ikk,1:ipt-1)=ak3(1:ipt-1)
			pchftc(i,ikk)=ak3(ipt)
			goto 1201
		endif
	enddo
	iph=iph+1
	idpatch(iph)=id3d(ix3,jy3,kz3)
	ikk=iex(iuk1)
	ikk=7-ikk
	iftpts(indirr(ix3,jy3,kz3),ikk,1)=-1           !?????????????? check later...
	icpatch(iph,ikk,:)=icolu3(:)
	akpatch(iph,ikk,1:ipt-1)=ak3(1:ipt-1)
	pchftc(iph,ikk)=ak3(ipt)
	1201 continue
	!print *,ix3,jy3,kz3,ikk,' fixed by patch2_1'
endif

call line2cube(ic,jc,kc,ixyz4)
if(irrpts(ic,jc,kc)==1 .and. iex(iuk2) .ne. 7) then
	do i=1,iph
		if (idpatch(i)==id3d(ix4,jy4,kz4)) then
			ikk=iex(iuk2)
			ikk=7-ikk
			iftpts(indirr(ix4,jy4,kz4),ikk,1)=-1
			icpatch(i,ikk,:)=icolu4(:)
			akpatch(i,ikk,1:ipt-1)=ak4(1:ipt-1)
			pchftc(i,ikk)=ak4(ipt)
			goto 1202
		endif
	enddo
	iph=iph+1
	idpatch(iph)=id3d(ix4,jy4,kz4)
	ikk=iex(iuk2)
	ikk=7-ikk
	iftpts(indirr(ix4,jy4,kz4),ikk,1)=-1           !?????????????? check later...
	icpatch(iph,ikk,:)=icolu4(:)
	akpatch(iph,ikk,1:ipt-1)=ak4(1:ipt-1)
	pchftc(iph,ikk)=ak4(ipt)
	1202 continue
	!print *,ix4,jy4,kz4,ikk,' fixed by patch2_2' 
endif

return
end


!-------------------------------------------------------------------------------
!
subroutine substitution(ix,jy,kz,inext)
!-----------------------------------------------
! four points 1,2,9,10 and geometry are shown below
!             \
!              \  
! --*---9---10--\--1---2---*---*------------
!               |
!     inside   /  outside
!-----------------------------------------------------
use pbedata
use comdata
implicit double precision(a-h,o-z)

idx=idre(1,inext)
idy=idre(2,inext)
idz=idre(3,inext)

if(sum(iftpts(indirr(ix+idx,jy+idy,kz+idz),:,1))==0)then
    iph=iph+1
	idpatch(iph)=indirr(ix+idx,jy+idy,kz+idz)
	icpatch(iph,7-inext,1)=id3d(ix,jy,kz)
	akpatch(iph,7-inext,1)=1.d0
    print *,'sub',ix+idx,jy+idy,kz+idz
!	print *,'iph=',iph,idpatch(iph),lexch(inext)
endif

if(sum(iftpts(indirr(ix,jy,kz),:,1))==0)then
    iph=iph+1
	idpatch(iph)=indirr(ix,jy,kz)
	icpatch(iph,inext,1)=id3d(ix+idx,jy+idy,kz+idz)
	akpatch(iph,inext,1)=1.d0
	print *,'sub',ix,jy,kz
!	print *,'iph=',iph,idpatch(iph),inext
endif

return
end


!------------------------------------------------------------
subroutine line2cube(ic,jc,kc,iline)
use comdata
implicit double precision(a-h,o-z)
!kc=mod(iline,nz)
!jc=mod(iline-mod(iline,nz),nz*ny)/nz+1
!ic=(iline-mod(iline,nz)-mod(iline-mod(iline,nz),nz*ny))/nz/ny+1		

ic=mod(iline,nx)
if (ic==0) then
	ic=nx
else
	continue
endif
jc=mod(iline-ic,nx*ny)/nx+1
kc=(iline-ic-mod(iline-ic,nx*ny))/nx/ny+1		
End

!-------------------------------------------------------------
function id3d(ic,jc,kc)
use comdata
implicit double precision(a-h,o-z)
!id3d=(ic-1)*ny*nz+(jc-1)*nz+kc
id3d=(kc-1)*ny*nx+(jc-1)*nx+ic
End

function iY2G(ii)
implicit double precision(a-h,o-z)
integer iarray(7)
iarray=(/3,7,4,2,5,1,6/)
iY2G=iarray(ii)
end

function iG2Y(ii)
implicit double precision(a-h,o-z)
integer iarray(7)
iarray=(/6,4,1,3,5,7,2/)
iG2Y=iarray(ii)
end


!--------------------------------------------------------------
subroutine readin_dist_angle
use pbedata
use comdata
use molecule
implicit double precision(a-h,o-z)

allocate(dist_angle(5,7,nirr))
dist_angle=0.d0

do i=1,nx
    do j=1,ny
	    do k=1,nz
		    if(irrpt(i,j,k)==1)then
			    lir=indirr(i,j,k)
				if(mcx(lir)==1)then    ! x+ direction
				    t3=dxl(lir)
					dist_angle(5,4,lir)=t3
	                dist_angle(1:4,4,lir)=clocal(:,1,lir)
				endif
				if(abs(irrpt(i-1,j,k))==1)then
				    if(mcx(indirr(i-1,j,k))==1)then  ! x-
						lirn=indirr(i-1,j,k)
						t1=-dxr(lirn)
						dist_angle(5,3,lir)=t1
						dist_angle(1:4,3,lir)=clocal(:,1,lirn)
					endif
				endif
				if(mcy(lir)==1)then   ! y+
				    t5=dyl(lir)
					dist_angle(5,5,lir)=t5
					dist_angle(1:4,5,lir)=clocal(:,2,lir)
				endif
				if(abs(irrpt(i,j-1,k))==1)then
					if(mcy(indirr(i,j-1,k))==1)then  ! y-
						lirn=indirr(i,j-1,k)
						t4=-dyr(lirn)
						dist_angle(5,2,lir)=t4
						dist_angle(1:4,2,lir)=clocal(:,2,lirn)
				    endif
				endif
				if(mcz(lir)==1)then  ! z+
				    t7=dzl(lir)
					dist_angle(5,6,lir)=t7
					dist_angle(1:4,6,lir)=clocal(:,3,lir)
                endif
				if(abs(irrpt(i,j,k-1))==1)then
					if(mcz(indirr(i,j,k-1))==1)then ! z-
						lirn=indirr(i,j,k-1)
						t6=-dzr(lirn)
						dist_angle(5,1,lir)=t6
						dist_angle(1:4,1,lir)=clocal(:,3,lirn)
					endif
				endif
			elseif(irrpt(i,j,k)==-1)then   ! on_interface
			    lir=indirr(i,j,k)
				do i_dir=1,7
					dist_angle(5,i_dir,lir)=0.d0
					dist_angle(1:4,i_dir,lir)=clocal(:,1,lir)
				enddo
				! print *,real(clocal(:,1,lir))
			endif
		enddo
	enddo
enddo

return
end

!------------------------------------------------------------------------------
subroutine getfictitious
use molecule
use comdata
use pbedata
implicit double precision(a-h,o-z)

!iph=0 ; itnew=0
!itrblini=itrbl
!do i=1,itrblini
!	call n2ijk(itps(i,1),i1,j1,k1)
!	inext=itps(i,2)
!   call settlecoefrp(i1,j1,k1,inext)   ! call patch 1 and patch 2
!enddo

!print *,'itnew= ', itnew
!print *,'iph= ',iph       
iphnow=iph
do i=1,itnew
    ix=itpsnew(i,1) 
    jy=itpsnew(i,2)
    kz=itpsnew(i,3)
    inext=itpsnew(i,4)
		   
    lpp=0  ;  lnn=0
    if(itpsnew(i,8)==1)then
        lpp=1
        do j=1,iphnow
	        if(idpatch(j)==id3d(ix,jy,kz)) lpp=0
	    enddo		 
	endif
    if(itpsnew(i,5)==1)then
        lnn=1
        in=ix+idre(1,inext)
        jn=jy+idre(2,inext)
        kn=kz+idre(3,inext)
        do j=1,iphnow
            if(idpatch(j)==id3d(in,jn,kn)) lnn=0
        enddo
    endif
    if(lpp+lnn.ne.0)then
        if(itpsnew(i,6)+itpsnew(i,7)==1.and.sum(itpsnew(i,5:8))<=3)then
			!print *, ix,jy,kz,inext,'need patch3'
			call replace3(ix,jy,kz,inext,itpsnew(i,5:8),lpp,lnn)
		elseif(itpsnew(i,6)+itpsnew(i,7)==0)then
		    print *,ix,jy,kz,'need replace4'
			call findfctin(ix,jy,kz,inext,4)
		else
		    call substitution(ix,jy,kz,inext)
		    print *,'sub used',itpsnew(i,5:8)
		endif
	endif
enddo

!print *,'number of patch applied iph=',iph

return
end


!-------------------------------------------------------------------------------
subroutine replace3(ix,jy,kz,inext,iun,lpp,lnn)
use molecule
use pbedata
use comdata
use ModInput
!---------------------------------------------------------------
!this subroutine give the representation of two fictitious points for inside irregular point.
!ipde(1:2) indicates which two partial derivatives include fictitious values
!ipde(3:4) indicates which two partial derivatives are given by real values
!ipde(5:6) indicates which two partial derivatives are to be cancelled
!aiums(2,2,ls-1) the coefficients of real values on the line of two irregualr point 
!                before solving (two equations,negative/positive side,# of points on each side)
!aiuas(2,2,ls-1,ls) the coefficients of real values for interpolation before solving
!                   (two equations, two partial derivatives, # of lines for each derivatives, # of points on each line)
!ack(2,2) the coefficients of fictitious
!ajums(2,4) the coefficients of jump conditions for two equations.(two equations, four jump conditions)
!addpt(3,2) the coef of point 5 and jump conditions for the nearby interface
!-------------------------------------------------------------
implicit double precision(a-h,o-z)
dimension dp(3,3),dw(2),iums(2,ls-1,3),iuas(2,ls-1,ls,3),coiums(2,ls-1,3),coiuas(2,ls-1,ls,3)&
          ,xfw(0:ls-2+(ls-1)/2),idn(3),weitp(0:ls-2+(ls-1)/2,0:1),weitn(0:ls-2+(ls-1)/2,0:1)&
		  ,winter(2,ls-1,0:ls-1,0:0),idnn(3,2),ipde(6),yfw(0:ls-1),pweit(2,0:ls-1,0:1),ajump(3),cpde(4)&
		  ,aiums(3,2,ls-1),aiuas(3,2,ls-1,ls),ajums(3,4),ack(3,3),indxack(3),aifs(2,2,0:(ls-1)/2-1)&
		  ,weitp2(0:ls-1,0:0),weitn2(0:ls-1,0:0),addjump(3),iun(4),icolu3(50),ak3(50),icolu4(50),ak4(50)&
		  ,ixx1(3),ixx2(3),intppts(2,6,3),rjp(4),intppts_temp(6,3),flux(3)                    ! new added

ack=0.d0 ; aiums=0.d0 ; aiuas=0.d0 ; ajums=0.d0

xx=x(ix)
yy=y(jy)
zz=z(kz)

idx=idre(1,inext)
idy=idre(2,inext)
idz=idre(3,inext)

inbx=ix+idx
inby=jy+idy
inbz=kz+idz

idn(:)=(/abs(idx),abs(idy),abs(idz)/) ! line of two irregular points parrallel to this axis
idnn(:,1)=(/abs(idz),abs(idx),abs(idy)/) ! other two axis
idnn(:,2)=(/abs(idy),abs(idz),abs(idx)/)

! to find u_o on the interface
! to find u_o on the interface
if(isf==0)then
	t1=min(0.d0,dx*idx+dy*idy+dz*idz)
	t2=max(0.d0,dx*idx+dy*idy+dz*idz)
	t1=t1+xx*abs(idx)+yy*abs(idy)+zz*abs(idz)
	t2=t2+xx*abs(idx)+yy*abs(idy)+zz*abs(idz)
	t=rtbis(t1,t2,1.d-10,xx,yy,zz,int(abs(idx)+2*abs(idy)+3*abs(idz)))
	xo=xx+(t-xx)*abs(idx)  ; yo=yy+(t-yy)*abs(idy) ; zo=zz+(t-zz)*abs(idz)
	ctht=costheta(xo,yo,zo)           
	stht=sintheta(xo,yo,zo)
	cphi=cosphi(xo,yo,zo)
	sphi=sinphi(xo,yo,zo)
else
	t=dist_angle(5,inext,indirr(ix,jy,kz))
	xo=xx+t*abs(idx)  ; yo=yy+t*abs(idy) ; zo=zz+t*abs(idz)

	sphi=dist_angle(1,inext,indirr(ix,jy,kz))
	cphi=dist_angle(2,inext,indirr(ix,jy,kz))
	ctht=dist_angle(3,inext,indirr(ix,jy,kz))
	stht=dist_angle(4,inext,indirr(ix,jy,kz))
endif

call coortrans(dp,sphi,cphi,stht,ctht)

btp=betan(xo,yo,zo)		!use my beta
btn=betap(xo,yo,zo)




! change the jump condtion for pb3 **********************
 
        xs=xo;	ys=yo;	zs=zo
        ik=inext
        
		if (ik==1 .or. ik==6) then
		    idn1=5
		elseif (ik==2 .or. ik==5) then    
		    idn1=3
		elseif (ik==3 .or. ik==4) then
		    idn1=1
		else
		    print *,'ik is wrong....'
		endif
		!print *,'ik= ',ik,' idn',idn1
		if (icg==0) then						! Directly set the jump condition
			call findjumps(xs,ys,zs,dP,rjp)		! Use this if the true solution known
		elseif (icg==2) then
                        call findRegJumps(xs,ys,zs,dP,rjp,eps0)
		endif


ww=-rjp(1)
vv=-rjp(2)
dw=-rjp(3:4)


! Note if for pbe, we need set ww,vv,dw as 0 here. 

do i=1,ls-1
    iums(1,i,:)=(/ix-(i-1)*idx,jy-(i-1)*idy,kz-(i-1)*idz/) ! the first main set is inside
	coiums(1,i,:)=(/x(iums(1,i,1)),y(iums(1,i,2)),z(iums(1,i,3))/) ! give the coordinate of those points
	iums(2,i,:)=(/ix+i*idx,jy+i*idy,kz+i*idz/)             ! the second main set is outside
	coiums(2,i,:)=(/x(iums(2,i,1)),y(iums(2,i,2)),z(iums(2,i,3))/)
enddo
do i=3,ls-1
    if(io(iums(1,i,1),iums(1,i,2),iums(1,i,3))==1)then
		print *,ix,jy,kz,inext,'replace3--> sub'
		call substitution(ix,jy,kz,inext)
    endif
	if(io(iums(2,i,1),iums(2,i,2),iums(2,i,3)).ne.1) write(*,*) 'inside.f90--not enough outside points'
enddo
in=iums(2,1,1) ; jn=iums(2,1,2) ; kn=iums(2,1,3)
ip=iums(1,1,1) ; jp=iums(1,1,2) ; kp=iums(1,1,3)

! ipde(i) could be (1:6), means (u^+_x,u^-_x,u^+_y,u^-_y,u^+_z,u^-_z)
! where i could be (1:6), means (two main derivatives, two auxiliary derivatives, two cancelled derivatives)
ipde(1)=idn(1)+3*idn(2)+5*idn(3)  
ipde(2)=2*idn(1)+4*idn(2)+6*idn(3) 


do i=1,2             ! search auxiliary points for partial derivative of other two direction
    call searchauxi(i,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,-1,ijust)
	!if (ijust ==-1) print *,i,ix,jy,kz,idx,idy,idz,idn,idnn,'lsf'
	if(ijust==1)then  ! derivative from negative side found
		ipde(2+i)=2*idnn(1,i)+4*idnn(2,i)+6*idnn(3,i) ! store derivative from negative side 
		ipde(4+i)=ipde(2+i)-1                         ! derivative from positive side will be cancelled
	else
		
		call searchauxi(i,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,1,ijust) 
		!print *,i,ix,jy,kz,idx,idy,idz,idn,idnn,'lsf'
		!print *,'ijust=',ijust
		if(ijust==1)then
			ipde(2+i)=idnn(1,i)+3*idnn(2,i)+5*idnn(3,i) ! store derivative from positive side
			ipde(4+i)=ipde(2+i)+1                       ! derivative from negative side will be cancelled
	    else             ! if searchauxi doesn't work, then try searchauxi4
		    
		    call searchauxi4(i,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,-1,ijust)
	        !print *,'patch 4 applied in patch 3'
	        !print *,ix,jy,kz,inext
	        !print *,'ijust=',ijust
	        !print *,iuas,real(coiuas)
	        if(ijust==1)then  ! derivative from negative side found
	            ipde(2+i)=2*idnn(1,i)+4*idnn(2,i)+6*idnn(3,i) ! store derivative from negative side 
		        ipde(4+i)=ipde(2+i)-1                         ! derivative from positive side will be cancelled
	        else
		        call searchauxi4(i,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,1,ijust) 
	            
	            if(ijust==1)then
			        ipde(2+i)=idnn(1,i)+3*idnn(2,i)+5*idnn(3,i) ! store derivative from positive side
			        ipde(4+i)=ipde(2+i)+1                       ! derivative from negative side will be cancelled
	            else
		            print *,ix,jy,kz,inext,'replace3 ->not enough point for partial derivatives'
                    call substitution(ix,jy,kz,inext)
			        return
  		        endif
            endif
 		endif
	endif

	do j=1,ls-1
	    do k=1,ls
			xfw(k-1)=sum(coiuas(i,j,k,:)*idn(:))
	    enddo
		call weights(xo*idn(1)+yo*idn(2)+zo*idn(3),xfw(0:ls-1),ls-1,ls-1,0,winter(i,j,:,:))
		yfw(j)=sum(coiuas(i,j,1,:)*idnn(:,i))
	enddo
	yfw(0)=xx*idnn(1,i)+yy*idnn(2,i)+zz*idnn(3,i)
	call weights(xo*idnn(1,i)+yo*idnn(2,i)+zo*idnn(3,i),yfw,ls-1,ls-1,1,pweit(i,:,:))
enddo

call simpjumpcon(ix,jy,kz,inext,dp,btp,btn,ipde,ajump,cpde)

ack=0.d0 ; aiuas=0.d0 ; aiums=0.d0 
ajums=0.d0 ; indxack=0.d0 ; aifs=0.d0
weitp=0.d0 ; weitn=0.d0  ; xfw=0.d0

! the weight of points f^+ and iums(1,:,:), xfw is the location for finding the weight 
xfw(0)=sum(coiums(2,1,:)*idn(:)) ! f^+ are the first ie points in um2 set
do i=1,ls-1
	xfw(i)=sum(coiums(1,i,:)*idn(:))
enddo
call weights(xo*idn(1)+yo*idn(2)+zo*idn(3),xfw(0:ls-1),ls-1,ls-1,1,weitp(0:ls-1,:))  

! the weight of points f^- and iums(2,:,:) 
xfw(0)=sum(coiums(1,1,:)*idn(:)) ! f^- are the first ie points in um1 set
do i=1,ls-1
	xfw(i)=sum(coiums(2,i,:)*idn(:))
enddo
call weights(xo*idn(1)+yo*idn(2)+zo*idn(3),xfw(0:ls-1),ls-1,ls-1,1,weitn(0:ls-1,:)) 

do j=1,ls-1
	aiums(2,1,j)=-cpde(1)*weitp(j,1)
	aiums(2,2,j)=-cpde(2)*weitn(j,1)
	aiums(1,1,j)=-weitp(j,0)
	aiums(1,2,j)=weitn(j,0)
enddo

aifs(2,1,0)=-cpde(1)*weitp(0,1)
aifs(2,2,0)=-cpde(2)*weitn(0,1)
aifs(1,1,0)=-weitp(0,0)
aifs(1,2,0)=weitn(0,0)

do i=3,4  ! the partial derivative on other two directions
	if(mod(ipde(i),2)==1)then  ! the partial derivative on positive side is represented.
		do j=1,ls-1		
			aiums(2,1,j)=aiums(2,1,j)-cpde(i)*pweit(i-2,0,1)*weitp(j,0)
		enddo
		aifs(2,1,0)=aifs(2,1,0)-cpde(i)*pweit(i-2,0,1)*weitp(0,0)
		!ack(2,1)=ack(2,1)+cpde(i)*pweit(i-2,0,1)*weitp(0,0)
	else   ! the partial derivative on negative side is represented
		do j=1,ls-1
			aiums(2,2,j)=aiums(2,2,j)-cpde(i)*pweit(i-2,0,1)*weitn(j,0)
		enddo
		aifs(2,2,0)=aifs(2,2,0)-cpde(i)*pweit(i-2,0,1)*weitn(0,0)
		!ack(2,2)=ack(2,2)+cpde(i)*pweit(i-2,0,1)*weitn(0,0)
	endif
enddo

ack(1,1:2)=-aifs(1,:,0)
ack(2,1:2)=-aifs(2,:,0)

ajums(2,2:4)=ajump(:)
ajums(2,1)=0.d0
ajums(1,2:4)=0.d0
ajums(1,1)=1.d0

do i=3,4
	do j=1,ls-1
		do k=1,ls
			aiuas(2,i-2,j,k)=-cpde(i)*pweit(i-2,j,1)*winter(i-2,j,k-1,0)
			aiuas(1,i-2,j,k)=0.d0
		enddo
	enddo
enddo
!---------------------------------------------------------
!the following add one more equation
!----------------------------------------------------------
if(iun(2)==1)then   ! the equation is about interface between 1- and 2+

	ix2=ix+2*idx ;  jy2=jy+2*idy ;  kz2=kz+2*idz
	x2=x(ix2) ;  y2=y(jy2) ; z2=z(kz2)	
	if(isf==0)then
		t1=min(0.d0,-dcel*idx-dcel*idy-dcel*idz)
		t2=max(0.d0,-dcel*idx-dcel*idy-dcel*idz)
		t1=t1+x2*abs(idx)+y2*abs(idy)+z2*abs(idz)
		t2=t2+x2*abs(idx)+y2*abs(idy)+z2*abs(idz)
		t=rtbis(t1,t2,1.d-10,x2,y2,z2,int(abs(idx)+2*abs(idy)+3*abs(idz)))
		t=t-x2*abs(idx)-y2*abs(idy)-z2*abs(idz)
		!pause 'patch3 for analytic surface, underconstruction'
	else
		t=dist_angle(5,7-inext,indirr(ix2,jy2,kz2))
    endif
	xo2=x2+t*abs(idx)  ; yo2=y2+t*abs(idy) ; zo2=z2+t*abs(idz)
    ! Adjust the jump condition
    IF (ibd .ne. 0) THEN 
        if (icg==2) then   !mibpb_newReg
            ujump2=-phi_star(xo2,yo2,zo2)
        else
            UJUMP2=0.D0                 !mibpb2, mibpb3
        endif
	ELSE 
       ujump2=-uj(xo2,yo2,zo2)
	ENDIF 

    ack(1,3)=-aiums(1,2,2)    ! obtain coef of f2 from rightside of eq1 and 2,
	ack(2,3)=-aiums(2,2,2)    ! aiums(eq#,outside,the second pt)

    xfw(0:ls-1)=(/-dcel,0.d0,dcel/)   

	call weights(dcel-abs(t),xfw(0:ls-1),ls-1,ls-1,0,weitn2)  ! in order (4,1,2)
	weitp2=weitn2
	    
    ack(3,:)=(/weitp2(1,0),-weitn2(0,0),-weitn2(2,0)/)    ! coef of (f1,f4,f2)
	aiums(:,2,2)=(/0.d0,0.d0,-weitp2(2,0)/)  ! coef of u2+
	aiums(3,2,1)=weitn2(1,0)                 ! coef of u1-
	aiums(3,1,1)=-weitp2(0,0)                ! coef of u4+
	addjump(:)=(/0.d0,0.d0,ujump2/)          ! ujump for interface 2

elseif(iun(3)==1)then  ! the interface is between 3 and 4

	if(isf==0)then
		t1=min(0.d0,-dcel*idx-dcel*idy-dcel*idz)
		t2=max(0.d0,-dcel*idx-dcel*idy-dcel*idz)
		t1=t1+xx*abs(idx)+yy*abs(idy)+zz*abs(idz)
		t2=t2+xx*abs(idx)+yy*abs(idy)+zz*abs(idz)
		t=rtbis(t1,t2,1.d-10,xx,yy,zz,int(abs(idx)+2*abs(idy)+3*abs(idz)))
		t=t-xx*abs(idx)-yy*abs(idy)-zz*abs(idz)
		!pause 'patch 3 for analytic surface, underconstruction'
	else
		t=dist_angle(5,7-inext,indirr(ix,jy,kz))
    endif
	xo2=xx+t*abs(idx)  ; yo2=yy+t*abs(idy) ; zo2=zz+t*abs(idz)
    
    IF (ibd .ne. 0) THEN
        if (icg==2) then   !mibpb_newReg
            ujump2=-phi_star(xo2,yo2,zo2)
        else
            UJUMP2=0.D0                 !mibpb2, mibpb3
        endif
    ELSE
       ujump2=-uj(xo2,yo2,zo2)
    ENDIF    

    ack(1,3)=-aiums(1,1,2)    ! obtain coef of f3 from rightside of eq1 and 2,
	ack(2,3)=-aiums(2,1,2)    ! aiums(eq#,outside,the second pt)
    xfw(0:ls-1)=(/-dcel,0.d0,dcel/)   

    call weights(-abs(t),xfw(0:ls-1),ls-1,ls-1,0,weitp2)  ! in order (3,4,1)
	weitn2=weitp2

    ack(3,:)=(/weitp2(2,0),-weitn2(1,0),weitp2(0,0)/)    ! coef of (f1,f4,f3)
  	aiums(:,1,2)=(/0.d0,0.d0,weitn2(0,0)/)  ! coef of u3-
   	aiums(3,1,1)=-weitp2(1,0)               ! coef of u4+
    aiums(3,2,1)=weitn2(2,0)                ! coef of u1-
   	addjump(:)=(/0.d0,0.d0,ujump2/)          ! ujump for interface 2
  
endif

!-----------------------------------------------------------
call ludcmp(ack,3,3,indxack,d)
do i=1,2
	do j=1,ls-1
		call lubksb(ack,3,3,indxack,aiums(:,i,j))
	enddo
	do j=1,ls-1
		do k=1,ls
			call lubksb(ack,3,3,indxack,aiuas(:,i,j,k))
		enddo
	enddo
enddo
call lubksb(ack,3,3,indxack,addjump(:))
do i=1,4
	call lubksb(ack,3,3,indxack,ajums(:,i))
enddo
ajums(:,1)=ww*ajums(:,1)
ajums(:,2)=vv*ajums(:,2)
ajums(:,3)=dw(1)*ajums(:,3)
ajums(:,4)=dw(2)*ajums(:,4)

icount=0

do i=1,2
	do j=1,ls-1
		icount=icount+1
		icolu3(icount)=id3d(iums(i,j,1),iums(i,j,2),iums(i,j,3))
		ak3(icount)=aiums(1,i,j)
		icolu4(icount)=icolu3(icount)
		ak4(icount)=aiums(2,i,j)
	enddo
enddo
do i=1,2
	do j=1,ls-1
		do k=1,ls
			icount=icount+1
			icolu3(icount)=id3d(iuas(i,j,k,1),iuas(i,j,k,2),iuas(i,j,k,3))
			ak3(icount)=aiuas(1,i,j,k)
			icolu4(icount)=icolu3(icount)
			ak4(icount)=aiuas(2,i,j,k)
		enddo
	enddo
enddo

icount=icount+1
icolu3(icount)=0
ak3(icount)=sum(ajums(1,:))+addjump(1)
icolu4(icount)=0
ak4(icount)=sum(ajums(2,:))+addjump(2)


ipt=icount ! add by weihua

if(lnn==1)then
do i=1,iph
    if(idpatch(i)==id3d(in,jn,kn))then
		iftpts(indirr(in,jn,kn),inext,1)=-1
		icpatch(i,inext,:)=icolu3(:)
	    akpatch(i,inext,1:ipt-1)=ak3(1:ipt-1)
		pchftc(i,inext)=ak3(ipt)
		goto 1205
	endif
enddo
iftpts(indirr(in,jn,kn),inext,1)=-1
iph=iph+1
idpatch(iph)=id3d(in,jn,kn)
!idpatchreal(iph,:)=(/in,jn,kn/)
icpatch(iph,inext,:)=icolu3(:)
akpatch(iph,inext,1:ipt-1)=ak3(1:ipt-1)
pchftc(iph,inext)=ak3(ipt)

1205 continue
!print *,in,jn,kn,inext,' patch3_1'
endif

if(lpp==1)then
do i=1,iph
    if(idpatch(i)==id3d(ix,jy,kz))then
	   iftpts(indirr(ix,jy,kz),7-inext,1)=-1
	   icpatch(i,7-inext,:)=icolu4(:)
	   akpatch(i,7-inext,1:ipt-1)=ak4(1:ipt-1)
	   pchftc(i,7-inext)=ak4(ipt)
	   !print *,ipt
	   !print *,icpatch(i,inext,:)
	   !print *,ix,jy,kz,in,jn,kn
       !print *,akpatch(i,inext,1:ipt-1)
       !print *,pchftc(i,inext)
	goto 1206
	endif
enddo
iftpts(indirr(ix,jy,kz),7-inext,1)=-1
iph=iph+1
idpatch(iph)=id3d(ix,jy,kz)
!idpatchreal(iph,:)=(/ix,jy,kz/)
icpatch(iph,7-inext,:)=icolu4(:)
akpatch(iph,7-inext,1:ipt-1)=ak4(1:ipt-1)
pchftc(iph,7-inext)=ak4(ipt)
!print *,ipt
!print *,icpatch(iph,inext,:)
!print *,akpatch(iph,inext,:)
!print *,pchftc(iph,inext)

1206 continue
!print *,ix,jy,kz,7-inext,' patch3_2'
endif

return
end

!------------------------------------------------------------------------------------------
!
subroutine simpjumpcon(ix,jy,kz,inext,dp,btp,btn,ipde,ajump,cpde)
!---------------------------------------------------------------
!this subroutine is to give equations to solve for fictitious values
!--------------------------------------------------------------

implicit double precision(a-h,o-z)
dimension dp(3,3),cjump(3,6),ipde(6),ajump(3),cpde(4)

cjump(1,:)=(/dp(1,1)*btp,-dp(1,1)*btn,dp(1,2)*btp,-dp(1,2)*btn,dp(1,3)*btp,-dp(1,3)*btn/)
cjump(2,:)=(/dp(2,1),-dp(2,1),dp(2,2),-dp(2,2),dp(2,3),-dp(2,3)/)
cjump(3,:)=(/dp(3,1),-dp(3,1),dp(3,2),-dp(3,2),dp(3,3),-dp(3,3)/)

i=ipde(5) ; j=ipde(6)
m1=ipde(1) ; m2=ipde(2)

summi=sum(abs(cjump(:,m1)*cjump(:,i)))
summj=sum(abs(cjump(:,m1)*cjump(:,j)))
sumij=sum(abs(cjump(:,i)*cjump(:,j)))

kk=0 ; ks=0
if(summi==0.d0.and.summj==0.d0)then
	k=1 
	do while(k.le.3.and.kk==0)
		if(cjump(k,i)==0.and.cjump(k,j)==0) kk=k
		k=k+1
	enddo
elseif(summi==0.d0.and.sumij==0.d0)then
	ileft=j
	k=1 
	do while(k.le.3.and.ks==0)
		if(cjump(k,m1)==0) ks=k
		k=k+1
	enddo
elseif(summj==0.d0.and.sumij==0.d0)then
	ileft=i
	k=1 
	do while(k.le.3.and.ks==0)
		if(cjump(k,m1)==0) ks=k
		k=k+1
	enddo
endif

if(kk==0.and.ks==0)then
    ajump(1)=cjump(2,i)*cjump(3,j)-cjump(2,j)*cjump(3,i)
	ajump(2)=cjump(3,i)*cjump(1,j)-cjump(3,j)*cjump(1,i)
	ajump(3)=cjump(1,i)*cjump(2,j)-cjump(1,j)*cjump(2,i)
	do l=1,4
	    m=ipde(l)
		cpde(l)=ajump(1)*cjump(1,m)+ajump(2)*cjump(2,m)+ajump(3)*cjump(3,m)
	enddo
elseif(kk.ne.0)then
    ajump=0.d0
	ajump(kk)=1.d0
	do l=1,4
	    m=ipde(l)
		cpde(l)=cjump(kk,m)
	enddo
else
    if(ks==1)then
        k1=2 ; k2=3
	elseif(ks==2)then
	    k1=3 ; k2=1
	else
	    k1=1 ; k2=2
	endif

	ajump(ks)=0.d0
	ajump(k1)=cjump(k2,ileft)
	ajump(k2)=-cjump(k1,ileft)
	do l=1,4
		m=ipde(l)
		cpde(l)=cjump(k2,ileft)*cjump(k1,m)-cjump(k1,ileft)*cjump(k2,m)
	enddo
endif

if(cpde(1)==0.or.cpde(2)==0)then
    print *,dp,ipde,ix,jy,kz
	write(*,*) 'inside_3d.f90 simjumpcon---> cpde(1) or cpde(2)==0'
endif

return
end


!---------------------------------------------------------------------------------------
function costheta(x,y,z)
implicit double precision(a-h,o-z)

psix=varphix(x,y,z)
psiy=varphiy(x,y,z)
psiz=varphiz(x,y,z)

psinorm=sqrt(psix**2+psiy**2+psiz**2)  !psinorm is |\grad psi|
cphi=psiz/psinorm
sphi=sqrt(1.d0-cphi**2)
psinorm2=psinorm*sphi                 !psinorm2 is |\grad psi2|

if(psinorm2<1.d-15)then                ! the normal direction perpendicular to xy 
    costheta=1.d0
else
	costheta=psix/psinorm2
endif

return
end
!
!------------------------------------------------------------------------------------
!
function sintheta(x,y,z)
implicit double precision(a-h,o-z)

psix=varphix(x,y,z)
psiy=varphiy(x,y,z)
psiz=varphiz(x,y,z)

psinorm=sqrt(psix**2+psiy**2+psiz**2)  !psinorm is |\grad psi|
cphi=psiz/psinorm
sphi=sqrt(1.d0-cphi**2)
psinorm2=psinorm*sphi                 !psinorm2 is |\grad psi2|

if(psinorm2<1.d-15)then                ! the normal direction perpendicular to xy 
    sintheta=0.d0
else
	sintheta=psiy/psinorm2
endif

return
end
!
!------------------------------------------------------------------------------------
!
function cosphi(x,y,z)
implicit double precision(a-h,o-z)

psix=varphix(x,y,z)
psiy=varphiy(x,y,z)
psiz=varphiz(x,y,z)

psinorm=sqrt(psix**2+psiy**2+psiz**2)  !psinorm is |\grad psi|
cosphi=psiz/psinorm

return
end
!
!------------------------------------------------------------------------------------
!
function sinphi(x,y,z)
implicit double precision(a-h,o-z)

psix=varphix(x,y,z)
psiy=varphiy(x,y,z)
psiz=varphiz(x,y,z)

psinorm=sqrt(psix**2+psiy**2+psiz**2)  !psinorm is |\grad psi|
cphi=psiz/psinorm
sinphi=sqrt(1.d0-cphi**2)

return
end


!-------------------------------------------------------------------
subroutine coortrans(dp,sphi,cphi,stht,ctht)
implicit double precision(a-h,o-z)
real*8,dimension(3,3):: dp

dp(1,1)=sphi*ctht ; dp(1,2)=sphi*stht ; dp(1,3)=cphi
dp(2,1)=-stht      ; dp(2,2)=ctht     ; dp(2,3)=0.d0
dp(3,1)=-cphi*ctht ; dp(3,2)=-cphi*stht ; dp(3,3)=sphi

return
end

!-------------------------------------------------------------------
subroutine searchauxi4(i4,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,iio,ijust)
use pbedata
use comdata
implicit double precision(a-h,o-z)
!-----------------------------------------------------------------------
! this subroutine is to search (ls-1)*ls auxiliary points for partial derivative on other two sides
! iio=1 means searching inside, iio=-1 means searching outside
!------------------------------------------------------------------------
dimension idn(3),idnn(3,2),iuas(2,ls-1,ls,3),coiuas(2,ls-1,ls,3),iuy(ls),iuy4(4),iuast(4,3)&
          ,iuasl(ls,ls,3),lid(3)
dimension iord(ls+1,-1:1),iordl(ls),iord5(5,-1:1),iordl4(4)

	iord=0
	iord(1,1)=-(ls-3)/2 ; iord(2,1)=-(ls-1)/2
	do i=1,(ls-1)/2
		iord(2*i+1,1)=iord(2*i-1,1)+1
		iord(2*i+2,1)=iord(2*i,1)-1
	enddo
	do i=1,ls+1
		iord(i,-1)=iord(i,1)-1
	enddo
	! define the order to determine the distribution of line of auxiliary points
	iordl=0
	iordl(1)=(ls+1)/2
	iordl(2)=iordl(1)-1
	iordl(3)=iordl(1)+1
	do i=2,(ls-1)/2
		iordl(2*i)=iordl(2*i-2)-1
		iordl(2*i+1)=iordl(2*i-1)+1
	enddo
	! define the order to search auiliary points on a line when 3 consective points
	! are not available, work for patch4
	iord5(:,1)=(/-1,0,-2,1,-3/)
	iord5(:,-1)=(/-2,-1,-3,0,-4/)
	! define the order to determine the distribution of line of quxiliary points
	! when 3 consective lines are not available, work for patch4
	iordl4(:)=(/2,1,3,4/)



    ijust=1
    j=1 ; lfind=0
	do while(j.le.ls.and.lfind<ls-1)  ! the distribution of lines of auxiliary points which determines pweit.
	    lfind=0 ; k=1 ; lpfind=ls
		do while(k.le.ls-1.and.lpfind==ls)  ! the # of lines of the auxiliary points
		! the condition lpfind==ls means points are found in previous line.                           
		    if(k<iordl(j))then
			    k1=k-iordl(j)
			else
			    k1=k-iordl(j)+1
			endif
			l=1 ; lpfind=0  ! lpfind=ls if the auxiliary points on one line were found
			do while(l.le.ls+1.and.lpfind==0) ! the distribution of auxiliary points on each line which determines winter.
		        do m=1,ls ! the # of auxiliary points on one line
				    ixau=ix+iord(l,idx)+(m-1)*idn(1)+k1*idnn(1,i4)
					jyau=jy+iord(l,idy)+(m-1)*idn(2)+k1*idnn(2,i4)
					kzau=kz+iord(l,idz)+(m-1)*idn(3)+k1*idnn(3,i4)
					iuas(i4,k,m,:)=(/ixau,jyau,kzau/)
					coiuas(i4,k,m,:)=(/x(ixau),y(jyau),z(kzau)/)
				    if(iuas(i4,k,m,1)>nx.or.iuas(i4,k,m,2)>ny.or.iuas(i4,k,m,3)>nz.or.minval(iuas(i4,k,m,:))<1)then
						ijust=-1
    				    if(iio==1) print *,'searchauxi4--> auxiliary point outside comput domain'
						return
					endif
					iuy(m)=-io(ixau,jyau,kzau)
				enddo
				if((iio==-1.and.maxval(iuy(1:ls))==-1).or.(iio==1.and.minval(iuy(1:ls)).ge.0))then
 					lpfind=ls ; lfind=lfind+1
				endif
				l=l+1
			enddo
			if(lpfind==0)then
			    l=1
				do while(l.le.5.and.lpfind==0)
			        do m=1,4 ! the # of auxiliary points on one line
					    ixau=ix+iord5(l,idx)+(m-1)*idn(1)+k1*idnn(1,i4)
						jyau=jy+iord5(l,idy)+(m-1)*idn(2)+k1*idnn(2,i4)
						kzau=kz+iord5(l,idz)+(m-1)*idn(3)+k1*idnn(3,i4)
						iuast(m,:)=(/ixau,jyau,kzau/)
						if(iuast(m,1)>nx.or.iuast(m,2)>ny.or.iuast(m,3)>nz.or.minval(iuast(m,:))<1)then
							ijust=-1
    						if(iio==1) print *,'searchauxi4 --> auxiliary point outside comput domain'
							return
						endif
						if(io(ixau,jyau,kzau)==1)then
						    iuy4(m)=-1
						else
						    iuy4(m)=1
						endif
					enddo
					if((iio==-1.and.sum(iuy4(1:4)).le.-2).or.(iio==1.and.sum(iuy4(1:4)).ge.2))then
 						lpfind=ls ; lfind=lfind+1
						m=1 ; mk=1
						do while(m.le.4.and.mk.le.3)
						    if(iuy4(m)==iio)then
								iuas(i4,k,mk,:)=iuast(m,:)
								coiuas(i4,k,mk,:)=(/x(iuast(m,1)),y(iuast(m,2)),z(iuast(m,3))/)
								mk=mk+1
							endif
							m=m+1
						enddo
					endif
					l=l+1
				enddo
			endif
			k=k+1
		enddo
		j=j+1 
	enddo

if(lfind<ls-1)then

    j=1 ; lfind=0
	do while(j.le.4.and.lfind<ls-1)  ! the distribution of lines of auxiliary points which determines pweit.
	    lid=0 ; lfind=0 ; k=1 
		do while(k.le.3)  ! the # of lines of the auxiliary points
		    if(k<iordl4(j))then
			    k1=k-iordl4(j)
			else
			    k1=k-iordl4(j)+1
			endif
			l=1 ; lpfind=0  ! lpfind=ls if the auxiliary points on one line were found
			do while(l.le.ls+1.and.lpfind==0) ! the distribution of auxiliary points on each line which determines winter.
		        do m=1,ls ! the # of auxiliary points on one line
				    ixau=ix+iord(l,idx)+(m-1)*idn(1)+k1*idnn(1,i4)
					jyau=jy+iord(l,idy)+(m-1)*idn(2)+k1*idnn(2,i4)
					kzau=kz+iord(l,idz)+(m-1)*idn(3)+k1*idnn(3,i4)
					iuasl(k,m,:)=(/ixau,jyau,kzau/)
				    if(iuasl(k,m,1)>nx.or.iuasl(k,m,2)>ny.or.iuasl(k,m,3)>nz.or.minval(iuasl(k,m,:))<1)then
						ijust=-1
    				    if(iio==1) print *,'searchauxi4--> auxiliary point outside comput domain'
						return
					endif
					iuy(m)=-io(ixau,jyau,kzau)
				enddo
				if((iio==-1.and.maxval(iuy(1:ls))==-1).or.(iio==1.and.minval(iuy(1:ls)).ge.0))then
 					lpfind=ls ; lfind=lfind+1 ; lid(k)=1
				endif
				l=l+1
			enddo
            k=k+1
		enddo
		j=j+1 
	enddo
    if(lfind.ge.ls-1)then
	    k=1 ; kreal=1
		do while(k.le.ls.and.kreal.le.ls-1)
		    if(lid(k)==1)then
			    do m=1,ls
				    iuas(i4,kreal,m,:)=iuasl(k,m,:)
					coiuas(i4,kreal,m,:)=(/x(iuasl(k,m,1)),y(iuasl(k,m,2)),z(iuasl(k,m,3))/)
			    enddo
				kreal=kreal+1
			endif
			k=k+1
		enddo
	else
	    ijust=-1
	endif
endif

return
end
!--------------------------------------------------------------------------------
subroutine searchauxi(i,ix,jy,kz,idx,idy,idz,idn,idnn,iuas,coiuas,iio,ijust)
use pbedata
use comdata
implicit double precision(a-h,o-z)
!-----------------------------------------------------------------------
! this subroutine is to search (ls-1)*ls auxiliary points for partial derivative on other two sides
! iio=1 means searching inside, iio=-1 means searching outside
!------------------------------------------------------------------------
dimension idn(3),idnn(3,2),iuas(2,ls-1,ls,3),coiuas(2,ls-1,ls,3),iuy(ls)
dimension iord(ls+1,-1:1),iordl(ls)

iord=0
iord(1,1)=-(ls-3)/2 ; iord(2,1)=-(ls-1)/2
do ii=1,(ls-1)/2
	iord(2*ii+1,1)=iord(2*ii-1,1)+1
	iord(2*ii+2,1)=iord(2*ii,1)-1
enddo

do ii=1,ls+1
	iord(ii,-1)=iord(ii,1)-1
enddo

! define the order to determine the distribution of line of auxiliary points
iordl=0
iordl(1)=(ls+1)/2
iordl(2)=iordl(1)-1
iordl(3)=iordl(1)+1
do ii=2,(ls-1)/2
	iordl(2*ii)=iordl(2*ii-2)-1
	iordl(2*ii+1)=iordl(2*ii-1)+1
enddo

    ijust=1
    j=1 ; lfind=0
	do while(j.le.ls.and.lfind<ls-1)  ! the distribution of lines of auxiliary points which determines pweit.
	    lfind=0 ; k=1 ; lpfind=ls
		do while(k.le.ls-1.and.lpfind==ls)  ! the # of lines of the auxiliary points
		! the condition lpfind==ls means points are found in previous line.                           
		    if(k<iordl(j))then
			    k1=k-iordl(j)
			else
			    k1=k-iordl(j)+1
			endif
			l=1 ; lpfind=0  ! lpfind=ls if the auxiliary points on one line were found
			do while(l.le.ls+1.and.lpfind==0) ! the distribution of auxiliary points on each line which determines winter.
		        do m=1,ls ! the # of auxiliary points on one line
				    ixau=ix+iord(l,idx)+(m-1)*idn(1)+k1*idnn(1,i)
					jyau=jy+iord(l,idy)+(m-1)*idn(2)+k1*idnn(2,i)
					kzau=kz+iord(l,idz)+(m-1)*idn(3)+k1*idnn(3,i)
					iuas(i,k,m,:)=(/ixau,jyau,kzau/)
					coiuas(i,k,m,:)=(/x(ixau),y(jyau),z(kzau)/)
				    if(iuas(i,k,m,1)>nx.or.iuas(i,k,m,2)>ny.or.iuas(i,k,m,3)>nz.or.minval(iuas(i,k,m,:))<1)then
						ijust=-1
    				    if(iio==1) print *,'inside_3d subroutin findficin--> auxiliary point outside comput domain'
						return
!    				    pause 'inside_3d subroutin findficin--> auxiliary point outisde'
					endif
					iuy(m)=-io(ixau,jyau,kzau)
				enddo
				if((iio==-1.and.maxval(iuy(1:ls))==-1).or.(iio==1.and.minval(iuy(1:ls)).ge.0))then
 					lpfind=ls ; lfind=lfind+1
				endif
				l=l+1
			enddo
			k=k+1
		enddo
		j=j+1 
	enddo

if(lfind<ls-1) ijust=-1

return
end
!
