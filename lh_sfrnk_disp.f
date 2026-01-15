      subroutine sfrnk_lh(
     & X_ar,Y_ar,T_ar,mass_ar,nbulk, 
     & npar_in,nperp_in,D_0)
      
c-----Sam Frank's modified version of the LH dispersion relation from
c     Bonoli/Englade 1986 which includes a thermal correction to the
c     episilon_|| term. This ensures proper ray bending in hot plasmas
c     for more information on this effect see Wright/Bertelli PPCF 2014.
c     This routine uses the zetbes routines from TORIC for calculation
c     of the thermal correction from the plasma dispersion function.
      implicit none
      include 'param.i'
      include 'eps.i'
c-----inputs
      integer nbulk
      real*8 X_ar(*), Y_ar(*), T_ar(*), mass_ar(*),
     & npar_in,nperp_in,r,z,phi,nosigma
      integer i_deriv
c     nbulk is a number of the plasma species
c     mass_ar(nbulk) - the mass of the given species (in electron mass) 
c     X_ar(nbulk) = (fp/f)**2
c     Y_ar(nbulk) = fc/f ! for electron Ye should has the sign opposite to Yi
c     T_ar(nbulk)=temperature in eV  
c     npar_in - parallel index of refraction n.
c     nperp_in - perpendicular index of refraction n.

c-----output
      real*8 D_0

c-----local
      real*8 c,k4,npar,nperp,npars,nperps,nperp4,nperp6,
     & eps_perp,eps_par,eps_xy,
     & vt,vt_dc,
     & P_0,P_2,P_4,P_6,
     & beta_p6,beta_p4,beta_p2,beta_p0,
     & frqncpl,frqncmn,df,hfrqnc,
     & xe0

      real*8 cz1, cz0, cdze
      
      integer i
      nosigma=0.d0
      beta_p6=1.d0
      beta_p4=1.d0
      beta_p2=1.d0
      beta_p0=1.d0

      c = 3.0d10                !speed of light
      k4 = 4.19d7               !constant for elec therm vel
      
      npar=npar_in
      nperp=nperp_in

      npars = npar**2
      nperps = nperp**2
      nperp4 = nperps**2
      nperp6 = nperp4*nperps

c-----calculation of thermal correction
      vt = k4*dsqrt(2.0d0*t_ar(1)/mass_ar(1))
      vt_dc = vt / c

      xe0 = 1.d0/(npar*vt_dc)
      call czetar(xe0,cz0,cz1)
      cdze = xe0*xe0*cz1
c-----electron terms
      eps_perp = 1.d0+X_ar(1)/y_ar(1)**2
      eps_par  = 1.d0-X_ar(1)*cdze
      eps_xy   = x_ar(1)/y_ar(1) 

c-----ion terms
      do i=2,nbulk
         eps_perp = eps_perp-x_ar(i)
         eps_par = eps_par-x_ar(i)
      enddo 
 
c-----coefficients of dispersion function
      P_0 = eps_par*((npars-eps_perp)**2 - eps_xy**2)
      P_2 = (eps_perp+eps_par)*(npars-eps_perp)+eps_xy**2
      P_4 = eps_perp
      P_6 = -(3.d0/8.d0)*(X_ar(1)/Y_ar(1)**4)*vt_dc**2
      

      do i=2,nbulk
         vt = k4*dsqrt(2.d0*t_ar(i)/mass_ar(i))
         vt_dc = vt/c
         P_6 = P_6 - (3.d0/2.d0)*X_ar(i)*vt_dc**2
      enddo
      P_6 = P_6*nosigma !zeros h.o.t.
c-----dispersion relation
      D_0 = P_6 * nperp6 + P_4 * nperp4 + P_2 * nperps + P_0

c-----dielectric tensor elements
      reps = dcmplx(0.d0,0.d0)
      reps(1,1) = dcmplx(eps_perp,0.d0)
      reps(1,2) = dcmplx(0.d0,eps_xy)
      reps(2,1) = -reps(1,2)
      reps(2,2) = reps(1,1)
      reps(3,3) = dcmplx(eps_par,0.d0)
      
      return
      end

     
      subroutine sfrnk_dervs(u,wf,dddz,dddr,dddphi,
     &                       dddcnz,dddcnr,dddcm,dddw)
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions_nml.i'
c-----input
      real*8 
     &u(6),
     &wf
c-----output
      real*8 
     &dddz,dddr,dddphi,
     &dddcnz,dddcnr,dddcm,
     &dddw
c-----externals
      real*8 b,cn,gamma1,x,y,tempe,
     &dxdr,dxdz,dxdphi,dydr,dydz,dydphi,dtempdr,dtempdz

c-----locals
      real*8  
     &z,r,phi,cnz,cnr,cm,
     &gam1,ds,dc,dcn,
     &xe0,
     &npar,nperp,npars,nperps,nperp4,nperp6,
     &eps_perp,eps_par,eps_xy,
     &k4,c,vt,vt_dc,P_0,P_2,P_4,P_6,D_0,
     &d_nperp_dr,d_nperp_dz,d_nperp_dphi,
     &d_npar_dr,d_npar_dz,d_npar_dphi, 
     &d_nperp_dcnr,d_nperp_dcnz,d_nperp_dcm,  
     &d_npar_dcnr,d_npar_dcnz,d_npar_dcm,  
     &d_nperp_dw,d_npar_dw,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

      real*8 :: cdze, cz0, cz1,nosigma
      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar
      real*8, dimension(1:nbulk) :: d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar 
      real*8, dimension(1:nbulk) :: d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar
      real*8, dimension(1:nbulk) :: d_t_dr_ar,d_t_dz_ar
      integer i

      nosigma=0.d0
      z=u(1)
      r=u(2)
      phi=u(3)
      cnz=u(4)
      cnr=u(5)
      cm=u(6)

      bmod=b(z,r,phi)  
      gam1=gamma1(z,r,phi,cnz,cnr,cm)
      ds=dsin(gam1)
      dc=dcos(gam1)
      dcn=cn(r,cnz,cnr,cm)
      npar=dcn*dc
      nperp=dcn*ds

      npars = npar**2
      nperps = nperp**2
      nperp4 = nperps**2
      nperp6 = nperp4*nperps
      
      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
         t_ar(i) = tempe(z,r,phi,i)*1.d3 !eV
         d_x_dr_ar(i)= dxdr(z,r,phi,i)
         d_x_dz_ar(i)= dxdz(z,r,phi,i)
         d_x_dphi_ar(i)=dxdphi(z,r,phi,i)
         d_y_dr_ar(i)= dydr(z,r,phi,i)
         d_y_dz_ar(i)= dydz(z,r,phi,i)
         d_y_dphi_ar(i)=dydphi(z,r,phi,i)
         d_t_dr_ar(i)=dtempdr(z,r,phi,i)*1.d3
         d_t_dz_ar(i)=dtempdz(z,r,phi,i)*1.d3
      enddo
      
      
c-----calculation of thermal correction
      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm v

      vt = k4*dsqrt(2.0d0*t_ar(1)/dmas(1))
      vt_dc = vt / c

      xe0 = 1/(npar*vt_dc)
      call czetar(xe0,cz0,cz1)
      cdze = xe0*xe0*cz1

      eps_perp = 1.d0+X_ar(1)/Y_ar(1)**2
      eps_par  = 1.d0-X_ar(1)*cdze
      eps_xy   = X_ar(1)/Y_ar(1)
      
      do i=2,nbulk
         eps_perp=eps_perp - x_ar(i)
         eps_par=eps_par - x_ar(i)
      enddo
      
      P_0=eps_par*((npars-eps_perp)**2 - eps_xy**2)
      P_2=(eps_perp+eps_par)*(npars-eps_perp)+eps_xy**2
      P_4=eps_perp
      P_6 = -(3.d0/8.d0)*(X_ar(1)/Y_ar(1)**4)*vt_dc**2

      do i=2,nbulk
         vt = k4*dsqrt(2.0d0*t_ar(i)/dmas(i))
         vt_dc = vt/c
         P_6  = P_6 - (3.d0/2.d0)*X_ar(i)*vt_dc**2
      enddo
      P_6 = P_6*nosigma !zeros h.o.t.
      call dnd(z,r,phi,cnz,cnr,cm,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw)

c-----calculate derivatives of eps elements wrt r,z,phi,cnr,cnz,cm
      call d_eps_d_rzphi_nrnznm_sfrnk(
     &wf,x_ar,y_ar,t_ar,xe0,cz0,cz1,npar,r,
     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar, 
     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar,
     &d_t_dr_ar,d_t_dz_ar,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)


c-----calculate derivatives of 
      call d_P0246_d_rzphi_nrnzcm(
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm, 
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)

      call d_P_6_d_rzphi_w(r,z,phi,wf,
     &d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw)

      dddr= d_P_6_dr*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dr+
     &      d_P_4_dr*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dr+
     &      d_P_2_dr*nperps +2.d0*P_2 * nperp*d_nperp_dr+
     &      d_P_0_dr

      dddz= d_P_6_dz*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dz+
     &      d_P_4_dz*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dz+
     &      d_P_2_dz*nperps +2.d0*P_2 * nperp*d_nperp_dz+
     &      d_P_0_dz 
    
      dddphi= d_P_6_dphi*nperp6  +6.d0*P_6 * nperp4*nperp*d_nperp_dphi+
     &        d_P_4_dphi*nperp4 +4.d0*P_4 * nperps*nperp*d_nperp_dphi+
     &        d_P_2_dphi*nperps +2.d0*P_2 * nperp*d_nperp_dphi+
     &        d_P_0_dphi

      dddcnr= d_P_6_dcnr*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dcnr+
     &      d_P_4_dcnr*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcnr+
     &      d_P_2_dcnr*nperps   +2.d0*P_2 * nperp*d_nperp_dcnr+
     &      d_P_0_dcnr

      dddcnz= d_P_6_dcnz*nperp6 +6.d0*P_6 * nperp4*nperp*d_nperp_dcnz+
     &      d_P_4_dcnz*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcnz+
     &      d_P_2_dcnz*nperps   +2.d0*P_2 * nperp*d_nperp_dcnz+
     &      d_P_0_dcnz

      dddcm= d_P_6_dcm*nperp6  +6.d0*P_6 * nperp4*nperp*d_nperp_dcm+
     &      d_P_4_dcm*nperp4   +4.d0*P_4 * nperps*nperp*d_nperp_dcm+
     &      d_P_2_dcm*nperps   +2.d0*P_2 * nperp*d_nperp_dcm+
     &      d_P_0_dcm
    
      dddw= d_P_6_dw* nperp6   +6.d0*P_6 * nperp4*nperp*d_nperp_dw+
     &      d_P_4_dw*nperp4    +4.d0*P_4 * nperps*nperp*d_nperp_dw+
     &      d_P_2_dw*nperps    +2.d0*P_2 * nperp*d_nperp_dw+
     &      d_P_0_dw

      return
      end

      subroutine d_eps_d_rzphi_nrnznm_sfrnk(
     &wf,x_ar,y_ar,t_ar,xe0,cz0,cz1,npar,r,
     &d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar, 
     &d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar, 
     &d_t_dr_ar,d_t_dz_ar,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)
c---------------------------------------------------------x
c     calculates derivatives of dielectric tensor
c     with respect to r,z,phi,cnr,cnz,cm
c------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'
      
c-----input
      real*8 wf,r,npar,  
     &x_ar(*),y_ar(*), t_ar(*), xe0, 
     &d_x_dr_ar(*),d_x_dz_ar(*),d_x_dphi_ar(*), 
     &d_y_dr_ar(*),d_y_dz_ar(*),d_y_dphi_ar(*),
     &d_t_dr_ar(*),d_t_dz_ar(*)
      real*8 cz0,cz1
c-----output
      real*8
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----locals
      integer i
      real*8, dimension(1:nbulk) :: d_eps_par_d_x_ar,d_eps_par_d_y_ar
      real*8, dimension(1:nbulk) :: d_eps_perp_d_x_ar,d_eps_perp_d_y_ar
      real*8, dimension(1:nbulk) :: d_eps_xy_d_x_ar,d_eps_xy_d_y_ar
      real*8, dimension(1:nbulk) :: d_eps_par_d_t_ar
      real*8, dimension(1:nbulk) :: d_eps_par_d_npar
      
      real*8 dxdw,dydw,cz0r,cz1r,cz2

      cz0r = cz0
      cz1r = cz1
      cz2 = -2*(cz0r+xe0*cz1r)
      do i=1,nbulk
         if (i.eq.1) then
c-----------electron term        
            d_eps_perp_d_x_ar(i)    =  1.d0/y_ar(i)**2  
            d_eps_par_d_x_ar(i)     =  -xe0**2*cz1r
            d_eps_xy_d_x_ar(i)      =  1.d0/y_ar(i)

            d_eps_perp_d_y_ar(i)    =  -2.d0*x_ar(i)/y_ar(i)**3   
            d_eps_par_d_y_ar(i)     =   0    
            d_eps_xy_d_y_ar(i)      =   -x_ar(i)/y_ar(i)**2

            d_eps_par_d_t_ar(i)     =   
     &      (x_ar(i)/T_ar(i))*xe0**2*cz1r
     &      +0.5*(x_ar(i)/T_ar(i))*xe0**3*cz2
       
            d_eps_par_d_npar(i)     = 
     &      2.d0*(x_ar(i)/npar)*xe0**2*cz1r     
     &      + (x_ar(i)/npar)*xe0**3*cz2

         else 
c-----------ion terms
            d_eps_perp_d_x_ar(i)    =   -1.d0
            d_eps_par_d_x_ar(i)     =   -1.d0
            d_eps_xy_d_x_ar(i)      =   0.d0

            d_eps_perp_d_y_ar(i)    =   0.d0 
            d_eps_par_d_y_ar(i)     =   0.d0
            d_eps_xy_d_y_ar(i)      =   0.d0

            d_eps_par_d_t_ar(i)     =   0.d0
            d_eps_par_d_npar(i)     =   0.d0
           
         endif
      enddo !nbulk

      d_eps_par_dr=0.d0
      d_eps_par_dz=0.d0
      d_eps_par_dphi=0.d0
      d_eps_perp_dr=0.d0
      d_eps_perp_dz=0.d0
      d_eps_perp_dphi=0.d0
      d_eps_xy_dr=0.d0
      d_eps_xy_dz=0.d0
      d_eps_xy_dphi=0.d0
      d_eps_par_dcnr=0.d0
      d_eps_par_dcnz=0.d0
      d_eps_par_dcm=0.d0
      d_eps_perp_dcnr=0.d0
      d_eps_perp_dcnz=0.d0
      d_eps_perp_dcm=0.d0
      d_eps_xy_dcnr=0.d0
      d_eps_xy_dcnz=0.d0
      d_eps_xy_dcm=0.d0
      d_eps_par_dw=0.d0
      d_eps_perp_dw=0.d0
      d_eps_xy_dw=0.d0

      do i=1,nbulk

         d_eps_par_dr=d_eps_par_dr + d_eps_par_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_par_d_y_ar(i)*d_y_dr_ar(i)+
     &                               d_eps_par_d_t_ar(i)*d_t_dr_ar(i)      

         d_eps_par_dz=d_eps_par_dz + d_eps_par_d_x_ar(i)*d_x_dz_ar(i)+
     &                               d_eps_par_d_y_ar(i)*d_y_dz_ar(i)+
     &                               d_eps_par_d_t_ar(i)*d_t_dz_ar(i)

         d_eps_par_dphi=d_eps_par_dphi+
     &                  d_eps_par_d_x_ar(i)*d_x_dphi_ar(i)+
     &                  d_eps_par_d_y_ar(i)*d_y_dphi_ar(i)


         d_eps_perp_dr=d_eps_perp_dr+d_eps_perp_d_x_ar(i)*d_x_dr_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dr_ar(i)

         d_eps_perp_dz=d_eps_perp_dz+d_eps_perp_d_x_ar(i)*d_x_dz_ar(i)+
     &                               d_eps_perp_d_y_ar(i)*d_y_dz_ar(i)

         d_eps_perp_dphi=d_eps_perp_dphi+
     &                   d_eps_perp_d_x_ar(i)*d_x_dphi_ar(i)+
     &                   d_eps_perp_d_y_ar(i)*d_y_dphi_ar(i)


         d_eps_xy_dr=d_eps_xy_dr+d_eps_xy_d_x_ar(i)*d_x_dr_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dr_ar(i)


         d_eps_xy_dz=d_eps_xy_dz+d_eps_xy_d_x_ar(i)*d_x_dz_ar(i)+
     &                           d_eps_xy_d_y_ar(i)*d_y_dz_ar(i)

         d_eps_xy_dphi=d_eps_xy_dphi+
     &                 d_eps_xy_d_x_ar(i)*d_x_dphi_ar(i)+
     &                 d_eps_xy_d_y_ar(i)*d_y_dphi_ar(i)

         d_eps_par_dcnr=0.d0
         d_eps_par_dcnz=0.d0
         d_eps_par_dcm=0.d0

         d_eps_perp_dcnr=0.d0        
         d_eps_perp_dcnz=0.d0
         d_eps_perp_dcm=0.d0
        
         d_eps_xy_dcnr=0.d0
         d_eps_xy_dcnz=0.d0
         d_eps_xy_dcm=0.d0
         
          
         dxdw=-2.d0*x_ar(i)/wf 
         dydw=-y_ar(i)/wf 
            
         d_eps_par_dw=d_eps_par_dw + d_eps_par_d_x_ar(i)*dxdw+
     &                               d_eps_par_d_y_ar(i)*dydw

         d_eps_perp_dw=d_eps_perp_dw+d_eps_perp_d_x_ar(i)*dxdw+
     &                               d_eps_perp_d_y_ar(i)*dydw

         d_eps_xy_dw=d_eps_xy_dw +   d_eps_xy_d_x_ar(i)*dxdw+
     &                               d_eps_xy_d_y_ar(i)*dydw
         
      enddo
      
      d_eps_par_dcnr = d_eps_par_d_npar(1)*br/bmod
      d_eps_par_dcnz = d_eps_par_d_npar(1)*bz/bmod
      d_eps_par_dcm  = d_eps_par_d_npar(1)*bphi/bmod/r

      return
      end

      subroutine eps_sfrnk(x_ar,y_ar,t_ar,eps_perp,eps_par,eps_xy) 
      implicit none
       include 'param.i'
       include 'one.i'
c------input
       real*8 x_ar(*),y_ar(*),t_ar(*)
c------output
       real*8  eps_perp,eps_par,eps_xy 
c------locals
       integer i
       real*8 vt,vt_dc,c,xe0,cdze,k4
       complex*8 cz1,cz0
c------externals
       real*8 dmas,npar
c-------------------------------------------------------------------
c      dielectric tesor calculations: eps_perp, eps_par, eps_xy
c      article formula 16(a-c)
c-------------------------------------------------------------------
c      thermal correction
       c = 3.d10
       k4 = 4.19d7
      vt = k4*dsqrt(2.0d0*t_ar(1))
      vt_dc = vt / c

      xe0 = 1/(npar*vt_dc)
      call czetar(xe0,cz0,cz1)
      cdze = xe0*xe0*REAL(cz1)

c      electron terms
       eps_perp=1.d0+X_ar(1)/y_ar(1)**2  
       eps_par=1.d0-x_ar(1)*cdze
       eps_xy=x_ar(1)/y_ar(1) 
c      ion terms
       do i=2,nbulk
         eps_perp=eps_perp-x_ar(i)
         eps_par=eps_par-x_ar(i)
       enddo

       return
       end
      subroutine d_P0246_d_rzphi_nrnzcm(
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw, 
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw,
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw)

c-----calculates derivatives of coefficients P0,P2,P4,P6
c     with respect to r,z,phi,nr,nz,cm
c------------------------------------------------------------
c     INPUT:
c     &npar,nperp,npars,
c     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
c     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
c     &d_npar_dz,d_npar_dr,d_npar_dphi,
c     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
c     &d_nperp_dw,d_npar_dw, 
c     &eps_perp,eps_par,eps_xy,
c     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
c     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
c     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
c     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
c     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
c     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
c     d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c     OUTPUT
c     d_P_6_dr,d_P_6_dz,d_P_6_dphi,
c     d_P_4_dr,d_P_4_dz,d_P_4_dphi,
c     d_P_2_dr,d_P_2_dz,d_P_2_dphi,
c     d_P_0_dr,d_P_0_dz,d_P_0_dphi,
c     d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
c     d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
c     d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
c     d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm,
c     d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
c-----------------------------------------------------
       implicit none
c-----input
      real*8 z,r,phi,cnz,cnr,cm ,
     &npar,nperp,npars,
     &d_nperp_dz,d_nperp_dr,d_nperp_dphi,
     &d_nperp_dcnz,d_nperp_dcnr,d_nperp_dcm,
     &d_npar_dz,d_npar_dr,d_npar_dphi,
     &d_npar_dcnz,d_npar_dcnr,d_npar_dcm,
     &d_nperp_dw,d_npar_dw,
     &eps_perp,eps_par,eps_xy,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dcnr,d_eps_par_dcnz,d_eps_par_dcm,
     &d_eps_perp_dcnr,d_eps_perp_dcnz,d_eps_perp_dcm,
     &d_eps_xy_dcnr,d_eps_xy_dcnz,d_eps_xy_dcm,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----output
      real*8
     & d_P_6_dr,d_P_6_dz,d_P_6_dphi,
     & d_P_4_dr,d_P_4_dz,d_P_4_dphi,
     & d_P_2_dr,d_P_2_dz,d_P_2_dphi,
     & d_P_0_dr,d_P_0_dz,d_P_0_dphi,
     & d_P_6_dcnr,d_P_6_dcnz,d_P_6_dcm,
     & d_P_4_dcnr,d_P_4_dcnz,d_P_4_dcm,
     & d_P_2_dcnr,d_P_2_dcnz,d_P_2_dcm,
     & d_P_0_dcnr,d_P_0_dcnz,d_P_0_dcm , 
     & d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw

c-----locals

      d_P_0_dr=d_eps_par_dr*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dr-d_eps_perp_dr)-
     &         2.d0*eps_xy*d_eps_xy_dr)

      d_P_0_dz=d_eps_par_dz*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dz-d_eps_perp_dz)-
     &         2.d0*eps_xy*d_eps_xy_dz)

      d_P_0_dphi=d_eps_par_dphi*((npars-eps_perp)**2 - eps_xy**2)+
     &eps_par*(2*(npars-eps_perp)*(2.d0*npar*d_npar_dphi-
     &                             d_eps_perp_dphi)-
     &         2.d0*eps_xy*d_eps_xy_dphi)

      d_P_2_dr=(d_eps_perp_dr+d_eps_par_dr)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dr-d_eps_perp_dr)+
     &                                    2.d0*eps_xy*d_eps_xy_dr

      d_P_2_dz=(d_eps_perp_dz+d_eps_par_dz)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dz-d_eps_perp_dz)+
     &                                    2.d0*eps_xy*d_eps_xy_dz

      d_P_2_dphi=(d_eps_perp_dphi+d_eps_par_dphi)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2*npar*d_npar_dphi-d_eps_perp_dphi)+
     &                                    2.d0*eps_xy*d_eps_xy_dphi

      d_P_4_dr=d_eps_perp_dr
      d_P_4_dz=d_eps_perp_dz
      d_P_4_dphi=d_eps_perp_dphi
c-----------------------------
      d_P_6_dr=0.d0
      d_P_6_dz=0.d0
      d_P_6_dphi=0.d0

c----------------------------------
      
      d_P_0_dcnr=d_eps_par_dcnr*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcnr-d_eps_perp_dcnr)-
     &          2.d0*eps_xy*d_eps_xy_dcnr)
     
      d_P_0_dcnz=d_eps_par_dcnz*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcnz-d_eps_perp_dcnz)-
     &          2.d0*eps_xy*d_eps_xy_dcnz)
     
      d_P_0_dcm=d_eps_par_dcm*((npars-eps_perp)**2 - eps_xy**2)+
     & eps_par*(2.d0*(npars-eps_perp)*
     &          (2.d0*npar*d_npar_dcm-d_eps_perp_dcm)-
     &          2.d0*eps_xy*d_eps_xy_dcm)
     
      d_P_2_dcnr=(d_eps_perp_dcnr+d_eps_par_dcnr)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcnr-d_eps_perp_dcnr)+
     &2.d0*eps_xy*d_eps_xy_dcnr

      d_P_2_dcnz=(d_eps_perp_dcnz+d_eps_par_dcnz)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcnz-d_eps_perp_dcnz)+
     &2.d0*eps_xy*d_eps_xy_dcnz

      d_P_2_dcm=(d_eps_perp_dcm+d_eps_par_dcm)*(npars-eps_perp)+
     &(eps_perp+eps_par)*(2.d0*npar*d_npar_dcm-d_eps_perp_dcm)+
     &2.d0*eps_xy*d_eps_xy_dcm
 
      d_P_4_dcnr=d_eps_perp_dcnr   
      d_P_4_dcnz=d_eps_perp_dcnz
      d_P_4_dcm=d_eps_perp_dcm


c-----------------------------
      d_P_6_dcnr=0.d0
      d_P_6_dcnz=0.d0
      d_P_6_dcm=0.d0 
c----------------------------------

      d_P_0_dw=d_eps_par_dw*((npars-eps_perp)**2 - eps_xy**2)+
     &          eps_par*(2.d0*(npars-eps_perp)*
     &         (2.d0*npar*d_npar_dw-d_eps_perp_dw)-
     &          2.d0*eps_xy*d_eps_xy_dw)

      d_P_2_dw=(d_eps_perp_dw+d_eps_par_dw)*(npars-eps_perp)+
     &         (eps_perp+eps_par)*(2.d0*npar*d_npar_dw-d_eps_perp_dw)+
     &         2.d0*eps_xy*d_eps_xy_dw


      d_P_4_dw=d_eps_perp_dw


      d_P_6_dw=0.d0

c----------------------------------

c      write(*,*)'d_P0246_d_rzphi_nrnzcm'
c      write(*,*)'npar,nperp,npars',npar,nperp,npars
c      write(*,*)'&d_nperp_dw,d_npar_dw',d_nperp_dw,d_npar_dw
c      write(*,*)'eps_perp,eps_par,eps_xy',eps_perp,eps_par,eps_xy
c      write(*,*)'d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw',
c     &           d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c      write(*,*)'d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw',
c     &           d_P_6_dw,d_P_4_dw,d_P_2_dw,d_P_0_dw
      return
      end

      subroutine d_P_6_d_rzphi_w(r,z,phi,wf,
     &d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw)
c-----------------------------------------------------------------
c     Calculate derivatives: d_P_0/d_r,d_P_0/d_z,d_P_0/d_phi,d_P_0/d_w
c-------------------------------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'
      include 'ions.i'

c-----input
      real*8 r,z,phi,wf
c-----output
      real*8 d_P_6_dr,d_P_6_dz,d_P_6_dphi,d_P_6_dw
c-----externals
      real*8 b,x,y,tempe,dxdr,dxdz,dxdphi,dydr,dydz,dydphi,
     &dtempdr,dtempdz
c-----locals
      integer i

      real*8 c,k4,vt,vt_dc,p,vt_dc_s,dv_dr,dv_dz,dv_dphi
     &       ,nosigma 
      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar,
     &        d_x_dr_ar,d_x_dz_ar,d_x_dphi_ar,
     &        d_y_dr_ar,d_y_dz_ar,d_y_dphi_ar,
     &        d_t_dr_ar,d_t_dz_ar,d_t_dphi_ar

      bmod=b(z,r,phi)

      nosigma=0.d0
      do i=1,nbulk
         x_ar(i) = x(z,r,phi,i)
         y_ar(i) = y(z,r,phi,i)
         t_ar(i) = tempe(z,r,phi,i)*1.d3 !eV
         d_x_dr_ar(i)= dxdr(z,r,phi,i)
         d_x_dz_ar(i)= dxdz(z,r,phi,i)
         d_x_dphi_ar(i)=dxdphi(z,r,phi,i)
         d_y_dr_ar(i)= dydr(z,r,phi,i)
         d_y_dz_ar(i)= dydz(z,r,phi,i)
         d_y_dphi_ar(i)=dydphi(z,r,phi,i) 
         d_t_dr_ar(i)= dtempdr(z,r,phi,i)*1.d3 
         d_t_dz_ar(i)= dtempdz(z,r,phi,i)*1.d3 
         d_t_dphi_ar(i)=0.d0
      enddo

      d_P_6_dr=0.d0
      d_P_6_dz=0.d0
      d_P_6_dphi=0.d0
      d_P_6_dw=0.d0

      c = 3.0d10                    !speed of light
      k4 = 4.19d7                   !constant for elec therm v
c-----electron terms
      vt =k4*dsqrt(2.0d0*t_ar(1)/dmas(1)) !vt= sqrt(2.0*kT[eV]/mass) 
                                             !cm/sec 
      vt_dc = vt / c
      vt_dc_s= vt_dc*vt_dc

      p=k4/dsqrt(2.d0*t_ar(1)*dmas(1))

      dv_dr=p* d_t_dr_ar(1)
      dv_dz=p* d_t_dz_ar(1)
      dv_dphi=0.d0

      d_P_6_dr=- (3.d0/8.d0)*(
     & (d_x_dr_ar(1)/y_ar(1)**4)*vt_dc_s -
     & (4.d0*x_ar(1)/y_ar(1)**5)*d_y_dr_ar(1)*vt_dc_s +
     &  (x_ar(1)/y_ar(1)**4)*(2.d0*vt_dc*(dv_dr/c)))

      d_P_6_dz=- (3.d0/8.d0)*(
     & (d_x_dz_ar(1)/y_ar(1)**4)*vt_dc_s -
     & (4.d0*x_ar(1)/y_ar(1)**5)*d_y_dz_ar(1)*vt_dc_s +
     &  (x_ar(1)/y_ar(1)**4)*(2.d0*vt_dc*(dv_dz/c)))

      d_P_6_dw=- (3.d0/8.d0)*
     &     (2.d0*x_ar(1)/y_ar(1)**4)*vt_dc_s/wf

c-----ion terms
      do i=2,nbulk
         vt =k4*dsqrt(2.0d0*t_ar(i)/dmas(i)) !vt= sqrt(2.0*kT[eV]/mass) 
                                                !cm/sec 
         vt_dc = vt / c
         vt_dc_s= vt_dc*vt_dc

         p=k4/dsqrt(2.d0*t_ar(i)*dmas(i))

         dv_dr=p* d_t_dr_ar(i)
         dv_dz=p* d_t_dz_ar(i)
         dv_dphi=0.d0

         d_P_6_dr = d_P_6_dr - (3.d0/2.d0)*(
     &                         d_x_dr_ar(i)*vt_dc_s +
     &                         x_ar(i)*(2.d0*vt_dc*(dv_dr/c)))

         
         d_P_6_dz = d_P_6_dz - (3.d0/2.d0)*(
     &                         d_x_dz_ar(i)*vt_dc_s +
     &                         x_ar(i)*(2.d0*vt_dc*(dv_dz/c)))

         d_P_6_dw = d_P_6_dw -(3.d0/2.d0)*
     &                         (-2.d0*x_ar(i)*vt_dc_s)/wf
      enddo

      d_P_6_dr = d_P_6_dr*nosigma
      d_P_6_dz = d_P_6_dz*nosigma
      d_P_6_dw = d_P_6_dw*nosigma
      
      return
      end

      
      subroutine sfrnk_numerical_deriv_of_eps(r,z,phi,cnr,cnz,cm,
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,    
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw)
c----------------------------------------------------------
c     calculates numerical derivatives of Bonoli eps
c     with respectr to r,z,phi,cnr,cnz,cm,
c-----------------------------------------------
      implicit none
      include 'param.i'
      include 'one.i'

c-----input
      real*8 r,z,phi,cnr,cnz,cm
c-----output
      real*8  
     &d_eps_par_dr,d_eps_par_dz,d_eps_par_dphi,
     &d_eps_perp_dr,d_eps_perp_dz,d_eps_perp_dphi,
     &d_eps_xy_dr,d_eps_xy_dz,d_eps_xy_dphi,
     &d_eps_par_dw,d_eps_perp_dw,d_eps_xy_dw
c-----externals
      real*8 b,x,y,tempe
c-----locals
      real*8 step,
     &rp,rm,zp,zm,phip,phim,cnrp,cnrm,cnzp,cnzm,cmp,cmm,
     &eps_perp_p,eps_par_p,eps_xy_p,
     &eps_perp_m,eps_par_m,eps_xy_m,
     &hw,hfrqnc,frqncpl,df,frqncmn

      integer i
      real*8, dimension(1:nbulk) :: x_ar,y_ar,t_ar,vp,wp


      step=1.d-6

      hw=step*1.0d0
      hfrqnc=hw
      hfrqnc=hfrqnc*frqncy

      rp=r+step
      bmod=b(z,rp,phi)
      do i=1,nbulk
         x_ar(i) = x(z,rp,phi,i)
         y_ar(i) = y(z,rp,phi,i)
         t_ar(i) = tempe(z,rp,phi,i)
      enddo

      call eps_sfrnk(x_ar,y_ar,t_ar,eps_perp_p,eps_par_p,eps_xy_p)

      rm=r-step
      bmod=b(z,rm,phi)
      do i=1,nbulk
         x_ar(i) = x(z,rm,phi,i)
         y_ar(i) = y(z,rm,phi,i)
         t_ar(i) = tempe(z,rm,phi,i)
      enddo
      call eps_sfrnk(x_ar,y_ar,t_ar,eps_perp_m,eps_par_m,eps_xy_m)

      d_eps_par_dr=(eps_par_p-eps_par_m)/(2.d0*step)    
      d_eps_perp_dr=(eps_perp_p-eps_perp_m)/(2.d0*step)
      d_eps_xy_dr=(eps_xy_p-eps_xy_m)/(2.d0*step)
c--------------------------------------------------------------
      zp=z+step
      bmod=b(zp,r,phi)
      do i=1,nbulk
         x_ar(i) = x(zp,r,phi,i)
         y_ar(i) = y(zp,r,phi,i)
         t_ar(i) = tempe(zp,r,phi,i)
      enddo
      call eps_sfrnk(x_ar,y_ar,eps_perp_p,eps_par_p,eps_xy_p)

      zm=z-step
      bmod=b(zm,r,phi)
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
         t_ar(i) = tempe(zm,r,phi,i)
      enddo
      call eps_sfrnk(x_ar,y_ar,t_ar,eps_perp_m,eps_par_m,eps_xy_m)
    
      d_eps_par_dz=(eps_par_p-eps_par_m)/(2.d0*step)    
      d_eps_perp_dz=(eps_perp_p-eps_perp_m)/(2.d0*step)
      d_eps_xy_dz=(eps_xy_p-eps_xy_m)/(2.d0*step)

      bmod=b(z,r,phi)
      do i=1,nbulk
	  vp(i)=v(i)
	  wp(i)=w(i)
      enddo

      frqncpl=frqncy+hfrqnc
      df=frqncy/frqncpl
      do i=1,nbulk
	  v(i)=vp(i)*df* df
	  w(i)=wp(i)*df
      enddo
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
         t_ar(i) = tempe(zm,r,phi,i)
      enddo
      call eps_sfrnk(x_ar,y_ar,t_ar,eps_perp_p,eps_par_p,eps_xy_p)
  
      frqncmn=frqncy-hfrqnc
      df=frqncy/frqncmn
      do i=1,nbulk
	  v(i)=vp(i)*df*df
	  w(i)=wp(i)*df
      enddo
      do i=1,nbulk
         x_ar(i) = x(zm,r,phi,i)
         y_ar(i) = y(zm,r,phi,i)
         t_ar(i) = tempe(zm,r,phi,i)
      enddo
      call eps_sfrnk(x_ar,y_ar,t_ar,eps_perp_m,eps_par_m,eps_xy_m)

      d_eps_par_dw=(eps_par_p-eps_par_m)/(2.d0*hfrqnc)  
      d_eps_perp_dw=(eps_perp_p-eps_perp_m)/(2.d0*hfrqnc)
      d_eps_xy_dw=(eps_xy_p-eps_xy_m)/(2.d0*hfrqnc)

      do i=1,nbulk
	  v(i)=vp(i)
	  w(i)=wp(i)
      enddo

      return
      end


      subroutine czetar(X,CZ0,CZ1)
C
C ==========================================================
C SFRNK copied from TORIC credit to M. Brambilla
C This subroutine evaluates the Plasma Dispersion Function
C 	and its first derivative, for real argument only
C       sfrnk modified to only output real part for speed.
C
C	Input:  X      - argument (real)
C
C	Output: CZ0    - Z(x)     (real)
C	        CZ1    - Z'(x)    (real)
C ==========================================================
C
      implicit none

      real*8, intent(out) :: cz0, cz1
      real*8, intent(in) :: x
      real*8 :: x2, zsum, terme, zn, term, tmult, z1
      integer :: n

C
C ==========================================================
C
      X2    = X*X

      IF(ABS(X).LE.5.7)  THEN
C ----------------------------------------------------------
C   Power series expansion
C
         ZSUM  = 1.d0
         TERME = 1.d0
         N     = 0

   10    N     = N+1
         ZN    = REAL(N)
         TERME = TERME*X2/ZN
         TERM  = TERME/(ZN+ZN+1.d0)
         ZSUM  = ZSUM + TERM
         IF(ABS(TERM/ZSUM).GT.1.E-14)  GOTO 10

         CZ0 = -2.d0*X*ZSUM*EXP(-X2)
         CZ1 = -2.d0*(1.d0 + X*CZ0)

      ELSE
C ----------------------------------------------------------
C  Asymptotic expansion
C
         N     = 0
         ZN    = 1.d0
         TERME = 1.d0
         TMULT = 0.5/X2
         ZSUM  = 1.d0

   20    TERME = ZN*TERME*TMULT
         N = N+1
         ZN = ZN+2.d0
         TERM = ZN*TERME
         ZSUM = ZSUM + TERM
         IF(ZN.LT.20.d0 .AND. ABS(TERM/ZSUM).GT.1.E-14)  GOTO 20

         Z1 = ZSUM/X2 

         CZ1 = Z1
         CZ0 = -(1.d0 + 0.5*CZ1)/X

      END IF

      return
      end subroutine czetar
