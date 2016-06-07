	program mc_shms_single

C+______________________________________________________________________________
!
! Monte-Carlo of SHMS spectrometer using uniform illumination.
!   This version uses TRANSPORT right-handed coordinate system.
C-______________________________________________________________________________

	implicit none

	include 'shms/struct_shms.inc'
	include 'spectrometers.inc'
	include 'constants.inc'

C HBOOK/NTUPLE common block and parameters.
	integer*4	pawc_size
	parameter	(pawc_size = 1000000)
	common		/pawc/ hbdata(pawc_size)
	integer*4	hbdata
	character*8	hut_nt_names(28)/
     >			'hsxfp', 'hsyfp', 'hsxpfp', 'hsypfp',
     >			'hsztari','hsytari', 'hsdeltai', 'hsyptari', 'hsxptari',
     >			'hsztar','hsytar', 'hsdelta', 'hsyptar', 'hsxptar', 
     >                  'hsxtari','yrast','xsnum','ysnum','xsieve'
     >                  ,'ysieve'
     >                  ,'evtype','normfac','xsn'
     >                  ,'theta','eprime','ecalc','Q2','W'
     >                  /
	real*4		hut(28)

	character*8	spec_nt_names(58)/
     >			's_hb1_x', 's_hb1_y','s_hb2_x', 's_hb2_y','s_hb3_x', 's_hb3_y','s_hb4_x', 's_hb4_y', 's_q1_x', 's_q1_y', ! 10
     >                  's_q2_x', 's_q2_y', 's_q3_x', 's_q3_y', !14
     >                  's_d1_x', 's_d1_y', 's_d1f_x', 's_d1f_y', !18
     >                  's_dme_x', 's_dme_y', 's_dm1_x', 's_dm1_y', !22
     >                  'S_dm2_x', 's_dm2_y', 'S_dm3_x','s_dm3_y', !26
     >                  'S_dm4_x', 's_dm4_y', 'S_dm5_x','s_dm5_y', !30
     >                  'S_dm6_x', 's_dm6_y', 'S_dm7_x','s_dm7_y', !34
     >                  'S_dmex_x', 's_dmex_y', 's_dex_x', 's_dex_y', !38
     >                  's_dc1_x', 's_dc1_y', 's_dc2_x', 's_dc2_y', !42
     >                  's_s1_x', 's_s1_y', 's_s2_x', 's_s2_y', !46
     >                  's_cal_x', 's_cal_y', 's_fcal_x', 's_fcal_y',
     >                  'sxfp', 'syfp', 
     >                  'sdelta', 'sxptar', 'syptar', 'sxcoll', 
     >                  'sycoll', 'sflag'/
	real*4          spec(58)
c
	real*8 xs_num,ys_num,xc_sieve,yc_sieve
	real*8 xsfr_num,ysfr_num,xc_frsieve,yc_frsieve
        logical use_front_sieve /.false./
c
        common /sieve_info/  xs_num,ys_num,xc_sieve,yc_sieve
     > ,xsfr_num,ysfr_num,xc_frsieve,yc_frsieve,use_front_sieve


C Local declarations.
	integer*4	i,
     >			chanin	/1/,
     >			chanout	/2/,
     >			n_trials,trial,
     >			tmp_int

	integer*4 Itrial                        ! TH - add this for gfortran: forces integer type cast
	logical*4	iss

	real*8 th_nsig_max                      ! TH - add this for gfortran
	parameter(th_nsig_max=3.0d0)            !max #/sigma for gaussian ran #s

C Event limits, topdrawer limits, physics quantities
	real*8 gen_lim(6)			!M.C. phase space limits.
	real*8 gen_lim_up(3)
	real*8 gen_lim_down(3)

	real*8 cut_dpp,cut_dth,cut_dph,cut_z	!cuts on reconstructed quantities
	real*8 th_ev,cos_ev,sin_ev		!cos and sin of event angle

	real*8 cos_ts,sin_ts			!cos and sin of spectrometer angle
	real*8 x,y,z,dpp,dxdz,dydz,t1,t2,t3,t4	!temporaries

	real*8 x_a,y_a,z_a,dydz_a,dif_a,dydz_aa,dif_aa   ! TH - for target aperture check
	real*8 musc_targ_len			!target length for multiple scattering
	real*8 m2				!particle mass squared.
	real*8 rad_len_cm			!conversion r.l. to cm for target
	real*8 pathlen				!path length through spectrometer.
	logical*4 ok_spec			!indicates whether event makes it in MC
	integer*4 hit_calo                      !flag for hitting the calorimeter

C Initial and reconstructed track quantities.
	real*8 dpp_init,dth_init,dph_init,xtar_init,ytar_init,ztar_init
	real*8 dpp_recon,dth_recon,dph_recon,ztar_recon,ytar_recon
	real*8 x_fp,y_fp,dx_fp,dy_fp		!at focal plane
	real*8 p_spec,th_spec			!spectrometer setting
	real*8 resmult

C Control flags (from input file)
	integer*4 p_flag			!particle identification
	logical*4 use_aer                       !Aerogel usage flag
	logical*4 ms_flag
	logical*4 wcs_flag
	logical*4 cer_flag
	logical*4 vac_flag

	common /hutflag/ cer_flag,vac_flag
C Hardwired control flags.
	logical*4 hut_ntuple	/.true./
        logical*4 spec_ntuple   /.false./
	logical*4 decay_flag	/.false./

	real*8	dpp_var(2),dth_var(2),dph_var(2),ztg_var(2)
	integer*8	stime,etime

	character*132	str_line
C Local  spectrometer varibales
	real*8 x_s,y_s,z_s
	real*8 dxdz_s,dydz_s,dpp_s
c carbon cross section
	integer doing_carbon /0/
	real*8 mass_tar,theta_pol,ebeam,eprime
        real*8 car_density /2.2/ ! g/cm3
	REAL*8 q2_vertex,W_vertex
	real event_type
	real*8 Q_E, N_A,lumin,ep_min,ep_max,domega,denergy
	real*8 sig_elastic,sig_inelastic
	real*8 cur,normfac,thick
	real*8 theta_recon,eprime_recon,eprime_calc
        PARAMETER (Q_E = 1.602d00)            !e- charge in uCoul (*1E-13)
        PARAMETER (N_A = 6.022d00)            !Avogadro's number (*1E+23)
	real*8 hbarcsq,sig_mott
        real*8 thrown_wt
c
C Function definitions.

	integer*4	last_char
	logical*4	rd_int,rd_real
	real*8          grnd,gauss1
	INTEGER irnd
	REAL rnd(99)
        integer      itime,ij
        character	timestring*30

        character*80 rawname, filename
	real*4  secnds,zero

	parameter(zero=0.0)

        integer iquest
        common/quest/iquest(100)

	save		!Remember it all!

C ================================ Executable Code =============================

C Initialize
C using SIMC unstructured version
C
	shmsSTOP_trials	= 0
	shmsSTOP_HB_in	= 0
        shmsSTOP_HB_men = 0
	shmsSTOP_HB_mex	= 0
	shmsSTOP_HB_out	= 0	
	shmsSTOP_targ_hor	= 0
	shmsSTOP_targ_vert	= 0
	shmsSTOP_targ_oct	= 0
	shmsSTOP_slit_hor	= 0
	shmsSTOP_slit_vert	= 0
	shmsSTOP_slit_oct	= 0
	shmsSTOP_Q1_in	= 0
	shmsSTOP_Q1_men	= 0
	shmsSTOP_Q1_mid	= 0
	shmsSTOP_Q1_mex	= 0
	shmsSTOP_Q1_out	= 0
	shmsSTOP_Q2_in	= 0
	shmsSTOP_Q2_men	= 0
	shmsSTOP_Q2_mid	= 0
	shmsSTOP_Q2_mex	= 0
	shmsSTOP_Q2_out	= 0
	shmsSTOP_Q3_in	= 0
	shmsSTOP_Q3_men	= 0
	shmsSTOP_Q3_mid	= 0
	shmsSTOP_Q3_mex	= 0
	shmsSTOP_Q3_out	= 0
c	shmsSTOP_Q3_out1	= 0
c	shmsSTOP_Q3_out2	= 0
c	shmsSTOP_Q3_out3	= 0
c	shmsSTOP_Q3_out4	= 0
c	shmsSTOP_Q3_out5	= 0
c	shmsSTOP_Q3_out6	= 0
	shmsSTOP_D1_in	= 0
        shmsSTOP_D1_flr = 0
	shmsSTOP_D1_men = 0
	shmsSTOP_D1_mid1 = 0
	shmsSTOP_D1_mid2 = 0
	shmsSTOP_D1_mid3 = 0
	shmsSTOP_D1_mid4 = 0
	shmsSTOP_D1_mid5 = 0
	shmsSTOP_D1_mid6 = 0
	shmsSTOP_D1_mid7 = 0
	shmsSTOP_D1_mex = 0
	shmsSTOP_D1_out	= 0
	shmsSTOP_BP_in  = 0
	shmsSTOP_BP_out = 0
	shmsSTOP_hut	= 0
	shmsSTOP_dc1	= 0
	shmsSTOP_dc2	= 0
	shmsSTOP_s1	= 0
	shmsSTOP_s2	= 0
	shmsSTOP_s3	= 0
	shmsSTOP_cal	= 0
	shmsSTOP_successes	= 0
	stop_id = 0

C Open setup file.

	write(*,*)'Enter input filename (assumed to be in infiles dir)'
	read(*,1968) rawname
 1968	format(a)
	filename = 'infiles/'//rawname(1:last_char(rawname))//'.inp'
	print *,filename,'opened'
	open(unit=chanin,status='old',file=filename)

C Initialize HBOOK/NTUPLE if used.
	if (hut_ntuple) then
	  call hlimit(pawc_size)
	  filename = 'worksim/'//rawname(1:last_char(rawname))//'.rzdat'

cmkj          iquest(10) = 256000
cmkj	  iquest(10) = 510000
! see for example
!   http://wwwasd.web.cern.ch/wwwasd/cgi-bin/listpawfaqs.pl/7
! the file size is limited to ~260M no matter how I change iquest !
cmkj	  call hropen(30,'HUT',filename,'NQ',4096,i) !CERNLIB
	  call hropen(30,'HUT',filename,'N',1024,i) !CERNLIB
 
	  if (i.ne.0) then
! TH - use "write" instead of "type" for gfortran. Change this everywhere below
!	    type *,'HROPEN error: istat = ',i
	    write(*,*),'HROPEN error: istat = ',i
	    stop
	  endif

	  call hbookn(1411,'HUT NTUPLE',28,'HUT',10000,hut_nt_names)
          if (spec_ntuple) then
           call hbookn(1412,'SPEC NTU',58,'HUT',10000,spec_nt_names)
          endif
	endif	   

C Open Output file.
	filename = 'outfiles/'//rawname(1:last_char(rawname))//'.out'
	open (unit=chanout,status='unknown',file=filename)

C Read in real*8's from setup file

	str_line = '!'

C Strip off header

	do while (str_line(1:1).eq.'!')
	  read (chanin,1001) str_line
	enddo

! Read data lines.

	write(*,*),str_line(1:last_char(str_line))
	iss = rd_int(str_line,n_trials)
	if (.not.iss) stop 'ERROR (ntrials) in setup!'

! Spectrometer momentum:
	read (chanin,1001) str_line
	iss = rd_real(str_line,p_spec)
	if (.not.iss) stop 'ERROR (Spec momentum) in setup!'

! Spectrometer angle:
	read (chanin,1001) str_line
	iss = rd_real(str_line,th_spec)
	if (.not.iss) stop 'ERROR (Spec theta) in setup!'
	th_spec = abs(th_spec) / degrad
	cos_ts = cos(th_spec)
	sin_ts = sin(th_spec)

! M.C. limits (half width's for dp,th,ph, full width's for x,y,z)
	do i=1,3
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_down(i) = gen_lim(i)
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	  gen_lim_up(i) = gen_lim(i)
	enddo

	do i = 4,6
	  read (chanin,1001) str_line
	  write(*,*),str_line(1:last_char(str_line))
	  iss = rd_real(str_line,gen_lim(i))
	  if (.not.iss) stop 'ERROR (M.C. limits) in setup!'
	enddo

! Cuts on reconstructed quantities
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dpp)) 
     > stop 'ERROR (CUT_DPP) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dth)) 
     > stop 'ERROR (CUT_DTH) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_dph)) 
     > stop 'ERROR (CUT_DPH) in setup!'

	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,cut_z)) 
     > stop 'ERROR (CUT_Z) in setup!'

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,rad_len_cm)) 
     > stop 'ERROR (RAD_LEN_CM) in setup!'

! read in flag for particle type.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,p_flag)) 
     > stop 'ERROR: p_flag in setup file!'

! Read in flag for aerogel usage in the HMS.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: use_aer in setup file!'
	if (tmp_int.eq.1) use_aer = .true.

! Read in flag for multiple scattering.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: ms_flag in setup file!'
	if (tmp_int.eq.1) ms_flag = .true.

! Read in flag for wire chamber smearing.
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: wcs_flag in setup file!'
	if (tmp_int.eq.1) wcs_flag = .true.

! Read in flag for 1st cerenkov stack 
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: cer_flag in setup file!'
	if (tmp_int.eq.1) cer_flag = .true.

! Read in flag for vacuum pipe or helium bag 
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: vac_flag in setup file!'
	if (tmp_int.eq.1) vac_flag = .true.


! Read in flag   doing carbon cross section for elastic , 4.4 excited state and QE scattering
	read (chanin,1001,end=888) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_int(str_line,tmp_int)) 
     > stop 'ERROR: carbon xn flag in setup file!'
	doing_carbon = tmp_int

! Read in radiation length of target material in cm
	read (chanin,1001) str_line
	write(*,*),str_line(1:last_char(str_line))
	if (.not.rd_real(str_line,ebeam)) 
     > stop 'ERROR beam energy in setup!'

 888	continue
C Set particle masses.
	m2 = me2			!default to electron
	if(p_flag.eq.0) then
	  m2 = me2
	else if(p_flag.eq.1) then
	  m2 = mp2
	else if(p_flag.eq.2) then
	  m2 = md2
	else if(p_flag.eq.3) then
	  m2 = mpi2
	else if(p_flag.eq.4) then
	  m2 = mk2
	endif

C------------------------------------------------------------------------------C
C                           Top of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

	stime = secnds(zero)

! TH - use "Itrial" instead of "trial" for gfortran. Somehow the stringlib.f
! function does not type cast string to integer otherwise.
          itime=time8()
   	  call ctime(itime,timestring)
	  call srand(itime)

	do Itrial = 1,n_trials
	  if(mod(Itrial,5000).eq.0) write(*,*)'event #: ',
     >Itrial,'       successes: ',shmsSTOP_successes


	  irnd=Itrial

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians, using numbers picked above.
C Units are cm.

! TH - use a double precision for random number generation here.
	  x = gauss1(th_nsig_max) * gen_lim(4) / 6.0	!beam width
	  y = gauss1(th_nsig_max) * gen_lim(5) / 6.0	!beam height
	  z = (rand() - 0.5) * gen_lim(6)		!along target


C Pick scattering angles and DPP from independent, uniform distributions.
C dxdz and dydz in HMS TRANSPORT coordinates.

	  dpp  = rand()*(gen_lim_up(1)-gen_lim_down(1))
     &             + gen_lim_down(1)
	  dydz = rand()*(gen_lim_up(2)-gen_lim_down(2))
     &          /1000.   + gen_lim_down(2)/1000.
	  dxdz = rand()*(gen_lim_up(3)-gen_lim_down(3))
     &          /1000.   + gen_lim_down(3)/1000.
c If doing carbon elastics and quasielastics
	    event_type=0.  
	  if (doing_carbon .ge. 1) then
	    mass_tar = 12.*931.5
            cur=20. ! microAmps
	    thick=car_density*gen_lim(6)  ! g/cm2
            ep_min = p_spec*(1.+0.01*gen_lim_down(1))
	    ep_max = p_spec*(1.+0.01*gen_lim_up(1))
	    domega = (gen_lim_up(3)-gen_lim_down(3))*(gen_lim_up(2)-gen_lim_down(2))/1000./1000.
	    denergy = ep_max-ep_min
	    lumin=thick*cur/12.*N_A/Q_E*1e+10 !per fm2 per sec at 20uA 
            theta_pol = acos( (cos_ts + dydz*sin_ts)
     +                        / sqrt( 1. + dxdz**2 + dydz**2 ) )
	    thrown_wt = 1.
	    if (doing_carbon .eq. 2) then
	     thrown_wt = 2.
	     if (grnd() .le. 0.5) then ! elastic carbon scattering
 	      eprime = mass_tar*ebeam/(ebeam*(1-cos(theta_pol))+mass_tar)
	      dpp = (eprime-p_spec)/p_spec*100.
	      event_type=1.  
	     else
	      event_type=3.  
	      eprime= p_spec*(1+0.01*dpp)
	      Q2_vertex= 4.0*ebeam*eprime*sin(theta_pol/2)**2
              W_vertex= 2.*938.27*(ebeam-eprime) + (938.27)**2 - Q2_vertex
              if ( W_vertex .gt. 0)  W_vertex = sqrt(W_vertex)
	      if (  W_vertex .le. 0 ) goto 500
	      if (  ebeam-eprime .le. 0 ) goto 500
	     endif
	    endif
	    if (doing_carbon .eq. 1) then
	     event_type=3.  
	     eprime= p_spec*(1+0.01*dpp)
	     Q2_vertex= 4.0*ebeam*eprime*sin(theta_pol/2)**2
             W_vertex= 2.*938.27*(ebeam-eprime) + (938.27)**2 - Q2_vertex
             if ( W_vertex .gt. 0)  W_vertex = sqrt(W_vertex)
	     if (  W_vertex .le. 0 ) goto 500
	     if (  ebeam-eprime .le. 0 ) goto 500
            endif	    
	  endif
c



C Transform from target to SHMS (TRANSPORT) coordinates.
C Version for a spectrometer on the left-hand side: (i.e. SHMS)
	  x_s    = -y
	  y_s    = x * cos_ts - z * sin_ts
	  z_s    = z * cos_ts + x * sin_ts

C Below assumes that HMS is on the right-hand side of the beam
C line (looking downstream).
!	  xs    = -y
!	  ys    = x * cos_ts + z * sin_ts
!	  zs    = z * cos_ts - x * sin_ts

	  dpp_s  = dpp
	  dxdz_s = dxdz
	  dydz_s = dydz

C Drift back to zs = 0, the plane through the target center
	  x_s = x_s - z_s * dxdz_s
	  y_s = y_s - z_s * dydz_s
	  z_s = 0.0

C Save init values for later.
	  xtar_init = x_s
	  ytar_init = y_s
	  ztar_init = z
	  dpp_init = dpp
	  dth_init = dydz_s*1000.		!mr
	  dph_init = dxdz_s*1000.		!mr

C Calculate multiple scattering length of target
	  cos_ev = (cos_ts+dydz_s*sin_ts)/sqrt(1+dydz_s**2+dxdz_s**2)
	  th_ev = acos(cos_ev)
	  sin_ev = sin(th_ev)

C Case 1 : extended target:
C   cryo LH2 (4.0 cm diamater = 2.0 radius)
C   Liquid + 2 x 0.13mm Al (X0=8.89cm) tuna can
C
cXZ	  if (abs(gen_lim(6)).gt.3.) then   
cXZ	     musc_targ_len = (2.**2.-(z*sin_ev)**2.)**0.5-z*cos_ev ! after scattering
cXZ                                                                   ! forward or backward
cXZ	     musc_targ_len = musc_targ_len +(gen_lim(6)/2. + z)    ! before scattering
cXZ	     musc_targ_len = musc_targ_len + 0.0127/8.89         ! entrance window
cXZ	     musc_targ_len = musc_targ_len 
cXZ     >	     + 0.0138*(2./(2.**2.-(z*sin_ev)**2.)**0.5)/8.89 ! downstream window
C Case 1 : cryo target (2.65 inches diamater --> 3.37 cm radius)
C   Liquid + 5 mil Al (X0=8.89cm) beer can
C DG - put back beer can stuff...
! TH - hardcoded values may cause trouble in gfortran
	  if (abs(gen_lim(6)).gt.3.) then
	     if (((gen_lim(6)/2. - z)/cos_ev .lt. 3.37/sin_ev).
     &           and.((gen_lim(6)/2. + z)/cos_ev .gt. 3.37/sin_ev)) then !forward or backward...
*		write(*,*)'forward/backward'
		musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev

c                write(6,*) "1:  ",musc_targ_len*rad_len_cm

		musc_targ_len = musc_targ_len + .005*2.54/8.89/cos_ev
  	    else						     !side
*	       write(*,*)'side'
	       musc_targ_len = 3.37/rad_len_cm/sin_ev
	       musc_targ_len = musc_targ_len + .005*2.54/8.89/sin_ev
	    endif
C Solid target
	  else
	     musc_targ_len = abs(gen_lim(6)/2. - z)/rad_len_cm/cos_ev
	  endif

C Scattering before magnets:  Approximate all scattering as occuring AT TARGET.
C  16 mil Al scattering chamber window (X0=8.89cm)
C  15 cm air (X0=30420cm)
C spectrometer entrance window
C  20 mil Al s (X0=8.89cm)

	  musc_targ_len = musc_targ_len + .016*2.54/8.89 +
     >          15./30420. +  .020*2.54/8.89
c
	  if (ms_flag ) call musc(m2,p_spec*(1.+dpp_s/100.),
     > musc_targ_len,dydz_s,dxdz_s)

!-----------------------------------------------------------------------------
! TH - START TARGET APERTURE TESTS
! ----------------------------------------------------------------------------

! Restore xs to values at pivot. 
!	   xs = x_transp
!	   ys = y_transp
	  x_a = 0
	  y_a = 2.99 !cm
	  z_a = 57.2 !cm

	  dydz_a = (y_a-ytar_init)/(z_a-ztar_init)
	  dydz_aa = atan(dydz_a)

! Check target aperture, at about 0.572 meter
! theta_a = lower limit of aperture window
! theta_s = scattering angle (=spectrometer angle + position)
! The difference between the scattering and the limiting angle of the
! window for a given central spectrometer angle.
	  dif_a = (th_spec*1000+dth_init-dydz_aa*1000)  ! mrad

! ----------------------------------------------------------------------------
	  call mc_shms(p_spec, th_spec, dpp_s, x_s, y_s, z, 
     >          dxdz_s, dydz_s,
     >		x_fp, dx_fp, y_fp, dy_fp, m2, spec,
     >		ms_flag, wcs_flag, decay_flag, resmult, xtar_init, ok_spec, 
     >          pathlen, 5,
     >          .false.)
          if (spec_ntuple) then
	     spec(58) = stop_id
c            if (ok_spec) spec(58) =1.
            call hfn(1412,spec)
          endif

	  if (ok_spec) then !Success, increment arrays
	    dpp_recon = dpp_s
            dth_recon = dydz_s*1000.			!mr
	    dph_recon = dxdz_s*1000.			!mr
	    ztar_recon = + y_s / sin_ts 
            ytar_recon = y_s
	    if ( doing_carbon .ge. 1) then
            theta_recon = acos( (cos_ts + dydz_s*sin_ts)
     +                        / sqrt( 1. + dxdz_s**2 + dydz_s**2 ) )
	    eprime_recon= p_spec*(1+0.01*dpp_s)
	    eprime_calc=  mass_tar*ebeam/(ebeam*(1-cos(theta_recon))+mass_tar)
	    endif
c
	    if (event_type .eq. 1) then
	       call calc_elastic_sig(theta_pol*180./3.14159,sig_elastic)
	       normfac=thrown_wt*sig_elastic*lumin*domega/n_trials
	    endif
	    if (event_type .eq. 3) then
	      call calc_inelastic_sig(theta_pol,ebeam-eprime,Q2_vertex,W_vertex,sig_inelastic)
                 hbarcsq=0.389379292d0
                 sig_mott = alpha**2 / (Q2_vertex/1.d6) / tan(theta_pol/2.0d0)**2 * eprime/ebeam
                 sig_mott = sig_mott * hbarcsq * 1.0e-4  !!xsec in fm2/MeV/str/nuc
		 sig_inelastic = sig_mott*sig_inelastic
	         normfac=thrown_wt*sig_inelastic*lumin*domega*denergy/n_trials
	    endif
c

C Output NTUPLE entry.

	    if (hut_ntuple) then
	      hut(1) = x_fp
	      hut(2) = y_fp
	      hut(3) = dx_fp
	      hut(4) = dy_fp
	      hut(5) = ztar_init
	      hut(6) = ytar_init
	      hut(7) = dpp_init
	      hut(8) = dth_init/1000.
	      hut(9) = dph_init/1000.
	      hut(10) = ztar_recon
	      hut(11) = ytar_recon
	      hut(12)= dpp_recon
	      hut(13)= dth_recon/1000.
	      hut(14)= dph_recon/1000.
	      hut(15)= xtar_init
	      hut(16)= y
	      hut(17)= xs_num
	      hut(18)= ys_num
	      hut(19)= xc_sieve
	      hut(20)= yc_sieve
              if (use_front_sieve) then
	      hut(17)= xsfr_num
	      hut(18)= ysfr_num
	      hut(19)= xc_frsieve
	      hut(20)= yc_frsieve
              endif
	      hut(21)=event_type
	      hut(22)=normfac
	      if (event_type .eq. 1) hut(23)=sig_elastic
	      if (event_type .eq. 3) hut(23)=sig_inelastic
	      hut(24)=theta_recon
	      hut(25)=eprime_recon
	      hut(26)=eprime_calc
	      hut(27)=Q2_vertex*1e-6
	      hut(28)=W_vertex*1e-3
	      call hfn(1411,hut)
	    endif


C Compute sums for calculating reconstruction variances.
	    dpp_var(1) = dpp_var(1) + (dpp_recon - dpp_init)
	    dth_var(1) = dth_var(1) + (dth_recon - dth_init)
	    dph_var(1) = dph_var(1) + (dph_recon - dph_init)
	    ztg_var(1) = ztg_var(1) + (ztar_recon - ztar_init)

	    dpp_var(2) = dpp_var(2) + (dpp_recon - dpp_init)**2
	    dth_var(2) = dth_var(2) + (dth_recon - dth_init)**2
	    dph_var(2) = dph_var(2) + (dph_recon - dph_init)**2
	    ztg_var(2) = ztg_var(2) + (ztar_recon - ztar_init)**2
	  endif			!Incremented the arrays

C We are done with this event, whether GOOD or BAD.
C Loop for remainder of trials.

500	  continue

	enddo				!End of M.C. loop

C------------------------------------------------------------------------------C
C                           End of Monte-Carlo loop                            C
C------------------------------------------------------------------------------C

C Close NTUPLE file.

	call hrout(1411,i,' ')
        if (spec_ntuple) call hrout(1412,i,' ')
	call hrend('HUT')

	write (chanout,1002)
	write (chanout,1003) p_spec,th_spec*degrad
        write (chanout,1004) (gen_lim(i),i=1,6)

	write (chanout,1005) n_trials

C Indicate where particles are lost in spectrometer.

	write (chanout,1015)
     >	shmsSTOP_targ_hor,shmsSTOP_targ_vert,shmsSTOP_targ_oct,
     >	shmsSTOP_slit_hor,shmsSTOP_slit_vert,shmsSTOP_slit_oct,
     >	shmsSTOP_HB_in,shmsSTOP_HB_men,shmsSTOP_HB_mex,
     >  shmsSTOP_HB_out,shmsSTOP_Q1_in,shmsSTOP_Q1_men,
     >  shmsSTOP_Q1_mid,shmsSTOP_Q1_mex,shmsSTOP_Q1_out,
     >	shmsSTOP_Q2_in,shmsSTOP_Q2_men,shmsSTOP_Q2_mid,
     >  shmsSTOP_Q2_mex,shmsSTOP_Q2_out,
     >  shmsSTOP_Q3_in,shmsSTOP_Q3_men,shmsSTOP_Q3_mid,
     >  shmsSTOP_Q3_mex,shmsSTOP_Q3_out,
     >	shmsSTOP_D1_in,shmsSTOP_D1_flr,shmsSTOP_D1_men,
     >  shmsSTOP_D1_mid1,shmsSTOP_D1_mid2,shmsSTOP_D1_mid3,
     >  shmsSTOP_D1_mid4,shmsSTOP_D1_mid5,shmsSTOP_D1_mid6,
     >  shmsSTOP_D1_mid7,shmsSTOP_D1_mex,shmsSTOP_D1_out,
     >  shmsSTOP_BP_in, shmsSTOP_BP_out

	write (chanout,1006)
     >	shmsSTOP_trials,shmsSTOP_hut,shmsSTOP_dc1,shmsSTOP_dc2,
     >  shmsSTOP_s1,shmsSTOP_s2,shmsSTOP_s3,shmsSTOP_cal,
     >  shmsSTOP_successes,shmsSTOP_successes

C Compute reconstruction resolutions.

	if (shmsSTOP_successes.eq.0) shmsSTOP_successes=1
	t1 = sqrt(max(0.,dpp_var(2)/shmsSTOP_successes 
     > - (dpp_var(1)/shmsSTOP_successes)**2))
	t2 = sqrt(max(0.,dth_var(2)/shmsSTOP_successes 
     > - (dth_var(1)/shmsSTOP_successes)**2))
	t3 = sqrt(max(0.,dph_var(2)/shmsSTOP_successes 
     > - (dph_var(1)/shmsSTOP_successes)**2))
	t4 = sqrt(max(0.,ztg_var(2)/shmsSTOP_successes 
     > - (ztg_var(1)/shmsSTOP_successes)**2))

	write (chanout,1011) dpp_var(1)/shmsSTOP_successes,t1,
     > dth_var(1)/shmsSTOP_successes,
     >		t2,dph_var(1)/shmsSTOP_successes,t3,
     > ztg_var(1)/shmsSTOP_successes,t4

	write(6,*) shmsSTOP_trials,' Trials',shmsSTOP_successes
     > ,' Successes'
	write (6,1011) dpp_var(1)/shmsSTOP_successes,t1,
     > dth_var(1)/shmsSTOP_successes,
     >		t2,dph_var(1)/shmsSTOP_successes,t3,
     > ztg_var(1)/shmsSTOP_successes,t4

C ALL done!

	stop ' '

C =============================== Format Statements ============================

1001	format(a)
1002	format('!',/,'! Uniform illumination Monte-Carlo results')
1003	format('!',/'! Spectrometer setting:',/,'!',/,
     >g11.5,' =  P  spect (MeV)',/,
     >g11.5,' =  TH spect (deg)')

1004	format('!',/'! Monte-Carlo limits:',/,'!',/,
     >g11.5,'= GEN_LIM(1) - DP/P   (half width,% )',/,
     >g11.5,'= GEN_LIM(2) - Theta  (half width,mr)',/,
     >g11.5,'= GEN_LIM(3) - Phi    (half width,mr)',/,
     >g11.5,'= GEN_LIM(4) - HORIZ (full width of 3 sigma cutoff,cm)',/,
     >g11.5,'= GEN_LIM(5) - VERT  (full width of 3 sigma cutoff,cm)',/,
     >g11.5,'= GEN_LIM(6) - Z      (Full width,cm)')

!inp     >	,/,
!inp     >	g18.8,' =  Hor. 1/2 gap size (cm)',/,
!inp     >	g18.8,' =  Vert. 1/2 gap size (cm)')

1005	format('!',/,'! Summary:',/,'!',/,
!     >	i,' Monte-Carlo trials:')
     > i11,' Monte-Carlo trials:')

1006	format(i11,' Initial Trials',/
     >i11,' Trials made it to the hut',/
     >i11,' Trial cut in dc1',/
     >i11,' Trial cut in dc2',/
     >i11,' Trial cut in s1',/
     >i11,' Trial cut in s2',/
     >i11,' Trial cut in s3',/
     >i11,' Trial cut in cal',/
     >i11,' Trials made it thru the detectors and were reconstructed',/
     >i11,' Trials passed all cuts and were histogrammed.',/
     >)

!1008	format(8i)
!1009	format(1x,i4,g,i)
!1010	format(a,i)
1011	format(
     >'DPP ave error, resolution = ',2g18.8,' %',/,
     >'DTH ave error, resolution = ',2g18.8,' mr',/,
     >'DPH ave error, resolution = ',2g18.8,' mr',/,
     >'ZTG ave error, resolution = ',2g18.8,' cm')

1012	format(1x,16i4)

1015	format(/,
     >i11,' stopped in the TARG APERT HOR',/
     >i11,' stopped in the TARG APERT VERT',/
     >i11,' stopped in the TARG APERT OCTAGON',/
     >i11,' stopped in the FIXED SLIT HOR',/
     >i11,' stopped in the FIXED SLIT VERT',/
     >i11,' stopped in the FIXED SLIT OCTAGON',/
     >i11,' stopped in HB ENTRANCE',/
     >i11,' stopped in HB MAG ENTRANCE',/
     >i11,' stopped in HB MAG EXIT',/
     >i11,' stopped in HB EXIT',/
     >i11,' stopped in Q1 ENTRANCE',/
     >i11,' stopped in Q1 MAG ENTRANCE',/
     >i11,' stopped in Q1 MIDPLANE',/
     >i11,' stopped in Q1 MAG EXIT',/
     >i11,' stopped in Q1 EXIT',/
     >i11,' stopped in Q2 ENTRANCE',/
     >i11,' stopped in Q2 MAG ENTRANCE',/
     >i11,' stopped in Q2 MIDPLANE',/
     >i11,' stopped in Q2 MAG EXIT',/
     >i11,' stopped in Q2 EXIT',/
     >i11,' stopped in Q3 ENTRANCE',/
     >i11,' stopped in Q3 MAG ENTRANCE',/
     >i11,' stopped in Q3 MIDPLANE',/
     >i11,' stopped in Q3 MAG EXIT',/
     >i11,' stopped in Q3 EXIT',/
     >i11,' stopped in D1 ENTRANCE',/
     >i11,' stopped in D1 FLARE',/
     >i11,' stopped in D1 MAG ENTRANCE',/
     >i11,' stopped in D1 MID-1',/
     >i11,' stopped in D1 MID-2',/
     >i11,' stopped in D1 MID-3',/
     >i11,' stopped in D1 MID-4',/
     >i11,' stopped in D1 MID-5',/
     >i11,' stopped in D1 MID-6',/
     >i11,' stopped in D1 MID-7',/
     >i11,' stopped in D1 MAG EXIT',/
     >i11,' stopped in D1 EXIT',/
     >i11,' stopped in BP ENTRANCE',/
     >i11,' stopped in BP EXIT',/
     > )

1100	format('!',79('-'),/,'! ',a,/,'!')
1200	format(/,'! ',a,' Coefficients',/,/,
     >(5(g18.8,','))
     >)
1300	format(/,'! ',a,' Coefficient uncertainties',/,/,
     >(5(g18.8,','))
     >)

	end

	subroutine calc_elastic_sig(theta_pol,sig_elastic)
	implicit none
	real*8 theta_pol   ! 
	real*8 sig_elastic  ! fm2/sr
	integer nang
	parameter (nang=200)
	real*8 sigma(nang),th_file(nang),frac
	integer nfile
	character*132 str_line
	logical found
	logical first /.true./
	integer nang_test
	real*8 q,q_eff,sig_mott,ratio,dsigde,dsigdth
	save
c       
	if ( first) then
	   write(*,*) ' opening file'
	   nfile=0
	   open(unit=23,status='old',file='carbon_elastic_xn.out')
	   str_line='#'
	   do while (str_line(1:1) .EQ. '#')
	      read ( 23,'(a132)') str_line
	      write(*,*) str_line
	      enddo	   
	    do while (str_line(1:4).ne.' -100')
	      nfile=nfile+1
	      if ( nfile .gt. nang) STOP
	      read(str_line
     >,'(1x,f6.3,2(1x,f5.3),2(2x,e11.5),1x,f8.3,1x,f9.3,1x,e13.5)') 
     >  th_file(nfile),q,q_eff,sig_mott,ratio,dsigde,dsigdth,sigma(nfile)
	      write(*,'(1x,f6.3,2(1x,f5.3),1x,2(1x,e11.5),1x,2(f9.3,1x),1x,e13.5)') 
     >  th_file(nfile),q,q_eff,sig_mott,ratio,dsigde,dsigdth,sigma(nfile)
	      read ( 23,'(a132)',end=100) str_line
	      enddo
	      close(unit=23)
	endif
c
 100	nang_test=1
	first=.false.
	found = .false.
	sig_elastic=-100.
	do while (nang_test .lt. nfile .and. .not. found)
	      frac= (theta_pol-th_file(nang_test))/(th_file(nang_test+1)-th_file(nang_test))
	   if (abs(frac) .lt. 1) then
	      sig_elastic=sigma(nang_test)+(sigma(nang_test+1)-sigma(nang_test))*frac
	      found = .true.
	      endif
	   nang_test = nang_test+ 1
	enddo
c
	return
	end
c
c
c
	subroutine calc_inelastic_sig(thr,e_nu,qsq,w,sigma_tot)
	implicit none
         real*8 Z,A,qsq,wsq,thr,e_nu,w
	 real*8 W1,W2,r09,Mp
	 real*8 F1,F2,siginel,sigma_qe,sigma_tot
c
	 wsq=w*w
	 siginel  = 0.
	 sigma_qe = 0.
	 Z = 6.
	 A = 12.
	 Mp =938.27
	 if ( w .gt. 1075.) then
         call F1F2IN09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2,R09)  !! Peter and Vahe new code 09	
               W1 = F1 / Mp
               W2 = F2 / e_nu
               siginel = ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
	       endif
          call F1F2QE09(Z, A, qsq/1.d6, wsq/1.d6, F1, F2)
               W1 = F1 / Mp
               W2 = F2 / e_nu
               sigma_qe = ( W2 + 2.0d0*W1*tan(thr/2.0d0)**2 )
c
        sigma_tot = siginel + sigma_qe
	return
	end
