*---------------------------------------------------------------------
      program interactive_data_set_create
*---------------------------------------------------------------------
*     JAT Jan 2005 [Version zzzz] Incorporates JLM potential as option
*     and can read BAB HF matter densities or other external file.
*     MSU/Betty Tsang July 2004 visit - for above. 
*---------------------------------------------------------------------
*     Revised March 2005 to include n and p densities in determining
*     alpha (isovector) rather than (N-Z)/A.
*     Modified: FSU (Kemper/Roeder) to (d,n) and (n,d) ('05)
*     Modified: MSU (Tsang) to include (3He,d) (June '05)
*     front5z.f is the front end - as per normal - but, if Vso=/0,
*     then this version requests input also for the bound state 
*     potential's spin-orbit geometry (at the end). You thus
*     can produce data sets with a specified rso and and aso
*     that are different to the central r0 and a0 - as previous. 
*
*     Numbering of reactions has reverted to earlier conventions
*     for (p,d) and (d,p) to allow earlier data set usage - MSU 2005
*---------------------------------------------------------------------
*     Version 6z - allows negative ireac values, for rotation
*     of amplitudes and then m-dependent cross sections - WNC 2006
*---------------------------------------------------------------------
*     Version 7z - uses the Bauge JLM parameterisation - BT 2006
*---------------------------------------------------------------------
*     Version 8z - has option to read the Sao Paulo potential for 3He
*     Helio Dias and Betty Tsang - got postponed until early 2008
*     Requires no changes to twofnr itself so can use twofnr7.f
*     Also a number of additional input error traps - courtesy of
*     United Flight 929 - 15 March 2008
*---------------------------------------------------------------------
*     Version 9z: Having found an old code for the Watanabe and finite
*     range adiabatic potentials from the Reid soft core potential
*     and deuteron - these options have been introduced -  March 2008.
*---------------------------------------------------------------------
*     Version 10: Includes the triton global potential of Li et al. 
*     For Jeff Thomas (Surrey) April/May 2009 and also corrections to
*     the Bauge JLM parameterisation (due to exchange with Pang and
*     he with Bauge, January 2011 - e-mail records these changes)
*---------------------------------------------------------------------
*     Version 11: James Benstead (Surrey/AWE)/JAT. Version includes 
*     (p,t) and (t,p) reactions as di-neutron transfer - May 2011.
*---------------------------------------------------------------------
*     Version 12: Global A=3 potential GDP08 -- Pang et al. Jan 2012
*---------------------------------------------------------------------
*     Version 13: Has (3He,alpha) option added. April 2012 
*     Includes a global apha-potential of Atzrott/Kumar (RIPL).
*---------------------------------------------------------------------
*     Version 14 just maintains counter with twofnr - no changes
*---------------------------------------------------------------------
*     Version 15: Has option to output the distorted waves - this uses
*     twofnr version 15 - and also outputs the partial wave S-matrices
*     This front has choice of deuteron wave functions for FR folding
*     and Johnson-Tandy adiabatic deuteron distorting potentials and 
*     includes the smaller additional terms in the folding expressions
*     Potentials and wave functions output to folded and deutwf: 2014
*---------------------------------------------------------------------
*     Koning-Delaroche global potential added as option for nucleon
*     potentials and adiabatic (from KD02 code, 2014)
*     Two different imaginary part shapes so potential is written. An
*     imaginary spin-orbit potential is included in the adiabatic
*     potential calculation for when needed - this very small in KD02.
*
*     Liang et al. 3He Global potential added. Two different imaginary 
*     part shapes so real and imaginary potentials are written (Pang)
*---------------------------------------------------------------------
*     Version 16: Changes for better structured reading wave functions 
*     Version 17: Option to read deuteron partial waves with rank-2
*     tensor structure included - June 2015
*---------------------------------------------------------------------
*     Version 18: The D0 and LEA ranges for different deuteron wave 
*     functions (RSC and AV18 wave function) choices are improved and
*     some tidying up is done, e.g. if negative separation energy is
*     encountered and the writing of files with jm-substatates sigma
*---------------------------------------------------------------------
*     Version 19: (3He,n) and (n,3He) options added - 9/2015
*     and also (3He,p) and (p,3He) - for JJV-D - January 2016
*     and also (t,alpha) - for JC - February 2017
*---------------------------------------------------------------------
*     Version 20: allows a radial sensitivity analysis for all cases.
*     Sensitivity test choice is signalled by a value of ktout(10)=8.
*     The bound state formfactor is set to zero beyond a given r_max
*     (if r_max positive) or for r < r_max (if r_max is negative).
*     r_max read into twofnr from fort.82 - March 2017 at TIT. 
*---------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      character fname*12,titf*46,title*58,form*12,mix*3
      character*30 ampnam(25)
      character*7 crea(15)
      integer ktout(10),ins
      real*8 a(8),b(5)
*---------------------------------------------------------------------
*     common blocks with the channel distorting potential parameters
*     uses the jlm array poti for imaginary part printout when KD02
*---------------------------------------------------------------------
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/che3/pspr(900),pspi(900),psps(900),ihesp
      common/cjlm/potr(900),poti(900),ijlm,ijlmp,ikd
      common/deut/iadia,iwat,ideutwf,ii,imso
      common/djlm/rlr,rli
      common/file/fname
      data crea,mix/'(p,d)  ','(d,p)  ','(n,d)  ','(d,n)  ','(d,t)  ',
     #        '(d,3He)','(3He,d)','(p,t)  ','(t,p)  ','(n,3He)',
     #        '(3He,n)','(p,3He)','(3He,p)','(3He,a)','(t,a)  ','mix'/
*---------------------------------------------------------------------
      print*,'====================================================== '
      print*,' Front end for generating TWOFNR transfer data sets    '
      print*,' front version 20: J.A. Tostevin March 2017         '
      print*,'====================================================== '
      print*,' Data set identifier (xxx in tran.xxx: max 12 chars)   '
      read '(a)',fname
      open(17,file='tran.'//fname,status='unknown')
      open(18,file='in.'//fname,status='unknown')
      write(18,'(a)') fname
      print '(a,a)','  >>>> ',fname
      print '(a,a)','  Output is to file:   ','tran.'//fname
      print '(a,a)','  Input saved to file: ','in.'//fname
      print*,'------------------------------------------------------ '
*---------------------------------------------------------------------
      nil=0
      unit=1.d0
*---------------------------------------------------------------------
*     version date
      iday=15
      imon=03
      iyear=2017
*---------------------------------------------------------------------
      iadia=0
      iwat=0
      ihesp=0
      form='(5e14.7)'
*---------------------------------------------------------------------
*     the following are used for the ranges with jlm/bauge potentials
*     also switch ikd for when the KD global potential used for n, p
*     ids is used for spin of n+p cluster transfer
*---------------------------------------------------------------------
      ikd=0
      ijlm=0
      rlr=0.d0
      rli=0.d0
      ids=10
      do i=1,10
       ktout(i)=0
      enddo
*---------------------------------------------------------------------
*     title line
*---------------------------------------------------------------------
      print*,' Enter title information (for info only, <= 46 chars) '
      read '(a)',titf
      write(18,'(a)') titf
      print '(a,a)','  >>>> ',titf
      title=titf//fname
*---------------------------------------------------------------------
*     menu of the available (automated input) transfer reactions
*---------------------------------------------------------------------
  932 print*,'------------------------------------------------------ '
      print*,' Reaction type: [ 1]  (p,d)       '
      print*,'                [ 2]  (d,p)       '
      print*,'                [ 3]  (n,d)       '
      print*,'                [ 4]  (d,n)       '
      print*,'                [ 5]  (d,t)       '
      print*,'                [ 6]  (d,3He)     '
      print*,'                [ 7]  (3He,d)     '
      print*,'                [ 8]  (p,t)       '
      print*,'                [ 9]  (t,p)       '
      print*,'                [10]  (n,3He)     '
      print*,'                [11]  (3He,n)     '
      print*,'                [12]  (p,3He)     '
      print*,'                [13]  (3He,p)     '
      print*,'                [14]  (3He,a)     '
      print*,'                [15]  (t,a)       '
      print*,'------------------------------------------------------ '
      print*
      print*,' Note: to store reaction transition amplitude, choose the'
      print*,' negative of the above (e.g. -4 for (d,n)) and then input' 
      print*,' Euler angles 0,0,0 when asked. For the cross sections   '
      print*,' in a rotated coordinate system use the negative reaction'
      print*,' type and input the appropriate three Euler angles       '
      print*,' ------ '
      print*,' Also, 5-->105 or -5-->-105 etc. for radial sensitivites '
*---------------------------------------------------------------------
*     read and check for a valid reaction option
      read*,ireaco
      ireac=abs(ireaco)
      irsign=ireaco/ireac
      if(ireac.gt.100) ireac=ireac-100
      if(ireac.lt.1.or.ireac.gt.15) goto 932
*     ireac now has the actual reaction choice (1-15) going forward
*---------------------------------------------------------------------
      write(18,*) ireaco
      print*,' >>>> ',ireaco
      print*
      print*,'----------------------------------------------------  '
      print'(a,a,a)','  ',crea(ireac),' reaction has been selected  '
      print*,'----------------------------------------------------  '
      print*
*---------------------------------------------------------------------
*     radial sensitivity analysis option - if abs(ireaco) > 100
      if(abs(ireaco).gt.100) then
       print*,'------------------------------------------------------'
       print*,' ktout(10)=8 if radial sensitivity analysis selected  ' 
       print*,' (twofnr reads maximum bound state radius on fort.82) '
       print*,' ----'
       print*,' cm angle for differential cross section sensitivity? '
       print*,' (twofnr will output sigma for this angle on fort.83) '
       print*,'------------------------------------------------------'
       read*,theta
       write(18,*) real(theta)
       print*,' >>>> ',real(theta),' degrees'
       ktout(10)=8
      endif
*---------------------------------------------------------------------
*     change for reading back scattering amplitudes after using
*     the rotation option. Here ktout(1)=8 to write the transition
*     amplitude (from twofnr) to an external file amp.xxx
*---------------------------------------------------------------------
      if(irsign.lt.0) then
       print*,'----------------------------------------------------  '
       print*,' You have requested to store the transition amplitude '
       print*,' or calculate the tansfer cross sections in a rotated '
       print*,' coordinate system. Now must specify the Euler angles '
       print*,' (alfa,beta,gama) in units of pi. E.g. (-1/2,-1/2,0)  '
       print*,' if z-axis is normal to the k_in x k_out plane. Input '
       print*,' 0, 0, 0 if no rotation is required                   '
       print*,'----------------------------------------------------  '
       read*,angal,angbe,angga
       print*,' Euler angles: '
       print*,' >>>> ',real(angal),real(angbe),real(angga)
       write(18,*) real(angal),real(angbe),real(angga)
       open(19,file='rot.'//fname,status='unknown')
       write(19,*) angal,angbe,angga
       print '(a,a)','  Euler angles written to file: ','rot.'//fname
       print '(a,a)','  amplitude will be written to: ','amp.'//fname
       close(19)
       ktout(1)=8
       print* 
      endif
*---------------------------------------------------------------------
*     changes for reading/printing the partial wave radial functions
*---------------------------------------------------------------------
  832 idut=0
      if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then
       print*,' Entrance channel (proton) distorted wave:  '
      else if(ireac.eq.3.or.ireac.eq.10) then
       print*,' Entrance channel (neutron) distorted wave: '
      else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +        .or.ireac.eq.14) then
       print*,' Entrance channel (3He) distorted wave:     '
      else if(ireac.eq.9.or.ireac.eq.15) then
       print*,' Entrance channel (triton) distorted wave:  '
      else
       print*,' Entrance channel (deuteron) distorted wave:'
       idut=1
      endif
      print*,' Options: '
      print*,'---------------------------------------------------- '
      print*,'    [0] Calculate - no print                 '
      print*,'    [1] Read      - no print                 '
      print*,'    [2] Read      -  print                   '
      print*,'    [3] Calculate -  print                   '
      if(idut.eq.1) then
       print*,'    [4] Read (with tensor) - no print       '
       print*,'    [5] Read (with tensor) - print          '
      endif
      print*,'---------------------------------------------------- '
      read*,ktout(3)
      if(ktout(3).lt.0.or.ktout(3).gt.5) goto 832
      write(18,*) ktout(3)
      print*,' >>>> ',ktout(3)
      if(ktout(3).eq.1.or.ktout(3).eq.2.or.ktout(3).gt.3) then
       print*,' wave functions will be read from fort.16   '
      endif
      print*
  833 idut=0
      if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
       idut=1
       print*,' Exit channel (deuteron) distorted wave:    '
      else if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
       print*,' Exit channel (proton) distorted wave:      '
      else if(ireac.eq.4.or.ireac.eq.11) then
       print*,' Exit channel (neutron) distorted wave:     '
      else if(ireac.eq.5.or.ireac.eq.8) then
       print*,' Exit channel (triton) distorted wave:      '
      else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then
       print*,' Exit channel (3He) distorted wave:         '
      else if(ireac.eq.14.or.ireac.eq.15) then
       print*,' Exit channel (alpha) distorted wave:       '
      endif
      print*,' Options: '
      print*,'---------------------------------------------------- '
      print*,'    [0] Calculate - no print                 '
      print*,'    [1] Read      - no print                 '
      print*,'    [2] Read      -  print                   '
      print*,'    [3] Calculate -  print                   '
      if(idut.eq.1) then
       print*,'    [4] Read (with tensor) - no print       '
       print*,'    [5] Read (with tensor) - print          '
      endif
      print*,'---------------------------------------------------- '
      read*,ktout(4)
      if(ktout(4).lt.0.or.ktout(4).gt.5) goto 833
      write(18,*) ktout(4)
      print*,' >>>> ',ktout(4)
      if(ktout(4).eq.1.or.ktout(4).eq.2.or.ktout(4).gt.3) then
       print*,' wave functions will be read from fort.17   '
      endif
      print*
      ins=1
*---------------------------------------------------------------------
*     start writing the tran.xxx data set
*---------------------------------------------------------------------
      write(17,101) (ktout(i),i=1,10),ins,iday,imon,iyear,title
  101 format(10i1,i2,i2,i2,i4,2x,a)
  102 format(0p8f10.4)
*---------------------------------------------------------------------
*     separation energies for A=2, 3 and 4 body light-ions (MeV)
*---------------------------------------------------------------------
*     deuteron --> n+p separation energy
      edeut=2.224573d0
*     triton binding energy  
      etwo =8.481821d0
*     triton --> deuteron+n separation energy
      etrit=etwo-edeut
*     3He --> deuteron+p separation energy
      eheli=7.718058d0-edeut
*     3He --> 2p+n separation energy
      etwop=7.718058d0
*     alpha --> 3he+n separation energy
      ealpa=20.57762d0
*     alpha --> t+p separation energy
      ealpb=19.81386d0
*---------------------------------------------------------------------
*     incident and outgoing reaction energies and masses, etc.
*---------------------------------------------------------------------
      m1=1
      if(ireac.gt.1) m1=2
      if(ireac.eq.3.or.ireac.eq.8.or.ireac.eq.10.or.ireac.eq.12) m1=1
      if(ireac.eq.7.or.ireac.eq.9.or.ireac.eq.11.or.ireac.gt.12) m1=3
      print*,' Laboratory incident energy per nucleon (MeV) '
      read*,energy
      write(18,*) real(energy)
      print*,' >>>> ',real(energy)
      energy=m1*energy
      print*, 'Total projectile energy =',real(energy)
  933 print*,' Target mass (a1) and charge (z1)    '
      read*,a1,z1
      if(z1.ge.a1.or.a1.le.0.d0.or.z1.lt.0.d0) go to 933
      write(18,*) real(a1),real(z1)
      print*,' >>>> ',real(a1),real(z1)
*---------------------------------------------------------------------
*     spins and charges
      s1=1.0
      zp1=1
      s2=0.5
      z2=z1
      zp2=1
*-----------------------------
      if(ireac.eq.1) then
       s1=0.5 
       a2=a1-1.0
       m2=2
       s2=1.0
*-----------------------------
      else if(ireac.eq.2) then
       a2=a1+1.0
       m2=1
*-----------------------------
      else if(ireac.eq.3) then
       zp1=0
       s1=0.5 
       a2=a1-1.0
       m2=2
       s2=1.0
       z2=z1-1.0
*-----------------------------
      else if(ireac.eq.4) then
       a2=a1+1.0
       m2=1
       zp2=0.0
       z2=z1+1.0
*-----------------------------
      else if(ireac.eq.5) then
       a2=a1-1.0
       m2=3
*-----------------------------
      else if(ireac.eq.6) then
       a2=a1-1.0
       m2=3
       z2=z1-1.0
       zp2=2.0
*-----------------------------
      else if(ireac.eq.7) then
       zp1=2
       s1=0.5 
       a2=a1+1.0
       m2=2
       s2=1.0
       z2=z1+1.0
*-----------------------------
      else if(ireac.eq.8) then
       s1=0.5
       a2=a1-2.0
       m2=3
*-----------------------------
      else if(ireac.eq.9) then
       s1=0.5
       a2=a1+2.0
       m2=1
*-----------------------------
      else if(ireac.eq.10) then
       zp1=0
       s1=0.5 
       a2=a1-2.0
       z2=z1-2
       m2=3
       zp2=2
       s2=0.5
*-----------------------------
      else if(ireac.eq.11) then
       zp1=2
       s1=0.5 
       a2=a1+2.0
       z2=z1+2
       m2=1
       zp2=0
       s2=0.5
*-----------------------------
      else if(ireac.eq.12) then
       zp1=1
       s1=0.5 
       a2=a1-2.0
       z2=z1-1
       m2=3
       zp2=2
       s2=0.5
*-----------------------------
      else if(ireac.eq.13) then
       zp1=2
       s1=0.5 
       a2=a1+2.0
       z2=z1+1
       m2=1
       zp2=1
       s2=0.5
*-----------------------------
      else if(ireac.eq.14) then
       zp1=2
       s1=0.5 
       a2=a1-1.0
       z2=z1
       m2=4
       zp2=2
       s2=0.0
*-----------------------------
      else if(ireac.eq.15) then
       zp1=1
       s1=0.5 
       a2=a1-1.0
       z2=z1-1.0
       m2=4
       zp2=2
       s2=0.0
      endif
*-----------------------------
      ia1=nint(a1)
      ia2=nint(a2)
      iz1=nint(z1)
      iz2=nint(z2)
*---------------------------------------------------------------------
*     energy and integration ranges line - card 1.0
*---------------------------------------------------------------------
      print*,'---------------------------------------------------- '
 873  print*,' Integration ranges: [1] use defaults   '
      print*,' (defaults: 0-30 fm in 0.10 fm steps)   '
      print*,'                     [2] specify values '
      print*,'---------------------------------------------------- '
      print*,' Note: by default these rmax and step apply to the   '
      print*,' entrance channel values. To specify that the values '
      print*,' apply to the exit channel, input instead -1 or -2   '
      read*,iiran
      irana=iiran
      iiran=abs(iiran)
      if(iiran.lt.1.or.iiran.gt.2) goto 873
      write(18,*) irana
      print*,' >>>> ',irana
 555  if(iiran.eq.1) then
       rmax=30.d0
       step1=0.10d0
       nrmax=nint(rmax/step1)
      else
       print*,' max integration radius and step '
       read*,rmax,step1
       nrmax=nint(rmax/step1)
       if(nrmax.le.900) then
        write(18,*) real(rmax),real(step1)
        print*,' >>>> ',real(rmax),real(step1)
       endif
      endif
*---------------------------------------------------------------------
      if(nrmax.gt.900) then
       print*,' Inputs require',nrmax,' radial steps but   '
       print*,' a maximum of 900 radial steps are allowed  '
       print*,' Revise your inputs accordingly:            '
       print*         
       goto 555
      endif  
      step2=step1*(a1/a2)
      if(irana.lt.0) then
       step2=step1
       step1=step2*(a2/a1)
       rmax=nrmax*step1
      endif
      print*,' integrations from 0 to',real(rmax),' fm'
      print77,'  in',nrmax,' radial steps of',real(step1),' fm'
   77 format(a,i5,a,f8.5,a)
      print*,' step length in outgoing channel =',real(step2)
      print*  
      a(1)=unit
      a(2)=nil
      a(3)=rmax
      a(4)=nil
      a(5)=nrmax
      a(6)=energy
      a(7)=nil
      a(8)=nil
      write(17,102) (a(i),i=1,7) 
      nr3max=nrmax+3
*---------------------------------------------------------------------
 874  print*,' number of partial waves [1] default (=40)   '
      print*,'                         [2] specify (<90)   '
      read*,ipw
      if(ipw.lt.1.or.ipw.gt.2) goto 874
      write(18,*) ipw
      print*,' >>>> ',ipw
      npw=40
      if(ipw.eq.2) then
 5503  print*,' input number of partial waves (<90) '
       read*,npw
       if(npw.gt.90) goto 5503
       write(18,*) npw
       print*,' >>>> ',npw
      endif
*---------------------------------------------------------------------
  408 print*,' Input the required centre of mass angles info:    '
      print*,' number of angles: step (degrees): starting value  '
      print*,' (maximum number of angles is 181)                 '
      print*,' (entering 0 0 0 will use 181  1.0  0.0 )          '
      read*,b(2),b(3),b(4)
      if(b(2).gt.181) then
       print*,' maximum number of angles is 181: please reenter  '
       go to 408
      endif
      write(18,*) real(b(2)),real(b(3)),real(b(4))
      if (abs(b(2)).lt.0.01d0) then
       b(2)=181
       b(3)=1.0
       b(4)=0.0
      endif
      b(5)=0.0
      print*,' >>>> ',real(b(2)),real(b(3)),real(b(4))
      if(ktout(10).eq.8) then
       theta=(theta-b(4))/b(3)+1
       b(5)=theta
      endif
*---------------------------------------------------------------------
*     transferred angular momenta line - card 2.2
*---------------------------------------------------------------------
      a(1)=2.2
      print*,'---------------------------------------------------- '
 911  if(ireac.eq.8.or.ireac.eq.9.or.ireac.eq.10.or.ireac.eq.11) then
       a(2)=0.0
       print*,' Enter quantum numbers L and J of transferred cluster'
       print*,' Uses simple di-nucleon model, so S = 0: enter L = J '
      else if(ireac.eq.12.or.ireac.eq.13) then
       print*,' Enter spin S of transferred n+p cluster (= 0 or 1)  '
       read*,a(2)
       ids=nint(a(2))
       if(ids.ne.0.and.ids.ne.1) then
        print*,' need S=0 or S=1 : reenter '
        go to 911
       endif
       write(18,*) real(a(2))
       print*,' >>>> ',real(a(2))
       print*,' Enter quantum numbers L and J of transferred cluster'
      else
       a(2)=0.5
       print*,' sp quantum numbers L and J of transferred nucleon   '
      endif
      read*,ltr,rjtr
      if(abs(rjtr-ltr).gt.a(2)) then
       print*,' need |J-S| < L < J+S : reenter '
       go to 911
      endif
      write(18,*) ltr,real(rjtr)
      print*,' >>>> ',ltr,real(rjtr)
 437  if(ireac.eq.8.or.ireac.eq.9.or.ireac.eq.10.or.ireac.eq.11
     *   .or.ireac.eq.12.or.ireac.eq.13)then
       print*,' number of nodes in cluster radial wave function     '
      else
       print*,' number of nodes in nucleon sp radial wave function  '
       print*,'---------------------------------------------------- '
       print*,' convention used here: the lowest state has nodes=0  '
       print*,' |+ 2|- 8|+  20|-  40|+     70|-    112|+       168| '
       print*,' | 0s| 0p|1s,0d|1p,0f|2s,1d,0g|2p,1f,0h|3s,2d,1g,0i| ' 
       print*,'---------------------------------------------------- '
      endif
      read*,nodes
      if(nodes.lt.0.or.nodes.gt.7) go to 437
      write(18,*) nodes
      print*,' >>>> ',nodes
      a(3)=ltr
      a(4)=rjtr
      write(17,102) (a(i),i=1,4) 
      print*
*---------------------------------------------------------------------
*     sort out separation energy and/or via the Q-value
*---------------------------------------------------------------------
  935 if(ireac.eq.1.or.ireac.eq.2.or.ireac.eq.5.or.ireac.eq.14) then
       print*,' specify : [1] neutron separation energy (>0 MeV)    '
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' specify : [1] two-neutron separation energy (>0 MeV)'
      else if(ireac.eq.10.or.ireac.eq.11) then
       print*,' specify : [1] two-proton separation energy (>0 MeV) '
      else if(ireac.eq.12.or.ireac.eq.13) then
       print*,' specify : [1] (n+p) separation energy (>0 MeV)      '
      else
       print*,' specify : [1] proton separation energy (>0 MeV)     '
      endif
      print*, ' or        [2] reaction Q-value (MeV)                '
      read*,ietyp
      if(ietyp.lt.1.or.ietyp.gt.2) go to 935
      write(18,*) ietyp
      print*,' >>>> ',ietyp
      if(ietyp.eq.1) then
       print*,' transferred particle separation energy (MeV: >0)    '
       read*,sn
       if(sn.le.0.0) then
        print*,' separation energy must be positive - reenter '
        goto 935
       endif
       write(18,*) real(sn)
       print*,' >>>> ',real(sn)
*---------------------------------------------------------------------
*     compute Q-value given the separation energy
*---------------------------------------------------------------------
       if(ireac.eq.1.or.ireac.eq.3) then
        qval=-sn+edeut
       else if(ireac.eq.2.or.ireac.eq.4) then
        qval=sn-edeut
       else if(ireac.eq.5) then
        qval=-sn+etrit
       else if(ireac.eq.6) then
        qval=-sn+eheli
       else if(ireac.eq.7) then
        qval=sn-eheli
       else if(ireac.eq.8) then
        qval=-sn+etwo
       else if(ireac.eq.9) then
        qval=sn-etwo
       else if(ireac.eq.10) then
        qval=-sn+etwop
       else if(ireac.eq.11) then
        qval=sn-etwop
       else if(ireac.eq.12) then
        if(ids.eq.0) qval=-sn+etwop
        if(ids.eq.1) qval=-sn+eheli
       else if(ireac.eq.13) then
        if(ids.eq.0) qval=sn-etwop
        if(ids.eq.1) qval=sn-eheli
       else if(ireac.eq.14) then
        qval=-sn+ealpa
       else if(ireac.eq.15) then
        qval=-sn+ealpb
       endif
       print*,'  Q-value is',real(qval),' MeV '
      else if(ietyp.eq.2) then
       print*,' reaction Q-value (MeV)        '
       read*,qval
       write(18,*) real(qval)
       print*,' >>>> ',real(qval)
*---------------------------------------------------------------------
*     compute separation energy given the Q-value
*---------------------------------------------------------------------
       if(ireac.eq.1.or.ireac.eq.3) then 
        sn=edeut-qval
       else if(ireac.eq.2.or.ireac.eq.4) then
        sn=edeut+qval
       else if(ireac.eq. 5) then
        sn=etrit-qval
       else if(ireac.eq. 6) then
        sn=eheli-qval
       else if(ireac.eq. 7) then
        sn=eheli+qval
       else if(ireac.eq. 8) then
        sn=etwo-qval
       else if(ireac.eq. 9) then
        sn=etwo+qval
       else if(ireac.eq.10) then
        sn=etwop-qval
       else if(ireac.eq.11) then
        sn=etwop+qval
       else if(ireac.eq.12) then
        if(ids.eq.0) sn=etwop-qval
        if(ids.eq.1) sn=eheli-qval
       else if(ireac.eq.13) then
        if(ids.eq.0) sn=etwop+qval
        if(ids.eq.1) sn=eheli+qval
       else if(ireac.eq.14) then
        sn=ealpa-qval
       else if(ireac.eq.15) then
        sn=ealpb-qval
       endif       
  476  print*,'  Separation energy is',real(sn),' MeV'
       if(sn.le.0.d0) then
        print*,'---------------------------------------------------- '
        print*,' with the input Q-value the state is particle unbound'
        print*,'    [1] reenter Q-value or separation energy?        '
        print*,'    [2] proceed and choose separation energy?        '
        print*,'    [3] exit                                         '
        print*,'---------------------------------------------------- '
        read*,ifix
        if(ifix.lt.1.or.ifix.gt.3) go to 476
        write(18,*) ifix
        print*,' >>>> ',ifix
        if(ifix.eq.1) goto 935
        if(ifix.eq.3) stop
        if(ifix.eq.2) then
  474    print*,' input chosen separation energy (>0):         '
         print*,' (note the Q-value remains the enetered value '
         print*,'  i.e. Q-value is',real(qval),' MeV)          '
         read*,sn
         if(sn.le.0.0) then
          print*,' separation energy must be positive - reenter '
          goto 474
         endif
         write(18,*) real(sn)
         print*,' >>>> ',real(sn)
        endif
       endif 
      endif
*---------------------------------------------------------------------
*     can now compute lab energy for final state potential (energy2)
*     print wavenumbers and look at likely ell mismatch of reaction
*---------------------------------------------------------------------
      print*,'====================================================== '
      ecm1=energy*a1/(a1+m1)
      print*,' entrance channel cm energy ',real(ecm1)
      fmu1=a1*m1/(a1+m1)
      fkay1=0.2195376d0*sqrt(fmu1*ecm1) 
      ecm2=ecm1+qval
      print*,' exit     channel cm energy ',real(ecm2)
      if(ecm2.lt.0.d0) then
       print*,' reaction is below threshold - stopping'
       stop
      else if(ecm2.lt.2.d0) then
       print*,' reaction is near threshold '
      endif
      fmu2=a2*m2/(a2+m2)
      fkay2=0.2195376d0*sqrt(fmu2*ecm2) 
      rad=1.2d0*(a1**0.3333333333d0)
      rl1=fkay1*rad
      rl2=fkay2*rad
      print*,' asymptotic wavenumbers and grazing angular momenta '
      print*,' kin  = ',real(fkay1),'  L(in ) = ',real(rl1)
      print*,' kout = ',real(fkay2),'  L(out) = ',real(rl2)
      rlmis=abs(rl1-rl2)
      print*,' so L mismatch is of order   ',real(rlmis),' hbar'
      print*,' from an estimated radius of ',real(rad),' fm'
      print*,' and a chosen L transfer of  ',ltr,' hbar'
      energy2=ecm2*(a2+m2)/a2
*---------------------------------------------------------------------
*     entrance channel partial waves/non-locality line - card 3.1
*---------------------------------------------------------------------
      print*,'====================================================== '
      if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then
       print*,' incident (proton) channel information           '
      else if(ireac.eq.3.or.ireac.eq.10) then
       print*,' incident (neutron) channel information          '
      else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +        .or.ireac.eq.14) then
       print*,' incident (3He) channel information              '
      else if(ireac.eq.9.or.ireac.eq.15) then
       print*,' incident (triton) channel information           '
      else
       print*,' incident (deuteron) channel information         '
      endif
      if(ktout(3).eq.1.or.ktout(3).eq.2.or.ktout(3).gt.3) then
        print*,'---------------------------------------------------- '
        print*,' wave functions will be read: chosen potential will  '
        print*,' be used when LEA finite-range correction is chosen  '
        print*,' and with any specified non-locality correction      '
        print*,'---------------------------------------------------- '
      endif
 876  print*,' nonlocality in incident channel [1] no           '
      print*,'                                 [2] yes          '      
      if(ireac.eq.2.or.ireac.eq.4) then
       print*,' It is recommended you do NOT include a non-locality  '
       print*,' with an adiabatic description of the deuteron channel'
      endif
      read*,inonloc
      if(inonloc.lt.1.or.inonloc.gt.2) go to 876
      write(18,*) inonloc
      print*,' >>>> ',inonloc
      a(5)=0.d0 
      if(inonloc.eq.2) then
       if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then
        print*,' input proton nonlocality range (~0.85 fm)      '
       else if(ireac.eq.3.or.ireac.eq.10) then
        print*,' input neutron nonlocality range (~0.85 fm)     '
       else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +         .or.ireac.eq.14) then
        print*,' input 3He nonlocality range (~0.20 fm)         '
       else if(ireac.eq.9.or.ireac.eq.15) then
        print*,' input triton nonlocality range (~0.20 fm)      '
       else
        print*,' input deuteron nonlocality range (~0.54 fm)    '
       endif
       read*,a(5)
       write(18,*) real(a(5))
       print*,' >>>> ',real(a(5))
      endif
      a(1)=3.1
      a(2)=nil
      a(3)=npw
      a(4)=nil
      write(17,102) (a(i),i=1,5) 
*---------------------------------------------------------------------
*     entrance channel masses/charges line - card 4.1
*---------------------------------------------------------------------
      a(1)=4.1
      a(2)=m1
      a(3)=a1
      a(4)=zp1
      a(5)=z1
      a(6)=s1
      print*,' target spin in incident channel '
      read*,spin1
      write(18,*) real(spin1)
      print*,' >>>> ',real(spin1)
      a(7)=spin1
      a(8)=0.0
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
  877 if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then
       print*,' incident (proton) channel potential           '
      else if(ireac.eq.3.or.ireac.eq.10) then
       print*,' incident (neutron) channel potential          '
      else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +        .or.ireac.eq.14) then
       print*,' incident (3He) channel potential              '
      else if(ireac.eq.9.or.ireac.eq.15) then
       print*,' incident (triton) channel potential           '
      else
       print*,' incident (deuteron) channel potential         '
      endif
      print*,'            [1] from those built in             '
      print*,'            [2] specify potential parameters    '
*---------------------------------------------------------------------
      read*,iopti
      if(iopti.lt.1.or.iopti.gt.2) go to 877
      write(18,*) iopti
      print*,' >>>> ',iopti
      print*,'----------------------------------------------------'
      print*,' initial potential at Elab=',real(energy),' MeV'
      if(iopti.eq.1) then
       if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then 
        call proton(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.3.or.ireac.eq.10) then 
        call neutron(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +         .or.ireac.eq.14) then
        call helium(energy,a1,z1,step1,nr3max)
       else if(ireac.eq.9.or.ireac.eq.15) then
        call triton(energy,a1,z1,step1,nr3max)
       else
        call deuteron(energy,a1,z1,ireac,inonloc)
       endif
      else
       print*,' Central terms '
       print*,' ------------- '
       print*,' Coulomb radius parameter '
       read*,rcd      
       write(18,*) real(rcd)
       print*,' >>>> ',real(rcd)      
       print*,' Real volume  : depth(>0), radius, diffuseness '
       read*,vd,rrd,ard
       write(18,*) real(vd),real(rrd),real(ard)
       print*,' >>>> ',real(vd),real(rrd),real(ard)
       print*,' Imag volume  : depth(>0), radius, diffuseness '
       read*,wd,rid,aid  
       write(18,*) real(wd),real(rid),real(aid)
       print*,' >>>> ',real(wd),real(rid),real(aid)  
*---------------------------------------------------------------------
*      print*,' Imag surface : depth(>0), radius, diffuseness '
*      read*,wisd,risid,aisid   
*---------------------------------------------------------------------
       print*,' Imag surface : depth(>0) '
       read*,wisd    
       write(18,*) real(wisd)
       print*,' >>>> ',real(wisd)    
       print*,' Spin-orbit terms '
       print*,' ---------------- '
       if(ireac.eq.1.or.ireac.eq.8.or.ireac.eq.12) then
        print*,' proton: coefficients are of L.sigma (~6.0 MeV)   '
       else if(ireac.eq.3.or.ireac.eq.10) then
        print*,' neutron: coefficients are of L.sigma (~6.0 MeV)  '
       else if(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +         .or.ireac.eq.14) then 
        print*,' 3He: coefficients are of L.sigma (~6.0 MeV)      '
       else if(ireac.eq.9.or.ireac.eq.15) then
        print*,' triton: coefficients are of L.sigma (~6.0 MeV)   '
       else 
        print*,' Careful of convention here: non-standard strength'
        print*,' deuteron: half coefficient of L.S (~3.0 MeV)     '
       endif
       print*,' Real s/orbit : depth(>0), radius, diffuseness     '
       read*,vsod,rsord,asord      
       write(18,*) real(vsod),real(rsord),real(asord)
       print*,' >>>> ',real(vsod),real(rsord),real(asord)      
       print*,' Imag s/orbit : depth(>0), radius, diffuseness     '
       read*,wsod,rsoid,asoid                                   
       write(18,*) real(wsod),real(rsoid),real(asoid)
       print*,' >>>> ',real(wsod),real(rsoid),real(asoid)               
      endif                               
*---------------------------------------------------------------------
*     potential line 1 - card 5.1
*---------------------------------------------------------------------
      a(1)=5.1
      a(2)=vd
      a(3)=(wd+wisd)
      a(4)=vsod
      a(5)=wsod
      a(6)=rrd
      a(7)=ard
      a(8)=rcd
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     potential line 2 - card 6.1
*---------------------------------------------------------------------
      a(1)=6.1
      a(2)=rsord
      a(3)=asord
      a(4)=rsoid
      a(5)=asoid
      write(17,102) (a(i),i=1,5)
*---------------------------------------------------------------------
*     potential line 3 - card 7.1
*---------------------------------------------------------------------
      a(1)=7.1
      if(abs(wd+wisd).gt.1.d-10) then
       a(2)=wisd/(wd+wisd)
      else
       a(2)=1.d0
      endif
      a(3)=rid
      a(4)=aid
      a(5)=nil
      a(6)=nil
      write(17,102) (a(i),i=1,6)
*---------------------------------------------------------------------
*     a(1)=8.1
*     a(2)=2.0
*     a(3)=nil
*     a(4)=nil
*     a(5)=nil
*     a(6)=wisd
*     a(7)=risid
*     a(8)=aisid
*     write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     exit channel partial waves/nonlocality line
*---------------------------------------------------------------------
      print*,'===================================================='
      if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
       print*,' outgoing (deuteron) channel information         '
      else if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
       print*,' outgoing (proton) channel information           '
      else if(ireac.eq.4.or.ireac.eq.11) then
       print*,' outgoing (neutron) channel information          '
      else if(ireac.eq.5.or.ireac.eq.8) then
       print*,' outgoing (triton) channel information           '
      else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then
       print*,' outgoing (3He) channel information              '
      else if(ireac.eq.14.or.ireac.eq.15) then
       print*,' outgoing (alpha) channel information            '
      endif
      if(ktout(4).eq.1.or.ktout(4).eq.2.or.ktout(4).gt.3) then
        print*,'---------------------------------------------------- '
        print*,' wave functions will be read: chosen potential will  '
        print*,' be used when LEA finite-range correction is chosen  '
        print*,' and with any specified non-locality correction      '
        print*,'---------------------------------------------------- '
      endif
 878  print*,' nonlocality in outgoing channel [1] no           '
      print*,'                                 [2] yes          '      
      if(ireac.eq.1.or.ireac.eq.3) then
       print*,' It is recommended you do NOT include a non-locality  '
       print*,' with an adiabatic description of the deuteron channel'
      endif
      read*,inonloc
      if(inonloc.lt.1.or.inonloc.gt.2) go to 878
      write(18,*) inonloc
      print*,' >>>> ',inonloc
      a(5)=0.d0 
      if(inonloc.eq.2) then
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
        print*,' input deuteron nonlocality range (~0.54 fm)    '
       else if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
        print*,' input proton nonlocality range (~0.85 fm)      '
       else if(ireac.eq.4.or.ireac.eq.11) then
        print*,' input neutron nonlocality range (~0.85 fm)     '
       else if(ireac.eq.5.or.ireac.eq.8) then
        print*,' input triton nonlocality range (~0.20 fm)      '
       else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then
        print*,' input 3He nonlocality range (~0.20 fm)         '
       else if(ireac.eq.14.or.ireac.eq.15) then
        print*,' input alpha nonlocality range (~0.20 fm)       '
       endif
       read*,a(5)
       write(18,*) real(a(5))
       print*,' >>>> ',real(a(5))
      endif
      a(1)=3.2
      a(2)=nil
      a(3)=npw
      a(4)=unit
      write(17,102) (a(i),i=1,5) 
*---------------------------------------------------------------------
*     exit channel masses/charges line - card 4.2
*---------------------------------------------------------------------
      a(1)=4.2
      a(2)=m2
      a(3)=a2
      a(4)=zp2
      a(5)=z2
      a(6)=s2
 879  print*,' target spin in outgoing channel '
      read*,spin2
      ispierr=0
*---------------------------------------------------------------------
*     check consistency of target and transferred angular momenta
*---------------------------------------------------------------------
      if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.5.or.ireac.eq.6.or.
     1 ireac.eq.8.or.ireac.eq.10.or.ireac.eq.12.or.ireac.eq.14
     2 .or.ireac.eq.15) then
       big=spin2+rjtr+0.1
       sma=abs(spin2-rjtr)-0.1
       if(spin1.gt.big.or.spin1.lt.sma) then
	ispierr=1
        print*,' ===================================================='
        print*,' input angular momenta are inconsistent:       '
        print*,real(spin2),' +',real(rjtr),' =',real(spin1),'  '
        print*,' one or more of the target and/or transferred  '
        print*,' nucleon angular momenta are wrong.            '
        print*,' ===================================================='
       endif
      else
       big=spin1+rjtr+0.1
       sma=abs(spin1-rjtr)-0.1
       if(spin2.gt.big.or.spin2.lt.sma) then
	ispierr=1
        print*,' ===================================================='
        print*,' input angular momenta are inconsistent:       '
        print*,real(spin1),' +',real(rjtr),' =',real(spin2),'  '
        print*,' one or more of the target and/or transferred  '
        print*,' nucleon angular momenta are wrong.            '
        print*,' ===================================================='
       endif
      endif
      if(ispierr.gt.0) then
      print*,' problem with spins: [1] re-enter outgoing target spin'
      print*,'                     [2] abort and start again        '
      read*,ispierr
      if(ispierr.eq.1) go to 879
      if(ispierr.eq.2) stop
      endif
*---------------------------------------------------------------------
      write(18,*) real(spin2)
      print*,' >>>> ',real(spin2)
      a(7)=spin2
      a(8)=qval
      write(17,102) (a(i),i=1,8) 
*---------------------------------------------------------------------
  880 if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then
       print*,' outgoing (deuteron) channel potential         '
      else if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
       print*,' outgoing (proton) channel potential           '
      else if(ireac.eq.4.or.ireac.eq.11) then
       print*,' outgoing (neutron) channel potential          '
      else if(ireac.eq.5.or.ireac.eq.8) then
       print*,' outgoing (triton) channel potential           '
      else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then
       print*,' outgoing (3He) channel potential              '
      else if(ireac.eq.14.or.ireac.eq.15) then
       print*,' outgoing (alpha) channel potential            '
      endif
      print*,'            [1] from those built in             '
      print*,'            [2] specify potential parameters    '
*---------------------------------------------------------------------
      read*,iopti
      if(iopti.lt.1.or.iopti.gt.2) go to 880
      write(18,*) iopti
      print*,' >>>> ',iopti
      print*,'----------------------------------------------------'
      print*,' final potential at Elab=',real(energy2),' MeV'
      if(iopti.eq.1) then
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then 
        call deuteron(energy2,a2,z2,ireac,inonloc)
       else if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
        call proton(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.4.or.ireac.eq.11) then
        call neutron(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.5.or.ireac.eq.8) then
        call triton(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then
        call helium(energy2,a2,z2,step2,nr3max)
       else if(ireac.eq.14.or.ireac.eq.15) then
        call alphap(energy2,a2,z2,step2,nr3max)
       endif
      else
       print*,' Central terms '
       print*,' ------------- '
       print*,' Coulomb radius parameter '
       read*,rcd      
       write(18,*) real(rcd)
       print*,' >>>> ',real(rcd)      
       print*,' Real volume  : depth(>0), radius, diffuseness '
       read*,vd,rrd,ard
       write(18,*) real(vd),real(rrd),real(ard)
       print*,' >>>> ',real(vd),real(rrd),real(ard)
       print*,' Imag volume  : depth(>0), radius, diffuseness '
       read*,wd,rid,aid  
       write(18,*) real(wd),real(rid),real(aid)
       print*,' >>>> ',real(wd),real(rid),real(aid)  
*---------------------------------------------------------------------
*      print*,' Imag surface : depth(>0), radius, diffuseness '
*      read*,wisd,risid,aisid   
*---------------------------------------------------------------------
       print*,' Imag surface : depth(>0) '
       read*,wisd   
       write(18,*) real(wisd)
       print*,' >>>> ',real(wisd)   
       if(ireac.ne.14.and.ireac.ne.15) then
        print*,' Spin-orbit terms '
        print*,' ---------------- '
        if(ireac.eq.2.or.ireac.eq.9.or.ireac.eq.13) then
         print*,' proton: coefficients are of L.sigma (~6.0 MeV)   '
        else if(ireac.eq.4.or.ireac.eq.11) then
         print*,' neutron: coefficients are of L.sigma (~6.0 MeV)  '
        else if(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12) then 
         print*,' 3He: coefficients are of L.sigma (~6.0 MeV)      '
        else if(ireac.eq.5.or.ireac.eq.8) then
         print*,' triton: coefficients are of L.sigma (~6.0 MeV)   '
        else 
         print*,' Careful of convention here: non-standard strength'
         print*,' deuteron: half coefficient of L.S (~3.0 MeV)     '
        endif
        print*,' Real s/orbit : depth(>0), radius, diffuseness      '
        read*,vsod,rsord,asord      
        write(18,*) real(vsod),real(rsord),real(asord)
        print*,' >>>> ',real(vsod),real(rsord),real(asord)      
        print*,' Imag s/orbit : depth(>0), radius, diffuseness      '
        read*,wsod,rsoid,asoid                                   
        write(18,*) real(wsod),real(rsoid),real(asoid)
        print*,' >>>> ',real(wsod),real(rsoid),real(asoid)   
       else
        vsod=0.d0
        wsod=0.d0
        rsord=1.d0
        asord=1.d0
        rsoid=1.d0
        asoid=1.d0
       endif
      endif
*---------------------------------------------------------------------
*     potential line 1 - card 5.2
*---------------------------------------------------------------------
      a(1)=5.2
      a(2)=vd
      a(3)=(wd+wisd)
      a(4)=vsod
      a(5)=wsod
      a(6)=rrd
      a(7)=ard
      a(8)=rcd
      write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     potential line 2 - card 6.2
*---------------------------------------------------------------------
      a(1)=6.2
      a(2)=rsord
      a(3)=asord
      a(4)=rsoid
      a(5)=asoid
      write(17,102) (a(i),i=1,5)
*---------------------------------------------------------------------
*     potential line 3 - card 7.2
*---------------------------------------------------------------------
      a(1)=7.2
      if(abs(wd+wisd).gt.1.d-10) then
       a(2)=wisd/(wd+wisd)
      else
       a(2)=1.d0
      endif
      a(3)=rid
      a(4)=aid
      a(5)=nil
      a(6)=nil
      write(17,102) (a(i),i=1,6)
*---------------------------------------------------------------------
*     a(1)=8.2
*     a(2)=2.0
*     a(3)=nil
*     a(4)=nil
*     a(5)=nil
*     a(6)=wisd
*     a(7)=risid
*     a(8)=aisid
*     write(17,102) (a(i),i=1,8)
*---------------------------------------------------------------------
*     angles line - card 9
*---------------------------------------------------------------------
      b(1)=9.0
      write(17,102) (b(i),i=1,5)
      write(17,102) 
*---------------------------------------------------------------------
*     write entrance channel nucleon JLM potentials if required
*     or write imaginary part if KD02 nucleon potential is used
*     or write Sao Paulo potential in the 3He case
*     or write Li et al. triton potential in the triton case
*---------------------------------------------------------------------
*     nucleon cases
*     for jlm potential
      if(ijlmp.eq.1.and.(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.8
     #  .or.ireac.eq.10.or.ireac.eq.12))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(potr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
      endif
*     for kd potential
      if(ikd.eq.1.and.(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.8
     #  .or.ireac.eq.10.or.ireac.eq.12))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
      endif
*---------------------------------------------------------------------
*     he3 cases
      if(ihesp.eq.1.and.(ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13
     +   .or.ireac.eq.14))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
      endif
*     triton cases
      if(ihesp.eq.1.and.(ireac.eq.9.or.ireac.eq.15))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(psps(ii),ii=1,nr3max)
      endif
*---------------------------------------------------------------------
*     calculate/write adiabatic or watanabe deuteron potentials 
*---------------------------------------------------------------------
      if(iadia.gt.0.or.iwat.gt.0) then
       print*,'----------------------------------------------------  '
       if(iadia.gt.0) then
        print*,' Now construct the Johnson-Tandy adiabatic potential:'
       else if (iwat.gt.0) then
        print*,' Now construct the Watanabe (folding model)potential:'
       endif
       print*,' select the nucleon optical potentials for folding   '
       print*,' - usually these should be consistent with that used '
       print*,' for the nucleon channel in (d,p), (n,d), ... etc.   '
       if(ireac.eq.1.or.ireac.eq.3.or.ireac.eq.7) then 
        call adiab(energy2,a2,z2,step2,nr3max)
       elseif(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.5.or.ireac.eq.6)then
        call adiab(energy ,a1,z1,step1,nr3max)
       endif
      endif
*---------------------------------------------------------------------
*     write exit channel nucleon JLM potentials if required
*     or write imaginary part if KD02 nucleon potential is used
*     or write Sao Paulo potential in the 3He case
*     or write Li et al. triton potential in the triton case
*---------------------------------------------------------------------
*     nucleon cases
*     fot jlm potential
      if(ijlmp.eq.1.and.(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.9
     #   .or.ireac.eq.11.or.ireac.eq.13))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(potr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
      endif
*     for kd potential
      if(ikd.eq.1.and.(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.9
     #   .or.ireac.eq.11.or.ireac.eq.13))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(poti(ii),ii=1,nr3max)
      endif
*---------------------------------------------------------------------
*     for 3he cases
      if(ihesp.eq.1.and.(ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
      endif
*     for triton cases
      if(ihesp.eq.1.and.(ireac.eq.5.or.ireac.eq.8))then
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspr(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(pspi(ii),ii=1,nr3max)
       write(17,'(a)') form
       write(17,'(5e14.7)')(psps(ii),ii=1,nr3max)
      endif
*---------------------------------------------------------------------
*     this ends the distorting potential inputs and outputs
*---------------------------------------------------------------------
*     following for zero-range or LEA light particle vertex
*     first formfactor line, D0^2 value needed - card 10
*---------------------------------------------------------------------
      print*,'----------------------------------------------------'
      a(1)=10.0
      a(2)=1.0
      a(3)=nil
      a(5)=nil
*---------------------------------------------------------------------
*     deuteron-nucleon overlap cases D0 assignments in a(4)
*---------------------------------------------------------------------
      if(ireac.lt.5) then
       if(ireac.eq.1.or.ireac.eq.2) then
        print*,' <p|d>  vertex constant D0 '
       else 
        print*,' <n|d>  vertex constant D0 '
       endif
 811   print*,'    [1] use default value -122.50 MeV fm^3/2 '
       print*,'    [2] use Reid SC value -125.19 MeV fm^3/2 '
       print*,'    [3] use AV18    value -126.11 MeV fm^3/2 '
       read*,idwf
       if(idwf.lt.1.or.idwf.gt.3) go to 811
       write(18,*) idwf
       print*,' >>>> ',idwf       
       a(4)=-122.5d0
       if(idwf.eq.2) a(4)=-125.19d0
       if(idwf.eq.3) a(4)=-126.11d0
       a(4)=a(4)*a(4)
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3 '
 881   print*,' use this value [1] yes '
       print*,'                [2] no  '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 881
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.5) then
       print*,' <d|t>  vertex constant D0 = -160.0 MeV fm^3/2  '
       print*,'            e.g. Phys. Rev. C 20 (1979) 1631    '
       a(4)=160.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3    '
 882   print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 882
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.6.or.ireac.eq.7) then
       print*,' <d|3He> vertex constant D0 = -160.0 MeV fm^3/2 '
       print*,'            e.g. Phys. Rev. C 20 (1979) 1631    '
       a(4)=160.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3    '
 883   print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 883
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' <p|t> vertex constant D0 = -469.0 MeV fm^3/2 '
       print*,'            e.g. Phys. Rev. C 4 (1971) 196    '
       a(4)=469.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3  '
 8882  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8882
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.10.or.ireac.eq.11) then
       print*,' <n|3He> vertex constant D0 = -469.0 MeV fm^3/2 '
       print*,'         e.g. Phys. Rev. C 4 (1971) 196    '
       a(4)=469.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3  '
 8889  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8889
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.12.or.ireac.eq.13) then
       print*,' <p|3He> vertex constant D0 = -469.0 MeV fm^3/2 '
       print*,'         e.g. Phys. Rev. C 4 (1971) 196    '
       a(4)=469.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3  '
 8669  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8669
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.14) then
 5559  print*,' <3He|a> vertex constant D0                       '
       print*,'   [1] -539.0 MeV fm^3/2 (from dwuck4 range)      '
       print*,'   [2] -455.0 MeV fm^3/2 (Barnwell,ZR+nonloc)     '
       print*,'   [3] -275.0 MeV fm^3/2 (Barnwell,ZR+LEA+nonloc) '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.3) go to 5559
       write(18,*) ianq
       print*,' >>>> ',ianq
       a(4)=539.0d0**2
       if(ianq.eq.2) a(4)=455.0d0**2
       if(ianq.eq.3) a(4)=275.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3    '
 8872  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8872
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
*---------------------------------------------------------------------
      else if(ireac.eq.15) then
 5779  print*,' <t|a> vertex constant D0                         '
       print*,'   [1] -539.0 MeV fm^3/2 (from dwuck4 range)      '
       print*,'   [2] -455.0 MeV fm^3/2 (Barnwell,ZR+nonloc)     '
       print*,'   [3] -275.0 MeV fm^3/2 (Barnwell,ZR+LEA+nonloc) '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.3) go to 5779
       write(18,*) ianq
       print*,' >>>> ',ianq
       a(4)=539.0d0**2
       if(ianq.eq.2) a(4)=455.0d0**2
       if(ianq.eq.3) a(4)=275.0d0**2
       print*,' this gives D0^2 = ',real(a(4)),' MeV^2 fm^3    '
 8772  print*,' use this default [1] yes '
       print*,'                  [2] no '
       read*,ianq
       if(ianq.lt.1.or.ianq.gt.2) go to 8772
       write(18,*) ianq
       print*,' >>>> ',ianq
       if(ianq.eq.2) then
        print*,' input D0^2 MeV^2 fm^3 '
        read*,a(4)
        write(18,*) a(4)
        print*,' >>>> ',a(4)
       endif
      endif
*---------------------------------------------------------------------
      write(17,33) (a(i),i=1,5)
  33  format(0p3f10.4,f10.2,3f10.4)
*---------------------------------------------------------------------
*     finite range correction: light particle vertex - card 10.01
*---------------------------------------------------------------------
      a(1)=10.01
      print*,'----------------------------------------------------'
      if(ireac.eq.1) then
       print*,' Treatment of finite range (fnrng) of <p|d> vertex  '
      else if(ireac.eq.3) then
       print*,' Treatment of finite range (fnrng) of <n|d> vertex  '
      else if(ireac.eq.2) then
       print*,' Treatment of finite range (fnrng) of <d|p> vertex  '
      else if(ireac.eq.4) then
       print*,' Treatment of finite range (fnrng) of <d|n> vertex  '
      else if(ireac.eq.5) then
       print*,' Treatment of finite range (fnrng) of <d|t> vertex  '
      else if(ireac.eq.6) then
       print*,' Treatment of finite range (fnrng) of <d|3He> vertex'
      else if(ireac.eq.7) then
       print*,' Treatment of finite range (fnrng) of <3He|d> vertex'
      else if(ireac.eq.8) then
       print*,' Treatment of finite range (fnrng) of <p|t> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
       print*,' (p,t) or (t,p) reactions.                    '
      else if(ireac.eq.9) then
       print*,' Treatment of finite range (fnrng) of <t|p> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
       print*,' (p,t) or (t,p) reactions.                    '
      else if(ireac.eq.10) then
       print*,' Treatment of finite range (fnrng) of <n|3He> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
       print*,' (p,t) or (t,p) reactions.                    '
      else if(ireac.eq.11) then
       print*,' Treatment of finite range (fnrng) of <3He|n> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
       print*,' (p,t) or (t,p) reactions.                    '
      else if(ireac.eq.12) then
       print*,' Treatment of finite range (fnrng) of <p|3He> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
      else if(ireac.eq.13) then
       print*,' Treatment of finite range (fnrng) of <3He|p> vertex  '
       print*,' Glendenning: Ch 9, Nuclear Spectroscopy and  '
       print*,' Reactions - recommends zero range for L<= 2  '
      else if(ireac.eq.14) then
       print*,' Treatment of finite range (fnrng) of <3He|a> vertex'
      else if(ireac.eq.15) then
       print*,' Treatment of finite range (fnrng) of <t|a>   vertex'
      endif
*---------------------------------------------------------------------
 884  print*,'          [1] zero-range (fnrng = 0)               '
      print*,'          [2] local-energy (default values)        '
*---------------------------------------------------------------------
*     if deuteron-nucleon vertex - default built-in values
*---------------------------------------------------------------------
      if(ireac.lt.5) then
       if(idwf.eq.1) print*,'              default fnrng=0.745712 fm'
       if(idwf.eq.2) print*,'              Reid SC fnrng=0.745650 fm'
       if(idwf.eq.3) print*,'              AV18    fnrng=0.759392 fm'
      endif
*---------------------------------------------------------------------
      print*,'          [3] local-energy (specify fnrng value)   '
      read*,izr
      if(izr.lt.1.or.izr.gt.3) go to 884
      write(18,*) izr
      print*,' >>>> ',izr
      print*,'---------------------------------------------------- '
*---------------------------------------------------------------------
      if(izr.eq.1) then
       a(2)=0.0
       print*,' Zero-range calculation (fnrng = 0) '
*---------------------------------------------------------------------
      else if(izr.eq.2) then
       if(ireac.lt.5) then
        if(idwf.eq.1) then
         a(2)=0.745712d0
         print*,' Default finite range factor, fnrng=0.745712 fm'
        else if(idwf.eq.2) then
         a(2)=0.745650d0
         print*,' Reid SC finite range factor, fnrng=0.745650 fm'
        else if(idwf.eq.3) then
         a(2)=0.759392d0
         print*,' AV18    finite range factor, fnrng=0.759392 fm'
        endif
       else if(ireac.eq.5) then
        a(2)=0.746269d0
        print*,' Hulthen finite range factor, fnrng=0.746269 fm'
        print*,'                 Nucl. Phys. A234 (1974) 301 '
       else if(ireac.eq.6.or.ireac.eq.7) then
        a(2)=0.746269d0
        print*,' Hulthen finite range factor, fnrng=0.746269 fm'
        print*,'                 Nucl. Phys. A234 (1974) 301 '
       else if(ireac.eq.8.or.ireac.eq.9) then
        a(2)=0.746269d0
        print*,' Hulthen finite range factor, fnrng=0.746269 fm'
        print*,' Uses same as (d,t), (d,3He), (3He,d) presently'
       else if(ireac.eq.10.or.ireac.eq.11) then
        a(2)=0.746269d0
        print*,' Hulthen finite range factor, fnrng=0.746269 fm'
        print*,' Uses same as (d,t), (d,3He), (3He,d) presently'
       else if(ireac.eq.12.or.ireac.eq.13) then
        a(2)=0.746269d0
        print*,' Hulthen finite range factor, fnrng=0.746269 fm'
        print*,' Uses same as (d,t), (d,3He), (3He,d) presently'
       else if(ireac.eq.14.or.ireac.eq.15) then
        a(2)=0.7d0
        print*,' Hulthen finite range factor, fnrng=0.70 fm  '
        print*,'                 from dwuck4 appendix        '
       endif
*---------------------------------------------------------------------
      else if(izr.eq.3) then     
       print*,' input local-energy vertex interaction range (fm)  '
       print*,' for conventions:  DelVecchio and Daehnick PRC     '
       print*,' Vol 6 (1972) p2095 (DVD)                          '
       print*,' [if +ve,  Hulthen, fnrng = R of DVD             ] ' 
       print*,' [fnrng = 1/beta of Nucl. Phys. A241 (1975)  36  ] '
       print*,' [if -ve, Gaussian, fnrng = 1/(2*epsilon) of DVD ] '
       read*,a(2)
       write(18,*) real(a(2))
       print*,' >>>> ',real(a(2))
      endif
      write(17,102) (a(i),i=1,2)
*---------------------------------------------------------------------
*     formfactor line
*---------------------------------------------------------------------
      a(1)=10.41
      a(2)=nodes
*     set the charges product correctly for stripping/pickup cases
      a(3)=nil
      if(ireac.eq.3.or.ireac.eq.6.or.ireac.eq.10.or.ireac.eq.12
     +   .or.ireac.eq.15) a(3)=z2
      if(ireac.eq.10) a(3)=2.d0*a(3)
      if(ireac.eq.4.or.ireac.eq.7.or.ireac.eq.11.or.ireac.eq.13)a(3)=z1
      if(ireac.eq.11) a(3)=2.d0*a(3)
      a(4)=sn
      a(5)=1.0
      if(ireac.eq.8.or.ireac.eq.9.or.ireac.eq.10.or.
     #    ireac.eq.11.or.ireac.eq.12.or.ireac.eq.13) a(5)=2.0
*     set the core mass correctly for stripping/pickup cases
      if(ireac.eq.2.or.ireac.eq.4.or.ireac.eq.7.or.ireac.eq.9
     #    .or.ireac.eq.11.or.ireac.eq.13) then
       a(6)=a1
      else
       a(6)=a2
      endif
      a(7)=nil
      write(17,102) (a(i),i=1,7)
*---------------------------------------------------------------------
*     formfactor potential line
*---------------------------------------------------------------------
      a(1)=10.42
      print*,'----------------------------------------------------'
      if(ireac.eq.1.or.ireac.eq.2.or.ireac.eq.5.or.ireac.eq.14) then
       print*,' neutron binding potential     '
      else if(ireac.eq.8.or.ireac.eq.9) then
       print*,' di-neutron binding potential  '
      else if(ireac.eq.10.or.ireac.eq.11) then
       print*,' di-proton binding potential   '
      else if(ireac.eq.12.or.ireac.eq.13) then
       print*,' (n+p) binding potential       '
      else
       print*,' proton binding potential      '
      endif
      print*,' radius and diffuseness (e.g. 1.25  0.65 fm)'
      read*,a(2),a(4)
      radius =a(2)
      diffuse=a(4)
      write(18,*) real(a(2)),real(a(4))
      print*,' >>>> ',real(a(2)),real(a(4))
      a(3)=a(2)
      if(ireac.eq.8.or.ireac.eq.9.or.ireac.eq.10.or.ireac.eq.11) then
       print*,' Spin-orbit: di-nucleon so input strength Vso = 0'
      else if(ireac.eq.12.or.ireac.eq.13) then
       if(ids.eq.0) print*,' Spin-orbit: S=0 so input  Vso = 0'
       if(ids.eq.1) then
        print*,' Spin-orbit: S=1 so input  Vso                '
        print*,' Spin-orbit: Vso is strength of L.S (~6.0 MeV)'
       endif
      else
       print*,' Spin-orbit: strength of l.sigma (~6.0 MeV)    '
      endif
      read*,a(5)  
      if(ids.eq.0) a(5)=0.d0
      if(ireac.eq. 8.or.ireac.eq. 9) a(5)=0.d0
      if(ireac.eq.10.or.ireac.eq.11) a(5)=0.d0
      vso=a(5)    
      write(18,*) real(a(5))
      print*,' >>>> ',real(a(5))      
      print*,' Bound state non-locality (0 usually)       '
      read*,a(6)      
      write(18,*) real(a(6))
      print*,' >>>> ',real(a(6))                              
      write(17,102) (a(i),i=1,6)
*---------------------------------------------------------------------
*     formfactor line (for distinct spin-orbit geometry if |vso|>0)
*---------------------------------------------------------------------
      if(abs(vso).gt.1.d-3) then
       a(1)=10.43
       print*,' Bound state spin-orbit radius parameter  '
       print*,' (if 0 entered, use the same geometry as  '
       print*,'  input for the real central interaction)  '
       read*,a(2)      
       write(18,*) real(a(2))
       print*,' >>>> ',real(a(2))      
       if(abs(a(2)).lt.0.01d0) then
	a(2)=radius
        a(3)=diffuse
	goto 887
       endif
       print*,' Bound state spin-orbit diffuseness parameter  '
       read*,a(3)      
       write(18,*) real(a(3))
       print*,' >>>> ',real(a(3))      
 887   write(17,102) (a(i),i=1,3)
      endif
      write(17,102) 
*---------------------------------------------------------------------
*     ends creation of tran.xxx data set - print summary information
*---------------------------------------------------------------------
      close(17)
      print*
      print*,'====================================================    '
      print'(a,a)','  tran.'//fname,'     dataset has been created:   '
      print 702,'  for the [',ia1,',',iz1,'] ',crea(ireac),'[',ia2,
     #          ',',iz2,'] reaction'
  702 format(a,i3,a,i2,a,a,a,i3,a,i2,a)
      print*,'====================================================    '
      print*,' when running twofnr the output files are as follows:   '
      print*
      print*,' Normal  kinematics observables                          '
      print*,' ------------------------------'
      print '(a,a)','  20.'//fname,' cm - with all spin observables   '
      print '(a,a)','  21.'//fname,' cm - sigma and Ay only           '
      print '(a,a)','  22.'//fname,' lab - light in - light out       '
      print*
      print*,' Inverse kinematics observables                         '
      print*,' ------------------------------ '
      print '(a,a)','  23.'//fname,' lab - heavy in - heavy out       '
      print '(a,a)','  24.'//fname,' lab - heavy in - light out       '
      if(ktout(1).eq.8) then
       print*
       print*,' Rotated frame cm observables                          '
       print*,' ---------------------------- '
       print*,' input Euler angles (alfa,beta,gama) in units of pi    '
       print'(3x,3f8.3)',real(angal),real(angbe),real(angga)
       print '(a,a)','  39.'//fname,' m-substates sigma_jm of transfer'
       print '(a,a)','  40.'//fname,' cm - sigma of rotated amplitudes' 
       print*
       print '(a,a)','  amp.'//fname,'transfer amplitude (for mixing) ' 
       print*,'----------------------------------------------------   '
       print*,' weighted linear combinations of such stored amp.xxx   '
       print*,' amplitudes and residual nucleus substate populations  '
       print*,' can be computed with a suitable additional data set   ' 
  814  print*,' Options:' 
       print*,'----------------------------------------------------   '
       print*,'    [1] do not create at this time ' 
       print*,'    [2] create a mixing data set   '
       read*,imix
       if(imix.lt.1.or.imix.gt.2) go to 814
       write(18,*) imix
       print*,' >>>> ',imix
       if(imix.eq.1) then
        print*,'====================================================  '
        stop
       endif
       ins=8
       ktout(1)=9
       ktout(3)=0
       ktout(4)=0
       fname=mix//fname
       title=titf//fname
       open(17,file='tran.'//fname,status='unknown')
       print '(a,a)','  Mix data is to file:  ','tran.'//fname
       write(17,101) (ktout(i),i=1,10),ins,iday,imon,iyear,title
  664  print*,' How many amplitudes are to be combined? '
       read*,namp
       if(namp.le.0) then
        print*,' need at least one amplitude '
        go to 664
       endif
       write(18,*) namp
       print*,' >>>> ',namp
       do ia=1,namp
        print'(a,i3)','  name, i.e. xxx (of amp.xxx) amplitude',ia
        read'(a)',ampnam(ia)
        write(18,'(a)') ampnam(ia)
        print'(a,i3)','  spectroscopic amplitude for amplitude',ia
        read*,specamp
        write(18,*) specamp
        print'(a,a12,f8.4)','  >>>> ',ampnam(ia),real(specamp)
        write(17,432) real(ia),specamp,ampnam(ia)
  432   format(f3.1,7x,f7.4,3x,a)
       enddo
       write(17,*)'                                                 '
       close(17)
       print*,'==================================================== '
       print'(a,a)','  tran.'//fname,' dataset has been created  '
      endif
      print*,'====================================================  '
      end
*---------------------------------------------------------------------
      subroutine proton(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      character nucleon,form*12
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/cjlm/potr(900),poti(900),ijlm,ijlmp,ikd
*---------------------------------------------------------------------
*     woods-saxon potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is coefficient of l.sigma for nucleons (ww~6MeV)
*---------------------------------------------------------------------
 887  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)       '
      print*,'            Phys Rev 182 (1969) 1190             '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57               '
      print*,' [3] Menet (30<E<60) see: ADNDT 17 (1976) p6     '
      print*,' [4] Perey (E<20MeV) see: ADNDT 17 (1976) p6     '
      print*,' [5] JLM microscopic optical potential           '
      print*,'            Bauge implementation, PRC 58, 1120   '
      print*,' [6] Koning-Delaroche (KD02) global potential    '
      print*,'            Nucl Phys A713 (2003) 231            '
      print*,'----------------------------------------------------'
*---------------------------------------------------------------------
      ijlmp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.6) go to 887
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=54.d0-0.32d0*ed+0.4d0*e/a13+24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.32d0
       api=0.51d0+0.7d0*an
       wpi=11.8d0-0.25d0*ed+12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*ed-2.7d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV proton energy ')
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=6.2d0
       rsord=1.01d0
       asord=0.75d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0
       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0
       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=1.73d0*z/rrc
       erp=e-ecpp
       vrp=v0+vt*((n-z)/a)+erp*ve
       rp=r0*a13+r00
       rpn=rp/a13
       ap=a0
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       awp=aw
       wsp=(ws0+wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       rcd=rcn
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       vd=vrp
       rrd=rpn
       ard=ap
       wd=wvp
       rid=rwpn
       aid=awp
       wisd=wsp
       risid=rwpn
       aisid=awp
*---------------------------------------------------------------------
*      CH86 parameters
       vsod=5.9d0
       rsord=(1.39d0*a13-1.43)/a13
       asord=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vsod=5.9d0
       rsord=(1.34d0*a13-1.20)/a13
       asord=0.63d0
*---------------------------------------------------------------------
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=49.9d0-0.22*e+26.4*(n-z)/a+0.4*z/a13
       wvp=1.2+0.09*e
       wsp=4.2-0.05*e+15.5*(n-z)/a
       if(wsp.lt.0.d0) wsp=0.d0
       awp=0.74d0-0.008*e+(n-z)/a
       rrd=1.16d0
       ard=0.75d0
       rid=1.37d0
       print 25,a
       print 101,z,e
   25  format(1h ,' Menet potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=6.04d0
       rsord=1.064d0
       asord=0.78d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=53.3d0-0.55*e+27.0*(n-z)/a+0.4*z/a13
       wvp=0.d0
       wsp=13.5d0
       rrd=1.25d0
       ard=0.65d0
       rid=1.25d0
       awp=0.47d0
       print 24,a
       print 101,z,e
   24  format(1h ,' Perey potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.50d0
       rsord=1.25d0
       asord=0.47d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.5) then
       ijlmp=1
       ijlm=1
       print 29,a1
   29  format(1h ,' JLM potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
       nucleon='p'
       call jlm(potr,poti,nucleon,a1,z1,energy,step,nrmax)
      endif
*---------------------------------------------------------------------
      if(inopt.eq.6) then
*      KD02 Koning-Delaroche for proton (Pang/KD02 code)
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*      defined in KD02
       rvKD=1.3039-0.4054*A**(-1./3.)
       avKD=0.6778-1.487e-4*A
       rwKD=rvKD
       awKD=avKD
       v4KD=7.0e-9
       w2KD=73.55+0.0795*A
       rvdKD=1.3424-0.01585*A**(1./3.)
       rwdKD=rvdKD
       vdKD=0.
       d2KD=0.0180+3.802e-3/(1.+exp((A-156.)/8.))
       d3KD=11.5
       vso1KD=5.922+0.0030*A
       vso2KD=0.0040
       rvsoKD=1.1854-0.647*A**(-1./3.)
       rwsoKD=rvsoKD
       avsoKD=0.59
       awsoKD=avsoKD
       wso1KD=-3.1
       wso2KD=160.
       efKD=-8.4075+0.01378*A
       v1KD=59.30+21.0*real(N-Z)/A-0.024*A
       v2KD=7.067e-3+4.23e-6*A
       v3KD=1.729e-5+1.136e-8*A
       w1KD=14.667+0.009629*A
       avdKD=0.5187+5.205e-4*A
       awdKD=avdKD
       d1KD=16.0+16.0*real(N-Z)/A
       rcKD=1.198+0.697*A**(-2./3.)+12.994*A**(-5./3.)
       fKD=E-efKD
*      Coulomb
       VcKD=1.73/rcKD*Z/(A**(1./3.))
       vcoulKD=VcKD*v1KD*(v2KD-2.*v3KD*fKD+3.*v4KD*fKD*fKD)
       vKD=v1KD*(1.-v2KD*fKD+v3KD*fKD**2-v4KD*fKD**3)+vcoulKD
       wKD=w1KD*fKD**2/(fKD**2+w2KD**2)
       vdKD=0.
       wdKD=d1KD*fKD**2*exp(-d2KD*fKD)/(fKD**2+d3KD**2)
       vsoKD=vso1KD*exp(-vso2KD*fKD)
       wsoKD=wso1KD*fKD**2/(fKD**2+wso2KD**2)
       vd  = vKD
       rrd = rvKD
       ard = avKD
       wd  = wKD
       rid = rrd
       aid = ard
       wisd  = wdKD
       risid = rwdKD
       aisid = awdKD
       print 60,a
       print 101,z,e
   60  format(1h ,' KD02 potential for a = ',f5.1)
       rcd=rcKD
       print*,' Coulomb radius parameter = ',real(rcd)
       print 1112
 1112  format(1h ,'    v      rv     av     w      rw     aw ')
       print 122, vd,rrd,ard,wd,rid,aid
       print 1312
 1312  format(1h ,'    ws     rws    aws  ')
       print 122, wisd,risid,aisid
  122  format(1h ,9f7.3)
*      KD02 spin-orbit 
       vsod=vsoKD
       rsord=rvsoKD
       asord=avsoKD
       wsod=wsoKD
       rsoid=rwsoKD
       asoid=awsoKD
       print 1129
 1129  format(1h ,'  vso    rso    aso    wso    rsoi   asoi ')
       print 122,vsod,rsord,asord,wsod,rsoid,asoid
       ikd=1
*---------------------------------------------------------------------
*      twofnr assumes the same geometry for the volume and surface
*      imaginary potentials. KD does not have the same geometry so
*      compute to print the imaginary part as two shapes involved
*---------------------------------------------------------------------
       a13=a**(1./3.)
       brid =rid  *a13
       brisd=risid*a13
       do ii=1,nrmax
        r=ii*step
        poti(ii)=ws(r,wd,brid,aid)+wsd(r,wisd,brisd,aisid)
       enddo
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine neutron(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      character nucleon,form*12
      common/cjlm/potr(900),poti(900),ijlm,ijlmp,ikd
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
*---------------------------------------------------------------------
*     potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is coefficient of l.sigma for nucleons (ww~6MeV)
*---------------------------------------------------------------------
 888  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)       '
      print*,'            Phys Rev 182 (1969) 1190             '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57               '
*     print*,' [3] Menet (30<E<60) see: ADNDT 17 (1976) p6     '
*     print*,' [4] Perey (E<20MeV) see: ADNDT 17 (1976) p6     '
      print*,' [3] JLM microscopic optical potential           '
      print*,'            Bauge implementation PRC 58, 1120    '
      print*,' [4] Koning-Delaroche global potential           '
      print*,'            Nucl Phys A713 (2003) 231            '
      print*,'----------------------------------------------------'
*---------------------------------------------------------------------
      ijlmp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.4) go to 888
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=56.3d0-0.32d0*ed-24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.26d0
       api=0.58d0
       wpi=13.0d0-0.25d0*ed-12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       wvpi=0.22d0*ed-1.56d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV neutron energy ')
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=6.2d0
       rsord=1.01d0
       asord=0.75d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0
       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0
       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=0.d0
       erp=e-ecpp
       vrp=v0-vt*((n-z)/a)+erp*ve
       rp=r0*a13+r00
       rpn=rp/a13
       ap=a0
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       awp=aw
       wsp=(ws0-wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       rcd=rcn
*      print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       vd=vrp
       rrd=rpn
       ard=ap
       wd=wvp
       rid=rwpn
       aid=awp
       wisd=wsp
       risid=rwpn
       aisid=awp
*---------------------------------------------------------------------
*      CH86 parameters
       vsod=5.9d0
       rsord=(1.39d0*a13-1.43)/a13
       asord=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vsod=5.9d0
       rsord=(1.34d0*a13-1.20)/a13
       asord=0.63d0
*---------------------------------------------------------------------
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.98) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=49.9d0-0.22*e+26.4*(n-z)/a+0.4*z/a13
       wvp=1.2+0.09*e
       wsp=4.2-0.05*e+15.5*(n-z)/a
       if(wsp.lt.0.d0) wsp=0.d0
       awp=0.74d0-0.008*e+(n-z)/a
       rrd=1.16d0
       ard=0.75d0
       rid=1.37d0
       print 25,a
       print 101,z,e
   25  format(1h ,' Menet potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=6.04d0
       rsord=1.064d0
       asord=0.78d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.99) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=53.3d0-0.55*e+27.0*(n-z)/a+0.4*z/a13
       wvp=0.d0
       wsp=13.5d0
       rrd=1.25d0
       ard=0.65d0
       rid=1.25d0
       awp=0.47d0
       print 24,a
       print 101,z,e
   24  format(1h ,' Perey potential for a = ',f5.1)
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
*      print 11
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.50d0
       rsord=1.25d0
       asord=0.47d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       ijlmp=1
       ijlm=1
       print 29,a1
   29  format(1h ,' JLM potential for a = ',f5.1)
       rcd=1.25d0
*      print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
       nucleon='n'
       call jlm(potr,poti,nucleon,a1,z1,energy,step,nrmax)
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
*      KD02 Koning-Delaroche for neutron (Pang code)
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
*      defined in KD02
       rvKD=1.3039-0.4054*A**(-1./3.)
       avKD=0.6778-1.487e-4*A
       rwKD=rvKD
       awKD=avKD
       v4KD=7.0e-9
       w2KD=73.55+0.0795*A
       rvdKD=1.3424-0.01585*A**(1./3.)
       rwdKD=rvdKD
       vdKD=0.
       d2KD=0.0180+3.802e-3/(1.+exp((A-156.)/8.))
       d3KD=11.5
       vso1KD=5.922+0.0030*A
       vso2KD=0.0040
       rvsoKD=1.1854-0.647*A**(-1./3.)
       rwsoKD=rvsoKD
       avsoKD=0.59
       awsoKD=avsoKD
       wso1KD=-3.1
       wso2KD=160.
       efKD=-11.2814+0.02646*A
       v1KD=59.30-21.0*real(N-Z)/A-0.024*A
       v2KD=7.228e-3-1.48e-6*A
       v3KD=1.994e-5-2.0e-8*A
       w1KD=12.195+0.0167*A
       d1KD=16.0-16.0*real(N-Z)/A
       avdKD=0.5446-1.656e-4*A
       awdKD=avdKD
       rcKD=1.
       fKD=E-efKD
*      Coulomb
       vcoulKD=0.
       vKD=v1KD*(1.-v2KD*fKD+v3KD*fKD**2-v4KD*fKD**3)+vcoulKD
       wKD=w1KD*fKD**2/(fKD**2+w2KD**2)
       vdKD=0.
       wdKD=d1KD*fKD**2*exp(-d2KD*fKD)/(fKD**2+d3KD**2)
       vsoKD=vso1KD*exp(-vso2KD*fKD)
       wsoKD=wso1KD*fKD**2/(fKD**2+wso2KD**2)
       vd  = vKD
       rrd = rvKD
       ard = avKD
       wd  = wKD
       rid = rrd
       aid = ard
       wisd  = wdKD
       risid = rwdKD
       aisid = awdKD
       print 60,a
       print 101,z,e
   60  format(1h ,' KD02 systematics for a = ',f5.1)
       rcd=rcKD
       print*,' Coulomb radius parameter = ',real(rcd)
       print 1112
 1112  format(1h ,'    v      rv     av     w      rw     aw ')
       print 122, vd,rrd,ard,wd,rid,aid
       print 1312
 1312  format(1h ,'    ws     rws    aws  ')
       print 122, wisd,risid,aisid
  122  format(1h ,9f7.3)
*      KD02 spin-orbit 
       vsod=vsoKD
       rsord=rvsoKD
       asord=avsoKD
       wsod=wsoKD
       rsoid=rwsoKD
       asoid=awsoKD
       print 1129
 1129  format(1h ,'  vso    rso    aso    wso    rsoi   asoi ')
       print 122,vsod,rsord,asord,wsod,rsoid,asoid
       ikd=1
*---------------------------------------------------------------------
*      twofnr assumes the same geometry for the volume and surface
*      imaginary potentials. KD does not have the same geometry so
*      compute to print the imaginary part as two shapes involved
*---------------------------------------------------------------------
       a13=a**(1./3.)
       brid =rid  *a13
       brisd=risid*a13
       do ii=1,nrmax
        r=ii*step
        poti(ii)=ws(r,wd,brid,aid)+wsd(r,wisd,brisd,aisid)
       enddo
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine deuteron(energy,a1,z1,ireac,inonloc)
      implicit real*8(a-h,o-z) 
      character fname*12
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/deut/iadia,iwat,ideutwf,ii,imso
      common/file/fname
*---------------------------------------------------------------------
 889  print*,'----------------------------------------------------'
      print*,' Optical potentials for DWBA                     '
      print*,'                                                 '
      print*,' [1] Lohr-Haeberli (A>40 8<E<13 MeV)             '
      print*,'             see: ADNDT 17 (1976) p6             '
      print*,' [2] Perey-Perey (12<E<25 MeV) no spin-orbit     '
      print*,'             see: ADNDT 17 (1976) p6             ' 
      print*,' [3] Daehnick Global (A>27 12<E<90 MeV)          '
      print*,'     Phys. Rev. C 21, 2253 (1980)                ' 
      print*,' [4] Watanabe folding model potential from       '
      print*,'     nucleon potentials                          ' 
      print*,'----------------------------------------------------'
      if(ireac.lt.5)then
       print*,' Adiabatic potentials for breakup               '
       print*,'                                                '
       print*,' [5] Zero range adiabatic potential             '
       print*,'     Johnson-Soper PRC 1 (1970) 976             '
       print*,' [6] Finite range adiabatic potential           '
       print*,'     Johnson-Tandy NPA 235 (1974) 56            '
       print*,'----------------------------------------------------'
      endif
*---------------------------------------------------------------------
      read*,inopt
      if(ireac.lt.5) then
       if(inopt.lt.1.or.inopt.gt.6) go to 889
      else
       if(inopt.lt.1.or.inopt.gt.4) go to 889
      endif
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=91.13+2.2*z/a13
       wvp=0.d0
       wsp=218.d0/a13/a13
       rrd=1.05d0
       ard=0.86d0
       rid=1.43d0
       awp=0.5d0+0.013*a13*a13
       print 20,a
   20  format(1h ,' LH deuteron potential for a = ',f5.1)
       print 101,z,e
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV deuteron energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
   12  format(1h ,7f7.3)
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=7.0d0/2.d0
       rsord=0.75d0
       asord=0.50d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       a=a1
       e=energy
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=81.1-0.22*e+2.0*z/a13
       wvp=0.d0
       wsp=14.4+0.24*e
       rrd=1.15d0
       ard=0.81d0
       rid=1.34d0
       awp=0.68d0
       print 21,a
   21  format(1h ,' P-P deuteron potential for a = ',f5.1)
       print 101,z,e
       rcd=1.15d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       vd=vrp
       wd=wvp
       aid=awp
       wisd=wsp
       risid=rid
       aisid=aid
       vsod=0.0d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       a=a1
       e=energy
       bet=-(e/100.d0)**2
       bet=exp(bet)
       z=z1
       n=nint(a-z)
       a13=a**0.3333333333d0
       vrp=88.5-0.26*e+0.88*z/a13
       wvp=(12.2+0.026*e)*(1.d0-bet)
       wsp=(12.2+0.026*e)*(bet)
       rrd=1.17d0
       ard=0.709d0+0.0017*e
       rid=1.325d0
       awp=0.53d0+0.07*a13
       awp=awp-0.04*exp(-((  8-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 20-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 28-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 50-n)/2.d0)**2)
       awp=awp-0.04*exp(-(( 82-n)/2.d0)**2)
       awp=awp-0.04*exp(-((126-n)/2.d0)**2)
       print 23,a
   23  format(1h ,' Daehnick deuteron potential for a = ',f5.1)
       print 101,z,e
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=Vrp
       wd=Wvp
       aid=awp
       wisd=Wsp
       risid=rid
       aisid=aid
       vsod=(7.33d0-0.029*e)/2.d0
       rsord=1.07d0
       asord=0.66d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 111
       print 12,vrp,rrd,ard,wsp,rid,awp,wvp
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
       iwat=1
       iadia=0
       print 33,a1
   33  format(1h ,' Watanabe deuteron potential for a = ',f5.1)
       print 101,z1,energy
       print*,' from one of the nucleon potentials:    '
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=1.d0
       rsoid=99.d0
       asoid=1.d0
       print*
 701   print*,' Deuteron wave function for Watanabe folding '
       print*,'----------------------------------------------------'
       print*,' [1] Reid soft-core wfn and interaction      '
       print*,' [2] Hulthen S+D-state wave function         '
       print*,' [3] Hulthen S-state only wave function      '
       print*,' [4] AV18 deuteron wave function             '
       print*,'----------------------------------------------------'
       read*,ideutwf
       if(ideutwf.lt.1.or.ideutwf.gt.4) go to 701
       open(77,file='deutwf.'//fname,status='unknown')
       write(18,*) ideutwf
       print*,' >>>> ',ideutwf
       if(ideutwf.eq.1) then
        print*,' Reid soft-core wfn and interaction      '
       else if(ideutwf.eq.2) then
        print*,' Hulthen S+D-state wave function         '
       else if(ideutwf.eq.3) then
        print*,' Hulthen S-state only wave function      '
       else if(ideutwf.eq.4) then
        print*,' AV18 deuteron wave function             '
       endif
      endif
*---------------------------------------------------------------------
      if(ireac.lt.5.and.inopt.gt.4) then
        if(inonloc.eq.2) then
         print*,'----------------------------------------------------'
         print*,' Conventional Perey-Buck type nonlocality was chosen'
         print*,' plus an adiabatic deuteron channel description     '
         print*,' This is not well founded - and NOT recommended     '
         print*,'----------------------------------------------------'
        endif
        iadia=inopt
       iw=0
       if(iadia.eq.5) then
        print 27,a1
   27   format(1h ,' ZR Adiabatic deuteron potential for a = ',f5.1)
       else 
        print 57,a1
   57   format(1h ,' FR Adiabatic deuteron potential for a = ',f5.1)
       endif
       print 101,z1,energy
       print*,' from one of the nucleon potentials:    '
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=1.d0
       rsoid=99.d0
       asoid=1.d0
       if(iadia.eq.6) then
        print*
 702    print*,' Deuteron wave function for adiabatic folding'
        print*,'----------------------------------------------------'
        print*,' [1] Reid soft-core wfn and interaction      '
        print*,' [2] Hulthen S+D-state wave function         '
        print*,' [3] Hulthen S-state only wave function      '
        print*,' [4] AV18 deuteron wave function             '
        print*,'----------------------------------------------------'
        read*,ideutwf
        if(ideutwf.lt.1.or.ideutwf.gt.4) go to 702
        open(77,file='deutwf.'//fname,status='unknown')
        write(18,*) ideutwf
        print*,' >>>> ',ideutwf
       if(ideutwf.eq.1) then
        print*,' Reid soft-core wfn and interaction      '
       else if(ideutwf.eq.2) then
        print*,' Hulthen S+D-state wave function         '
       else if(ideutwf.eq.3) then
        print*,' Hulthen S-state only wave function      '
       else if(ideutwf.eq.4) then
        print*,' AV18 deuteron wave function             '
       endif
       endif
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine adiab(energy,a,z,step,nrmax)
      implicit real*8(a-h,o-z)
      real*8 prn(900),pin(900),prp(900),pip(900),pson(900),psop(900)
      real*8 vreal(900),vimag(900),vspin(900),vspii(900),pgrid(900)
      character form*12,nucleon,fname*12
      real*8 psin(900),psip(900)
      common/folded/prn,pin,pson,psin,prp,pip,psop,psip,pgrid
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/cjlm/potr(900),poti(900),ijlm,ijlmp,ikd
      common/fromfold/vreal,vimag,vspin,vspii
      common/deut/iadia,iwat,ideutwf,ii,imso
      common/djlm/rlr,rli
      common/file/fname
*---------------------------------------------------------------------
*     woods-saxon potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is coefficient of l.sigma for nucleons (ww~6MeV)
*---------------------------------------------------------------------
*     we are here if iadia > 0 or iwat > 0
*---------------------------------------------------------------------
      open(36,file='folded.'//fname,status='unknown')
      do ii=1,900
       prn(ii)= 0.d0
       pin(ii)= 0.d0
       pson(ii)=0.d0
       psin(ii)=0.d0
       prp(ii)= 0.d0
       pip(ii)= 0.d0
       psop(ii)=0.d0
       psip(ii)=0.d0
       vreal(ii)=0.d0
       vimag(ii)=0.d0
       vspin(ii)=0.d0
       vspii(ii)=0.d0
      enddo
      form='(5e14.7)'
      imso=0
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (A>40 20<E<50 MeV)       '
      print*,'            Phys Rev 182 (1969) 1190             '
      print*,' [2] Chapel-Hill 89 Global set (A>40 E>10 MeV)   '
      print*,'            Phys Rep 201 (1991) 57               '
      print*,' [3] JLM microscopic optical potential           '
      print*,'            Bauge implementation PRC 58, 1120    '
      print*,' [4] Koning-Delaroche global potential           '
      print*,'            Nucl Phys A713 (2003) 231            '
      print*,'----------------------------------------------------'
*---------------------------------------------------------------------
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.4) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
*---------------------------------------------------------------------
      if(inopt.eq.1) then
*---------------------------------------------------------------------
*      use half the deuteron energy
*---------------------------------------------------------------------
       ed=energy/2.d0
       e=z
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=54.d0-0.32d0*ed+0.4d0*e/a13+24.d0*an
       vnr=56.3d0-0.32d0*ed-24.d0*an
       rr=1.17d0
       ar=0.75d0
       rpi=1.32d0
       rni=1.26d0
       api=0.51d0+0.7d0*an
       ani=0.58d0
       wpi=11.8d0-0.25d0*ed+12.d0*an
       wni=13.0d0-0.25d0*ed-12.d0*an
       if(wpi.lt.0.d0) wpi=0.d0
       if(wni.lt.0.d0) wni=0.d0
       wvpi=0.22d0*ed-2.7d0
       wvni=0.22d0*ed-1.56d0
       if(wvpi.lt.0.d0) wvpi=0.d0
       if(wvni.lt.0.d0) wvni=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV nucleon energy ')
       print 11
       print 111
   11  format(1h ,'   proton ')
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
       print 19
       print 191
       print*
   19  format(1h ,'   neutron')
  191  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vnr,rr,ar,wni,rni,ani,wvni
   12  format(1h ,7f7.3)
       vso=6.2d0
       rso=1.01d0
       aso=0.75d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vso,rso,aso
*---------------------------------------------------------------------
       rr =rr *a13
       rpi=rpi*a13
       rni=rni*a13
       rso=rso*a13
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        prn(ii)=ws(r,vnr,rr,ar)
        pin(ii)=ws(r,wvni,rni,ani)+wsd(r,wni,rni,ani)
        pson(ii)=wso(r,vso,rso,aso)
        prp(ii)=ws(r,vpr,rr,ar)
        pip(ii)=ws(r,wvpi,rpi,api)+wsd(r,wpi,rpi,api)
        psop(ii)=pson(ii)
*       if zero-range adiabatic then construct coincidence potentials
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
*        mean coincidence n and p spin-orbit form factor
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif 
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspii(ii),ii=1,nrmax)
        print*,'potential is written to tran and folded '
        do ii=1,nrmax
         r=ii*step
         write(36,1077) r,vreal(ii),vimag(ii),vspin(ii),vspii(ii)
        enddo
       endif
      endif 
*---------------------------------------------------------------------
      if(inopt.eq.2) then
       e=energy/2.d0
       n=nint(a-z)
*---------------------------------------------------------------------
*      CH86 parameters
       v0=    52.9d0
       vt=    13.d0
       ve=   -0.3d0
       r0=    1.25d0
       r00=  -0.24d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   10.d0
       wve0=  35.d0
       wvew=  15.d0
       rw=    1.32d0
       rw0=  -0.41d0
       aw=    0.72d0
       ws0=   9.d0
       wst=   14.d0
       wse0=  29.d0
       wsew=  23.d0
*---------------------------------------------------------------------
*      CH89 parameters
       v0=    52.9d0
       vt=    13.1d0
       ve=   -0.299d0
       r0=    1.25d0
       r00=  -0.225d0
       a0=    0.69d0
       rc=    1.24d0
       rc0=   0.12d0
       wv0=   7.8d0 
       wve0=  35.d0
       wvew=  16.d0
       ws0=   10.d0
       wst=   18.d0
       wse0=  36.d0
       wsew=  37.d0
       rw=    1.33d0
       rw0=  -0.42d0
       aw=    0.69d0
*---------------------------------------------------------------------
       a13=a**(1.d0/3.d0)
       rrc=rc*a13+rc0
       rcn=rrc/a13
       ecpp=1.73d0*z/rrc
       ecnn=0.d0
       erp=e-ecpp
       ern=e-ecnn
       vrp=v0+vt*((n-z)/a)+erp*ve
       vrn=v0-vt*((n-z)/a)+ern*ve
       rp=r0*a13+r00
       rpn=rp/a13
       rn=rpn
       ap=a0
       an=ap
       wvp=wv0/(1.d0+exp((wve0-erp)/wvew))
       wvn=wv0/(1.d0+exp((wve0-ern)/wvew))
       if(wvp.lt.0.d0) wvp=0.d0
       if(wvn.lt.0.d0) wvn=0.d0
       rwp=rw*a13+rw0
       rwpn=rwp/a13
       rwn=rwpn
       awp=aw
       awn=awp
       wsp=(ws0+wst*((n-z)/a))/(1.d0+exp((erp-wse0)/wsew))
       wsn=(ws0-wst*((n-z)/a))/(1.d0+exp((ern-wse0)/wsew))
       if(wsp.lt.0.d0) wsp=0.d0
       if(wsn.lt.0.d0) wsn=0.d0
       print 20,a
       print 101,z,e
   20  format(1h ,' Chapel Hill 89 potentials for a = ',f5.1)
       print 11
       print 111
       print 12,vrp,rpn,ap,wsp,rwpn,awp,wvp
       print *
       print 19
       print 191
       print 12,vrn,rn,an,wsn,rwn,awn,wvn
*---------------------------------------------------------------------
*      CH86 parameters
       vso=5.9d0
       rso=(1.39d0*a13-1.43)/a13
       aso=0.65d0
*---------------------------------------------------------------------
*      CH89 parameters
       vso=5.9d0
       rso=(1.34d0*a13-1.20)/a13
       aso=0.63d0
*---------------------------------------------------------------------
       print 112
       print 12,vso,rso,aso
*---------------------------------------------------------------------
       rr =rpn*a13
       rpi=rwpn*a13
       rni=rwn*a13
       rso=rso*a13
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        prn(ii)=ws(r,vrn,rr,ap)
        pin(ii)=ws(r,wvn,rni,awn)+wsd(r,wsn,rni,awn)
        pson(ii)=wso(r,vso,rso,aso)
        prp(ii)=ws(r,vrp,rr,ap)
        pip(ii)=ws(r,wvp,rpi,awp)+wsd(r,wsp,rpi,awp)
        psop(ii)=pson(ii)
*       if zero-range adiabatic then construct coincidence potentials
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
*        mean coincidence n and p spin-orbit form factor
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspii(ii),ii=1,nrmax)
        print*,'potential is written to tran and folded '
        do ii=1,nrmax
         r=ii*step
         write(36,1077) r,vreal(ii),vimag(ii),vspin(ii),vspii(ii)
        enddo
       endif
      endif
*---------------------------------------------------------------------
      if(inopt.eq.3) then
       if(ijlm.eq.0) ijlm=1
       ed=energy/2.d0
       nucleon='n'
       call jlm(prn,pin,nucleon,a,z,ed,step,nrmax)
       nucleon='p'
       call jlm(prp,pip,nucleon,a,z,ed,step,nrmax)
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
*       no spin-orbit currently in jlm option
        pson(ii)=0.d0
        psop(ii)=pson(ii)
*       if zero-range adiabatic then construct coincidence potentials
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
*        mean coincidence n and p spin-orbit form factor
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
        endif 
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspii(ii),ii=1,nrmax)
        print*,'potential is written to tran and folded '
        do ii=1,nrmax
         r=ii*step
         write(36,1077) r,vreal(ii),vimag(ii),vspin(ii),vspii(ii)
        enddo
       endif
      endif
*---------------------------------------------------------------------
      if(inopt.eq.4) then
*      KD02 systematics
       e=energy/2.d0
       n=nint(a-z)
       imso=1  
*---------------------------------------------------------------------
*      first the proton
       rvKD=1.3039-0.4054*A**(-1./3.)
       avKD=0.6778-1.487e-4*A
       rwKD=rvKD
       awKD=avKD
       v4KD=7.0e-9
       w2KD=73.55+0.0795*A
       rvdKD=1.3424-0.01585*A**(1./3.)
       rwdKD=rvdKD
       vdKD=0.
       d2KD=0.0180+3.802e-3/(1.+exp((A-156.)/8.))
       d3KD=11.5
       vso1KD=5.922+0.0030*A
       vso2KD=0.0040
       rvsoKD=1.1854-0.647*A**(-1./3.)
       rwsoKD=rvsoKD
       avsoKD=0.59
       awsoKD=avsoKD
       wso1KD=-3.1
       wso2KD=160.
       efKD=-8.4075+0.01378*A
       v1KD=59.30+21.0*real(N-Z)/A-0.024*A
       v2KD=7.067e-3+4.23e-6*A
       v3KD=1.729e-5+1.136e-8*A
       w1KD=14.667+0.009629*A
       avdKD=0.5187+5.205e-4*A
       awdKD=avdKD
       d1KD=16.0+16.0*real(N-Z)/A
       rcKD=1.198+0.697*A**(-2./3.)+12.994*A**(-5./3.)
*      proton energy-dependence part
       fKD=E-efKD
*      Coulomb
       VcKD=1.73/rcKD*Z/(A**(1./3.))
       vcoulKD=VcKD*v1KD*(v2KD-2.*v3KD*fKD+3.*v4KD*fKD*fKD)
       vKD=v1KD*(1.-v2KD*fKD+v3KD*fKD**2-v4KD*fKD**3)+vcoulKD
       wKD=w1KD*fKD**2/(fKD**2+w2KD**2)
       vdKD=0.
       wdKD=d1KD*fKD**2*exp(-d2KD*fKD)/(fKD**2+d3KD**2)
       vsoKD=vso1KD*exp(-vso2KD*fKD)
       wsoKD=wso1KD*fKD**2/(fKD**2+wso2KD**2)
       vkdp  = vKD
       rvkdp = rvKD
       avkdp = avKD
       wkdp  = wKD
       rwkdp = rvkdp
       awkdp = avkdp
       wskdp  = wdKD
       rwskdp = rwdKD
       awskdp = awdKD
       vsokdp = vsoKD
       wsokdp = wsoKD
       rsokdp = rvsoKD
       asokdp = avsoKD
*      neutron energy-dependence part
       efKD=-11.2814+0.02646*A
       v1KD=59.30-21.0*real(N-Z)/A-0.024*A
       v2KD=7.228e-3-1.48e-6*A
       v3KD=1.994e-5-2.0e-8*A
       w1KD=12.195+0.0167*A
       d1KD=16.0-16.0*real(N-Z)/A
       avdKD=0.5446-1.656e-4*A
       awdKD=avdKD
       rcKD=0.
       fKD=E-efKD
*      Coulomb
       vcoulKD=0.
       vKD=v1KD*(1.-v2KD*fKD+v3KD*fKD**2-v4KD*fKD**3)+vcoulKD
       wKD=w1KD*fKD**2/(fKD**2+w2KD**2)
       vdKD=0.
       wdKD=d1KD*fKD**2*exp(-d2KD*fKD)/(fKD**2+d3KD**2)
       vsoKD=vso1KD*exp(-vso2KD*fKD)
       wsoKD=wso1KD*fKD**2/(fKD**2+wso2KD**2)
       vkdn  = vKD
       rvkdn = rvKD
       avkdn = avKD
       wkdn  = wKD
       rwkdn = rvkdn
       awkdn = avkdn
       wskdn  = wdKD
       rwskdn = rwdKD
       awskdn = awdKD
       vsokdn = vsoKD
       wsokdn = wsoKD
       rsokdn = rvsoKD
       asokdn = avsoKD
*      print potential parameters
       print 202,a
       print 101,z,e
  202  format(1h ,' KD02 systematics for a = ',f5.1)
*--------------------------------------------------------------------
       print 11
       print 1112
 1112  format(1h ,'    v      rv     av     w      rw     aw ')
       print 122, vkdp,rvkdp,avkdp,wkdp,rwkdp,awkdp
       print 1312
 1312  format(1h ,'    ws     rws    aws  ')
       print 122, wskdp,rwskdp,awskdp
       print 1129
 1129  format(1h ,'  vso    rso    aso    wso    rsoi   asoi ')
       print 12,vsokdp,rsokdp,asokdp,wsokdp,rsokdp,asokdp
*--------------------------------------------------------------------
       print *
       print 19
       print 1112
       print 122, vkdn,rvkdn,avkdn,wkdn,rwkdn,awkdn
       print 1312
       print 122, wskdn,rwskdn,awskdn
       print 1129
       print 12,vsokdn,rsokdn,asokdn,wsokdn,rsokdn,asokdn
  122  format(1h ,9f7.3)
*---------------------------------------------------------------------
       a13 = a**(1./3.)
       rrvp =rvkdp*a13
       rrwp =rwkdp*a13
       rrwsp=rwskdp*a13
       rrvn =rvkdn*a13
       rrwn =rwkdn*a13
       rrwsn=rwskdn*a13
       rrsop=rsokdp*a13
       rrson=rsokdn*a13
*      for /folded/prn,pin,pson,prp,pip,psop
       do ii=1,nrmax
        r=ii*step
        pgrid(ii)=r
        prn(ii)=ws(r,vkdn,rrvn,avkdn)
        pin(ii)=ws(r,wkdn,rrwn,awkdn)+wsd(r,wskdn,rrwsn,awskdn)
        pson(ii)=wso(r,vsokdn,rrson,asokdn)
*       imaginary spin-orbit included
        psin(ii)=wso(r,wsokdn,rrson,asokdn)
        prp(ii)=ws(r,vkdp,rrvp,avkdp)
        pip(ii)=ws(r,wkdp,rrvp,awkdp)+wsd(r,wskdp,rrwsp,awkdp)
        psop(ii)=wso(r,vsokdp,rrsop,asokdp)
*       imaginary spin-orbit included
        psip(ii)=wso(r,wsokdp,rrsop,asokdp)
*       if zero-range adiabatic then construct coincidence potentials
        if(iadia.eq.5) then
         vreal(ii)=prn(ii)+prp(ii)
         vimag(ii)=pin(ii)+pip(ii)
*        mean coincidence n and p spin-orbit form factor
         vspin(ii)=(pson(ii)+psop(ii))/2.d0
         vspii(ii)=(psin(ii)+psip(ii))/2.d0
        endif
       enddo
       if(iadia.eq.5) then
        print*,'printout is in',nrmax,' steps of',real(step)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
        write(17,'(a)') form
        write(17,'(5e14.7)')(vspii(ii),ii=1,nrmax)
        print*,'potential is written to tran and folded '
        do ii=1,nrmax
         r=ii*step
         write(36,1077) r,vreal(ii),vimag(ii),vspin(ii),vspii(ii)
        enddo
       endif
      endif
 1077 format(2x,f5.2,4(2x,d12.5))
*---------------------------------------------------------------------
*     if finite range folding to be done, nucleon potentials to
*     be used have already been calculated above. call folder
*---------------------------------------------------------------------
      if(iwat.eq.1.or.iadia.eq.6) then
       call folder(step,nrmax)
*      print the returned folded potentials
       print*,'printout is in',nrmax,' steps of',real(step)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vreal(ii),ii=1,nrmax)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vimag(ii),ii=1,nrmax)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vspin(ii),ii=1,nrmax)
       write(17,'(a)') form
       write(17,'(5e14.7)')(vspii(ii),ii=1,nrmax)
       print*,'the potential was written to tran and folded '
      endif 
*---------------------------------------------------------------------
      return
      end
*---------------------------------------------------------------------
      subroutine folder(step,nrmax)
      implicit real*8(a-h,o-z)
      real*8 prn(900),pin(900),prp(900),pip(900),pson(900),psop(900)
      real*8 pgrid(900),br(12),wri(200),xri(200),sint(10),pint(13)
      real*8 g(8,400),vtr(8,900),funct(13,200),psin(900),psip(900)
      common/fromfold/vreal(900),vimag(900),vspin(900),vspii(900)
      common/avstuff/uavs(1501),vavs(1501),uavd(1501),vavd(1501)
      common/folded/prn,pin,pson,psin,prp,pip,psop,psip,pgrid
      common/deut/iadia,iwat,ideutwf,k,imso
      common/avstuf2/uasp(1501),uadp(1501)
*     parameter(ndim=900,ndata=201,ndatar=200,mmts=96,sr2=sqrt(2.d0))
      parameter(ndim=900,ndata=201,mmts=96,sr2=sqrt(2.d0))
      data vtr/7200*0.d0/
c     ----------------------------------------------------------------
c     icalc=1 watanabe : icalc=2 Johnson and Tandy adiabatic 
c     ----------------------------------------------------------------
      icalc=1
      if(iadia.gt.0) icalc=2
      ndatar=180
c     ----------------------------------------------------------------
c     gauss quadrature points - number of and weights arrays
c     if av18 choice then set up the wave function arrays needed
c     ----------------------------------------------------------------
      call gauss(-1.d0,1.d0,mmts,xri,wri)
      if(ideutwf.eq.4) call av18wf
c     ------------------------------------------------------------------
c     small r - internal ndata,sep  big R - external ndatar,step
c     adjust ndatar according to the chosen reaction step length
c     ------------------------------------------------------------------
      ndatar=nint(ndatar*0.1d0/step)
      if(ndatar.gt.nrmax) ndatar=nrmax
      sep=0.1d0
      if(icalc.eq.2) sep=0.05d0
      do 1001 k=1,ndatar
      rcap=k*step
      j21=2
      do ii=1,10
       sint(ii)=0.0
      enddo
c     ----------------------------------------------------------------
      do 1002 j=1,ndata
      r=(j-1)*sep
      if(j.eq.1) r=1.d-3
      rr1=rcap*rcap+0.25d0*r*r
      rr2=rcap*r
      do ii=1,13
       pint(ii)=0.0
      enddo
c     ----------------------------------------------------------------
c     loop on angular gauss quadrature points - cosine(theta)
      do i=1,mmts
       um=xri(i)
       p2=0.5d0*(3.d0*um*um-1.d0)
c     ----------------------------------------------------------------
c      will use function terp(r,fun,rgrid,npts,ndim) to interpolate
c     ----------------------------------------------------------------
       r1=sqrt(rr1+rr2*um)
       gr= terp(r1,prn ,pgrid,nrmax,ndim)
       gi= terp(r1,pin ,pgrid,nrmax,ndim)
       gso=terp(r1,pson,pgrid,nrmax,ndim)
       r2=sqrt(rr1-rr2*um)
       gr= gr +terp(r2,prp ,pgrid,nrmax,ndim)
       gi= gi +terp(r2,pip ,pgrid,nrmax,ndim)
       gso=gso+terp(r2,psop,pgrid,nrmax,ndim)
       gsi=0.d0
       if(imso.gt.0) then
        gsi=terp(r1,psin,pgrid,nrmax,ndim)
        gsi=gsi+terp(r2,psip,pgrid,nrmax,ndim)
       endif
c     ----------------------------------------------------------------
c      angle weighted integrands of the different potential terms
c     ----------------------------------------------------------------
       funct(1,i)= gr
       funct(2,i)= gi
       funct(3,i)= 0.d0 
       funct(4,i)= gr*p2
       funct(5,i)= gi*p2
       funct(6,i)= gso
       funct(7,i)= gso*um
       funct(8,i)= gso*p2
       funct(9,i)= gso*um*abs(1.d0-um*um)
       funct(10,i)=gsi
       funct(11,i)=gsi*um
       funct(12,i)=gsi*p2
       funct(13,i)=gsi*um*abs(1.d0-um*um)
c     ----------------------------------------------------------------
       do ii=1,13
        pint(ii)=pint(ii)+funct(ii,i)*wri(i)
       enddo
      enddo
c     ----------------------------------------------------------------
c     ends loop over angles gauss quadrature points if given (r,R)
c     ----------------------------------------------------------------
      call fact(r,br,icalc)
c     central  
c     with added nucleon spin-orbit potential contributions
      g(1,j)=br(1)*pint(1)-1.5d0*br(6)*pint(6)
      g(1,j)=g(1,j)-3.d0*(rcap/r)*br(6)*pint(7)
      g(1,j)=g(1,j)/2.d0
c     added spin-orbit terms in imaginary part
      g(2,j)=br(1)*pint(2)-1.5d0*br(6)*pint(10)
      g(2,j)=g(2,j)-3.d0*(rcap/r)*br(6)*pint(11)
      g(2,j)=g(2,j)/2.d0
c     tensor tr  
c     with added nucleon spin-orbit potential contributions
      g(4,j)=3.d0/sr2*br(2)*pint(4)
      g(4,j)=g(4,j)-4.5d0*(rcap/r)*br(7)*pint(7)
      g(4,j)=g(4,j)+4.50*rcap*br(8)*pint(9)
      g(4,j)=g(4,j)-4.50*br(9)*pint(8)
      g(4,j)=g(4,j)/2.d0
c     with added nucleon spin-orbit potential contributions
      g(5,j)=3.d0/sr2*br(2)*pint(5)
      g(5,j)=g(5,j)-4.5d0*(rcap/r)*br(7)*pint(11)
      g(5,j)=g(5,j)+4.50*rcap*br(8)*pint(13)
      g(5,j)=g(5,j)-4.50*br(9)*pint(12)
      g(5,j)=g(5,j)/2.d0
c     spin-orbit  
      g(6,j)=(br(5)-br(6)/2.0)*pint(6)/2.d0
c     g(6,j)=g(6,j)+br(3)*(r/(2.0*rcap))*pint(7)/2.d0
c     g(6,j)=g(6,j)-br(4)*pint(8)/2.0
c     overall L.S factor of one-half in Keaton + Armstrong
      g(6,j)=g(6,j)/2.d0
c     imaginary spin-orbit, use(g(3,*)  
      g(3,j)=0.d0
      if(imso.gt.0) then
       g(3,j)=(br(5)-br(6)/2.0)*pint(10)/2.d0
       g(3,j)=g(3,j)/2.d0
      endif
c     adiabatic model denominator 
      g(7,j)=br(1)
c     s-state part of the denominator to calculate prob_s and prob_d
      g(8,j)=br(5)
c     ----------------------------------------------------------------
      if(j-j21.le.0) go to 1002
      j1=j-1
      j2=j-2
c     ----------------------------------------------------------------
      do ii=1,8
       sint(ii)=sint(ii)+sep*(g(ii,j2)+4.d0*g(ii,j1)+g(ii,j))/3.0
      enddo
c     ----------------------------------------------------------------
      j21=j+1
 1002 continue
c     ----------------------------------------------------------------
      do ii=1,6
       vtr(ii,k)=sint(ii)/sint(7)
      enddo
c     ----------------------------------------------------------------
c     potentials to usual arrays for writing to tran data set
c     ----------------------------------------------------------------
      vreal(k)=vtr(1,k)
      vimag(k)=vtr(2,k)
      vspin(k)=vtr(6,k)
      vspii(k)=vtr(3,k)
c     ----------------------------------------------------------------
c     write out the folded potentials in the order real central, imag
c     central, real and imag spin orbit (coefficient of L.S), real TR 
c     tensor and imag TR tensor. The tensor interactions are not used 
c     by twofnr at this time
c     ----------------------------------------------------------------
      write(36,17) rcap,vreal(k),vimag(k),vspin(k),vspii(k),
     #                vtr(4,k),vtr(5,k)
  17  format(2x,f5.2,6(2x,d12.5))
 1001 continue
c     ----------------------------------------------------------------
c     print the s- and d-state probabilities in the Watanabe or the
c     Johnson and Tandy (JT) approaches. For Watanabe these are the 
c     wfn probabilities, for JT adiabatic, those of <phi|Vnp|phi>
c     ----------------------------------------------------------------
      probs=100*sint(8)/sint(7)
      print*,'----------------------------------------------------'
      if(ideutwf.eq.1) then
       print*,' Reid soft-core wfn and interaction      '
      else if(ideutwf.eq.2) then
       print*,' Hulthen S+D-state wave function         '
      else if(ideutwf.eq.3) then
       print*,' Hulthen S-state only wave function      '
      else if(ideutwf.eq.4) then
       print*,' AV18 deuteron wave function             '
      endif
      if(icalc.eq.1) then
       print*,' Deuteron s and d-state probabilities (%): '
      else 
       print*,' Johnson-Tandy s and d-state probabilities (%): '
      endif
      print*,' Ps = ',real(probs),' Pd = ',real(100.0-probs)
      if(icalc.eq.2) then
       print*,' Johnson-Tandy denominator: '
       print*,' <phi|Vnp|phi> = ',real(sint(7)),' MeV'
      endif
      print*,'----------------------------------------------------'
      return
      end
c-----------------------------------------------------------------------
      subroutine fact(r,br,icalc)
      implicit real*8(a-h,o-z)
      real*8 br(12)
      common/avstuff/uavs(1501),vavs(1501),uavd(1501),vavd(1501)
      common/avstuf2/uasp(1501),uadp(1501)
      common/deut/iadia,iwat,ideutwf,k,imso
      parameter (h=0.01d0,twoh=2.d0*h,sr2=sqrt(2.d0),sr8=2.d0*sr2)
      z=0.d0
      if(r.le.1.d-3) goto 1
c     --------------------------------------------------
      if(ideutwf.eq.1)then
c     Reid soft core case
       u =us(r)
       w =ud(r)
       up=(us(r+h)-us(r-h))/twoh
       wp=(ud(r+h)-ud(r-h))/twoh
       if(icalc.ne.1)then
        v0=vs(r)
        v2=vd(r)
       endif
      endif
c     --------------------------------------------------
      if(ideutwf.eq.2)then
c     Hulthen s+d case
       u =uhs(r)
       w =uhd(r)
       up=(uhs(r+h)-uhs(r-h))/twoh
       wp=(uhd(r+h)-uhd(r-h))/twoh
       if(icalc.ne.1)then
        v0=vhs(r)
        v2=vhd(r)
       endif
      endif
c     --------------------------------------------------
      if(ideutwf.eq.3)then
c     Hulthen s-only case
       u =uhp(r)
       w =0.d0
       up=(uhp(r+h)-uhp(r-h))/twoh
       wp=0.d0
       if(icalc.ne.1)then
        v0=vhp(r)
        v2=0.d0
       endif
      endif
c     --------------------------------------------------
      if(ideutwf.eq.4) then
c     AV18 case
       irele=nint(r*100)+1
       if(irele.gt.1501) go to 1
       u =uavs(irele)
       w= uavd(irele)
       up=uasp(irele)
       wp=uadp(irele)
       if(icalc.ne.1)then
        v0=vavs(irele)
        v2=vavd(irele)
       endif
      endif
c     --------------------------------------------------
      if(icalc.eq.1) then
       v0=u
       v2=w
      endif
c     --------------------------------------------------
      if(k.eq.1) write(77,'(f7.2,4d18.9)') r,u,w,v0,v2
c     --------------------------------------------------
      br(1)=(v0*u+v2*w)
      br(2)=v0*w+v2*u-v2*w/sr2
      br(3)=v0*u-v2*w-(v0*w+v2*u)/sr8
      br(4)=(v0*w+v2*u)/sr8+w*v2/2.d0
      br(5)=v0*u
      br(6)=v2*w
      br(7)=sr2*v0*w-v2*w
      br(8)=((v2*up-v0*wp)-(v2*u-v0*w)/r)/sr2
      br(8)=br(8)+(sr2*v0*w-2.d0*v2*w)/r
      br(9)=br(7)/2.d0
c     --------------------------------------------------
      return
    1 do jt=1,12
       br(jt)=0.d0
      enddo
      if(k.eq.1) then
       write(77,'(f7.2,4d18.9)') z,z,z,z,z
      endif
      return
      end 
c-----------------------------------------------------------------------
      real*8 function us(r)
*----------------------------------------------------------------------
*     Reid soft core s-state wavefunction
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension x(33),y(33),v(33)
      data x,y/1.000d-2,4.125d-2,7.250d-2,1.350d-1,1.975d-1,
     1         2.600d-1,3.225d-1,3.850d-1,4.475d-1,5.100d-1,
     2         5.725d-1,6.350d-1,6.975d-1,7.600d-1,8.850d-1,
     3         1.010d-0,1.135d-0,1.260d-0,1.385d-0,1.510d-0,
     4         1.760d-0,2.010d-0,2.510d-0,3.010d-0,3.510d-0,
     5         4.010d-0,4.510d-0,5.010d-0,5.510d-0,6.010d-0,
     6         7.010d-0,8.010d-0,9.010d-0,
     7         0.0000d-0,3.3373d-5,2.3901d-4,2.7621d-3,1.2737d-2,
     8         3.6062d-2,7.5359d-2,1.2847d-1,1.8993d-1,2.5349d-1,
     9         3.1390d-1,3.6770d-1,4.1317d-1,4.4992d-1,4.9953d-1,
     a         5.2406d-1,5.3166d-1,5.2864d-1,5.1926d-1,5.0621d-1,
     b         4.7505d-1,4.4200d-1,3.7864d-1,3.2249d-1,2.7399d-1,
     c         2.3251d-1,1.9719d-1,1.6718d-1,1.4172d-1,1.2012d-1,
     d         8.6290d-2,6.1983d-2,4.4523d-2/
      data v/2.9751d-4,2.6127d-3,1.2335d-2,8.2951d-2,2.5326d-1,
     1  5.0052d-1,7.5072d-1,9.3349d-1,1.0162d-0,1.0034d-0,
     2  9.2042d-1,7.9674d-1,6.5744d-1,5.1985d-1,2.8494d-1,
     3  1.1865d-1,1.1397d-2,-5.4090d-2,-9.2523d-2,-1.1418d-1,
     4  -1.3097d-1,-1.3193d-1,-1.2004d-1,-1.0453d-1,-8.9709d-2,
     5  -7.6510d-2,-6.5054d-2,-5.5229d-2,-4.6850d-2,-3.9726d-2,
     6  -2.8547d-2,-2.0508d-2,-1.4731d-2/
      z=0.7*r
      if(z-0.01)10,10,20
   10 us=0.0d0
      return
   20 if(z-9.01) 40,40,30
   30 us=0.87758d0*exp(-.33087d0*z)
      return
   40 do 100 i=1,33
      if(x(i)-z) 70,80,100
   70 if(x(i+1)-z) 100,60,50
  100 continue
   80 us=y(i)
      return
   60 us=y(i+1)
      return
   50 h=x(i+1)-x(i)
      p=(z-x(i))/h
      y1=y(i)
      y2=y(i+1)
      v1=v(i)
      v2=v(i+1)
      a1=y1
      a2=h*v1
      a3=3.d0*(y2-y1)-h*(2.d0*v1+v2)
      a4=2.d0*(y1-y2)+h*(v1+v2)
      us=((a4*p+a3)*p+a2)*p+a1
      return
      end
c-----------------------------------------------------------------------
      real*8 function ud(r)
*----------------------------------------------------------------------
*     Reid soft core d-state wavefunction
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      dimension x(33),y(33),v(33)
      data x,y/1.000d-2,4.125d-2,7.250d-2,1.350d-1,1.975d-1,
     1         2.600d-1,3.225d-1,3.850d-1,4.475d-1,5.100d-1,
     2         5.725d-1,6.350d-1,6.975d-1,7.600d-1,8.850d-1,
     3         1.010d-0,1.135d-0,1.260d-0,1.385d-0,1.510d-0,
     4         1.760d-0,2.010d-0,2.510d-0,3.010d-0,3.510d-0,
     5         4.010d-0,4.510d-0,5.010d-0,5.510d-0,6.010d-0,
     6         7.010d-0,8.010d-0,9.010d-0,
     8         0.0000d-0,1.0850d-5,8.4073d-5,1.0369d-3,4.9642d-3,
     9         1.4446d-2,3.0795d-2,5.3157d-2,7.8995d-2,1.0525d-1,
     a         1.2933d-1,1.4958d-1,1.6529d-1,1.7645d-1,1.8710d-1,
     b         1.8654d-1,1.7946d-1,1.6910d-1,1.5742d-1,1.4553d-1,
     c         1.2314d-1,1.0373d-1,7.3859d-2,5.3293d-2,3.9077d-2,
     d         2.9115d-2,2.2016d-2,1.6871d-2,1.3079d-2,1.0243d-2,
     f         6.4412d-3,4.1575d-3,2.7363d-3/
      data v/7.4202d-5,8.9321d-4,4.4755d-3,3.1982d-2,1.0121d-1,
     1  2.0596d-1,3.1477d-1,3.9371d-1,4.2466d-1,4.0842d-1,
     2  3.5772d-1,2.8848d-1,2.1428d-1,1.4427d-1,3.3294d-2,
     3  -3.5899d-2,-7.3151d-2,-9.0091d-2,-9.5299d-2,-9.4158d-2,
     4  -8.3973d-2,-7.1282d-2,-4.9327d-2,-3.3940d-2,-2.3612d-2,
     5  -1.6691d-2,-1.2002d-2,-8.7768d-3,-6.5201d-3,-4.9136d-3
     6  ,-2.8956d-3,-1.7758d-3,-1.1227d-3/
      z=0.7*r
      if(z-0.01)10,10,20
   10 ud=0.0d0
      return
   20 if(z-9.01) 40,40,30
   30 z=0.33087d0*z
      w=1.d0/z
      ud=0.023013d0*exp(-z)*(1.d0+3.d0*w*(1.d0+w))
      return
   40 do 100 i=1,33
      if(x(i)-z) 70,80,100
   70 if(x(i+1)-z) 100,60,50
  100 continue
   80 ud=y(i)
      return
   60 ud=y(i+1)
      return
   50 h=x(i+1)-x(i)
      p=(z-x(i))/h
      y1=y(i)
      y2=y(i+1)
      v1=v(i)
      v2=v(i+1)
      a1=y1
      a2=h*v1
      a3=3.d0*(y2-y1)-h*(2.d0*v1+v2)
      a4=2.d0*(y1-y2)+h*(v1+v2)
      ud=((a4*p+a3)*p+a2)*p+a1
      return
      end
c-----------------------------------------------------------------------
      real*8 function vs(r)
*----------------------------------------------------------------------
*     Reid soft core s-state of Vnp*wavefunction
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      if(r.eq.0.0) r=1.0d-10
      vs=(vc(r)*us(r)+2.82842712*vt(r)*ud(r))
      return
      end
c-----------------------------------------------------------------------
      real*8 function vd(r)
*----------------------------------------------------------------------
*     Reid soft core d-state of Vnp*wavefunction
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      if(r) 10,10,20
   10 vd=0.0
      return
   20 vd=((vc(r)-2.*vt(r)-3.*vls(r))*ud(r)+2.82842712*vt(r)*us(r))
      return
      end
c-----------------------------------------------------------------------
      real*8 function vc(r)
      implicit real*8(a-h,o-z)
      x=0.7*r
      e=exp(-x)
      vc=e*(-10.463+e*(105.468+e*e*(-3187.8+9924.3*e*e)))/x
      return
      end
c-----------------------------------------------------------------------
      real*8 function vt(r)
      implicit real*8(a-h,o-z)
      x=0.7*r
      y=1.0/x
      e=exp(-x)
      vt=e*y*(-10.463*(1.+3.*y*(1.+y)-3.*y*(4.+y)*e**3)+
     1e**3*(351.77-1673.5*e*e))
      return
      end
c-----------------------------------------------------------------------
      real*8 function vls(r)
      implicit real*8(a-h,o-z)
      x=0.7*r
      e=exp(-x)
      vls=e**4*(708.91-2713.1*e*e)/x
      return
      end
*----------------------------------------------------------------------
      real*8 function uhs(r)
*----------------------------------------------------------------------
*     hulthen s-state wavefunction beta ~ 6*alpha (r space)
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (alpha=0.2314549d0, beta=1.3887294d0)
      parameter (rnorm=0.8474405d0)
      uhs=(exp(-alpha*r)-exp(-beta*r))*rnorm
      return
      end
*-----------------------------------------------------------------------
      real*8 function uhd(r)
*-----------------------------------------------------------------------
*     hulthen d-state wfunction gamma(t)~4alpha (r space) x~d/s ratio
*-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (alpha=0.2314549d0, t=4.d0*alpha)  
      parameter (rnorm=0.8474405d0, x=0.026d0)
      uhd=0.d0
      if(r.gt.1.d-12) then
      exg=1.d0-exp(-t*r)     
      exgsq=exg*exg  
      alrd=alpha*r
      alrdsq=alrd*alrd   
      uhd=(1.d0+3.d0*exg/alrd+3.d0*exgsq/alrdsq)*exgsq
      uhd=uhd*x*exp(-alrd)*rnorm
      endif
      return
      end
*------------------------------------------------------------------------
      real*8 function vhs(r)
*------------------------------------------------------------------------
*     hulthen s-state wavefunction * vnp
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (rmn=939.565d0,rmp=938.272d0, rmu=rmn*rmp/(rmn+rmp)) 
      parameter (ed=-2.22452d0, h=0.01d0, hc=197.32696010352811d0)
      parameter (alpha=0.2314549d0, beta=1.3887294d0)
      cons=-hc*hc/2.d0/rmu
      twel=12.d0*h*h
      d2=16.d0*(uhs(r+h)+uhs(r-h))-30.d0*uhs(r)-uhs(r+2*h)-uhs(r-2*h)
      vhs=ed*uhs(r)-cons*d2/twel 
      return
      end 
*------------------------------------------------------------------------
      real*8 function vhd(r)
*------------------------------------------------------------------------
*     hulthen d-state wavefunction * vnp
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z) 
      parameter (rmn=939.565d0,rmp=938.272d0, rmu=rmn*rmp/(rmn+rmp)) 
      parameter (ed=-2.22452d0, h=0.01d0, hc=197.32696010352811d0)
      parameter (alpha=0.2314549d0, beta=1.3887294d0)
      vhd=0.d0
      if(r.gt.1.d-12) then
      cons=-hc*hc/2.d0/rmu
      twel=12.d0*h*h
      d2=16.d0*(uhd(r+h)+uhd(r-h))-30.d0*uhd(r)-uhd(r+2*h)-uhd(r-2*h)
      vhd=ed*uhd(r)-cons*(d2/twel-6.d0*uhd(r)/r/r)
      endif
      return
      end 
*------------------------------------------------------------------------
      real*8 function uhp(r)
*----------------------------------------------------------------------
*     hulthen s-state wavefunction beta~6*alpha (r space)
*----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (alpha=0.2314549d0, beta=1.3887294d0)
      parameter (rnorm=0.8818672d0)
      uhp=(exp(-alpha*r)-exp(-beta*r))*rnorm
      return
      end
*-----------------------------------------------------------------------
      real*8 function vhp(r)
*------------------------------------------------------------------------
*     hulthen s-state wavefunction * vnp
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (rmn=939.565d0,rmp=938.272d0, rmu=rmn*rmp/(rmn+rmp)) 
      parameter (ed=-2.22452d0, h=0.01d0, hc=197.32696010352811d0)  
      parameter (alpha=0.2314549d0, beta=1.3887294d0)   
      cons=-hc*hc/2.d0/rmu
      twel=12.d0*h*h
      d2=16.d0*(uhp(r+h)+uhp(r-h))-30.d0*uhp(r)-uhp(r+2*h)-uhp(r-2*h)
      vhp=ed*uhp(r)-cons*d2/twel 
      return
      end
*---------------------------------------------------------------------
      subroutine alphap(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      dimension aji(100), e(100)
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
*---------------------------------------------------------------------
*     potential formfactor statement functions
*---------------------------------------------------------------------
*     ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
*     wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
*---------------------------------------------------------------------
      print*,'----------------------------------------------------'
 890  print*,' [1] Atzrott/Kumar et al. global alpha potential '
      print*,'       PRC 53 (1996) 1336 and NP A776 (2006) 105 '
      print*,'----------------------------------------------------'
*---------------------------------------------------------------------
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.1) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      a=a1
      a13=a**0.3333333333d0
      if(inopt.eq.1) then
*---------------------------------------------------------------------
      amass=a1
      az=z1
*     program alphaop
*     dimension aji(100), e(100)
c     calculation of real and imaginary potentials for alpha -
c     nucleus  systems
c     please send your comments to  kailas@magnum.barc.ernet.in
c     modified on 30.8.2002
c     real and imaginary potentials starting from volume integral
c     systematics, potential and slope of potential at radius
c     close to strong absorption radius
c     ref. mahaux,ngo and satchler nucl.,phys.a 446,1986,354
c     for dispersion correction
*     open(unit=5, file='input')
*     open(unit=6, file='output')
      pi=4.*atan(1.)
      ex=1.0
      elab=energy
      at13=amass**(1./3.)
*---------------------------------------------------------------------
c     read mass of target, z of target, 1st excited state of target(in mev)
c     and e lab (in mev)
*     read(5, 99001)amass, az, ex, elab
*9001 format(f5.0, f5.0, 2f10.3)
*     write(6, 99002)amass, az, ex, elab
*9002 format('amass ', f5.1, 5x, 'az ', f5.1, 5x, 'ex ', f10.3, '(mev)', 
*    &       5x, 'elab', f10.3, '(mev)')
*---------------------------------------------------------------------
c     r2.4 systematics for real and imaginary parts of potential
c     ar, ai are the diffuseness parameter values
      r24r=1.35*at13 + 2.55
      r24i=1.35*at13 + 2.14
      ar=0.76
      ai=0.60
c     eref is the reference energy for normalisation in making dispersion
c     relation calculation
      eref=140.01
      erefc=140.01*amass/(amass + 4.)
c     jref is the volume integral at the reference energy and is
c     calculated using the empirical relation given below
      jref=(224. - 0.98*erefc/amass**0.184 + 2.57*az/at13)
     &       *(1. + (2.05/at13))
      jref=-jref
c     imaginary part systematics is taken from a.shridahar et al.
c     phys. rev. c30,1760 (1984)
      ajii=32.8*(1. + 7.1/at13)
c     calculation of real and imaginary potential radii by fitting
c     r 2.4 and volume integral values at e=90 mev
      const=4.*pi/4./amass/3.
      const1=pi*ar
      const2=pi*ai
      ec=90.*amass/(amass + 4.)
c     bar is the coulomb barrier
      bar=1.44*2.*az/1.5/(4.**(1./3.) + amass**(1./3.))
c     alpha and beta are the energy coefficients for imaginary part
      alpha=0.05
      beta=0.309
      ebar=bar + ex
      jrref=(224. - 0.98*ec/amass**0.184 + 2.57*az/at13)
     &        *(1. + (2.05/at13))
      jiief=ajii*(1. - exp( - (90.-ebar)*alpha))
      j=0
      ror=1.0
 100  rr=ror*at13
      vor=jrref/const/rr**3./(1 + (const1/rr)**2.)
      r24rc=rr + ar*log((vor - 2.4)/2.4)
      diff=r24rc - r24r
      j=j + 1
      if(j.le.200)then
         if(abs(diff).gt.0.05)then
            ror=ror + 0.005
            goto 100
         endif
      endif
      roi=1.2
      j=0
 200  ri=roi*at13
      voi=jiief/const/ri**3./(1. + (const2/ri)**2.)
      r24ic=ri + ai*log((voi - 2.4)/2.4)
      diff=r24ic - r24i
      j=j + 1
      if(j.le.200)then
         if(abs(diff).gt.0.05)then
            roi=roi + 0.005
            goto 200
         endif
      endif
      ecm=elab*amass/(amass + 4.)
      ed=-elab + ebar
      edd=elab - ebar - 3.
      if(edd.le.0)then
         ajie=ajii*0.055*exp( - beta*ed)
      else
         ajie=ajii*(1. - exp( - (elab-ebar)*alpha))
      endif
c     matching of two forms of imaginary potentials is done
c     at ebar + 3 mev.
      n=0
      ee=ex
      emax=ebar + 3.
 300  n=n + 1
      e(n)=ee
      if(e(n).le.emax)then
         aji(n)=ajii*0.055*exp( - beta*( - e(n) + ebar))
         ee=ee + 2.
         goto 300
      else
 350     if(e(n).le.140.)then
            aji(n)=ajii*(1. - exp( - (e(n)-ebar)*alpha))
            n=n + 1
            ee=ee + 4.
            e(n)=ee
            goto 350
         else
            no=n
c           do 189 i=1,no-1
c           write(6,188) e(i), aji(i)
c189        continue
c188        format ( 'e',f10.3,'aji',f10.3)
c           assume the dispersion correction  to vanish at eref.
            dvef=0.
            dverf=0.
            do j=1, no - 3
               dj=e(j + 1) - e(j)
               an1=(elab - e(j))*log(abs((elab-e(j))/dj))
               an2=(elab - e(j + 1))*log(abs((elab-e(j+1))/dj))
               dve=((aji(j+1) - aji(j))/dj)*(an1 - an2)/pi
               anr1=(eref - e(j))*log(abs((eref-e(j))/dj))
               anr2=(eref - e(j + 1))*log(abs((eref-e(j+1))/dj))
               dver=((aji(j+1) - aji(j))/dj)*(anr1 - anr2)/pi
               dvef=dvef + dve
               dverf=dverf + dver
            enddo
            v0=jref - dverf
            ajre=dvef + v0
            ror=rr/at13
            roi=ri/at13
            vor=ajre/const/rr**3./(1 + (const1/rr)**2.)
            voi=ajie/const/ri**3./(1. + (const2/ri)**2.)
            vor=-vor
            ajre=-ajre
*     -------------------------------------------------------------
*     comment the original output statements
*     -------------------------------------------------------------
*           write(6, 99003)
*9003       format(2x, 'elab(mev) ', 2x, 'real v.i.(mev fm3) ', 2x, 
*    &             'img v.i.(mev  fm3)')
*           write(6, 99004)elab, ajre, ajie
*9004       format(f10.1, 9x, f10.1, 9x, f10.2)
*           write(6, 99005)
*9005       format(2x, 'elab(mev) ', 2x, 'vor(mev)  ', 2x, 'ror(fm)   ', 
*    &             2x, 'ar(fm)   ', 'voi(mev)  ', 2x, 'roi(fm)   ', 
*    &             'ai(fm)    ')
*           write(6, 99006)elab, vor, ror, ar, voi, roi, ai
*9006       format(f10.1, 2x, f10.3, 2x, f8.3, 2x, f8.3, 2x, f10.3, 2x, 
*    &             f8.3, 2x, f8.3)
*     -------------------------------------------------------------
*           print*,elab, vor, ror, ar, voi, roi, ai
         endif
       endif
*     -------------------------------------------------------------
*     convert to usual parameters
*     -------------------------------------------------------------
       vpr=vor
       rr=ror
       wvpi=voi
       wpi=0.d0
       rpi=roi
       api=ai
*     -------------------------------------------------------------
       print 10,a
   10  format(1h ,' alpha particle potentials for a=',f5.1)
       print 15,z1,energy
   15  format(1h ,' z = ',f4.1,' at ',f5.1,' mev alpha energy ')
       rcd=ror
       print*,' coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=0.0
       wsod=0.d0
       rsord=1.0d0
       asord=1.0d0
       rsoid=1.d0
       asoid=1.d0
      endif
      return
      end
*---------------------------------------------------------------------
      subroutine triton(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/che3/pspr(900),pspi(900),psps(900),ihesp
*---------------------------------------------------------------------
*     woods-saxon potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
      wso(r,ww,w0,wa)= wsd(r,ww,w0,wa)/(2.d0*wa*r)
*---------------------------------------------------------------------
*     latter is the coefficient of l.sigma 
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (not well determined)    '
      print*,'                      see: ADNDT 17 (1976) p6    ' 
      print*,'     Potential B2 of Nucl Phys A190 (1972) 437   '
      print*,' [2] X. Li et al. Global potential E<40 MeV      '    
      print*,'                      see: NPA 789 (2007) 103    '
      print*,' [3] D.Y. Pang et al. GDP08                      '
      print*,'                      see: PRC 79 (2009) 024615  '
      print*,'----------------------------------------------------'
*---------------------------------------------------------------------
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.3) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      ihesp=0
      a=a1
      ed=energy
      e=z1
      an=(a-2.d0*e)/a
      a13=a**0.3333333333d0
      z13=e**0.3333333333d0
      if(inopt.eq.1) then
       vpr=165.d0-0.17d0*ed-6.4d0*an
       rr=1.20d0
       ar=0.72d0
       rpi=1.40d0
       api=0.84d0
       wvpi=46.0d0-0.33d0*ed-110.d0*an
       if(wvpi.lt.0.d0) wvpi=0.d0
       wpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees 3h potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV triton energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=2.5d0
       rsord=1.20d0
       asord=0.72d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
      if(inopt.eq.2) then
*---------------------------------------------------------------------
*     X Li et al. Global potential NPA 2007
*---------------------------------------------------------------------
       ihesp=1
       do i=1,900
        pspr(i)=0.d0
        pspi(i)=0.d0
        psps(i)=0.d0
       enddo
       print 29,a1
   29  format(1h ,' Li et al. triton potential for a = ',f5.1)
       rcd=1.4219d0
       print*,' Coulomb radius parameter = ',real(rcd)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=1.d0
       rsord=99.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
*---------------------------------------------------------------------
       vlir=137.6d0-0.1456d0*ed+0.0436d0*ed*ed+4.3751d0*an
       vlir=vlir+1.0474d0*e/a13
       wlis=37.06d0-0.6451d0*ed-47.19d0*an
       wliv=7.383d0+0.5025d0*ed-0.0097d0*ed*ed
       vliso=1.9029d0
       rlir= (1.12010d0-0.1504d0/a13)
       rliv= (1.32020d0-0.1776d0/a13)
       rlis= (1.25100d0-0.4622d0/a13)
       rliso=(0.46991d0+0.1294d0/a13)
       alir= (0.68330d0+0.01910d0*a13)
       aliv= (1.11900d0+0.01913d0*a13)
       alis= (0.81140d0+0.01159d0*a13)
       aliso=(0.35450d0-0.05220d0*a13)
       print 117
  117  format(1h ,'    vr     ro     ao     ws     ris     ais  ')
       print 12,vlir,rlir,alir,wlis,rlis,alis
       print 119
  119  format(1h ,'    wv     riv     aiv   ')
       print 12,wliv,rliv,aliv
       print 118
  118  format(1h ,'    vso    rso    aso ')
       print 12,vliso,rliso,aliso
       rlir= rlir*a13
       rliv= rliv*a13
       rlis= rlis*a13
       rliso=rliso*a13
       print*,' printout is in',nrmax,' steps of',real(step)
*---------------------------------------------------------------------
       do ii=1,nrmax
        r=ii*step
        pspr(ii)=ws(r,vlir,rlir,alir)
        pspi(ii)=ws(r,wliv,rliv,aliv)+wsd(r,wlis,rlis,alis)
        psps(ii)=wso(r,vliso,rliso,aliso)
       enddo
      endif
      if(inopt.eq.3) then
*---------------------------------------------------------------------
*     GDP08 potential for 3H dypang
*---------------------------------------------------------------------
       a=a1
       ed=energy
       e=z1
       ap=3.0d0
       zp=1.0d0
       call gdp08(ap,zp,a1,z1,energy)
       print 710,a
       print 701,e,ed
  710  format(1h ,' DYPang et al 3H potentials for a = ',f5.1)
  701  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV triton energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 711
  711  format(1h ,'    vr     ro     ao     wv     riv    aiv ')
       print 713,vd,rrd,ard,wd,rid,aid
       print 715
  715  format(1h ,'    ws     ris    ais    ')
       print 713,wisd,risid,aisid
  713  format(1h ,7f7.3)
       print 712
  712  format(1h ,'    vso    rso    aso ')
       print 713,vsod,rsord,asord
       print 714
  714  format(1h ,'    wso   wrso   waso ')
       print 713,wsod,rsoid,asoid
      endif
*---------------------------------------------------------------------
      return
      end
*---------------------------------------------------------------------
      subroutine helium(energy,a1,z1,step,nrmax)
      implicit real*8(a-h,o-z)
      character spfile*30,guff*10
      real*8 funr(1500),funi(1500),rgrid(1500)
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/che3/pspr(900),pspi(900),psps(900),ihesp
      character form*12
*---------------------------------------------------------------------
*     potential formfactor statement functions
*---------------------------------------------------------------------
      ws (r,vv,r0,aa)=-vv/(1.d0+dexp((r-r0)/aa))
      wsd(r,ww,w0,wa)=-ww*4.d0*dexp((r-w0)/wa)/(1.+dexp((r-w0)/wa))**2
*---------------------------------------------------------------------
 890  print*,' [1] Bechetti-Greenlees (not well determined)    '
      print*,'                      see: ADNDT 17 (1976) p6    ' 
      print*,'     Potential B2 of Nucl Phys A190 (1972) 437   '
      print*,' [2] Read Sao Paulo potential from file          '
      print*,'                      see: Reference             ' 
      print*,' [3] D.Y. Pang et al. GDP08                      '
      print*,'                      see: PRC 79 (2009) 024615  '
      print*,' [4] C.T. Liang et al.                           '
      print*,'                      see: JPG 36 (2009) 085104  '
      print*,'----------------------------------------------------'
      ndim=1500
*---------------------------------------------------------------------
      ihesp=0
      read*,inopt
      if(inopt.lt.1.or.inopt.gt.4) go to 890
      write(18,*) inopt
      print*,' >>>> ',inopt
      if(inopt.eq.1) then
       a=a1
       ed=energy
       e=z1
       an=(a-2.d0*e)/a
       a13=a**0.3333333333d0
       vpr=151.9d0-0.17d0*ed+50.d0*an
       rr=1.20d0
       ar=0.72d0
       rpi=1.40d0
       api=0.88d0
       wvpi=41.7d0-0.33d0*ed+44.d0*an
       if(wvpi.lt.0.d0) wvpi=0.d0
       wpi=0.d0
       print 10,a
       print 101,e,ed
   10  format(1h ,' bechetti greenlees 3he potentials for a = ',f5.1)
  101  format(1h ,' z = ',f4.1,' at ',f6.2,' MeV 3He  energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 111
  111  format(1h ,'    vr     ro     ao     ws     ri     ai     wv ')
       print 12,vpr,rr,ar,wpi,rpi,api,wvpi
   12  format(1h ,7f7.3)
       vd=vpr
       rrd=rr
       ard=ar
       wd=wvpi
       rid=rpi
       aid=api
       wisd=wpi
       risid=rpi
       aisid=api
       vsod=2.5d0
       rsord=1.20d0
       asord=0.72d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       print 112
  112  format(1h ,'    vso    rso    aso ')
       print 12,vsod,rsord,asord
      endif
      if(inopt.eq.2) then
*---------------------------------------------------------------------
*     the externally read Sao Paulo potential option requested
*---------------------------------------------------------------------
       ihesp=1
       print 29,a1
   29  format(1h ,' Sao Paulo 3He potential for a = ',f5.1)
       rcd=1.25d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print*,' printout is in',nrmax,' steps of',real(step)
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
       vsod=0.d0
       rsord=1.d0
       asord=1.d0
       wsod=0.d0
       rsoid=1.d0
       asoid=1.d0
       form='(5e14.7)'
*---------------------------------------------------------------------
       print*,' filename with Sao Paulo potential ?  '
       read '(a)',spfile
       open(22,file=spfile,status='unknown')
       print '(a,a)','  >>>> ',spfile
       write(18,'(a)') spfile
*---------------------------------------------------------------------
*      read things from the external file on channel 22
*      format is to be provided here (JAT 14/3/2008)
*---------------------------------------------------------------------
       read(22,'(a)') guff
       read(22,*) npts,spstep,spstart
       do ii=1,npts
	rgrid(ii)=spstart+(ii-1)*spstep
        read(22,*) funr(ii),funi(ii)
*       write(31,*)rgrid(ii),funr(ii),funi(ii)
       enddo
       rmax=rgrid(npts)
*---------------------------------------------------------------------
*      interpolate the potential to the internal working grid
*      c/i/so (nrmax points with step) first radial point at step
*---------------------------------------------------------------------
       do ii=1,nrmax
        r=ii*step
	if(r.lt.rmax) then
         pspr(ii)=terp(r,funr,rgrid,npts,ndim)
         pspi(ii)=terp(r,funi,rgrid,npts,ndim)
	else
    	 pspr(ii)=0.d0
	 pspi(ii)=0.d0
	endif
*       diagnostic prints
*       write(30,*) r,pspr(ii),pspi(ii)
       enddo
      endif
      if(inopt.eq.3) then
*---------------------------------------------------------------------
*     GDP08 potential for 3He  dypang
*---------------------------------------------------------------------
       a=a1
       ed=energy
       e=z1
       ap=3.0d0
       zp=2.0d0
       call gdp08(ap,zp,a1,z1,energy)
       print 710,a
       print 701,e,ed
  710  format(1h ,' DYPang et al 3He potentials for a = ',f5.1)
  701  format(1h ,' z = ',f4.1,' at ',f5.1,' MeV helion energy ')
       rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 711
  711  format(1h ,'    vr     ro     ao     wv     riv    aiv ')
       print 713,vd,rrd,ard,wd,rid,aid
       print 715
  715  format(1h ,'    ws     ris    ais    ')
       print 713,wisd,risid,aisid
  713  format(1h ,7f7.3)
       print 712
  712  format(1h ,'    vso    rso    aso ')
       print 713,vsod,rsord,asord
       print 714
  714  format(1h ,'    wso   wrso   waso ')
       print 713,wsod,rsoid,asoid
      endif
      if(inopt.eq.4) then
*---------------------------------------------------------------------
*     3He systematic potential of Liang Chun-Tian et al., JPG-2009
*---------------------------------------------------------------------
       ihesp=1
       do i=1,900
        pspr(i)=0.d0
        pspi(i)=0.d0
       enddo
       a=a1
       a13=a**0.3333333333d0
       ed=energy
       e=z1
       ap=3.0d0
       zp=2.0d0
       call liangchuntian(ap,zp,a1,z1,energy)
       rrrd=rrd*a13
       rrid=rid*a13
       rrisid=risid*a13
       do ii=1,nrmax
        r=ii*step
        pspr(ii)=ws(r,vd,rrrd,ard)
        pspi(ii)=ws(r,wd,rrid,aid)+wsd(r,wisd,rrisid,aisid)
       enddo
       print 720,a
       print 701,e,ed
  720  format(1h ,' CTLiang et al. 3He potentials for a = ',f5.1)
*      rcd=1.30d0
       print*,' Coulomb radius parameter = ',real(rcd)
       print 711
       print 713,vd,rrd,ard,wd,rid,aid
       print 715
       print 713,wisd,risid,aisid
       print 712
       print 713,vsod,rsord,asord
       print 714
       print 713,wsod,rsoid,asoid
       vd=1.d0
       rrd=99.d0
       ard=1.d0
       wd=1.d0
       rid=99.d0
       aid=1.d0
       wisd=0.d0
       risid=rid
       aisid=1.d0
      endif
*---------------------------------------------------------------------
      return
      end
*---------------------------------------------------------------------
*  3He systematic pot of Liang Chun-Tian et al., JPG 36 (2009) 085104
*---------------------------------------------------------------------
      subroutine liangchuntian(ap,zp,at,zt,Einc)
      implicit real*8(a-h,o-z)
      real*8 AT, ZT, AP, ZP, at13, Einc
      real*8 vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      real*8 vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      at13 = AT**(1./3.)
      vd = 118.36-0.2071*Einc + 6.3961e-5*Einc*Einc
      vd = vd+26.001*(AT-2.*ZT)/AT+0.5668*ZT/at13
      rrd= 1.1657 + 0.0401/at13
      ard= 0.6641 + 0.0305*at13
      wd = -6.8871 + 0.3115*Einc -6.8096e-4*Einc*Einc
      rid= 1.4022+ 0.0418/at13
      aid= 0.7732+ 0.0219*at13
      wisd = 20.119-0.1626*Einc-5.4067*(AT-2.*ZT)/AT+1.2087*at13
      risid= 1.1802+ 0.0587/at13
      aisid= 0.6292+0.0657*at13
      vsod = 2.0491d0+9.9804e-3*at13
      rsord= 0.7211+0.0586/at13
      asord= 0.7643 + 0.0535*at13
      wsod = -1.1591d0
      rsoid= rsord
      asoid= asord
      rcd  = 1.289d0
      return
      end
*-------------------------------------------------------------
      subroutine gdp08(ap,zp,at,zt,Einc)
      implicit real*8(a-h,o-z)
      real*8 Wv, rrwv, aawv, Ws, rrws, aaws, V, rv, av
      real*8 VVso, rvso, avso, WWso, rwso, awso
      real*8 V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0
      real*8 WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT
      real*8 WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST
      real*8 VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP, RC0,RCA, RCAP, RC
      real*8 WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP
      real*8 AT, ZT, AP, ZP, Einc, EC, Epot, Ecm, varpsilon
      real*8 MDA,MDB,MDC,VMD, PI, AT13, zero
      common/pot1/vd,rrd,ard,wd,rid,aid,wisd,risid,aisid
      common/pot2/vsod,rsord,asord,wsod,rsoid,asoid,rcd
      common/iunit/iunit
      common/GDP08Para/
     &        V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0,
     &        WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT,
     &        WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST,
     &        VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP,
     &        WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP,
     &        RC0, RCA, RCAP
      PI=4.0d0*datan(1.0d0)
      call initGDP08
c     real part
      V0=118.25d0
      VE=-0.12512d0
      R0=1.3007
      R0A =-0.4816
      A0 =0.8148
      VC=1.0
c     Wv
      WV0=38.481
      RWV=1.3120
      RWVA=-0.1290
      AWV=0.8399
      WVE0=156.09
      WVEW=52.442
c     Ws
      WS0=35.037
      WST=34.181
      RWS=1.3120
      RWSA=-0.1290
      AWS=0.8399
      WSE0=30.755
      WSEW=106.36
c     Coulomb
      RC0=1.238
      RCA=0.116
c     isovector parameter
      if(ZP.EQ.2.) then
        varpsilon=(AT-2.0d0*ZT)/AT
      elseif(ZP.EQ.1.) then
        varpsilon=(2.0d0*ZT-AT)/AT
      else
        varpsilon=0.0d0
      endif
c     Coulomb correction, asumming alpha=1.0 as in CH89
      if(RCAP.EQ.0.0) RCAP=-1.0d0/3.0d0
      Ecm=AT*Einc/(AT+AP)
      RC=RC0+RCA*AT**RCAP
      if(VC .EQ. 1.0) then
        EC =(1.73d0/RC)*ZT*ZP/AT**(1.d0/3.d0)
      else
        EC=0.0d0
      endif
      Epot=Einc-EC
      AT13=AT**(1.0d0/3.0d0)
      zero=0.0d0
c     real part
      V =V0 + VE*Epot + (VT+VTE*Epot)*varpsilon
      rv=R0 + (R0A+R0AE*Epot)/AT13
      av=A0
c     volume imag
      if(WVE0 .LE. 0.0) then
        Wv =(WV0+WVE*Epot) + (WVT+WVTE*Epot)*varpsilon
      else
        Wv =(WV0+WVT*varpsilon)/(1.0d0+dexp(-(Epot-WVE0)/WVEW))
      endif
      rrwv= RWV + RWVA/AT13
      aawv= AWV + AWVT*varpsilon
c     surface imag
      if(WSE0 .LE. 0.0) then
        Ws=(WS0+WSE*Epot) + (WST+WSTE*Epot)*varpsilon
      else
        Ws=(WS0+WST*varpsilon)/(1.0d0+dexp( (Epot-WSE0)/WSEW))
      endif
      rrws= RWS + RWSA/AT13
      aaws= AWS + AWST*varpsilon
c     spin-orbit real
      VVso=VSO + VSOE*Epot + VSOA*AT**VSAP
      rvso=RSO + RSOA*AT13
      avso=ASO
c     spin-orbit imag
      WWso=WSO + WSOE*Epot + WSOA*AT**WSAP
      rwso=RWO + RWOA*AT13
      awso=AWO
c     Coulomb radius
      RC=RC0 +RCA*AT**RCAP
      vd  =V
      rrd =RV
      ard =AV
      wd  =Wv
      rid =rrwv
      aid =aawv
      wisd=Ws
      risid= rrws
      aisid= aaws
      vsod=VVso
      rsord= rvso
      asord= avso
      wsod=WWso
      rsoid= rwso
      asoid= awso
      return
      end
*-------------------------------------------------------------
      subroutine initGDP08
      implicit none
      real*8 V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0
      real*8 WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT
      real*8 WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST
      real*8 VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP
      real*8 WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP
      real*8 RC0,RCA, RCAP
      common/GDP08Para/
     &        V0, VE,  VT, VTE, R0,  R0A, R0AP, A0, VC, R0AE,MDAE,MDA0,
     &        WV0,WVE, WVT,WVTE,RWV, RWVA,RWVP, AWV,WVC,WVE0,WVEW,AWVT,
     &        WS0,WSE, WST,WSTE,RWS, RWSA,RWSP, AWS,WSC,WSE0,WSEW,AWST,
     &        VSO,VSOE,RSO,RSOA,RSOP,ASO, VSOA, VSAP,
     &        WSO,WSOE,RWO,RWOA,RWOP,AWO, WSOA, WSAP,
     &        RC0, RCA, RCAP
c     ------------------
      V0  = 0.0d0
      VE  = 0.0d0
      VT  = 0.0d0
      VTE = 0.0d0
      R0  = 0.0d0
      R0A = 0.0d0
      R0Ap= 0.0d0
      A0  = 0.0d0
      VC  = 1.0d0
      R0AE= 0.0d0
      MDAE= 0.0d0
      MDA0= 0.0d0
c     ------------------
      WV0 = 0.0d0
      WVE = 0.0d0
      WVT = 0.0d0
      WVTE= 0.0d0
      RWV = 0.0d0
      RWVA= 0.0d0
      RWVP= 0.0d0
      AWV = 0.0d0
      WVC = 1.0d0
      WVE0= 0.0d0
      WVEW= 0.0d0
      AWVT= 0.0d0
c     ------------------
      WS0 = 0.0d0
      WSE = 0.0d0
      WST = 0.0d0
      WSTE= 0.0d0
      RWS = 0.0d0
      RWSA= 0.0d0
      RWSP= 0.0d0
      AWS = 0.0d0
      WSC = 1.0d0
      WSE0= 0.0d0
      WSEW= 0.0d0
      AWST= 0.0d0
c     ------------------
      VSO = 0.0d0
      VSOE= 0.0d0
      RSO = 1.2d0
      RSOA= 0.0d0
      RSOP= 0.0d0
      ASO = 0.6d0
      VSOA= 0.0d0
      VSAp= 0.0d0
c     ------------------
      WSO = 0.0d0
      WSOE= 0.0d0
      RWO = 1.2d0
      RWOA= 0.0d0
      RWOP= 0.0d0
      AWO = 0.6d0
      WSOA= 0.0d0
      WSAp= 0.0d0
c     ------------------
      RC0 = 0.0d0
      RCA = 0.0d0
      RCAP= 0.0d0
      return
      end
*---------------------------------------------------------------------
      subroutine jlm(potr,poti,nucleon,aa,zz,Elab,step,nrmax)
*     ------------------------------------------------------------
*     JLM local density approximation (LDA) N+A potentials          
*     JAT February 1999 - tested March 1999.
*     ------------------------------------------------------------
      implicit real*8(a-h,o-z)
      common/density/rad(500),rhon(500),rhop(500),rhom(500)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4),nn
      real*8 xrt(200),wrt(200),xmu(200),wmu(200)
      real*8 potr(900),poti(900)
*     real*8 gbden(401)
      character dname*12,guff*4
*     character fmt*20
      character nucleon,ans,nucl*12
      common/arrays/a,b,c,d,f
      common /dcom/ rho0, rad0, diff, beta2, b2bar
      common/djlm/rlr,rli
      pi=4.d0*datan(1.d0)
      do i=1,900
       potr(i)=0.d0
       poti(i)=0.d0
      enddo
      print*,' ---------------------------------------------------'
      print*,' JLM local density approximation nucleon potentials'
      print*
*     ------------------------------------------------------------
*     print*,' Output file trailer for potentials:'
*     read '(a)',fname
*     open(17,file='jlm.'//fname,status='unknown')
*     open(22,file='density',status='unknown')
*     open(27,file='jlmf.'//fname,status='unknown')
*     print '(a,a)','  Output is to file: ','jlm.'//fname
*     ------------------------------------------------------------
*  16 print*,' neutron (n) or proton (p) optical potential'  
*     read '(a)',nucleon
*     if(nucleon.ne.'n'.and.nucleon.ne.'p') goto 16
      nuc=-1
      nucl='proton'
      if(nucleon.eq.'n') then
       nuc=1
       nucl='neutron'
      endif
      print '(a,2x,a)','  Optical potential for ',nucl
*     ------------------------------------------------------------
*     print*,' Target mass A and charge Z '
*     read*,aa,zz
*     print9,aa,zz
      nn=aa-zz
      aa13=aa**(1.d0/3.d0)
      if(nuc.eq.-1) rc=1.123d0*aa13+2.35d0/aa13-2.07d0/aa
*     print*,' Lab Energy (MeV) '
*     read*,E
*     read*,Elab
*     print9,Elab
*     compute cm energy
      E=Elab*aa/(aa+1.d0) 
*     print'(a,f10.4)',' incident cm energy (MeV) =  ',E
*     ------------------------------------------------------------
*     print*,' potential at radii: rmin,rmax,rstep '
*     read*,rmin,rmax,rstep
*     ------------------------------------------------------------
      rmin=0.d0
      rmax=nrmax*step
      rstep=step
   14 print*,' Uses the JLM parameterisation of Bauge '
*     ------------------------------------------------------------
*     print9,rmin,rmax,rstep
    9 format(3f10.4)
*  14 print*,' Choose JLM potential parameter set '
*     print*,'   1) High energy set (>15 MeV) '
*     print*,'   2) Low energy set  (<15 MeV) '
*     read*,iset
*     print*,iset
*     if(iset.ne.1.and.iset.ne.2) goto 14
      iset=1
*     if(Elab.lt.15.d0) iset=2 
*  11 print*,' Version of local density approximation '
*     print*,'   1) LDA (rx=rt)'
*     print*,'   2) LDA (rx=rp)'
*     print*,'   3) LDA (Mid-point) ******** '
*     read*,ilda
*     print*,ilda
*     if(ilda.lt.1.or.ilda.gt.3) goto 11
      print*,' Uses the mid-point LDA prescription '
      ilda=3
   17 print*,' Choose assumed target density     '
      print*,'   1) Negele (Fermi) form          '
      print*,'   2) Specify rms radius           '
      print*,'   3) Oscillator form              '
      print*,'   4) 3pF three parameter Fermi    '
      print*,'   5) read Alex Brown HF densities '
*     print*,'   6) read other density           '
      read*,icho
      if(icho.lt.1.or.icho.gt.5) go to 17
      print*,' >>>> ',icho
      write(18,*) icho
      if(abs(icho).gt.5) goto 17
      if(icho.eq.2) then
        print*,' ------------------------------------------------------'
        print*,' input rms radius (fm) '
        read*,rms
        print*,' >>>> ',rms
	write(18,*) rms
        print*,' ------------------------------------------------------'
   18   print*,'   1) Gaussian density '
        print*,'   2) Woods-Saxon density'
        read*,irms
        if(irms.lt.1.or.irms.gt.2) go to 18
        print*,' >>>> ',irms
	write(18,*)irms
        if(irms.ne.1.and.irms.ne.2) goto 18
        if(irms.eq.2) then
         diff=0.54d0
         print*,'  default diffuseness (0.54 fm) (y/n) '
         read '(a)',ans
         print '(a,a)',' >>>> ',ans
         write(18,'(a)') ans
         if(ans.eq.'n') then
          print*,'  enter diffuseness '
          read*,diff
          print*,' >>>> ',diff
	  write(18,*)diff
         endif
         call finder(aa,rms,rho0,rad0,diff)
         iden=3
         print*,' -----------------------------------------------------'
         print*,' iden = 3: WooSax:  rho0 = ',rho0
         print*,' rad0 = ',rad0,' diff = ',diff
        else
         iden=2
         gamma=sqrt(2.d0/3.d0)*rms
         rho0=aa/(sqrt(pi)*gamma)**3
         print*,' -----------------------------------------------------'
         print*,' iden = 2: Gauss :  gamma = ',gamma
         print*,' rho0 = ',rho0 
        endif       
       else if(icho.eq.1) then
*      Negele Woods-Saxon parameters
        iden=1
        diff=0.54d0
        rad0=(0.978d0+0.0206d0*aa13)*aa13
        rho0=3.d0*aa/(4.d0*pi*rad0**3*(1.d0+(pi*diff/rad0)**2))
	vol=aa
	call volrms(vol,rms)
        print*,' -----------------------------------------------------'
        print*,' iden = 1: Negele:  rho0 = ',rho0
        print*,' rad0 = ',rad0,' diff = ',diff
	print*,' rms mass radius =',rms
      else if(icho.eq.3) then
        print*,' Input a and alfa values in '
        print*,' (1+alfa*[r/a]**2)exp(-r**2/a**2) '
        read*,rad0,diff
        print*,' >>>> ',rad0,diff
        write(18,*) rad0,diff
        iden=4
        call finder2(aa,rms,rho0,rad0,diff)
        print*,' -----------------------------------------------------'
        print*,' iden = 4: Oscill:  rho0 = ',rho0
        print*,' a    = ',rad0,' alfa = ',diff
        print*,' rho0 = ',rho0,' rms  = ',rms 
      else if(icho.eq.4) then
        print*,' Input rad, diff and w values in '
        print*,' (1+w*[r/rad]**2)/(1+exp(r-rad)/diff)  '
        read*,rad0,diff,www
        print*,' >>>> ',rad0,diff,www
        write(18,*) rad0,diff,www
        iden=5
        call finder3(aa,rms,rho0,rad0,diff,www)
        print*,' -----------------------------------------------------'
        print*,' iden = 5: 3pFermi:  rho0 = ',rho0
        print*,' rad  = ',rad0,' diff = ',diff
        print*,'  w =   ',www
        print*,' rho0 = ',rho0,' rms  = ',rms 
*     ------------------------------------------------------------
*       print to file for reading by smat 
*       drin=0.05d0
*       nval=401
*       istart=0
*       fmt='(4d19.8)'
*       write(22,10) fmt,drin,nval,istart
*  10   format(a20,f10.4,2i10)
*       do iir=1,nval
*        rr=(iir-1)*drin
*        gbden(iir)=oscden(rr,rho0,rad0,diff)
*        write(19,*) rr,gbden(iir)
*       enddo
*       write(22,'(4(d19.8))') (gbden(iir),iir=1,nval)
*     ------------------------------------------------------------
      else if(icho.eq.5.or.icho.eq.6) then
*       density is read from file
        iden=6
        ijax=500
        print*,' using read matter density '
        print*,' -----------------------------------------------------'
        print*,' iden = 6: read matter density  '
        if(icho.eq.5) print*,' f(r)=Hartree-Fock density'
        if(icho.eq.6) print*,' f(r)=external density'
        print*
        print*,' file with required density is'
        read '(a)',dname
        print '(a)','  >>>> '//dname
        write(18,'(a)') dname
        open(19,file=dname,status='unknown')
        print*
*       following for reading of Hartree-Fock densities
        if(icho.eq.5)then
         rewind 19
         ico=0
         read(19,'(a)') guff
         read(19,'(a)') guff
         read(19,'(a)') guff
         do ijj=1,1000
          read(19,*,end=707) rad(ijj),rhop(ijj),rhon(ijj),rhom(ijj)
          ico=ico+1
         enddo
  707    continue
         close(19)
         ival=ico
         drin=rad(2)-rad(1)
         print*,' number of read HF radii ',ival
         print*,' with step ',drin
*     ------------------------------------------------------------
*        do ijj=1,ival
*         write(30,*) rad(ijj),rhom(ijj)
*        enddo
*     ------------------------------------------------------------
        endif
        if(icho.eq.6)then
         rewind 19
         ico=0
         do ijj=1,1000
          read(19,*,end=708) rad(ijj),rhom(ijj)
          ico=ico+1
         enddo
  708    continue
         close(19)
         ival=ico
         drin=rad(2)-rad(1)
         print*,' number of read radii ',ival
         print*,' with step ',drin
*     ------------------------------------------------------------
*        do ijj=1,ival
*         write(31,*) rad(ijj),rhom(ijj)
*        enddo
*     ------------------------------------------------------------
        endif
      endif
      print*,' -----------------------------------------------------'
*     ------------------------------------------------------------
*     look up potential parameters
      call assign(iset)
      alpha=nuc*(nn-zz)/aa
*     ------------------------------------------------------------
*     print*,' Maximum NN relative separation (rptmax) for folding'
*     read*,rptmax
*     ------------------------------------------------------------
      rptmax=4.d0
*     ------------------------------------------------------------
*     print9,rptmax
*     print*,' Quadrature points: radial and cos(theta) integrals'
*     read*,mqr,mqmu
*     ------------------------------------------------------------
      mqr=32
      mqmu=48
*     ------------------------------------------------------------
*     print*,mqr,mqmu
      zero=0.d0
      one=1.d0
*     for integration over radius of target density 
      call gauss(zero,rptmax,mqr,xrt,wrt)
*     for integration over cos(theta) 
      call gauss(-one,one,mqmu,xmu,wmu)
*     ------------------------------------------------------------
*     Gaussian effective interaction strength (for t=1.0 fm)
*     and for unit integrated strength
*     print*,' real and imaginary Gaussian folding ranges t'
*     read*,tr,ti
*     print9,tr,ti
*     ------------------------------------------------------------
      tr=1.d0
      ti=1.d0
*     ------------------------------------------------------------
      t2r=tr*tr
      t2i=ti*ti
      rmagr=1.d0/(sqrt(pi)*tr)**3
      rmagi=1.d0/(sqrt(pi)*ti)**3
*     ------------------------------------------------------------
      print*,' real and imaginary potential scalings lambda '
      if(rlr.lt.1.d-3.and.rli.lt.1.d-3) then
       print*,' systematics usually suggest  1.0  0.8 '
       read*,rlr,rli
       print*,' >>>> ',rlr,rli
       write(18,*) rlr,rli
       print*
      else
       print*,' using  ',rlr,rli
       print*
      endif 
*     ------------------------------------------------------------
      i=0
      do rp=rstep,rmax+0.01d0,rstep
      i=i+1
      potr(i)=0.d0
      poti(i)=0.d0
      rp2=rp*rp
      do imqr=1,mqr
       rpt=xrt(imqr)
       rpt2=rpt*rpt
       gg3r=g3(rmagr,t2r,rpt2)*rpt2
       gg3i=g3(rmagi,t2i,rpt2)*rpt2
       summr=0.d0
       summi=0.d0
        do imqmu=1,mqmu      
        bval=rp*rpt*xmu(imqmu)
        rt2=rpt2+rp2-2.d0*bval
        rt=sqrt(rt2)     
        if(ilda.eq.1) then
         rx=rt
        elseif(ilda.eq.2) then
         rx=rp
        elseif(ilda.eq.3) then
         rx=sqrt(rpt2/4.d0+rp2-bval)
        endif
       if(iden.eq.1) then
        rhot=woosax(rt,rho0,rad0,diff)
        rhox=woosax(rx,rho0,rad0,diff)
       elseif(iden.eq.2) then
        rhot=gauden(rt,rho0,gamma)
        rhox=gauden(rx,rho0,gamma)
       elseif(iden.eq.3) then
        rhot=woosax(rt,rho0,rad0,diff)
        rhox=woosax(rx,rho0,rad0,diff)
       elseif(iden.eq.4) then
        rhot=oscden(rt,rho0,rad0,diff)
        rhox=oscden(rx,rho0,rad0,diff)
       elseif(iden.eq.5) then
        rhot=fermi3(rt,rho0,rad0,diff,www)
        rhox=fermi3(rx,rho0,rad0,diff,www)
       elseif(iden.eq.6) then
        rhot=1.d-15
        if(rt.lt.rad(ival)) rhot=terp(rt,rhom,rad,ival,ijax)
        rhox=1.d-15
        if(rx.lt.rad(ival)) then
           rhoxn=terp(rx,rhon,rad,ival,ijax)
           rhoxp=terp(rx,rhop,rad,ival,ijax)
           rhox=rhoxn+rhoxp
           alpha=nuc*(rhoxn-rhoxp)/rhox
*          alpha=nuc*(nn-zz)/aa
        endif
       endif
*       if proton take care of Coulomb interaction
        Ecal=E
        if(nuc.eq.-1) then
         if(rx.ge.rc) then
          vc=1.44d0*zz/rx
         else
          vc=0.72d0*zz/rc*(3.d0-(rx/rc)**2) 
         endif
         Ecal=E-vc
        endif
        call pots(Ecal,rhox,iset,V0,W0,V1,W1,rmtilde)
        summr=summr+wmu(imqmu)*rhot*(V0+alpha*V1)/rhox
        summi=summi+wmu(imqmu)*rhot*rmtilde*(W0+alpha*W1)/rhox
       enddo
       potr(i)=potr(i)+wrt(imqr)*summr*gg3r
       poti(i)=poti(i)+wrt(imqr)*summi*gg3i
      enddo
      potr(i)=2.d0*pi*potr(i)*rlr
      poti(i)=2.d0*pi*poti(i)*rli
*     print*,i
*     write(24,*) rp,potr(i),poti(i)
      enddo     
*     ------------------------------------------------------------
*     call print1(17,aa,zz,nuc,Elab)
*     call print2(17)
*     write(17,'(5(e14.7))') (potr(i),i=2,301) 
*     call print2(17)
*     write(17,'(5(e14.7))') (poti(i),i=2,301)
*     ------------------------------------------------------------
*     for cupid
*     ------------------------------------------------------------
*     write(18,'(a)')' 31.       1.0   1.20    0.61    0.       0.0  '
*     write(18,'(a)')'150.0   0.0'
*     write(18,'(5(e16.7))') (potr(i),i=2,151) 
*     write(18,'(a)')' 31.       0.0   1.20    0.61    0.       1.0  '
*     write(18,'(a)')'150.0   1.0'
*     write(18,'(5(e16.7))') (poti(i),i=2,151) 
*     write(18,'(a)') ' 0'
*     ------------------------------------------------------------
*     for fresco
*     ------------------------------------------------------------
*     write(27,'(a)')'201   0.1    0.0'
*     write(27,'(5(e16.7))') (potr(i),i=1,201) 
*     write(27,'(a)')'201   0.1    0.0'
*     write(27,'(5(e16.7))') (poti(i),i=1,201) 
*     ------------------------------------------------------------
      return
      end
*     ------------------------------------------------------------
      subroutine pots(E,rho,iset,V0,W0,V1,W1,rmtilde)
      implicit real*8(a-h,o-z)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4),ImN
      common/arrays/a,b,c,d,f
      ef=fermi(E,rho,iset)  
      V0=0.d0
      ReN=0.d0
      rmtilde=0.d0
      do i=1,3
       do j=1,3
        con=rho**i*E**(j-1)
        V0=V0+a(i,j)*con
        ReN=ReN+b(i,j)*con
        rmtilde=rmtilde+c(i,j)*con
       enddo
      enddo
      rmtilde=1.d0-rmtilde
      V0P=0.d0
      do i=1,3
       do j=2,3
        con=(j-1)*rho**i*E**(j-2)
        V0P=V0P+a(i,j)*con
       enddo
      enddo
      rmstar=1.d0-V0P  
      W0=0.d0
      ImN=0.d0
      do i=1,4
       do j=1,4
        con=rho**i*E**(j-1)
        W0=W0+d(i,j)*con
        ImN=ImN+f(i,j)*con
       enddo
      enddo
c----------------------------------------------------------
*     dd=100.d0
*     if(iset.eq.1) dd=600.d0
*     dd=625.d0
*     changed according to Pang/Bauge 2011
c----------------------------------------------------------
      dd=126.25d0
      W0=W0/(1.d0+dd/(E-ef)**2)
      ImN=ImN/(1.d0+1.d0/(E-ef))
      rmbar=rmstar/rmtilde
      V1=rmtilde*ReN
      W1=ImN/rmbar
      return
      end
*     ------------------------------------------------------------
      subroutine assign(iset)
      implicit real*8(a-h,o-z)
      real*8 a(3,3),b(3,3),c(3,3),d(4,4),f(4,4)
      common/arrays/a,b,c,d,f
c     ------------------
      a(1,1)=-0.9740d+3
      a(1,2)= 0.1126d+2
      a(1,3)=-0.4250d-1
      a(2,1)= 0.7097d+4
      a(2,2)=-0.1257d+3
      a(2,3)= 0.5853d+0
      a(3,1)=-0.1953d+5
      a(3,2)= 0.4180d+3
      a(3,3)=-0.2054d+1
c     ------------------     
      b(1,1)= 0.3601d+3
      b(1,2)=-0.5224d+1
      b(1,3)= 0.2051d-1
      b(2,1)=-0.2691d+4
      b(2,2)= 0.5130d+2
      b(2,3)=-0.2470d+0 
      b(3,1)= 0.7733d+4
      b(3,2)=-0.1717d+3
      b(3,3)= 0.8846d+0
c     ------------------
      c(1,1)= 0.4557d+1
      c(1,2)=-0.5291d-2
      c(1,3)= 0.6108d-5
      c(2,1)=-0.2051d+1
      c(2,2)=-0.4906d+0
      c(2,3)= 0.1812d-2
      c(3,1)=-0.6509d+2
      c(3,2)= 0.3095d+1   
      c(3,3)=-0.1190d-1  
      go to 100
c     ------------------
*     if(iset.eq.1) then
      d(1,1)=-0.1483d+4
      d(1,2)= 0.3718d+2
      d(1,3)=-0.3549d+0
      d(1,4)= 0.1119d-2
      d(2,1)= 0.2988d+5
      d(2,2)=-0.9318d+3
      d(2,3)= 0.9591d+1
      d(2,4)=-0.3160d-1
      d(3,1)=-0.2128d+6
      d(3,2)= 0.7209d+4
      d(3,3)=-0.7752d+2
      d(3,4)= 0.2611d+0
      d(4,1)= 0.5125d+6
      d(4,2)=-0.1796d+5
      d(4,3)= 0.1980d+3
      d(4,4)=-0.6753d+0  
c     ------------------            
      f(1,1)= 0.5461d+3
      f(1,2)=-0.1120d+2
      f(1,3)= 0.1065d+0
      f(1,4)=-0.3541d-3
      f(2,1)=-0.8471d+4
      f(2,2)= 0.2300d+3
      f(2,3)=-0.2439d+1
      f(2,4)= 0.8544d-2
      f(3,1)= 0.5172d+5
      f(3,2)=-0.1520d+4
      f(3,3)= 0.1717d+2
      f(3,4)=-0.6211d-1
      f(4,1)=-0.1140d+6
      f(4,2)= 0.3543d+4
      f(4,3)=-0.4169d+2
      f(4,4)= 0.1537d+0
*     else
*     These from JLM
      d(1,1)=-0.5138d+3
      d(1,2)=-0.2985d+2
      d(1,3)= 0.1452d+1
      d(1,4)= 0.9428d-1
      d(2,1)= 0.9078d+4
      d(2,2)= 0.5757d+3
      d(2,3)=-0.3435d+2
      d(2,4)=-0.2310d+1
      d(3,1)=-0.6192d+5
      d(3,2)=-0.4155d+4
      d(3,3)= 0.2657d+3
      d(3,4)= 0.1882d+2
      d(4,1)= 0.1516d+6
      d(4,2)= 0.1037d+5
      d(4,3)=-0.6748d+3
      d(4,4)=-0.5014d+2 
c     ------------------            
      f(1,1)= 0.6597d+3
      f(1,2)= 0.4509d+1
      f(1,3)=-0.2383d+1
      f(1,4)=-0.4324d-1
      f(2,1)=-0.1263d+5
      f(2,2)=-0.6572d+2
      f(2,3)= 0.5866d+2
      f(2,4)= 0.1348d+1
      f(3,1)= 0.9428d+5
      f(3,2)= 0.5972d+3
      f(3,3)=-0.4923d+3
      f(3,4)=-0.1295d+2
      f(4,1)=-0.2453d+6
      f(4,2)=-0.1800d+4
      f(4,3)= 0.1358d+4
      f(4,4)= 0.3836d+2
*     These from Bauge (old, Bauge paper)
      d(1,1)=-0.6599d+3
      d(1,2)= 0.1077d+2
      d(1,3)=-0.7886d-1
      d(1,4)= 0.1875d-3
      d(2,1)= 0.1144d+5
      d(2,2)=-0.2908d+3
      d(2,3)= 0.2443d+1
      d(2,4)=-0.6203d-2
      d(3,1)=-0.7451d+5
      d(3,2)= 0.2207d+4
      d(3,3)=-0.1993d+2
      d(3,4)= 0.5175d-1
      d(4,1)= 0.1761d+6
      d(4,2)=-0.5458d+4
      d(4,3)= 0.5113d+2
      d(4,4)=-0.1339d+0   
c     ------------------           
      f(1,1)= 0.4596d+3
      f(1,2)=-0.6440d+1
      f(1,3)= 0.4040d-1
      f(1,4)=-0.9009d-4
      f(2,1)=-0.7693d+4
      f(2,2)= 0.1464d+3
      f(2,3)=-0.1025d+1
      f(2,4)= 0.2337d-2
      f(3,1)= 0.5525d+5
      f(3,2)=-0.1112d+4
      f(3,3)= 0.7967d+1
      f(3,4)=-0.1802d-1
      f(4,1)=-0.1437d+6
      f(4,2)= 0.3038d+4
      f(4,3)=-0.2220d+2
      f(4,4)= 0.5026d-1
  100 continue
c     ------------------
c     These from Bauge (from MOM code, Pang 2011)
      d(1,1)=-0.65986d+3
      d(1,2)= 0.10768d+2
      d(1,3)=-0.78863d-1
      d(1,4)= 0.18755d-3
      d(2,1)= 0.11437d+5
      d(2,2)=-0.29076d+3
      d(2,3)= 0.24430d+1
      d(2,4)=-0.62028d-2
      d(3,1)=-0.74505d+5
      d(3,2)= 0.22068d+4
      d(3,3)=-0.19926d+2
      d(3,4)= 0.51754d-1
      d(4,1)= 0.17609d+6
      d(4,2)=-0.54579d+4
      d(4,3)= 0.51127d+2
      d(4,4)=-0.13386d+0
c     ------------------
      f(1,1)= 0.45959d+3
      f(1,2)=-0.64399d+1
      f(1,3)= 0.40403d-1
      f(1,4)=-0.90086d-4
      f(2,1)=-0.76929d+4
      f(2,2)= 0.14639d+3
      f(2,3)=-0.10244d+1
      f(2,4)= 0.23367d-2
      f(3,1)= 0.55250d+5
      f(3,2)=-0.11121d+4
      f(3,3)= 0.79667d+1
      f(3,4)=-0.18008d-1
      f(4,1)=-0.14373d+6
      f(4,2)= 0.30382d+4
      f(4,3)=-0.22202d+2
      f(4,4)= 0.50258d-1
c     ------------------
*     endif
      return  
      end       
c--------------------------------------------------------------- 
      real*8 function fermi(E,rho,iset)
      implicit real*8(a-h,o-z)
*     Bauge modification PRC 58 page 1120
*     E0=10.d0
*     change accirding to Pang/Bauge 2011
      E0=9.d0
      ae=2.d0
      fermih=rho*(-510.8d0+3222.d0*rho-6250.d0*rho*rho)
      fermil=-22.d0-rho*(298.52d0-3760.23d0*rho+
     +       12345.82d0*rho*rho)
      fwt=1.d0/(1.d0+exp((E-E0)/ae))
      fermi=fwt*fermil+(1.d0-fwt)*fermih
      return
      end
c--------------------------------------------------------------- 
      real*8 function woosax(r,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      con=(r-rad0)/diff
      woosax=rho0/(1.d0+exp(con))
      return
      end
c--------------------------------------------------------------- 
      real*8 function fermi3(r,rho0,rad0,diff,www)
      implicit real*8(a-h,o-z)
      con=(r-rad0)/diff
      con1=(r/rad0)**2
      fermi3=rho0*(1+www*con1)/(1.d0+exp(con))
      return
      end
c--------------------------------------------------------------- 
      subroutine finder3(aa,rms,rho0,rad0,diff,www)
      implicit real*8(a-h,o-z)
      rho0=1.d0
      call vol3pf(vol,rms,rho0,rad0,diff,www)
      rho0=aa/vol
      call vol3pf(vol,rms,rho0,rad0,diff,www)
      print 24, vol, rms
 24   format('  Volume: ',f12.7,' rms ',f12.7)
      return
      end 
c--------------------------------------------------------------- 
      real*8 function oscden(r,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      con=(r/rad0)**2
      oscden=rho0*(1+diff*con)*exp(-con)
      return
      end
c--------------------------------------------------------------- 
      real*8 function gauden(r,rho0,gamma)
      implicit real*8(a-h,o-z)
      gauden=rho0*exp(-r*r/gamma/gamma)
      return
      end
c--------------------------------------------------------------- 
      real*8 function g3(rmag,t2,r2)
      implicit real*8(a-h,o-z)
      g3=rmag*exp(-r2/t2)
      return
      end
c--------------------------------------------------------------- 
      subroutine print1(icc,aa,zz,nuc,E)
      implicit real*8(a-h,o-z)
      zn=1.0
      if(nuc.eq.1) zn=0.0
      write(icc,'(f4.1,a,2f5.1)') zn,'  1.0',zz,aa
      write(icc,'(f5.1)') E
      write(icc,1) '0.005  16000                                      '
      write(icc,1) '0.  30.  0.1                                      '   
      write(icc,1) '1.25                                              '
    1 format(a,a,a)
      return
      end
c--------------------------------------------------------------- 
      subroutine print2(icc)
      implicit real*8(a-h,o-z)
      write(icc,1) 'read                                              '  
      write(icc,1) '  300       0.100    ((5e14.7))                   '
    1 format(a,a,a)
      return
      end
c--------------------------------------------------------------- 
      subroutine finder(aa,rms,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      common /dcom/ v0, r0, a0, beta2, b2bar
      pi=4.d0*datan(1.d0) 
c     Woods-Saxon for quadrupole deformed nucleus (beta4 = 0)
      a0=diff
      beta2=0.d0
      rwant=rms
      v0 = 0.13269d0
      r0 = 1.1d0*aa**(1.d0/3.d0) 
      b2bar = r0 * sqrt(5.d0/(4.d0*pi)) * beta2
      call volrms(vol,rms)
      step = 0.05d0
      r0 = r0 + step
      err = rwant -rms
      icount =0
 10   call volrms(vol,rms2)
      err2 = rwant -rms2 
      if (abs(err2).gt.1.d-5) then
       icount = icount+1
       write(6,*) 'radius',icount, rms2
       rms =  rms2
       grad = (err2-err)/step
       err = err2
       step = - err/grad
       r0 = r0 + step
       go to 10
      end if
      call volrms(vol,rms)
      step = 0.05d0
      v0 = v0 + step
      err = aa - vol
      icount =0
 20   call volrms(vol2,rms)
      err2 = aa - vol2
      if (abs(err2).gt.1.d-8) then
       icount = icount+1
       write(6,*) 'depth',icount, vol2
       vol = vol2 
       grad = (err2-err)/step
       err = err2
       step = - err/grad
       v0 = v0 + step
       go to 20
      end if
      call volrms(vol,rms)
      print*
      print 24, vol, rms
 24   format('Volume: ',f12.7,' rms ',f12.7)
      print 25, v0, r0, a0, beta2
 25   format('depth: ',f12.7,' radius ',f12.7,'\n',
     +       'diffuse: ',f12.7,' defm ',f12.7)
      print*
      rho0=v0
      rad0=r0
      return
      end 
c--------------------------------------------------------------- 
      real*8 function den0(r)
      implicit real*8(a-h,o-z)
      common /dcom/ v0, r0, a0, beta2, b2bar
      e0 = 1.d0/(1.d0+exp((r-r0)/a0))
      e1 = e0*(1.d0-e0)/a0
      e2 = e1*(1.d0 - 2.d0*e0)/a0
      e3 = e1*(1.d0-6.d0*e0+6.d0*e0*e0)/(a0*a0)
      den0 = v0*( e0 + b2bar*b2bar*e2/10.d0 +       
     +       b2bar*b2bar*b2bar*e3/105.d0 )
      return
      end 
c-------------------------------------------------------------
      subroutine volrms(vol,rms)
      implicit real*8(a-h,o-z)
      pi=4.d0*datan(1.d0)
      h = 0.1d0
      rmax = 20.d0
      nval = int(rmax/h) + 1
      if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
      sum = 0.d0 
      sumr2 = 0.d0 
      do i=0, nval
       r = i*h
       r2 = r*r
       v = den0(r)
       r2vs = r2*v*simpfac(i,nval)
       sum = sum + r2vs
       sumr2 = sumr2 + r2*r2vs
      end do  
      sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
      sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
      vol = sum
      rms = sqrt(sumr2/sum)
      return
      end 
c-----------------------------------------------------------
      real*8 function simpfac(i,nstep) 
c     convert degrees to radians
      implicit none 
      integer i,nstep
      real*8 sfact 
      if ((i.eq.0).or.(i.eq.nstep)) then
       sfact= 1.d0 
      else if (mod(i,3).eq.0) then
       sfact= 2.d0
      else
       sfact=3.d0
      end if
      simpfac =sfact
      return
      end 
c-----------------------------------------------------------
      subroutine gauss(a,b,npoint,xri,wri)
      implicit real*8(a-h,o-z)
      real*8 xg(200),wg(200),xri(200),wri(200)
      call setmgl(npoint,xg,wg)
      do j=1,npoint
       xri(j) = (a+b)/2.d0 + (b-a)/2.d0*xg(j)
       wri(j) = (b-a)/2.d0*wg(j)
      enddo
      return
      end
c-----------------------------------------------------------
      subroutine setmgl( n, points, weight )
      implicit real*8  ( a-h, o-z )
      real*8 points(200), weight(200)
      real*8 poin16(300)
      pi = 4.d0*atan(1.d0)
      if (n.gt.200)  write (1,50)
50    format(' setmlg call with too many points')
      m = ( n + 1 ) / 2
      e1 = n * ( n + 1 )
      do 1 i = 1, m
      t = ( 4*i - 1 ) * pi / ( 4*n + 2 )
      x0 = ( 1.d0 - ( 1.d0 - 1.d0/n ) / ( 8.d0*n*n ) ) * cos(t)
      pkm1 = 1.d0
      pk = x0
      do 3 k = 2, n
      t1 = x0 * pk
      pkp1 = t1 - pkm1 - ( t1-pkm1 )/k + t1
      pkm1 = pk
      pk = pkp1
3     continue
      den = 1.d0 - x0*x0
      d1 = n * ( pkm1 - x0*pk )
      dpn = d1 / den
      d2pn = ( 2.d0*x0*dpn - e1*pk ) / den
      d3pn = ( 4.d0*x0*d2pn + (2.d0-e1)*dpn ) / den
      d4pn = ( 6.d0*x0*d3pn + (6.d0-e1)*d2pn ) / den
      u = pk / dpn
      v = d2pn / dpn
      h = -u * ( 1.d0 + 0.5d0*u*(v+u*(v*v-u*d3pn/(3.d0*dpn))))
      p = pk + h*(dpn+0.5d0*h*(d2pn+h/3.d0*(d3pn+0.25d0*h*d4pn)))
      dp = dpn + h*(d2pn+0.5d0*h*(d3pn+h*d4pn/3.d0))
      h = h - p / dp
      poin16(i) = x0 + h
      fx = d1 - h*e1*(pk+0.5d0*h*(dpn+h/3.d0* 
     1     (d2pn+0.25d0*h*(d3pn+0.2d0*h*d4pn))))
      weight(i) = 2.d0 * ( 1.d0 - poin16(i)*poin16(i)) / (fx*fx)
1     continue
      if ( m + m .gt. n ) poin16(m) = 0.d0
      do 10 i = n/2 + 1, n
      poin16(i) = poin16( n + 1 - i )
      weight(i) = weight( n + 1 - i )
      poin16( n + 1 - i ) = -poin16( n + 1 - i )
10    continue
      do 30 i=1,n
 30   points(i)=poin16(i)
      return
      end
c--------------------------------------------------------------- 
      subroutine finder2(aa,rms,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      pi=4.d0*datan(1.d0) 
      rho0=1.d0
      call volosc(vol,rms,rho0,rad0,diff)
      rho0=aa/vol
      call volosc(vol,rms,rho0,rad0,diff)
      print 24, vol, rms
   24 format('  Volume: ',f12.7,' rms ',f12.7)
      return
      end 
c-------------------------------------------------------------
      subroutine volosc(vol,rms,rho0,rad0,diff)
      implicit real*8(a-h,o-z)
      pi=4.d0*datan(1.d0)
      h = 0.1d0
      rmax = 20.d0
      nval = int(rmax/h) + 1
      if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
      sum = 0.d0 
      sumr2 = 0.d0 
      do i=0, nval
       r = i*h
       r2 = r*r
       v=oscden(r,rho0,rad0,diff)
       r2vs = r2*v*simpfac(i,nval)
       sum = sum + r2vs
       sumr2 = sumr2 + r2*r2vs
      end do  
      sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
      sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
      vol = sum
      rms = sqrt(sumr2/sum)
      return
      end 
c-------------------------------------------------------------
      subroutine vol3pf(vol,rms,rho0,rad0,diff,www)
      implicit real*8(a-h,o-z)
      pi=4.d0*datan(1.d0)
      h = 0.1d0
      rmax = 20.d0
      nval = int(rmax/h) + 1
      if (mod(nval,3).ne.0) nval = nval+3 - mod(nval,3)
      sum = 0.d0 
      sumr2 = 0.d0 
      do i=0, nval
       r = i*h
       r2 = r*r
       v=fermi3(r,rho0,rad0,diff,www)
       r2vs = r2*v*simpfac(i,nval)
       sum = sum + r2vs
       sumr2 = sumr2 + r2*r2vs
      end do  
      sum = 4.d0 * pi * sum * (3.d0*h/8.d0)
      sumr2 = 4.d0 * pi * sumr2 * (3.d0*h/8.d0)
      vol = sum
      rms = sqrt(sumr2/sum)
      return
      end 
c-------------------------------------------------------------
      real*8 function terp(r,fun,rgrid,npts,ndim)
c------------------------------------------------------------------------------
c     this function calculates, by interpolation, the value of a real
c     function at an arbitrary point r, when the value of the function 
c     (stored in array fun) is known on a grid of points. The values of the
c     npts points at which the function is known is stored in array rgrid.
c     ndim is the externally defined dimensions of the arrays fun and rgrid
c     JAT routine of some vintage.
c------------------------------------------------------------------------------
      implicit double precision(a-h,o-z)
      real*8 fun(ndim),y1,y2,y3,y4,y5,y6
      double precision rgrid(ndim)
      do 30 k=1,npts
      nst=0
      if(rgrid(k).lt.r) goto 30
      nst=max0(k-3,1)
      goto 33
   30 continue
   33 if(nst.gt.npts-5) nst=npts-5
      x1=rgrid(nst+0)
      x2=rgrid(nst+1)
      x3=rgrid(nst+2)
      x4=rgrid(nst+3)
      x5=rgrid(nst+4)
      x6=rgrid(nst+5)
      y1=fun(nst+0)
      y2=fun(nst+1)
      y3=fun(nst+2)
      y4=fun(nst+3)
      y5=fun(nst+4)
      y6=fun(nst+5)
      pii1=(x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)*(x1-x6)
      pii2=(x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)*(x2-x6)
      pii3=(x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)*(x3-x6)
      pii4=(x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)*(x4-x6)
      pii5=(x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)*(x5-x6)
      pii6=(x6-x1)*(x6-x2)*(x6-x3)*(x6-x4)*(x6-x5)
  777 xd1=r-x1
      xd2=r-x2
      xd3=r-x3
      xd4=r-x4
      xd5=r-x5
      xd6=r-x6
      pi1= xd2*xd3*xd4*xd5*xd6
      pi2= xd1*xd3*xd4*xd5*xd6
      pi3= xd1*xd2*xd4*xd5*xd6
      pi4= xd1*xd2*xd3*xd5*xd6
      pi5= xd1*xd2*xd3*xd4*xd6
      pi6= xd1*xd2*xd3*xd4*xd5
      terp=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4+
     + y5*pi5/pii5+y6*pi6/pii6
      return
      end
*------------------------------------------------------------------------
      subroutine av18wf
*------------------------------------------------------------------------
*     av18 wavefunction r-space - called once to set up function arrays
*     from routine folder if this is the wavefunction choice: ideutwf=4
*------------------------------------------------------------------------
      implicit real*8(a-h,o-z)
      parameter (rmn=939.565d0,rmp=938.272d0, rmu=rmn*rmp/(rmn+rmp)) 
      parameter (ed=-2.22452d0, h=0.01d0, hc=197.32696010352811d0)
      common/avstuff/uavs(1501),vavs(1501),uavd(1501),vavd(1501)
      common/avstuf2/uasp(1501),uadp(1501)
c     -------------------------------------------------------------------
      data vavs/1501*0.0d0/,vavd/1501*0.0d0/
c     s-state radial wave function: r=0 to 15 fm, steps of 0.01 fm
      data uavs/0.0d0,
     + 0.7920652E-03, 0.1588732E-02, 0.2394609E-02, 0.3214325E-02,
     + 0.4052529E-02, 0.4913901E-02, 0.5803156E-02, 0.6725043E-02,
     + 0.7684354E-02, 0.8685919E-02, 0.9734610E-02, 0.1083533E-01,
     + 0.1199303E-01, 0.1321268E-01, 0.1449927E-01, 0.1585780E-01,
     + 0.1729329E-01, 0.1881074E-01, 0.2041512E-01, 0.2211138E-01,
     + 0.2390440E-01, 0.2579898E-01, 0.2779984E-01, 0.2991158E-01,
     + 0.3213867E-01, 0.3448538E-01, 0.3695584E-01, 0.3955395E-01,
     + 0.4228335E-01, 0.4514745E-01, 0.4814936E-01, 0.5129187E-01,
     + 0.5457743E-01, 0.5800815E-01, 0.6158573E-01, 0.6531150E-01,
     + 0.6918633E-01, 0.7321068E-01, 0.7738455E-01, 0.8170749E-01,
     + 0.8617855E-01, 0.9079633E-01, 0.9555893E-01, 0.1004640E+00,
     + 0.1055086E+00, 0.1106896E+00, 0.1160030E+00, 0.1214447E+00,
     + 0.1270100E+00, 0.1326939E+00, 0.1384908E+00, 0.1443949E+00,
     + 0.1504000E+00, 0.1564997E+00, 0.1626871E+00, 0.1689552E+00,
     + 0.1752968E+00, 0.1817044E+00, 0.1881704E+00, 0.1946871E+00,
     + 0.2012467E+00, 0.2078415E+00, 0.2144635E+00, 0.2211049E+00,
     + 0.2277580E+00, 0.2344151E+00, 0.2410687E+00, 0.2477113E+00,
     + 0.2543357E+00, 0.2609349E+00, 0.2675021E+00, 0.2740307E+00,
     + 0.2805143E+00, 0.2869471E+00, 0.2933231E+00, 0.2996370E+00,
     + 0.3058836E+00, 0.3120581E+00, 0.3181560E+00, 0.3241731E+00,
     + 0.3301056E+00, 0.3359498E+00, 0.3417027E+00, 0.3473612E+00,
     + 0.3529229E+00, 0.3583853E+00, 0.3637465E+00, 0.3690047E+00,
     + 0.3741586E+00, 0.3792069E+00, 0.3841486E+00, 0.3889832E+00,
     + 0.3937100E+00, 0.3983290E+00, 0.4028399E+00, 0.4072429E+00,
     + 0.4115383E+00, 0.4157267E+00, 0.4198085E+00, 0.4237846E+00,
     + 0.4276559E+00, 0.4314232E+00, 0.4350879E+00, 0.4386510E+00,
     + 0.4421138E+00, 0.4454777E+00, 0.4487441E+00, 0.4519145E+00,
     + 0.4549904E+00, 0.4579734E+00, 0.4608651E+00, 0.4636671E+00,
     + 0.4663810E+00, 0.4690086E+00, 0.4715515E+00, 0.4740113E+00,
     + 0.4763899E+00, 0.4786888E+00, 0.4809097E+00, 0.4830544E+00,
     + 0.4851244E+00, 0.4871213E+00, 0.4890469E+00, 0.4909027E+00,
     + 0.4926903E+00, 0.4944113E+00, 0.4960671E+00, 0.4976594E+00,
     + 0.4991895E+00, 0.5006589E+00, 0.5020691E+00, 0.5034215E+00,
     + 0.5047174E+00, 0.5059581E+00, 0.5071451E+00, 0.5082795E+00,
     + 0.5093627E+00, 0.5103958E+00, 0.5113801E+00, 0.5123168E+00,
     + 0.5132069E+00, 0.5140517E+00, 0.5148521E+00, 0.5156093E+00,
     + 0.5163244E+00, 0.5169982E+00, 0.5176319E+00, 0.5182264E+00,
     + 0.5187826E+00, 0.5193014E+00, 0.5197838E+00, 0.5202307E+00,
     + 0.5206428E+00, 0.5210210E+00, 0.5213661E+00, 0.5216790E+00,
     + 0.5219604E+00, 0.5222110E+00, 0.5224316E+00, 0.5226229E+00,
     + 0.5227856E+00, 0.5229203E+00, 0.5230279E+00, 0.5231088E+00,
     + 0.5231638E+00, 0.5231934E+00, 0.5231983E+00, 0.5231790E+00,
     + 0.5231362E+00, 0.5230704E+00, 0.5229821E+00, 0.5228719E+00,
     + 0.5227403E+00, 0.5225878E+00, 0.5224149E+00, 0.5222221E+00,
     + 0.5220100E+00, 0.5217788E+00, 0.5215292E+00, 0.5212615E+00,
     + 0.5209762E+00, 0.5206737E+00, 0.5203544E+00, 0.5200188E+00,
     + 0.5196671E+00, 0.5192999E+00, 0.5189175E+00, 0.5185202E+00,
     + 0.5181085E+00, 0.5176826E+00, 0.5172430E+00, 0.5167899E+00,
     + 0.5163237E+00, 0.5158448E+00, 0.5153533E+00, 0.5148497E+00,
     + 0.5143342E+00, 0.5138072E+00, 0.5132689E+00, 0.5127196E+00,
     + 0.5121596E+00, 0.5115891E+00, 0.5110085E+00, 0.5104179E+00,
     + 0.5098177E+00, 0.5092080E+00, 0.5085892E+00, 0.5079614E+00,
     + 0.5073249E+00, 0.5066799E+00, 0.5060266E+00, 0.5053653E+00,
     + 0.5046962E+00, 0.5040194E+00, 0.5033352E+00, 0.5026438E+00,
     + 0.5019453E+00, 0.5012400E+00, 0.5005280E+00, 0.4998096E+00,
     + 0.4990848E+00, 0.4983539E+00, 0.4976171E+00, 0.4968744E+00,
     + 0.4961261E+00, 0.4953724E+00, 0.4946133E+00, 0.4938491E+00,
     + 0.4930799E+00, 0.4923058E+00, 0.4915269E+00, 0.4907435E+00,
     + 0.4899556E+00, 0.4891634E+00, 0.4883670E+00, 0.4875666E+00,
     + 0.4867622E+00, 0.4859540E+00, 0.4851421E+00, 0.4843266E+00,
     + 0.4835076E+00, 0.4826853E+00, 0.4818598E+00, 0.4810311E+00,
     + 0.4801994E+00, 0.4793648E+00, 0.4785273E+00, 0.4776871E+00,
     + 0.4768443E+00, 0.4759989E+00, 0.4751511E+00, 0.4743009E+00,
     + 0.4734484E+00, 0.4725938E+00, 0.4717370E+00, 0.4708782E+00,
     + 0.4700175E+00, 0.4691549E+00, 0.4682906E+00, 0.4674245E+00,
     + 0.4665568E+00, 0.4656875E+00, 0.4648167E+00, 0.4639445E+00,
     + 0.4630710E+00, 0.4621961E+00, 0.4613200E+00, 0.4604428E+00,
     + 0.4595644E+00, 0.4586850E+00, 0.4578047E+00, 0.4569233E+00,
     + 0.4560412E+00, 0.4551582E+00, 0.4542744E+00, 0.4533900E+00,
     + 0.4525049E+00, 0.4516191E+00, 0.4507329E+00, 0.4498461E+00,
     + 0.4489588E+00, 0.4480712E+00, 0.4471832E+00, 0.4462948E+00,
     + 0.4454062E+00, 0.4445173E+00, 0.4436282E+00, 0.4427390E+00,
     + 0.4418497E+00, 0.4409603E+00, 0.4400708E+00, 0.4391813E+00,
     + 0.4382919E+00, 0.4374025E+00, 0.4365133E+00, 0.4356242E+00,
     + 0.4347352E+00, 0.4338465E+00, 0.4329580E+00, 0.4320697E+00,
     + 0.4311818E+00, 0.4302942E+00, 0.4294069E+00, 0.4285200E+00,
     + 0.4276336E+00, 0.4267475E+00, 0.4258620E+00, 0.4249769E+00,
     + 0.4240923E+00, 0.4232083E+00, 0.4223249E+00, 0.4214420E+00,
     + 0.4205598E+00, 0.4196782E+00, 0.4187972E+00, 0.4179170E+00,
     + 0.4170374E+00, 0.4161586E+00, 0.4152805E+00, 0.4144031E+00,
     + 0.4135266E+00, 0.4126508E+00, 0.4117759E+00, 0.4109018E+00,
     + 0.4100285E+00, 0.4091561E+00, 0.4082847E+00, 0.4074141E+00,
     + 0.4065444E+00, 0.4056757E+00, 0.4048079E+00, 0.4039411E+00,
     + 0.4030753E+00, 0.4022105E+00, 0.4013467E+00, 0.4004839E+00,
     + 0.3996222E+00, 0.3987615E+00, 0.3979018E+00, 0.3970433E+00,
     + 0.3961858E+00, 0.3953294E+00, 0.3944742E+00, 0.3936201E+00,
     + 0.3927671E+00, 0.3919152E+00, 0.3910645E+00, 0.3902150E+00,
     + 0.3893666E+00, 0.3885195E+00, 0.3876735E+00, 0.3868287E+00,
     + 0.3859852E+00, 0.3851428E+00, 0.3843017E+00, 0.3834619E+00,
     + 0.3826232E+00, 0.3817859E+00, 0.3809498E+00, 0.3801150E+00,
     + 0.3792814E+00, 0.3784491E+00, 0.3776182E+00, 0.3767885E+00,
     + 0.3759601E+00, 0.3751331E+00, 0.3743073E+00, 0.3734829E+00,
     + 0.3726598E+00, 0.3718380E+00, 0.3710176E+00, 0.3701986E+00,
     + 0.3693809E+00, 0.3685645E+00, 0.3677495E+00, 0.3669359E+00,
     + 0.3661236E+00, 0.3653127E+00, 0.3645032E+00, 0.3636951E+00,
     + 0.3628884E+00, 0.3620830E+00, 0.3612791E+00, 0.3604766E+00,
     + 0.3596754E+00, 0.3588757E+00, 0.3580774E+00, 0.3572805E+00,
     + 0.3564850E+00, 0.3556909E+00, 0.3548983E+00, 0.3541070E+00,
     + 0.3533172E+00, 0.3525289E+00, 0.3517420E+00, 0.3509565E+00,
     + 0.3501724E+00, 0.3493898E+00, 0.3486086E+00, 0.3478289E+00,
     + 0.3470506E+00, 0.3462738E+00, 0.3454984E+00, 0.3447244E+00,
     + 0.3439520E+00, 0.3431809E+00, 0.3424114E+00, 0.3416433E+00,
     + 0.3408766E+00, 0.3401114E+00, 0.3393477E+00, 0.3385854E+00,
     + 0.3378246E+00, 0.3370653E+00, 0.3363074E+00, 0.3355510E+00,
     + 0.3347960E+00, 0.3340425E+00, 0.3332905E+00, 0.3325400E+00,
     + 0.3317909E+00, 0.3310433E+00, 0.3302971E+00, 0.3295525E+00,
     + 0.3288093E+00, 0.3280675E+00, 0.3273273E+00, 0.3265885E+00,
     + 0.3258512E+00, 0.3251153E+00, 0.3243809E+00, 0.3236480E+00,
     + 0.3229166E+00, 0.3221866E+00, 0.3214581E+00, 0.3207310E+00,
     + 0.3200055E+00, 0.3192814E+00, 0.3185587E+00, 0.3178376E+00,
     + 0.3171179E+00, 0.3163996E+00, 0.3156828E+00, 0.3149675E+00,
     + 0.3142537E+00, 0.3135413E+00, 0.3128304E+00, 0.3121209E+00,
     + 0.3114129E+00, 0.3107064E+00, 0.3100013E+00, 0.3092976E+00,
     + 0.3085955E+00, 0.3078948E+00, 0.3071955E+00, 0.3064977E+00,
     + 0.3058013E+00, 0.3051064E+00, 0.3044130E+00, 0.3037209E+00,
     + 0.3030304E+00, 0.3023413E+00, 0.3016536E+00, 0.3009674E+00,
     + 0.3002826E+00, 0.2995992E+00, 0.2989173E+00, 0.2982368E+00,
     + 0.2975578E+00, 0.2968802E+00, 0.2962040E+00, 0.2955293E+00,
     + 0.2948560E+00, 0.2941841E+00, 0.2935137E+00, 0.2928446E+00,
     + 0.2921770E+00, 0.2915108E+00, 0.2908461E+00, 0.2901827E+00,
     + 0.2895208E+00, 0.2888603E+00, 0.2882012E+00, 0.2875435E+00,
     + 0.2868872E+00, 0.2862324E+00, 0.2855789E+00, 0.2849269E+00,
     + 0.2842762E+00, 0.2836269E+00, 0.2829791E+00, 0.2823326E+00,
     + 0.2816876E+00, 0.2810439E+00, 0.2804016E+00, 0.2797607E+00,
     + 0.2791212E+00, 0.2784831E+00, 0.2778464E+00, 0.2772110E+00,
     + 0.2765770E+00, 0.2759444E+00, 0.2753132E+00, 0.2746834E+00,
     + 0.2740549E+00, 0.2734278E+00, 0.2728020E+00, 0.2721776E+00,
     + 0.2715546E+00, 0.2709330E+00, 0.2703127E+00, 0.2696937E+00,
     + 0.2690762E+00, 0.2684599E+00, 0.2678450E+00, 0.2672315E+00,
     + 0.2666193E+00, 0.2660085E+00, 0.2653990E+00, 0.2647908E+00,
     + 0.2641840E+00, 0.2635785E+00, 0.2629743E+00, 0.2623715E+00,
     + 0.2617700E+00, 0.2611698E+00, 0.2605709E+00, 0.2599734E+00,
     + 0.2593772E+00, 0.2587823E+00, 0.2581887E+00, 0.2575965E+00,
     + 0.2570055E+00, 0.2564158E+00, 0.2558275E+00, 0.2552405E+00,
     + 0.2546547E+00, 0.2540703E+00, 0.2534871E+00, 0.2529053E+00,
     + 0.2523247E+00, 0.2517454E+00, 0.2511674E+00, 0.2505907E+00,
     + 0.2500153E+00, 0.2494412E+00, 0.2488683E+00, 0.2482967E+00,
     + 0.2477264E+00, 0.2471574E+00, 0.2465896E+00, 0.2460231E+00,
     + 0.2454579E+00, 0.2448939E+00, 0.2443312E+00, 0.2437697E+00,
     + 0.2432095E+00, 0.2426505E+00, 0.2420928E+00, 0.2415364E+00,
     + 0.2409812E+00, 0.2404272E+00, 0.2398745E+00, 0.2393230E+00,
     + 0.2387727E+00, 0.2382237E+00, 0.2376759E+00, 0.2371293E+00,
     + 0.2365840E+00, 0.2360398E+00, 0.2354970E+00, 0.2349553E+00,
     + 0.2344148E+00, 0.2338756E+00, 0.2333375E+00, 0.2328007E+00,
     + 0.2322651E+00, 0.2317307E+00, 0.2311975E+00, 0.2306655E+00,
     + 0.2301346E+00, 0.2296050E+00, 0.2290766E+00, 0.2285494E+00,
     + 0.2280233E+00, 0.2274985E+00, 0.2269748E+00, 0.2264523E+00,
     + 0.2259310E+00, 0.2254108E+00, 0.2248919E+00, 0.2243741E+00,
     + 0.2238574E+00, 0.2233420E+00, 0.2228277E+00, 0.2223146E+00,
     + 0.2218026E+00, 0.2212918E+00, 0.2207821E+00, 0.2202736E+00,
     + 0.2197662E+00, 0.2192600E+00, 0.2187550E+00, 0.2182511E+00,
     + 0.2177483E+00, 0.2172466E+00, 0.2167461E+00, 0.2162468E+00,
     + 0.2157485E+00, 0.2152514E+00, 0.2147554E+00, 0.2142606E+00,
     + 0.2137669E+00, 0.2132742E+00, 0.2127827E+00, 0.2122924E+00,
     + 0.2118031E+00, 0.2113149E+00, 0.2108279E+00, 0.2103419E+00,
     + 0.2098571E+00, 0.2093734E+00, 0.2088907E+00, 0.2084092E+00,
     + 0.2079287E+00, 0.2074494E+00, 0.2069711E+00, 0.2064939E+00,
     + 0.2060178E+00, 0.2055428E+00, 0.2050689E+00, 0.2045961E+00,
     + 0.2041243E+00, 0.2036536E+00, 0.2031840E+00, 0.2027154E+00,
     + 0.2022479E+00, 0.2017815E+00, 0.2013161E+00, 0.2008518E+00,
     + 0.2003886E+00, 0.1999264E+00, 0.1994652E+00, 0.1990052E+00,
     + 0.1985461E+00, 0.1980881E+00, 0.1976312E+00, 0.1971753E+00,
     + 0.1967204E+00, 0.1962666E+00, 0.1958138E+00, 0.1953620E+00,
     + 0.1949113E+00, 0.1944616E+00, 0.1940130E+00, 0.1935653E+00,
     + 0.1931187E+00, 0.1926731E+00, 0.1922285E+00, 0.1917849E+00,
     + 0.1913424E+00, 0.1909008E+00, 0.1904603E+00, 0.1900208E+00,
     + 0.1895822E+00, 0.1891447E+00, 0.1887082E+00, 0.1882727E+00,
     + 0.1878381E+00, 0.1874046E+00, 0.1869721E+00, 0.1865405E+00,
     + 0.1861100E+00, 0.1856804E+00, 0.1852518E+00, 0.1848242E+00,
     + 0.1843975E+00, 0.1839719E+00, 0.1835472E+00, 0.1831235E+00,
     + 0.1827007E+00, 0.1822789E+00, 0.1818581E+00, 0.1814383E+00,
     + 0.1810194E+00, 0.1806015E+00, 0.1801845E+00, 0.1797685E+00,
     + 0.1793534E+00, 0.1789393E+00, 0.1785262E+00, 0.1781140E+00,
     + 0.1777027E+00, 0.1772924E+00, 0.1768830E+00, 0.1764745E+00,
     + 0.1760670E+00, 0.1756604E+00, 0.1752548E+00, 0.1748501E+00,
     + 0.1744463E+00, 0.1740435E+00, 0.1736415E+00, 0.1732405E+00,
     + 0.1728404E+00, 0.1724412E+00, 0.1720430E+00, 0.1716456E+00,
     + 0.1712492E+00, 0.1708537E+00, 0.1704591E+00, 0.1700654E+00,
     + 0.1696726E+00, 0.1692807E+00, 0.1688897E+00, 0.1684996E+00,
     + 0.1681103E+00, 0.1677220E+00, 0.1673346E+00, 0.1669481E+00,
     + 0.1665624E+00, 0.1661777E+00, 0.1657938E+00, 0.1654108E+00,
     + 0.1650287E+00, 0.1646474E+00, 0.1642671E+00, 0.1638876E+00,
     + 0.1635090E+00, 0.1631312E+00, 0.1627544E+00, 0.1623784E+00,
     + 0.1620032E+00, 0.1616289E+00, 0.1612555E+00, 0.1608829E+00,
     + 0.1605112E+00, 0.1601404E+00, 0.1597704E+00, 0.1594012E+00,
     + 0.1590329E+00, 0.1586655E+00, 0.1582989E+00, 0.1579331E+00,
     + 0.1575682E+00, 0.1572041E+00, 0.1568409E+00, 0.1564784E+00,
     + 0.1561169E+00, 0.1557561E+00, 0.1553962E+00, 0.1550371E+00,
     + 0.1546789E+00, 0.1543214E+00, 0.1539648E+00, 0.1536090E+00,
     + 0.1532540E+00, 0.1528999E+00, 0.1525465E+00, 0.1521940E+00,
     + 0.1518423E+00, 0.1514914E+00, 0.1511413E+00, 0.1507920E+00,
     + 0.1504435E+00, 0.1500958E+00, 0.1497489E+00, 0.1494028E+00,
     + 0.1490576E+00, 0.1487131E+00, 0.1483694E+00, 0.1480264E+00,
     + 0.1476843E+00, 0.1473430E+00, 0.1470024E+00, 0.1466627E+00,
     + 0.1463237E+00, 0.1459855E+00, 0.1456481E+00, 0.1453114E+00,
     + 0.1449756E+00, 0.1446405E+00, 0.1443062E+00, 0.1439726E+00,
     + 0.1436398E+00, 0.1433078E+00, 0.1429766E+00, 0.1426461E+00,
     + 0.1423163E+00, 0.1419874E+00, 0.1416592E+00, 0.1413317E+00,
     + 0.1410050E+00, 0.1406791E+00, 0.1403539E+00, 0.1400294E+00,
     + 0.1397058E+00, 0.1393828E+00, 0.1390606E+00, 0.1387391E+00,
     + 0.1384184E+00, 0.1380984E+00, 0.1377792E+00, 0.1374607E+00,
     + 0.1371429E+00, 0.1368259E+00, 0.1365096E+00, 0.1361940E+00,
     + 0.1358791E+00, 0.1355650E+00, 0.1352516E+00, 0.1349389E+00,
     + 0.1346270E+00, 0.1343157E+00, 0.1340052E+00, 0.1336954E+00,
     + 0.1333863E+00, 0.1330779E+00, 0.1327703E+00, 0.1324633E+00,
     + 0.1321571E+00, 0.1318515E+00, 0.1315467E+00, 0.1312425E+00,
     + 0.1309391E+00, 0.1306364E+00, 0.1303343E+00, 0.1300330E+00,
     + 0.1297324E+00, 0.1294324E+00, 0.1291332E+00, 0.1288346E+00,
     + 0.1285367E+00, 0.1282395E+00, 0.1279430E+00, 0.1276472E+00,
     + 0.1273521E+00, 0.1270576E+00, 0.1267638E+00, 0.1264707E+00,
     + 0.1261783E+00, 0.1258866E+00, 0.1255955E+00, 0.1253051E+00,
     + 0.1250154E+00, 0.1247263E+00, 0.1244379E+00, 0.1241502E+00,
     + 0.1238631E+00, 0.1235767E+00, 0.1232910E+00, 0.1230059E+00,
     + 0.1227214E+00, 0.1224377E+00, 0.1221546E+00, 0.1218721E+00,
     + 0.1215903E+00, 0.1213091E+00, 0.1210286E+00, 0.1207488E+00,
     + 0.1204696E+00, 0.1201910E+00, 0.1199131E+00, 0.1196358E+00,
     + 0.1193591E+00, 0.1190831E+00, 0.1188078E+00, 0.1185330E+00,
     + 0.1182589E+00, 0.1179855E+00, 0.1177126E+00, 0.1174404E+00,
     + 0.1171688E+00, 0.1168979E+00, 0.1166276E+00, 0.1163579E+00,
     + 0.1160888E+00, 0.1158203E+00, 0.1155525E+00, 0.1152853E+00,
     + 0.1150187E+00, 0.1147527E+00, 0.1144873E+00, 0.1142226E+00,
     + 0.1139584E+00, 0.1136949E+00, 0.1134320E+00, 0.1131696E+00,
     + 0.1129079E+00, 0.1126468E+00, 0.1123863E+00, 0.1121264E+00,
     + 0.1118671E+00, 0.1116084E+00, 0.1113503E+00, 0.1110928E+00,
     + 0.1108359E+00, 0.1105795E+00, 0.1103238E+00, 0.1100687E+00,
     + 0.1098141E+00, 0.1095602E+00, 0.1093068E+00, 0.1090540E+00,
     + 0.1088018E+00, 0.1085502E+00, 0.1082991E+00, 0.1080487E+00,
     + 0.1077988E+00, 0.1075495E+00, 0.1073007E+00, 0.1070526E+00,
     + 0.1068050E+00, 0.1065580E+00, 0.1063115E+00, 0.1060657E+00,
     + 0.1058204E+00, 0.1055756E+00, 0.1053315E+00, 0.1050879E+00,
     + 0.1048448E+00, 0.1046023E+00, 0.1043604E+00, 0.1041190E+00,
     + 0.1038782E+00, 0.1036380E+00, 0.1033983E+00, 0.1031592E+00,
     + 0.1029206E+00, 0.1026825E+00, 0.1024451E+00, 0.1022081E+00,
     + 0.1019717E+00, 0.1017359E+00, 0.1015006E+00, 0.1012658E+00,
     + 0.1010316E+00, 0.1007979E+00, 0.1005648E+00, 0.1003322E+00,
     + 0.1001002E+00, 0.9986865E-01, 0.9963767E-01, 0.9940722E-01,
     + 0.9917730E-01, 0.9894792E-01, 0.9871906E-01, 0.9849073E-01,
     + 0.9826293E-01, 0.9803566E-01, 0.9780892E-01, 0.9758269E-01,
     + 0.9735699E-01, 0.9713181E-01, 0.9690716E-01, 0.9668302E-01,
     + 0.9645940E-01, 0.9623629E-01, 0.9601370E-01, 0.9579163E-01,
     + 0.9557007E-01, 0.9534902E-01, 0.9512848E-01, 0.9490846E-01,
     + 0.9468894E-01, 0.9446993E-01, 0.9425142E-01, 0.9403342E-01,
     + 0.9381592E-01, 0.9359893E-01, 0.9338244E-01, 0.9316645E-01,
     + 0.9295095E-01, 0.9273596E-01, 0.9252146E-01, 0.9230746E-01,
     + 0.9209395E-01, 0.9188094E-01, 0.9166842E-01, 0.9145639E-01,
     + 0.9124485E-01, 0.9103380E-01, 0.9082324E-01, 0.9061316E-01,
     + 0.9040357E-01, 0.9019447E-01, 0.8998584E-01, 0.8977770E-01,
     + 0.8957004E-01, 0.8936287E-01, 0.8915617E-01, 0.8894994E-01,
     + 0.8874420E-01, 0.8853893E-01, 0.8833413E-01, 0.8812981E-01,
     + 0.8792596E-01, 0.8772259E-01, 0.8751968E-01, 0.8731724E-01,
     + 0.8711527E-01, 0.8691377E-01, 0.8671273E-01, 0.8651215E-01,
     + 0.8631205E-01, 0.8611240E-01, 0.8591321E-01, 0.8571449E-01,
     + 0.8551622E-01, 0.8531842E-01, 0.8512107E-01, 0.8492418E-01,
     + 0.8472774E-01, 0.8453175E-01, 0.8433622E-01, 0.8414114E-01,
     + 0.8394652E-01, 0.8375234E-01, 0.8355861E-01, 0.8336533E-01,
     + 0.8317250E-01, 0.8298011E-01, 0.8278817E-01, 0.8259667E-01,
     + 0.8240561E-01, 0.8221500E-01, 0.8202482E-01, 0.8183509E-01,
     + 0.8164579E-01, 0.8145694E-01, 0.8126852E-01, 0.8108053E-01,
     + 0.8089298E-01, 0.8070586E-01, 0.8051918E-01, 0.8033293E-01,
     + 0.8014710E-01, 0.7996171E-01, 0.7977675E-01, 0.7959221E-01,
     + 0.7940810E-01, 0.7922442E-01, 0.7904116E-01, 0.7885833E-01,
     + 0.7867591E-01, 0.7849392E-01, 0.7831235E-01, 0.7813120E-01,
     + 0.7795047E-01, 0.7777016E-01, 0.7759026E-01, 0.7741078E-01,
     + 0.7723172E-01, 0.7705307E-01, 0.7687483E-01, 0.7669701E-01,
     + 0.7651959E-01, 0.7634259E-01, 0.7616599E-01, 0.7598981E-01,
     + 0.7581403E-01, 0.7563865E-01, 0.7546369E-01, 0.7528913E-01,
     + 0.7511497E-01, 0.7494121E-01, 0.7476786E-01, 0.7459490E-01,
     + 0.7442235E-01, 0.7425020E-01, 0.7407844E-01, 0.7390708E-01,
     + 0.7373612E-01, 0.7356555E-01, 0.7339538E-01, 0.7322560E-01,
     + 0.7305621E-01, 0.7288722E-01, 0.7271862E-01, 0.7255040E-01,
     + 0.7238258E-01, 0.7221514E-01, 0.7204809E-01, 0.7188143E-01,
     + 0.7171515E-01, 0.7154925E-01, 0.7138374E-01, 0.7121862E-01,
     + 0.7105387E-01, 0.7088951E-01, 0.7072552E-01, 0.7056192E-01,
     + 0.7039869E-01, 0.7023584E-01, 0.7007337E-01, 0.6991128E-01,
     + 0.6974955E-01, 0.6958821E-01, 0.6942723E-01, 0.6926663E-01,
     + 0.6910640E-01, 0.6894654E-01, 0.6878705E-01, 0.6862793E-01,
     + 0.6846917E-01, 0.6831079E-01, 0.6815277E-01, 0.6799511E-01,
     + 0.6783782E-01, 0.6768090E-01, 0.6752433E-01, 0.6736813E-01,
     + 0.6721229E-01, 0.6705681E-01, 0.6690169E-01, 0.6674693E-01,
     + 0.6659253E-01, 0.6643848E-01, 0.6628479E-01, 0.6613146E-01,
     + 0.6597848E-01, 0.6582585E-01, 0.6567358E-01, 0.6552166E-01,
     + 0.6537009E-01, 0.6521887E-01, 0.6506800E-01, 0.6491748E-01,
     + 0.6476731E-01, 0.6461749E-01, 0.6446801E-01, 0.6431888E-01,
     + 0.6417009E-01, 0.6402165E-01, 0.6387355E-01, 0.6372579E-01,
     + 0.6357837E-01, 0.6343130E-01, 0.6328456E-01, 0.6313817E-01,
     + 0.6299211E-01, 0.6284639E-01, 0.6270101E-01, 0.6255597E-01,
     + 0.6241126E-01, 0.6226688E-01, 0.6212284E-01, 0.6197913E-01,
     + 0.6183575E-01, 0.6169271E-01, 0.6155000E-01, 0.6140761E-01,
     + 0.6126556E-01, 0.6112384E-01, 0.6098244E-01, 0.6084137E-01,
     + 0.6070062E-01, 0.6056020E-01, 0.6042011E-01, 0.6028034E-01,
     + 0.6014089E-01, 0.6000177E-01, 0.5986297E-01, 0.5972449E-01,
     + 0.5958633E-01, 0.5944849E-01, 0.5931096E-01, 0.5917376E-01,
     + 0.5903687E-01, 0.5890030E-01, 0.5876405E-01, 0.5862811E-01,
     + 0.5849248E-01, 0.5835717E-01, 0.5822217E-01, 0.5808749E-01,
     + 0.5795311E-01, 0.5781905E-01, 0.5768530E-01, 0.5755185E-01,
     + 0.5741872E-01, 0.5728589E-01, 0.5715337E-01, 0.5702115E-01,
     + 0.5688925E-01, 0.5675764E-01, 0.5662635E-01, 0.5649535E-01,
     + 0.5636466E-01, 0.5623427E-01, 0.5610418E-01, 0.5597440E-01,
     + 0.5584491E-01, 0.5571572E-01, 0.5558683E-01, 0.5545824E-01,
     + 0.5532995E-01, 0.5520195E-01, 0.5507425E-01, 0.5494685E-01,
     + 0.5481974E-01, 0.5469292E-01, 0.5456640E-01, 0.5444017E-01,
     + 0.5431423E-01, 0.5418859E-01, 0.5406323E-01, 0.5393817E-01,
     + 0.5381339E-01, 0.5368890E-01, 0.5356470E-01, 0.5344079E-01,
     + 0.5331716E-01, 0.5319382E-01, 0.5307077E-01, 0.5294800E-01,
     + 0.5282551E-01, 0.5270331E-01, 0.5258139E-01, 0.5245975E-01,
     + 0.5233839E-01, 0.5221732E-01, 0.5209652E-01, 0.5197600E-01,
     + 0.5185576E-01, 0.5173580E-01, 0.5161612E-01, 0.5149672E-01,
     + 0.5137759E-01, 0.5125873E-01, 0.5114016E-01, 0.5102185E-01,
     + 0.5090382E-01, 0.5078606E-01, 0.5066858E-01, 0.5055136E-01,
     + 0.5043442E-01, 0.5031775E-01, 0.5020135E-01, 0.5008521E-01,
     + 0.4996935E-01, 0.4985375E-01, 0.4973842E-01, 0.4962336E-01,
     + 0.4950857E-01, 0.4939404E-01, 0.4927977E-01, 0.4916577E-01,
     + 0.4905203E-01, 0.4893856E-01, 0.4882535E-01, 0.4871240E-01,
     + 0.4859971E-01, 0.4848728E-01, 0.4837511E-01, 0.4826320E-01,
     + 0.4815155E-01, 0.4804016E-01, 0.4792903E-01, 0.4781815E-01,
     + 0.4770753E-01, 0.4759717E-01, 0.4748706E-01, 0.4737721E-01,
     + 0.4726761E-01, 0.4715826E-01, 0.4704917E-01, 0.4694032E-01,
     + 0.4683173E-01, 0.4672340E-01, 0.4661531E-01, 0.4650747E-01,
     + 0.4639988E-01, 0.4629254E-01, 0.4618545E-01, 0.4607861E-01,
     + 0.4597201E-01, 0.4586566E-01, 0.4575956E-01, 0.4565370E-01,
     + 0.4554809E-01, 0.4544272E-01, 0.4533760E-01, 0.4523271E-01,
     + 0.4512807E-01, 0.4502368E-01, 0.4491952E-01, 0.4481561E-01,
     + 0.4471193E-01, 0.4460850E-01, 0.4450530E-01, 0.4440234E-01,
     + 0.4429963E-01, 0.4419715E-01, 0.4409490E-01, 0.4399289E-01,
     + 0.4389112E-01, 0.4378959E-01, 0.4368829E-01, 0.4358722E-01,
     + 0.4348639E-01, 0.4338579E-01, 0.4328542E-01, 0.4318528E-01,
     + 0.4308538E-01, 0.4298571E-01, 0.4288627E-01, 0.4278706E-01,
     + 0.4268808E-01, 0.4258932E-01, 0.4249080E-01, 0.4239250E-01,
     + 0.4229443E-01, 0.4219659E-01, 0.4209897E-01, 0.4200158E-01,
     + 0.4190442E-01, 0.4180748E-01, 0.4171076E-01, 0.4161427E-01,
     + 0.4151800E-01, 0.4142195E-01, 0.4132613E-01, 0.4123053E-01,
     + 0.4113515E-01, 0.4103999E-01, 0.4094505E-01, 0.4085033E-01,
     + 0.4075582E-01, 0.4066154E-01, 0.4056748E-01, 0.4047363E-01,
     + 0.4038000E-01, 0.4028658E-01, 0.4019339E-01, 0.4010040E-01,
     + 0.4000764E-01, 0.3991508E-01, 0.3982275E-01, 0.3973062E-01,
     + 0.3963871E-01, 0.3954701E-01, 0.3945552E-01, 0.3936425E-01,
     + 0.3927319E-01, 0.3918233E-01, 0.3909169E-01, 0.3900126E-01,
     + 0.3891103E-01, 0.3882102E-01, 0.3873121E-01, 0.3864161E-01,
     + 0.3855222E-01, 0.3846303E-01, 0.3837405E-01, 0.3828528E-01,
     + 0.3819671E-01, 0.3810835E-01, 0.3802019E-01, 0.3793223E-01,
     + 0.3784448E-01, 0.3775693E-01, 0.3766959E-01, 0.3758244E-01,
     + 0.3749550E-01, 0.3740876E-01, 0.3732222E-01, 0.3723588E-01,
     + 0.3714974E-01, 0.3706380E-01, 0.3697806E-01, 0.3689251E-01,
     + 0.3680717E-01, 0.3672202E-01, 0.3663707E-01, 0.3655231E-01,
     + 0.3646775E-01, 0.3638339E-01, 0.3629922E-01, 0.3621525E-01,
     + 0.3613147E-01, 0.3604788E-01, 0.3596449E-01, 0.3588129E-01,
     + 0.3579828E-01, 0.3571547E-01, 0.3563284E-01, 0.3555041E-01,
     + 0.3546817E-01, 0.3538612E-01, 0.3530426E-01, 0.3522259E-01,
     + 0.3514110E-01, 0.3505981E-01, 0.3497870E-01, 0.3489778E-01,
     + 0.3481705E-01, 0.3473651E-01, 0.3465615E-01, 0.3457598E-01,
     + 0.3449599E-01, 0.3441619E-01, 0.3433657E-01, 0.3425714E-01,
     + 0.3417789E-01, 0.3409882E-01, 0.3401994E-01, 0.3394123E-01,
     + 0.3386272E-01, 0.3378438E-01, 0.3370622E-01, 0.3362825E-01,
     + 0.3355045E-01, 0.3347284E-01, 0.3339540E-01, 0.3331815E-01,
     + 0.3324107E-01, 0.3316417E-01, 0.3308745E-01, 0.3301090E-01,
     + 0.3293454E-01, 0.3285835E-01, 0.3278233E-01, 0.3270650E-01,
     + 0.3263083E-01, 0.3255535E-01, 0.3248003E-01, 0.3240489E-01,
     + 0.3232993E-01, 0.3225514E-01, 0.3218052E-01, 0.3210607E-01,
     + 0.3203180E-01, 0.3195770E-01, 0.3188377E-01, 0.3181001E-01,
     + 0.3173642E-01, 0.3166300E-01, 0.3158975E-01, 0.3151667E-01,
     + 0.3144376E-01, 0.3137102E-01, 0.3129845E-01, 0.3122605E-01,
     + 0.3115381E-01, 0.3108174E-01, 0.3100983E-01, 0.3093810E-01,
     + 0.3086652E-01, 0.3079512E-01, 0.3072388E-01, 0.3065280E-01,
     + 0.3058189E-01, 0.3051114E-01, 0.3044056E-01, 0.3037014E-01,
     + 0.3029988E-01, 0.3022978E-01, 0.3015985E-01, 0.3009008E-01,
     + 0.3002047E-01, 0.2995102E-01, 0.2988173E-01, 0.2981261E-01,
     + 0.2974364E-01, 0.2967483E-01, 0.2960618E-01, 0.2953769E-01,
     + 0.2946936E-01, 0.2940118E-01, 0.2933317E-01, 0.2926531E-01,
     + 0.2919761E-01, 0.2913006E-01, 0.2906267E-01, 0.2899544E-01,
     + 0.2892836E-01, 0.2886144E-01, 0.2879467E-01, 0.2872806E-01,
     + 0.2866160E-01, 0.2859529E-01, 0.2852914E-01, 0.2846314E-01,
     + 0.2839730E-01, 0.2833160E-01, 0.2826606E-01, 0.2820067E-01,
     + 0.2813543E-01, 0.2807034E-01, 0.2800541E-01, 0.2794062E-01,
     + 0.2787598E-01, 0.2781149E-01, 0.2774716E-01, 0.2768297E-01,
     + 0.2761892E-01, 0.2755503E-01, 0.2749129E-01, 0.2742769E-01/
c     -------------------------------------------------------------------
c     d-state radial wave function: r=0 to 15 fm, steps of 0.01 fm
      data uavd/0.0d0,
     + 0.5404087E-06, 0.4249015E-05, 0.1410349E-04, 0.3290019E-04,
     + 0.6328127E-04, 0.1077597E-03, 0.1687417E-03, 0.2485474E-03,
     + 0.3494285E-03, 0.4735840E-03, 0.6231737E-03, 0.8003291E-03,
     + 0.1007162E-02, 0.1245774E-02, 0.1518254E-02, 0.1826691E-02,
     + 0.2173166E-02, 0.2559756E-02, 0.2988529E-02, 0.3461540E-02,
     + 0.3980826E-02, 0.4548396E-02, 0.5166226E-02, 0.5836246E-02,
     + 0.6560331E-02, 0.7340290E-02, 0.8177850E-02, 0.9074650E-02,
     + 0.1003222E-01, 0.1105198E-01, 0.1213520E-01, 0.1328303E-01,
     + 0.1449645E-01, 0.1577626E-01, 0.1712311E-01, 0.1853744E-01,
     + 0.2001948E-01, 0.2156928E-01, 0.2318665E-01, 0.2487119E-01,
     + 0.2662228E-01, 0.2843906E-01, 0.3032046E-01, 0.3226518E-01,
     + 0.3427167E-01, 0.3633819E-01, 0.3846279E-01, 0.4064327E-01,
     + 0.4287728E-01, 0.4516226E-01, 0.4749548E-01, 0.4987402E-01,
     + 0.5229486E-01, 0.5475479E-01, 0.5725052E-01, 0.5977863E-01,
     + 0.6233563E-01, 0.6491795E-01, 0.6752197E-01, 0.7014403E-01,
     + 0.7278045E-01, 0.7542753E-01, 0.7808161E-01, 0.8073905E-01,
     + 0.8339623E-01, 0.8604960E-01, 0.8869570E-01, 0.9133112E-01,
     + 0.9395255E-01, 0.9655680E-01, 0.9914078E-01, 0.1017015E+00,
     + 0.1042362E+00, 0.1067420E+00, 0.1092165E+00, 0.1116573E+00,
     + 0.1140620E+00, 0.1164284E+00, 0.1187548E+00, 0.1210391E+00,
     + 0.1232798E+00, 0.1254753E+00, 0.1276242E+00, 0.1297253E+00,
     + 0.1317774E+00, 0.1337797E+00, 0.1357313E+00, 0.1376314E+00,
     + 0.1394796E+00, 0.1412754E+00, 0.1430184E+00, 0.1447084E+00,
     + 0.1463452E+00, 0.1479289E+00, 0.1494595E+00, 0.1509372E+00,
     + 0.1523621E+00, 0.1537347E+00, 0.1550552E+00, 0.1563242E+00,
     + 0.1575420E+00, 0.1587094E+00, 0.1598268E+00, 0.1608950E+00,
     + 0.1619146E+00, 0.1628864E+00, 0.1638111E+00, 0.1646896E+00,
     + 0.1655226E+00, 0.1663110E+00, 0.1670556E+00, 0.1677574E+00,
     + 0.1684172E+00, 0.1690358E+00, 0.1696143E+00, 0.1701536E+00,
     + 0.1706544E+00, 0.1711178E+00, 0.1715447E+00, 0.1719359E+00,
     + 0.1722924E+00, 0.1726150E+00, 0.1729047E+00, 0.1731623E+00,
     + 0.1733887E+00, 0.1735848E+00, 0.1737514E+00, 0.1738894E+00,
     + 0.1739995E+00, 0.1740827E+00, 0.1741396E+00, 0.1741711E+00,
     + 0.1741779E+00, 0.1741609E+00, 0.1741207E+00, 0.1740580E+00,
     + 0.1739736E+00, 0.1738682E+00, 0.1737424E+00, 0.1735970E+00,
     + 0.1734325E+00, 0.1732496E+00, 0.1730489E+00, 0.1728310E+00,
     + 0.1725965E+00, 0.1723459E+00, 0.1720799E+00, 0.1717990E+00,
     + 0.1715036E+00, 0.1711944E+00, 0.1708718E+00, 0.1705363E+00,
     + 0.1701883E+00, 0.1698284E+00, 0.1694569E+00, 0.1690744E+00,
     + 0.1686812E+00, 0.1682778E+00, 0.1678645E+00, 0.1674417E+00,
     + 0.1670098E+00, 0.1665692E+00, 0.1661203E+00, 0.1656633E+00,
     + 0.1651986E+00, 0.1647266E+00, 0.1642475E+00, 0.1637616E+00,
     + 0.1632693E+00, 0.1627708E+00, 0.1622664E+00, 0.1617564E+00,
     + 0.1612410E+00, 0.1607205E+00, 0.1601952E+00, 0.1596652E+00,
     + 0.1591308E+00, 0.1585923E+00, 0.1580497E+00, 0.1575035E+00,
     + 0.1569536E+00, 0.1564005E+00, 0.1558441E+00, 0.1552848E+00,
     + 0.1547227E+00, 0.1541580E+00, 0.1535908E+00, 0.1530213E+00,
     + 0.1524497E+00, 0.1518760E+00, 0.1513006E+00, 0.1507234E+00,
     + 0.1501447E+00, 0.1495645E+00, 0.1489830E+00, 0.1484004E+00,
     + 0.1478167E+00, 0.1472321E+00, 0.1466466E+00, 0.1460605E+00,
     + 0.1454737E+00, 0.1448864E+00, 0.1442987E+00, 0.1437107E+00,
     + 0.1431225E+00, 0.1425341E+00, 0.1419457E+00, 0.1413574E+00,
     + 0.1407691E+00, 0.1401811E+00, 0.1395933E+00, 0.1390059E+00,
     + 0.1384189E+00, 0.1378324E+00, 0.1372465E+00, 0.1366611E+00,
     + 0.1360765E+00, 0.1354926E+00, 0.1349095E+00, 0.1343272E+00,
     + 0.1337459E+00, 0.1331655E+00, 0.1325861E+00, 0.1320078E+00,
     + 0.1314305E+00, 0.1308544E+00, 0.1302796E+00, 0.1297059E+00,
     + 0.1291335E+00, 0.1285624E+00, 0.1279927E+00, 0.1274244E+00,
     + 0.1268574E+00, 0.1262920E+00, 0.1257280E+00, 0.1251655E+00,
     + 0.1246046E+00, 0.1240452E+00, 0.1234875E+00, 0.1229314E+00,
     + 0.1223769E+00, 0.1218241E+00, 0.1212730E+00, 0.1207236E+00,
     + 0.1201760E+00, 0.1196301E+00, 0.1190860E+00, 0.1185438E+00,
     + 0.1180033E+00, 0.1174647E+00, 0.1169279E+00, 0.1163930E+00,
     + 0.1158600E+00, 0.1153288E+00, 0.1147996E+00, 0.1142723E+00,
     + 0.1137470E+00, 0.1132236E+00, 0.1127021E+00, 0.1121826E+00,
     + 0.1116651E+00, 0.1111496E+00, 0.1106361E+00, 0.1101246E+00,
     + 0.1096151E+00, 0.1091076E+00, 0.1086021E+00, 0.1080987E+00,
     + 0.1075973E+00, 0.1070979E+00, 0.1066006E+00, 0.1061053E+00,
     + 0.1056121E+00, 0.1051209E+00, 0.1046318E+00, 0.1041448E+00,
     + 0.1036598E+00, 0.1031769E+00, 0.1026961E+00, 0.1022173E+00,
     + 0.1017406E+00, 0.1012660E+00, 0.1007935E+00, 0.1003230E+00,
     + 0.9985458E-01, 0.9938824E-01, 0.9892397E-01, 0.9846177E-01,
     + 0.9800164E-01, 0.9754358E-01, 0.9708758E-01, 0.9663365E-01,
     + 0.9618177E-01, 0.9573196E-01, 0.9528420E-01, 0.9483849E-01,
     + 0.9439483E-01, 0.9395322E-01, 0.9351365E-01, 0.9307613E-01,
     + 0.9264064E-01, 0.9220718E-01, 0.9177575E-01, 0.9134634E-01,
     + 0.9091895E-01, 0.9049358E-01, 0.9007022E-01, 0.8964887E-01,
     + 0.8922951E-01, 0.8881215E-01, 0.8839679E-01, 0.8798341E-01,
     + 0.8757201E-01, 0.8716258E-01, 0.8675512E-01, 0.8634963E-01,
     + 0.8594609E-01, 0.8554451E-01, 0.8514487E-01, 0.8474717E-01,
     + 0.8435141E-01, 0.8395757E-01, 0.8356566E-01, 0.8317566E-01,
     + 0.8278757E-01, 0.8240138E-01, 0.8201708E-01, 0.8163468E-01,
     + 0.8125415E-01, 0.8087551E-01, 0.8049873E-01, 0.8012381E-01,
     + 0.7975075E-01, 0.7937953E-01, 0.7901016E-01, 0.7864262E-01,
     + 0.7827691E-01, 0.7791302E-01, 0.7755094E-01, 0.7719066E-01,
     + 0.7683219E-01, 0.7647550E-01, 0.7612060E-01, 0.7576747E-01,
     + 0.7541612E-01, 0.7506652E-01, 0.7471868E-01, 0.7437259E-01,
     + 0.7402824E-01, 0.7368562E-01, 0.7334472E-01, 0.7300554E-01,
     + 0.7266807E-01, 0.7233231E-01, 0.7199824E-01, 0.7166585E-01,
     + 0.7133515E-01, 0.7100612E-01, 0.7067875E-01, 0.7035304E-01,
     + 0.7002898E-01, 0.6970656E-01, 0.6938578E-01, 0.6906662E-01,
     + 0.6874909E-01, 0.6843316E-01, 0.6811884E-01, 0.6780612E-01,
     + 0.6749499E-01, 0.6718544E-01, 0.6687746E-01, 0.6657106E-01,
     + 0.6626621E-01, 0.6596291E-01, 0.6566116E-01, 0.6536095E-01,
     + 0.6506226E-01, 0.6476510E-01, 0.6446945E-01, 0.6417532E-01,
     + 0.6388268E-01, 0.6359153E-01, 0.6330187E-01, 0.6301369E-01,
     + 0.6272698E-01, 0.6244173E-01, 0.6215794E-01, 0.6187560E-01,
     + 0.6159470E-01, 0.6131523E-01, 0.6103720E-01, 0.6076058E-01,
     + 0.6048537E-01, 0.6021157E-01, 0.5993917E-01, 0.5966816E-01,
     + 0.5939854E-01, 0.5913029E-01, 0.5886341E-01, 0.5859790E-01,
     + 0.5833374E-01, 0.5807094E-01, 0.5780947E-01, 0.5754934E-01,
     + 0.5729054E-01, 0.5703306E-01, 0.5677689E-01, 0.5652203E-01,
     + 0.5626848E-01, 0.5601622E-01, 0.5576524E-01, 0.5551555E-01,
     + 0.5526713E-01, 0.5501998E-01, 0.5477408E-01, 0.5452945E-01,
     + 0.5428606E-01, 0.5404391E-01, 0.5380299E-01, 0.5356331E-01,
     + 0.5332484E-01, 0.5308759E-01, 0.5285155E-01, 0.5261671E-01,
     + 0.5238306E-01, 0.5215060E-01, 0.5191933E-01, 0.5168924E-01,
     + 0.5146031E-01, 0.5123255E-01, 0.5100594E-01, 0.5078049E-01,
     + 0.5055618E-01, 0.5033301E-01, 0.5011098E-01, 0.4989007E-01,
     + 0.4967028E-01, 0.4945161E-01, 0.4923405E-01, 0.4901759E-01,
     + 0.4880222E-01, 0.4858795E-01, 0.4837476E-01, 0.4816266E-01,
     + 0.4795162E-01, 0.4774166E-01, 0.4753275E-01, 0.4732490E-01,
     + 0.4711810E-01, 0.4691235E-01, 0.4670763E-01, 0.4650395E-01,
     + 0.4630130E-01, 0.4609967E-01, 0.4589905E-01, 0.4569945E-01,
     + 0.4550085E-01, 0.4530325E-01, 0.4510664E-01, 0.4491103E-01,
     + 0.4471640E-01, 0.4452274E-01, 0.4433007E-01, 0.4413835E-01,
     + 0.4394760E-01, 0.4375781E-01, 0.4356897E-01, 0.4338108E-01,
     + 0.4319412E-01, 0.4300811E-01, 0.4282302E-01, 0.4263887E-01,
     + 0.4245563E-01, 0.4227331E-01, 0.4209190E-01, 0.4191139E-01,
     + 0.4173179E-01, 0.4155308E-01, 0.4137527E-01, 0.4119834E-01,
     + 0.4102229E-01, 0.4084712E-01, 0.4067282E-01, 0.4049939E-01,
     + 0.4032682E-01, 0.4015511E-01, 0.3998425E-01, 0.3981424E-01,
     + 0.3964507E-01, 0.3947674E-01, 0.3930925E-01, 0.3914258E-01,
     + 0.3897674E-01, 0.3881172E-01, 0.3864752E-01, 0.3848413E-01,
     + 0.3832154E-01, 0.3815976E-01, 0.3799878E-01, 0.3783859E-01,
     + 0.3767919E-01, 0.3752058E-01, 0.3736275E-01, 0.3720569E-01,
     + 0.3704941E-01, 0.3689389E-01, 0.3673914E-01, 0.3658515E-01,
     + 0.3643191E-01, 0.3627943E-01, 0.3612769E-01, 0.3597669E-01,
     + 0.3582644E-01, 0.3567692E-01, 0.3552813E-01, 0.3538007E-01,
     + 0.3523273E-01, 0.3508611E-01, 0.3494021E-01, 0.3479501E-01,
     + 0.3465053E-01, 0.3450674E-01, 0.3436366E-01, 0.3422128E-01,
     + 0.3407958E-01, 0.3393857E-01, 0.3379825E-01, 0.3365861E-01,
     + 0.3351965E-01, 0.3338136E-01, 0.3324374E-01, 0.3310678E-01,
     + 0.3297049E-01, 0.3283485E-01, 0.3269987E-01, 0.3256555E-01,
     + 0.3243187E-01, 0.3229883E-01, 0.3216644E-01, 0.3203468E-01,
     + 0.3190356E-01, 0.3177307E-01, 0.3164320E-01, 0.3151396E-01,
     + 0.3138534E-01, 0.3125733E-01, 0.3112994E-01, 0.3100316E-01,
     + 0.3087699E-01, 0.3075142E-01, 0.3062645E-01, 0.3050208E-01,
     + 0.3037830E-01, 0.3025511E-01, 0.3013251E-01, 0.3001049E-01,
     + 0.2988905E-01, 0.2976820E-01, 0.2964791E-01, 0.2952820E-01,
     + 0.2940906E-01, 0.2929048E-01, 0.2917247E-01, 0.2905501E-01,
     + 0.2893811E-01, 0.2882176E-01, 0.2870597E-01, 0.2859072E-01,
     + 0.2847602E-01, 0.2836185E-01, 0.2824823E-01, 0.2813514E-01,
     + 0.2802258E-01, 0.2791056E-01, 0.2779906E-01, 0.2768809E-01,
     + 0.2757763E-01, 0.2746770E-01, 0.2735828E-01, 0.2724937E-01,
     + 0.2714098E-01, 0.2703309E-01, 0.2692571E-01, 0.2681882E-01,
     + 0.2671244E-01, 0.2660656E-01, 0.2650117E-01, 0.2639627E-01,
     + 0.2629186E-01, 0.2618793E-01, 0.2608449E-01, 0.2598153E-01,
     + 0.2587905E-01, 0.2577704E-01, 0.2567551E-01, 0.2557445E-01,
     + 0.2547386E-01, 0.2537373E-01, 0.2527407E-01, 0.2517486E-01,
     + 0.2507612E-01, 0.2497783E-01, 0.2487999E-01, 0.2478261E-01,
     + 0.2468567E-01, 0.2458919E-01, 0.2449314E-01, 0.2439754E-01,
     + 0.2430237E-01, 0.2420765E-01, 0.2411335E-01, 0.2401949E-01,
     + 0.2392606E-01, 0.2383306E-01, 0.2374048E-01, 0.2364833E-01,
     + 0.2355660E-01, 0.2346528E-01, 0.2337439E-01, 0.2328391E-01,
     + 0.2319384E-01, 0.2310417E-01, 0.2301492E-01, 0.2292608E-01,
     + 0.2283763E-01, 0.2274959E-01, 0.2266195E-01, 0.2257471E-01,
     + 0.2248786E-01, 0.2240140E-01, 0.2231534E-01, 0.2222966E-01,
     + 0.2214437E-01, 0.2205947E-01, 0.2197495E-01, 0.2189081E-01,
     + 0.2180705E-01, 0.2172366E-01, 0.2164065E-01, 0.2155802E-01,
     + 0.2147575E-01, 0.2139386E-01, 0.2131233E-01, 0.2123117E-01,
     + 0.2115037E-01, 0.2106993E-01, 0.2098985E-01, 0.2091013E-01,
     + 0.2083077E-01, 0.2075176E-01, 0.2067310E-01, 0.2059480E-01,
     + 0.2051684E-01, 0.2043923E-01, 0.2036196E-01, 0.2028504E-01,
     + 0.2020846E-01, 0.2013222E-01, 0.2005631E-01, 0.1998075E-01,
     + 0.1990552E-01, 0.1983062E-01, 0.1975605E-01, 0.1968181E-01,
     + 0.1960790E-01, 0.1953431E-01, 0.1946105E-01, 0.1938811E-01,
     + 0.1931550E-01, 0.1924320E-01, 0.1917122E-01, 0.1909956E-01,
     + 0.1902821E-01, 0.1895717E-01, 0.1888645E-01, 0.1881604E-01,
     + 0.1874593E-01, 0.1867613E-01, 0.1860664E-01, 0.1853744E-01,
     + 0.1846856E-01, 0.1839997E-01, 0.1833168E-01, 0.1826369E-01,
     + 0.1819599E-01, 0.1812859E-01, 0.1806148E-01, 0.1799466E-01,
     + 0.1792813E-01, 0.1786189E-01, 0.1779594E-01, 0.1773028E-01,
     + 0.1766489E-01, 0.1759979E-01, 0.1753497E-01, 0.1747044E-01,
     + 0.1740618E-01, 0.1734219E-01, 0.1727849E-01, 0.1721505E-01,
     + 0.1715189E-01, 0.1708900E-01, 0.1702639E-01, 0.1696404E-01,
     + 0.1690195E-01, 0.1684014E-01, 0.1677859E-01, 0.1671730E-01,
     + 0.1665627E-01, 0.1659551E-01, 0.1653500E-01, 0.1647476E-01,
     + 0.1641477E-01, 0.1635503E-01, 0.1629555E-01, 0.1623633E-01,
     + 0.1617735E-01, 0.1611863E-01, 0.1606016E-01, 0.1600193E-01,
     + 0.1594395E-01, 0.1588622E-01, 0.1582873E-01, 0.1577148E-01,
     + 0.1571448E-01, 0.1565772E-01, 0.1560119E-01, 0.1554491E-01,
     + 0.1548886E-01, 0.1543305E-01, 0.1537747E-01, 0.1532213E-01,
     + 0.1526702E-01, 0.1521215E-01, 0.1515750E-01, 0.1510308E-01,
     + 0.1504889E-01, 0.1499493E-01, 0.1494119E-01, 0.1488768E-01,
     + 0.1483439E-01, 0.1478133E-01, 0.1472848E-01, 0.1467586E-01,
     + 0.1462346E-01, 0.1457127E-01, 0.1451930E-01, 0.1446755E-01,
     + 0.1441602E-01, 0.1436469E-01, 0.1431358E-01, 0.1426269E-01,
     + 0.1421200E-01, 0.1416152E-01, 0.1411126E-01, 0.1406120E-01,
     + 0.1401134E-01, 0.1396170E-01, 0.1391226E-01, 0.1386302E-01,
     + 0.1381399E-01, 0.1376515E-01, 0.1371652E-01, 0.1366809E-01,
     + 0.1361986E-01, 0.1357183E-01, 0.1352399E-01, 0.1347635E-01,
     + 0.1342890E-01, 0.1338165E-01, 0.1333460E-01, 0.1328773E-01,
     + 0.1324106E-01, 0.1319458E-01, 0.1314829E-01, 0.1310218E-01,
     + 0.1305627E-01, 0.1301054E-01, 0.1296500E-01, 0.1291964E-01,
     + 0.1287447E-01, 0.1282948E-01, 0.1278467E-01, 0.1274005E-01,
     + 0.1269561E-01, 0.1265134E-01, 0.1260726E-01, 0.1256335E-01,
     + 0.1251963E-01, 0.1247607E-01, 0.1243270E-01, 0.1238950E-01,
     + 0.1234647E-01, 0.1230362E-01, 0.1226094E-01, 0.1221843E-01,
     + 0.1217610E-01, 0.1213393E-01, 0.1209193E-01, 0.1205010E-01,
     + 0.1200844E-01, 0.1196695E-01, 0.1192562E-01, 0.1188446E-01,
     + 0.1184346E-01, 0.1180262E-01, 0.1176195E-01, 0.1172144E-01,
     + 0.1168110E-01, 0.1164091E-01, 0.1160089E-01, 0.1156102E-01,
     + 0.1152131E-01, 0.1148176E-01, 0.1144237E-01, 0.1140313E-01,
     + 0.1136405E-01, 0.1132513E-01, 0.1128635E-01, 0.1124774E-01,
     + 0.1120927E-01, 0.1117096E-01, 0.1113280E-01, 0.1109479E-01,
     + 0.1105692E-01, 0.1101921E-01, 0.1098165E-01, 0.1094424E-01,
     + 0.1090697E-01, 0.1086985E-01, 0.1083287E-01, 0.1079604E-01,
     + 0.1075936E-01, 0.1072282E-01, 0.1068642E-01, 0.1065016E-01,
     + 0.1061405E-01, 0.1057808E-01, 0.1054225E-01, 0.1050655E-01,
     + 0.1047100E-01, 0.1043559E-01, 0.1040031E-01, 0.1036518E-01,
     + 0.1033018E-01, 0.1029531E-01, 0.1026058E-01, 0.1022599E-01,
     + 0.1019153E-01, 0.1015720E-01, 0.1012301E-01, 0.1008895E-01,
     + 0.1005502E-01, 0.1002123E-01, 0.9987562E-02, 0.9954027E-02,
     + 0.9920621E-02, 0.9887344E-02, 0.9854196E-02, 0.9821175E-02,
     + 0.9788282E-02, 0.9755516E-02, 0.9722875E-02, 0.9690361E-02,
     + 0.9657971E-02, 0.9625706E-02, 0.9593564E-02, 0.9561546E-02,
     + 0.9529651E-02, 0.9497878E-02, 0.9466227E-02, 0.9434697E-02,
     + 0.9403288E-02, 0.9371998E-02, 0.9340829E-02, 0.9309778E-02,
     + 0.9278846E-02, 0.9248031E-02, 0.9217334E-02, 0.9186754E-02,
     + 0.9156291E-02, 0.9125943E-02, 0.9095711E-02, 0.9065593E-02,
     + 0.9035590E-02, 0.9005701E-02, 0.8975925E-02, 0.8946262E-02,
     + 0.8916711E-02, 0.8887272E-02, 0.8857945E-02, 0.8828728E-02,
     + 0.8799622E-02, 0.8770625E-02, 0.8741738E-02, 0.8712960E-02,
     + 0.8684290E-02, 0.8655729E-02, 0.8627275E-02, 0.8598927E-02,
     + 0.8570687E-02, 0.8542552E-02, 0.8514523E-02, 0.8486600E-02,
     + 0.8458781E-02, 0.8431066E-02, 0.8403455E-02, 0.8375947E-02,
     + 0.8348543E-02, 0.8321240E-02, 0.8294040E-02, 0.8266941E-02,
     + 0.8239944E-02, 0.8213047E-02, 0.8186250E-02, 0.8159553E-02,
     + 0.8132956E-02, 0.8106457E-02, 0.8080057E-02, 0.8053755E-02,
     + 0.8027551E-02, 0.8001444E-02, 0.7975434E-02, 0.7949520E-02,
     + 0.7923702E-02, 0.7897980E-02, 0.7872353E-02, 0.7846821E-02,
     + 0.7821383E-02, 0.7796039E-02, 0.7770789E-02, 0.7745632E-02,
     + 0.7720567E-02, 0.7695595E-02, 0.7670715E-02, 0.7645926E-02,
     + 0.7621229E-02, 0.7596622E-02, 0.7572106E-02, 0.7547680E-02,
     + 0.7523343E-02, 0.7499096E-02, 0.7474938E-02, 0.7450868E-02,
     + 0.7426886E-02, 0.7402992E-02, 0.7379185E-02, 0.7355465E-02,
     + 0.7331832E-02, 0.7308285E-02, 0.7284824E-02, 0.7261449E-02,
     + 0.7238159E-02, 0.7214953E-02, 0.7191832E-02, 0.7168795E-02,
     + 0.7145842E-02, 0.7122972E-02, 0.7100186E-02, 0.7077482E-02,
     + 0.7054860E-02, 0.7032320E-02, 0.7009862E-02, 0.6987486E-02,
     + 0.6965190E-02, 0.6942975E-02, 0.6920840E-02, 0.6898785E-02,
     + 0.6876810E-02, 0.6854914E-02, 0.6833097E-02, 0.6811359E-02,
     + 0.6789699E-02, 0.6768117E-02, 0.6746613E-02, 0.6725186E-02,
     + 0.6703835E-02, 0.6682562E-02, 0.6661365E-02, 0.6640244E-02,
     + 0.6619199E-02, 0.6598229E-02, 0.6577335E-02, 0.6556515E-02,
     + 0.6535769E-02, 0.6515098E-02, 0.6494501E-02, 0.6473977E-02,
     + 0.6453527E-02, 0.6433149E-02, 0.6412844E-02, 0.6392612E-02,
     + 0.6372452E-02, 0.6352363E-02, 0.6332346E-02, 0.6312400E-02,
     + 0.6292525E-02, 0.6272720E-02, 0.6252986E-02, 0.6233322E-02,
     + 0.6213728E-02, 0.6194203E-02, 0.6174747E-02, 0.6155360E-02,
     + 0.6136042E-02, 0.6116792E-02, 0.6097610E-02, 0.6078496E-02,
     + 0.6059449E-02, 0.6040470E-02, 0.6021558E-02, 0.6002712E-02,
     + 0.5983932E-02, 0.5965219E-02, 0.5946572E-02, 0.5927990E-02,
     + 0.5909473E-02, 0.5891022E-02, 0.5872635E-02, 0.5854313E-02,
     + 0.5836056E-02, 0.5817862E-02, 0.5799732E-02, 0.5781665E-02,
     + 0.5763662E-02, 0.5745722E-02, 0.5727844E-02, 0.5710029E-02,
     + 0.5692277E-02, 0.5674586E-02, 0.5656957E-02, 0.5639389E-02,
     + 0.5621883E-02, 0.5604438E-02, 0.5587053E-02, 0.5569729E-02,
     + 0.5552465E-02, 0.5535262E-02, 0.5518118E-02, 0.5501033E-02,
     + 0.5484008E-02, 0.5467042E-02, 0.5450135E-02, 0.5433286E-02,
     + 0.5416496E-02, 0.5399764E-02, 0.5383090E-02, 0.5366474E-02,
     + 0.5349914E-02, 0.5333413E-02, 0.5316968E-02, 0.5300580E-02,
     + 0.5284248E-02, 0.5267973E-02, 0.5251754E-02, 0.5235591E-02,
     + 0.5219484E-02, 0.5203431E-02, 0.5187435E-02, 0.5171493E-02,
     + 0.5155606E-02, 0.5139773E-02, 0.5123995E-02, 0.5108271E-02,
     + 0.5092601E-02, 0.5076984E-02, 0.5061422E-02, 0.5045912E-02,
     + 0.5030455E-02, 0.5015052E-02, 0.4999701E-02, 0.4984402E-02,
     + 0.4969156E-02, 0.4953962E-02, 0.4938820E-02, 0.4923729E-02,
     + 0.4908690E-02, 0.4893702E-02, 0.4878765E-02, 0.4863879E-02,
     + 0.4849043E-02, 0.4834258E-02, 0.4819523E-02, 0.4804839E-02,
     + 0.4790204E-02, 0.4775619E-02, 0.4761083E-02, 0.4746597E-02,
     + 0.4732159E-02, 0.4717771E-02, 0.4703431E-02, 0.4689140E-02,
     + 0.4674898E-02, 0.4660703E-02, 0.4646556E-02, 0.4632457E-02,
     + 0.4618406E-02, 0.4604402E-02, 0.4590446E-02, 0.4576536E-02,
     + 0.4562673E-02, 0.4548857E-02, 0.4535088E-02, 0.4521365E-02,
     + 0.4507688E-02, 0.4494057E-02, 0.4480472E-02, 0.4466932E-02,
     + 0.4453438E-02, 0.4439989E-02, 0.4426585E-02, 0.4413226E-02,
     + 0.4399912E-02, 0.4386643E-02, 0.4373418E-02, 0.4360237E-02,
     + 0.4347100E-02, 0.4334008E-02, 0.4320959E-02, 0.4307953E-02,
     + 0.4294991E-02, 0.4282072E-02, 0.4269197E-02, 0.4256364E-02,
     + 0.4243574E-02, 0.4230826E-02, 0.4218121E-02, 0.4205459E-02,
     + 0.4192838E-02, 0.4180259E-02, 0.4167723E-02, 0.4155227E-02,
     + 0.4142774E-02, 0.4130361E-02, 0.4117990E-02, 0.4105660E-02,
     + 0.4093370E-02, 0.4081122E-02, 0.4068914E-02, 0.4056746E-02,
     + 0.4044619E-02, 0.4032531E-02, 0.4020484E-02, 0.4008476E-02,
     + 0.3996508E-02, 0.3984580E-02, 0.3972691E-02, 0.3960841E-02,
     + 0.3949030E-02, 0.3937258E-02, 0.3925525E-02, 0.3913831E-02,
     + 0.3902175E-02, 0.3890557E-02, 0.3878977E-02, 0.3867436E-02,
     + 0.3855932E-02, 0.3844467E-02, 0.3833038E-02, 0.3821648E-02,
     + 0.3810294E-02, 0.3798978E-02, 0.3787699E-02, 0.3776457E-02,
     + 0.3765252E-02, 0.3754083E-02, 0.3742951E-02, 0.3731855E-02,
     + 0.3720796E-02, 0.3709773E-02, 0.3698785E-02, 0.3687834E-02,
     + 0.3676918E-02, 0.3666038E-02, 0.3655193E-02, 0.3644384E-02,
     + 0.3633609E-02, 0.3622870E-02, 0.3612166E-02, 0.3601497E-02,
     + 0.3590862E-02, 0.3580262E-02, 0.3569696E-02, 0.3559165E-02,
     + 0.3548667E-02, 0.3538204E-02, 0.3527775E-02, 0.3517379E-02,
     + 0.3507018E-02, 0.3496689E-02, 0.3486395E-02, 0.3476133E-02,
     + 0.3465905E-02, 0.3455710E-02, 0.3445547E-02, 0.3435418E-02,
     + 0.3425321E-02, 0.3415257E-02, 0.3405225E-02, 0.3395226E-02,
     + 0.3385259E-02, 0.3375324E-02, 0.3365421E-02, 0.3355550E-02,
     + 0.3345711E-02, 0.3335903E-02, 0.3326127E-02, 0.3316382E-02,
     + 0.3306669E-02, 0.3296987E-02, 0.3287335E-02, 0.3277715E-02,
     + 0.3268126E-02, 0.3258568E-02, 0.3249040E-02, 0.3239542E-02,
     + 0.3230075E-02, 0.3220639E-02, 0.3211232E-02, 0.3201856E-02,
     + 0.3192510E-02, 0.3183193E-02, 0.3173907E-02, 0.3164650E-02,
     + 0.3155422E-02, 0.3146224E-02, 0.3137055E-02, 0.3127916E-02,
     + 0.3118806E-02, 0.3109724E-02, 0.3100672E-02, 0.3091648E-02,
     + 0.3082654E-02, 0.3073687E-02, 0.3064750E-02, 0.3055840E-02,
     + 0.3046960E-02, 0.3038107E-02, 0.3029282E-02, 0.3020485E-02,
     + 0.3011717E-02, 0.3002976E-02, 0.2994263E-02, 0.2985577E-02,
     + 0.2976919E-02, 0.2968288E-02, 0.2959685E-02, 0.2951109E-02,
     + 0.2942560E-02, 0.2934038E-02, 0.2925543E-02, 0.2917075E-02,
     + 0.2908633E-02, 0.2900218E-02, 0.2891830E-02, 0.2883468E-02,
     + 0.2875133E-02, 0.2866824E-02, 0.2858541E-02, 0.2850284E-02,
     + 0.2842053E-02, 0.2833848E-02, 0.2825669E-02, 0.2817515E-02,
     + 0.2809387E-02, 0.2801285E-02, 0.2793208E-02, 0.2785157E-02,
     + 0.2777131E-02, 0.2769130E-02, 0.2761154E-02, 0.2753203E-02,
     + 0.2745277E-02, 0.2737376E-02, 0.2729500E-02, 0.2721648E-02,
     + 0.2713821E-02, 0.2706018E-02, 0.2698240E-02, 0.2690486E-02,
     + 0.2682756E-02, 0.2675051E-02, 0.2667369E-02, 0.2659712E-02,
     + 0.2652078E-02, 0.2644469E-02, 0.2636883E-02, 0.2629320E-02,
     + 0.2621781E-02, 0.2614266E-02, 0.2606774E-02, 0.2599305E-02,
     + 0.2591860E-02, 0.2584438E-02, 0.2577039E-02, 0.2569663E-02,
     + 0.2562309E-02, 0.2554979E-02, 0.2547672E-02, 0.2540387E-02,
     + 0.2533124E-02, 0.2525885E-02, 0.2518667E-02, 0.2511472E-02,
     + 0.2504300E-02, 0.2497149E-02, 0.2490021E-02, 0.2482915E-02,
     + 0.2475831E-02, 0.2468768E-02, 0.2461728E-02, 0.2454709E-02,
     + 0.2447712E-02, 0.2440737E-02, 0.2433783E-02, 0.2426850E-02,
     + 0.2419939E-02, 0.2413050E-02, 0.2406181E-02, 0.2399334E-02,
     + 0.2392508E-02, 0.2385703E-02, 0.2378919E-02, 0.2372155E-02,
     + 0.2365413E-02, 0.2358691E-02, 0.2351990E-02, 0.2345310E-02,
     + 0.2338650E-02, 0.2332010E-02, 0.2325391E-02, 0.2318792E-02,
     + 0.2312214E-02, 0.2305656E-02, 0.2299117E-02, 0.2292599E-02,
     + 0.2286101E-02, 0.2279623E-02, 0.2273164E-02, 0.2266726E-02,
     + 0.2260307E-02, 0.2253908E-02, 0.2247528E-02, 0.2241168E-02,
     + 0.2234827E-02, 0.2228506E-02, 0.2222203E-02, 0.2215921E-02,
     + 0.2209657E-02, 0.2203413E-02, 0.2197187E-02, 0.2190981E-02,
     + 0.2184793E-02, 0.2178624E-02, 0.2172474E-02, 0.2166343E-02,
     + 0.2160231E-02, 0.2154137E-02, 0.2148061E-02, 0.2142005E-02,
     + 0.2135966E-02, 0.2129946E-02, 0.2123944E-02, 0.2117961E-02,
     + 0.2111995E-02, 0.2106048E-02, 0.2100119E-02, 0.2094208E-02,
     + 0.2088314E-02, 0.2082439E-02, 0.2076581E-02, 0.2070741E-02,
     + 0.2064919E-02, 0.2059115E-02, 0.2053328E-02, 0.2047558E-02,
     + 0.2041806E-02, 0.2036072E-02, 0.2030354E-02, 0.2024654E-02,
     + 0.2018972E-02, 0.2013306E-02, 0.2007657E-02, 0.2002026E-02,
     + 0.1996412E-02, 0.1990814E-02, 0.1985233E-02, 0.1979670E-02,
     + 0.1974123E-02, 0.1968592E-02, 0.1963079E-02, 0.1957582E-02,
     + 0.1952101E-02, 0.1946637E-02, 0.1941190E-02, 0.1935758E-02,
     + 0.1930344E-02, 0.1924945E-02, 0.1919563E-02, 0.1914196E-02,
     + 0.1908846E-02, 0.1903512E-02, 0.1898194E-02, 0.1892892E-02,
     + 0.1887606E-02, 0.1882336E-02, 0.1877082E-02, 0.1871843E-02,
     + 0.1866620E-02, 0.1861412E-02, 0.1856221E-02, 0.1851044E-02,
     + 0.1845884E-02, 0.1840738E-02, 0.1835608E-02, 0.1830494E-02,
     + 0.1825394E-02, 0.1820310E-02, 0.1815241E-02, 0.1810188E-02,
     + 0.1805149E-02, 0.1800125E-02, 0.1795117E-02, 0.1790123E-02,
     + 0.1785144E-02, 0.1780180E-02, 0.1775231E-02, 0.1770297E-02,
     + 0.1765377E-02, 0.1760472E-02, 0.1755581E-02, 0.1750705E-02,
     + 0.1745844E-02, 0.1740997E-02, 0.1736164E-02, 0.1731346E-02,
     + 0.1726542E-02, 0.1721753E-02, 0.1716977E-02, 0.1712216E-02,
     + 0.1707469E-02, 0.1702736E-02, 0.1698017E-02, 0.1693312E-02,
     + 0.1688621E-02, 0.1683944E-02, 0.1679281E-02, 0.1674631E-02,
     + 0.1669996E-02, 0.1665374E-02, 0.1660765E-02, 0.1656171E-02,
     + 0.1651590E-02, 0.1647022E-02, 0.1642468E-02, 0.1637928E-02,
     + 0.1633401E-02, 0.1628887E-02, 0.1624387E-02, 0.1619900E-02,
     + 0.1615426E-02, 0.1610965E-02, 0.1606518E-02, 0.1602083E-02,
     + 0.1597662E-02, 0.1593254E-02, 0.1588858E-02, 0.1584476E-02,
     + 0.1580107E-02, 0.1575750E-02, 0.1571406E-02, 0.1567075E-02,
     + 0.1562757E-02, 0.1558452E-02, 0.1554159E-02, 0.1549879E-02,
     + 0.1545611E-02, 0.1541356E-02, 0.1537113E-02, 0.1532883E-02,
     + 0.1528665E-02, 0.1524460E-02, 0.1520267E-02, 0.1516086E-02,
     + 0.1511918E-02, 0.1507762E-02, 0.1503618E-02, 0.1499486E-02,
     + 0.1495366E-02, 0.1491258E-02, 0.1487162E-02, 0.1483079E-02,
     + 0.1479007E-02, 0.1474947E-02, 0.1470899E-02, 0.1466863E-02,
     + 0.1462839E-02, 0.1458826E-02, 0.1454825E-02, 0.1450836E-02/
c     ------------------------------------------------------------
c     derivative of s-state wave function
      data uasp/0.7937903E-01,
     + 0.7937903E-01, 0.8006949E-01, 0.8122170E-01, 0.8283769E-01,
     + 0.8492012E-01, 0.8747218E-01, 0.9049745E-01, 0.9399978E-01,
     + 0.9798318E-01, 0.1024517E+00, 0.1074091E+00, 0.1128592E+00,
     + 0.1188050E+00, 0.1252492E+00, 0.1321935E+00, 0.1396387E+00,
     + 0.1475846E+00, 0.1560296E+00, 0.1649706E+00, 0.1744030E+00,
     + 0.1843203E+00, 0.1947140E+00, 0.2055736E+00, 0.2168864E+00,
     + 0.2286373E+00, 0.2408088E+00, 0.2533809E+00, 0.2663312E+00,
     + 0.2796346E+00, 0.2932635E+00, 0.3071879E+00, 0.3213751E+00,
     + 0.3357902E+00, 0.3503961E+00, 0.3651535E+00, 0.3800212E+00,
     + 0.3949563E+00, 0.4099142E+00, 0.4248492E+00, 0.4397145E+00,
     + 0.4544627E+00, 0.4690457E+00, 0.4834156E+00, 0.4975242E+00,
     + 0.5113243E+00, 0.5247691E+00, 0.5378130E+00, 0.5504120E+00,
     + 0.5625236E+00, 0.5741074E+00, 0.5851253E+00, 0.5955415E+00,
     + 0.6053231E+00, 0.6144402E+00, 0.6228659E+00, 0.6305766E+00,
     + 0.6375519E+00, 0.6437751E+00, 0.6492330E+00, 0.6539157E+00,
     + 0.6578171E+00, 0.6609344E+00, 0.6632684E+00, 0.6648231E+00,
     + 0.6656059E+00, 0.6656271E+00, 0.6649001E+00, 0.6634410E+00,
     + 0.6612683E+00, 0.6584031E+00, 0.6548686E+00, 0.6506899E+00,
     + 0.6458935E+00, 0.6405079E+00, 0.6345624E+00, 0.6280874E+00,
     + 0.6211142E+00, 0.6136745E+00, 0.6058004E+00, 0.5975243E+00,
     + 0.5888782E+00, 0.5798943E+00, 0.5706041E+00, 0.5610387E+00,
     + 0.5512287E+00, 0.5412036E+00, 0.5309922E+00, 0.5206223E+00,
     + 0.5101206E+00, 0.4995128E+00, 0.4888233E+00, 0.4780752E+00,
     + 0.4672906E+00, 0.4564902E+00, 0.4456934E+00, 0.4349184E+00,
     + 0.4241820E+00, 0.4135001E+00, 0.4028870E+00, 0.3923560E+00,
     + 0.3819192E+00, 0.3715877E+00, 0.3613713E+00, 0.3512791E+00,
     + 0.3413188E+00, 0.3314975E+00, 0.3218214E+00, 0.3122956E+00,
     + 0.3029247E+00, 0.2937124E+00, 0.2846618E+00, 0.2757754E+00,
     + 0.2670550E+00, 0.2585018E+00, 0.2501168E+00, 0.2419002E+00,
     + 0.2338519E+00, 0.2259716E+00, 0.2182583E+00, 0.2107111E+00,
     + 0.2033284E+00, 0.1961088E+00, 0.1890503E+00, 0.1821509E+00,
     + 0.1754085E+00, 0.1688206E+00, 0.1623848E+00, 0.1560986E+00,
     + 0.1499593E+00, 0.1439641E+00, 0.1381104E+00, 0.1323953E+00,
     + 0.1268161E+00, 0.1213697E+00, 0.1160535E+00, 0.1108645E+00,
     + 0.1058000E+00, 0.1008571E+00, 0.9603301E-01, 0.9132500E-01,
     + 0.8673035E-01, 0.8224636E-01, 0.7787041E-01, 0.7359989E-01,
     + 0.6943226E-01, 0.6536502E-01, 0.6139572E-01, 0.5752197E-01,
     + 0.5374143E-01, 0.5005180E-01, 0.4645086E-01, 0.4293641E-01,
     + 0.3950634E-01, 0.3615856E-01, 0.3289106E-01, 0.2970186E-01,
     + 0.2658905E-01, 0.2355075E-01, 0.2058514E-01, 0.1769046E-01,
     + 0.1486498E-01, 0.1210703E-01, 0.9414958E-02, 0.6787192E-02,
     + 0.4222181E-02, 0.1718422E-02,-0.7255480E-03,-0.3111152E-02,
     +-0.5439776E-02,-0.7712768E-02,-0.9931439E-02,-0.1209707E-01,
     +-0.1421089E-01,-0.1627413E-01,-0.1828795E-01,-0.2025350E-01,
     +-0.2217191E-01,-0.2404425E-01,-0.2587159E-01,-0.2765496E-01,
     +-0.2939536E-01,-0.3109379E-01,-0.3275118E-01,-0.3436847E-01,
     +-0.3594658E-01,-0.3748638E-01,-0.3898874E-01,-0.4045451E-01,
     +-0.4188450E-01,-0.4327951E-01,-0.4464033E-01,-0.4596772E-01,
     +-0.4726243E-01,-0.4852518E-01,-0.4975668E-01,-0.5095763E-01,
     +-0.5212870E-01,-0.5327056E-01,-0.5438385E-01,-0.5546920E-01,
     +-0.5652722E-01,-0.5755853E-01,-0.5856370E-01,-0.5954331E-01,
     +-0.6049793E-01,-0.6142810E-01,-0.6233436E-01,-0.6321723E-01,
     +-0.6407723E-01,-0.6491485E-01,-0.6573059E-01,-0.6652493E-01,
     +-0.6729833E-01,-0.6805125E-01,-0.6878414E-01,-0.6949743E-01,
     +-0.7019156E-01,-0.7086694E-01,-0.7152398E-01,-0.7216308E-01,
     +-0.7278464E-01,-0.7338903E-01,-0.7397663E-01,-0.7454780E-01,
     +-0.7510291E-01,-0.7564231E-01,-0.7616633E-01,-0.7667531E-01,
     +-0.7716959E-01,-0.7764947E-01,-0.7811527E-01,-0.7856731E-01,
     +-0.7900587E-01,-0.7943126E-01,-0.7984376E-01,-0.8024365E-01,
     +-0.8063120E-01,-0.8100669E-01,-0.8137038E-01,-0.8172251E-01,
     +-0.8206336E-01,-0.8239315E-01,-0.8271214E-01,-0.8302055E-01,
     +-0.8331862E-01,-0.8360657E-01,-0.8388463E-01,-0.8415300E-01,
     +-0.8441190E-01,-0.8466153E-01,-0.8490210E-01,-0.8513381E-01,
     +-0.8535684E-01,-0.8557138E-01,-0.8577763E-01,-0.8597576E-01,
     +-0.8616594E-01,-0.8634835E-01,-0.8652316E-01,-0.8669054E-01,
     +-0.8685064E-01,-0.8700363E-01,-0.8714965E-01,-0.8728887E-01,
     +-0.8742143E-01,-0.8754747E-01,-0.8766713E-01,-0.8778056E-01,
     +-0.8788789E-01,-0.8798925E-01,-0.8808478E-01,-0.8817459E-01,
     +-0.8825882E-01,-0.8833759E-01,-0.8841101E-01,-0.8847920E-01,
     +-0.8854227E-01,-0.8860034E-01,-0.8865352E-01,-0.8870191E-01,
     +-0.8874562E-01,-0.8878474E-01,-0.8881939E-01,-0.8884965E-01,
     +-0.8887562E-01,-0.8889740E-01,-0.8891508E-01,-0.8892874E-01,
     +-0.8893847E-01,-0.8894437E-01,-0.8894651E-01,-0.8894498E-01,
     +-0.8893986E-01,-0.8893122E-01,-0.8891914E-01,-0.8890370E-01,
     +-0.8888498E-01,-0.8886303E-01,-0.8883795E-01,-0.8880978E-01,
     +-0.8877861E-01,-0.8874450E-01,-0.8870751E-01,-0.8866770E-01,
     +-0.8862515E-01,-0.8857990E-01,-0.8853202E-01,-0.8848157E-01,
     +-0.8842861E-01,-0.8837318E-01,-0.8831535E-01,-0.8825516E-01,
     +-0.8819268E-01,-0.8812795E-01,-0.8806102E-01,-0.8799194E-01,
     +-0.8792077E-01,-0.8784754E-01,-0.8777230E-01,-0.8769510E-01,
     +-0.8761598E-01,-0.8753499E-01,-0.8745217E-01,-0.8736756E-01,
     +-0.8728120E-01,-0.8719313E-01,-0.8710339E-01,-0.8701201E-01,
     +-0.8691905E-01,-0.8682452E-01,-0.8672847E-01,-0.8663094E-01,
     +-0.8653195E-01,-0.8643155E-01,-0.8632976E-01,-0.8622662E-01,
     +-0.8612215E-01,-0.8601640E-01,-0.8590939E-01,-0.8580115E-01,
     +-0.8569170E-01,-0.8558109E-01,-0.8546933E-01,-0.8535646E-01,
     +-0.8524250E-01,-0.8512747E-01,-0.8501141E-01,-0.8489434E-01,
     +-0.8477628E-01,-0.8465726E-01,-0.8453730E-01,-0.8441642E-01,
     +-0.8429466E-01,-0.8417202E-01,-0.8404854E-01,-0.8392423E-01,
     +-0.8379912E-01,-0.8367322E-01,-0.8354656E-01,-0.8341916E-01,
     +-0.8329103E-01,-0.8316220E-01,-0.8303268E-01,-0.8290250E-01,
     +-0.8277167E-01,-0.8264020E-01,-0.8250812E-01,-0.8237544E-01,
     +-0.8224218E-01,-0.8210836E-01,-0.8197399E-01,-0.8183908E-01,
     +-0.8170366E-01,-0.8156774E-01,-0.8143133E-01,-0.8129444E-01,
     +-0.8115710E-01,-0.8101931E-01,-0.8088109E-01,-0.8074246E-01,
     +-0.8060342E-01,-0.8046399E-01,-0.8032418E-01,-0.8018400E-01,
     +-0.8004347E-01,-0.7990260E-01,-0.7976139E-01,-0.7961987E-01,
     +-0.7947804E-01,-0.7933592E-01,-0.7919350E-01,-0.7905082E-01,
     +-0.7890786E-01,-0.7876466E-01,-0.7862121E-01,-0.7847752E-01,
     +-0.7833361E-01,-0.7818948E-01,-0.7804515E-01,-0.7790062E-01,
     +-0.7775590E-01,-0.7761101E-01,-0.7746594E-01,-0.7732071E-01,
     +-0.7717532E-01,-0.7702979E-01,-0.7688412E-01,-0.7673832E-01,
     +-0.7659240E-01,-0.7644636E-01,-0.7630022E-01,-0.7615397E-01,
     +-0.7600764E-01,-0.7586121E-01,-0.7571471E-01,-0.7556813E-01,
     +-0.7542148E-01,-0.7527478E-01,-0.7512802E-01,-0.7498122E-01,
     +-0.7483437E-01,-0.7468749E-01,-0.7454058E-01,-0.7439365E-01,
     +-0.7424669E-01,-0.7409973E-01,-0.7395276E-01,-0.7380578E-01,
     +-0.7365881E-01,-0.7351185E-01,-0.7336490E-01,-0.7321797E-01,
     +-0.7307106E-01,-0.7292418E-01,-0.7277734E-01,-0.7263053E-01,
     +-0.7248376E-01,-0.7233704E-01,-0.7219037E-01,-0.7204375E-01,
     +-0.7189719E-01,-0.7175070E-01,-0.7160427E-01,-0.7145791E-01,
     +-0.7131163E-01,-0.7116542E-01,-0.7101930E-01,-0.7087326E-01,
     +-0.7072731E-01,-0.7058145E-01,-0.7043569E-01,-0.7029003E-01,
     +-0.7014447E-01,-0.6999902E-01,-0.6985367E-01,-0.6970844E-01,
     +-0.6956332E-01,-0.6941832E-01,-0.6927345E-01,-0.6912869E-01,
     +-0.6898406E-01,-0.6883956E-01,-0.6869520E-01,-0.6855096E-01,
     +-0.6840687E-01,-0.6826291E-01,-0.6811910E-01,-0.6797543E-01,
     +-0.6783190E-01,-0.6768853E-01,-0.6754531E-01,-0.6740224E-01,
     +-0.6725932E-01,-0.6711657E-01,-0.6697397E-01,-0.6683154E-01,
     +-0.6668927E-01,-0.6654717E-01,-0.6640524E-01,-0.6626347E-01,
     +-0.6612188E-01,-0.6598046E-01,-0.6583922E-01,-0.6569815E-01,
     +-0.6555727E-01,-0.6541656E-01,-0.6527603E-01,-0.6513569E-01,
     +-0.6499554E-01,-0.6485557E-01,-0.6471579E-01,-0.6457620E-01,
     +-0.6443680E-01,-0.6429759E-01,-0.6415858E-01,-0.6401976E-01,
     +-0.6388114E-01,-0.6374272E-01,-0.6360449E-01,-0.6346647E-01,
     +-0.6332865E-01,-0.6319103E-01,-0.6305361E-01,-0.6291640E-01,
     +-0.6277939E-01,-0.6264260E-01,-0.6250601E-01,-0.6236963E-01,
     +-0.6223345E-01,-0.6209750E-01,-0.6196175E-01,-0.6182621E-01,
     +-0.6169089E-01,-0.6155578E-01,-0.6142089E-01,-0.6128622E-01,
     +-0.6115176E-01,-0.6101752E-01,-0.6088350E-01,-0.6074970E-01,
     +-0.6061612E-01,-0.6048276E-01,-0.6034962E-01,-0.6021670E-01,
     +-0.6008400E-01,-0.5995153E-01,-0.5981929E-01,-0.5968726E-01,
     +-0.5955546E-01,-0.5942389E-01,-0.5929255E-01,-0.5916143E-01,
     +-0.5903054E-01,-0.5889987E-01,-0.5876944E-01,-0.5863923E-01,
     +-0.5850925E-01,-0.5837951E-01,-0.5824999E-01,-0.5812070E-01,
     +-0.5799164E-01,-0.5786282E-01,-0.5773422E-01,-0.5760586E-01,
     +-0.5747773E-01,-0.5734983E-01,-0.5722217E-01,-0.5709474E-01,
     +-0.5696754E-01,-0.5684057E-01,-0.5671384E-01,-0.5658734E-01,
     +-0.5646108E-01,-0.5633505E-01,-0.5620926E-01,-0.5608370E-01,
     +-0.5595838E-01,-0.5583329E-01,-0.5570844E-01,-0.5558382E-01,
     +-0.5545944E-01,-0.5533530E-01,-0.5521139E-01,-0.5508772E-01,
     +-0.5496428E-01,-0.5484108E-01,-0.5471812E-01,-0.5459539E-01,
     +-0.5447290E-01,-0.5435065E-01,-0.5422864E-01,-0.5410686E-01,
     +-0.5398531E-01,-0.5386401E-01,-0.5374294E-01,-0.5362211E-01,
     +-0.5350151E-01,-0.5338116E-01,-0.5326103E-01,-0.5314115E-01,
     +-0.5302150E-01,-0.5290209E-01,-0.5278292E-01,-0.5266398E-01,
     +-0.5254528E-01,-0.5242682E-01,-0.5230859E-01,-0.5219060E-01,
     +-0.5207285E-01,-0.5195533E-01,-0.5183805E-01,-0.5172100E-01,
     +-0.5160420E-01,-0.5148762E-01,-0.5137129E-01,-0.5125518E-01,
     +-0.5113932E-01,-0.5102369E-01,-0.5090829E-01,-0.5079313E-01,
     +-0.5067821E-01,-0.5056352E-01,-0.5044906E-01,-0.5033484E-01,
     +-0.5022086E-01,-0.5010711E-01,-0.4999359E-01,-0.4988030E-01,
     +-0.4976725E-01,-0.4965444E-01,-0.4954186E-01,-0.4942951E-01,
     +-0.4931739E-01,-0.4920551E-01,-0.4909385E-01,-0.4898243E-01,
     +-0.4887125E-01,-0.4876029E-01,-0.4864957E-01,-0.4853908E-01,
     +-0.4842882E-01,-0.4831879E-01,-0.4820899E-01,-0.4809942E-01,
     +-0.4799008E-01,-0.4788097E-01,-0.4777210E-01,-0.4766345E-01,
     +-0.4755503E-01,-0.4744684E-01,-0.4733888E-01,-0.4723114E-01,
     +-0.4712364E-01,-0.4701636E-01,-0.4690931E-01,-0.4680249E-01,
     +-0.4669590E-01,-0.4658953E-01,-0.4648339E-01,-0.4637747E-01,
     +-0.4627178E-01,-0.4616632E-01,-0.4606108E-01,-0.4595607E-01,
     +-0.4585128E-01,-0.4574672E-01,-0.4564238E-01,-0.4553826E-01,
     +-0.4543437E-01,-0.4533071E-01,-0.4522726E-01,-0.4512404E-01,
     +-0.4502104E-01,-0.4491826E-01,-0.4481570E-01,-0.4471337E-01,
     +-0.4461126E-01,-0.4450936E-01,-0.4440769E-01,-0.4430624E-01,
     +-0.4420501E-01,-0.4410399E-01,-0.4400320E-01,-0.4390262E-01,
     +-0.4380227E-01,-0.4370213E-01,-0.4360221E-01,-0.4350251E-01,
     +-0.4340302E-01,-0.4330375E-01,-0.4320470E-01,-0.4310586E-01,
     +-0.4300724E-01,-0.4290884E-01,-0.4281064E-01,-0.4271267E-01,
     +-0.4261491E-01,-0.4251736E-01,-0.4242003E-01,-0.4232291E-01,
     +-0.4222600E-01,-0.4212931E-01,-0.4203282E-01,-0.4193655E-01,
     +-0.4184050E-01,-0.4174465E-01,-0.4164901E-01,-0.4155359E-01,
     +-0.4145837E-01,-0.4136336E-01,-0.4126857E-01,-0.4117398E-01,
     +-0.4107960E-01,-0.4098543E-01,-0.4089147E-01,-0.4079771E-01,
     +-0.4070417E-01,-0.4061082E-01,-0.4051769E-01,-0.4042476E-01,
     +-0.4033204E-01,-0.4023952E-01,-0.4014721E-01,-0.4005510E-01,
     +-0.3996320E-01,-0.3987150E-01,-0.3978001E-01,-0.3968871E-01,
     +-0.3959762E-01,-0.3950674E-01,-0.3941605E-01,-0.3932557E-01,
     +-0.3923529E-01,-0.3914521E-01,-0.3905533E-01,-0.3896565E-01,
     +-0.3887617E-01,-0.3878689E-01,-0.3869780E-01,-0.3860892E-01,
     +-0.3852024E-01,-0.3843175E-01,-0.3834346E-01,-0.3825537E-01,
     +-0.3816747E-01,-0.3807978E-01,-0.3799227E-01,-0.3790497E-01,
     +-0.3781786E-01,-0.3773094E-01,-0.3764422E-01,-0.3755769E-01,
     +-0.3747136E-01,-0.3738522E-01,-0.3729927E-01,-0.3721352E-01,
     +-0.3712796E-01,-0.3704259E-01,-0.3695741E-01,-0.3687242E-01,
     +-0.3678763E-01,-0.3670302E-01,-0.3661861E-01,-0.3653438E-01,
     +-0.3645034E-01,-0.3636650E-01,-0.3628284E-01,-0.3619937E-01,
     +-0.3611609E-01,-0.3603299E-01,-0.3595008E-01,-0.3586736E-01,
     +-0.3578483E-01,-0.3570248E-01,-0.3562031E-01,-0.3553834E-01,
     +-0.3545654E-01,-0.3537494E-01,-0.3529351E-01,-0.3521227E-01,
     +-0.3513121E-01,-0.3505034E-01,-0.3496965E-01,-0.3488914E-01,
     +-0.3480881E-01,-0.3472867E-01,-0.3464870E-01,-0.3456892E-01,
     +-0.3448931E-01,-0.3440989E-01,-0.3433065E-01,-0.3425158E-01,
     +-0.3417270E-01,-0.3409399E-01,-0.3401546E-01,-0.3393711E-01,
     +-0.3385894E-01,-0.3378094E-01,-0.3370312E-01,-0.3362548E-01,
     +-0.3354801E-01,-0.3347072E-01,-0.3339360E-01,-0.3331666E-01,
     +-0.3323990E-01,-0.3316330E-01,-0.3308688E-01,-0.3301064E-01,
     +-0.3293457E-01,-0.3285867E-01,-0.3278294E-01,-0.3270738E-01,
     +-0.3263200E-01,-0.3255679E-01,-0.3248175E-01,-0.3240688E-01,
     +-0.3233218E-01,-0.3225765E-01,-0.3218328E-01,-0.3210909E-01,
     +-0.3203507E-01,-0.3196121E-01,-0.3188753E-01,-0.3181401E-01,
     +-0.3174066E-01,-0.3166747E-01,-0.3159445E-01,-0.3152160E-01,
     +-0.3144892E-01,-0.3137639E-01,-0.3130404E-01,-0.3123185E-01,
     +-0.3115982E-01,-0.3108796E-01,-0.3101626E-01,-0.3094473E-01,
     +-0.3087336E-01,-0.3080215E-01,-0.3073110E-01,-0.3066022E-01,
     +-0.3058949E-01,-0.3051893E-01,-0.3044853E-01,-0.3037829E-01,
     +-0.3030821E-01,-0.3023829E-01,-0.3016853E-01,-0.3009893E-01,
     +-0.3002949E-01,-0.2996021E-01,-0.2989108E-01,-0.2982212E-01,
     +-0.2975331E-01,-0.2968466E-01,-0.2961616E-01,-0.2954782E-01,
     +-0.2947964E-01,-0.2941161E-01,-0.2934374E-01,-0.2927603E-01,
     +-0.2920846E-01,-0.2914106E-01,-0.2907380E-01,-0.2900671E-01,
     +-0.2893976E-01,-0.2887297E-01,-0.2880633E-01,-0.2873984E-01,
     +-0.2867351E-01,-0.2860732E-01,-0.2854129E-01,-0.2847541E-01,
     +-0.2840968E-01,-0.2834410E-01,-0.2827867E-01,-0.2821340E-01,
     +-0.2814827E-01,-0.2808329E-01,-0.2801845E-01,-0.2795377E-01,
     +-0.2788924E-01,-0.2782485E-01,-0.2776061E-01,-0.2769652E-01,
     +-0.2763257E-01,-0.2756877E-01,-0.2750512E-01,-0.2744161E-01,
     +-0.2737825E-01,-0.2731504E-01,-0.2725196E-01,-0.2718904E-01,
     +-0.2712626E-01,-0.2706362E-01,-0.2700112E-01,-0.2693877E-01,
     +-0.2687656E-01,-0.2681450E-01,-0.2675258E-01,-0.2669079E-01,
     +-0.2662916E-01,-0.2656766E-01,-0.2650630E-01,-0.2644508E-01,
     +-0.2638401E-01,-0.2632307E-01,-0.2626228E-01,-0.2620162E-01,
     +-0.2614111E-01,-0.2608073E-01,-0.2602049E-01,-0.2596039E-01,
     +-0.2590043E-01,-0.2584060E-01,-0.2578091E-01,-0.2572136E-01,
     +-0.2566195E-01,-0.2560267E-01,-0.2554353E-01,-0.2548453E-01,
     +-0.2542566E-01,-0.2536692E-01,-0.2530832E-01,-0.2524986E-01,
     +-0.2519153E-01,-0.2513333E-01,-0.2507527E-01,-0.2501734E-01,
     +-0.2495955E-01,-0.2490188E-01,-0.2484435E-01,-0.2478695E-01,
     +-0.2472969E-01,-0.2467255E-01,-0.2461555E-01,-0.2455868E-01,
     +-0.2450194E-01,-0.2444533E-01,-0.2438885E-01,-0.2433249E-01,
     +-0.2427627E-01,-0.2422018E-01,-0.2416422E-01,-0.2410838E-01,
     +-0.2405268E-01,-0.2399710E-01,-0.2394165E-01,-0.2388633E-01,
     +-0.2383113E-01,-0.2377607E-01,-0.2372113E-01,-0.2366631E-01,
     +-0.2361162E-01,-0.2355706E-01,-0.2350262E-01,-0.2344831E-01,
     +-0.2339413E-01,-0.2334006E-01,-0.2328613E-01,-0.2323231E-01,
     +-0.2317862E-01,-0.2312506E-01,-0.2307161E-01,-0.2301830E-01,
     +-0.2296510E-01,-0.2291202E-01,-0.2285907E-01,-0.2280624E-01,
     +-0.2275353E-01,-0.2270095E-01,-0.2264848E-01,-0.2259614E-01,
     +-0.2254391E-01,-0.2249181E-01,-0.2243982E-01,-0.2238796E-01,
     +-0.2233621E-01,-0.2228459E-01,-0.2223308E-01,-0.2218169E-01,
     +-0.2213042E-01,-0.2207927E-01,-0.2202824E-01,-0.2197732E-01,
     +-0.2192652E-01,-0.2187584E-01,-0.2182528E-01,-0.2177483E-01,
     +-0.2172450E-01,-0.2167428E-01,-0.2162418E-01,-0.2157419E-01,
     +-0.2152432E-01,-0.2147457E-01,-0.2142493E-01,-0.2137540E-01,
     +-0.2132599E-01,-0.2127669E-01,-0.2122751E-01,-0.2117844E-01,
     +-0.2112948E-01,-0.2108064E-01,-0.2103191E-01,-0.2098329E-01,
     +-0.2093478E-01,-0.2088638E-01,-0.2083810E-01,-0.2078993E-01,
     +-0.2074186E-01,-0.2069391E-01,-0.2064607E-01,-0.2059834E-01,
     +-0.2055072E-01,-0.2050321E-01,-0.2045581E-01,-0.2040852E-01,
     +-0.2036134E-01,-0.2031427E-01,-0.2026730E-01,-0.2022044E-01,
     +-0.2017370E-01,-0.2012706E-01,-0.2008052E-01,-0.2003410E-01,
     +-0.1998778E-01,-0.1994157E-01,-0.1989546E-01,-0.1984946E-01,
     +-0.1980357E-01,-0.1975779E-01,-0.1971210E-01,-0.1966653E-01,
     +-0.1962106E-01,-0.1957569E-01,-0.1953043E-01,-0.1948528E-01,
     +-0.1944022E-01,-0.1939528E-01,-0.1935043E-01,-0.1930569E-01,
     +-0.1926105E-01,-0.1921652E-01,-0.1917209E-01,-0.1912776E-01,
     +-0.1908353E-01,-0.1903940E-01,-0.1899538E-01,-0.1895146E-01,
     +-0.1890764E-01,-0.1886392E-01,-0.1882030E-01,-0.1877678E-01,
     +-0.1873337E-01,-0.1869005E-01,-0.1864683E-01,-0.1860372E-01,
     +-0.1856070E-01,-0.1851778E-01,-0.1847496E-01,-0.1843224E-01,
     +-0.1838962E-01,-0.1834710E-01,-0.1830467E-01,-0.1826234E-01,
     +-0.1822012E-01,-0.1817798E-01,-0.1813595E-01,-0.1809401E-01,
     +-0.1805217E-01,-0.1801043E-01,-0.1796878E-01,-0.1792723E-01,
     +-0.1788577E-01,-0.1784441E-01,-0.1780315E-01,-0.1776198E-01,
     +-0.1772090E-01,-0.1767992E-01,-0.1763904E-01,-0.1759825E-01,
     +-0.1755755E-01,-0.1751695E-01,-0.1747644E-01,-0.1743603E-01,
     +-0.1739571E-01,-0.1735548E-01,-0.1731534E-01,-0.1727530E-01,
     +-0.1723535E-01,-0.1719549E-01,-0.1715573E-01,-0.1711605E-01,
     +-0.1707647E-01,-0.1703698E-01,-0.1699758E-01,-0.1695827E-01,
     +-0.1691905E-01,-0.1687993E-01,-0.1684089E-01,-0.1680194E-01,
     +-0.1676309E-01,-0.1672432E-01,-0.1668564E-01,-0.1664705E-01,
     +-0.1660856E-01,-0.1657015E-01,-0.1653182E-01,-0.1649359E-01,
     +-0.1645545E-01,-0.1641739E-01,-0.1637942E-01,-0.1634154E-01,
     +-0.1630375E-01,-0.1626604E-01,-0.1622842E-01,-0.1619089E-01,
     +-0.1615345E-01,-0.1611609E-01,-0.1607882E-01,-0.1604163E-01,
     +-0.1600453E-01,-0.1596752E-01,-0.1593059E-01,-0.1589374E-01,
     +-0.1585698E-01,-0.1582031E-01,-0.1578372E-01,-0.1574722E-01,
     +-0.1571080E-01,-0.1567446E-01,-0.1563821E-01,-0.1560204E-01,
     +-0.1556596E-01,-0.1552996E-01,-0.1549404E-01,-0.1545820E-01,
     +-0.1542245E-01,-0.1538678E-01,-0.1535119E-01,-0.1531569E-01,
     +-0.1528027E-01,-0.1524493E-01,-0.1520967E-01,-0.1517449E-01,
     +-0.1513939E-01,-0.1510438E-01,-0.1506944E-01,-0.1503459E-01,
     +-0.1499982E-01,-0.1496512E-01,-0.1493051E-01,-0.1489598E-01,
     +-0.1486152E-01,-0.1482715E-01,-0.1479286E-01,-0.1475864E-01,
     +-0.1472451E-01,-0.1469045E-01,-0.1465647E-01,-0.1462257E-01,
     +-0.1458875E-01,-0.1455501E-01,-0.1452135E-01,-0.1448776E-01,
     +-0.1445425E-01,-0.1442082E-01,-0.1438746E-01,-0.1435419E-01,
     +-0.1432098E-01,-0.1428786E-01,-0.1425481E-01,-0.1422184E-01,
     +-0.1418895E-01,-0.1415613E-01,-0.1412339E-01,-0.1409072E-01,
     +-0.1405813E-01,-0.1402561E-01,-0.1399317E-01,-0.1396080E-01,
     +-0.1392851E-01,-0.1389630E-01,-0.1386415E-01,-0.1383209E-01,
     +-0.1380009E-01,-0.1376817E-01,-0.1373633E-01,-0.1370456E-01,
     +-0.1367286E-01,-0.1364123E-01,-0.1360968E-01,-0.1357820E-01,
     +-0.1354679E-01,-0.1351546E-01,-0.1348420E-01,-0.1345301E-01,
     +-0.1342189E-01,-0.1339084E-01,-0.1335987E-01,-0.1332897E-01,
     +-0.1329814E-01,-0.1326738E-01,-0.1323669E-01,-0.1320607E-01,
     +-0.1317553E-01,-0.1314505E-01,-0.1311464E-01,-0.1308431E-01,
     +-0.1305404E-01,-0.1302385E-01,-0.1299372E-01,-0.1296367E-01,
     +-0.1293368E-01,-0.1290377E-01,-0.1287392E-01,-0.1284414E-01,
     +-0.1281443E-01,-0.1278479E-01,-0.1275522E-01,-0.1272571E-01,
     +-0.1269628E-01,-0.1266691E-01,-0.1263761E-01,-0.1260838E-01,
     +-0.1257921E-01,-0.1255012E-01,-0.1252109E-01,-0.1249212E-01,
     +-0.1246323E-01,-0.1243440E-01,-0.1240564E-01,-0.1237694E-01,
     +-0.1234831E-01,-0.1231975E-01,-0.1229125E-01,-0.1226282E-01,
     +-0.1223445E-01,-0.1220615E-01,-0.1217792E-01,-0.1214975E-01,
     +-0.1212165E-01,-0.1209361E-01,-0.1206563E-01,-0.1203772E-01,
     +-0.1200988E-01,-0.1198210E-01,-0.1195438E-01,-0.1192673E-01,
     +-0.1189914E-01,-0.1187161E-01,-0.1184415E-01,-0.1181676E-01,
     +-0.1178942E-01,-0.1176215E-01,-0.1173494E-01,-0.1170780E-01,
     +-0.1168072E-01,-0.1165370E-01,-0.1162674E-01,-0.1159984E-01,
     +-0.1157301E-01,-0.1154624E-01,-0.1151953E-01,-0.1149289E-01,
     +-0.1146630E-01,-0.1143978E-01,-0.1141332E-01,-0.1138691E-01,
     +-0.1136057E-01,-0.1133429E-01,-0.1130808E-01,-0.1128192E-01,
     +-0.1125582E-01,-0.1122978E-01,-0.1120381E-01,-0.1117789E-01,
     +-0.1115203E-01,-0.1112624E-01,-0.1110050E-01,-0.1107482E-01,
     +-0.1104920E-01,-0.1102364E-01,-0.1099814E-01,-0.1097270E-01,
     +-0.1094732E-01,-0.1092200E-01,-0.1089673E-01,-0.1087153E-01,
     +-0.1084638E-01,-0.1082129E-01,-0.1079626E-01,-0.1077128E-01,
     +-0.1074636E-01,-0.1072151E-01,-0.1069670E-01,-0.1067196E-01,
     +-0.1064727E-01,-0.1062264E-01,-0.1059807E-01,-0.1057356E-01,
     +-0.1054910E-01,-0.1052469E-01,-0.1050035E-01,-0.1047606E-01,
     +-0.1045182E-01,-0.1042765E-01,-0.1040353E-01,-0.1037946E-01,
     +-0.1035545E-01,-0.1033149E-01,-0.1030760E-01,-0.1028375E-01,
     +-0.1025996E-01,-0.1023623E-01,-0.1021255E-01,-0.1018893E-01,
     +-0.1016536E-01,-0.1014184E-01,-0.1011838E-01,-0.1009497E-01,
     +-0.1007162E-01,-0.1004832E-01,-0.1002508E-01,-0.1000189E-01,
     +-0.9978751E-02,-0.9955668E-02,-0.9932638E-02,-0.9909661E-02,
     +-0.9886737E-02,-0.9863866E-02,-0.9841049E-02,-0.9818284E-02,
     +-0.9795571E-02,-0.9772912E-02,-0.9750304E-02,-0.9727749E-02,
     +-0.9705246E-02,-0.9682795E-02,-0.9660396E-02,-0.9638049E-02,
     +-0.9615754E-02,-0.9593510E-02,-0.9571317E-02,-0.9549176E-02,
     +-0.9527086E-02,-0.9505047E-02,-0.9483060E-02,-0.9461123E-02,
     +-0.9439236E-02,-0.9417401E-02,-0.9395616E-02,-0.9373881E-02,
     +-0.9352196E-02,-0.9330562E-02,-0.9308978E-02,-0.9287444E-02,
     +-0.9265959E-02,-0.9244524E-02,-0.9223139E-02,-0.9201803E-02,
     +-0.9180517E-02,-0.9159279E-02,-0.9138091E-02,-0.9116952E-02,
     +-0.9095862E-02,-0.9074821E-02,-0.9053828E-02,-0.9032884E-02,
     +-0.9011988E-02,-0.8991140E-02,-0.8970341E-02,-0.8949590E-02,
     +-0.8928887E-02,-0.8908232E-02,-0.8887624E-02,-0.8867065E-02,
     +-0.8846552E-02,-0.8826088E-02,-0.8805670E-02,-0.8785300E-02,
     +-0.8764977E-02,-0.8744701E-02,-0.8724472E-02,-0.8704289E-02,
     +-0.8684153E-02,-0.8664064E-02,-0.8644022E-02,-0.8624025E-02,
     +-0.8604075E-02,-0.8584171E-02,-0.8564314E-02,-0.8544502E-02,
     +-0.8524736E-02,-0.8505015E-02,-0.8485340E-02,-0.8465711E-02,
     +-0.8446127E-02,-0.8426589E-02,-0.8407095E-02,-0.8387647E-02,
     +-0.8368244E-02,-0.8348885E-02,-0.8329572E-02,-0.8310303E-02,
     +-0.8291078E-02,-0.8271899E-02,-0.8252763E-02,-0.8233672E-02,
     +-0.8214625E-02,-0.8195622E-02,-0.8176662E-02,-0.8157747E-02,
     +-0.8138876E-02,-0.8120048E-02,-0.8101263E-02,-0.8082523E-02,
     +-0.8063825E-02,-0.8045171E-02,-0.8026560E-02,-0.8007992E-02,
     +-0.7989467E-02,-0.7970984E-02,-0.7952545E-02,-0.7934148E-02,
     +-0.7915794E-02,-0.7897482E-02,-0.7879212E-02,-0.7860985E-02,
     +-0.7842800E-02,-0.7824657E-02,-0.7806556E-02,-0.7788497E-02,
     +-0.7770480E-02,-0.7752504E-02,-0.7734570E-02,-0.7716677E-02,
     +-0.7698826E-02,-0.7681016E-02,-0.7663247E-02,-0.7645520E-02,
     +-0.7627833E-02,-0.7610187E-02,-0.7592582E-02,-0.7575018E-02,
     +-0.7557495E-02,-0.7540012E-02,-0.7522569E-02,-0.7505167E-02,
     +-0.7487805E-02,-0.7470483E-02,-0.7453201E-02,-0.7435959E-02,
     +-0.7418758E-02,-0.7401595E-02,-0.7384473E-02,-0.7367390E-02,
     +-0.7350347E-02,-0.7333343E-02,-0.7316379E-02,-0.7299453E-02,
     +-0.7282567E-02,-0.7265720E-02,-0.7248912E-02,-0.7232143E-02,
     +-0.7215412E-02,-0.7198721E-02,-0.7182067E-02,-0.7165453E-02,
     +-0.7148877E-02,-0.7132339E-02,-0.7115839E-02,-0.7099378E-02,
     +-0.7082955E-02,-0.7066569E-02,-0.7050222E-02,-0.7033912E-02,
     +-0.7017640E-02,-0.7001406E-02,-0.6985209E-02,-0.6969050E-02,
     +-0.6952928E-02,-0.6936844E-02,-0.6920796E-02,-0.6904786E-02,
     +-0.6888813E-02,-0.6872877E-02,-0.6856977E-02,-0.6841115E-02,
     +-0.6825289E-02,-0.6809500E-02,-0.6793747E-02,-0.6778031E-02,
     +-0.6762351E-02,-0.6746707E-02,-0.6731100E-02,-0.6715528E-02,
     +-0.6699993E-02,-0.6684493E-02,-0.6669030E-02,-0.6653602E-02,
     +-0.6638210E-02,-0.6622853E-02,-0.6607532E-02,-0.6592247E-02,
     +-0.6576997E-02,-0.6561782E-02,-0.6546602E-02,-0.6531457E-02,
     +-0.6516348E-02,-0.6501273E-02,-0.6486233E-02,-0.6471228E-02,
     +-0.6456258E-02,-0.6441322E-02,-0.6426421E-02,-0.6411555E-02,
     +-0.6396722E-02,-0.6381925E-02,-0.6367161E-02,-0.6352432E-02/
c     ------------------------------------------------------------
c     derivative of d-state wave function
      data uadp/0.1739843E-03,
     + 0.1739843E-03, 0.6421329E-03, 0.1398629E-02, 0.2426725E-02,
     + 0.3712279E-02, 0.5243528E-02, 0.7010854E-02, 0.9006557E-02,
     + 0.1122462E-01, 0.1366046E-01, 0.1631075E-01, 0.1917313E-01,
     + 0.2224603E-01, 0.2552848E-01, 0.2901982E-01, 0.3271960E-01,
     + 0.3662735E-01, 0.4074237E-01, 0.4506362E-01, 0.4958954E-01,
     + 0.5431787E-01, 0.5924557E-01, 0.6436868E-01, 0.6968221E-01,
     + 0.7518004E-01, 0.8085488E-01, 0.8669817E-01, 0.9270007E-01,
     + 0.9884943E-01, 0.1051338E+00, 0.1115393E+00, 0.1180511E+00,
     + 0.1246526E+00, 0.1313266E+00, 0.1380545E+00, 0.1448167E+00,
     + 0.1515929E+00, 0.1583619E+00, 0.1651019E+00, 0.1717907E+00,
     + 0.1784057E+00, 0.1849242E+00, 0.1913234E+00, 0.1975808E+00,
     + 0.2036742E+00, 0.2095817E+00, 0.2152823E+00, 0.2207558E+00,
     + 0.2259827E+00, 0.2309449E+00, 0.2356252E+00, 0.2400079E+00,
     + 0.2440788E+00, 0.2478250E+00, 0.2512353E+00, 0.2543000E+00,
     + 0.2570112E+00, 0.2593627E+00, 0.2613498E+00, 0.2629696E+00,
     + 0.2642210E+00, 0.2651041E+00, 0.2656211E+00, 0.2657754E+00,
     + 0.2655719E+00, 0.2650168E+00, 0.2641177E+00, 0.2628834E+00,
     + 0.2613236E+00, 0.2594492E+00, 0.2572718E+00, 0.2548040E+00,
     + 0.2520588E+00, 0.2490499E+00, 0.2457915E+00, 0.2422981E+00,
     + 0.2385844E+00, 0.2346653E+00, 0.2305559E+00, 0.2262712E+00,
     + 0.2218261E+00, 0.2172353E+00, 0.2125134E+00, 0.2076746E+00,
     + 0.2027330E+00, 0.1977021E+00, 0.1925950E+00, 0.1874243E+00,
     + 0.1822024E+00, 0.1769409E+00, 0.1716509E+00, 0.1663432E+00,
     + 0.1610277E+00, 0.1557141E+00, 0.1504112E+00, 0.1451275E+00,
     + 0.1398708E+00, 0.1346485E+00, 0.1294674E+00, 0.1243337E+00,
     + 0.1192531E+00, 0.1142310E+00, 0.1092721E+00, 0.1043807E+00,
     + 0.9956078E-01, 0.9481580E-01, 0.9014882E-01, 0.8556255E-01,
     + 0.8105933E-01, 0.7664115E-01, 0.7230970E-01, 0.6806637E-01,
     + 0.6391224E-01, 0.5984817E-01, 0.5587475E-01, 0.5199234E-01,
     + 0.4820110E-01, 0.4450100E-01, 0.4089182E-01, 0.3737319E-01,
     + 0.3394457E-01, 0.3060531E-01, 0.2735461E-01, 0.2419159E-01,
     + 0.2111524E-01, 0.1812448E-01, 0.1521815E-01, 0.1239500E-01,
     + 0.9653732E-02, 0.6993009E-02, 0.4411429E-02, 0.1907558E-02,
     +-0.5200681E-03,-0.2872942E-02,-0.5152576E-02,-0.7360498E-02,
     +-0.9498243E-02,-0.1156736E-01,-0.1356939E-01,-0.1550588E-01,
     +-0.1737836E-01,-0.1918838E-01,-0.2093746E-01,-0.2262710E-01,
     +-0.2425880E-01,-0.2583404E-01,-0.2735428E-01,-0.2882095E-01,
     +-0.3023549E-01,-0.3159928E-01,-0.3291369E-01,-0.3418009E-01,
     +-0.3539979E-01,-0.3657409E-01,-0.3770427E-01,-0.3879157E-01,
     +-0.3983722E-01,-0.4084241E-01,-0.4180830E-01,-0.4273604E-01,
     +-0.4362675E-01,-0.4448150E-01,-0.4530136E-01,-0.4608737E-01,
     +-0.4684054E-01,-0.4756184E-01,-0.4825223E-01,-0.4891266E-01,
     +-0.4954402E-01,-0.5014720E-01,-0.5072306E-01,-0.5127244E-01,
     +-0.5179614E-01,-0.5229496E-01,-0.5276967E-01,-0.5322101E-01,
     +-0.5364970E-01,-0.5405646E-01,-0.5444197E-01,-0.5480689E-01,
     +-0.5515186E-01,-0.5547752E-01,-0.5578447E-01,-0.5607329E-01,
     +-0.5634457E-01,-0.5659885E-01,-0.5683668E-01,-0.5705857E-01,
     +-0.5726505E-01,-0.5745659E-01,-0.5763367E-01,-0.5779676E-01,
     +-0.5794630E-01,-0.5808273E-01,-0.5820647E-01,-0.5831794E-01,
     +-0.5841752E-01,-0.5850560E-01,-0.5858255E-01,-0.5864874E-01,
     +-0.5870451E-01,-0.5875021E-01,-0.5878617E-01,-0.5881269E-01,
     +-0.5883011E-01,-0.5883870E-01,-0.5883877E-01,-0.5883059E-01,
     +-0.5881445E-01,-0.5879060E-01,-0.5875930E-01,-0.5872081E-01,
     +-0.5867536E-01,-0.5862319E-01,-0.5856452E-01,-0.5849958E-01,
     +-0.5842858E-01,-0.5835172E-01,-0.5826922E-01,-0.5818126E-01,
     +-0.5808803E-01,-0.5798971E-01,-0.5788650E-01,-0.5777854E-01,
     +-0.5766603E-01,-0.5754910E-01,-0.5742794E-01,-0.5730268E-01,
     +-0.5717348E-01,-0.5704047E-01,-0.5690381E-01,-0.5676362E-01,
     +-0.5662004E-01,-0.5647319E-01,-0.5632320E-01,-0.5617019E-01,
     +-0.5601427E-01,-0.5585556E-01,-0.5569417E-01,-0.5553020E-01,
     +-0.5536376E-01,-0.5519495E-01,-0.5502387E-01,-0.5485061E-01,
     +-0.5467526E-01,-0.5449791E-01,-0.5431865E-01,-0.5413756E-01,
     +-0.5395473E-01,-0.5377023E-01,-0.5358414E-01,-0.5339654E-01,
     +-0.5320750E-01,-0.5301708E-01,-0.5282536E-01,-0.5263240E-01,
     +-0.5243827E-01,-0.5224303E-01,-0.5204675E-01,-0.5184947E-01,
     +-0.5165125E-01,-0.5145216E-01,-0.5125225E-01,-0.5105157E-01,
     +-0.5085016E-01,-0.5064809E-01,-0.5044539E-01,-0.5024212E-01,
     +-0.5003831E-01,-0.4983402E-01,-0.4962928E-01,-0.4942414E-01,
     +-0.4921863E-01,-0.4901280E-01,-0.4880669E-01,-0.4860032E-01,
     +-0.4839374E-01,-0.4818698E-01,-0.4798007E-01,-0.4777305E-01,
     +-0.4756595E-01,-0.4735879E-01,-0.4715161E-01,-0.4694444E-01,
     +-0.4673730E-01,-0.4653023E-01,-0.4632324E-01,-0.4611636E-01,
     +-0.4590962E-01,-0.4570304E-01,-0.4549665E-01,-0.4529046E-01,
     +-0.4508450E-01,-0.4487879E-01,-0.4467335E-01,-0.4446820E-01,
     +-0.4426336E-01,-0.4405885E-01,-0.4385468E-01,-0.4365088E-01,
     +-0.4344745E-01,-0.4324442E-01,-0.4304180E-01,-0.4283960E-01,
     +-0.4263785E-01,-0.4243655E-01,-0.4223572E-01,-0.4203537E-01,
     +-0.4183552E-01,-0.4163617E-01,-0.4143735E-01,-0.4123905E-01,
     +-0.4104130E-01,-0.4084409E-01,-0.4064745E-01,-0.4045138E-01,
     +-0.4025590E-01,-0.4006100E-01,-0.3986670E-01,-0.3967302E-01,
     +-0.3947994E-01,-0.3928750E-01,-0.3909568E-01,-0.3890451E-01,
     +-0.3871398E-01,-0.3852410E-01,-0.3833488E-01,-0.3814633E-01,
     +-0.3795845E-01,-0.3777124E-01,-0.3758472E-01,-0.3739888E-01,
     +-0.3721373E-01,-0.3702929E-01,-0.3684554E-01,-0.3666249E-01,
     +-0.3648016E-01,-0.3629854E-01,-0.3611763E-01,-0.3593745E-01,
     +-0.3575798E-01,-0.3557925E-01,-0.3540124E-01,-0.3522396E-01,
     +-0.3504741E-01,-0.3487160E-01,-0.3469653E-01,-0.3452220E-01,
     +-0.3434861E-01,-0.3417576E-01,-0.3400365E-01,-0.3383230E-01,
     +-0.3366168E-01,-0.3349182E-01,-0.3332271E-01,-0.3315434E-01,
     +-0.3298673E-01,-0.3281987E-01,-0.3265376E-01,-0.3248840E-01,
     +-0.3232379E-01,-0.3215994E-01,-0.3199684E-01,-0.3183449E-01,
     +-0.3167289E-01,-0.3151204E-01,-0.3135195E-01,-0.3119260E-01,
     +-0.3103401E-01,-0.3087617E-01,-0.3071907E-01,-0.3056272E-01,
     +-0.3040712E-01,-0.3025227E-01,-0.3009816E-01,-0.2994480E-01,
     +-0.2979218E-01,-0.2964030E-01,-0.2948916E-01,-0.2933876E-01,
     +-0.2918910E-01,-0.2904018E-01,-0.2889199E-01,-0.2874453E-01,
     +-0.2859781E-01,-0.2845181E-01,-0.2830654E-01,-0.2816200E-01,
     +-0.2801819E-01,-0.2787509E-01,-0.2773272E-01,-0.2759107E-01,
     +-0.2745013E-01,-0.2730991E-01,-0.2717040E-01,-0.2703160E-01,
     +-0.2689351E-01,-0.2675612E-01,-0.2661944E-01,-0.2648346E-01,
     +-0.2634818E-01,-0.2621360E-01,-0.2607972E-01,-0.2594652E-01,
     +-0.2581402E-01,-0.2568220E-01,-0.2555107E-01,-0.2542062E-01,
     +-0.2529085E-01,-0.2516176E-01,-0.2503334E-01,-0.2490560E-01,
     +-0.2477852E-01,-0.2465212E-01,-0.2452638E-01,-0.2440130E-01,
     +-0.2427688E-01,-0.2415312E-01,-0.2403001E-01,-0.2390755E-01,
     +-0.2378575E-01,-0.2366459E-01,-0.2354407E-01,-0.2342420E-01,
     +-0.2330496E-01,-0.2318636E-01,-0.2306839E-01,-0.2295105E-01,
     +-0.2283434E-01,-0.2271825E-01,-0.2260279E-01,-0.2248794E-01,
     +-0.2237371E-01,-0.2226009E-01,-0.2214709E-01,-0.2203469E-01,
     +-0.2192290E-01,-0.2181171E-01,-0.2170112E-01,-0.2159112E-01,
     +-0.2148172E-01,-0.2137291E-01,-0.2126469E-01,-0.2115706E-01,
     +-0.2105000E-01,-0.2094353E-01,-0.2083763E-01,-0.2073231E-01,
     +-0.2062756E-01,-0.2052338E-01,-0.2041977E-01,-0.2031672E-01,
     +-0.2021423E-01,-0.2011229E-01,-0.2001091E-01,-0.1991009E-01,
     +-0.1980981E-01,-0.1971008E-01,-0.1961089E-01,-0.1951225E-01,
     +-0.1941414E-01,-0.1931657E-01,-0.1921953E-01,-0.1912302E-01,
     +-0.1902704E-01,-0.1893159E-01,-0.1883665E-01,-0.1874224E-01,
     +-0.1864834E-01,-0.1855495E-01,-0.1846208E-01,-0.1836971E-01,
     +-0.1827785E-01,-0.1818650E-01,-0.1809564E-01,-0.1800529E-01,
     +-0.1791542E-01,-0.1782605E-01,-0.1773718E-01,-0.1764878E-01,
     +-0.1756088E-01,-0.1747345E-01,-0.1738650E-01,-0.1730004E-01,
     +-0.1721404E-01,-0.1712852E-01,-0.1704347E-01,-0.1695888E-01,
     +-0.1687476E-01,-0.1679110E-01,-0.1670790E-01,-0.1662515E-01,
     +-0.1654286E-01,-0.1646103E-01,-0.1637964E-01,-0.1629870E-01,
     +-0.1621820E-01,-0.1613814E-01,-0.1605853E-01,-0.1597935E-01,
     +-0.1590060E-01,-0.1582229E-01,-0.1574441E-01,-0.1566695E-01,
     +-0.1558992E-01,-0.1551331E-01,-0.1543713E-01,-0.1536136E-01,
     +-0.1528600E-01,-0.1521106E-01,-0.1513653E-01,-0.1506241E-01,
     +-0.1498870E-01,-0.1491539E-01,-0.1484248E-01,-0.1476997E-01,
     +-0.1469786E-01,-0.1462614E-01,-0.1455482E-01,-0.1448389E-01,
     +-0.1441334E-01,-0.1434318E-01,-0.1427341E-01,-0.1420402E-01,
     +-0.1413500E-01,-0.1406637E-01,-0.1399811E-01,-0.1393022E-01,
     +-0.1386271E-01,-0.1379556E-01,-0.1372878E-01,-0.1366237E-01,
     +-0.1359631E-01,-0.1353062E-01,-0.1346529E-01,-0.1340031E-01,
     +-0.1333569E-01,-0.1327142E-01,-0.1320750E-01,-0.1314393E-01,
     +-0.1308071E-01,-0.1301783E-01,-0.1295529E-01,-0.1289310E-01,
     +-0.1283124E-01,-0.1276972E-01,-0.1270853E-01,-0.1264768E-01,
     +-0.1258715E-01,-0.1252696E-01,-0.1246709E-01,-0.1240755E-01,
     +-0.1234833E-01,-0.1228943E-01,-0.1223086E-01,-0.1217260E-01,
     +-0.1211465E-01,-0.1205702E-01,-0.1199970E-01,-0.1194270E-01,
     +-0.1188600E-01,-0.1182960E-01,-0.1177352E-01,-0.1171773E-01,
     +-0.1166225E-01,-0.1160707E-01,-0.1155218E-01,-0.1149759E-01,
     +-0.1144330E-01,-0.1138930E-01,-0.1133559E-01,-0.1128217E-01,
     +-0.1122904E-01,-0.1117619E-01,-0.1112363E-01,-0.1107135E-01,
     +-0.1101935E-01,-0.1096763E-01,-0.1091619E-01,-0.1086503E-01,
     +-0.1081414E-01,-0.1076352E-01,-0.1071318E-01,-0.1066311E-01,
     +-0.1061330E-01,-0.1056376E-01,-0.1051449E-01,-0.1046547E-01,
     +-0.1041673E-01,-0.1036824E-01,-0.1032001E-01,-0.1027204E-01,
     +-0.1022432E-01,-0.1017686E-01,-0.1012966E-01,-0.1008270E-01,
     +-0.1003600E-01,-0.9989541E-02,-0.9943333E-02,-0.9897370E-02,
     +-0.9851651E-02,-0.9806176E-02,-0.9760942E-02,-0.9715948E-02,
     +-0.9671194E-02,-0.9626676E-02,-0.9582395E-02,-0.9538349E-02,
     +-0.9494536E-02,-0.9450954E-02,-0.9407604E-02,-0.9364483E-02,
     +-0.9321589E-02,-0.9278923E-02,-0.9236481E-02,-0.9194264E-02,
     +-0.9152269E-02,-0.9110496E-02,-0.9068943E-02,-0.9027608E-02,
     +-0.8986492E-02,-0.8945591E-02,-0.8904906E-02,-0.8864435E-02,
     +-0.8824176E-02,-0.8784128E-02,-0.8744291E-02,-0.8704662E-02,
     +-0.8665241E-02,-0.8626027E-02,-0.8587018E-02,-0.8548213E-02,
     +-0.8509611E-02,-0.8471210E-02,-0.8433010E-02,-0.8395010E-02,
     +-0.8357208E-02,-0.8319603E-02,-0.8282194E-02,-0.8244980E-02,
     +-0.8207959E-02,-0.8171131E-02,-0.8134495E-02,-0.8098049E-02,
     +-0.8061792E-02,-0.8025724E-02,-0.7989842E-02,-0.7954147E-02,
     +-0.7918636E-02,-0.7883310E-02,-0.7848166E-02,-0.7813204E-02,
     +-0.7778423E-02,-0.7743821E-02,-0.7709398E-02,-0.7675153E-02,
     +-0.7641084E-02,-0.7607191E-02,-0.7573472E-02,-0.7539927E-02,
     +-0.7506554E-02,-0.7473354E-02,-0.7440323E-02,-0.7407462E-02,
     +-0.7374770E-02,-0.7342246E-02,-0.7309888E-02,-0.7277696E-02,
     +-0.7245669E-02,-0.7213805E-02,-0.7182105E-02,-0.7150566E-02,
     +-0.7119189E-02,-0.7087972E-02,-0.7056914E-02,-0.7026014E-02,
     +-0.6995271E-02,-0.6964686E-02,-0.6934255E-02,-0.6903980E-02,
     +-0.6873858E-02,-0.6843889E-02,-0.6814072E-02,-0.6784407E-02,
     +-0.6754892E-02,-0.6725526E-02,-0.6696309E-02,-0.6667240E-02,
     +-0.6638318E-02,-0.6609542E-02,-0.6580911E-02,-0.6552424E-02,
     +-0.6524082E-02,-0.6495882E-02,-0.6467824E-02,-0.6439907E-02,
     +-0.6412131E-02,-0.6384494E-02,-0.6356997E-02,-0.6329637E-02,
     +-0.6302414E-02,-0.6275328E-02,-0.6248378E-02,-0.6221562E-02,
     +-0.6194881E-02,-0.6168333E-02,-0.6141918E-02,-0.6115635E-02,
     +-0.6089483E-02,-0.6063461E-02,-0.6037569E-02,-0.6011806E-02,
     +-0.5986172E-02,-0.5960664E-02,-0.5935284E-02,-0.5910029E-02,
     +-0.5884900E-02,-0.5859896E-02,-0.5835016E-02,-0.5810258E-02,
     +-0.5785624E-02,-0.5761111E-02,-0.5736720E-02,-0.5712449E-02,
     +-0.5688298E-02,-0.5664266E-02,-0.5640353E-02,-0.5616557E-02,
     +-0.5592879E-02,-0.5569317E-02,-0.5545871E-02,-0.5522541E-02,
     +-0.5499325E-02,-0.5476222E-02,-0.5453234E-02,-0.5430358E-02,
     +-0.5407594E-02,-0.5384941E-02,-0.5362399E-02,-0.5339968E-02,
     +-0.5317646E-02,-0.5295433E-02,-0.5273328E-02,-0.5251331E-02,
     +-0.5229442E-02,-0.5207659E-02,-0.5185981E-02,-0.5164410E-02,
     +-0.5142943E-02,-0.5121580E-02,-0.5100321E-02,-0.5079165E-02,
     +-0.5058112E-02,-0.5037161E-02,-0.5016311E-02,-0.4995561E-02,
     +-0.4974913E-02,-0.4954363E-02,-0.4933913E-02,-0.4913562E-02,
     +-0.4893309E-02,-0.4873153E-02,-0.4853094E-02,-0.4833132E-02,
     +-0.4813266E-02,-0.4793495E-02,-0.4773819E-02,-0.4754237E-02,
     +-0.4734750E-02,-0.4715355E-02,-0.4696054E-02,-0.4676844E-02,
     +-0.4657727E-02,-0.4638701E-02,-0.4619766E-02,-0.4600921E-02,
     +-0.4582166E-02,-0.4563500E-02,-0.4544923E-02,-0.4526434E-02,
     +-0.4508034E-02,-0.4489721E-02,-0.4471494E-02,-0.4453354E-02,
     +-0.4435301E-02,-0.4417332E-02,-0.4399449E-02,-0.4381651E-02,
     +-0.4363936E-02,-0.4346305E-02,-0.4328758E-02,-0.4311293E-02,
     +-0.4293911E-02,-0.4276610E-02,-0.4259391E-02,-0.4242253E-02,
     +-0.4225196E-02,-0.4208218E-02,-0.4191321E-02,-0.4174502E-02,
     +-0.4157763E-02,-0.4141102E-02,-0.4124519E-02,-0.4108013E-02,
     +-0.4091585E-02,-0.4075234E-02,-0.4058958E-02,-0.4042759E-02,
     +-0.4026636E-02,-0.4010587E-02,-0.3994613E-02,-0.3978714E-02,
     +-0.3962888E-02,-0.3947136E-02,-0.3931457E-02,-0.3915851E-02,
     +-0.3900317E-02,-0.3884856E-02,-0.3869465E-02,-0.3854146E-02,
     +-0.3838898E-02,-0.3823721E-02,-0.3808613E-02,-0.3793575E-02,
     +-0.3778607E-02,-0.3763707E-02,-0.3748876E-02,-0.3734114E-02,
     +-0.3719419E-02,-0.3704792E-02,-0.3690231E-02,-0.3675738E-02,
     +-0.3661311E-02,-0.3646951E-02,-0.3632656E-02,-0.3618426E-02,
     +-0.3604262E-02,-0.3590162E-02,-0.3576127E-02,-0.3562155E-02,
     +-0.3548248E-02,-0.3534403E-02,-0.3520622E-02,-0.3506904E-02,
     +-0.3493248E-02,-0.3479654E-02,-0.3466122E-02,-0.3452651E-02,
     +-0.3439241E-02,-0.3425892E-02,-0.3412603E-02,-0.3399375E-02,
     +-0.3386206E-02,-0.3373097E-02,-0.3360047E-02,-0.3347056E-02,
     +-0.3334124E-02,-0.3321250E-02,-0.3308434E-02,-0.3295675E-02,
     +-0.3282974E-02,-0.3270330E-02,-0.3257743E-02,-0.3245212E-02,
     +-0.3232738E-02,-0.3220319E-02,-0.3207956E-02,-0.3195649E-02,
     +-0.3183396E-02,-0.3171198E-02,-0.3159055E-02,-0.3146966E-02,
     +-0.3134930E-02,-0.3122949E-02,-0.3111021E-02,-0.3099145E-02,
     +-0.3087323E-02,-0.3075553E-02,-0.3063836E-02,-0.3052170E-02,
     +-0.3040557E-02,-0.3028994E-02,-0.3017483E-02,-0.3006023E-02,
     +-0.2994614E-02,-0.2983255E-02,-0.2971946E-02,-0.2960687E-02,
     +-0.2949478E-02,-0.2938318E-02,-0.2927207E-02,-0.2916146E-02,
     +-0.2905133E-02,-0.2894168E-02,-0.2883251E-02,-0.2872383E-02,
     +-0.2861562E-02,-0.2850789E-02,-0.2840062E-02,-0.2829383E-02,
     +-0.2818751E-02,-0.2808164E-02,-0.2797625E-02,-0.2787131E-02,
     +-0.2776683E-02,-0.2766280E-02,-0.2755923E-02,-0.2745611E-02,
     +-0.2735344E-02,-0.2725122E-02,-0.2714943E-02,-0.2704809E-02,
     +-0.2694719E-02,-0.2684673E-02,-0.2674670E-02,-0.2664711E-02,
     +-0.2654795E-02,-0.2644921E-02,-0.2635091E-02,-0.2625302E-02,
     +-0.2615556E-02,-0.2605852E-02,-0.2596190E-02,-0.2586570E-02,
     +-0.2576990E-02,-0.2567453E-02,-0.2557955E-02,-0.2548499E-02,
     +-0.2539084E-02,-0.2529708E-02,-0.2520373E-02,-0.2511078E-02,
     +-0.2501823E-02,-0.2492608E-02,-0.2483431E-02,-0.2474295E-02,
     +-0.2465196E-02,-0.2456137E-02,-0.2447117E-02,-0.2438135E-02,
     +-0.2429191E-02,-0.2420285E-02,-0.2411417E-02,-0.2402587E-02,
     +-0.2393794E-02,-0.2385038E-02,-0.2376320E-02,-0.2367639E-02,
     +-0.2358994E-02,-0.2350386E-02,-0.2341814E-02,-0.2333279E-02,
     +-0.2324779E-02,-0.2316316E-02,-0.2307888E-02,-0.2299496E-02,
     +-0.2291139E-02,-0.2282817E-02,-0.2274531E-02,-0.2266279E-02,
     +-0.2258061E-02,-0.2249878E-02,-0.2241730E-02,-0.2233616E-02,
     +-0.2225536E-02,-0.2217489E-02,-0.2209476E-02,-0.2201497E-02,
     +-0.2193551E-02,-0.2185639E-02,-0.2177759E-02,-0.2169912E-02,
     +-0.2162098E-02,-0.2154317E-02,-0.2146567E-02,-0.2138850E-02,
     +-0.2131166E-02,-0.2123512E-02,-0.2115891E-02,-0.2108302E-02,
     +-0.2100744E-02,-0.2093217E-02,-0.2085721E-02,-0.2078256E-02,
     +-0.2070822E-02,-0.2063419E-02,-0.2056047E-02,-0.2048705E-02,
     +-0.2041393E-02,-0.2034111E-02,-0.2026859E-02,-0.2019638E-02,
     +-0.2012446E-02,-0.2005283E-02,-0.1998150E-02,-0.1991046E-02,
     +-0.1983972E-02,-0.1976926E-02,-0.1969909E-02,-0.1962921E-02,
     +-0.1955961E-02,-0.1949031E-02,-0.1942128E-02,-0.1935253E-02,
     +-0.1928407E-02,-0.1921588E-02,-0.1914797E-02,-0.1908034E-02,
     +-0.1901299E-02,-0.1894591E-02,-0.1887910E-02,-0.1881256E-02,
     +-0.1874629E-02,-0.1868030E-02,-0.1861457E-02,-0.1854910E-02,
     +-0.1848390E-02,-0.1841897E-02,-0.1835429E-02,-0.1828988E-02,
     +-0.1822573E-02,-0.1816184E-02,-0.1809820E-02,-0.1803483E-02,
     +-0.1797170E-02,-0.1790883E-02,-0.1784622E-02,-0.1778386E-02,
     +-0.1772174E-02,-0.1765988E-02,-0.1759827E-02,-0.1753690E-02,
     +-0.1747578E-02,-0.1741490E-02,-0.1735427E-02,-0.1729388E-02,
     +-0.1723373E-02,-0.1717382E-02,-0.1711415E-02,-0.1705472E-02,
     +-0.1699552E-02,-0.1693656E-02,-0.1687784E-02,-0.1681935E-02,
     +-0.1676109E-02,-0.1670307E-02,-0.1664527E-02,-0.1658771E-02,
     +-0.1653037E-02,-0.1647326E-02,-0.1641638E-02,-0.1635972E-02,
     +-0.1630329E-02,-0.1624708E-02,-0.1619109E-02,-0.1613532E-02,
     +-0.1607977E-02,-0.1602445E-02,-0.1596934E-02,-0.1591445E-02,
     +-0.1585977E-02,-0.1580531E-02,-0.1575107E-02,-0.1569704E-02,
     +-0.1564322E-02,-0.1558961E-02,-0.1553621E-02,-0.1548302E-02,
     +-0.1543004E-02,-0.1537727E-02,-0.1532470E-02,-0.1527234E-02,
     +-0.1522019E-02,-0.1516823E-02,-0.1511649E-02,-0.1506494E-02,
     +-0.1501359E-02,-0.1496245E-02,-0.1491150E-02,-0.1486075E-02,
     +-0.1481020E-02,-0.1475985E-02,-0.1470969E-02,-0.1465973E-02,
     +-0.1460996E-02,-0.1456038E-02,-0.1451099E-02,-0.1446180E-02,
     +-0.1441280E-02,-0.1436398E-02,-0.1431536E-02,-0.1426692E-02,
     +-0.1421867E-02,-0.1417061E-02,-0.1412273E-02,-0.1407504E-02,
     +-0.1402753E-02,-0.1398020E-02,-0.1393306E-02,-0.1388610E-02,
     +-0.1383932E-02,-0.1379272E-02,-0.1374630E-02,-0.1370005E-02,
     +-0.1365399E-02,-0.1360809E-02,-0.1356238E-02,-0.1351684E-02,
     +-0.1347147E-02,-0.1342628E-02,-0.1338126E-02,-0.1333641E-02,
     +-0.1329174E-02,-0.1324723E-02,-0.1320289E-02,-0.1315873E-02,
     +-0.1311473E-02,-0.1307090E-02,-0.1302723E-02,-0.1298374E-02,
     +-0.1294040E-02,-0.1289723E-02,-0.1285423E-02,-0.1281139E-02,
     +-0.1276871E-02,-0.1272619E-02,-0.1268383E-02,-0.1264164E-02,
     +-0.1259960E-02,-0.1255772E-02,-0.1251601E-02,-0.1247445E-02,
     +-0.1243304E-02,-0.1239179E-02,-0.1235070E-02,-0.1230976E-02,
     +-0.1226897E-02,-0.1222834E-02,-0.1218786E-02,-0.1214754E-02,
     +-0.1210736E-02,-0.1206734E-02,-0.1202747E-02,-0.1198774E-02,
     +-0.1194817E-02,-0.1190874E-02,-0.1186947E-02,-0.1183033E-02,
     +-0.1179135E-02,-0.1175250E-02,-0.1171381E-02,-0.1167526E-02,
     +-0.1163685E-02,-0.1159858E-02,-0.1156046E-02,-0.1152248E-02,
     +-0.1148464E-02,-0.1144695E-02,-0.1140939E-02,-0.1137197E-02,
     +-0.1133469E-02,-0.1129756E-02,-0.1126055E-02,-0.1122369E-02,
     +-0.1118696E-02,-0.1115037E-02,-0.1111391E-02,-0.1107759E-02,
     +-0.1104140E-02,-0.1100534E-02,-0.1096942E-02,-0.1093363E-02,
     +-0.1089797E-02,-0.1086245E-02,-0.1082706E-02,-0.1079179E-02,
     +-0.1075666E-02,-0.1072165E-02,-0.1068678E-02,-0.1065203E-02,
     +-0.1061741E-02,-0.1058292E-02,-0.1054855E-02,-0.1051431E-02,
     +-0.1048020E-02,-0.1044621E-02,-0.1041234E-02,-0.1037860E-02,
     +-0.1034498E-02,-0.1031149E-02,-0.1027812E-02,-0.1024487E-02,
     +-0.1021174E-02,-0.1017873E-02,-0.1014585E-02,-0.1011308E-02,
     +-0.1008043E-02,-0.1004790E-02,-0.1001549E-02,-0.9983194E-03,
     +-0.9951020E-03,-0.9918960E-03,-0.9887015E-03,-0.9855187E-03,
     +-0.9823475E-03,-0.9791875E-03,-0.9760393E-03,-0.9729024E-03,
     +-0.9697769E-03,-0.9666628E-03,-0.9635598E-03,-0.9604679E-03,
     +-0.9573872E-03,-0.9543177E-03,-0.9512593E-03,-0.9482121E-03,
     +-0.9451756E-03,-0.9421500E-03,-0.9391355E-03,-0.9361315E-03,
     +-0.9331388E-03,-0.9301561E-03,-0.9271845E-03,-0.9242238E-03,
     +-0.9212739E-03,-0.9183342E-03,-0.9154051E-03,-0.9124864E-03,
     +-0.9095785E-03,-0.9066808E-03,-0.9037934E-03,-0.9009163E-03,
     +-0.8980496E-03,-0.8951934E-03,-0.8923472E-03,-0.8895112E-03,
     +-0.8866855E-03,-0.8838693E-03,-0.8810633E-03,-0.8782674E-03,
     +-0.8754817E-03,-0.8727060E-03,-0.8699397E-03,-0.8671835E-03,
     +-0.8644371E-03,-0.8616999E-03,-0.8589730E-03,-0.8562552E-03,
     +-0.8535475E-03,-0.8508497E-03,-0.8481610E-03,-0.8454818E-03,
     +-0.8428125E-03,-0.8401528E-03,-0.8375017E-03,-0.8348603E-03,
     +-0.8322283E-03,-0.8296055E-03,-0.8269918E-03,-0.8243868E-03,
     +-0.8217918E-03,-0.8192057E-03,-0.8166287E-03,-0.8140609E-03,
     +-0.8115017E-03,-0.8089517E-03,-0.8064108E-03,-0.8038786E-03,
     +-0.8013558E-03,-0.7988420E-03,-0.7963368E-03,-0.7938399E-03,
     +-0.7913518E-03,-0.7888729E-03,-0.7864023E-03,-0.7839404E-03,
     +-0.7814873E-03,-0.7790429E-03,-0.7766066E-03,-0.7741790E-03,
     +-0.7717595E-03,-0.7693488E-03,-0.7669462E-03,-0.7645518E-03,
     +-0.7621662E-03,-0.7597889E-03,-0.7574198E-03,-0.7550587E-03,
     +-0.7527059E-03,-0.7503618E-03,-0.7480255E-03,-0.7456972E-03,
     +-0.7433770E-03,-0.7410649E-03,-0.7387607E-03,-0.7364645E-03,
     +-0.7341760E-03,-0.7318952E-03,-0.7296228E-03,-0.7273580E-03,
     +-0.7251014E-03,-0.7228524E-03,-0.7206111E-03,-0.7183777E-03,
     +-0.7161522E-03,-0.7139341E-03,-0.7117235E-03,-0.7095204E-03,
     +-0.7073251E-03,-0.7051374E-03,-0.7029571E-03,-0.7007845E-03,
     +-0.6986192E-03,-0.6964610E-03,-0.6943110E-03,-0.6921680E-03,
     +-0.6900323E-03,-0.6879040E-03,-0.6857830E-03,-0.6836689E-03,
     +-0.6815627E-03,-0.6794634E-03,-0.6773721E-03,-0.6752873E-03,
     +-0.6732099E-03,-0.6711397E-03,-0.6690759E-03,-0.6670193E-03,
     +-0.6649695E-03,-0.6629269E-03,-0.6608911E-03,-0.6588621E-03,
     +-0.6568397E-03,-0.6548248E-03,-0.6528165E-03,-0.6508145E-03,
     +-0.6488202E-03,-0.6468324E-03,-0.6448512E-03,-0.6428769E-03,
     +-0.6409092E-03,-0.6389487E-03,-0.6369945E-03,-0.6350470E-03,
     +-0.6331059E-03,-0.6311713E-03,-0.6292435E-03,-0.6273218E-03,
     +-0.6254072E-03,-0.6234993E-03,-0.6215974E-03,-0.6197018E-03,
     +-0.6178126E-03,-0.6159298E-03,-0.6140532E-03,-0.6121831E-03,
     +-0.6103191E-03,-0.6084618E-03,-0.6066107E-03,-0.6047655E-03,
     +-0.6029264E-03,-0.6010939E-03,-0.5992672E-03,-0.5974467E-03,
     +-0.5956327E-03,-0.5938245E-03,-0.5920221E-03,-0.5902253E-03,
     +-0.5884351E-03,-0.5866511E-03,-0.5848729E-03,-0.5831011E-03,
     +-0.5813356E-03,-0.5795756E-03,-0.5778212E-03,-0.5760725E-03,
     +-0.5743302E-03,-0.5725932E-03,-0.5708623E-03,-0.5691371E-03,
     +-0.5674179E-03,-0.5657043E-03,-0.5639960E-03,-0.5622941E-03,
     +-0.5605974E-03,-0.5589058E-03,-0.5572205E-03,-0.5555405E-03,
     +-0.5538665E-03,-0.5521971E-03,-0.5505334E-03,-0.5488758E-03,
     +-0.5472240E-03,-0.5455782E-03,-0.5439365E-03,-0.5423010E-03,
     +-0.5406705E-03,-0.5390449E-03,-0.5374256E-03,-0.5358109E-03,
     +-0.5342017E-03,-0.5325983E-03,-0.5309996E-03,-0.5294063E-03,
     +-0.5278181E-03,-0.5262365E-03,-0.5246597E-03,-0.5230875E-03,
     +-0.5215202E-03,-0.5199585E-03,-0.5184012E-03,-0.5168499E-03,
     +-0.5153034E-03,-0.5137622E-03,-0.5122257E-03,-0.5106950E-03,
     +-0.5091690E-03,-0.5076480E-03,-0.5061315E-03,-0.5046204E-03,
     +-0.5031140E-03,-0.5016131E-03,-0.5001169E-03,-0.4986261E-03,
     +-0.4971397E-03,-0.4956578E-03,-0.4941807E-03,-0.4927087E-03,
     +-0.4912411E-03,-0.4897788E-03,-0.4883211E-03,-0.4868690E-03,
     +-0.4854206E-03,-0.4839762E-03,-0.4825383E-03,-0.4811047E-03,
     +-0.4796756E-03,-0.4782512E-03,-0.4768315E-03,-0.4754157E-03,
     +-0.4740039E-03,-0.4725967E-03,-0.4711940E-03,-0.4697970E-03,
     +-0.4684047E-03,-0.4670166E-03,-0.4656327E-03,-0.4642533E-03,
     +-0.4628787E-03,-0.4615085E-03,-0.4601421E-03,-0.4587810E-03,
     +-0.4574244E-03,-0.4560721E-03,-0.4547236E-03,-0.4533790E-03,
     +-0.4520393E-03,-0.4507038E-03,-0.4493726E-03,-0.4480453E-03,
     +-0.4467234E-03,-0.4454049E-03,-0.4440913E-03,-0.4427819E-03,
     +-0.4414773E-03,-0.4401762E-03,-0.4388790E-03,-0.4375870E-03,
     +-0.4362983E-03,-0.4350139E-03,-0.4337336E-03,-0.4324572E-03,
     +-0.4311850E-03,-0.4299171E-03,-0.4286534E-03,-0.4273936E-03,
     +-0.4261375E-03,-0.4248857E-03,-0.4236378E-03,-0.4223936E-03,
     +-0.4211535E-03,-0.4199176E-03,-0.4186854E-03,-0.4174572E-03,
     +-0.4162333E-03,-0.4150126E-03,-0.4137957E-03,-0.4125828E-03,
     +-0.4113745E-03,-0.4101702E-03,-0.4089697E-03,-0.4077730E-03,
     +-0.4065796E-03,-0.4053908E-03,-0.4042051E-03,-0.4030234E-03,
     +-0.4018450E-03,-0.4006710E-03,-0.3995011E-03,-0.3983342E-03/
c     -------------------------------------------------------------------
      r=h
      twoh=2.d0*h
      cons=hc*hc/2.d0/rmu
      do i=2,1500
       ds=(uasp(i+1)-uasp(i-1))/twoh
       dd=(uadp(i+1)-uadp(i-1))/twoh
       vavs(i)=ed*uavs(i)+cons*(ds)
       vavd(i)=ed*uavd(i)+cons*(dd-6.d0*uavd(i)/r/r)
       r=r+h
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine saxon(r,rt,deet,g)
      implicit real*8(a-h,o-z)
      g=1.d0/(1.0+exp((r-rt)*deet))
      return
      end
c-----------------------------------------------------------------------
      subroutine dsaxon(r,rt,deet,dg)
      implicit real*8(a-h,o-z)
      gg=exp((r-rt)*deet)
      g=(1.d0+gg)
      dg=4.0*gg/g/g
      return
      end
c-----------------------------------------------------------------------
      subroutine thomas(r,rt,deet,dg)
      implicit real*8(a-h,o-z)
      dg=0.0
      if(r.eq.0.0) return
      gg=exp((r-rt)*deet)
      g=1+gg
      dg=2.0*gg*deet/r/g/g
      return
      end
