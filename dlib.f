      subroutine dlib(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi) 

************************************************************************
************************************************************************
* Stefan Kretzer (kretzer@pa.msu.edu) Nov 2001:                        *
*                                                                      *
* Library of FF parametrizations updated to be posted on               *
* http://www.pv.infn.it/~radici/FFdatabase/                            *
* site maintained by Marco Radici (radici@pv.infn.it)                  *
* and Rainer Jakob (rainer@theorie.physik.uni-wuppertal.de)            *
* within the EU network:                                               *
* "Hadronic Physics with High Energy Electromagnetic Probes"           *
*                                                                      *
* Also posted on that site are the individual FF sets included in this *
* library as distributed by the corresponding authors.                 *
************************************************************************
* K, KKP and BFGW fragmentation functions for                          *
* iparton -> ihadron                                                   *
*                                                                      *
* Subroutine "dlib" returns D(z,Q2) [NOT z*D(z,Q2)] as an array        *
* "dff(iparton)" of dimension dff(-5:5)                                *
*                                                                      *
* iparton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar      *
*                                                                      *
* [Q2] = GeV^2                                                         *
*                                                                      *
* ffset = 1,2,3 means K, KKP, BFGW                                     *
************************************************************************
*        ALEXEI PROKUDIN prokudin@to.infn.it                           *
*        making the same interface FF fDSS as for DLIB                 *
************************************************************************
* ffset = 4 means                                                      *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                     *
*  D.de Florian, R.Sassot, M.Stratmann   hep-ph/0703242)               *
*    Phys.Rev.D.75:114010,2007                                         *
*                                                                      *
************************************************************************
*        ALEXEI PROKUDIN prokudin@to.infn.it                           *
*        making the same interface FF AKK as for DLIB                  *
*                    26/03/2008                                        *
* ffset = 5 means                                                      *
c----------------------------------------------------------------------
*
c     AKK ROUTINES 2008
*
c----------------------------------------------------------------------*
*
*      SUBROUTINE AKK(IH,Z,Q,DH)*
*
* AKK Update: Improvements from New Theoretical Input and Experimental Data.
* S. Albino, B.A. Kniehl, G. Kramer . Mar 2008. 
* e-Print: arXiv:0803.2768 [hep-ph]
* Details are to be found in the corresponding references:             *
* S.Kretzer, Phys.Rev.D62, 054001 (2000)                               *
* B.A.Kniehl, G.Kramer, B.Potter, Nucl.Phys.B582, 514 (2000)           *
* L.Bourhis, M.Fontannaz, J.P.Guillet, M.Werlen,Eur.Phys.J.C19,89(2001)*
*                                                                      *
* fforder=0,1 is LO, NLO(MSbar)                                        *
* no LO set for BFGW                                                   *
*                                                                      *
* ihadron=1,2,3,4,5 is pi,K,h,p,n                                      *
* no pi,K sets for BFGW                                                *
* no p,n sets for K, BFGW                                              *
* no n set for fDSS                                                    * AP
*                                                                      *
* icharge=0,1,2,3 is 0,+,-,+&-                                         *
* note: 3 = +&- = charge sum (NOT average) in this library             *
* for the neutral (icharge=0) particles K^0, n:                        *
* icp = 1,2,3 chooses between particle, anti-particle or sum of both   *
* icp is inactive for pi,h,p and charged Kaons                         *
*                                                                      *
* Following KKP, FFs into {pi^0;(anti-)K^0} and into {pi^+/-;K^+/-}    *
* are related by isospin.                                              *
*                                                                      *
* ipi = 1,2,3 is a flag for BFGW; inactive for K, KKP:                 *
* ipi = 1: best fit     (formula (8))  of BFGW                         *
* ipi = 2: large Ng set (formula (9))  of BFGW                         *
* ipi = 3: low Ng set   (formula (10)) of BFGW                         *
************************************************************************
* If options (ihadron, icharge etc.) are chosen which do not exist for *
* a given parametrization (ffset) or if a flag is chosen outside its   *
* range as defined above then the code returns a corresponding warning *
* and stops.                                                           *
************************************************************************
************************************************************************

      implicit double precision ( a - z )
      dimension dff(-5:5)  
      dimension uqff(2), dqff(2), sqff(2), cqff(2), bqff(2)
      dimension dh(0:10)
      integer ffset, fforder, iparton, ipi, ipisk, ifini, old1, 
     >        iset, ihadron, icharge, old2, old3, ih, isu2, icp
      integer old1_deflorian, old2_deflorian, old3_deflorian ! AP
      integer IC,IO ! AP
      real*8 BB,CB ! AP
      real*8 Q ! AP
      character*12 readme
      COMMON / FRAGINI / IFINI
      COMMON / FRAGINI_DEFLORIAN / IFINI_DEFLORIAN
      common / order / old1
      common / old / old2, old3
      common / order_deflorian / old1_deflorian ! AP
      common / old_deflorian / old2_deflorian, old3_deflorian !AP
      data readme/ '; see MANUAL' /
      ipisk = ipi  
      
      if ( .not. ( (fforder .eq. 0).or.(fforder .eq. 1) ) ) then
      write(6,*) 'choose fforder=0,1 for LO,NLO'//readme
      stop
      endif      

      if ( (ihadron .gt. 3 ) .and. (.not. (ffset .eq. 2)) ) then
      write(6,*) 'no p,n FFs (ihadron=4,5) for K, BFGW (ffset=1,3)'
     > //readme
      stop
      endif

      if ( ihadron .lt. 1 .or. ihadron .gt. 5 ) then 
      write(6,*)
     > 'set ihadron=1,2,3,4,5 for pi,K,h,p,n'//readme
      stop
      endif    

      if ( (icharge .eq. 0 ) .and. 
     >   (  ihadron .eq. 3 .or. ihadron .eq. 4 ) ) then
       write(6,*) 'incl hadrons (ihadron=3) or protons (ihadron=4)'
     > //' must be charged (icharge=1,2,3)'//readme
      stop
      endif 

      if ( icharge .lt. 0 .or. icharge .gt. 3 ) then
      write(6,*) 'set icharge=0,1,2,3 for 0,+,-,+&-'//readme
      stop
      endif    

      if (
     > ( (ihadron .eq. 2 .or. ihadron .eq. 5) .and. icharge .eq. 0 )
     >  .and. ( icp .lt. 1 .or. icp .gt. 3 )  ) then
      write(6,*) 'for neutral (icharge=0) K,n (ihadron=2,5) choose'
     >//' particle, anti-particle or sum (icp=1,2,3)'//readme 
      stop
      endif

      if ( ffset .eq. 1 ) then

      norm = 1.d0
      isu2 = 0
 
      if ( .not. ( fforder  .eq. old1 ) .or.
     >     .not. ( ihadron  .eq. old2 ) .or.
     >     .not. ( icharge  .eq. old3 )
     >   ) ifini = 0
      old1 = fforder
      old2 = ihadron
      old3 = icharge

      iset = 2*ihadron + fforder -1
      if ( ihadron .eq. 1 .and. icharge .eq. 0 ) then
      norm = 0.5d0
      icharge = 3
      endif
      if ( ihadron .eq. 2 .and. icharge .eq. 0 ) then
      icharge = icp
      isu2 = 1
      endif
 
      call PKHFF(ISET,ICHARGE,Z,Q2,uqff,dqff,sqff,cqff,bqff,gff)

      dff(-5     ) = bqff(2) * norm 
      dff(-4     ) = cqff(2) * norm
      dff(-3     ) = sqff(2) * norm
      dff(-2+isu2) = dqff(2) * norm
      dff(-1-isu2) = uqff(2) * norm
      dff( 0     ) = gff     * norm
      dff( 1+isu2) = uqff(1) * norm
      dff( 2-isu2) = dqff(1) * norm
      dff( 3     ) = sqff(1) * norm
      dff( 4     ) = cqff(1) * norm
      dff( 5     ) = bqff(1) * norm
c      print*,icp,isu2
      return

      elseif ( ffset .eq. 2 ) then

      norm = 2.d0

      qs = dsqrt(Q2)
      if ( icharge .eq. 3 ) then
      if ( ihadron .eq. 5 ) then
      write(6,*)
     > 'set icharge=0 for neutrons (ihadron=5)'//readme
      stop
      endif
      ih = ihadron
      if ( ihadron .eq. 3) then
      ih = 7
      norm = 1.d0
      endif
      elseif ( icharge .eq. 0 ) then
       if ( ihadron .eq. 1 ) then
       ih = 5
       norm = 1.d0
       elseif ( ihadron .eq. 2 ) then
       if ( .not. icp .eq. 3 ) then 
       write(6,*)
     > 'only sum of particle & anti-particle (icp=3)'
     > //'for KKP K^0 set'//readme
       stop
       endif 
       ih = 3
       elseif ( ihadron .eq. 5 ) then
       if ( .not. icp .eq. 3 ) then
       write(6,*)
     > 'only sum of particle & anti-particle (icp=3)'
     > //'for KKP neutron set'//readme
       stop
       endif 
       ih = 6
       else
       write(6,*) 
     >'only neutral (icharge=0) pions (ihadron=1), Kaons (ihadron=2)'
     > //' and neutrons (ihadron=5)'//readme
       stop
       endif
      else
      write(6,*) 
     >'icharge must be 0 (neutral) or 3 (charge sum) for KKP (ffset=2)'
     >//readme
      stop
      endif

      call kkp(ih,fforder,z,qs,dh)
      dff(-5) = dh(10) * norm
      dff(-4) = dh( 8) * norm
      dff(-3) = dh( 6) * norm
      dff(-2) = dh( 4) * norm
      dff(-1) = dh( 2) * norm
      dff( 0) = dh( 0) * norm
      dff( 1) = dh( 1) * norm
      dff( 2) = dh( 3) * norm
      dff( 3) = dh( 5) * norm
      dff( 4) = dh( 7) * norm
      dff( 5) = dh( 9) * norm
      return

      elseif ( ffset .eq. 3) then

      if (fforder .eq. 0) then
      write(6,*) 'only NLO (fforder=1) set for BFGW (ffset=3)'//readme
      stop
      endif

      if (.not. ihadron .eq. 3) then
      write(6,*) 'only incl hadron set (ihadron=3) for BFGW (ffset=3)'
     > //readme
      stop
      endif      

      if (.not. icharge .eq. 3) then
      write(6,*) 'only charge sum (icharge=3) for BFGW (ffset=3)'
     > //readme
      stop
      endif

      if ( ipisk .lt. 1 .or. ipisk .gt. 3 ) then
      write(6,*)
     > 'set ipi=1,2,3 for central,large,small Ng set of BFGW (ffset=3)'
     > //readme
      stop
      endif

      call fonfra(z,ipisk,Q2,xdup,xdubp,xddp,xddbp,xdsp,xdcp
     > ,xdbp,xdbbp,xdgp)
      dff(-5) = xdbbp / z
      dff(-4) = xdcp  / z
      dff(-3) = xdsp  / z
      dff(-2) = xddbp / z
      dff(-1) = xdubp / z
      dff( 0) = xdgp  / z
      dff( 1) = xdup  / z
      dff( 2) = xddp  / z
      dff( 3) = xdsp  / z
      dff( 4) = xdcp  / z
      dff( 5) = xdbp  / z
      return

c     Alexei Prokudin 18/07/2007
      elseif ( ffset .eq. 4) then


      if ( .not. ( fforder  .eq. old1_deflorian ) .or.
     >     .not. ( ihadron  .eq. old2_deflorian ) .or.
     >     .not. ( icharge  .eq. old3_deflorian )
     >   ) ifini_deflorian = 0
      old1_deflorian = fforder
      old2_deflorian = ihadron
      old3_deflorian = icharge

      norm = 1.d0 

c* ihadron=1,2,3,4,5 is pi,K,h,p,n                                  *
c*  IH = hadron type    1: PION                                     *
c*                      2: KAON                                     *
c*                      3: PROTON                                   *
c*                      4: CHARGED HADRONS                          *
      if ( ihadron .eq. 1 ) then
         IH = 1
      elseif ( ihadron .eq. 2 ) then
         IH = 2
      elseif ( ihadron .eq. 3 ) then
         IH = 4
      elseif ( ihadron .eq. 4 ) then
         IH = 3
      elseif ( ihadron .eq. 5 ) then
         IH = 4
      write(6,*) 'only pi,K,h,p (ihadron=1,2,3,4) for fDSS (ffset=4)'
     > //readme
      stop
      endif      

c* icharge=0,1,2,3 is 0,+,-,+&-                                     *
c* note: 3 = +&- = charge sum (NOT average) in this library         *
c*  IC = Hadron Charge  0: 0 (as average of + and -)                *
c*                      1: +                                        *
c*                     -1: -                                        *
      if (icharge .eq. 0) then
         IC = 0
      elseif (icharge .eq. 1) then
         IC = 1
      elseif (icharge .eq. 2) then
         IC = -1
      elseif (icharge .eq. 3) then
         IC = 0
         norm = 2.d0 ! Transform from charge average to charge sum ....
      endif

c* fforder=0,1 is LO, NLO(MSbar)                                    *
c*  IO= Order           0: LO                                       *
c*                      1: NLO                                      *
      IO = fforder

      call fDSS (IH,IC,IO, z, Q2, U, UB, D, DB, S, SB, C, B, GL)
      dff(-5) = norm * B  / z
      dff(-4) = norm * C  / z
      dff(-3) = norm * SB / z
      dff(-2) = norm * DB / z
      dff(-1) = norm * UB / z
      dff( 0) = norm * GL / z
      dff( 1) = norm * U  / z
      dff( 2) = norm * D  / z
      dff( 3) = norm * S  / z
      dff( 4) = norm * C  / z
      dff( 5) = norm * B  / z
      return

c     Alexei Prokudin 26/03/2008
      elseif ( ffset .eq. 5) then ! AKK08

      norm = 1.d0 

      Q = DSQRT(Q2)

C IH hadron species
C***************************************************************
C*************** CHARGE-SIGN UNIDENTIFIED FFS ******************
C***************************************************************
C  1  pi^+ + pi^-
C  2  K^+ + K^-
C  3  p + anti-p
C  4  K0short
C  5  Lambda+anti-Lambda
C  6  pi^0, calculated from (pi^+ + pi^-)/2
C  7  K0short, calculated from (K^+ + K^-)/2 with u<->d
C  8  n + anti-n, calculated from p + anti-p with u<->d
C  9  h^+ + h^- (= pi^+ + pi^- + K^+ + K^- + p + anti-p)
C***************************************************************
C************** CHARGE-SIGN ASYMMETRY FFS **********************
C***************************************************************
C 10  h^+ - h^- (= pi^+ - pi^- + K^+ - K^- + p - anti-p)
C 11  pi^+ - pi^-
C 12  K^+ - K^-
C 13  p - anti-p
C ALEXEI PROKUDIN 26/03/2008
C 14  pi^+
C 15  pi^-
C 16  K^+
C 17  K^-
C 18  h^+
c 19  h^-
C END  ALEXEI PROKUDIN 26/03/2008
C DH(i): fragmentation function of parton i
C 0  1    2    3    4    5    6    7    8    9   10
C g  u  u-bar  d  d-bar  s  s-bar  c  c-bar  b  b-bar
C Output:

      if ( ihadron .eq. 1 ) then ! pions

         if (icharge .eq. 0) then ! 0

            IH = 11
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)

         elseif (icharge .eq. 1) then ! +
            IH = 14
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         elseif (icharge .eq. 2) then ! -
            IC = 15
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         elseif (icharge .eq. 3) then ! pi^+ + pi^-
            IH = 1
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         endif


      elseif ( ihadron .eq. 2 ) then ! K

         if (icharge .eq. 0) then ! 0
            IH = 4 ! K0short
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         elseif (icharge .eq. 1) then ! +
            IH = 16

            IH = 2
            CALL AKK(IH,z,Q,DH)

C DH(i): fragmentation function of parton i
C 0  1    2    3    4    5    6    7    8    9   10
C g  u  u-bar  d  d-bar  s  s-bar  c  c-bar  b  b-bar
C Output:

C DH(i): fragmentation function of parton i
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)

            IH = 12
            CALL AKK(IH,z,Q,DH)

            GL= (GL + DH(0) )/2.d0        
            U = (U  + DH(1) )/2.d0        
            UB= (UB + DH(2) )/2.d0        
            D = (D  + DH(3) )/2.d0        
            DB= (DB + DH(4) )/2.d0        
            S = (S  + DH(5) )/2.d0        
            SB= (SB + DH(6) )/2.d0        
            C = (C  + DH(7) )/2.d0        
            CB= (CB + DH(8) )/2.d0        
            B = (B  + DH(9) )/2.d0        
            BB= (BB + DH(10))/2.d0  ! PI^+


         elseif (icharge .eq. 2) then ! -
            IH = 17
             CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
        elseif (icharge .eq. 3) then ! h+ + h-
            IH = 12
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         endif
      
      elseif ( ihadron .eq. 3 ) then ! h

         if (icharge .eq. 0) then ! 0
            IH = 9 ! h+ + h-
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         norm = 0.5d0 ! Transform to charge average ....
         elseif (icharge .eq. 1) then ! +
            IH = 18
            IH = 2
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)

            IH = 12
            CALL AKK(IH,z,Q,DH)

            GL= (GL + DH(0) )/2.d0        
            U = (U  + DH(1) )/2.d0        
            UB= (UB + DH(2) )/2.d0        
            D = (D  + DH(3) )/2.d0        
            DB= (DB + DH(4) )/2.d0        
            S = (S  + DH(5) )/2.d0        
            SB= (SB + DH(6) )/2.d0        
            C = (C  + DH(7) )/2.d0        
            CB= (CB + DH(8) )/2.d0        
            B = (B  + DH(9) )/2.d0        
            BB= (BB + DH(10))/2.d0  ! K^+


         elseif (icharge .eq. 2) then ! -
            IH = 19
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         elseif (icharge .eq. 3) then ! h+ + h-
            IH = 10
            CALL AKK(IH,z,Q,DH)
            GL= DH(0)
            U = DH(1)
            UB= DH(2)
            D = DH(3)
            DB= DH(4)
            S = DH(5)
            SB= DH(6)
            C = DH(7)
            CB= DH(8)
            B = DH(9)
            BB= DH(10)
         endif
               
      elseif ( ihadron .gt. 3 ) then ! p
         IH = 3
         IH = 4
      write(6,*) 'only pi,K,h (ihadron=1,2,3) for AKK08 (ffset=5)'
     > //readme
      stop
      endif
      




      dff(-5) = norm * BB 
      dff(-4) = norm * CB 
      dff(-3) = norm * SB
      dff(-2) = norm * DB
      dff(-1) = norm * UB
      dff( 0) = norm * GL
      dff( 1) = norm * U 
      dff( 2) = norm * D 
      dff( 3) = norm * S 
      dff( 4) = norm * C 
      dff( 5) = norm * B 
 

      else
 
      write(6,*) 'choose ffset=1,2,3,4 for K,KKP,BFGW, FDSS'//readme
      stop      

      endif

      return

      end

