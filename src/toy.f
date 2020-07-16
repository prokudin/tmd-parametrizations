      program toy
      implicit double precision ( a-z )
      dimension dff(-5:5)
      integer ffset, fforder, iparton, ipi, ifini, old1, 
     >        iset, ihadron, icharge, old2, old3, ih, icp
      integer n
      data ipi / 1 /  


      ffset    = 4 ! DSS
      fforder  = 0 ! LO
      ihadron  = 1 ! PION
      icharge  = 3 ! +&- = charge sum
      icp      = 1

      z  = 0.2d0
      Q2 = 1.d4

      call dlib(z,Q2,dff,ffset,fforder,ihadron,icharge,icp,ipi) 
      write(6,*) dff

      su3sing = 0.d0
      do 10 n = -3, 3
      if (n .eq. 0) goto 10
      su3sing = su3sing + dff(n)
 10   continue
      write(6,*) su3sing

      end
