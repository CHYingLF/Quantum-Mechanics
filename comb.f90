        Program Comb
          !combine 2 energy profile and normalize diffusion coordinates
          !C. Ying, 06/15/2021
          real, dimension(6,3) :: p12
          real, dimension(5,3) :: p23
          real, dimension(10,2) :: p13
          integer :: i, j

          OPEN(1,FILE='../p1_2_p2/neb1/peptide_neb.dat')
            do i=1,6
              READ(1,*) p12(i,:)
            enddo
          CLOSE(1)

          OPEN(1,FILE='../p2_2_p3/neb1/peptide_neb.dat')
            do i=1,5
              READ(1,*) p23(i,:)
            enddo
          CLOSE(1)
  
!        write(*,*) p23(5,1)
          p23(:,1)=p23(:,1)+1.0
          p23(:,2)=p23(:,2)+p12(6,2)

          p13(1:6,1:2)=p12(:,1:2)
          p13(7:10,1:2)=p23(2:5,1:2)
          p13(:,1)=p13(:,1)/2.0

          OPEN(1,FILE='ra_clo4_p13.dat')
            do i=1,10
               write(1,*) p13(i,:)
            enddo
          CLOSE(1)

          

        end program
