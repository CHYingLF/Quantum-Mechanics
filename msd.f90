      program msdcal
      implicit none
!--------------------------------------------------------------------------
!      This code calculate msd of Lithium
!      Before run, update Nt, Natoms, NLi,Nstart,Nend,
!      Vol(volume),Temp,filename
!--------------------------------------------------------------------------
      integer,parameter :: Nt=580,Natoms=230,NLi=230,Nstart=300,Nend=Nt
      real::Nstep(Nt),t(Nt),r(Nt,Natoms,3),rLi(Nt,3),msd(Nt/2),D(Nt/2),&
           &msd2(Nt)
      real::kB=1.38064852e-23,Nd=1,e=1.60217662e-19,Z=1,Temp=300.0,sigma
      real::Vol=2.003762E+04,dt
      integer:: i,j,k,tao,ntao
     
      open(90,file='peptide.pos')  
        do i=1,Nt
          read(90,*)Nstep(i),t(i)
          do j=1,Natoms
            read(90,*)r(i,j,1:3)
          enddo
        enddo
      close(90)
      dt=t(2)-t(1)

      

!-----Lithium trajectory, convert from a.u. to cm
      do i=1,Nt
        rLi(i,1:3)=r(i,NLi,1:3)*5.29177e-11*1e2
      enddo      
!-----msd in cm^2
      do tao=1,(Nend-Nstart)/2 
      msd(tao)=0.0
      ntao=0
      do i=Nstart,Nend
        j=i+tao
        if(j.le.Nend) then
          ntao=ntao+1
          msd(tao)=msd(tao)+sum(abs(rLi(j,1:3)-rLi(i,1:3)))**2
        endif
      enddo !end i
      msd(tao)=msd(tao)/real(ntao)
      enddo !end tao

!-----msd2 in angstrom^2
      do i=1,Nt
        msd2(i)=sum((rLi(i,1:3)-rLi(1,1:3))**2)*(1e8)**2
      enddo

!-----duffusion coefficient, cm^2/s
      do i=1,(Nend-Nstart)/2
        D(i)=msd(i)/(i*dt*1e-12)/3.0/2.0
      enddo
!-----ionic conductivity
      Vol=Vol*(5.29177e-11*1e2)**3  !in cm^3
      Nd=Nd/Vol
      sigma=Nd*(Z*e)**2/(kB*1e4)/Temp*D(1)
!-----output, msd in angstorm^2, D in cm^2/s 
      write(*,*)'diffusion coefficient=',D(1)
      write(*,*)'ionic conductivity',sigma
      open(unit=90,file='msd.dat',form='formatted',status='unknown') 
        do i=1,(Nend-Nstart)/2
          write(90,*)i*dt,msd(i)*(1e10)**2*(1e-2)**2,D(i),log(D(i))
        enddo
      close(90)

      open(unit=90,file='msd2.dat',form='formatted',status='unknown')
        do i=1,Nt
          write(90,*)i*dt,msd2(i)
        enddo
      close(90)

      end program
