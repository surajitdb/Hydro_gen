c
      subroutine covar(filecov,idim,ilrf,nn1,jfin,indec,dx,dy,
     %sclx,scly,c0,omega2,cov,cov1,cov2)  


c
c----- This subroutine compute the covariance matrix
c      or the semivariogram matrix when power law semivariograms
c      are employed. 
c      It uses the function covy(x,y) which gives the
c      covariance function or the semivariogram for the lag (x,y)
c
c      For kriging systems defined by the semivariogram the kriging
c      coefficients are computed with reference to a pseudo covariance
c      function defined as: C(x_x,x_j)=A-gamma(x_i,x_j), where A
c      is a positive constant greater than the greatest semivariogram used
c      in the kriging system (Journell and Huijbregts, Mining Geostatistics
c      pg 306).
c           This subroutine computes the variogram and the main program
c      computes the pseudo covariance employed in the kriging system
c

      
      implicit double precision(a-g,o-z)  
      
      parameter(ilevmx=4)
    
      dimension cov(idim,*),cov1(3,4,ilevmx),cov2(10,ilevmx)
      character*30 filecov
      external covy
      data zero,two/0.0d0,2.0d0/
      
      

      
      if(ilrf.gt.ilevmx) then
          write(*,*)'increases the parameter ilevmx in the files'
          write(*,*)'covariance.f and hydro_gen.fto the following'
          write(*,*)'value:',ilrf
          stop
      end if

      if(indec.eq.0) then
c
c---- the covariance function is read from a file ---
c
      open(88,file=filecov)
     
      do i=1,nn1+1
      read(88,*)(cov(i,j),j=1,jfin)
      end do  
      
      do kk=1,ilrf
c----- level # 1 ------
       do i=1,2
      read(88,*)(cov1(i,j,kk),j=1,2)
       end do
       read(88,*)(cov1(3,j,kk),j=1,4)
c
c----- level # 2 -------
c
         do i=1,10
         read(88,*)cov2(i,kk)
           end do 
      end do

       close(88,status='keep')

      
       else

      do i=1,nn1+1
        y=dble(i-1)*dy/scly
      do j=1,jfin
        x=dble(j-1)*dx/sclx
        cov(i,j)=covy(x,y,c0,omega2,sclx,indec) 
      end do
      end do

        dx1=dx/sclx
        dy1=dy/scly  


        
        do kk=1,ilrf
        
c
c------ level # 1------
c
        do i=1,2
        y=dble(i-1)*dy1
        do j=1,2
        x=dble(j-1)*dx1
        cov1(i,j,kk)=covy(x,y,c0,omega2,sclx,indec)
        end do
        end do
        cov1(3,1,kk)=covy(-dx1/two,-dy1/two,c0,omega2,sclx,indec)
        cov1(3,2,kk)=covy(dx1/two,-dy1/two,c0,omega2,sclx,indec)
        cov1(3,3,kk)=covy(-dx1/two,dy1/two,c0,omega2,sclx,indec)
        cov1(3,4,kk)=covy(dx1/two,dy1/two,c0,omega2,sclx,indec)
c
c---- level # 2
c

      cov2(1,kk)=covy(zero,zero,c0,omega2,sclx,indec)
      cov2(2,kk)=covy(dx1/two,-dy1/two,c0,omega2,sclx,indec)
      cov2(3,kk)=covy(dx1,zero,c0,omega2,sclx,indec)
      cov2(4,kk)=covy(dx1/two,dy1/two,c0,omega2,sclx,indec)
      cov2(5,kk)=covy(zero,dy1,c0,omega2,sclx,indec)
      cov2(6,kk)=covy(-dx1/two,dy1/two,c0,omega2,sclx,indec)
      cov2(7,kk)=covy(-dx1/two,zero,c0,omega2,sclx,indec)
      cov2(8,kk)=covy(zero,-dy1/two,c0,omega2,sclx,indec)
      cov2(9,kk)=covy(dx1/two,zero,c0,omega2,sclx,indec)
      cov2(10,kk)=covy(zero,dy1/two,c0,omega2,sclx,indec)
                
       dx1=dx1/two
       dy1=dy1/two
       end do
       

       end if

c      write(*,*)'------- covariance matrix -----'
c      do i=1,nn1+1
c      write(*,888)(cov(i,j),j=1,jfin)
c      end do
c         write(*,*)'level # 1'
c         do i=1,2
c       write(*,*)(cov1(i,j),j=1,2)
c         end do
c       write(*,*)(cov1(3,j),j=1,4)
c
c        write(*,*)'level # 2'
c         do i=1,10
c          write(*,*)cov2(i)
c         end do 
c888   format(10f8.5)
      return
      end
 
