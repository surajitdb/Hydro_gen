        program hydro_gen 
c=============================================================================
c
c  H   H Y        Y DDD     RRRRR      OOOO      GGGGGGG  EEEEEE N         N
c  H   H  Y      Y  D  D    R    R    O    O     G        E      N N       N
c  H   H   Y    Y   D   D   R     R  O      O    G        E      N  N      N
c  H   H    Y  Y    D    D  R     R  O      O    G        E      N   N     N
c  HHHHH     Y      D    D  R    R   O      O    G        EEEE   N    N    N
c  H   H     Y      D    D  R  R     O      O    G   GGG  E      N     N   N
c  H   H     Y      D   D   R   R    O      O    G     G  E      N      N  N
c  H   H     Y      D  D    R    R    O    O     G     G  E      H       N N
c  H   H     Y      DDD     R     R    OOO  ____ GGGGGGG  EEEEEE H         N
c
c============================================================================
c
c Generation of normally ditributed random fields with a given covariance
c function
c                               BY:
c
c               ALBERTO BELLIN^1 AND YORAM RUBIN^2
c
c       1:   Dipartimento di Ingegneria Civile ed Ambientale
c            Universita' di Trento
c            via Mesiano, 77, I-38050 Trento, Italy
c            phone:  +39 461 882620
c            fax:    +39 461 882672
c            e-mail: Alberto.Bellin@unitn.it
c
c       2:   Department of Civil Engineering
c            University of California, Berkeley
c            Berkeley, CA 94720, USA
c            phone:  +1 510 642 2282
c            fax:    +1 510 642 7476
c            e-mail: rubin@ce.berkeley.edu
c
c=============================================================================
c
c Software developed by Alberto Bellin  
c Universita' di Trento,
c Dipartimento di Ingegneria Civile ed Ambientale,
c 38050-I Mesiano di Povo, TRENTO
c
c e-mail:   Alberto.Bellin@unitn.it 
c
c Copyright: Alberto Bellin and Yoram Rubin.
c
c
c
c Summary:  This is the main program of the package aimed at the generation
c           of two dimensional random fields with assigned covariance
c           function.
c
c
c
c Package Version: 2.0 April 1997
c
c
c please cite the following paper in papers or reports that use
c the present software:
c
c Bellin A., Y. Rubin, Hydro_gen: A new random field generator
c for correlated properties, Stochastic Hydrology and Hydraulics,
c 10(4), 1996.
c
c
c Permission is  granted to anyone to use and modify this packages provided
c that:
c i) the authors are acknowledged by citing the abofe referenced paper;
c ii) the use in any kind of research or job will be cited in the relative
c papers or reports;
c iii) the use of the package is under the user responsability
c NO WARRANTY is given concerning bugs and errors.
c iv) The use or distribution must be free of charge.                 
c  v) The package uses the following libraries:
c     a) LINPACK by J. J. Dongarra, J. R. Bunch, C. B. Moler 
c        e G.W. Stewart, for the linear system solution
c     b) BLAS, for linear algebra
c     d) RANLIB by Barry W. Brown and James Lovato,
c        Department of Biomathematics, Box 237
c        the University of Texas, M.D. Anderson Cancer Center
c        1515 Holcombe Boulevard, Huston, TX 77030, for the generation
c        of independent normally distributed random numbers.
c     e) Numerical Recipes by W. H. Press, B. P. Flannery, S. A.
c        Teukolsky, W. T. Vetterling, for the function computing
c        the Bessel Function
c  
c Copyright conditions of the above referenced libraries are 
c extended to hydro_gen  
c
c Bug reports and hints are welcomed to the following e-mail address:
c Alberto.Bellin@unitn.it
c----------------------------------------------------------------------------
c              parameter dimensions description:
c
c  igrid1 => max number of grid points in y direction
c  igrid2 => max number of grid points in x direction
c            BOTH ARE FOR THE COARSE GRID
c  icond  => max number of kriging points
c  iptmx  => max dimension for the vectors storing the kriging
c            coefficients
c            the coefficients are stored in the following order:
c            side #1: kriging area exiting form the field border
c                      x=0
c            side #2: kriging area exiting from the field border
c                     x=lx
c            side #3: kriging area exiting from the field border
c                     y=0
c            corners: kriging area exiting from the two couple
c                     of borders y=0 ,x=0 or y=0, x=lx.
c
c  iptmxcv=> max dimension for the  vectors containing the
c            conditional variances relatives to points referenced
c            in the previous entry.
c            the coefficients are stored in the following order:
c            side #1: kriging area exiting form the field border
c                      x=0
c            side #2: kriging area exiting from the field border
c                     x=lx
c            side #3: kriging area exiting from the field border
c                     y=0
c            corners: kriging area exiting from the two couple
c                     of border y=0 ,x=0 or y=0, x=lx.
c
c  idim   => max dimension of the vector storing the covariance
c            function to be reproduced
c
c
c
c
c
c---------------------------------------------------------------------
c       cond(i,j) -> generated field value at position i,j
c
c
c       u1sid(I) -> vector of the kriging coefficients
c                   which are used in the Monte Carlo generations
c
c       vcxsid(i) -> vector of  the conditional variances
c
c
c--------------------------------------------------------------------------

       
      implicit double precision(a-g,o-z)
       
c      INCLUDE 'dimension.blk' 
      parameter(idim=50,iptmx=5295780,iptmxcv=4536,igrid1=1005,
     %igrid2=1005,ilevmx=4,icond=1670,iwork=icond*(icond+1)/2)
 
 
          

        character*30 filecov,filecoef
        character*12 covartype(5)
        character*1 gset,answ
	  character*12 fileout
	  
c        character*30 file1

        real gennor
c      

      dimension u1sid(iptmx),vcxsid(iptmxcv)  
      dimension u1l1(4,ilevmx),u1l2(4,ilevmx),vcxl1(ilevmx),
     % vcxl2(ilevmx)
      dimension r(igrid1*igrid2)
      dimension cond(igrid1,igrid2)
      DIMENSION u1(icond),cov1(3,4,ilevmx),cov2(10,ilevmx),
     %cov(idim,idim),n(2),cov1a(3,4),cov2a(10)

      real*8 lx,ly
         
      data zero,one/0.0d0,1.0d0/
      data covartype/"Exponential ","Gaussian ","Whittle ",
     %"Mizell (B)","Self-Similar"/

8888  format(///,/,
     %'H   H Y       Y DDD     RRR      OO    GGGGGG  EEEEE N      N',
     %/,'H   H  Y     Y  D   D   R  R    O  O   G       E     N N    N',
     %/,'H   H   Y   Y   D    D  R   R  O    O  G       E     N  N   N',
     %/,'HHHHH    Y Y    D    D  R  R   O    O  G       EEE   N   N  N',
     %/,'H   H     Y     D   D   R R    O    O  G  GGG  E     N    N N',
     %/,'H   H     Y     D  D    R  R    O  O   G    G  E     N     NN',
     %/,'H   H     Y     DDD     R   R    OO __ GGGGGG  EEEEE N      N',
     %////,12x,' BY: ALBERTO BELLIN and YORAM RUBIN '
     %,////)

      write(*,8888)


c
c---------------------------------------------------------------
c     COND(I,J): log-conductivity of position (I*DX,J*DY)]
c
c--------------------------------------------------------------------
c
887       format(4f15.8)
c
c------ Set the seed number for the Gaussian random generation--first
c       inquire whether the use would like to set the seed manually
c       or using the clock:
c
c      write(*,*)'Enter ''y'' to set the seed number manually, or'
c      write(*,*)'enter any other letter to have the clock set the seed:'
c      read(*,'(A1)') gset
c      if (gset .eq. 'y') then
c          write(*,*)'Enter the seed number as an integer'
c          write(*,*)'between 2 and 2147483398:'
c          read(*,*)nseed
c         nseed1 = nseed-1
c      else
c
c  This part of the code is written for Unix sustems.  DOS users
c  must modify the following instructions according
c  to the system specifications:
c
c 1001        call system('date >date.skratch')
c             open(88,file='date.skratch',status='unknown')
c             read(88,336)date
c             read(date(12:13),337)ihour
c             read(date(15:16),337)minute
c             read(date(18:19),337)isec
c337          format(i2)
c336          format(a19)
c             close(88,status='delete')
c             nseed=ihour*100000+minute*1000+isec*10+1
c             nseed1=nseed-1
c         end if

      write(*,*)'Enter the seed number as an integer'
      write(*,*)'between 2 and 2147483398:'
      read(*,*)nseed
         nseed1 = nseed-1

         call setall(nseed,nseed1)         
c
       write(*,*)'SEEDS:',nseed,nseed1
       write(*,*)
c
c-------   INPUT ------------------
c      DX: grid dimension along x
c      DY: grid dimension along y
c      lx & ly: field dimensions along x and y respectively
c      sigy: log-conductivity variance
c      COND10: mean log-conductivity
c      xsp,ysp: dimensions of the area used for the conditioning
c
c-------------------------------------------------------------
        write(*,*)'Enter final grid spacing dx,dy'
        read(*,*)dx,dy
        write(*,'(2x,A,2F6.2)')'dx and dy',dx,dy
        write(*,*)'Enter number of Monte Carlo realizations'
        read(*,*)np
        write(*,*)'# of monte carlo iterations:',np
        write(*,*)'Enter field dimensions'
        read(*,*)lx,ly
        write(*,'(2X,A,2F8.2)')'field dimensions:',lx,ly
        write(*,*)'To read the interpolation "kriging" coefficients'
        write(*,*)'from a file enter [1]; if not,'
        write(*,*)'enter any other integer:'
        read(*,*)imark
        write(*,*)'Enter the covariance type:'
        write(*,*)'itype=0 ==> discrete covariance function (from file)'
        write(*,*)'itype=1 ==> exponential'
        write(*,*)'itype=2 ==> Gaussian'
        write(*,*)'itype=3 ==> Whittle'
        write(*,*)'itype=4 ==> Mizell (B)'
        write(*,*)'itype=5 ==> power law semivariogram'
        write(*,*)'(self similar field)'
        read(*,*)itype

        write(*,101)covartype(itype)
101     format('Covariance type:',a11)

        write(*,*)'covariance of type:',itype 
        
        if(itype.ne.5) then          
           write(*,'(2X,A,F8.2)')'Enter the variance'
           read(*,*)sigy
           write(*,*)'variance:',sigy
        end if
      
         if(itype.eq.5) then
          write(*,*)'do you want to force the mean to be constant in'
          write(*,*)'each realization? [Y/N]'
          read(*,'(A1)')answ
          if(answ.eq.'Y') then
             write(*,*)'enter the mean:'
             read(*,*)fixcond
             cond10=0.0d0
             write(*,'(2X,A,F8.2)')'mean:',fixcond
          end if                                 
        else
          write(*,*)'Enter the mean '
          read(*,*)cond10
          write(*,'(2X,A,F8.2)')'mean:',cond10        
        end if
        
        
        if(imark.ne.1) then
c
        if(itype.eq.1) then
         write(*,*)'Enter integral scales in x and y directions:'
         read(*,*)sclx,scly
         write(*,'(2X,A,2F8.2)')'integral scales:',sclx,scly
         xspsug = DBLE(3.0)*sclx
         yspsug = DBLE(4.0)*scly
         xspasug = sclx
        end if
c
        if(itype.eq.2) then
         write(*,*)'Enter correlation lengths in x and y directions:'
         read(*,*)sclx,scly
         write(*,'(2X,A,2F8.2)')'correlation lengths:',sclx,scly
         xspsug = DBLE(5.0)*sclx
         yspsug = DBLE(6.0)*scly
         xspasug = sclx
        end if
c
      if(itype.eq.3.or.itype.eq.4) then
         write(*,*)'Enter the integral scale (isotropic)'
         read(*,*)sclx
         write(*,'(2X,A,F8.2)')'integral scale (isotropic):',sclx
         scly=sclx
        if (itype .eq. 3) then
           xspsug = DBLE(3.0)*sclx
           yspsug = DBLE(4.0)*scly
           xspasug = sclx
        else
           xspsug = DBLE(5.0)*sclx
           yspsug = DBLE(6.0)*scly
           xspasug = sclx
        end if
      end if

        if(itype.eq.5) then
          write(*,*)'Enter the reference scales in x and y directions:'
          read(*,*)sclx,scly
          write(*,'(2X,A,2F8.2)')'reference  scales:',sclx,scly
          write(*,*)'enter the cofficients a and beta of the power'
          write(*,*)'semivariogram ( g=a r^beta)'
          read(*,*)c,beta
        end if                                                
                                             
        if(itype.ne.5) then
                                       
        write(*,*)'search neighborhood dimensions:'
        write(*,*)'--------           -'
        write(*,*)'|      |           |'
        write(*,*)'|      --------   ysp '
        write(*,*)'|             |    |'
        write(*,*)'|             |    |'
        write(*,*)'---------------    -'
        write(*,*)'|--xsp-|-xspa-|'
        write(*,*)
        write(*,*)'Suggested values:' 
        write(*,'(2X,A,F8.3)')'xsp = ',xspsug
        write(*,'(2X,A,F8.3)')'ysp =',yspsug
        write(*,'(2X,A,F8.3)')'xspa =',xspasug
        write(*,*)
c
        if (xspsug+xspasug .GT. lx-3) then
           write(*,*)'WARNING:'
           write(*,*)'Suggested xsp/xspa is too large.' 
           write(*,*)'Enter xsp and xspa such that their sum'
           write(*,*)'is smaller than:',lx-3
           write(*,*)'or, stop program and begin again with an'
           write(*,*)'x dimension larger than:',xspsug+xspasug+3
        end if
c
        if (yspsug .GT. ly-3) then
           write(*,*)'WARNING:'
           write(*,*)'Suggested ysp is too large.' 
           write(*,*)'Enter an ysp smaller than:',ly-3
           write(*,*)'or, stop program and begin again with a'
           write(*,*)'y dimension larger than:',yspsug+3
        end if
        write(*,*)

        write(*,*)'Enter xsp, ysp'
        read(*,*)xsp,ysp
        write(*,'(2X,A,2F8.3)')'xsp and ysp:',xsp,ysp
        write(*,*)'Enter xspa'
        read(*,*)xspa
        write(*,'(2X,A,F8.3)')'xspa:',xspa
        end if

      end if
      
        write(*,*)'Enter the file name for the file storing the'
        write(*,*)'interpolation coefficients:'
        read(*,10)filecoef  
        open(88,file=filecoef,form='unformatted',status='unknown')

      if(itype.le.4) then
        if(imark.eq.1) then
           read(88)xsp,ysp
           read(88)xspa 
        else
           write(88)xsp,ysp
           write(88)xspa
        end if        
      end if
      
      if(itype.eq.5) then
         if(imark.eq.1) then
           read(88)c,beta
         else
           write(88)c,beta
         end if
      end if
               
        write(*,*)'enter the number of refinement levels:'
        write(*,*)'[0] ==> No refinement'
        write(*,*)'[n] ==> refinement at n levels'
        read(*,*)ilevref

c      write(*,*)'Enter the name  of the output file for the'
c      write(*,*)'replicates [max 30 characters]:'
c      read(*,10)file1
       write(*,*)'Enter the format of the output file:'
      write(*,*)'type 1 for 3 columns x y z'
      write(*,*)'type 0 for the matrix format (only z values'
      write(*,*)'on a regular grid)'
       read(*,*)iformat

10     format(a30)


c
c------------ open output files----------------


c          open(18,file=file1,status='unknown')
          open(20,file='stats.out',status='unknown')
c
c ---------------------------------------------------------
c     SD1: S.D. of the conductivity
c----------------------------------------------------------
        SD1=DSQRT(SIGY)
c
2000    format('particle#: ',i5)
c
c-----------------------------------------------------------
c       N1: # of grid points along y
c       N2: # of grid points along x
c       NN1: # of previous points along y used for the conditioning
c       NN2: # of previous points along x used for the conditioning
c       NN2A: # of following points along x used for the conditioning
c                       


      if(ilevref.eq.0) then
c          iref=2
          idiv=1
      else
c          iref=1
          idiv=2**ilevref
      end if

          ddx=dx*dble(idiv)
          ddy=dy*dble(idiv)
          
      n(1)=int((ly+dy/10.)/dy)+1
      n(2)=int((lx+dx/10.)/dx)+1  
       
      kin=idiv  
      n1in=1
      do kk=1,ilevref
         kin=kin/2
         n1in=n1in+kin
      end do  
      n2in=n1in
                  
      n1fin=n1in+n(1)-1
      n2fin=n2in+n(2)-1
      
      nx=n2fin-n2in+1
      ny=n1fin-n1in+1   
      nxy=nx*ny
            
      if(ilevref.gt.0) then
         n1tot=n1fin+idiv+mod(n(1),2)
         n2tot=n2fin+idiv+mod(n(2),2)
         n1=int(n1tot/idiv)+1
         n2=int(n2tot/idiv)+1
      else
         n1tot=n(1)
         n2tot=n(2)
         n1=n1tot
         n2=n2tot
      end if
      

c      ngen=nx*ny
c      if(ilevref.gt.0) then
c      ngen=ngen+ilevref*(nx+ny+2*ilevref-2)+(n1+n2-2)*2  
c      end if
       ngen=((n1-1)*idiv+1)*((n2-1)*idiv+1)
      
      
      write(*,*)'ngen:',ngen


c        igrid=2*igrid1
c        N2=INT(LX/DX)+1
c        N1=INT(LY/DY)+1
c
c        n2d2=int((2*n2-1)/2)
c        n1d2=int((2*n1-1)/2)
c
c            ilast=2*n1-1
c            jlast=2*n2-1
c        if(iref.eq.2) then
c         isrt=0
c          else
c         isrt=1
c          end if

c        if(iref.eq.1) then
c        n1n2=(2*n1-1)*(2*n2-1)
c        n1n2=n1n2-4*(n1+n2-2)
c        else
c        n1n2=n1*n2
c        end if
c        rn1n2=dble(n1n2)

c        NN2=int(xsp/dx)
c        NN1=int(ysp/dy)
c        NN2A=int(xspa/dx)
c        jfin=nn2+nn2a+1
c 
c
c------- grid spacing for the coarse grid:
c

 

        if(itype.ne.5) then
          nn2=int(xsp/ddx)
          nn2a=int(xspa/ddx)
          nn1=int(ysp/ddy)
        else
c
c---- Power law semivariogram------
c          
          NN2=n2-1
          NN1=n1-1
          NN2A=0
        end if 

        jfin=nn2+nn2a+1
                  
      if(itype.eq.5) then   
          betad2=beta/dble(2.)
          fd=dble(3)-betad2
       
          write(*,*)'the semivariogram is computed with reference to'
          write(*,*)'the following parameters:'
          write(*,*)'C=',c
          write(*,*)'beta=',beta
          write(*,*)'FRACTAL DIMENSION:',fd
          write(*,*)'the employed semivariogram is stored in the file:'
          write(*,*)'semivariog.out'
          open(8,file='semivariog.out')
          spmx=dsqrt(lx*lx+ly*ly)
          nlgs=int(spmx/ddx)+1
          do i=1,nlgs
          xlgs=(i-1)*ddx
          semivar=c*xlgs**beta
          write(8,*)xlgs,semivar
          end do
          close(8,status='keep')
       end if
                         
        ierror=0
        if(n1tot.gt.igrid1) then
           write(*,*)'DIMENSIONING ERROR: the parameter igrid1'
           write(*,*)'in the file dimension.blk is not properly'
           write(*,*)'assigned. Suggested value:',n1tot+1
           ierror=ierror+1
        end if
c
        if(n2tot.gt.igrid2) then
           write(*,*)'DIMENSIONING ERROR: the parameter igrid2'
           write(*,*)'in the file dimesnion.blk is not properly'
           write(*,*)'assigned. Suggested value:',n2tot+1
           ierror=ierror+1
        end if
c 
        ijcount=(nn2+nn2a+1)*nn1+nn2
        if(ijcount.gt.icond) then
           write(*,*)'DIMENSIONING ERROR: the parameter icond'
           write(*,*)'in the file dimension.blk is not properly'
           write(*,*)'assigned. Suggested value:'
     %     ,ijcount
           ierror=ierror+1
        end if
c
        if(ierror.gt.0) stop

c         write(18,101)covartype(itype)
c         write(18,102)np
c102     format('This file stores',i4,'independent replicate(s)')
c        
c        if(iformat.eq.0) then
c          write(18,*)'the data are stored in the matrix format with:'
c          write(18,67)nx,ny        
c        else
c          write(18,*)'the data are stored in column format (x y z)'
c          write(*,*)'total number of data:',nx*ny
c
c        end if


c67     format(i3,1x,'lines and',1x,i3,1x,'columns')
c          write(18,*)'grid size:',dx,dy
c      write(18,*)'---------------------------------------------------'
      write(20,*)'replicate N.           mean        variance'
      write(20,*)'-------------------------------------------'
c
c------ compute the covariance matrix --
c
c      itype=0   discrete covariance function (it is read from a file)
c      itype=1   exponential covariance function (defined by a function)

c
       if(itype.eq.0) then
       write(*,*)'the covariance function is read from the file:',
     %  filecov
       read(*,442)filecov
       end if
442    format(a30)


c-----------------------------------------------
c
      if(imark.ne.1) then
c
      
      ierror=max(nn1+1,jfin)
      if(ierror.gt.idim) then
         write(*,*)'DIMENSIONING ERROR: the parameter idim'
         write(*,*)'in the file hydro_gen.f is not properly'
         write(*,*)'assigned. Suggested value:',ierror
         stop
      end if

      nvectpos=0
      kkcv=0
      ierror = 0

c
      call covar(filecov,idim,ilevref,nn1,jfin,itype,ddx,ddy,
     %sclx,scly,c,betad2,cov,cov1,cov2)

      end if
      
      if(itype.eq.5) then
c
c----  self-similar RF---------
c
      if(imark.ne.1) then
c
c---- compute the kriging coefficients
c                               
      write(*,'(2x,A)')'line   nvectpos   kkcv'
      jfin=n2 
      i=1
      do j=2,n2      
              ijcountm1=j-1
              call krig(0,j,idim,j,itype,cov,u1,vcx) 
              write(88)(u1(kcount),kcount=1,ijcountm1)
              write(88)dsqrt(vcx)
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
                u1sid(nvectpos)=u1(kcount)
             end do             
 
             kkcv=kkcv+1
             vcxsid(kkcv)=dsqrt(vcx)       
      end do 
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv  

      do i=2,n1            
         do j=1,n2  
             ijcount=(i-1)*n2+j  
             ijcountm1=ijcount-1
c            write(*,*)i,j
             call krig(0,ijcount,idim,n2,itype,cov,u1,vcx) 
             write(88)(u1(kcount),kcount=1,ijcountm1)
             write(88)dsqrt(vcx)
             do kcount=1,ijcountm1
                nvectpos=nvectpos+1
                u1sid(nvectpos)=u1(kcount)
             end do    
             kkcv=kkcv+1
             vcxsid(kkcv)=dsqrt(vcx)

         end do 
c         write(*,*)'line:',i,'done!'
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv
      end do   

      else
c
c----- read the kriging coefficients
c
       i=1 
      write(*,'(2x,A)')'line   nvectpos   kkcv'       
       do j=2,n2
          kkcv=kkcv+1
          read(88)(u1sid(k+nvectpos),k=1,j-1)
          read(88)vcxsid(kkcv)
          nvectpos=nvectpos+j-1
       end do

      write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv     
       do i=2,n1
          do j=1,n2
             ijcount=(i-1)*n2+j-1
             kkcv=kkcv+1
             read(88)(u1sid(k+nvectpos),k=1,ijcount)
             read(88)vcxsid(kkcv)
             nvectpos=nvectpos+ijcount
          end do
c        write(*,*)'line:',i,'done!'                    
         write(*,'(2x,i3,i8,i8)')i,nvectpos,kkcv          
       end do

      end if           
      else
c
c
c----- compute the kriging coefficients for kriging area that exits
c      from the side x=0.0
c

       write(*,*)'side #1'


      do icol=1,nn2
      ijcount=(nn2a+icol)*nn1+icol

c---------- compute the kriging coefficients -------
c
      kkcv=kkcv+1
      if(imark.ne.1) then
      jfin=nn2a+icol
      call krig(0,ijcount,idim,jfin,itype,cov,u1,vcxsid(kkcv))      
      write(88)(u1(j),j=1,ijcount-1)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv) 

      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do
c
c----- kriging coefficients for kriging area that exits
c      from the side x=lx
c
c------- iptkr2  pointer for the kriging coefficients of side #2
c                the vector position from iptkr2 to iptkr3-1 contain
c                the kriging coefficients of side #2.
c                The positions between 1 and iptkr2 contain
c                the coefficients of side #1
c        iptcv2  pointer for the conditional variances of side #2
c                the vector position from iptcv2 to iptcv3-1 contain
c                the conditional variances of side #2.
c                The positions between 1 and iptcv2 contain
c                the conditional variances of side #1
c      write(88)nvectpos,kkcv
       iptkr2=nvectpos
       iptcv2=kkcv

       write(*,*)'# of vector positions for the side #1:',iptkr2,kkcv
c
      write(*,*)'side #2'

      do icol=1,nn2a
      jfin=nn2+nn2a+1-icol
      ijcount=jfin*nn1+nn2+1
c------- compute the kriging coefficients--------------
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,idim,jfin,itype,cov,u1,vcxsid(kkcv))
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)      
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do               
 
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do
c
c------- iptkr3  pointer for the kriging coefficients of side #3
c                the vector position from iptkr3 to iptkrc-1 contain
c                the kriging coefficients of side #3.
c        iptcv3  pointer for the conditional variances of side #3
c                the vector position from iptcv3 to iptcvc-1 contain
c                the conditional variances of side #3.
c      write(88)nvectpos,kkcv     
      iptkr3=nvectpos
      iptcv3=kkcv

       write(*,*)'# of vector positions for the side #2:',nvectpos,kkcv

c
c------ kriging coefficients along the line y=0.0
c
        write(*,*)'side #3'

      jjfin=nn2+nn2a+1
      jfin1=nn2+1
                do iline=1,nn1
      if(iline.eq.1)  then
             jfin=jfin1
          else
             jfin=jjfin
          end if

        ijcount=jfin*(iline-1)+jfin1
c----- compute the kriging coefficients ------
       kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,idim,jfin,itype,cov,u1,vcxsid(kkcv))
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do  
   
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if
      end do

c
c------- iptkrc  pointer for the kriging coefficients of corners
c                the vector position from iptkrc to the iptkrs contain
c                the kriging coefficients of corners.
c        iptcvc  pointer for the conditional variances of side #3
c                the vector position from iptcvc to the end contain
c                the conditional variances of corners.
       
c       write(88)nvectpos,kkcv   
       iptkrc=nvectpos
       iptcvc=kkcv

        write(*,*)'# of vector positions for the side #3:',nvectpos,kkcv
c
c----- compute the remaining kriging coefficients----
c
        write(*,*)'corners:'


        DO 111 I=1,NN1
        jstart=1
        jend=nn2

112     DO 110 J=jstart,jend
c
c
c
c----- select the points used for the kriging--------
       IINIZ=I-NN1
       JINIZ=J-NN2
       idec=1
        jfin=j+nn2a
         if(jiniz.ge.1.and.jfin.le.n2)  goto 110
       iiniz=max(iiniz,1)
       jiniz=max(jiniz,1)
       jfin=min(jfin,n2)
        jfin1=j
        if(i.eq.1) jfin=j
c
c----- select the points used for the kriging--------
c
        jjfin=jfin-jiniz+1
        jjfin1=jfin1-jiniz+1
        inum=I-IINIZ
        ijcount=jjfin*inum+jjfin1

                if(ijcount.eq.1)  go to 110

c
c------------- compute the kriging coefficients -------
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,idim,jjfin,itype,cov,u1,vcxsid(kkcv))
      write(88)(u1(kcount),kcount=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do kcount=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(kcount)
      end do
      else
      read(88)(u1sid(nvectpos+kcount),kcount=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      end if

110         continue
      if(jstart.eq.1) then
      jstart=n2-nn2a+1
      jend=n2
        go to 112
          end if
111   continue
c
c------- iptkrs  pointer for the kriging coefficients of the larger
c                kriging area
c                the vector position from iptkrs to the end of
c                vectors contain the kriging coefficients of corners.
c      write(88)nvectpos
      iptkrs=nvectpos
       write(*,*)'# of vector positions for the corners:',
     & nvectpos,kkcv
c
c----------- KRIGING COEFFICIENTS FOR THE LARGER CONDITIONING AREA
c
      jfin=nn2+nn2a+1
      ijcount=nn1*jfin+nn2+1
      ijmax=ijcount

c------- compute the kriging coefficients-------
c      write(*,*)'compute the kriging coefficients for the larger area'
      kkcv=kkcv+1
      if(imark.ne.1) then
      call krig(0,ijcount,idim,jfin,itype,cov,u1,vcxsid(kkcv))
      write(88)(u1(j),j=1,ijcount-1)
      vcxsid(kkcv)=dsqrt(one-vcxsid(kkcv))
      write(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)
      do k=1,ijcount-1
         nvectpos=nvectpos+1
         u1sid(nvectpos)=u1(k)
      end do   
      iptcvst=kkcv
      else
      read(88)(u1sid(nvectpos+j),j=1,ijcount-1)
      nvectpos=nvectpos+ijcount-1
      read(88)vcxsid(kkcv)
      vcxsid(kkcv)=sd1*vcxsid(kkcv)   
      iptcvst=kkcv
      end if 

      end if       

c88       format(6i8)
       write(*,*)'total # of vector positions:',nvectpos,kkcv
c
       ierror=0
       if(iptmx.lt.nvectpos) then
        write(*,*)'ERROR IN DIMENSIONING THE VECTOR STORING'
        write(*,*)'THE INTERPOLATION COEFFICIENTS:'
        write(*,*)
        write(*,*)'ACTION: set the parameter iptmx to:',nvectpos
        write(*,*)'in the file hydro_gen.f and recompile the'
        write(*,*)'package. As an alternative, reduce the dimensions'
        write(*,*)'of the larger search neighborhood or increase the'
        write(*,*)'coarse grid spacing.'
        ierror=ierror+1
       end if
c
       if(iptmxcv.lt.kkcv) then
        write(*,*)'ERROR IN DIMENSIONING THE VECTOR STORING'
        write(*,*)'THE CONDITIONAL VARIANCES:'
        write(*,*)
        write(*,*)'ACTION: set the parameter iptmxcv to:',kkcv
        write(*,*)'in the file hydrogen.dim and compile again the'
        write(*,*)'package. As an alternative reduce the dimensions'
        write(*,*)'of the larger search neighborhood or increase the'
        write(*,*)'coarse grid spacings'
        ierror=ierror+1
       end if
c
       if (ierror .gt. 0) stop



c
c---------- compute the kriging coefficients for
c           the local interpolation
c
c-------- Interpolation level # 1  
c
      ijcount=5
      iptl1=nvectpos
      iptcv1=kkcv    
      

       do kk=1,ilevref  
      
      if(imark.ne.1) then                                  
      
c------- compute the kriging coefficients for the level #1        
          jfin=2   
          do ii=1,4
             u1l1(ii,kk)=dble(0.0)
             u1l2(ii,kk)=dble(0.0)
          end do         
  

      do ii=1,3
         do jj=1,4
            cov1a(ii,jj)=cov1(ii,jj,kk) 
         end do
      end do
      
c
      call krig(1,ijcount,3,jfin,itype,cov1a,u1,vcxl1(kk))
      do ii=1,4
         u1l1(ii,kk)=u1(ii)
      end do   
      write(88)(u1(ii),ii=1,4)
      
      if(itype.ne.5) then
         vcxl1(kk)=dsqrt(one-vcxl1(kk))
         write(88)vcxl1(kk)
         vcxl1(kk)=sd1*vcxl1(kk)   
      else
      vcxl1(kk)=dsqrt(vcxl1(kk))
      write(88)vcxl1(kk)   
      end if
      
         

c      write(*,*)'level # 1:'
c      write(*,*)(u1l1(k),k=1,4)
c      write(*,*)dsqrt(vcxl1)

c------- compute the kriging coefficients for the level #2 

      do ii=1,10
         cov2a(ii)=cov2(ii,kk)
      end do
c
          call coefl2(itype,cov2a,u1,vcxl2(kk))       
c
      do ii=1,4
         u1l2(ii,kk)=u1(ii)
      end do   
      write(88)(u1(ii),ii=1,4)
      
      if(itype.ne.5) then
         vcxl2(kk)=dsqrt(one-vcxl2(kk))
         write(88)vcxl2(kk)
         vcxl2(kk)=sd1*vcxl2(kk)   
      else
      vcxl2(kk)=dsqrt(vcxl2(kk))
      write(88)vcxl2(kk)         
      end if
      
      else
c      
c---- read the kriging coefficients for the refinement levels
c
c
c------ level # 1
c      
         read(88)(u1l1(ii,kk),ii=1,4)
         read(88)vcxl1(kk)
         if(itype.ne.5) vcxl1(kk)=sd1*vcxl1(kk)
c
c------- level # 2
c
         read(88)(u1l2(ii,kk),ii=1,4)
         read(88)vcxl2(kk)
         if(itype.ne.5) vcxl2(kk)=sd1*vcxl2(kk)
      end if
      
                        
      end do
      close(88)

c
c----------------START MONTE CARLO SIMULATION--------------
c
c  NP  Number of Monte Carlo Simulations.
c  i    row indicator
c  j    column indicator
c
c----- pointers -------
c
c      zone       krig. coef.    cond. variance
c ----------------------------------------------
c |   side #1 |     ip1s1     |    ipnod1s1    |
c |   side #2 |     ip1s2     |    ipnod1s2    |
c |   side #3 |     ips3      |    iptcv3      |
c |   corners |     ip1       |    ipnod1      |
c |larger area|     iptkrs    |     ---        |
c-----------------------------------------------
c
c
      WRITE(*,*)
      WRITE(*,*)' kriging coefficients computed; starting Monte Carlo..'

        
        iprint=0
        DO 200 JJ=1,NP
c
c----- generate the set of normally distributed random numbers---
c

        do k=1,ngen
           r(k)=dble(gennor(0.e0,1.e0))
        end do


c
c
c--------------Initialize  RF to zero.--------
c
       do i=1,n1tot
       do j=1,n2tot
       cond(i,j)=zero
       end do
       end do
 
c
c----- coarse grid ----
c
c    
        irv=1
        ivectpos=0
        kkcv=0  
        
        if(itype.eq.5) then
     
c----- Fractal Field
        i=1
        do j=1,n2
           if(j.eq.1) then
              cond(1,1)=r(irv)            
           else
              kkcv=kkcv+1
              CALL comb(igrid1,idiv,i,j,1,1,j,irv,vcxsid(kkcv),
     1 u1sid(ivectpos+1),r,cond)    
           ivectpos=ivectpos+j-1
           end if
        end do

 
        do i=2,n1
           do j=1,n2
             kkcv=kkcv+1
             CALL comb(igrid1,idiv,i,j,1,1,n2,irv,vcxsid(kkcv),
     1          u1sid(ivectpos+1),r,cond)     
             ivectpos=ivectpos+(i-1)*n2+j-1
           end do
        end do   


       else         
c      
c----- regular random field
c
c
c---- set the pointers for the kriging ceofficients and the conditional
c     variances at the filed corners

         ip1=iptkrc
         ipnod1=iptcvc
         ips3=iptkr3

c
c------cycle through each row and column of the field-------
c
c 

      do 120 i=1,n1


      ip1s1=0
      ip1s2=iptkr2
      ipnod1s1=0
      ipnod1s2=iptcv2


c
           do 121 j=1,n2
c
c   if we are at the first node, then we have no kriging to do,
c   create a random number.
c
             if(i.eq.1.and.j.eq.1) then
               cond(1,1)=sd1*r(irv)
               go to 121
             end if

             IINIZ=I-NN1
             JINIZ=J-NN2
             jfin=j+nn2a

             IF(iiniz.lt.1) then
c
c  that is, are we along side #3?
c
               IF(jiniz.ge.1.and.jfin.le.n2) then
c
c  are we along side #3 but not on the corner?
c
c----- compute the velocities along the side # 3
c    jiniz
c      |
c      v
c      --------------------* (i,j)       Note:  ijcount is the number of
c      |   |   |   |   |   |                    kriging coeffs for
c      |---|---|---|---|---|---|---| _          this position on
c      |   |   |   |   |   |   |   | ^          side #3
c -----|---|---|---|---|---|---|---| nn1-----------
c      |   |   |   |   |   |   |   | v     edge of random field.
c   -> ----------------------------- -
c iiniz                            ^
c      |<-----nn2--------->| nn2a  | jfin
c

                 iiniz=1
c  loop through all kriging points to extract kriging coefficients
c

      ij=(jfin-jiniz+1)*(i-1)+nn2
      ij3=ij


      call comb(igrid1,idiv,i,j,iiniz,jiniz,jfin,irv,vcxsid(i+iptcv3),
     %u1sid(ips3+1),r,cond)
c                write(*,*)'side #3',i,j,ips3
c
               ELSE
c
c  we are on a corner.
c
                 iiniz=1
                 jiniz=max(jiniz,1)
                 jfin=min(jfin,n2)
                 ij=(i-iiniz)*(jfin-jiniz+1)+j-jiniz

                  ipnod1=ipnod1+1
c             sc1=dsqrt(sigy-sigy*vcxkr(ipnod1))
c  loop over all kriging coefficients for the corner to extract
      call comb(igrid1,idiv,i,j,iiniz,jiniz,jfin,irv,vcxsid(ipnod1),
     %u1sid(ip1+1),r,cond)
        ip1=ip1+ij
               END IF
             ELSE
c
c  we are one side #1
c
               IF(jiniz.lt.1.and.jfin.le.n2) then
                 jiniz=1
                 ij=(jfin-jiniz+1)*nn1+j-1


                 ipnod1s1=ipnod1s1+1
c             sc1=dsqrt(sigy-sigy*vcxsid1(ipnod1s1))
      call comb(igrid1,idiv,i,j,iiniz,jiniz,jfin,irv,vcxsid(ipnod1s1),
     %u1sid(ip1s1+1),r,cond)
        ip1s1=ip1s1+ij
                 go to 121
               END IF
C
c  we are on side #2
c
               IF(jfin.gt.n2) then
                 jfin=n2
                 ij=(jfin-jiniz+1)*nn1+j-jiniz

                 ipnod1s2=ipnod1s2+1
c             sc1=dsqrt(sigy-sigy*vcxsid2(ipnod1s2))
      call comb(igrid1,idiv,i,j,iiniz,jiniz,jfin,irv,vcxsid(ipnod1s2),
     %u1sid(ip1s2+1),r,cond)
          ip1s2=ip1s2+ij
c      nvectpos2=nvectpos2+(jfin-jiniz+1)*nn1+j-jiniz
c                write(*,*)'side #2:',i,j,nvectpos2,ipnod1s2
                 go to  121
               END IF
c
c  we are in the greater kriging area
c
               IF(jiniz.ge.1.and.jfin.le.n2) then
               ij=ijmax-1


      call comb(igrid1,idiv,i,j,iiniz,jiniz,jfin,irv,vcxsid(iptcvst),
     %u1sid(iptkrs+1),r,cond)
c                write(*,*)'steady:',i,j,ijcount
               END IF
             END IF
c
c  Use extracted kriging coefficients to determine value of conductivity
c
121    continue

      ips3=ips3+ij3
120   continue

      end if
 
      if(ilevref.gt.0) then  

c
c---- multistage grid refinement----
c       
c
c---- modifications: August 22, 1997
c
      n1l1=n1tot-idiv+1
      n2l1=n2tot-idiv+1  
      n1finl2=n1tot
      n2finl2=n2tot
      
c
c------------------------------
c
      jump1=idiv  
      jump2=idiv/2
      inl1=1 
      inl2=idiv/2 +1
      infl2=inl2-jump2
      
      
      do kk=1,ilevref  
 

c
c------- first level refinement -----


       do 150 i1=inl1,n1l1,jump1
       i0=i1+jump1/2
 
       do 150 j1=inl1,n2l1,jump1   
       j0=j1+jump1/2  

       
       ijcount=0
       
       do i2=i1,i1+jump1,jump1

       do j2=j1,j1+jump1,jump1
       ijcount=ijcount+1 

       cond(i0,j0)=cond(i0,j0)+u1l1(ijcount,kk)*cond(i2,j2)
       end do
       end do 

c
c   Add the fluctuation
c
       irv=irv+1
       cond(i0,j0)=cond(i0,j0)+vcxl1(kk)*r(irv)
150    continue
       

c       
c---- second level refinement ----
c       
       
       do 151 i1=inl2,n1finl2,jump1
       do 151 j1=inl2,n2finl2,jump1  
       

       i1a=i1
       j1a=j1-jump2
       if(j1a.ne.infl2) then   
 

              
       cond(i1a,j1a)=cond(i1a,j1a)+
     %cond(i1a,j1a-jump2)*u1l2(1,kk)+
     %cond(i1a,j1a+jump2)*u1l2(3,kk)+
     %cond(i1a-jump2,j1a)*u1l2(2,kk)+
     %cond(i1a+jump2,j1a)*u1l2(4,kk)
c
c  Add the fluctuation
c
       irv=irv+1
       cond(i1a,j1a)=cond(i1a,j1a)+vcxl2(kk)*r(irv)
       end if

       j1a=j1
       i1a=i1-jump2

       if(i1a.ne.infl2) then  


       cond(i1a,j1a)=cond(i1a,j1a)+
     %cond(i1a,j1a-jump2)*u1l2(1,kk)+
     %cond(i1a,j1a+jump2)*u1l2(3,kk)+
     %cond(i1a-jump2,j1a)*u1l2(2,kk)+
     %cond(i1a+jump2,j1a)*u1l2(4,kk)
c
c  Add the fluctuation
c
       irv=irv+1
       cond(i1a,j1a)=cond(i1a,j1a)+vcxl2(kk)*r(irv)
       end if
151    continue 
 
        jump1=jump1/2
       inl1=inl1+jump1 
       jump2=jump2/2
       infl2=inl2
       inl2=inl2+jump2    
       n1finl2=n1finl2-jump2
       n2finl2=n2finl2-jump2 
       
       end do   
    
       if(irv.gt.ngen) then  
           write(*,*)'ngen=',ngen
           write(*,*)'code error: increase to',irv,'the parameter ngen'
           stop
       end if

      end if
c
c----- add the mean ----
c

      sumy=dble(0.0)
      varlogy=dble(0.0)
      condm=dble(0.0)
      condv=dble(0.0)                     
c      
c*******************************************************
c
c        STATISTICS
c
c*******************************************************
c
c---------- spatial statistics------------
      
      do i=n1in,n1fin
         do j=n2in,n2fin
            condm=condm+cond(i,j)
            condv=condv+cond(i,j)*cond(i,j)
            cond(i,j)=cond(i,j)+cond10
         end do
      end do    
      
      condm=condm/dble(nxy)
      condv=condv/dble(nxy-1) - condm*condm 
      condm=condm+cond10  
      
      if(itype.eq.5.and.answ.eq.'Y') then
      
      do i=n1in,n1fin
         do j=n2in,n2fin
            cond(i,j)=cond(i,j)-condm+fixcond
         end do
      end do                                 
      condm=fixcond
      end if
               
      
             write(20,9999)jj,condm,condv
9999         format(3x,i5,7x,2f15.8)      
c
c------ print the field -----
c
c		
	fileout='real0000.dat'

333   format(I1)
334   format(I2)
335   format(I3)
336   format(I4)
	
		
	if(jj.lt.10) then
	   write(fileout(8:8),333)jj
	end if

	if(jj.ge.10.and.jj.lt.100) then
	   write(fileout(7:8),334)jj
	end if

	if(jj.ge.100.and.jj.lt.1000) then
	   write(fileout(6:8),335)jj
	end if  

	if(jj.ge.1000.and.jj.lt.10000) then
	   write(fileout(5:8),336)jj
	end if  

      open(18,file=fileout)


c          write(18,*)'replicate #',jj

            if(iformat.eq.1) then   

             do i=n1in,n1fin
               do j=n2in,n2fin
           xprint=(j-1)*dx
           yprint=(i-1)*dy
           write(18,*)xprint,yprint,cond(i,j)
               end do
             end do
            else
             do i=n1in,n1fin
           write(18,866)(cond(i,j),j=n2in,n2fin)
             end do
             end if


866        format(1f20.10)
          close(18,status='keep')
200       continue


       close(20,status='keep')

        WRITE(*,*)' PROGRAM COMPLETE.'
        STOP
        END



