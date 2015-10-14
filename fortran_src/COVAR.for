c file covar.f
c
c-------------------------------------------------------------------
c Software developed by Alberto Bellin  
c Universita' di Trento,
c Dipartimento di Ingegneria Civile ed Ambientale,
c 38050-I Mesiano di Povo, TRENTO
c
c e-mail:   Alberto.Bellin@ing.unitn.it 
c
c Copyright: Alberto Bellin and Yoram Rubin.
c
c
c
c Summary:  This function provides different types of covariance
c           functions, namely: the exponential, the Gaussian (bell shape)
c           the Whittle (see Mizell et al. WRR 18(4), 1053-1067, 1982),
c           the Mizell B (see Mizell et al. WRR 18(4), 1053-1067, 1982). 
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
c Alberto.Bellin@ing.unitn.it
c----------------------------------------------------------------------------
c

      real*8 function covy(x,y,c0,omega2,scl,indec)
                                    
c
c----- exponential covariance function ----
c

      implicit double precision(a-g,o-z)

      
      external bessk0,bessk1 
      data pi/3.141592654d0/

c
      r2=x*x+y*y
c#1
      if(r2.gt.0.0) then

         r=dsqrt(r2)
c
c-----  exponential covariance function -----
c
         if(indec.eq.1) covy=dexp(-r) 
c
c-----  Gaussian covariance function -----
c

         if(indec.eq.2) covy=dexp(-r2)

       if(indec.eq.3)  then
c
c-----  Whittle covariance function 
c       (Mizell et al. WRR 18(4), 1053-1067, 1982) -----
c
       r=r*pi/(dble(2.0)*scl)
       covy=r*bessk1(r)
              
       end if 



        if(indec.eq.4) then 
c
c----- Mizell covariance function (Type B) -----
c     (Mizell et al. WRR 18(4), 1053-1067, 1982) -----
c

         r=dble(3.)*pi*r/(dble(16.0)*scl)
         covy=(dble(1.0)+r*r/dble(8.))*r*bessk1(r)-r*r*bessk0(r)

           end if
c#1       
      else
         covy=dble(1.0)
c#1
      end if
      
      if(indec.eq.5) then
         if(r.gt.1.0d-10) then
            covy=c0*r2**omega2
         else
            covy=0.0d0
         end if  
      end if   
            
      return
      end
   
