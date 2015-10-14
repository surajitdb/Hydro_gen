c program  coefcy.f
c
c
      subroutine krig(ilevel,n,idm,jfin,itype,cov,u1,vcx)


      IMPLICIT DOUBLE PRECISION(A-G,O-Z)
 
c      include 'dimension.blk' 
      parameter(icond=1670,iwork=icond*(icond+1)/2)
    

      DIMENSION ap(iwork),GAMA(icond),U1(*),cov(idm,*)
        integer iwk(iwork)

        data zero,one/0.0d0,1.0d0/



C
C
C-- compute the kriging matrix ----
c

	N1=N-1
        k=0
	DO  j=1,N1
            i2=int((j-1)/jfin)
              do i=1,j
                 i1=int((I-1)/jfin)
                 ilag=i2-i1+1
                 jlag=iabs(j-i-(ilag-1)*jfin)+1
                 k=k+1
                 ap(k)=cov(ilag,jlag)
              end do
        end do


       
c
C
c---- compute the known term vector 
c
      if(ilevel.eq.0) then

          i2=int((n-1)/jfin)
          DO  I=1,N1
          i1=int((i-1)/jfin)
          ilag=i2-i1+1
          jlag=iabs(n-i-(ilag-1)*jfin)+1
          GAMA(i)=cov(ilag,jlag)
          u1(i)=gama(i)
          end do
       else

          do i=1,n1
          gama(i)=cov(3,i)
          u1(i)=gama(i)
          end do
       end if


      
       ntot=n1

       if(itype.eq.5) then
       do i=1,n1
          k=k+1
          ap(k)=one
       end do
       ap(k+1)=zero
       gama(n)=one
       u1(n)=gama(n)
       ntot=n
       end if
c
c--------------------------------------------------------------------
c         factors a double precision symmetric matrix stored in
c         packed form by elimination with symmetric pivoting.
c----------------------------------------------------------------------
c
      info=0
      call dspfa(ap,ntot,iwk,info)

       if(info.ne.0) then
           write(*,*)'error in matrix factorization'
           stop
        end if
c
c---------------------------------------------------------------------
c     solves the double precision symmetric system
c     a * x = b
c     using the factors computed by dspfa.
c---------------------------------------------------------------------

      call dspsl(ap,ntot,iwk,u1)

      VCX=dble(0.)


      DO I=1,ntot
      VCX=VCX+U1(i)*GAMA(i)
        end do
c       write(*,*)nd,vcx
c
C
	RETURN
	END
