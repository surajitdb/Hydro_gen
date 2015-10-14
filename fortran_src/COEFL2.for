c program coef.f
c
c
        subroutine coefl2(itype,cov,u1,vcx)




C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
 
        parameter(icond=5,iwork=icond*(icond+1)/2)

       DIMENSION ap(iwork),GAMA(icond),U1(*),cov(10)
       integer iwk(iwork)

        data zero,one/0.0d0,1.0d0/
C
C
C
C
C-- compute the kriging matrix ----
c
      n1=4
      n=5
      ap(1)=cov(1)
      ap(2)=cov(2)
      ap(3)=cov(1)
      ap(4)=cov(3)
      ap(5)=cov(4)
      ap(6)=cov(1)
      ap(7)=cov(4)
      ap(8)=cov(5)
      ap(9)=cov(6)
      ap(10)=cov(1)
 
      ij=7
      do i=1,4
      gama(i)=cov(ij)
      u1(i)=gama(i)
      ij=ij+1
      end do


       nd=n1

       if(itype.eq.5) then
       do i=11,14
       ap(i)=one
       end do
       ap(15)=zero
       gama(n)=one
       u1(n)=gama(n)
       nd=n

       end if
c
c--------------------------------------------------------------------
c         factors a double precision symmetric matrix stored in
c         packed form by elimination with symmetric pivoting.
c----------------------------------------------------------------------
c
      info=0
      call dspfa(ap,nd,iwk,info)

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

      call dspsl(ap,nd,iwk,u1)


        VCX=dble(0.)

        DO I=1,nd
        VCX=VCX+U1(i)*GAMA(i)
        end do

c         write(*,*)'conditional variance:',vcx

	RETURN
	END

