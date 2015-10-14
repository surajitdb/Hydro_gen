c  file comb.f
c
c
c-------------------------------------------------------------------
c Author: Alberto Bellin  
c Universita' di Trento,
c Dipartimento di Ingegneria Civile ed Ambientale,
c 38050-I Mesiano di Povo, TRENTO

c Internet: alberto@itnca1.ing.unitn.it 
c           Alberto.Bellin@ing.unitn.it 


c Summary: It performs the linear combinations of selected field values
c          to compute the conditional mean and generates the random
c          deviate. 


c Package Version: 2.0  April 1997 

c Copyright: Alberto Bellin and Yoram Rubin.
c
c please cite the following report in papers or reports that use
c the present software:
c
c 1) Yoram Rubin and Alberto Bellin,  Hydro_gen: A New Random Field
c    Generator for Correlated Properties, Geotecnical Engineering Report
c    No. UCB/GT/94-04
c
c Permission is  granted to anyone to use and modify this packages provided
c that:
c i) the author is acknowledged;
c ii) the use in any kind of research or job will be cited in the relative
c papers
c or reports;
c iii) anyone who uses the packages is under his/her own responsibility.
c NO WARRANTY is given it is free of bugs and errors.
c iv) The use or distribution must be free of charge.

c Bug reports and hints are welcomed to the above e-mail address.
c----------------------------------------------------------------------------


      subroutine comb(idm,idiv,i,j,iiniz,jiniz,jfin,irv,sc1,u1,r,cond)


      implicit double precision(a-g,o-z)

      dimension cond(idm,*),u1(*),r(*)  
      
             ii0=idiv*i-idiv+1
             jj0=idiv*j-idiv+1 
 
      
             k=0
              jjfin=jfin
             do  ii=iiniz,i
               if(ii.eq.i) jjfin=j-1
               do  j1=jiniz,jjfin
                 k=k+1
                 iic=idiv*ii-idiv+1
                 j1c=idiv*j1-idiv+1 

 
                 cond(ii0,jj0)=cond(ii0,jj0)+cond(iic,j1c)*u1(k)
               end do
             end do

              irv=irv+1
             cond(ii0,jj0)=cond(ii0,jj0)+sc1*r(irv)
      return
      end
