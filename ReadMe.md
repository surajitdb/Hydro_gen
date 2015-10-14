# Hydro_gen

Hydro_gen is a computer code for generating two-dimensional space random functions with an assigned covariance structure. The original code is written in Ansi Fortran 77, however this repository will collect also the porting of the original code in the OMS 3.0 framework.

## Authors

               ALBERTO BELLIN^1 AND YORAM RUBIN^2

       1:   Dipartimento di Ingegneria Civile ed Ambientale
            Universita' di Trento
            via Mesiano, 77, I-38050 Trento, Italy
            phone:  +39 461 882620
            fax:    +39 461 882672
            e-mail: Alberto.Bellin@ing.unitn.it

       2:   Department of Civil Engineering
            University of California, Berkeley
            Berkeley, CA 94720, USA
            phone:  +1 510 642 2282
            fax:    +1 510 642 7476
            e-mail: rubin@ce.berkeley.edu


Software developed by Alberto Bellin  
Universita' di Trento,
Dipartimento di Ingegneria Civile ed Ambientale,
38050-I Mesiano di Povo, TRENTO

Official web-page [Hydrogen: a random field generator](http://www.ing.unitn.it/~bellin/frames/hydrogen.php)

## Copyright

Copyright: Alberto Bellin and Yoram Rubin.

please cite the following paper in papers or reports that use the present software:

> Bellin A., Y. Rubin, **Hydro_gen: A new random field generator for correlated properties**, Stochastic Hydrology and Hydraulics, 10(4), 1996.



Permission is  granted to anyone to use and modify this packages provided that:

1. the authors are acknowledged by citing the abofe referenced paper;
2. the use in any kind of research or job will be cited in the relative papers or reports;
3. the use of the package is under the user responsability NO WARRANTY is given concerning bugs and errors.
4. the use or distribution must be free of charge.
5. the package uses the following libraries:
    * LINPACK by J. J. Dongarra, J. R. Bunch, C. B. Moler e G.W. Stewart, for the linear system solution
    * BLAS, for linear algebra
    * RANLIB by Barry W. Brown and James Lovato, Department of Biomathematics, Box 237 the University of Texas, M.D. Anderson Cancer Center 1515 Holcombe Boulevard, Huston, TX 77030, for the generation of independent normally distributed random numbers.
    * Numerical Recipes by W. H. Press, B. P. Flannery, S. A. Teukolsky, W. T. Vetterling, for the function computing the Bessel Function

Copyright conditions of the above referenced libraries are extended to hydro_gen.

Bug reports and hints are welcomed to the following e-mail address:
Alberto.Bellin@ing.unitn.it