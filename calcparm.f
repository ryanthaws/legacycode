      subroutine CALCPARM(j, parent) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      parameter (P=2, QN=2)
      include 'params.f'
      dimension parent(nparmax, 200)

c      do i=1, 18
c      write(6, *) parent(i, j), i, j  
c      enddo    
      Zs=parent(7, j)
      Zp=parent(8, j)
      Gpp=parent(12, j)
      Gp2=parent(13 ,j)
      Hsp=parent(14, j)
      Gss=parent(10, j)
      
C         
         P2=P**2
         P4=P**4
         HPP=0.5D0*(GPP-GP2)
         DD=(2.D0*QN+1)*(4.D0*ZS*ZP)**(QN+0.5D0)/(ZS+ZP)
     1   **(2.D0*QN+2)/SQRT(3.D0)
         QQ=SQRT((4.D0*QN*QN+6.D0*QN+2.D0)/20.D0)/ZP
C     CALCULATE ADDITIVE TERMS, IN ATOMIC UNITS.
         GDD1= sign((P2*abs(HSP)/(27.21* 4.*DD**2))**(1./3.), HSP)
         GQQ= (P4*HPP/(27.21d0*48.d0*QQ**4))
         GQQ= sign(abs(GQQ)**(1.d0/5.d0), GQQ)
         D1=GDD1
         D2=GDD1+0.04
         Q1=GQQ
         Q2=GQQ+0.04
         del=10.d0
         DO while (del.gt.0.0001d0)
            DF=D2-D1
            HSP1= 2.*D1 - 2./SQRT(4.*DD**2+1./D1**2)
            HSP2= 2.*D2 - 2./SQRT(4.*DD**2+1./D2**2)
            HSP1= HSP1/P2
            HSP2= HSP2/P2
            D3= D1 + DF*(HSP/27.21-HSP1)/(HSP2-HSP1)
            D1= D2
            D2= D3
            del=abs(d1-d2)
c            write(6, *) d2
         enddo
         del=10.d0
         DO while (del.gt.0.0001d0)
            QF=Q2-Q1
            HPP1= 4.*Q1 - 8./SQRT(4.*QQ**2+1./Q1**2)
     1            + 4./SQRT(8.*QQ**2+1./Q1**2)
            HPP2= 4.*Q2 - 8./SQRT(4.*QQ**2+1./Q2**2)
     1            + 4./SQRT(8.*QQ**2+1./Q2**2)
            HPP1= HPP1/P4
            HPP2= HPP2/P4
            Q3= Q1 + QF*(HPP/27.21-HPP1)/(HPP2-HPP1)
            Q1= Q2
            Q2= Q3
c            write(6, *) q2
            del=abs(q1-q2)            
         enddo
         am=Gss/27.21
         ad=d2
         aq=q2
c         write(6, *) dd, qq, am, ad, aq
         open(10, file="pm3.prm")
         write(10, *) '0'
         do i=1, 2
           write(10, *) parent(i, j)
         enddo
         write(10, *) '7'
         do i=3, 20
           write(10, *) parent(i, j)
         enddo
         write(10, *) dd
         write(10, *) qq
         write(10, *) am 
         write(10, *) ad
         write(10, *) aq
         close(10)
         end
