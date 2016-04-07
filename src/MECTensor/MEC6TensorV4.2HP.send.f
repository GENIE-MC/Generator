ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  2ph hadronic tensor
c  as in Phys.Rev. C83 (2011) 045501
c    Takes q in the z direction
c  INPUT: from terminal
c  OUTPUT: File 'HadTensor.dat' with a grid of data
c                htfull(5,120,120)
c                1st index (1-5) components of the tensor (REAL)
c                2nd q0
c                3rd q
c USES: startnuclei --initialization
c       inicia      --initialization
c       HadTensor(q0,q,HTpart) returns HTpart(5), the 
c                 5 components of the Hadronic tensor
c
c COMMENTS: in this example, the head program generates a grid
c           with a binning of 10 MeV for both q0 and q 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Version 6 of 2ph nu,nubar 
c
c uses symmetry properties of Had. Tensor. See below!
c See e.g. Eq. (6) of Phys. Rev. C 73, 025504 (2006)
c
c   for PARALLEL openmp. It compiles with "gfortran -fopenmp"
ccccccccccccccccccccccccccccccccccccccccccccccccccccc

       program main
       implicit complex*16 (c)
       implicit real*8 (a,b,d-h,o-z)
       character*16 name
       dimension htfull(5,120,120)
       dimension htpart(5)
       common /nucleus/rho0,a,th,na
       common /options/iopt,isospin

c Chooses Nucleus
c Nuclei are 12C, 16O or  40Ca
       write(6,*)'12C, 16O or 40Ca? Choose 12,16,28,40,56,112,208'
       read(5,*)na
c options
       write(6,*)'options: Full=1,Only delta=2,Nodelta=3'
       read(5,*)iopt
c options isospin
       write(6,*)'Isospin options: Full=1',
     &        ' Only Equal Nucleons in Final State=2'
       read(5,*)isospin
c
       name='HadTensor.dat'
c
       call startnuclei        !! selects the nucleus and initializes 
       call inicia             !! some defs
c
       do iq0=1,120
         write(6,*)iq0
! parallel
        nchunk=1
!$OMP PARALLEL SHARED(NTHREADS,nCHUNK,iq0,q0,htfull)
!$OMP&  firstPRIVATE(ilorentz,iq,q,htpart)
!$OMP  DO SCHEDULE(DYNAMIC,nCHUNK)  
       do iq=1,120
        q0=0.01*iq0
        q=0.01*iq
c only q > q0. Leptonic tensor implies that
         if(iq.ge.iq0)then
          call HadTensor(q0,q,HTpart)
           do ilorentz=1,5
           htfull(ilorentz,iq0,iq)=htpart(ilorentz)
           enddo
         else
           do ilorentz=1,5
           htfull(ilorentz,iq0,iq)=0.d0
           enddo
         endif
       enddo
!$OMP END  DO
!$OMP END PARALLEL                 
       enddo
        open(unit=27,file=name)
        write(27,*)htfull
       end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c uses symmetry properties of Had. Tensor. for the case where
c LAB system and q=(q0,0,0,qz)
c only a few (5) components are needed
c
c See e.g. Eq. (6) of Phys. Rev. C 73, 025504 (2006)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine HadTensor(q0,qin,HTpart)
       implicit complex*16 (c)
       implicit real*8 (a,b,d-h,o-z)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      real*8 c3a,c4a,c5a,c6a,c3v,c4v,c5v,c6v
      dimension qd(0:3),htpart(5)
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
!$OMP THREADPRIVATE(/formfactors/)
       common /options/iopt,isospin

c
         qd(0)=q0
         qd(1)=0.
         qd(2)=0.
         qd(3)=qin
c        
         c3aff=c3a(qd)
         c4aff=c4a(qd)
         c5aff=c5a(qd)
         c6aff=c6a(qd)
         c3vff=c3v(qd)
         c4vff=c4v(qd)
         c5vff=c5v(qd)
         c6vff=c6v(qd)
         if(iopt.eq.3)then
c no Delta
          c3aff=0.d0
          c4aff=0.d0
          c5aff=0.d0
          c6aff=0.d0
          c3vff=0.d0
          c4vff=0.d0
          c5vff=0.d0
          c6vff=0.d0
         endif

       do icounter=1,5
c  
        if(icounter.eq.1)then   !00 component
         mu=0
         nu=0
        elseif(icounter.eq.2)then  !03 component
         mu=0
         nu=3
        elseif(icounter.eq.3)then  !11 component
         mu=1
         nu=1
        elseif(icounter.eq.4)then  !12 component
         mu=1
         nu=2
        elseif(icounter.eq.5)then  !33 component
         mu=3
         nu=3
        endif
            if(iopt.eq.1.or.iopt.eq.3)then
	     call HadronicTensorDelta(mu,nu,qd,cWD)
             call HadronicTensor2body(mu,nu,qd,cW2)
             call HadronicTensorDifferentBubbles(mu,nu,qd,cWDiff)
             call HadronicTensorRho(mu,nu,qd,cWrho)
             cWtot=cWD+cW2+cWDiff+cWrho
            else if(iopt.eq.2)then
c only Delta
	     call HadronicTensorDelta(mu,nu,qd,cWD)
             cWtot=cWD
            endif
c
           htpart(icounter)=real(cwtot)     
           if(icounter.eq.4)htpart(icounter)=aimag(cwtot)     

      enddo
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc matrix initialization
cc
      subroutine inicia
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
c
      dimension cga(4,4,0:3),cga5(4,4)
      dimension cIdent(4,4),g(0:3,0:3)
      dimension cMdeltaMatrix(4,4),cMidentity(4,4)

      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
   
      data ci/(0.d0,1.d0)/

      data cga/1.d0,0.d0,0.d0,0.d0,   
     &         0.d0,1.d0,0.d0,0.d0,   
     &         0.d0,0.d0,-1.d0,0.d0, 
     &         0.d0,0.d0,0.d0,-1.d0,
     &
     &         0.d0,0.d0,0.d0,-1.d0,
     &         0.d0,0.d0,-1.d0,0.d0,
     &         0.d0,1.d0,0.d0,0.d0,
     &         1.d0,0.d0,0.d0,0.d0,
     &     
     &         0.d0,0.d0,0.d0,(0.d0,-1.d0),
     &         0.d0,0.d0,(0.d0,1.d0),0.d0,
     &         0.d0,(0.d0,1.d0),0.d0,0.d0,
     &         (0.d0,-1.d0),0.d0,0.d0,0.d0,
     &     
     &         0.d0,0.d0,-1.d0,0.d0,
     &         0.d0,0.d0,0.d0,1.d0,
     &         1.d0,0.d0,0.d0,0.d0,
     &         0.d0,-1.d0,0.d0,0.d0/

      data cga5/0.d0,0.d0,1.d0,0.d0,
     &          0.d0,0.d0,0.d0,1.d0,
     &          1.d0,0.d0,0.d0,0.d0,
     &          0.d0,1.d0,0.d0,0.d0/

      data cIdent/1.d0,0.d0,0.d0,0.d0,
     &            0.d0,1.d0,0.d0,0.d0,
     &            0.d0,0.d0,1.d0,0.d0,
     &            0.d0,0.d0,0.d0,1.d0/

      data cMidentity/0.94d0,0.d0,0.d0,0.d0,
     &                0.d0,0.94d0,0.d0,0.d0,
     &                0.d0,0.d0,0.94d0,0.d0,
     &                0.d0,0.d0,0.d0,0.94d0/

      data cMdeltaMatrix/1.232d0,0.d0,0.d0,0.d0,
     &                   0.d0,1.232d0,0.d0,0.d0,
     &                   0.d0,0.d0,1.232d0,0.d0,
     &                   0.d0,0.d0,0.d0,1.232d0/

      data g/1.d0,0.d0,0.d0,0.d0,
     &       0.d0,-1.d0,0.d0,0.d0,
     &       0.d0,0.d0,-1.d0,0.d0,
     &       0.d0,0.d0,0.d0,-1.d0/

 
c init. gammamu*gammanu matrix for Delta
      call gammunus
c       
      return
      end


****************************************************************
****************************************************************
****************************************************************
c The matrices must be dimensioned in the main program
c This is for calculating the matricial product of two matrices cA,cB
      subroutine MatrixProduct(cA,cB,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cA(4,4),cB(4,4),cC(4,4)

      do i=1,4
         do j=1,4
            csuma=0.d0
            do k=1,4
               csuma=csuma+cA(i,k)*cB(k,j)
            enddo
            cC(i,j)=csuma
         enddo
      enddo
      return
      end
c The matrices must be dimensioned in the main program
c This is for calculating the matricial product of matrices cA,gamma5
      subroutine MatrixProductg5(cA,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cA(4,4),cC(4,4)

      do i=1,4
         do j=1,2
          cc(i,j)=cA(i,j+2)  
          cc(i,j+2)=cA(i,j)  
        enddo
      enddo
      return
      end
c The matrices must be dimensioned in the main program
c This is for calculating the matricial product of matrices gamma5,ca
      subroutine MatrixProductg5x(cA,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cA(4,4),cC(4,4)

      do i=1,2
         do j=1,4
          cc(i,j)=cA(i+2,j)  
          cc(i+2,j)=cA(i,j)  
        enddo
      enddo
      return
      end
****************************************************************
****************************************************************
****************************************************************
c The matrices must be dimensioned in the main program
c This subroutine is for calculating the matricial sum
      subroutine MatrixSum(cA,cB,cC)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cA(4,4),cB(4,4),cC(4,4)

      do i=1,4
         do j=1,4
            cC(i,j)=cA(i,j)+cB(i,j)
         enddo
      enddo

      return
      end
*****************************************************************
*****************************************************************
*****************************************************************
c The matrices must be dimensioned in the main program
c This subroutine is for calculating the matricial substraction
      subroutine MatrixResta(cA,cB,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cA(4,4),cB(4,4),cC(4,4)

      do i=1,4
         do j=1,4
            cC(i,j)=cA(i,j)-cB(i,j)
         enddo
      enddo

      return
      end
*****************************************************************
*****************************************************************
*****************************************************************
c The matrices must be dimensioned in the main program
c This subroutine is for calculating the product of a real scalar 
c by a matrix
      subroutine ScalarProductMatrix(a,cB,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cB(4,4),cC(4,4)

      do i=1,4
         do j=1,4
            cC(i,j)=a*cB(i,j)
         enddo
      enddo

      return
      end
*****************************************************************
*****************************************************************
*****************************************************************
c The matrices must be dimensioned in the main program
c This subroutine is for calculating the product of a complex scalar
c by a matrix
      subroutine ComplexScalarProductMatrix(ca,cB,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cB(4,4),cC(4,4)

      do i=1,4
         do j=1,4
            cC(i,j)=ca*cB(i,j)
         enddo
      enddo

      return
      end
*****************************************************************
*****************************************************************
*****************************************************************
c This is subroutine Dagger to calculate the hermitic conjugate of a
c matrix
      subroutine Dagger(cA,cAdagger)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cA(4,4),cAdagger(4,4)

      do i=1,4
         do j=1,4
            cAdagger(i,j)=dconjg(cA(j,i))
         enddo
      enddo

      return
      end

********************************************************************
c  use Eliecer subroutine slash
      subroutine slash(p,upsl)         
      implicit real*8 (a-h,o-t,v-z)
      implicit complex*16 (u)
      dimension p(0:3)
      dimension upsl(4,4)
       uy=(0.d0,1.d0)
      
      upsl(1,1)=p(0)
      upsl(1,2)=0
      upsl(1,3)=-p(3)
      upsl(1,4)=uy*p(2)-p(1)
      upsl(2,1)=0
      upsl(2,2)=p(0)
      upsl(2,3)=-uy*p(2)-p(1)
      upsl(2,4)=p(3)
      upsl(3,1)=p(3)
      upsl(3,2)=-uy*p(2)+p(1)
      upsl(3,3)=-p(0)
      upsl(3,4)=0
      upsl(4,1)=uy*p(2)+p(1)
      upsl(4,2)=-p(3)
      upsl(4,3)=0
      upsl(4,4)=-p(0)

      return
      end




********************************************************************
c This is subroutine LorentzScalarProduct for calculating the
c scalar product of two four-vectors
      subroutine LorentzScalarProduct(p,q,SP)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),g(0:3,0:3)
      common /metrictensor/ g

      SP=0.d0
      
      do mu=0,3
         nu=mu
            SP=SP+g(mu,nu)*p(mu)*q(nu)
      enddo

      return
      end
********************************************************************
********************************************************************
********************************************************************
c This is subroutine SumFourVectors to calculate the sum of two 
c four-vectors
      subroutine SumFourVectors(p,q,pplusq)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),pplusq(0:3)

      do mu=0,3
         pplusq(mu)=p(mu)+q(mu)
      enddo

      return
      end
*********************************************************************
*********************************************************************
*********************************************************************
c This is subroutine RestaFourVectors to calculate the substraction
c of two four-vectors
      subroutine RestaFourVectors(cp,cq,cpminusq)
      implicit real*8 (a-h,o-z)
      dimension cp(0:3),cq(0:3),cpminusq(0:3)

      do mu=0,3
         cpminusq(mu)=cp(mu)-cq(mu)
      enddo

      return
      end

*********************************************************************
*********************************************************************
*********************************************************************
c These are the vector form factors for the N-Delta transition matrix
c element. They depend on q^2. 'cq' is the four-momentum
      real*8 function C3V(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)
cparameter Mv in GeV	
      parameter(xMv=0.84d0)

c First we call the subroutine LorentzScalarProduct to calculate q^2
      call LorentzScalarProduct(cq,cq,cq2)

      C3V=(2.13d0/(1.d0-cq2/xMv**2)**2)*(1.d0/(1.d0-cq2/(4.d0*
     &     xMv**2)))

      return
      end
*********************************************************************
*********************************************************************
*********************************************************************
      real*8 function C4V(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)
cparameter Mv in GeV	
      parameter(xMv=0.84d0)

c First we call the subroutine LorentzScalarProduct to calculate q^2
      call LorentzScalarProduct(cq,cq,cq2)

      C4V=(-1.51d0/(1.d0-cq2/xMv**2)**2)*(1.d0/(1.d0-cq2/(4.d0*
     &     xMv**2)))

      return
      end
********************************************************************
********************************************************************
********************************************************************
      real*8 function C5V(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)
cparameter Mv in GeV	
      parameter(xMv=0.84d0)

c First we call the subroutine LorentzScalarProduct to calculate q^2
      call LorentzScalarProduct(cq,cq,cq2)

      C5V=(0.48d0/(1.d0-cq2/xMv**2)**2)*(1.d0/(1.d0-cq2/(0.776d0*
     &     xMv**2)))

      return
      end
**********************************************************************
**********************************************************************
**********************************************************************
      real*8 function C6V(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)

      C6V=0.d0

      return
      end
***********************************************************************
***********************************************************************
***********************************************************************
c These are the axial form factors for the N-Delta transition matrix
c element. They depend on q^2. 'cq' is the four-momentum

      real*8 function C3A(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)

      C3A=0.d0

      return
      end
*************************************************************************
*************************************************************************
*************************************************************************
      real*8 function C4A(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)

      C4A=-1.d0/4.d0*C5A(cq)

      return
      end
***********************************************************************
***********************************************************************
***********************************************************************
      real*8 function C5A(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)
cparameter MADelta in GeV
      parameter(xMADelta=0.93d0)

c call to the subroutine LorentzScalarProduct to calculate q^2
      call LorentzScalarProduct(cq,cq,cq2)

      C5A=(1.00d0/(1.d0-cq2/xMADelta**2)**2)*(1.d0/(1.d0-cq2/(3.d0*
     &     xMADelta**2)))

      return
      end
************************************************************************
************************************************************************
************************************************************************
      real*8 function C6A(cq)
      implicit real*8 (a-h,o-z)
      dimension cq(0:3)
cparameters M(nucleon mass) and mpion in GeV
      parameter(xM=0.94d0,xmpion=0.13957d0)

c call to the subroutine LorentzScalarProduct to calculate q^2
      call LorentzScalarProduct(cq,cq,cq2)

      C6A=C5A(cq)*xM**2/(xmpion**2-cq2)

      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
c  define the subroutine to calculate the vector
c part of the W+n->Delta+ vertex. We don't include Cabibbo angle in this
c definition. We must supply two Lorentz indices, the four-momentum 'p'
c of the neutron, the four momentum 'q' of the W+. And we receive as
c output variable a matrix 'cV' which depends on the two Lorentz indices.

      subroutine Vector(ialpha,mu,p,q,pdelta,gamu,cqslashx,cV)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      dimension p(0:3),q(0:3),cV(4,4),cpdeltaslash(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension pDelta(0:3),cqslash(4,4),cAUX1(4,4),cAUX2(4,4)
      dimension cAUX3(4,4),cAUX4(4,4),cAUX5(4,4),cAUX6(4,4)
      dimension cAUX7(4,4),cMdeltaMatrix(4,4),cMidentity(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
!$OMP THREADPRIVATE(/formfactors/)
      parameter(xM=0.94d0)
      dimension cqslashx(4,4)

c First we must define pDelta=p+q
c      call SumFourVectors(p,q,pDelta)
c      call slash(q,cqslash)
	
c work on the first term 
c	(g^{\alpha\mu}qslash-q^{\alpha}gamma^{\mu})
      xconstant=C3Vff/xM
      xconstant3=q(ialpha)*xconstant
      if(ialpha.eq.mu)then
       xconstant=gamu*xconstant
       call ScalarProductMatrix(xconstant,cqslashx,cAUX1)
       call ScalarProductMatrix(xconstant3,cga(1,1,mu),cAUX2)
       call MatrixResta(cAUX1,cAUX2,cAUX4)
      else
        call ScalarProductMatrix(-xconstant3,cga(1,1,mu),cAUX4)
      endif
c So we have the first term stored in the matrix CAUX4. We cannot
c touch this matrix.

c Let's work on the second term
c	g^{\alpha\mu}q.pDelta-q^{\alpha}pDelta^{\mu}
      call LorentzScalarProduct(q,pDelta,qpDelta)
      xconstant0=(C4Vff/xM**2)*(gamu*qpDelta-q(ialpha)*
     &     pDelta(mu))

c Let's work on the third term
c	g^{\alpha\mu}q.p-q^{\alpha}p^{\mu}
      call LorentzScalarProduct(q,p,qp)
      xconstant1=(C5Vff/xM**2)*(gamu*qp-q(ialpha)*p(mu))

c Let's work on the fourth term
c	g^{\alpha\mu}
      xconstant2=C6Vff*g(ialpha,mu)
      call ScalarProductMatrix(xconstant0+xconstant1+xconstant2,
     &     cIdent,cAUX7)
c So we have the fourth term stored in matrix CAUX7. We cannot touch
c this matrix.

c must sum all the matricial terms
      call MatrixSum(cAUX4,cAUX7,cAUX3)
c We have all the matricial terms summed stored in cAUX3.
c have to multiply this cAUX3 by gamma5
      call MatrixProductg5(cAUX3,cv)

      return
      end
*********************************************************************
*********************************************************************
*********************************************************************
c  define the subroutine to calculate the axial
c part of the W+n->Delta+ vertex. We don't include Cabibbo angle in this
c definition. We must supply two Lorentz indices, the four-momentum 'p'
c of the neutron, the four momentum 'q' of the W+. And we receive as
c output variable a matrix 'cAxial' which depends on the two Lorentz 
c indices.

      subroutine Axial(ialpha,mu,p,q,pdelta,gamu,cqslashx,cAxial)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      dimension p(0:3),q(0:3),cAxial(4,4),cpdeltaslash(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension pDelta(0:3),cqslash(4,4),cAUX1(4,4),cAUX2(4,4)
      dimension cAUX3(4,4),cAUX4(4,4),cAUX5(4,4),cAUX6(4,4)
      dimension cAUX7(4,4),cMdeltaMatrix(4,4),cMidentity(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
!$OMP THREADPRIVATE(/formfactors/)
      parameter(xM=0.94d0)
      dimension cqslashx(4,4)

c First we must define pDelta=p+q
c      call SumFourVectors(p,q,pDelta)
c      call slash(q,cqslash)
c work on the first term 
c	(g^{\alpha\mu}qslash-q^{\alpha}gamma^{\mu})
      xconstant=C3Aff/xM
      xconstant3=q(ialpha)*xconstant
      if(ialpha.eq.mu)then
       xconstant=gamu*xconstant
       call ScalarProductMatrix(xconstant,cqslashx,cAUX1)
       call ScalarProductMatrix(xconstant3,cga(1,1,mu),cAUX2)
       call MatrixResta(cAUX1,cAUX2,cAUX4)
      else
        call ScalarProductMatrix(-xconstant3,cga(1,1,mu),cAUX4)
      endif
c So we have the first term stored in the matrix CAUX4. We cannot
c touch this matrix.

c Let's work on the second term
c	g^{\alpha\mu}q.pDelta-q^{\alpha}pDelta^{\mu}
      call LorentzScalarProduct(q,pDelta,qpDelta)
      xconstant0=(C4Aff/xM**2)*(gamu*qpDelta-q(ialpha)*
     &     pDelta(mu))
c Let's work on the third term
c	g^{\alpha\mu}
      xconstant1=C5Aff*gamu
c Let's work on the fourth term
c	q^{\alpha}q^{\mu}
      xconstant2=(C6Aff/xM**2)*q(ialpha)*q(mu)
      call ScalarProductMatrix(xconstant0+xconstant1+xconstant2,
     & cIdent,cAUX7)
c So we have the fourth term stored in matrix CAUX7. We cannot touch
c this matrix.

c must sum all the matricial terms
      call MatrixSum(cAUX4,cAUX7,cAxial)

      return
      end
*********************************************************************
*********************************************************************
*********************************************************************
c  define the subroutine Vertex to calculate 
c the sum of the vector and axial part of the W+n-->Delta+ vertex.
c It is called \Gamma^{\alpha\mu}(p,q) in the paper by Eliecer Hernández,
c Juan Nieves and Manolo Valverde (arXiv:hep-ph/0701149v2)
c We must supply two Lorentz indices, the four-momentum 'p'
c of the neutron, the four momentum 'q' of the W+. And we receive as
c output variable a matrix 'cGamma' which depends on the two Lorentz 
c indices.
      subroutine Vertex(ialpha,mu,p,q,cqslashx,cGamma)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),cGamma(4,4)
      dimension cVector(4,4),cAxial(4,4)

c
      common /metrictensor/ g
      dimension pdelta(0:3),cqslashx(4,4),g(0:3,0:3)
c
    
c First we must define pDelta=p+q
      call SumFourVectors(p,q,pDelta)
c      call slash(q,cqslashx)
      gamu=g(ialpha,mu)
c
    
c first we call to subroutine Vector to calculate the vector part
      call Vector(ialpha,mu,p,q,pdelta,gamu,cqslashx,cVector)

c call to subroutine Axial to calculate the axial part
      call Axial(ialpha,mu,p,q,pdelta,gamu,cqslashx,cAxial)

c must sum both contributions
      call MatrixSum(cVector,cAxial,
     &     cGamma)

      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
c  define the subroutine to calculate the vector
c part of the W+n->Delta+ hermitic conjugate vertex.
c If the vector part of the W+n->Delta+ vertex is V^{\alpha\mu}, then
c the vector part of the hermitic conjugate vertex is equal to
c gamma^0*V^{\alpha\mu\dagger}*gamma^0.
c We don't include Cabibbo angle in this definition.
c We must supply two Lorentz indices, the four-momentum 'p'
c of the neutron, the four momentum 'q' of the W+. And we receive as
c output variable a matrix 'cV' which depends on the two Lorentz indices

      subroutine VectorHermiticConjugate(ialpha,mu,p,q,pdelta,
     &                                gamu,cqslashx,cV)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      dimension p(0:3),q(0:3),cV(4,4),cpdeltaslash(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension pDelta(0:3),cqslash(4,4),cAUX1(4,4),cAUX2(4,4)
      dimension cAUX3(4,4),cAUX4(4,4),cAUX5(4,4),cAUX6(4,4)
      dimension cAUX7(4,4),cMdeltaMatrix(4,4),cMidentity(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
!$OMP THREADPRIVATE(/formfactors/)
      parameter(xM=0.94d0,cosCabibbo=0.974d0)
      dimension cqslashx(4,4)

c First we must define pDelta=p+q
c      call SumFourVectors(p,q,pDelta)
c      call slash(q,cqslash)
c work on the first term 
c	(g^{\alpha\mu}qslash-q^{\alpha}gamma^{\mu})
      xconstant=C3Vff/xM
      xconstant3=q(ialpha)*xconstant
      if(ialpha.eq.mu)then
       xconstant=gamu*xconstant
       call ScalarProductMatrix(xconstant,cqslashx,cAUX1)
       call ScalarProductMatrix(xconstant3,cga(1,1,mu),cAUX2)
       call MatrixResta(cAUX1,cAUX2,cAUX4)
      else
        call ScalarProductMatrix(-xconstant3,cga(1,1,mu),cAUX4)
      endif
c So we have the first term stored in the matrix CAUX4. We cannot
c touch this matrix.

c Let's work on the second term
c	g^{\alpha\mu}q.pDelta-q^{\alpha}pDelta^{\mu}
      call LorentzScalarProduct(q,pDelta,qpDelta)
      xconstant0=-(C4Vff/xM**2)*(gamu*qpDelta-q(ialpha)*
     &     pDelta(mu))

c Let's work on the third term
c	g^{\alpha\mu}q.p-q^{\alpha}p^{\mu}
      call LorentzScalarProduct(q,p,qp)
      xconstant1=-(C5Vff/xM**2)*(gamu*qp-q(ialpha)*p(mu))

c Let's work on the fourth term
c	g^{\alpha\mu}
      xconstant2=-C6Vff*gamu
      call ScalarProductMatrix(xconstant0+xconstant1+xconstant2,
     &    cIdent,cAUX7)
c So we have the fourth term stored in matrix CAUX7. We cannot touch
c this matrix.

      call MatrixSum(cAUX4,cAUX7,cAUX3)
c We have all the matricial terms summed stored in cAUX3.
c have to multiply this cAUX3 by gamma5
      call MatrixProductg5(cAUX3,cv)

      return
      end
********************************************************************
********************************************************************
********************************************************************
c  define the subroutine to calculate the axial
c part of the W+n->Delta+ hermitic conjugate vertex.
c If the axial part of the W+n->Delta+ vertex is A^{\alpha\mu}, then
c the axial part of the hermitic conjugate vertex is equal to
c gamma^0*A^{\alpha\mu\dagger}*gamma^0.
c We don't include Cabibbo angle in this definition.
c We must supply two Lorentz indices, the four-momentum 'p'
c of the neutron, the four momentum 'q' of the W+. And we receive as
c output variable a matrix 'cAxial' which depends on the two Lorentz 
c indices
      subroutine AxialHermiticConjugate(ialpha,mu,p,q,pdelta,
     &             gamu,cqslashx,cAxial)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),cAxial(4,4)
c
      dimension pdelta(0:3),cqslashx(4,4)
c

c As the expression is the same than A^{\alpha\mu}, we can safely
c call to the subroutine Axial

      call Axial(ialpha,mu,p,q,pdelta,gamu,cqslashx,cAxial)

      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
c  define the subroutine to calculate the hermitic
c conjugate vertex, i.e, the vertex Delta+ -->W+ + n
c If the normal vertex is called \Gamma^{\alpha\mu}(p,q) then the 
c hermitic conjugate vertex is gamma^0*\Gamma^{\alpha\mu\dagger}*gamma^0
c we must supply two Lorentz indices, the four-momentum 'p' of the
c neutron, the four-momentum 'q' of the W+. And we receive as output
c variable a matrix 'cGamma' which depends on the two Lorentz indices
      subroutine VertexHermiticConjugate(ialpha,mu,p,q,cGamma)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),cGamma(4,4)
      dimension cVector(4,4),cAxial(4,4)
c
      common /metrictensor/ g
      dimension pdelta(0:3),cqslashx(4,4),g(0:3,0:3)
c
c First we must define pDelta=p+q
      call SumFourVectors(p,q,pDelta)
      call slash(q,cqslashx)
      gamu=g(ialpha,mu)
c      
c first we call to subroutine VectorHermiticConjugate
c to calculate the vector part
      call VectorHermiticConjugate(ialpha,mu,p,q,pdelta,
     & gamu,cqslashx,cVector)

c call to subroutine AxialHermiticConjugate
c to calculate the axial part
      call AxialHermiticConjugate(ialpha,mu,p,q,pdelta,
     & gamu,cqslashx,cAxial)

c must sum both contributions
      call MatrixSum(cVector,cAxial,
     &     cGamma)

      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
c  define the subroutine DeltaProjector which is the
c spin 3/2 projection operator. It is a matrix which depends on two 
c Lorentz indices. We must supply two Lorentz indices (mu,nu), the 
c four-momentum of the delta (pdelta),and we obtain as output 
c variable a matrix 'cProjector' which depends on the two Lorentz
c indices.
      subroutine DeltaProjector(mu,nu,pdelta,caux1,cProjector)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension pdelta(0:3),cProjector(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cqslash(4,4),cMidentity(4,4)
      dimension cAUX1(4,4),cAUX2(4,4),cAUX3(4,4),cAUX4(4,4)
      dimension cAUX5(4,4),cAUX6(4,4),cAUX7(4,4),cAUX8(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /matrizqslash/cMdeltaMatrix,cMidentity
      parameter(xMdelta=1.232d0)

c (pdeltaslash+Mdelta) is stored in caux1

c must construct the second multiplicative term. First the part
c which is proportional to the Identity Matrix.
c This part is (-g^{\mu\nu}+2/3*pdelta^\mu*pdelta^\nu/Mdelta**2)
      xconstant=-g(mu,nu)+2.d0/3.d0*pdelta(mu)*pdelta(nu)/xMdelta**2
cpart (1/3*gamma^\mu*gamma^\nu)
      call gammunu3(mu,nu,caux4)   
c This part is stored in cAUX4, we cannot touch this matrix

cpart (-1/3(pdelta^\mu*gamma^\nu-pdelta^\nu*gamma^\mu)/Mdelta)
      xconstant2=-1.d0/(3.d0*xMdelta)
      xcmu=pdelta(mu)*xconstant2
      xcnu=pdelta(nu)*xconstant2

      do i=1,4
      do j=1,4
      cAUX5(i,j)=xcmu*cga(i,j,nu)-xcnu*cga(i,j,mu)
     &  +caux4(i,j)+xconstant*cIdent(i,j)
      enddo
      enddo

c have to multiply (pdeltaslash+Mdelta) by cAUX5 in order to
c obtain cProjector
      call MatrixProduct(cAUX1,cAUX5,cProjector)

      return
      end


*******************************************************************
c First we must define the Callen Lambda function, which is defined as
c lambda(x,y,z)=x**2+y**2+z**2-2*x*y-2*x*z-2*y*z

      real*8 function xlambdaCallen(x,y,z)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      xlambdaCallen=x**2+y**2+z**2-2.d0*x*y-2.d0*x*z-2.d0*y*z

      return
      end
**********************************************************************
**********************************************************************
**********************************************************************
c  define the delta resonance width in its rest 
c frame. It depends on s=pdelta^2

      real*8 function DeltaWidth(s)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      parameter(pi=3.141592653589793d0,fstar=2.14d0,xmpion=0.13957d0,
     &     xM=0.94d0)

      if((s.lt.0.d0).or.((dsqrt(s)-xM-xmpion).lt.0.d0))then
         DeltaWidth=0.d0
      else
         DeltaWidth=1.d0/(6.d0*pi)*(fstar/xmpion)**2*xM/
     &        dsqrt(s)*(dsqrt(xlambdaCallen(s,xmpion**2,xM**2))/
     &        (2.d0*dsqrt(s)))**3
      endif
	
      return
      end
********************************************************************
********************************************************************
********************************************************************
      subroutine startnuclei
      IMPLICIT REAL*8 (A-H,O-Z)
       common /nucleus/rho0,a,th,na
c
       dimension xd(2000),fd(2000)
       pi=3.141592653589793d0
c       
       if(na.eq.12)then
c          a=2.5d0
          a=1.789d0
          th=0.55d0
        elseif(na.eq.16)then
          a=2.34d0
          th=0.55d0 
       elseif(na.eq.40)then
c       a=3.51d0
c       th=.56d0
          a=3.61d0
          th=0.56d0 
       elseif(na.eq.56)then
          a=4.106d0
          th=0.519d0
       elseif(na.eq.28)then
          a=2.86d0
          th=0.569d0
       elseif(na.eq.208)then
          a=6.89d0
          th=0.549d0
       elseif(na.eq.112)then
          a=4.608d0
          th=0.532d0
       else
          stop' sorry no data'
       endif
       
c normalize density
       rlim=5.d0*a
       np=2
       rho0=1.d0
       call DSG20R(0.d0,rlim,Np,Xd,NPp)
       do i=1,npp
          r=xd(i)
          fd(i)=r**2*density(r)
       enddo
       call DRG20R(0.d0,rlim,Np,Fd,RES)
       res=res*4.d0*pi
       rho0=dble(na)/res

       return
       end
*******************************************************************
*******************************************************************
*******************************************************************
      real*8 function Density(r)
      IMPLICIT REAL*8 (A-H,O-Z)
      common /nucleus/rho0,a,th,na
       density=rho0/(1.d0+dexp((r-a)/th))
      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
ccccccccccccccccccccccccccccc INTEGRATION cccccccccc

      SUBROUTINE DSG20R(A,B,N,X,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Y(10),X(2000)
      DATA Y/.9931285991d0,.9639719272d0,.9122344282d0,.8391169718d0,
     F .7463319064d0,.6360536807d0,.5108670019d0,.3737060887d0,
     F .2277858511d0,.0765265211d0/
      NP=20*N
      XN= N
      DINT=B - A
      DINT = DINT/XN
      DELT=DINT*0.5D0
      ORIG=A-DELT
      I1=-20
      DO 1 I=1,N
      ORIG=ORIG+DINT
      DORIG=ORIG+ORIG
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
      X(J1)=ORIG-DELT*Y(J)
 2    X(J2)=DORIG-X(J1)
 1    CONTINUE
      RETURN
      END
*********************************************************************
*********************************************************************
*********************************************************************
      SUBROUTINE DRG20R(A,B,N,CF,CRES)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION W(10),CF(2000)
      DATA W/.0176140071d0,.0406014298d0,.0626720483d0,.0832767415d0,
     F .1019301198d0,.1181945319d0,.1316886384d0,.1420961093d0,
     f .1491729864d0,.1527533871d0/
      cr=0.d0
      I1=-20
      DO 1 I=1,N
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
 2    CR=CR+W(J)*(CF(J1)+CF(J2))
 1    CONTINUE
      CRES=CR*0.5D0*(B-A)/DBLE(N)
      RETURN
      END       
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE DRG20C(A,B,N,CF,CRES)
      IMPLICIT REAL*8 (A,B,D-H,O-Z)
      IMPLICIT complex*16(C)
      DIMENSION W(10),CF(2000)
      DATA W/.0176140071d0,.0406014298d0,.0626720483d0,.0832767415d0,
     F .1019301198d0,.1181945319d0,.1316886384d0,.1420961093d0,
     f .1491729864d0,.1527533871d0/
      CR=(0.D0,0.D0)
      I1=-20
      DO 1 I=1,N
      I1=I1+20
      I2=I1+21
      DO 2 J=1,10
      J1=I1+J
      J2=I2-J
 2    CR=CR+W(J)*(CF(J1)+CF(J2))
 1    CONTINUE
      CRES=CR*0.5D0*(B-A)/DBLE(N)
      RETURN
      END
***********************************************************************
***********************************************************************
***********************************************************************
c  define the main subroutine of the program, the 
c subroutine to calculate the hadronic tensor, which has four 
c integrations. One integration in the radius, another one in the
c modulus of three momentum and two angular integrals.
c To calculate the hadronic tensor, we must supply two Lorentz indices,
c (mu,nu) and the four-momentum of the W+ (cq). And we get as an 
c output variable the hadronic tensor 'cW' which is only a complex 
c number (not a matrix).
      subroutine HadronicTensorDelta(mu,nu,q,cW)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      real*8 c3a,c4a,c5a,c6a,c3v,c4v,c5v,c6v
      dimension q(0:3),p(0:3),pdelta(0:3),pdeltakin(0:3)
      dimension rd(2000),pd(2000),thetad(2000),phid(2000)
      dimension crrd(2000),crpd(2000),crthetad(2000),crphid(2000)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),cpdeltaslash(4,4)
      dimension cqslash(4,4),cMdeltaMatrix(4,4),cMidentity(4,4)
      common /nucleus/rho0,a,th,na
      common /matrices/ cga,cga5,cIdent
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
!$OMP THREADPRIVATE(/formfactors/)
      common /options/iopt,isospin
      parameter(pi=3.141592653589793d0,xM=0.94d0,xMdelta=1.232d0,
     &     hbarc=0.19733d0,cMdelta=1.232d0,cM=0.94d0)

c CONSTANTS OF INTEREST
      xCabibbo=0.974d0

c precision in the integrals
      nr=1
      np=1
      ntheta=1
      nphi=1

c upper limits for the integrations
      rlimsup=5.d0*a
      rlimsupGeV=rlimsup/hbarc
      philimsup=2.d0*pi

c call to subroutines that make the partitions
      call DSG20R(0.d0,rlimsupGeV,nr,rd,nrr)
      call DSG20R(0.d0,pi,ntheta,thetad,nthetatheta)
      call DSG20R(0.d0,philimsup,nphi,phid,nphiphi)

c do loops to evaluate the function in the points of the partitions
      do ir=1,nrr
c rfm is the radius in fermi, rGeV is the radius in GeV^(-1).
c xkfermi is the Fermi momentum at 'r' in fm^(-1)
c xkfGeV is the Fermi momentum at 'r' in GeV
c         rfm=rd(ir)
c         rGeV=rfm/hbarc
         rGeV=rd(ir)
         rfm=rGeV*hbarc
         rhofm=density(rfm)
         xkfermi=(3.d0*pi**2/2.d0*rhofm)**(1.d0/3.d0)
         xkfGeV=xkfermi*hbarc

c make the partition in pd from 0 to kFermi
         call DSG20R(0.d0,xkfGeV,np,pd,npp)

         do ip=1,npp
            pm=pd(ip)
            energy=dsqrt(pm**2+xM**2)
c            cp(0)=energy
           p(0)=energy
	   p0kin=energy-xkfgev**2/(2.d0*xm)

            do itheta=1,nthetatheta
               theta=thetad(itheta)
               p(3)=pm*dcos(theta)

               do iphi=1,nphiphi
                  phi=phid(iphi)
                  p(1)=pm*dsin(theta)*dcos(phi)
                  p(2)=pm*dsin(theta)*dsin(phi)
c can define pdelta=p+q
                  call SumFourVectors(p,q,pdelta)
		  do ii=1,3
		  pdeltakin(ii)=pdelta(ii)
		  enddo
		  pdeltakin(0)=q(0)+p0kin
                  call LorentzScalarProduct(pdeltakin,pdeltakin,pdelta2)
                  s=pdelta2
                  if(s.lt.0.d0)then
                     s=0.d0
                  endif
                  W=dsqrt(s)
                  call hadronMathNacho(mu,nu,p,q,cTensorA)
                 
                  crphid(iphi)=cTensorA*(0.d0*PauliDeltaWidth(s,rhofm)/
     &                 2.d0-xImSelfenergyDelta2(s,rhofm))/((W+xMdelta)*
     &                 ((W-xMdelta-ReSelfenergyDelta(rhofm))**2+(
     &                 PauliDeltaWidth(s,rhofm)/2.d0-
     &                 xImSelfenergyDelta(s,rhofm))**2))
c                  write(*,*)crphid(iphi)
               enddo
               
               call DRG20C(0.d0,philimsup,nphi,crphid,cres)
               crthetad(itheta)=cres*dsin(theta)
c               write(*,*)crthetad(itheta)
            enddo
            call DRG20C(0.d0,pi,ntheta,crthetad,cres)
            crpd(ip)=cres*pm**2/energy
         enddo
c         pause
         call DRG20C(0.d0,xkfGeV,np,crpd,cres)
         crrd(ir)=cres*rGeV**2
c            write(*,*)ir,crrd(ir)
      enddo
      call DRG20C(0.d0,rlimsupGeV,nr,crrd,cres)
c This is for symmetric nuclear matter and for the processes
c p + W+ --> Delta++ and n + W+ --> Delta+ (3 + 1 == 4)
      cW=cres*4.d0*2.d0*xCabibbo**2/(2.d0*pi)**3
c bNEW isospin	
        if(isospin.eq.2)then		
         cW=cW*5./6.		
        endif		
c eNEW isospin

c If one only wants to see the pi+ in the final state, then one has to
c separate the Deltawidth in the numerator in 
c Deltawidth(Delta+ --> p + pi0) + Deltawidth(Delta+ --> n + pi+)
c and only consider the DeltaWidth to neutron + pion+, which is 1/3 of
c the total Width. This gives us a factor (3 + 1/3) instead of giving us
c the factor (3+1). So, the final result is

c	cW=cres*(3.d0+1.d0/3.d0)*2.d0/(2.d0*pi)**3

      return
      end
******************************************************************
******************************************************************
******************************************************************
      real*8 function epsfunNacho(mu,nu,p,q)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),g(0:3,0:3)
      common /metrictensor/ g
      
      if (mu.eq.nu)then
         epsfunNacho=0.d0
      else
         suma=0.d0
         do ialpha=0,3
            do ibeta=0,3
c               do ialphaprime=0,3
c                  do ibetaprime=0,3
               ialphaprime=ialpha
               ibetaprime=ibeta
               suma=suma+iepsi(mu,nu,ialpha,ibeta)*g(ialpha,
     &              ialphaprime)*g(ibeta,ibetaprime)*
     &              p(ialphaprime)*q(ibetaprime)
c                  enddo
c               enddo
            enddo
         enddo
         epsfunNacho=suma
      endif

      return
      end

*******************************************************************
*******************************************************************
*******************************************************************
      subroutine hadronmathNacho(mu,nu,p,q,cA)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),g(0:3,0:3)
      real*8 c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,c6vff
      common /metrictensor/ g
      common/formfactors/c3aff,c4aff,c5aff,c6aff,c3vff,c4vff,c5vff,
     &     c6vff
!$OMP THREADPRIVATE(/formfactors/)

      epsmunupq=epsfunNacho(mu,nu,p,q)
      call LorentzScalarProduct(p,p,p2)
      call LorentzScalarProduct(q,q,q2)
      call LorentzScalarProduct(p,q,pdotq)
      xM=0.94d0
      xMdelta=1.232d0
      gmunu=g(mu,nu)
      pmu=p(mu)
      pnu=p(nu)
      qmu=q(mu)
      qnu=q(nu)

      cA=(-4*((0,4)*c3aff*c3vff*epsmunupq*pdotq**2*xM**2 + 
     -      2*c3aff**2*gmunu*p2*pdotq**2*xM**2 + 
     -      2*c3vff**2*gmunu*p2*pdotq**2*xM**2 + 
     -      2*c3aff**2*gmunu*pdotq**3*xM**2 + 
     -      2*c3vff**2*gmunu*pdotq**3*xM**2 + 
     -      (0,8)*c3aff*c3vff*epsmunupq*pdotq*q2*xM**2 + 
     -      4*c3aff**2*gmunu*p2*pdotq*q2*xM**2 + 
     -      4*c3vff**2*gmunu*p2*pdotq*q2*xM**2 + 
     -      4*c3aff**2*gmunu*pdotq**2*q2*xM**2 + 
     -      4*c3vff**2*gmunu*pdotq**2*q2*xM**2 + 
     -     2*c3aff**2*p2*pmu*pnu*q2*xM**2 + 
     -      2*c3vff**2*p2*pmu*pnu*q2*xM**2 + 
     -     2*c3aff**2*pdotq*pmu*pnu*q2*xM**2 + 
     -      2*c3vff**2*pdotq*pmu*pnu*q2*xM**2 + 
     -      (0,4)*c3aff*c3vff*epsmunupq*q2**2*xM**2 + 
     -      2*c3aff**2*gmunu*p2*q2**2*xM**2 + 
     -     2*c3vff**2*gmunu*p2*q2**2*xM**2 + 
     -      2*c3aff**2*gmunu*pdotq*q2**2*xM**2 + 
     -      2*c3vff**2*gmunu*pdotq*q2**2*xM**2 - 
     -      2*c3aff**2*p2*pdotq*pnu*qmu*xM**2 - 
     -     2*c3vff**2*p2*pdotq*pnu*qmu*xM**2 - 
     -      2*c3aff**2*pdotq**2*pnu*qmu*xM**2 - 
     -     2*c3vff**2*pdotq**2*pnu*qmu*xM**2 - 
     -      2*c3aff**2*p2*pdotq*pmu*qnu*xM**2 - 
     -     2*c3vff**2*p2*pdotq*pmu*qnu*xM**2 - 
     -      2*c3aff**2*pdotq**2*pmu*qnu*xM**2 - 
     -     2*c3vff**2*pdotq**2*pmu*qnu*xM**2 - 
     -      4*c3aff**2*p2*pdotq*qmu*qnu*xM**2 - 
     -     4*c3vff**2*p2*pdotq*qmu*qnu*xM**2 - 
     -      4*c3aff**2*pdotq**2*qmu*qnu*xM**2 - 
     -     4*c3vff**2*pdotq**2*qmu*qnu*xM**2 - 
     -      2*c3aff**2*p2*q2*qmu*qnu*xM**2 - 
     -     2*c3vff**2*p2*q2*qmu*qnu*xM**2 - 
     -      2*c3aff**2*pdotq*q2*qmu*qnu*xM**2 - 
     -     2*c3vff**2*pdotq*q2*qmu*qnu*xM**2 - 
     -      2*c5aff**2*p2*pmu*pnu*xM**4 - 2*c6vff**2*p2*pmu*pnu*xM**4 - 
     -      2*c5aff**2*pdotq*pmu*pnu*xM**4 - 
     -     2*c6vff**2*pdotq*pmu*pnu*xM**4 - 
     -      2*c5aff**2*p2*pnu*qmu*xM**4 - 2*c6vff**2*p2*pnu*qmu*xM**4 - 
     -      2*c5aff**2*pdotq*pnu*qmu*xM**4 - 
     -     2*c6vff**2*pdotq*pnu*qmu*xM**4 - 
     -      2*c5aff**2*p2*pmu*qnu*xM**4 - 2*c6vff**2*p2*pmu*qnu*xM**4 - 
     -      2*c5aff**2*pdotq*pmu*qnu*xM**4 - 
     -     2*c6vff**2*pdotq*pmu*qnu*xM**4 - 
     -      2*c5aff**2*p2*qmu*qnu*xM**4 - 2*c6vff**2*p2*qmu*qnu*xM**4 - 
     -      2*c5aff**2*pdotq*qmu*qnu*xM**4 - 
     -     2*c6vff**2*pdotq*qmu*qnu*xM**4 + 
     -      (0,2)*c3vff*c4aff*epsmunupq*pdotq**2*xM*xMdelta + 
     -      (0,2)*c3aff*c4vff*epsmunupq*pdotq**2*xM*xMdelta + 
     -      2*c3aff*c4aff*gmunu*p2*pdotq**2*xM*xMdelta + 
     -      2*c3vff*c4vff*gmunu*p2*pdotq**2*xM*xMdelta + 
     -      2*c3aff*c4aff*gmunu*pdotq**3*xM*xMdelta + 
     -      2*c3vff*c4vff*gmunu*pdotq**3*xM*xMdelta + 
     -      (0,4)*c3vff*c4aff*epsmunupq*pdotq*q2*xM*xMdelta + 
     -      (0,4)*c3aff*c4vff*epsmunupq*pdotq*q2*xM*xMdelta + 
     -      4*c3aff*c4aff*gmunu*p2*pdotq*q2*xM*xMdelta + 
     -      4*c3vff*c4vff*gmunu*p2*pdotq*q2*xM*xMdelta + 
     -      4*c3aff*c4aff*gmunu*pdotq**2*q2*xM*xMdelta + 
     -      4*c3vff*c4vff*gmunu*pdotq**2*q2*xM*xMdelta + 
     -      2*c3aff*c4aff*p2*pmu*pnu*q2*xM*xMdelta + 
     -      2*c3vff*c4vff*p2*pmu*pnu*q2*xM*xMdelta + 
     -      2*c3aff*c4aff*pdotq*pmu*pnu*q2*xM*xMdelta + 
     -      2*c3vff*c4vff*pdotq*pmu*pnu*q2*xM*xMdelta + 
     -      (0,2)*c3vff*c4aff*epsmunupq*q2**2*xM*xMdelta + 
     -      (0,2)*c3aff*c4vff*epsmunupq*q2**2*xM*xMdelta + 
     -      2*c3aff*c4aff*gmunu*p2*q2**2*xM*xMdelta + 
     -      2*c3vff*c4vff*gmunu*p2*q2**2*xM*xMdelta + 
     -      2*c3aff*c4aff*gmunu*pdotq*q2**2*xM*xMdelta + 
     -      2*c3vff*c4vff*gmunu*pdotq*q2**2*xM*xMdelta - 
     -      2*c3aff*c4aff*p2*pdotq*pnu*qmu*xM*xMdelta - 
     -      2*c3vff*c4vff*p2*pdotq*pnu*qmu*xM*xMdelta - 
     -      2*c3aff*c4aff*pdotq**2*pnu*qmu*xM*xMdelta - 
     -      2*c3vff*c4vff*pdotq**2*pnu*qmu*xM*xMdelta - 
     -      2*c3aff*c4aff*p2*pdotq*pmu*qnu*xM*xMdelta - 
     -      2*c3vff*c4vff*p2*pdotq*pmu*qnu*xM*xMdelta - 
     -      2*c3aff*c4aff*pdotq**2*pmu*qnu*xM*xMdelta - 
     -      2*c3vff*c4vff*pdotq**2*pmu*qnu*xM*xMdelta - 
     -      4*c3aff*c4aff*p2*pdotq*qmu*qnu*xM*xMdelta - 
     -      4*c3vff*c4vff*p2*pdotq*qmu*qnu*xM*xMdelta - 
     -      4*c3aff*c4aff*pdotq**2*qmu*qnu*xM*xMdelta - 
     -      4*c3vff*c4vff*pdotq**2*qmu*qnu*xM*xMdelta - 
     -      2*c3aff*c4aff*p2*q2*qmu*qnu*xM*xMdelta - 
     -      2*c3vff*c4vff*p2*q2*qmu*qnu*xM*xMdelta - 
     -      2*c3aff*c4aff*pdotq*q2*qmu*qnu*xM*xMdelta - 
     -      2*c3vff*c4vff*pdotq*q2*qmu*qnu*xM*xMdelta + 
     -      (0,2)*c3vff*c5aff*epsmunupq*pdotq*xM**3*xMdelta + 
     -      (0,2)*c3aff*c6vff*epsmunupq*pdotq*xM**3*xMdelta + 
     -      2*c3aff*c5aff*gmunu*p2*pdotq*xM**3*xMdelta + 
     -      2*c3vff*c6vff*gmunu*p2*pdotq*xM**3*xMdelta + 
     -      2*c3aff*c5aff*gmunu*pdotq**2*xM**3*xMdelta + 
     -      2*c3vff*c6vff*gmunu*pdotq**2*xM**3*xMdelta + 
     -      (0,2)*c3vff*c5aff*epsmunupq*q2*xM**3*xMdelta + 
     -      (0,2)*c3aff*c6vff*epsmunupq*q2*xM**3*xMdelta + 
     -      2*c3aff*c5aff*gmunu*p2*q2*xM**3*xMdelta + 
     -      2*c3vff*c6vff*gmunu*p2*q2*xM**3*xMdelta + 
     -      2*c3aff*c5aff*gmunu*pdotq*q2*xM**3*xMdelta + 
     -      2*c3vff*c6vff*gmunu*pdotq*q2*xM**3*xMdelta + 
     -      2*c3aff*c5aff*pmu*pnu*q2*xM**3*xMdelta + 
     -      2*c3vff*c6vff*pmu*pnu*q2*xM**3*xMdelta - 
     -      c3aff*c5aff*p2*pnu*qmu*xM**3*xMdelta - 
     -      c3vff*c6vff*p2*pnu*qmu*xM**3*xMdelta - 
     -      2*c3aff*c5aff*pdotq*pnu*qmu*xM**3*xMdelta - 
     -      2*c3vff*c6vff*pdotq*pnu*qmu*xM**3*xMdelta + 
     -      c3aff*c5aff*pnu*q2*qmu*xM**3*xMdelta + 
     -      c3vff*c6vff*pnu*q2*qmu*xM**3*xMdelta - 
     -      c3aff*c5aff*p2*pmu*qnu*xM**3*xMdelta - 
     -      c3vff*c6vff*p2*pmu*qnu*xM**3*xMdelta - 
     -      2*c3aff*c5aff*pdotq*pmu*qnu*xM**3*xMdelta - 
     -      2*c3vff*c6vff*pdotq*pmu*qnu*xM**3*xMdelta + 
     -      c3aff*c5aff*pmu*q2*qnu*xM**3*xMdelta + 
     -      c3vff*c6vff*pmu*q2*qnu*xM**3*xMdelta - 
     -      2*c3aff*c5aff*p2*qmu*qnu*xM**3*xMdelta - 
     -      2*c3vff*c6vff*p2*qmu*qnu*xM**3*xMdelta - 
     -      4*c3aff*c5aff*pdotq*qmu*qnu*xM**3*xMdelta - 
     -      4*c3vff*c6vff*pdotq*qmu*qnu*xM**3*xMdelta - 
     -      2*c5aff**2*pmu*pnu*xM**5*xMdelta + 
     -     2*c6vff**2*pmu*pnu*xM**5*xMdelta - 
     -      2*c5aff**2*pnu*qmu*xM**5*xMdelta + 
     -     2*c6vff**2*pnu*qmu*xM**5*xMdelta - 
     -      2*c5aff**2*pmu*qnu*xM**5*xMdelta + 
     -     2*c6vff**2*pmu*qnu*xM**5*xMdelta - 
     -      2*c5aff**2*qmu*qnu*xM**5*xMdelta + 
     -     2*c6vff**2*qmu*qnu*xM**5*xMdelta - 
     -      (0,2)*c4aff*c4vff*epsmunupq*pdotq**2*xMdelta**2 + 
     -      2*c4aff**2*gmunu*p2*pdotq**2*xMdelta**2 + 
     -      2*c4vff**2*gmunu*p2*pdotq**2*xMdelta**2 + 
     -      2*c4aff**2*gmunu*pdotq**3*xMdelta**2 + 
     -      2*c4vff**2*gmunu*pdotq**3*xMdelta**2 - 
     -      (0,4)*c4aff*c4vff*epsmunupq*pdotq*q2*xMdelta**2 + 
     -      4*c4aff**2*gmunu*p2*pdotq*q2*xMdelta**2 + 
     -      4*c4vff**2*gmunu*p2*pdotq*q2*xMdelta**2 + 
     -      4*c4aff**2*gmunu*pdotq**2*q2*xMdelta**2 + 
     -      4*c4vff**2*gmunu*pdotq**2*q2*xMdelta**2 + 
     -      2*c4aff**2*p2*pmu*pnu*q2*xMdelta**2 + 
     -      2*c4vff**2*p2*pmu*pnu*q2*xMdelta**2 + 
     -      2*c4aff**2*pdotq*pmu*pnu*q2*xMdelta**2 + 
     -      2*c4vff**2*pdotq*pmu*pnu*q2*xMdelta**2 - 
     -      (0,2)*c4aff*c4vff*epsmunupq*q2**2*xMdelta**2 + 
     -      2*c4aff**2*gmunu*p2*q2**2*xMdelta**2 + 
     -      2*c4vff**2*gmunu*p2*q2**2*xMdelta**2 + 
     -      2*c4aff**2*gmunu*pdotq*q2**2*xMdelta**2 + 
     -      2*c4vff**2*gmunu*pdotq*q2**2*xMdelta**2 - 
     -      2*c4aff**2*p2*pdotq*pnu*qmu*xMdelta**2 - 
     -      2*c4vff**2*p2*pdotq*pnu*qmu*xMdelta**2 - 
     -      2*c4aff**2*pdotq**2*pnu*qmu*xMdelta**2 - 
     -      2*c4vff**2*pdotq**2*pnu*qmu*xMdelta**2 - 
     -      2*c4aff**2*p2*pdotq*pmu*qnu*xMdelta**2 - 
     -      2*c4vff**2*p2*pdotq*pmu*qnu*xMdelta**2 - 
     -      2*c4aff**2*pdotq**2*pmu*qnu*xMdelta**2 - 
     -      2*c4vff**2*pdotq**2*pmu*qnu*xMdelta**2 - 
     -      4*c4aff**2*p2*pdotq*qmu*qnu*xMdelta**2 - 
     -      4*c4vff**2*p2*pdotq*qmu*qnu*xMdelta**2 - 
     -      4*c4aff**2*pdotq**2*qmu*qnu*xMdelta**2 - 
     -      4*c4vff**2*pdotq**2*qmu*qnu*xMdelta**2 - 
     -      2*c4aff**2*p2*q2*qmu*qnu*xMdelta**2 - 
     -      2*c4vff**2*p2*q2*qmu*qnu*xMdelta**2 - 
     -      2*c4aff**2*pdotq*q2*qmu*qnu*xMdelta**2 - 
     -      2*c4vff**2*pdotq*q2*qmu*qnu*xMdelta**2 - 
     -      (0,8)*c3aff*c3vff*epsmunupq*pdotq*xM**2*xMdelta**2 - 
     -      (0,4)*c3vff*c4aff*epsmunupq*pdotq*xM**2*xMdelta**2 + 
     -      (0,4)*c3aff*c4vff*epsmunupq*pdotq*xM**2*xMdelta**2 - 
     -      (0,2)*c4vff*c5aff*epsmunupq*pdotq*xM**2*xMdelta**2 - 
     -      (0,2)*c4aff*c6vff*epsmunupq*pdotq*xM**2*xMdelta**2 + 
     -      4*c4aff*c5aff*gmunu*p2*pdotq*xM**2*xMdelta**2 + 
     -      4*c4vff*c6vff*gmunu*p2*pdotq*xM**2*xMdelta**2 + 
     -      2*c3aff**2*gmunu*pdotq**2*xM**2*xMdelta**2 + 
     -      2*c3vff**2*gmunu*pdotq**2*xM**2*xMdelta**2 + 
     -      4*c3aff*c4aff*gmunu*pdotq**2*xM**2*xMdelta**2 - 
     -      4*c3vff*c4vff*gmunu*pdotq**2*xM**2*xMdelta**2 + 
     -      4*c4aff*c5aff*gmunu*pdotq**2*xM**2*xMdelta**2 + 
     -      4*c4vff*c6vff*gmunu*pdotq**2*xM**2*xMdelta**2 - 
     -      (0,6)*c3aff*c3vff*epsmunupq*q2*xM**2*xMdelta**2 - 
     -      (0,4)*c3vff*c4aff*epsmunupq*q2*xM**2*xMdelta**2 + 
     -      (0,4)*c3aff*c4vff*epsmunupq*q2*xM**2*xMdelta**2 - 
     -      (0,2)*c4vff*c5aff*epsmunupq*q2*xM**2*xMdelta**2 - 
     -      (0,2)*c4aff*c6vff*epsmunupq*q2*xM**2*xMdelta**2 - 
     -      2*c3aff**2*gmunu*p2*q2*xM**2*xMdelta**2 - 
     -      2*c3vff**2*gmunu*p2*q2*xM**2*xMdelta**2 + 
     -      4*c4aff*c5aff*gmunu*p2*q2*xM**2*xMdelta**2 + 
     -      4*c4vff*c6vff*gmunu*p2*q2*xM**2*xMdelta**2 + 
     -      8*c3aff*c4aff*gmunu*pdotq*q2*xM**2*xMdelta**2 - 
     -      8*c3vff*c4vff*gmunu*pdotq*q2*xM**2*xMdelta**2 + 
     -      4*c4aff*c5aff*gmunu*pdotq*q2*xM**2*xMdelta**2 + 
     -      4*c4vff*c6vff*gmunu*pdotq*q2*xM**2*xMdelta**2 + 
     -      2*c3aff**2*pmu*pnu*q2*xM**2*xMdelta**2 + 
     -      2*c3vff**2*pmu*pnu*q2*xM**2*xMdelta**2 + 
     -      4*c3aff*c4aff*pmu*pnu*q2*xM**2*xMdelta**2 - 
     -      4*c3vff*c4vff*pmu*pnu*q2*xM**2*xMdelta**2 + 
     -      4*c3aff*c4aff*gmunu*q2**2*xM**2*xMdelta**2 - 
     -      4*c3vff*c4vff*gmunu*q2**2*xM**2*xMdelta**2 - 
     -      2*c4aff*c5aff*p2*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c4vff*c6vff*p2*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c3aff**2*pdotq*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c3vff**2*pdotq*pnu*qmu*xM**2*xMdelta**2 - 
     -      4*c3aff*c4aff*pdotq*pnu*qmu*xM**2*xMdelta**2 + 
     -      4*c3vff*c4vff*pdotq*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c4aff*c5aff*pdotq*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c4vff*c6vff*pdotq*pnu*qmu*xM**2*xMdelta**2 - 
     -      2*c4aff*c5aff*p2*pmu*qnu*xM**2*xMdelta**2 - 
     -      2*c4vff*c6vff*p2*pmu*qnu*xM**2*xMdelta**2 - 
     -      2*c3aff**2*pdotq*pmu*qnu*xM**2*xMdelta**2 - 
     -      2*c3vff**2*pdotq*pmu*qnu*xM**2*xMdelta**2 - 
     -      4*c3aff*c4aff*pdotq*pmu*qnu*xM**2*xMdelta**2 + 
     -      4*c3vff*c4vff*pdotq*pmu*qnu*xM**2*xMdelta**2 - 
     -      2*c4aff*c5aff*pdotq*pmu*qnu*xM**2*xMdelta**2 - 
     -      2*c4vff*c6vff*pdotq*pmu*qnu*xM**2*xMdelta**2 + 
     -      2*c3aff**2*p2*qmu*qnu*xM**2*xMdelta**2 + 
     -      2*c3vff**2*p2*qmu*qnu*xM**2*xMdelta**2 - 
     -      4*c4aff*c5aff*p2*qmu*qnu*xM**2*xMdelta**2 - 
     -      4*c4vff*c6vff*p2*qmu*qnu*xM**2*xMdelta**2 - 
     -      8*c3aff*c4aff*pdotq*qmu*qnu*xM**2*xMdelta**2 + 
     -      8*c3vff*c4vff*pdotq*qmu*qnu*xM**2*xMdelta**2 - 
     -      4*c4aff*c5aff*pdotq*qmu*qnu*xM**2*xMdelta**2 - 
     -      4*c4vff*c6vff*pdotq*qmu*qnu*xM**2*xMdelta**2 - 
     -      4*c3aff*c4aff*q2*qmu*qnu*xM**2*xMdelta**2 + 
     -      4*c3vff*c4vff*q2*qmu*qnu*xM**2*xMdelta**2 - 
     -      (0,4)*c3vff*c5aff*epsmunupq*xM**4*xMdelta**2 + 
     -      (0,4)*c3aff*c6vff*epsmunupq*xM**4*xMdelta**2 - 
     -      (0,2)*c5aff*c6vff*epsmunupq*xM**4*xMdelta**2 + 
     -      2*c5aff**2*gmunu*p2*xM**4*xMdelta**2 + 
     -      2*c6vff**2*gmunu*p2*xM**4*xMdelta**2 + 
     -      4*c3aff*c5aff*gmunu*pdotq*xM**4*xMdelta**2 + 
     -      2*c5aff**2*gmunu*pdotq*xM**4*xMdelta**2 - 
     -      4*c3vff*c6vff*gmunu*pdotq*xM**4*xMdelta**2 + 
     -      2*c6vff**2*gmunu*pdotq*xM**4*xMdelta**2 + 
     -      4*c3aff*c5aff*gmunu*q2*xM**4*xMdelta**2 - 
     -      4*c3vff*c6vff*gmunu*q2*xM**4*xMdelta**2 - 
     -      2*c3aff*c5aff*pnu*qmu*xM**4*xMdelta**2 + 
     -      2*c3vff*c6vff*pnu*qmu*xM**4*xMdelta**2 - 
     -      2*c3aff*c5aff*pmu*qnu*xM**4*xMdelta**2 + 
     -      2*c3vff*c6vff*pmu*qnu*xM**4*xMdelta**2 - 
     -      4*c3aff*c5aff*qmu*qnu*xM**4*xMdelta**2 + 
     -      4*c3vff*c6vff*qmu*qnu*xM**4*xMdelta**2 - 
     -      (0,4)*c3vff*c4aff*epsmunupq*pdotq*xM*xMdelta**3 - 
     -      (0,4)*c3aff*c4vff*epsmunupq*pdotq*xM*xMdelta**3 + 
     -      2*c3aff*c4aff*gmunu*pdotq**2*xM*xMdelta**3 + 
     -      2*c4aff**2*gmunu*pdotq**2*xM*xMdelta**3 + 
     -      2*c3vff*c4vff*gmunu*pdotq**2*xM*xMdelta**3 - 
     -      2*c4vff**2*gmunu*pdotq**2*xM*xMdelta**3 - 
     -      (0,4)*c3vff*c4aff*epsmunupq*q2*xM*xMdelta**3 - 
     -      (0,4)*c3aff*c4vff*epsmunupq*q2*xM*xMdelta**3 + 
     -      2*c3aff*c4aff*gmunu*pdotq*q2*xM*xMdelta**3 + 
     -      4*c4aff**2*gmunu*pdotq*q2*xM*xMdelta**3 + 
     -      2*c3vff*c4vff*gmunu*pdotq*q2*xM*xMdelta**3 - 
     -      4*c4vff**2*gmunu*pdotq*q2*xM*xMdelta**3 + 
     -      2*c3aff*c4aff*pmu*pnu*q2*xM*xMdelta**3 + 
     -      2*c4aff**2*pmu*pnu*q2*xM*xMdelta**3 + 
     -      2*c3vff*c4vff*pmu*pnu*q2*xM*xMdelta**3 - 
     -      2*c4vff**2*pmu*pnu*q2*xM*xMdelta**3 + 
     -      2*c4aff**2*gmunu*q2**2*xM*xMdelta**3 - 
     -      2*c4vff**2*gmunu*q2**2*xM*xMdelta**3 - 
     -      2*c3aff*c4aff*pdotq*pnu*qmu*xM*xMdelta**3 - 
     -      2*c4aff**2*pdotq*pnu*qmu*xM*xMdelta**3 - 
     -      2*c3vff*c4vff*pdotq*pnu*qmu*xM*xMdelta**3 + 
     -      2*c4vff**2*pdotq*pnu*qmu*xM*xMdelta**3 - 
     -      2*c3aff*c4aff*pdotq*pmu*qnu*xM*xMdelta**3 - 
     -      2*c4aff**2*pdotq*pmu*qnu*xM*xMdelta**3 - 
     -      2*c3vff*c4vff*pdotq*pmu*qnu*xM*xMdelta**3 + 
     -      2*c4vff**2*pdotq*pmu*qnu*xM*xMdelta**3 - 
     -      2*c3aff*c4aff*pdotq*qmu*qnu*xM*xMdelta**3 - 
     -      4*c4aff**2*pdotq*qmu*qnu*xM*xMdelta**3 - 
     -      2*c3vff*c4vff*pdotq*qmu*qnu*xM*xMdelta**3 + 
     -      4*c4vff**2*pdotq*qmu*qnu*xM*xMdelta**3 - 
     -      2*c4aff**2*q2*qmu*qnu*xM*xMdelta**3 + 
     -      2*c4vff**2*q2*qmu*qnu*xM*xMdelta**3 - 
     -      (0,4)*c3vff*c5aff*epsmunupq*xM**3*xMdelta**3 - 
     -      (0,4)*c3aff*c6vff*epsmunupq*xM**3*xMdelta**3 + 
     -      2*c3aff*c5aff*gmunu*pdotq*xM**3*xMdelta**3 + 
     -      4*c4aff*c5aff*gmunu*pdotq*xM**3*xMdelta**3 + 
     -      2*c3vff*c6vff*gmunu*pdotq*xM**3*xMdelta**3 - 
     -      4*c4vff*c6vff*gmunu*pdotq*xM**3*xMdelta**3 + 
     -      2*c3aff**2*gmunu*q2*xM**3*xMdelta**3 - 
     -      2*c3vff**2*gmunu*q2*xM**3*xMdelta**3 + 
     -      4*c4aff*c5aff*gmunu*q2*xM**3*xMdelta**3 - 
     -      4*c4vff*c6vff*gmunu*q2*xM**3*xMdelta**3 - 
     -      c3aff*c5aff*pnu*qmu*xM**3*xMdelta**3 - 
     -      2*c4aff*c5aff*pnu*qmu*xM**3*xMdelta**3 - 
     -      c3vff*c6vff*pnu*qmu*xM**3*xMdelta**3 + 
     -      2*c4vff*c6vff*pnu*qmu*xM**3*xMdelta**3 - 
     -      c3aff*c5aff*pmu*qnu*xM**3*xMdelta**3 - 
     -      2*c4aff*c5aff*pmu*qnu*xM**3*xMdelta**3 - 
     -      c3vff*c6vff*pmu*qnu*xM**3*xMdelta**3 + 
     -      2*c4vff*c6vff*pmu*qnu*xM**3*xMdelta**3 - 
     -      2*c3aff**2*qmu*qnu*xM**3*xMdelta**3 + 
     -      2*c3vff**2*qmu*qnu*xM**3*xMdelta**3 - 
     -      4*c4aff*c5aff*qmu*qnu*xM**3*xMdelta**3 + 
     -      4*c4vff*c6vff*qmu*qnu*xM**3*xMdelta**3 + 
     -      2*c5aff**2*gmunu*xM**5*xMdelta**3 - 
     -     2*c6vff**2*gmunu*xM**5*xMdelta**3 - 
     -      2*c6aff**2*qmu*qnu*(p2 + pdotq + xM*xMdelta)*
     -       (pdotq**2 + 2*pdotq*q2 + q2*(q2 - xMdelta**2)) - 
     -      2*c5vff**2*(p2 + pdotq - xM*xMdelta)*
     -       (pmu*(pnu*q2 - pdotq*qnu)*(q2 - xMdelta**2) + 
     -         pdotq*(pnu*qmu*(-q2 + xMdelta**2) + 
     -            pdotq*(qmu*qnu - gmunu*xMdelta**2))) - 
     -      c6aff*((pnu*q2*qmu + pmu*q2*qnu - 2*pdotq*qmu*qnu)*xMdelta*
     -          (2*c4aff*xMdelta*(p2 + pdotq + xM*xMdelta) + 
     -            c3aff*xM*(p2 - q2 + 2*xM*xMdelta + xMdelta**2)) + 
     -         2*c5aff*xM**2*(p2 + pdotq + xM*xMdelta)*
     -          (pnu*q2*qmu + pdotq*(pnu*qmu + (pmu + 2*qmu)*qnu) + 
     -            qnu*(pmu*q2 + 2*q2*qmu - 2*qmu*xMdelta**2))) + 
     -      2*c5vff*(c6vff*xM**2*(p2 + pdotq - xM*xMdelta)*
     -          (pnu*qmu*(q2 - xMdelta**2) - 
     -            pdotq*(pnu*qmu + 2*qmu*qnu - 2*gmunu*xMdelta**2) + 
     -            pmu*(2*pnu*q2 - qnu*(pdotq - q2 + xMdelta**2))) + 
     -         xMdelta*((0,1)*c3aff*epsmunupq*pdotq*xM*
     -             (pdotq + q2 + 2*(xM - xMdelta)*xMdelta) - 
     -            xMdelta*((0,1)*c4aff*epsmunupq*pdotq*(pdotq + q2) + 
     -               (0,1)*c5aff*epsmunupq*pdotq*xM**2 - 
     -               2*c4vff*(gmunu*pdotq*(pdotq + q2) - 
     -     pdotq*qmu*(pnu + qnu) + 
     -     pmu*(pnu*q2 - pdotq*qnu))*(p2 + pdotq - xM*xMdelta)) + 
     -  c3vff*xM*(pdotq*pmu*pnu*q2 - pmu*pnu*q2**2 - pdotq**2*pnu*qmu + 
     -    pdotq*pnu*q2*qmu - pdotq**2*pmu*qnu + pdotq*pmu*q2*qnu - 
     -               2*pdotq**2*qmu*qnu - 
     - p2*(pdotq*qmu*(pnu + qnu) + pmu*(-(pnu*q2) + pdotq*qnu)) - 
     -     2*pmu*pnu*q2*xM*xMdelta + 2*pdotq*pnu*qmu*xM*xMdelta + 
     -  2*pdotq*pmu*qnu*xM*xMdelta + 2*pdotq*qmu*qnu*xM*xMdelta + 
     -   pmu*pnu*q2*xMdelta**2 - pdotq*pnu*qmu*xMdelta**2 - 
     -      pdotq*pmu*qnu*xMdelta**2 + 
     - gmunu*pdotq*(pdotq**2 + p2*(pdotq + q2) - 2*q2*xM*xMdelta + 
     -                  pdotq*(q2 - 2*xM*xMdelta + xMdelta**2)))))))/
     -  (3.*xM**4*xMdelta**2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
         subroutine gammunus
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension cres(4,4)
      common/gammunu3common/cresmul(4,4,16)
      dimension cAUX1(4,4)
c
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
c (1/3*gamma^\mu*gamma^\nu)
       do mu=0,3
       do nu=0,3
        index=mu+4*nu+1  
        call MatrixProduct(cga(1,1,mu),cga(1,1,nu),cAUX1) 
        xconstant=1.d0/3.d0
        call ScalarProductMatrix(xconstant,cAUX1,cres)
       do i=1,4
       do j=1,4
        cresmul(j,i,index)=cres(j,i)
       enddo
       enddo
c
       enddo
       enddo
       return
       end


      subroutine gammunu3(mu,nu,cres)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension cres(4,4)
      common/gammunu3common/cresmul(4,4,16)
      dimension cAUX1(4,4)
c

       index=mu+4*nu+1  
       do i=1,4
       do j=1,4
        cres(j,i)=cresmul(j,i,index)
       enddo
       enddo     
   
      return 
      end       
      
c For the Levi-Civita symbol in four dimensions	
	integer function iepsi(ia,ib,ic,id)
	iepsi=0
	if(ia.eq.0)then
	 if(ib.eq.1)then
	   if(ic.eq.2.and.id.eq.3)iepsi=1
	   if(ic.eq.3.and.id.eq.2)iepsi=-1
	 elseif(ib.eq.2)then	 
	   if(ic.eq.3.and.id.eq.1)iepsi=1
	   if(ic.eq.1.and.id.eq.3)iepsi=-1
	 elseif(ib.eq.3)then
	   if(ic.eq.1.and.id.eq.2)iepsi=1
	   if(ic.eq.2.and.id.eq.1)iepsi=-1	 
	 endif
	else if(ia.eq.1)then
	 if(ib.eq.0)then
	   if(ic.eq.3.and.id.eq.2)iepsi=1
	   if(ic.eq.2.and.id.eq.3)iepsi=-1
	 elseif(ib.eq.2)then	 
	   if(ic.eq.0.and.id.eq.3)iepsi=1
	   if(ic.eq.3.and.id.eq.0)iepsi=-1
	 elseif(ib.eq.3)then
	   if(ic.eq.2.and.id.eq.0)iepsi=1
	   if(ic.eq.0.and.id.eq.2)iepsi=-1	 
	 endif
	else if(ia.eq.2)then
	 if(ib.eq.0)then
	   if(ic.eq.1.and.id.eq.3)iepsi=1
	   if(ic.eq.3.and.id.eq.1)iepsi=-1
	 elseif(ib.eq.1)then	 
	   if(ic.eq.3.and.id.eq.0)iepsi=1
	   if(ic.eq.0.and.id.eq.3)iepsi=-1
	 elseif(ib.eq.3)then
	   if(ic.eq.0.and.id.eq.1)iepsi=1
	   if(ic.eq.1.and.id.eq.0)iepsi=-1	 
	 endif
	else if(ia.eq.3)then
	 if(ib.eq.0)then
	   if(ic.eq.2.and.id.eq.1)iepsi=1
	   if(ic.eq.1.and.id.eq.2)iepsi=-1
	 elseif(ib.eq.1)then	 
	   if(ic.eq.0.and.id.eq.2)iepsi=1
	   if(ic.eq.2.and.id.eq.0)iepsi=-1
	 elseif(ib.eq.2)then
	   if(ic.eq.1.and.id.eq.0)iepsi=1
	   if(ic.eq.0.and.id.eq.1)iepsi=-1	 
	 endif
	endif
	return
	end
***********************************************************************	
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccc    CURRENTS      ccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  define the first current, the contribution 
c coming from the direct Delta(1232) pole term for the 2body diagrams. 
c Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCdelta, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c rho is the density in fm^(-3)
c
c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors
      subroutine CurrentDeltaPole2body(mu,xkpion,p,q,rho,cjota)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),pdelta(0:3)
      dimension cjota(4,4),g(0:3,0:3),cga5(4,4),cIdent(4,4)
      dimension cga(4,4,0:3),cProjector(4,4),cGamma(4,4)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cAUX(4,4),cAUX1(4,4)
      dimension cAUX2(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
c new
      dimension cdproj(0:3,0:3,4,4)
      common/deltapole2/cdproj
!$OMP THREADPRIVATE(/deltapole2/)
c new
      parameter(fstar=2.14d0,xmpion=0.13957d0,xCabibbo=0.974d0,
     &     xMdelta=1.232d0)
      xCdelta=1.d0
c      xmpion=0.13957d0
      gprime=0.63d0

c First we  define the constant that goes in front of
c the matrix
      call SumFourVectors(p,q,pdelta)
      call LorentzScalarProduct(pdelta,pdelta,pdelta2)
      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
      xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)
      ff=Fpi(xkpion)

      s=pdelta2
      if(s.lt.0.d0)then
         s=0.d0
      endif
      W=dsqrt(s)
      constant=ci*xCdelta*fstar/xmpion*dsqrt(3.d0)*xCabibbo*
     &     (1.d0+gprime*(xkpion2-xmpion**2)/(ff**2*xmodkpi**2))/((W+
     &     xMdelta)*(W-xMdelta-ReSelfenergyDelta(rho)+ci*(
     &     PauliDeltaWidth(s,rho)/2.d0-xImSelfenergyDelta(s,rho))))

c have to construct (pdeltaslash+Mdelta) to give this matrix to
c the subroutine DeltaProjector
c      call slash(pdelta,cpdeltaslash)
c      call MatrixSum(cpdeltaslash,cMdeltaMatrix,cpdMd)

c initialize a matrix to zero
      do i=1,4
         do j=1,4
            cAUX(i,j)=0.d0
         enddo
      enddo

c construct the matrix before multiplying by constant

       call slash(q,cqslash)
       do ibeta=0,3
       do i=1,4
         do j=1,4
            cAUX1(i,j)=0.d0
         enddo
       enddo
           call Vertex(ibeta,mu,p,q,cqslash,cGamma)
         do ialpha=0,3
            ialphaprima=ialpha
            ibetaprima=ibeta
c new
                  do ii=1,4
                  do jj=1,4
                  cprojector(ii,jj)=cdproj(ialpha,ibeta,ii,jj)
                   enddo
                  enddo
c new
c            call MatrixProduct(cProjector,cGamma,cAUX1)
            xcons=xkpion(ialpha)*g(ialpha,ialphaprima)*
     &           g(ibetaprima,ibeta)
            call ScalarProductMatrix(xcons,cProjector,cAUX2)
            call MatrixSum(cAUX1,cAUX2,cAUX1)
         enddo
            call MatrixProduct(caux1,cGamma,cAUX2)
            call MatrixSum(cAUX2,cAUX,cAUX)
      enddo
c have the matricial sum stored in cAUX. And we have to multiply
c this matrix by the complex constant stored in constant
      call ComplexScalarProductMatrix(constant,cAUX,cjota)
	
      return	
      end

c******************************************************************
c  define the second current, the contribution 
c coming from the crossed Delta(1232) pole term. Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCcrossdelta, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c
c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors
c******************************************************************
      subroutine CurrentCrossedDeltaPole2body(mu,xkpion,p,q,cjota)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),pdelta(0:3)
      dimension cjota(4,4),g(0:3,0:3),cga5(4,4),cIdent(4,4)
      dimension cga(4,4,0:3),cProjector(4,4),cGamma(4,4)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cAUX(4,4),cAUX1(4,4)
      dimension cAUX2(4,4),qprima(0:3),pprima(0:3)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/caux1sc/caux1s(4,4,0:3)
!$OMP THREADPRIVATE(/caux1sc/)

      parameter(fstar=2.14d0,xmpion=0.13957d0,xCabibbo=0.974d0,
     &     xMdelta=1.232d0)
      xCdelta=1.d0
c      xmpion=0.13957d0
      gprime=0.63d0
      ff=Fpi(xkpion)
c First we  define the constant that goes in front of
c the matrix
      call RestaFourVectors(p,xkpion,pdelta)
      call LorentzScalarProduct(pdelta,pdelta,pdelta2)
      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
      xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)

      constant=ci*xCdelta*fstar/xmpion*1.d0/dsqrt(3.d0)*xCabibbo*
     &     (1.d0+gprime*(xkpion2-xmpion**2)/(ff**2*xmodkpi**2))/
     &     (pdelta2-xMdelta**2+ci*xMdelta*DeltaWidth(pdelta2))

c have to construct (pdeltaslash+Mdelta) to give this matrix to
c the subroutine DeltaProjector
      call slash(pdelta,cpdeltaslash)
      call MatrixSum(cpdeltaslash,cMdeltaMatrix,cpdMd)

c initialize a matrix to zero
      do i=1,4
         do j=1,4
            cAUX(i,j)=0.d0
         enddo
      enddo
c need to define the four-vectors pprima and qprima=-q
      call SumFourVectors(q,pdelta,pprima)
      do i=0,3
         qprima(i)=-q(i)
      enddo
      

c construct the matrix before multiplying by constant
      do ialpha=0,3
            call VertexHermiticConjugate(ialpha,mu,pprima,
     &           qprima,cGamma)
         do i=1,4
           do j=1,4
            cAUX1(i,j)=0.d0
          enddo
        enddo
 
         do ibeta=0,3
            ialphaprima=ialpha
            ibetaprima=ibeta
           call DeltaProjector(ialphaprima,ibetaprima,pdelta,
     &           cpdMd,cProjector)
c            call MatrixProduct(cGamma,cProjector,cAUX1)
            xcons=xkpion(ibeta)*g(ialpha,ialphaprima)*
     &           g(ibetaprima,ibeta)
            call ScalarProductMatrix(xcons,cprojector,cAUX2)
            call MatrixSum(cAUX1,cAUX2,cAUX1)
         enddo

         do i=1,4
           do j=1,4
            cAUX1s(i,j,ialpha)=cAUX1(i,j)
          enddo
        enddo
 
           call MatrixProduct(cGamma,caux1,cAUX2)
           call MatrixSum(cAUX2,cAUX,cAUX)
      enddo
c have the matricial sum stored in cAUX. And we have to multiply
c this matrix by the complex constant stored in constant
      call ComplexScalarProductMatrix(constant,cAUX,cjota)
	
      return	
      end
c
c
      subroutine CurrentCrossedDeltaPole2body2(mu,xkpion,p,q,cjota)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),pdelta(0:3)
      dimension cjota(4,4),g(0:3,0:3),cga5(4,4),cIdent(4,4)
      dimension cga(4,4,0:3),cProjector(4,4),cGamma(4,4)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cAUX(4,4),cAUX1(4,4)
      dimension cAUX2(4,4),qprima(0:3),pprima(0:3)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common/caux1sc/caux1s(4,4,0:3)
!$OMP THREADPRIVATE(/caux1sc/)

      parameter(fstar=2.14d0,xmpion=0.13957d0,xCabibbo=0.974d0,
     &     xMdelta=1.232d0)
      xCdelta=1.d0
c      xmpion=0.13957d0
      gprime=0.63d0
      ff=Fpi(xkpion)
c First we  define the constant that goes in front of
c the matrix
      call RestaFourVectors(p,xkpion,pdelta)
      call LorentzScalarProduct(pdelta,pdelta,pdelta2)
      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
      xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)

      constant=ci*xCdelta*fstar/xmpion*1.d0/dsqrt(3.d0)*xCabibbo*
     &     (1.d0+gprime*(xkpion2-xmpion**2)/(ff**2*xmodkpi**2))/
     &     (pdelta2-xMdelta**2+ci*xMdelta*DeltaWidth(pdelta2))

c have to construct (pdeltaslash+Mdelta) to give this matrix to
c the subroutine DeltaProjector
      call slash(pdelta,cpdeltaslash)
      call MatrixSum(cpdeltaslash,cMdeltaMatrix,cpdMd)

c initialize a matrix to zero
      do i=1,4
         do j=1,4
            cAUX(i,j)=0.d0
         enddo
      enddo
c need to define the four-vectors pprima and qprima=-q
      call SumFourVectors(q,pdelta,pprima)
      do i=0,3
         qprima(i)=-q(i)
      enddo
      

c construct the matrix before multiplying by constant
      do ialpha=0,3
            call VertexHermiticConjugate(ialpha,mu,pprima,
     &           qprima,cGamma)
 
         do i=1,4
           do j=1,4
            cAUX1(i,j)=cAUX1s(i,j,ialpha)
          enddo
        enddo
 
           call MatrixProduct(cGamma,caux1,cAUX2)
           call MatrixSum(cAUX2,cAUX,cAUX)
      enddo
c have the matricial sum stored in cAUX. And we have to multiply
c this matrix by the complex constant stored in constant
      call ComplexScalarProductMatrix(constant,cAUX,cjota)
	
      return	
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc FORM FACTORS FOR THE WEAK TRANSITION NUCLEON-NUCLEON cccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c First the electric form factor G_Eproton(q^2)
	real*8 function G_Eproton(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xMD=0.843d0)

c	call LorentzScalarProduct(q,q,q2)

	G_Eproton=(1.d0/(1.d0-q2/xMD**2))**2

	return
	end

cmagnetic form factor for the proton G_Mproton(q^2)
	real*8 function G_Mproton(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xmagnetonproton=2.792847d0)

	G_Mproton=xmagnetonproton*G_Eproton(q2)

	return
	end

cmagnetic form factor for the neutron G_Mneutron(q^2)
	real*8 function G_Mneutron(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xmagnetonneutron=-1.913043d0)

	G_Mneutron=xmagnetonneutron*G_Eproton(q2)

	return
	end

celectric form factor for the neutron G_Eneutron(q^2)
	real*8 function G_Eneutron(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xmagnetonneutron=-1.913043d0,xlambdan=5.6d0,
     &	xM=0.94d0)

c	call LorentzScalarProduct(q,q,q2)
	
	tau=-q2/(4.d0*xM**2)

	G_Eneutron=-xmagnetonneutron*tau/(1.d0+xlambdan*tau)*
     &	G_Eproton(q2)

	return
	end

cform factor F1p(q^2)
	real*8 function F1p(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xM=0.94d0)

c	call LorentzScalarProduct(q,q,q2)

	tau=-q2/(4.d0*xM**2)

	F1p=(G_Eproton(q2)+tau*G_Mproton(q2))/(1.d0+tau)

	return
	end

cform factor F1n(q^2)
	real*8 function F1n(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xM=0.94d0)

c	call LorentzScalarProduct(q,q,q2)

	tau=-q2/(4.d0*xM**2)

	F1n=(G_Eneutron(q2)+tau*G_Mneutron(q2))/(1.d0+tau)

	return
	end

cform factor mupF2p(q^2)
	real*8 function xmupF2p(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xM=0.94d0)

c	call LorentzScalarProduct(q,q,q2)

	tau=-q2/(4.d0*xM**2)

	xmupF2p=(G_Mproton(q2)-G_Eproton(q2))/(1.d0+tau)

	return
	end

cform factor munF2n(q^2)
	real*8 function xmunF2n(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(xM=0.94d0)

c	call LorentzScalarProduct(q,q,q2)

	tau=-q2/(4.d0*xM**2)

	xmunF2n=(G_Mneutron(q2)-G_Eneutron(q2))/(1.d0+tau)

	return
	end

cvector form factor F1V(q^2)
	real*8 function F1V(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)

	F1V=1.d0/2.d0*(F1p(q2)-F1n(q2))

	return
	end

cvector form factor muVF2V(q^2)
	real*8 function xmuVF2V(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)

	xmuVF2V=1.d0/2.d0*(xmupF2p(q2)-xmunF2n(q2))

	return
	end

caxial form factor G_A(q^2)
	real*8 function G_A(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)
	parameter(gA=1.26d0,xMA=1.05d0)

c	call LorentzScalarProduct(q,q,q2)
	
	G_A=gA/(1.d0-q2/xMA**2)**2

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc NUCLEON CURRENTS ccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c First the VECTOR NUCLEON CURRENT.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Input variables
c 'ialpha' is a Lorentz index for the current
c 'q' is the four-momentum of the W+

c Output variable
c cVector is the corresponding matrix

      subroutine VectorN(ialpha,q,cqslash,cVector)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3),cVector(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4)
c Auxiliary matrices
      dimension cAUX1(4,4),cqslash(4,4),cAUX2(4,4),cAUX3(4,4)
      dimension cAUX4(4,4),cAUX5(4,4)
      parameter(xM=0.94d0)
      common /matrices/ cga,cga5,cIdent
      common /complexi/ ci

      call LorentzScalarProduct(q,q,q2)
c Let's work on the first term 2*F1V(q2)*gamma^\alpha
      xconstant1=2.d0*F1V(q2)
      call ScalarProductMatrix(xconstant1,cga(1,1,ialpha),cAUX1)
c the first term is stored in cAUX1. We cannot overwrite this 
c matrix

c Let's work on the second term of the current
      constant2=2.d0*ci*xmuVF2V(q2)/(2.d0*xM)*ci/2.d0
c  call slash(q,cqslash)
      call MatrixProduct(cga(1,1,ialpha),cqslash,cAUX2)
      call MatrixProduct(cqslash,cga(1,1,ialpha),cAUX3)
      call MatrixResta(cAUX2,cAUX3,cAUX4)
      call ComplexScalarProductMatrix(constant2,cAUX4,cAUX5)
c have the second term stored in cAUX5, we cannot overwrite
c this matrix

c must sum both contributions
      call MatrixSum(cAUX1,cAUX5,cVector)

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Secondly, the AXIAL-VECTOR NUCLEON CURRENT.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Input variables
c 'ialpha' is a Lorentz index for the current
c 'q' is the four-momentum of the W+

c Output variable
c cAxial is the corresponding matrix

      subroutine AxialN(ialpha,q,cqslash,cAxial)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3),cAxial(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4)
c Auxiliary matrices
      dimension cqslash(4,4),cAUX1(4,4),cAUX2(4,4),cAUX3(4,4)
      parameter(xmpion=0.13957d0)
      common /matrices/ cga,cga5,cIdent

      call LorentzScalarProduct(q,q,q2)

c The first term is only gamma^\alpha
c So, let's work on the second term
c The second term is q^\alpha/(mpion^2-q^2)*qslash
c      call slash(q,cqslash)
      xconstant=q(ialpha)/(xmpion**2-q2)
      call ScalarProductMatrix(xconstant,cqslash,cAUX1)
c have the second term stored in cAUX1. We cannot overwrite
c this matrix

c have to sum the first term with the second term
      call MatrixSum(cga(1,1,ialpha),cAUX1,cAUX2)
c have the sum stored in cAUX2. This matrix cannot be 
c overwritten

c multiply by the axial form factor G_A(q^2)
      xcons=G_A(q2)
      call ScalarProductMatrix(xcons,cAUX2,cAUX3)
c must multiply cAUX3 by gamma5 in order to obtain
c cAxial
      call MatrixProductg5(cAUX3,cAxial)

      return
      end
******************************************************************
*
******************************************************************
c  define the third current, the contribution 
c coming from the Nucleon pole term for the 2body absorption.
c With the factor correction of gprime
c Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCNP, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four-momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c rho is the density in fm^(-3)

c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors
      subroutine CurrentNucleonPole2body(mu,xkpion,p,q,rho,cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),cAUX(4,4),cjota(4,4)
      dimension pplusq(0:3),ckpionslash(4,4),cpplusqslash(4,4)
      dimension cAUX1(4,4),cAUX2(4,4),cqslash(4,4),cpdeltaslash(4,4)
      dimension cMdeltaMatrix(4,4),cMidentity(4,4),cAUX3(4,4)
      dimension cVector(4,4),cAxial(4,4),cAUX4(4,4)
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
      parameter(gA=1.26d0,fpion=0.093d0,xCabibbo=0.974d0,xM=0.94d0)

cfactor that goes in front of the matrix
c First we must define p+q (pplusq)
      call SumFourVectors(p,q,pplusq)
      call LorentzScalarProduct(pplusq,pplusq,pplusq2)
      xmpion=0.13957d0
      rho0=0.17d0
      width=0.020d0             !!! width in GeV ~ 20 MeV
      hbarc=0.19733d0
      gprime=0.63d0
      ff=Fpi(xkpion)
c	if(pplusq2.lt.(xmpion+xM)**2)then
c	   pplusq2=(xmpion+xM)**2
c	endif
c      s=pplusq2
c      if (s .lt. 0.d0)then
c         s=0.d0
c      endif
c      W=dsqrt(s)
c      dk0=(pplusq(0)-xM)/hbarc
c      dk=dsqrt(pplusq(1)**2+pplusq(2)**2+pplusq(3)**2)/hbarc
c      xImSigma=DIMASM(dk0,dk)*hbarc
      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
      xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)

      constant=-ci*gA/(dsqrt(2.d0)*fpion)*xCabibbo/(pplusq2-xM**2+ci*
     &     xM*width)*(1.d0+gprime*(xkpion2-xmpion**2)/
     &     (ff**2*xmodkpi**2))

c need to define cpplusqslash=cpslash+cqslash
c and ckpionslash

      call slash(xkpion,ckpionslash)
      call slash(pplusq,cpplusqslash)

      call MatrixProductg5(ckpionslash,cAUX1)
c this first result is stored in cAUX1, we cannot overwrite this matrix
      call MatrixSum(cpplusqslash,cMidentity,cAUX2)
c this second result is stored in cAUX2. must multiply cAUX1 and
c cAUX2

      call MatrixProduct(cAUX1,cAUX2,cAUX3)
c This third result is stored in cAUX3. We cannot overwrite this matrix

c must call VectorN and AxialN
      call slash(q,cqslash)
      call VectorN(mu,q,cqslash,cVector)
      call AxialN(mu,q,cqslash,cAxial)

      call MatrixResta(cVector,cAxial,cAUX4)
c Finally, we must multiply cAUX3 by cAUX4
      call MatrixProduct(cAUX3,cAUX4,cAUX)

c And obtain the current
      call ComplexScalarProductMatrix(constant,cAUX,cjota)

      return
      end
******************************************************************

******************************************************************
c  define the fourth current,the contribution 
c coming from the Crossed Nucleon pole term for the 2Body Absorption.
c Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCCrossedNP, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four-momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c
c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors
	subroutine CurrentCrossedNucleonPole2body(mu,xkpion,p,q,cjota)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	dimension xkpion(0:3),p(0:3),q(0:3),cAUX(4,4),cjota(4,4)
	dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
	dimension cMidentity(4,4),pminuskpion(0:3),cVector(4,4)
	dimension cAxial(4,4),cAUX1(4,4),cpminuskpionslash(4,4)
	dimension cAUX2(4,4),cAUX3(4,4),ckpionslash(4,4),cAUX4(4,4)
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
	parameter(gA=1.26d0,fpion=0.093d0,xCabibbo=0.974d0,xM=0.94d0)
c Constants of interest
        gprime=0.63d0
        xmpion=0.13957d0 !!! 139.57 MeV
        ff=Fpi(xkpion)
c We need pprime-q, where pprime is the four-momentum of the
c outgoing nucleon.
c Applying global conservation of four-momenta, we have
c pprime-q == p-xkpion
	call RestaFourVectors(p,xkpion,pminuskpion)
	call LorentzScalarProduct(pminuskpion,pminuskpion,
     &	pminuskpion2)
c  complex constant that goes in front of the matrix
        call LorentzScalarProduct(xkpion,xkpion,xkpion2)
        xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)

	constant=-ci*gA/(dsqrt(2.d0)*fpion)*xCabibbo/(pminuskpion2-
     &	xM**2)*(1.d0+gprime*(xkpion2-xmpion**2)/
     &       (ff**2*xmodkpi**2))

c the matrices. First we must call to VectorN and AxialN
        call slash(q,cqslash)
	call VectorN(mu,q,cqslash,cVector)
	call AxialN(mu,q,cqslash,cAxial)

c  must substract these two matrices.
	call MatrixResta(cVector,cAxial,cAUX1)
c This first result is stored in cAUX1. We cannot overwrite this matrix

cnumerator of the propagator of the nucleon
	call slash(pminuskpion,cpminuskpionslash)
	call MatrixSum(cpminuskpionslash,cMidentity,cAUX2)
c this second result is stored in cAUX2. We cannot overwrite this matrix

c must multiply cAUX1 and cAUX2 to obtain cAUX3
	call MatrixProduct(cAUX1,cAUX2,cAUX3)
c We cannot overwrite cAUX3.

c kpionslash
	call slash(xkpion,ckpionslash)

c must multiply cAUX3 by ckpionslash
	call MatrixProduct(cAUX3,ckpionslash,cAUX4)
c We cannot overwrite cAUX4

c And finally we must multiply by gamma5
	call MatrixProductg5(cAUX4,cAUX)

c And finally, to obtain cjota, we must multiply the complex constant
c by cAUX
	call ComplexScalarProductMatrix(constant,cAUX,cjota)

	return
	end

************************************************************************
************************************************************************
************************************************************************
c before defining the contact term current, we must define form 
c factors.
c First the contribution proportional to gA in the CT (contact term)
c diagram, which is the one due to the vector weak transition.
	real*8 function FVCT(q2)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension q(0:3)

	FVCT=2.d0*F1V(q2)

	return
	end

c Secondly, the form factor Frho(t)
	real*8 function Frho(t)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
c	dimension pprime(0:3)
	parameter(xmrho=0.7758d0)

c	call LorentzScalarProduct(pprime,pprime,t)

	Frho=1.d0/(1.d0-t/xmrho**2)

	return
	end

******************************************************************
******************************************************************
******************************************************************
c  define the fifth current,the contribution 
c coming from the contact term. Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCContactTerm, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four-momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c
c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors
	subroutine CurrentContactTerm(mu,xkpion,p,q,cjota)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	dimension xkpion(0:3),p(0:3),q(0:3)
	dimension cAUX(4,4),cjota(4,4)
	dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4)
	dimension cAUX1(4,4),qminuskpion(0:3),cAUX2(4,4),cAUX3(4,4)
	common /matrices/ cga,cga5,cIdent
      	common /complexi/ ci
	parameter(gA=1.26d0,fpion=0.093d0,xCabibbo=0.974d0)

        call LorentzScalarProduct(q,q,q2)
c First the complex constant that goes in front of the matrix
	cconstant=-ci/(dsqrt(2.d0)*fpion)*xCabibbo

cterm that goes with gamma5
	xconstant1=gA*FVCT(q2)
c must multiply gamma^\mu by gamma5
	call MatrixProductg5(cga(1,1,mu),cAUX1)
c must multiply the real scalar xconstant1 by cAUX1
	call ScalarProductMatrix(xconstant1,cAUX1,cAUX2)
cfirst term of the matricial substraction is stored in cAUX2


c need to calculate the four-vector (q-xkpion)
	call RestaFourVectors(q,xkpion,qminuskpion)
        call LorentzScalarProduct(qminuskpion,qminuskpion,t)
cform factor Frho((q-kpion)^2)
	xconstant2=Frho(t)
c must multiply this real constant by gamma^\mu
	call ScalarProductMatrix(xconstant2,cga(1,1,mu),cAUX3)
csecond term is stored in cAUX3

c have to substract cAUX2-cAUX3 in order to obtain cAUX
	call MatrixResta(cAUX2,cAUX3,cAUX)

c In order to obtain cjota we must multiply cAUX by the complex 
c constant
	call ComplexScalarProductMatrix(cconstant,cAUX,cjota)

	return
	end

******************************************************************
******************************************************************
******************************************************************
c  define the sixth current,the contribution 
c coming from the pion pole term. Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCPionPole, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four-momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c
c Output variable:
c cjota is the matrix for the current
c cAUX is only the matrix, without the constant factors

      subroutine CurrentPionPole(mu,xkpion,p,q,cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),cAUX(4,4),cjota(4,4)
      dimension qminuskpion(0:3)
      common /complexi/ ci
      parameter(fpion=0.093d0,xCabibbo=0.974d0,xmpion=0.13957d0)

c First we need (q-kpion) and q2
      call LorentzScalarProduct(q,q,q2)
      call RestaFourVectors(q,xkpion,qminuskpion)
      call LorentzScalarProduct(qminuskpion,qminuskpion,t)
c can define the complex constant that goes in front of the
c matrix
	cconstant=-ci*Frho(t)*xCabibbo/(dsqrt(2.d0)*fpion)*
     &	q(mu)/(q2-xmpion**2)

cmatrix, which is simply qslash
	call slash(q,cAUX)

ccurrent cjota
	call ComplexScalarProductMatrix(cconstant,cAUX,cjota)

	return
	end
******************************************************************
******************************************************************
******************************************************************
c before programming the pion-in-flight current, we must define
c the form factor FPF(q), which is equal to FVCT(q)
      real*8 function FPF(q2)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
c      dimension q(0:3)

      FPF=FVCT(q2)

      return
      end
******************************************************************
******************************************************************
******************************************************************
c  define the seventh current,the contribution 
c coming from the pion-in-flight term. Note that without
c the spinors ubar and u multiplying from the left and from the right
c respectively, the current is a matrix.

cccccccccccccccccccccccc CAREFULLY cccccccccccccccccccccc
c It is lacking the isospin factor xCPionInFlight, which will be external
c to the subroutine and will be multiplying the cjota
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Input variables:
c mu is the Lorentz index of the current
c xkpion is the four-momentum of the outgoing pion
c p is the four-momentum of the incoming nucleon
c q is the four-momentum of the boson W+
c
c Output variable:
c cjota is the matrix for the current

      subroutine CurrentPionInFlight(mu,xkpion,p,q,cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),p(0:3),q(0:3),cjota(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4)
      dimension xkpionminusq(0:3)
      common /matrices/ cga,cga5,cIdent
      common /complexi/ ci
      parameter(gA=1.26d0,fpion=0.093d0,xCabibbo=0.974d0,
     &     xmpion=0.13957d0,xM=0.94d0)

      call LorentzScalarProduct(q,q,q2)
c We need the four-vector (xkpion-q)
      call RestaFourVectors(xkpion,q,xkpionminusq)
c need its squared
      call LorentzScalarProduct(xkpionminusq,xkpionminusq,xkpionminusq2)

ccomplex constant that goes in front of the matrix
      cconstant=-ci*FPF(q2)*gA/(dsqrt(2.d0)*fpion)*xCabibbo*
     &     (2.d0*xkpion(mu)-q(mu))/(xkpionminusq2-xmpion**2)*2.d0*xM

c must multiply this complex scalar by matrix gamma5
      call ComplexScalarProductMatrix(cconstant,cga5,cjota)

      return
      end

*******************************************************************
*******************************************************************
*******************************************************************
c define the subroutine SumCurrents, to optimize
c the sum of the seven currents with the proper isospin factors
c INPUT VARIABLES
c 'xc1','xc2',...,'xc7' will be the isospin factors
c 'cj1','cj2',...,'cj7' will be the currents (matrices)

c OUTPUT VARIABLES
c 'cJota' will be the total current (matrix)

      subroutine SumCurrents(xc1,xc2,xc3,xc4,xc5,xc6,xc7,cj1,cj2,cj3,
     &     cj4,cj5,cj6,cj7,cJota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cj1(4,4),cj2(4,4),cj3(4,4),cj4(4,4),cj5(4,4),cj6(4,4)
      dimension cj7(4,4),cJota(4,4)

      do i=1,4
         do j=1,4
            cJota(i,j)=xc1*cj1(i,j)+xc2*cj2(i,j)+xc3*cj3(i,j)+
     &           xc4*cj4(i,j)+xc5*cj5(i,j)+xc6*cj6(i,j)+xc7*cj7(i,j)
         enddo
      enddo
      
      return
      end

*******************************************************************

*******************************************************************
*******************************************************************
*******************************************************************
c  define the subroutine TensorABackground2body, which
c is the trace of the fermionic loop with the currents in every vertex


c For getting the pi+ and pi0 cross section, one must take into account the 
c currents J^\mu(p --> p + pion+) (and its hermitic conjugate),
c J^\mu(n --> n + pion+) (and its hermitic conjugate)
c and J^\mu(n --> p + pion0) (and its hermitic conjugate)
c So, there are three tensors, the one coming from p + W+ --> p + pi+
c the other one coming from n + W+ --> n + pi+ and
c the other one coming from n + W+ --> p + pi0

c INPUT VARIABLES
c 'mu' is the Lorentz index for the current in one vertex
c 'nu' is the other Lorentz index for the current in the hermitic
c conjugate vertex.

c 'xkpion' is the momentum of the outgoing pi+
c 'p' is the momentum of the incoming nucleon (proton or neutron)
c 'q' is the momentum of the incoming W+

c 'rho' is the density in fm^(-3) for the CurrentDeltaPole
c and for the CurrentNucleonPole2body


c OUTPUT VARIABLES
c 'cAproton' is the complex number which is the tensorial component
c for the given two Lorentz indices and the momenta configuration 
c for the reaction in proton (p + W+ --> p + pi+)

c 'cAneutronpiplus' is the complex number which is the tensorial component
c for the given two Lorentz indices and the momenta configuration
c for the reaction in neutron going to neutron + pi+

c 'cAneutronpi0' is the complex number which is the tensorial component
c for the given two Lorentz indices and the momenta configuration for
c the reaction in neutron going to proton + pi0

      subroutine TensorABackground2body(mu,nu,xkpion,pnum,pden,
     &     q,rho,cAproton,cAneutronpiplus,cAneutronpi0)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3),pnum(0:3),q(0:3),pplusq(0:3),pprime(0:3)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpslash(4,4),cpprimeslash(4,4)
      dimension cpM(4,4),cpprimeM(4,4)
      dimension cjota1(4,4),cjota2(4,4),cjota3(4,4),cjota4(4,4)
      dimension cjota5(4,4),cjota6(4,4),cjota7(4,4),cAUX1p(4,4)
      dimension cjota1neutron(4,4),cjota2neutron(4,4),cAUX1n(4,4)
      dimension cjota1bis(4,4),cjota2bis(4,4),cjota3bis(4,4)
      dimension cjota4bis(4,4),cjota5bis(4,4),cjota6bis(4,4)
      dimension cjota7bis(4,4),cAUX1pbis(4,4),cJotanudaggerp(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),cAUX2pbis(4,4)
      dimension cJtildenup(4,4),cjota1bisneutron(4,4)
      dimension cjota2bisneutron(4,4)
      dimension cAUX1nbis(4,4),cJotanudaggern(4,4),cAUX2nbis(4,4)
      dimension cJtildenun(4,4),cAUX3(4,4),cAUX4(4,4),cAUX5(4,4)
      dimension cjota1bisdaggerp(4,4),cjota1tildenup(4,4)
      dimension cjota1bisdaggern(4,4),cjota1tildenun(4,4)
      dimension cAUXpizero(4,4),cAUXpi0(4,4),cjota1pi0(4,4)
      dimension cJtildenupi0(4,4),cjota1tildepi0(4,4)
      dimension pden(0:3),cjota1p(4,4),cjota1bisp(4,4),cJotaproton(4,4)
      dimension cJotabisproton(4,4),cAUX2p(4,4),cJotatildeproton(4,4)
      dimension cjota1tildeproton(4,4),cJotaneutron(4,4)
      dimension cJotabisneutron(4,4),cAUX2n(4,4),cJotatildeneutron(4,4)
      dimension cjota1tildeneutron(4,4),cjota3tildeneutron(4,4)
      dimension cJotapi0(4,4),cJotatildepi0(4,4),cjota3proton(4,4)
      dimension cjota3pi0(4,4),cjota3tildeproton(4,4)
      dimension cjota3tildepi0(4,4)
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common /matrices/ cga,cga5,cIdent

c First we need pprime = p + q - kpion
      call SumFourVectors(pnum,q,pplusq)
      call RestaFourVectors(pplusq,xkpion,pprime)
c     need pslash and pprimeslash
      call slash(pnum,cpslash)
      call slash(pprime,cpprimeslash)
c need (pslash+M) and (pprimeslash+M)
      call MatrixSum(cpslash,cMidentity,cpM)
      call MatrixSum(cpprimeslash,cMidentity,cpprimeM)
c have stored (pslash+M) in cpM
c and (pprimeslash+M) in cpprimeM. We cannot overwrite these two 
c matrices.


c  put the isospin factors, for having them in mind
c First the isospin factors for p --> p + pi+
      xCdeltaproton=1.d0
      xCCrossedDeltaproton=1.d0
      xCNucleonPoleproton=0.d0
      xCCrossedNucleonPoleproton=1.d0
      xCcontactTermproton=1.d0
      xCPionPoleproton=1.d0
      xCPionInFlightproton=1.d0

c call all the currents in the vertex 'mu'
      call CurrentDeltaPole2body(mu,xkpion,pden,q,rho,cjota1)
      call CurrentCrossedDeltaPole2body(mu,xkpion,pden,q,cjota2)
      call CurrentNucleonPole2body(mu,xkpion,pden,q,rho,cjota3)
      call CurrentCrossedNucleonPole2body(mu,xkpion,pden,q,cjota4)
      call CurrentContactTerm(mu,xkpion,pden,q,cjota5)
      call CurrentPionPole(mu,xkpion,pden,q,cjota6)
      call CurrentPionInFlight(mu,xkpion,pden,q,cjota7)

c And we call all the currents in the vertex 'nu'
              if(nu.ne.mu)then
      call CurrentDeltaPole2body(nu,xkpion,pden,q,rho,cjota1bis)
      call CurrentCrossedDeltaPole2body2(nu,xkpion,pden,q,cjota2bis)
      call CurrentNucleonPole2body(nu,xkpion,pden,q,rho,cjota3bis)
      call CurrentCrossedNucleonPole2body(nu,xkpion,pden,q,cjota4bis)
      call CurrentContactTerm(nu,xkpion,pden,q,cjota5bis)
      call CurrentPionPole(nu,xkpion,pden,q,cjota6bis)
      call CurrentPionInFlight(nu,xkpion,pden,q,cjota7bis)
             else
              do iii=1,4
              do jjj=1,4
              cjota1bis(iii,jjj)=cjota1(iii,jjj)
              cjota2bis(iii,jjj)=cjota2(iii,jjj)
              cjota3bis(iii,jjj)=cjota3(iii,jjj)
              cjota4bis(iii,jjj)=cjota4(iii,jjj)
              cjota5bis(iii,jjj)=cjota5(iii,jjj)
              cjota6bis(iii,jjj)=cjota6(iii,jjj)
              cjota7bis(iii,jjj)=cjota7(iii,jjj)
              enddo
              enddo
              endif     


c IMPORTANT CURRENTS FOR THE SUBSTRACTION
c cjota1 and cjota1bis for the substraction of DeltaPoleXDeltaPole
c cjota3 and cjota3bis for the substraction of NPtimesNP


c Firstly, for protons, the current J^\mu for the process p --> p + pi+
      call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &     xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &     xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &     cjota1,cjota2,cjota3,cjota4,cjota5,cjota6,cjota7,
     &     cJotaproton)
c cJotaproton is the sum of all jotas for the process p --> p + pi+

c the same for the vertex 'nu' and we will store the matrix in
c cJotabisproton (which is J^\nu)
      call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &     xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &     xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &     cjota1bis,cjota2bis,cjota3bis,cjota4bis,cjota5bis,cjota6bis,
     &     cjota7bis,cJotabisproton)

c call Dagger and multiply by gamma0 from the left and from the
c right in order to obtain cJotatildeproton
      call Dagger(cJotabisproton,cAUX1p)
c      call MatrixProduct(cga(1,1,0),cAUX1p,cAUX2p)
c      call MatrixProduct(cAUX2p,cga(1,1,0),cJotatildeproton)
      call MatrixProductg0Mg0(cAUX1p,cJotatildeproton)

csame for the current cjota1bis for the substraction of
c DeltaPoleXDeltaPole
      call Dagger(cjota1bis,cAUX1p)
c      call MatrixProduct(cga(1,1,0),cAUX1p,cAUX2p)
c      call MatrixProduct(cAUX2p,cga(1,1,0),cjota1tildeproton)
      call MatrixProductg0Mg0(cAUX1p,cjota1tildeproton)

c  don't calculate the substraction of NPtimesNP because for 
c protons, cjota3bis goes multiplied by zero (see xCNucleonPoleproton)

c So, we must calculate the trace of [J^\mu*(cpM)*Jtilde^\nu*(cpprimeM)]
      call MatrixProduct(cJotaproton,cpM,cAUX3)
      call MatrixProduct(cAUX3,cJotatildeproton,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAprotontotal)
      call MatrixPt(cAUX4,cpprimeM,cAprotontotal)


c And must calculate the same for the term DeltaPoleXDeltaPole
      call MatrixProduct(cjota1,cpM,cAUX3)
      call MatrixProduct(cAUX3,cjota1tildeproton,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAproton11)
       call MatrixPt(cAUX4,cpprimeM,cAproton11)

      cAproton=cAprotontotal-cAproton11


c Secondly, for neutrons, the current J^\mu for the process n-->n + pi+
c Secondly, the isospin factors for n --> n + pi+
      xCdeltaneutron=1.d0/3.d0
      xCCrossedDeltaneutron=3.d0
      xCNucleonPoleneutron=1.d0
      xCCrossedNucleonPoleneutron=0.d0
      xCContactTermneutron=-1.d0
      xCPionPoleneutron=-1.d0
      xCPionInFlightneutron=-1.d0

      call ScalarProductMatrix(xCdeltaneutron,cjota1,cjota1neutron)
      call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &     xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &     xCcontactTermneutron,xCPionPoleneutron,xCPionInFlightneutron,
     &     cjota1,cjota2,cjota3,cjota4,cjota5,cjota6,cjota7,
     &     cJotaneutron)
c  in cJotaneutron we have stored the matrix J^\mu for n --> n + pi+

c the same for the vertex 'nu' and we will store the matrix in
c cJotabisneutron (which is J^\nu)
      call ScalarProductMatrix(xCdeltaneutron,cjota1bis,
     &     cjota1bisneutron)
      call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &     xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &     xCcontactTermneutron,xCPionPoleneutron,xCPionInFlightneutron,
     &     cjota1bis,cjota2bis,cjota3bis,cjota4bis,cjota5bis,cjota6bis,
     &     cjota7bis,cJotabisneutron)
c call Dagger and multiply by gamma0 from the left and from the
c right in order to obtain cJotatildeneutron
      call Dagger(cJotabisneutron,cAUX1n)
c      call MatrixProduct(cga(1,1,0),cAUX1n,cAUX2n)
c      call MatrixProduct(cAUX2n,cga(1,1,0),cJotatildeneutron)
      call MatrixProductg0Mg0(cAUX1n,cJotatildeneutron)

csame for the current cjota1bisneutron for the substraction
c of DeltaPoleXDeltaPole
      call Dagger(cjota1bisneutron,cAUX1n)
c      call MatrixProduct(cga(1,1,0),cAUX1n,cAUX2n)
c      call MatrixProduct(cAUX2n,cga(1,1,0),cjota1tildeneutron)
      call MatrixProductg0Mg0(cAUX1n,cjota1tildeneutron)


csame for the current cjota3bis for the substraction of
c NPtimesNP
      call Dagger(cjota3bis,cAUX1n)
c      call MatrixProduct(cga(1,1,0),cAUX1n,cAUX2n)
c      call MatrixProduct(cAUX2n,cga(1,1,0),cjota3tildeneutron)
      call MatrixProductg0Mg0(cAUX1n,cjota3tildeneutron)

c So, we must calculate the trace of [J^\mu*(cpM)*Jtilde^\nu*(cpprimeM)]
      call MatrixProduct(cJotaneutron,cpM,cAUX3)
      call MatrixProduct(cAUX3,cJotatildeneutron,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutrontotal)
      call MatrixPt(cAUX4,cpprimeM,cAneutrontotal)

c Andsame for the term DeltaPoleXDeltaPole
      call MatrixProduct(cjota1neutron,cpM,cAUX3)
      call MatrixProduct(cAUX3,cjota1tildeneutron,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutron11)
       call MatrixPt(cAUX4,cpprimeM,cAneutron11)

c Andsame for the term NPtimesNP
      call MatrixProduct(cjota3,cpM,cAUX3)
      call MatrixProduct(cAUX3,cjota3tildeneutron,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutron33)
      call MatrixPt(cAUX4,cpprimeM,cAneutron33)

      cAneutronpiplus=cAneutrontotal-cAneutron11-cAneutron33


c Thirdly, for the process n + W+ --> p + pi0
      xconstante=-1.d0/dsqrt(2.d0)
c We need J^\mu for this process
      call MatrixResta(cJotaproton,cJotaneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cJotapi0)
c We also need Jtilde^\nu for this process
      call MatrixResta(cJotatildeproton,cJotatildeneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cJotatildepi0)

c We need also j1^\mu
      call MatrixResta(cjota1,cjota1neutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cjota1pi0)
c We also need j1tilde^\nu
      call MatrixResta(cjota1tildeproton,cjota1tildeneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cjota1tildepi0)

c We need also j3^\mu
      call ScalarProductMatrix(xCNucleonPoleproton,cjota3,cjota3proton)
      call MatrixResta(cjota3proton,cjota3,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cjota3pi0)
c We also need j3tilde^\nu
      call ScalarProductMatrix(xCNucleonPoleproton,cjota3bis,
     &     cjota3tildeproton)
      call MatrixResta(cjota3tildeproton,cjota3tildeneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cjota3tildepi0)

ctraces
      call MatrixProduct(cJotapi0,cpM,cAUX3)
      call MatrixProduct(cAUX3,cJotatildepi0,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutronpi0total)
      call MatrixPt(cAUX4,cpprimeM,cAneutronpi0total)

      call MatrixProduct(cjota1pi0,cpM,cAUX3)
      call MatrixProduct(cAUX3,cjota1tildepi0,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutronpi011)
      call MatrixPt(cAUX4,cpprimeM,cAneutronpi011)

      call MatrixProduct(cjota3pi0,cpM,cAUX3)
      call MatrixProduct(cAUX3,cjota3tildepi0,cAUX4)
c      call MatrixProduct(cAUX4,cpprimeM,cAUX5)
c      call trace(cAUX5,cAneutronpi033)
      call MatrixPt(cAUX4,cpprimeM,cAneutronpi033)

      cAneutronpi0=cAneutronpi0total-cAneutronpi011-cAneutronpi033


      return
      end

*********************************************************************
*********************************************************************
*********************************************************************
c INPUT VARIABLES
c mu is a Lorentz index
c nu is the other Lorentz index
c p is the momentum of the first bubble
c p1 is the momentum of the second bubble
c q is the momentum of the W+
c xkpion is the momentum of the propagating first pion
c rho is the density in fm^(-3) that we will need for some currents

c OUTPUT VARIABLES 
c cAneutronpiplus is the complex tensorial component for the given Lorentz
c indices and momenta configuration which involves the currents
c J(n-->p + pi0) and J(n-->n + pi+)

c cAprotonpiplus is the complex tensorial component for the given Lorentz
c indices and momenta configuration which involves the currents
c J(n-->p + pi0) and J(p-->p + pi+)

      subroutine TensorADifferentBubbles(mu,nu,p,p1,q,xkpion,
     &     rho,cAneutronpiplus,cAprotonpiplus)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),p1(0:3),q(0:3),xkpion(0:3),pplusq(0:3)
      dimension pprime(0:3),qminuskpion(0:3),p1prime(0:3)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cga(4,4,0:3),cga5(4,4),cIdent(4,4)
      dimension cpslash(4,4),cqminuskpionslash(4,4),cpprimeslash(4,4)
      dimension ckpionslash(4,4),cp1slash(4,4),cp1primeslash(4,4)
      dimension cpM(4,4),cpprimeM(4,4),cp1M(4,4),cp1primeM(4,4)
      dimension cqminuskpiong5(4,4),ckpiong5(4,4),cjota1mu(4,4)
      dimension cjota2mu(4,4),cjota3mu(4,4),cjota4mu(4,4),cjota5mu(4,4)
      dimension cjota6mu(4,4),cjota7mu(4,4),cjota1nu(4,4),cjota2nu(4,4)
      dimension cjota3nu(4,4),cjota4nu(4,4),cjota5nu(4,4),cjota6nu(4,4)
      dimension cjota7nu(4,4),cJmuproton(4,4),cJnuproton(4,4)
      dimension cjota1muneutron(4,4),cjota2muneutron(4,4)
      dimension cJmuneutron(4,4),cjota1nuneutron(4,4),cAUXpizero(4,4)
      dimension cjota2nuneutron(4,4),cJnuneutron(4,4),cJmupi0(4,4)
      dimension cJnupi0(4,4),cJnupi0dagger(4,4),cJtildenupi0(4,4)
      dimension cjota1mup1(4,4),cjota2mup1(4,4),cjota3mup1(4,4)
      dimension cjota4mup1(4,4),cjota5mup1(4,4),cjota6mup1(4,4)
      dimension cjota7mup1(4,4)
      dimension cjota1nup1(4,4),cjota2nup1(4,4),cjota3nup1(4,4)
      dimension cjota4nup1(4,4),cjota5nup1(4,4),cjota6nup1(4,4)
      dimension cjota7nup1(4,4),cJmuprotonp1(4,4),cJnuprotonp1(4,4)
      dimension cJnuprotonp1dagger(4,4),cAUX(4,4),cJtildenuprotonp1(4,4)
      dimension cjota1mup1neutron(4,4),cjota2mup1neutron(4,4)
      dimension cJmuneutronp1(4,4),cjota1nup1neutron(4,4)
      dimension cjota2nup1neutron(4,4),cJnuneutronp1(4,4)
      dimension cJnuneutronp1dagger(4,4),cJtildenuneutronp1(4,4)
      dimension cAUX1(4,4),cAUX2(4,4)
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common /matrices/ cga,cga5,cIdent
      common /options/iopt,isospin



c first we need pprime=p+q-xkpion
      call SumFourVectors(p,q,pplusq)
      call RestaFourVectors(pplusq,xkpion,pprime)
c secondly, we need qminuskpion=q-xkpion
      call RestaFourVectors(q,xkpion,qminuskpion)
c thirdly, we need p1prime=p1+xkpion
      call SumFourVectors(p1,xkpion,p1prime)

c  need to call slash for some combination of four-momenta
      call slash(p,cpslash)
      call slash(qminuskpion,cqminuskpionslash)
      call slash(pprime,cpprimeslash)
      call slash(xkpion,ckpionslash)
      call slash(p1,cp1slash)
      call slash(p1prime,cp1primeslash)

c need (pslash+M),(pprimeslash+M),(p1slash+M),(p1primeslash+M)
      call MatrixSum(cpslash,cMidentity,cpM)
      call MatrixSum(cpprimeslash,cMidentity,cpprimeM)
      call MatrixSum(cp1slash,cMidentity,cp1M)
      call MatrixSum(cp1primeslash,cMidentity,cp1primeM)

c need to multiply some slash by gamma5, for instance 
c cqminuskpionslash and ckpionslash
      call MatrixProductg5(cqminuskpionslash,cqminuskpiong5)
      call MatrixProductg5(ckpionslash,ckpiong5)

c firstly, we will need all the currents j1,j2,... evaluated at 
c momenta configuration (p,q,xkpion)

c  put the isospin factors, for having them in mind
c First the isospin factors for p --> p + pi+
      xCdeltaproton=1.d0
      xCCrossedDeltaproton=1.d0
      xCNucleonPoleproton=0.d0
      xCCrossedNucleonPoleproton=1.d0
      xCcontactTermproton=1.d0
      xCPionPoleproton=1.d0
      xCPionInFlightproton=1.d0

c Secondly, the isospin factors for n --> n + pi+
      xCdeltaneutron=1.d0/3.d0
      xCCrossedDeltaneutron=3.d0
      xCNucleonPoleneutron=1.d0
      xCCrossedNucleonPoleneutron=0.d0
      xCContactTermneutron=-1.d0
      xCPionPoleneutron=-1.d0
      xCPionInFlightneutron=-1.d0


c call all the currents with indices mu and nu and momenta
c configuration (p,q,xkpion)
      call CurrentDeltaPole2body(mu,xkpion,p,q,rho,cjota1mu)
      call CurrentCrossedDeltaPole2body(mu,xkpion,p,q,cjota2mu)
      call CurrentNucleonPole2body(mu,xkpion,p,q,rho,cjota3mu)
      call CurrentCrossedNucleonPole2body(mu,xkpion,p,q,cjota4mu)
      call CurrentContactTerm(mu,xkpion,p,q,cjota5mu)
      call CurrentPionPole(mu,xkpion,p,q,cjota6mu)
      call CurrentPionInFlight(mu,xkpion,p,q,cjota7mu)

c have to construct J^\mu(p,q,xkpion) for the process
c p --> p + pi+. This matrix will be stored in cJmuproton
      call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &     xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &     xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &     cjota1mu,cjota2mu,cjota3mu,cjota4mu,cjota5mu,cjota6mu,
     &     cjota7mu,cJmuproton)

c have to construct J^\mu(p,q,xkpion) for the process 
c n --> n + pi+. This matrix will be stored in cJmuneutron
      call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &     xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &     xCcontactTermneutron,xCPionPoleneutron,xCPionInFlightneutron,
     &     cjota1mu,cjota2mu,cjota3mu,cjota4mu,cjota5mu,cjota6mu,
     &     cjota7mu,cJmuneutron)



      if(nu .eq. mu) then  !! We don't need to call the currents again
c     because we have them in the cjota1mu, cjota2mu,...
         do i=1,4
            do j=1,4
               cJnuproton(i,j)=cJmuproton(i,j)
               cJnuneutron(i,j)=cJmuneutron(i,j)
            enddo
         enddo
         
      else  !!! We need to call the currents again with index nu

         call CurrentDeltaPole2body(nu,xkpion,p,q,rho,cjota1nu)
         call CurrentCrossedDeltaPole2body2(nu,xkpion,p,q,cjota2nu)
         call CurrentNucleonPole2body(nu,xkpion,p,q,rho,cjota3nu)
         call CurrentCrossedNucleonPole2body(nu,xkpion,p,q,cjota4nu)
         call CurrentContactTerm(nu,xkpion,p,q,cjota5nu)
         call CurrentPionPole(nu,xkpion,p,q,cjota6nu)
         call CurrentPionInFlight(nu,xkpion,p,q,cjota7nu)
         
c have to construct J^\nu(p,q,xkpion) for the process
c p --> p + pi+. This matrix will be stored in cJnuproton
         call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &        xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &        xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &        cjota1nu,cjota2nu,cjota3nu,cjota4nu,cjota5nu,cjota6nu,
     &        cjota7nu,cJnuproton)

c have to construct J^\nu(p,q,xkpion) for the process
c n --> n + pi+. This matrix will be stored in cJnuneutron
         call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &        xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &        xCcontactTermneutron,xCPionPoleneutron,
     &        xCPionInFlightneutron,cjota1nu,cjota2nu,cjota3nu,cjota4nu,
     &        cjota5nu,cjota6nu,cjota7nu,cJnuneutron)

      endif


c can calculate J^\mu(p,q,xkpion) for the process n --> p + pi0
c and J^\nu(p,q,xkpion) for the same process. Let's remember that
c <p pi0|J|n>=-1/Sqrt[2]*(<p pi+|J|p>-<n pi+|J|n>). THESE ARE THE 
c CURRENTS WHICH APPEAR IN THE TENSOR A(mu,nu)

      xconstante=-1.d0/dsqrt(2.d0)
      call MatrixResta(cJmuproton,cJmuneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cJmupi0)

      call MatrixResta(cJnuproton,cJnuneutron,cAUXpizero)
      call ScalarProductMatrix(xconstante,cAUXpizero,cJnupi0)
c But for the index 'nu', we need Jtilde=gamma0*J^\dagger*gamma0
      call Dagger(cJnupi0,cJnupi0dagger)
c      call MatrixProduct(cga(1,1,0),cJnupi0dagger,cAUXpizero)
c      call MatrixProduct(cAUXpizero,cga(1,1,0),cJtildenupi0)
      call MatrixProductg0Mg0(cJnupi0dagger,cJtildenupi0)


c secondly, we will need all the currents j1,j2,... evaluated at 
c momenta configuration (p1,q,q-kpion)

c call all the currents with indices mu and nu and momenta
c configuration (p1,q,qminuskpion)
      call CurrentDeltaPole2body(mu,qminuskpion,p1,q,rho,cjota1mup1)
      call CurrentCrossedDeltaPole2body(mu,qminuskpion,p1,q,cjota2mup1)
      call CurrentNucleonPole2body(mu,qminuskpion,p1,q,rho,cjota3mup1)
      call CurrentCrossedNucleonPole2body(mu,qminuskpion,p1,q,
     &     cjota4mup1)
      call CurrentContactTerm(mu,qminuskpion,p1,q,cjota5mup1)
      call CurrentPionPole(mu,qminuskpion,p1,q,cjota6mup1)
      call CurrentPionInFlight(mu,qminuskpion,p1,q,cjota7mup1)


c have to construct J^\mu(p1,q,q-kpion) for the process
c p --> p + pi+. This matrix will be stored in cJmuprotonp1
      call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &     xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &     xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &     cjota1mup1,cjota2mup1,cjota3mup1,cjota4mup1,cjota5mup1,
     &     cjota6mup1,cjota7mup1,cJmuprotonp1)

c have to construct J^\mu(p1,q,q-kpion) for the process
c n --> n + pi+. This matrix will be stored in cJmuneutronp1
      call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &     xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &     xCcontactTermneutron,xCPionPoleneutron,xCPionInFlightneutron,
     &     cjota1mup1,cjota2mup1,cjota3mup1,cjota4mup1,cjota5mup1,
     &     cjota6mup1,cjota7mup1,cJmuneutronp1)

      if(nu .eq. mu)then !! We don't need to call all the currents again
c because we already have them
         do i=1,4
            do j=1,4
               cJnuprotonp1(i,j)=cJmuprotonp1(i,j)
               cJnuneutronp1(i,j)=cJmuneutronp1(i,j)
            enddo
         enddo

      else !!! We need to call the currents with index 'nu'
         
         call CurrentDeltaPole2body(nu,qminuskpion,p1,q,rho,cjota1nup1)
         call CurrentCrossedDeltaPole2body2(nu,qminuskpion,p1,q,
     &        cjota2nup1)
         call CurrentNucleonPole2body(nu,qminuskpion,p1,q,rho,
     &        cjota3nup1)
         call CurrentCrossedNucleonPole2body(nu,qminuskpion,p1,q,
     &        cjota4nup1)
         call CurrentContactTerm(nu,qminuskpion,p1,q,cjota5nup1)
         call CurrentPionPole(nu,qminuskpion,p1,q,cjota6nup1)
         call CurrentPionInFlight(nu,qminuskpion,p1,q,cjota7nup1)

c have to construct J^\nu(p1,q,q-kpion) for the process
c p --> p + pi+. This matrix will be stored in cJnuprotonp1
         call SumCurrents(xCdeltaproton,xCCrossedDeltaproton,
     &        xCNucleonPoleproton,xCCrossedNucleonPoleproton,
     &        xCcontactTermproton,xCPionPoleproton,xCPionInFlightproton,
     &        cjota1nup1,cjota2nup1,cjota3nup1,cjota4nup1,cjota5nup1,
     &        cjota6nup1,cjota7nup1,cJnuprotonp1)

c  have to construct J^\nu(p1,q,q-kpion) for the process
c n --> n + pi+. This matrix will be stored in cJnuneutronp1
         call SumCurrents(xCdeltaneutron,xCCrossedDeltaneutron,
     &        xCNucleonPoleneutron,xCCrossedNucleonPoleneutron,
     &        xCcontactTermneutron,xCPionPoleneutron,
     &        xCPionInFlightneutron,cjota1nup1,cjota2nup1,cjota3nup1,
     &        cjota4nup1,cjota5nup1,cjota6nup1,cjota7nup1,cJnuneutronp1)
      endif


c But I wanted J^\nu(p1,q,q-kpion) for the process p --> p + pi+
c to obtain cJtildenuprotonp1
      call Dagger(cJnuprotonp1,cJnuprotonp1dagger)
c      call MatrixProduct(cga(1,1,0),cJnuprotonp1dagger,cAUX)
c      call MatrixProduct(cAUX,cga(1,1,0),cJtildenuprotonp1)
      call MatrixProductg0Mg0(cJnuprotonp1dagger,cJtildenuprotonp1)

c SO, THE MATRICES WHICH APPEAR IN A(mu,nu) ARE cJmuprotonp1 AND
c cJtildenuprotonp1

c But I wanted J^\nu(p1,q,q-kpion) for the process n --> n + pi+
c to obtain Jtilde^\nu(p1,q,q-kpion)
      call Dagger(cJnuneutronp1,cJnuneutronp1dagger)
c      call MatrixProduct(cga(1,1,0),cJnuneutronp1dagger,cAUX)
c      call MatrixProduct(cAUX,cga(1,1,0),cJtildenuneutronp1)
      call MatrixProductg0Mg0(cJnuneutronp1dagger,cJtildenuneutronp1)



c SO, THE MATRICES WHICH APPEAR IN A(mu,nu) ARE cJmuneutronp1 AND
c cJtildenuneutronp1


cTRACES, FIRST FOR A(mu,nu) n --> p + pi0 and n --> n + pi+
      call MatrixProduct(cJmupi0,cpM,cAUX)
      call MatrixProduct(cAUX,cqminuskpiong5,cAUX1)
c      call MatrixProduct(cAUX1,cpprimeM,cAUX2)
c      call trace(cAUX2,ctraza1)
      call MatrixPt(cAUX1,cpprimeM,ctraza1)

      call MatrixProduct(ckpiong5,cp1M,cAUX)
      call MatrixProduct(cAUX,cJtildenuneutronp1,cAUX1)
c      call MatrixProduct(cAUX1,cp1primeM,cAUX2)
c      call trace(cAUX2,ctraza2)
      call MatrixPt(cAUX1,cp1primeM,ctraza2)
      
      call MatrixProduct(cJmuneutronp1,cp1M,cAUX)
      call MatrixProduct(cAUX,ckpiong5,cAUX1)
c      call MatrixProduct(cAUX1,cp1primeM,cAUX2)
c      call trace(cAUX2,ctraza3)
      call MatrixPt(cAUX1,cp1primeM,ctraza3)

      call MatrixProduct(cqminuskpiong5,cpM,cAUX)
      call MatrixProduct(cAUX,cJtildenupi0,cAUX1)
c      call MatrixProduct(cAUX1,cpprimeM,cAUX2)
c      call trace(cAUX2,ctraza4)
      call MatrixPt(cAUX1,cpprimeM,ctraza4)
c
      cAneutronpiplus=ctraza1*ctraza2+ctraza3*ctraza4

c bNEW ISOSPIN
         if(isospin.eq.2)then
          cAneutronpiplus=0.
         endif
c eNEW ISOSPIN

c Secondly, for A(mu,nu) n --> p + pi0 and p --> p + pi+
c the first trace is the same than ctraza1

csecond trace
      call MatrixProduct(ckpiong5,cp1M,cAUX)
      call MatrixProduct(cAUX,cJtildenuprotonp1,cAUX1)
c      call MatrixProduct(cAUX1,cp1primeM,cAUX2)
c      call trace(cAUX2,ctraza2)
      call MatrixPt(cAUX1,cp1primeM,ctraza2)

cthird trace
      call MatrixProduct(cJmuprotonp1,cp1M,cAUX)
      call MatrixProduct(cAUX,ckpiong5,cAUX1)
c      call MatrixProduct(cAUX1,cp1primeM,cAUX2)
c      call trace(cAUX2,ctraza3)
      call MatrixPt(cAUX1,cp1primeM,ctraza3)

c the fourth trace is the same than ctraza4

      cAprotonpiplus=ctraza1*ctraza2+ctraza3*ctraza4


      return
      end



*******************************************************************
c  define Manolo's version of the imaginary part
c of the Lindhard function

c INPUT VARIABLES
c 'q' is the four-momentum
c xkFh is the local Fermi momentum of the hole
c xkFp is the local Fermi momentum of the particle

c OUTPUT VARIABLES
c the function AImUr(q,xkfh,xkfp)

      real*8 function AImUr(q,xkfh,xkfp)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3)
      parameter(xM=0.94d0,pi=3.141592653589793d0)

      AImUr=0.d0
c First, we need q2, so we must call LorentzScalarProduct
      call LorentzScalarProduct(q,q,q2)

      if(q2.gt.0.d0)return
      if(q(0).lt.0.d0)return
c We also need the Fermi energy corresponding to Fermi momentum for
c holes and particles
      Efh=dsqrt(xM**2+xkfh**2)
      Efp=dsqrt(xM**2+xkfp**2)

c We also need the modulus of q (three-vector)
      xmodq=dsqrt(q(1)**2+q(2)**2+q(3)**2)

c And we also need epsilonRparticle, which is the maximum of
c three quantities
      epsilonRp=max(xM,Efp-q(0),(-q(0)+xmodq*dsqrt(1.d0-4.d0*xM**2/q2))/
     &     2.d0)

      if ((Efh-Efp+q(0)).lt.0.d0)return
      if ((Efh-epsilonRp).lt.0.d0)return

      AImUr=-xM**2/(2.d0*pi*xmodq)*(Efh-epsilonRp)

      return
      end
*********************************************************************

*******************************************************************
*******************************************************************
*******************************************************************
c first of all we  define the free pion propagator
      real*8 function Dpion(xkpion)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3)
      parameter(xmpion=0.13957d0)

      gprime=0.d0
      
      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
c      xmodkpi=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)!!modulus of kpion


      Dpion=1.d0/(xkpion2-xmpion**2)
      
      return
      end
*******************************************************************
*******************************************************************
*******************************************************************
cform factor for pions off-shell Fpi(q) from the work by
c A. Gil et al
      real*8 function Fpi(q)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3)
c constants of interest in GeV
      xlambdapi=1.2d0 !!! ~1200 MeV
      xmpion=0.13957d0 !!! 139.57 MeV

      call LorentzScalarProduct(q,q,q2)

      Fpi=(xlambdapi**2-xmpion**2)/(xlambdapi**2-q2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cform factor for rho mesons off-shell FormFactorRho(k)
      real*8 function FormFactorRho(xkrho)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkrho(0:3)
c constants of interest in GeV
      xlambdarho=2.5d0 !!! 2500 MeV
      xmrho=0.77d0     !!! 770 MeV

      call LorentzScalarProduct(xkrho,xkrho,xkrho2)
      
      FormFactorRho=(xlambdarho**2-xmrho**2)/(xlambdarho**2-xkrho2)
      return
      end
*******************************************************************
*******************************************************************
c  define the subroutine HadronicTensor2body,
c for the absorption of a W+ by 2 nucleons W+NN --> NN
      subroutine HadronicTensor2body(mu,nu,q,cW2body)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3),xkpi0d(2000),phid(2000),thetad(2000)
      dimension xmodkpid(2000),rd(2000),crkpi0d(2000),crphid(2000)
      dimension crthetad(2000),crmodkpid(2000),crrd(2000)
      dimension dcosid(2000),crdcosid(2000)
      dimension pnum(0:3),xkpion(0:3),qminuskpion(0:3),vectora(1:3)
      dimension pden(0:3)
       common /nucleus/rho0,a,th,na
       common/densidad/dro,dmdf   !!!in fermis
       common/datos/dpi,dm,hbb,hbbi
       common /options/iopt,isospin

c new
      common/matrizqslash/cMdeltaMatrix,cMidentity
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cProjector(4,4)
c new
c new
      dimension pdeltax(0:3),cdproj(0:3,0:3,4,4)
      dimension cpdeltaslashx(4,4)
      common/deltapole2/cdproj
!$OMP THREADPRIVATE(/deltapole2/)
c      common/deltapole2x/cdprojx(0:3,0:3,4,4)
c new




c constants of interest for the subroutine
      pi=3.141592653589793d0
      xM=0.94d0
      gA=1.26d0
      fpion=0.093d0
      hbarc=0.19733d0
      xmpion=0.13957d0
c constants for DIMASM
      DM=4.7559925D0   !!! Nucleon mass in fermis
      DPI=3.141592654D0
      hbb=197.33d0
      hbbi=1.d0/hbb
      xmodq=dsqrt(q(1)**2+q(2)**2+q(3)**2)
      

c if q(0) is negative, then the step function that goes in front of
c the integrals is zero and the Hadronic Tensor too.
      if (q(0).lt.0.d0)then
         cW2body=0.d0
         return
      endif

c precision in the integrals
      nr=1
      ncosid=2
      nmodkpi=1
      ntheta=1
      nphi=1
      nkpi0=1
c
c Upper limits for the integrations
      rlimfermis=5.d0*a
      rlimGeV=rlimfermis/hbarc !!! in GeV^(-1)

      xcutoff=2.0d0  !!!cut-off for the modulus of kpion in GeV

c CALL to subroutines to make the partitions

c      call DSG20R(0.d0,rlimGeV,nr,rd,nrr)
      call DSG20R(-1.d0,1.d0,ncosid,dcosid,ncosidcosid)
      call DSG20R(0.d0,xcutoff,nmodkpi,xmodkpid,nmodkpimodkpi)
c      call DSG20R(0.d0,pi,ntheta,thetad,nthetatheta)
c      call DSG20R(0.d0,2.d0*pi,nphi,phid,nphiphi)
c      call DSG20R(0.d0,q(0),nkpi0,xkpi0d,nkpi0kpi0)
      call STGS8(0.d0,rlimGeV,rd)
c      call STGS8(-1.d0,1.d0,dcosid)
c      call STGS8(0.d0,xcutoff,xmodkpid)
      call STGS8(0.d0,pi,thetad)
      call STGS8(0.d0,2.d0*pi,phid)
      call STGS8(0.d0,q(0),xkpi0d)

      do ir=1,8
         rGeV=rd(ir)
         rfermi=rGeV*hbarc
         rhofermi=density(rfermi) !!! density in fm^(-3)
         dro=rhofermi

         xkfermi=(3.d0*pi**2/2.d0*rhofermi)**(1.d0/3.d0) !!! momentum of Fermi in fm^(-1)
         dmdf=xkfermi
         xkFGeV=xkfermi*hbarc  !!! Momentum of Fermi in GeV
c can calculate the averaged <pmom>
         pmom_averaged=dsqrt(3.d0/5.d0)*xkFGeV
c         pmom_averaged=0.d0
c And the energy that corresponds to this averaged <pmom>
	 energyfermi=dsqrt(xM**2+xkFGeV**2)-xM
         pnum(0)=dsqrt(xM**2+pmom_averaged**2)-energyfermi*0.d0
         pden(0)=dsqrt(xM**2+pmom_averaged**2)-energyfermi

         do icosid=1,ncosidcosid  ! 8
            dcoseno=dcosid(icosid)
            seno=dsqrt(1.d0-dcoseno**2)
            pnum(1)=pmom_averaged*(dcoseno*q(1)/xmodq)
            pnum(2)=pmom_averaged*(dcoseno*q(2)/xmodq+seno)
            pnum(3)=pmom_averaged*(dcoseno*q(3)/xmodq)

            pden(1)=pmom_averaged*(dcoseno*q(1)/xmodq)
            pden(2)=pmom_averaged*(dcoseno*q(2)/xmodq+seno)
            pden(3)=pmom_averaged*(dcoseno*q(3)/xmodq)


cccc  test change
      call SumFourVectors(pden,q,pDeltax)
      call slash(pdeltax,cpdeltaslashx)
      call MatrixSum(cpdeltaslashx,cMdeltaMatrix,cpdMd)
           do i11=0,3
           do i22=0,3
           call DeltaProjector(i11,i22,pdeltax,
     &           cpdMd,cProjector)
              do ii=1,4
              do jj=1,4
           cdproj(i11,i22,ii,jj)=cProjector(ii,jj)
              enddo
              enddo
           enddo
           enddo
cccc  end test change


            do imodkpi=1,nmodkpimodkpi   !8
               xmodkpi=xmodkpid(imodkpi) !!! Modulus of kpi

               do itheta=1,8
                  theta=thetad(itheta)
                  xkpion(3)=xmodkpi*dcos(theta)

                  do iphi=1,8
                     phi=phid(iphi)
                     xkpion(1)=xmodkpi*dsin(theta)*dcos(phi)
                     xkpion(2)=xmodkpi*dsin(theta)*dsin(phi)

                     do ikpi0=1,8
                        xkpion(0)=xkpi0d(ikpi0)
                     
                        call LorentzScalarProduct(xkpion,xkpion,xkpion2)
                        call PolarizationPion(xkpion,rhofermi,cpolariz)
                        xImDpion=xImPolarizationPion(xkpion,rhofermi)/
     &                       (zabs(xkpion2-xmpion**2-cpolariz))**2
                        formfactor=Fpi(xkpion)
                        call RestaFourVectors(q,xkpion,qminuskpion)
c                     xImLind1=AImUr(xkpion,xkFGeV,xkFGeV)
                        xImLind2=AImUr(qminuskpion,xkFGeV,xkFGeV)
c                        call CrossProduct(q,xkpion,vectora)

c                        do jj=1,3
c                           pnum(jj)=pmom_averaged*vectora(jj)
c                           pden(jj)=pmom_averaged*vectora(jj)
c                        enddo

      if((abs(xImLind2).gt.1.d-20).and.(abs(xImDpion).gt.1.d-20))then
                           call TensorABackground2body(mu,nu,xkpion,
     &                          pnum,pden,q,rhofermi,
     &                    cAproton,cAneutronpiplus,cAneutronpi0)
                        else
                           cAproton=0.d0
                           cAneutronpiplus=0.d0
                           cAneutronpi0=0.d0
                        endif
c                     crkpi0d(ikpi0)=xkpion2*propagator**2*xImLind1*
c     &      xImLind2*formfactor**4*(cAproton+cAneutronpiplus+
c     &      cAneutronpi0)

C bNEW ISOSPIN
       if(isospin.eq.2)then
        cAproton=caproton  ! leads to pp for neutrino
        cAneutronpiplus=0.  ! leads to pn  
        cAneutronpi0=cAneutronpi0/2.  ! w+ n --> pi0 +p ; and pi0 absorbed
c                                    ! half of times by p
       endif
C eNEW ISOSPIN



                        crkpi0d(ikpi0)=xImDpion*formfactor**2*xImLind2*
     &                       (cAproton+cAneutronpiplus+cAneutronpi0)
                     enddo
c                     call DRG20C(0.d0,q(0),nkpi0,crkpi0d,cres)
                     call REGS8C(0.d0,q(0),crkpi0d,cres)
                     crphid(iphi)=cres
                  enddo
c                  call DRG20C(0.d0,2.d0*pi,nphi,crphid,cres)
                  call REGS8C(0.d0,2.d0*pi,crphid,cres)
                  crthetad(itheta)=cres*dsin(theta)
               enddo
c              call DRG20C(0.d0,pi,ntheta,crthetad,cres)
               call REGS8C(0.d0,pi,crthetad,cres)
               crmodkpid(imodkpi)=cres*xmodkpi**2
            enddo
            call DRG20C(0.d0,xcutoff,nmodkpi,crmodkpid,cres)
c            call REGS8C(0.d0,xcutoff,crmodkpid,cres)
            crdcosid(icosid)=cres/2.d0
         enddo
         call DRG20C(-1.d0,1.d0,ncosid,crdcosid,cres)
c         call REGS8C(-1.d0,1.d0,crdcosid,cres)
         crrd(ir)=cres*rGeV**2
c
      enddo
c      call DRG20C(0.d0,rlimGeV,nr,crrd,cresultado)
      call REGS8C(0.d0,rlimGeV,crrd,cresultado)
     
c      factor=-1.d0/xM**2*(gA/(2.d0*fpion))**2*4.d0*pi/(2.d0*pi)**5
       factor=1.d0/(2.d0*xM**2)*4.d0*pi/(2.d0*pi)**5

       cW2body=factor*cresultado

       return
       end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
crho meson propagator. xk is given in GeV
      real*8 function Drho(xk)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xk(0:3)

      xmrho=0.77d0 !!! in GeV

      call LorentzScalarProduct(xk,xk,xk2)
      
      Drho=1.d0/(xk2-xmrho**2)
      return
      end
*******************************************************************
*******************************************************************
c  define the subroutine HadronicTensorRho,
c for the absorption of a W+ by 2 nucleons W+NN --> NN via rho exchange
      subroutine HadronicTensorRho(mu,nu,q,cW2body)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3),xkpi0d(2000),phid(2000),thetad(2000)
      dimension xmodkpid(2000),rd(2000),crkpi0d(2000),crphid(2000)
      dimension crthetad(2000),crmodkpid(2000),crrd(2000)
      dimension dcosid(2000),crdcosid(2000)
      dimension pnum(0:3),xkpion(0:3),qminuskpion(0:3),vectora(1:3)
      dimension pden(0:3)
      common /nucleus/rho0,a,th,na
      common/densidad/dro,dmdf  !!!in fermis
      common/datos/dpi,dm,hbb,hbbi
      common /options/iopt,isospin

c new
      common/matrizqslash/cMdeltaMatrix,cMidentity
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cProjector(4,4)
c new
c new
      dimension pdeltax(0:3),cdproj(0:3,0:3,4,4)
      dimension cpdeltaslashx(4,4)
      common/deltapole/cdproj
      common/deltapolex/cdprojx(0:3,0:3,4,4)
!$OMP THREADPRIVATE(/deltapole/,/deltapolex/)
c new

c constants of interest for the subroutine
      pi=3.141592653589793d0
      xM=0.94d0
      gA=1.26d0
      fpion=0.093d0
      hbarc=0.19733d0
      xcrho=2.d0
c constants for DIMASM
      DM=4.7559925D0   !!! Nucleon mass in fermis
      DPI=3.141592654D0
      hbb=197.33d0
      hbbi=1.d0/hbb
      xmodq=dsqrt(q(1)**2+q(2)**2+q(3)**2)
      

c if q(0) is negative, then the step function that goes in front of
c the integrals is zero and the Hadronic Tensor too.
      if (q(0).lt.0.d0)then
         cW2body=0.d0
         return
      endif

c precision in the integrals
      nr=1
      ncosid=2
      nmodkpi=1
      ntheta=1
      nphi=1
      nkpi0=1
c
c Upper limits for the integrations
      rlimfermis=5.d0*a
      rlimGeV=rlimfermis/hbarc !!! in GeV^(-1)

      xcutoff=2.0d0  !!!cut-off for the modulus of kpion in GeV

c CALL to subroutines to make the partitions

c      call DSG20R(0.d0,rlimGeV,nr,rd,nrr)
      call DSG20R(-1.d0,1.d0,ncosid,dcosid,ncosidcosid)
      call DSG20R(0.d0,xcutoff,nmodkpi,xmodkpid,nmodkpimodkpi)
c      call DSG20R(0.d0,pi,ntheta,thetad,nthetatheta)
c      call DSG20R(0.d0,2.d0*pi,nphi,phid,nphiphi)
c      call DSG20R(0.d0,q(0),nkpi0,xkpi0d,nkpi0kpi0)
      call STGS8(0.d0,rlimGeV,rd)
c      call STGS8(-1.d0,1.d0,dcosid)
c      call STGS8(0.d0,xcutoff,xmodkpid)
      call STGS8(0.d0,pi,thetad)
      call STGS8(0.d0,2.d0*pi,phid)
      call STGS8(0.d0,q(0),xkpi0d)

      do ir=1,8
         rGeV=rd(ir)
         rfermi=rGeV*hbarc
         rhofermi=density(rfermi) !!! density in fm^(-3)
         dro=rhofermi

         xkfermi=(3.d0*pi**2/2.d0*rhofermi)**(1.d0/3.d0) !!! momentum of Fermi in fm^(-1)
         dmdf=xkfermi
         xkFGeV=xkfermi*hbarc  !!! Momentum of Fermi in GeV
c can calculate the averaged <pmom>
         pmom_averaged=dsqrt(3.d0/5.d0)*xkFGeV
c         pmom_averaged=0.d0
c And the energy that corresponds to this averaged <pmom>
	 energyfermi=dsqrt(xM**2+xkFGeV**2)-xM
c         pnum(0)=dsqrt(xM**2+pmom_averaged**2)-energyfermi*0.d0
         pden(0)=dsqrt(xM**2+pmom_averaged**2)-energyfermi

         do icosid=1,ncosidcosid !8
            dcoseno=dcosid(icosid)
            seno=dsqrt(1.d0-dcoseno**2)
c            pnum(1)=pmom_averaged*(dcoseno*q(1)/xmodq)
c            pnum(2)=pmom_averaged*(dcoseno*q(2)/xmodq+seno)
c            pnum(3)=pmom_averaged*(dcoseno*q(3)/xmodq)

            pden(1)=pmom_averaged*(dcoseno*q(1)/xmodq)
            pden(2)=pmom_averaged*(dcoseno*q(2)/xmodq+seno)
            pden(3)=pmom_averaged*(dcoseno*q(3)/xmodq)

cccc  test change
      call SumFourVectors(pden,q,pDeltax)
      call slash(pdeltax,cpdeltaslashx)
      call MatrixSum(cpdeltaslashx,cMdeltaMatrix,cpdMd)
           do i11=0,3
           do i22=0,3
           call DeltaProjector(i11,i22,pdeltax,
     &           cpdMd,cProjector)
              do ii=1,4
              do jj=1,4
           cdproj(i11,i22,ii,jj)=cProjector(ii,jj)
              enddo
              enddo
           enddo
           enddo
cccc  end test change
           
            do imodkpi=1,nmodkpimodkpi  !8
               xmodkpi=xmodkpid(imodkpi) !!! Modulus of kpi
              do ikpi0=1,8
                        xkpion(0)=xkpi0d(ikpi0)
                       xk0fermis=xkpion(0)/hbarc !! in fm^(-1)
c must work in fermis in order to call to ULIND
                        xmodkfermis=xmodkpi/hbarc !! in fm^(-1)
                        call ULIND3(xk0fermis,xmodkfermis,rhofermi,
     &                       CUnuc)
                    if(abs(dimag(CUnuc)).gt.1.d-15)then
                        call ULIND2(xk0fermis,xmodkfermis,rhofermi,
     &                       CUfuncompleta)
cc Conversion of the Lindhard functions to GeV^2
                        CUfuncompletaGeV=CUfuncompleta*hbarc**2 !!!GeV²
                        CUnucGeV=CUnuc*hbarc**2    !!! in GeV²
               




               do itheta=1,8
                  theta=thetad(itheta)
                  xkpion(3)=xmodkpi*dcos(theta)

                  do iphi=1,8
                     phi=phid(iphi)
                     xkpion(1)=xmodkpi*dsin(theta)*dcos(phi)
                     xkpion(2)=xmodkpi*dsin(theta)*dsin(phi)


                   
                        propagator=Drho(xkpion)
                        formfactor=FormFactorRho(xkpion)
cc Conversion of the Lindhard functions to GeV^2
                        CUfuncompletaGeV=CUfuncompleta*hbarc**2 !!!GeV²
                        CUnucGeV=CUnuc*hbarc**2    !!! in GeV²

                        xImUlambda=dimag(CUnucGeV)/zabs(1.d0-
     &                       CUfuncompletaGeV*Vt(xkpion))**2

                        call RestaFourVectors(q,xkpion,qminuskpion)
c                     xImLind1=AImUr(xkpion,xkFGeV,xkFGeV)
                        xImLind2=AImUr(qminuskpion,xkFGeV,xkFGeV)


      if((abs(xImLind2).gt.1.d-20).and.(abs(xImUlambda).gt.1.d-20))then
cccc  test change
      call RestaFourVectors(pden,xkpion,pDeltax)
      call slash(pdeltax,cpdeltaslashx)
      call MatrixSum(cpdeltaslashx,cMdeltaMatrix,cpdMd)
           do i11=0,3
           do i22=0,3
           call DeltaProjector(i11,i22,pdeltax,
     &           cpdMd,cProjector)
              do ii=1,4
              do jj=1,4
           cdprojx(i11,i22,ii,jj)=cProjector(ii,jj)
              enddo
              enddo
           enddo
           enddo
cccc  end test change
 

                           call TensorArho(mu,nu,pden,q,xkpion,
     &                    cAproton,cAneutronpiplus,cAneutronpi0)
                        else
                           cAproton=0.d0
                           cAneutronpiplus=0.d0
                           cAneutronpi0=0.d0
                        endif

C bNEW ISOSPIN
       if(isospin.eq.2)then
       cAproton=caproton  ! leads to pp for neutrino
       cAneutronpiplus=0.  ! leads to pn  
       cAneutronpi0=cAneutronpi0/2.  ! w+ n --> pi0 +p ; and pi0 absorbed
c                                    ! half of times by p
       endif
C eNEW ISOSPIN



                        crphid(iphi)=formfactor**4*propagator**2*
     &                       xImLind2*xcrho*(gA/(2.d0*fpion))**2*
     &                       xmodkpi**2*xImUlambda*
     &                       (cAproton+cAneutronpiplus+cAneutronpi0)
                     enddo
c                  call DRG20C(0.d0,2.d0*pi,nphi,crphid,cres)
                  call REGS8C(0.d0,2.d0*pi,crphid,cres)
                  crthetad(itheta)=cres*dsin(theta)
               enddo
c               call DRG20C(0.d0,pi,ntheta,crthetad,cres)
               call REGS8C(0.d0,pi,crthetad,cres)
               crkpi0d(ikpi0)=cres
               else
               crkpi0d(ikpi0)=0.d0
               endif
               enddo
c                     call DRG20C(0.d0,q(0),nkpi0,crkpi0d,cres)
                     call REGS8C(0.d0,q(0),crkpi0d,cres)
               crmodkpid(imodkpi)=cres*xmodkpi**2
            enddo
            call DRG20C(0.d0,xcutoff,nmodkpi,crmodkpid,cres)
c            call REGS8C(0.d0,xcutoff,crmodkpid,cres)
            crdcosid(icosid)=cres/2.d0
         enddo
         call DRG20C(-1.d0,1.d0,ncosid,crdcosid,cres)
c         call REGS8C(-1.d0,1.d0,crdcosid,cres)
         crrd(ir)=cres*rGeV**2
c
      enddo
c      call DRG20C(0.d0,rlimGeV,nr,crrd,cresultado)
      call REGS8C(0.d0,rlimGeV,crrd,cresultado)
     
c      factor=-1.d0/xM**2*(gA/(2.d0*fpion))**2*4.d0*pi/(2.d0*pi)**5
       factor=-1.d0/(2.d0*xM**2)*4.d0*pi/(2.d0*pi)**5

       cW2body=factor*cresultado

       return
       end
      
********************************************************************
********************************************************************
********************************************************************
c  define the subroutine that calculates the 
c Hadronic tensor for the coupling of the W+ to different bubbles
      subroutine HadronicTensorDifferentBubbles(mu,nu,q,cW)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension q(0:3),p(0:3),p1(0:3),xkpion(0:3),qminuskpion(0:3)
      dimension xkpi0d(2000),phid(2000),thetad(2000),xmodkpid(2000)
      dimension rd(2000),crrd(2000),crmodkpid(2000),crthetad(2000)
      dimension crphid(2000),crkpi0d(2000),dcosid(2000),crdcosid(2000)
      common /nucleus/rho0,a,th,na
c new
      common/matrizqslash/cMdeltaMatrix,cMidentity
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cProjector(4,4)
c new
c new
      dimension pdeltax(0:3),cdproj(0:3,0:3,4,4)
      dimension cpdeltaslashx(4,4)
      common/deltapole2/cdproj
!$OMP THREADPRIVATE(/deltapole2/)
c      common/deltapole2x/cdprojx(0:3,0:3,4,4)
c new




c constants of interest for the subroutine
      pi=3.141592653589793d0
      xM=0.94d0
      gA=1.26d0
      fpion=0.093d0
      hbarc=0.19733d0  !!! in GeV*fm
      xmodq=dsqrt(q(1)**2+q(2)**2+q(3)**2)

c if q(0) is negative, then the step function that goes in front of
c the integrals is zero and the Hadronic Tensor too.
      if (q(0).lt.0.d0)then
         cW=0.d0
         return
      endif

c precision in the integrals
      nr=1
      ncosid=2
      nmodkpi=1
      ntheta=1
      nphi=1
      nkpi0=1 

c Upper limits for the integrations
      rlimfermis=5.d0*a
      rlimGeV=rlimfermis/hbarc !!! in GeV^(-1)

      xcutoff=2.0d0  !!!cut-off for the modulus of kpion in GeV 

c CALL to subroutines to make the partitions
c      call DSG20R(0.d0,rlimGeV,nr,rd,nrr)
      call DSG20R(-1.d0,1.d0,ncosid,dcosid,ncosidcosid)
      call DSG20R(0.d0,xcutoff,nmodkpi,xmodkpid,nmodkpimodkpi)
c      call DSG20R(0.d0,pi,ntheta,thetad,nthetatheta)
c      call DSG20R(0.d0,2.d0*pi,nphi,phid,nphiphi)
c      call DSG20R(0.d0,q(0),nkpi0,xkpi0d,nkpi0kpi0)
      call STGS8(0.d0,rlimGeV,rd)
c      call STGS8(-1.d0,1.d0,dcosid)
c      call STGS8(0.d0,xcutoff,xmodkpid)
      call STGS8(0.d0,pi,thetad)
      call STGS8(0.d0,2.d0*pi,phid)
      call STGS8(0.d0,q(0),xkpi0d)

      do ir=1,8
         rGeV=rd(ir)
         rfermi=rGeV*hbarc
         rhofermi=density(rfermi) !!! density in fm^(-3)

         xkfermi=(3.d0*pi**2/2.d0*rhofermi)**(1.d0/3.d0) !!! momentum of Fermi in fm^(-1)
         xkFGeV=xkfermi*hbarc !!! Momentum of Fermi in GeV

c can calculate the averaged <pmom>
         pmom_averaged=dsqrt(3.d0/5.d0)*xkFGeV !!! in GeV

c And the energy that corresponds to this averaged <pmom>
c can calculate p and p1, assuming that the nucleons are at rest
         energyfermi=xkFGeV**2/(2.d0*xM)
c         p(0)=xM-energyfermi
c         p1(0)=xM-energyfermi
         p(0)=dsqrt(xM**2+pmom_averaged**2)-energyfermi
         p1(0)=p(0)

         do icosid=1,ncosidcosid !8
            dcoseno=dcosid(icosid)
            seno=dsqrt(1.d0-dcoseno**2)
            p(1)=pmom_averaged*(dcoseno*q(1)/xmodq)
            p(2)=pmom_averaged*(dcoseno*q(2)/xmodq+seno)
            p(3)=pmom_averaged*(dcoseno*q(3)/xmodq)
            do jj=1,3
               p1(jj)=p(jj)
            enddo
            
ccc  test change
      call SumFourVectors(p,q,pDeltax)
      call slash(pdeltax,cpdeltaslashx)
      call MatrixSum(cpdeltaslashx,cMdeltaMatrix,cpdMd)
           do i11=0,3
           do i22=0,3
           call DeltaProjector(i11,i22,pdeltax,
     &           cpdMd,cProjector)
              do ii=1,4
              do jj=1,4
           cdproj(i11,i22,ii,jj)=cProjector(ii,jj)
              enddo
              enddo
           enddo
           enddo
cccc  end test change



            do imodkpi=1,nmodkpimodkpi   !8
               xmodkpi=xmodkpid(imodkpi) !!! Modulus of kpi

               do itheta=1,8
                  theta=thetad(itheta)
                  xkpion(3)=xmodkpi*dcos(theta)

                  do iphi=1,8
                     phi=phid(iphi)
                     xkpion(1)=xmodkpi*dsin(theta)*dcos(phi)
                     xkpion(2)=xmodkpi*dsin(theta)*dsin(phi)

                     do ikpi0=1,8
                        xkpion(0)=xkpi0d(ikpi0)

                        call RestaFourVectors(q,xkpion,qminuskpion)
                        propagator1=Dpion(xkpion)
                        propagator2=Dpion(qminuskpion)
                        xImLind1=AImUr(qminuskpion,xkFGeV,xkFGeV)
                        xImLind2=AImUr(xkpion,xkFGeV,xkFGeV)
                        formfactor1=Fpi(xkpion)
                        formfactor2=Fpi(qminuskpion)

      if((abs(xImLind1).gt.1.d-20).and.(abs(xImLind2).gt.1.d-20))then
         call TensorADifferentBubbles(mu,nu,p,p1,q,xkpion,rhofermi,
     &        cAneutronpiplus,cAprotonpiplus)

      else
         cAneutronpiplus=0.d0
         cAprotonpiplus=0.d0
      endif
                        crkpi0d(ikpi0)=propagator1*propagator2*
     &               xImLind1*xImLind2*formfactor1**2*formfactor2**2*
     &               (cAneutronpiplus-cAprotonpiplus)
                     enddo
c                     call DRG20C(0.d0,q(0),nkpi0,crkpi0d,cres)
                     call REGS8C(0.d0,q(0),crkpi0d,cres)
                     crphid(iphi)=cres
                  enddo
c                  call DRG20C(0.d0,2.d0*pi,nphi,crphid,cres)
                  call REGS8C(0.d0,2.d0*pi,crphid,cres)
                  crthetad(itheta)=cres*dsin(theta)
               enddo
c               call DRG20C(0.d0,pi,ntheta,crthetad,cres)
               call REGS8C(0.d0,pi,crthetad,cres)
               crmodkpid(imodkpi)=cres*xmodkpi**2
            enddo
            call DRG20C(0.d0,xcutoff,nmodkpi,crmodkpid,cres)
c            call REGS8C(0.d0,xcutoff,crmodkpid,cres)
            crdcosid(icosid)=cres/2.d0
         enddo
         call DRG20C(-1.d0,1.d0,ncosid,crdcosid,cres)
c         call REGS8C(-1.d0,1.d0,crdcosid,cres)
         crrd(ir)=cres*rGeV**2
      enddo
c      call DRG20C(0.d0,rlimGeV,nr,crrd,cresultado)
      call REGS8C(0.d0,rlimGeV,crrd,cresultado)

      factor=-1.d0*(gA/(2.d0*fpion))**2*(4.d0*dsqrt(2.d0)/(64.d0*
     &     xM**4))*(4.d0*pi/(2.d0*pi)**5)

      cW=factor*cresultado

      return
      end

**********************************************************

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccc MEDIUM MODIFICATIONS!!!!!!!!CCCCCCCCCCCCCCCCCCCCCCCCCCC
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Tpion is the kinetic energy of the pions in the nucleon
c rest frame. It must be given in GeV.
c a1,a2,a3 are coefficients for the parabolic form for C(Tpion)
c as given in eq 4.5 NPA 468,631 (1987)

c For Cq,Ca2 and Ca3, xcoefficient is returned in MeV for a1,a2 and a3
c given in Table 2 of the above reference NPA 468, 631 (1987).

c For alpha,beta and gamma (which are dimensionless), xcoefficient is 
c returned dimensionless for a1,a2 and a3 given in Table 2 of the above
c reference
	real*8 function xcoefficient(a1,a2,a3,Tpion)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(xmu=0.13957d0)

c Definition of adimensional variable x=Tpion/mu
	x=Tpion/xmu

	xcoefficient=a1*x**2+a2*x+a3

	return
	end

*******************************************************************
*******************************************************************
*******************************************************************
c  define Heaviside(x), the step function, which
c is equal to 1 if x>=0 and equal to 0 if x<0

	real*8 function Heaviside(x)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)

	if (x.ge.0.d0)then
	   Heaviside=1.d0
	else
	   Heaviside=0.d0
	endif

	return
	end
*****************************************************************
*****************************************************************
*****************************************************************
c  program the real part of the Delta Selfenergy
c In the local density approximation, it depends on the density at the
c point 'r'. Rho is the density and must be given in fm^(-3).The result
c is given in GeV because xconstant=40 MeV=0.04 GeV is also given in
c GeV.
	real*8 function ReSelfenergyDelta(rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(xconstant=0.04d0,xrho0=0.17d0)

	ReSelfenergyDelta=xconstant*rho/xrho0*0.d0
	
	return
	end
******************************************************************
******************************************************************
******************************************************************
c  define the first integral angular of the appendix
c B in reference NPA 554, 554 (1993). It depends on q (the CM-momentum
c of the pion in which the Delta decays when it has an invariant mass 
c squared 's'). And it also depends on density via the Fermi momentum.
c 'q' is given in GeV and rho in fm^(-3)

	real*8 function xIntegral1(q,rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(hbarc=0.19733d0,pi=3.141592653589793d0)

	xkfermi=(3.d0*pi**2/2.d0*rho)**(1.d0/3.d0)
	xkFGeV=xkfermi*hbarc
c definition of qtilde, which is dimensionless
	qtilde=q/xkFGeV

	xIntegral1=1.d0+Heaviside(qtilde-1.d0)*(-1.d0/(5.d0*qtilde**2)
     &    )+Heaviside(1.d0-qtilde)*(qtilde-1.d0-qtilde**3/5.d0)

	return
	end
*********************************************************************
*********************************************************************
*********************************************************************
c  define the second angular integral of the 
c appendix B in reference NPA 554, 554 (1993). It depends on 'q'
c (the CM-momentum of the pion in which the Delta decays when it has
c an invariant mass squared 's'). And it also depends on density via 
c the Fermi momentum. 'q' is given in GeV and 'rho' in fm^(-3).
	real*8 function xIntegral2(q,rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(hbarc=0.19733d0,pi=3.141592653589793d0)

	xkfermi=(3.d0*pi**2/2.d0*rho)**(1.d0/3.d0)
	xkFGeV=xkfermi*hbarc

c definition of qtilde, which is dimensionless
	qtilde=q/xkFGeV

	xIntegral2=1.d0+Heaviside(qtilde-1.d0)*(-3.d0/5.d0*1.d0/
     &	qtilde**2-4.d0/21.d0*1.d0/qtilde**6+18.d0/35.d0*1.d0/
     &	qtilde**4)+Heaviside(1.d0-qtilde)*(-1.d0+33.d0/35.d0*
     &	qtilde-23.d0/105.d0*qtilde**3)

	return
	end
***********************************************************************
***********************************************************************
***********************************************************************
c  define the Pauli Delta Width, which depends on 
c 's' (the invariant mass squared of the Delta Resonance) and 'rho' 
c (the density in fm^(-3)) because the angular integrals of equation
c 18 in reference PRC 75, 055501 (2007) depend on rho. The Pauli
c Delta width is returned in GeV because the free delta width is in 
c GeV too.
	real*8 function PauliDeltaWidth(s,rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(xmpion=0.13957d0,xM=0.94d0)
	
        xlambCall=xlambdaCallen(s,xmpion**2,xM**2)
	if ((s.gt.0.d0).and.(xlambCall.gt.0.d0))then
	    q=dsqrt(xlambCall)/(2.d0*dsqrt(s))
	    PauliDeltaWidth=DeltaWidth(s)*(xIntegral1(q,rho)+
     &	    xIntegral2(q,rho))/2.d0
	else
	    PauliDeltaWidth=0.d0
	endif

	return
	end
***********************************************************************
***********************************************************************
***********************************************************************
c  parametrize the imaginary part of the Delta
c selfenergy. It depends on rho (density in fm^(-3)) explicitly and on
c 's' (the invariant mass squared of the delta resonance) implicitly
c via the kinetic energy of the pion (Tpion). The result is given in
c GeV because the constants are in GeV.
	real*8 function xImSelfenergyDelta(s,rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(xmpion=0.13957d0,xM=0.94d0,xrho0=0.17d0)
c Coefficients for eq 4.5 for an analytical interpolation of 
c ImSelfenergyDelta. Everything in GeV for Cq,Ca2 and Ca3.
c dimensionless coefficients for alpha,beta and gamma
	a1Cq=-0.00519d0
	a2Cq=0.01535d0
	a3Cq=0.00206d0

	a1Ca2=0.00106d0
	a2Ca2=-0.00664d0
	a3Ca2=0.02266d0

	a1Ca3=-0.01346d0
	a2Ca3=0.04617d0
	a3Ca3=-0.02034d0

	a1alpha=0.382d0
	a2alpha=-1.322d0
	a3alpha=1.466d0

	a1beta=-0.038d0
	a2beta=0.204d0
	a3beta=0.613d0

c Tpion is in GeV
	Tpion=(s-xmpion**2-xM**2)/(2.d0*xM)-xmpion

	if ((Tpion.gt.0.090d0).and.(Tpion.lt.0.315d0)) then
	    xCq=xcoefficient(a1Cq,a2Cq,a3Cq,Tpion)
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,Tpion)
	    xCa3=xcoefficient(a1Ca3,a2Ca3,a3Ca3,Tpion)
	    alpha=xcoefficient(a1alpha,a2alpha,a3alpha,Tpion)
	    beta=xcoefficient(a1beta,a2beta,a3beta,Tpion)
	    gamma=2.d0*beta

c Definition of the imaginary part of the Delta Selfenergy
	    xImSelfenergyDelta=-xCq*(rho/xrho0)**alpha-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**gamma

	 elseif((Tpion.ge.0.315d0)) then
	    xCq=xcoefficient(a1Cq,a2Cq,a3Cq,0.315d0)
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,0.315d0)
	    xCa3=xcoefficient(a1Ca3,a2Ca3,a3Ca3,0.315d0)
	    alpha=xcoefficient(a1alpha,a2alpha,a3alpha,0.315d0)
	    beta=xcoefficient(a1beta,a2beta,a3beta,0.315d0)
	    gamma=2.d0*beta

c Definition of the imaginary part of the Delta Selfenergy
	    xImSelfenergyDelta=-xCq*(rho/xrho0)**alpha-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**gamma

	 elseif((Tpion.gt.0.d0).and.(Tpion.le.0.090d0))then
	    xCq=Tpion/xmpion*(20.2d0-8.58d0*Tpion/xmpion+
     &      0.702d0*(Tpion/xmpion)**2)*1.d-3
	    alphaQ=1.d0+Tpion/xmpion*(-0.309d0-0.315d0*Tpion/xmpion+
     &      0.151d0*(Tpion/xmpion)**2)
	    xCa3=Tpion*1.d3/85.d0*3.7d0*1.d-3
	    alphaA3=1.d0+Tpion/xmpion*(0.984d0-0.512d0*Tpion/xmpion+
     &      0.1d0*(Tpion/xmpion)**2)
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,Tpion)
	    beta=xcoefficient(a1beta,a2beta,a3beta,Tpion)

ccccccc 
	    xImSelfenergyDelta=-xCq*(rho/xrho0)**alphaQ-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**alphaA3

	 else
            if(s.gt.xM**2)then
               xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,0.d0)
               beta=xcoefficient(a1beta,a2beta,a3beta,0.d0)
               xImSelfenergyDelta=-xCa2*(rho/xrho0)**beta
            else
               xImSelfenergyDelta=0.d0
            endif
         endif

	return
	end
***********************************************************************
***********************************************************************
***********************************************************************
c  parametrize the imaginary part of the Delta
c selfenergy. It depends on rho (density in fm^(-3)) explicitly and on
c 's' (the invariant mass squared of the delta resonance) implicitly
c via the kinetic energy of the pion (Tpion). The result is given in
c GeV because the constants are in GeV.
	real*8 function xImSelfenergyDelta2(s,rho)
	implicit complex*16 (c)
	implicit real*8 (a,b,d-h,o-z)
	parameter(xmpion=0.13957d0,xM=0.94d0,xrho0=0.17d0)
c Coefficients for eq 4.5 for an analytical interpolation of 
c ImSelfenergyDelta. Everything in GeV for Cq,Ca2 and Ca3.
c dimensionless coefficients for alpha,beta and gamma
	a1Cq=-0.00519d0
	a2Cq=0.01535d0
	a3Cq=0.00206d0

	a1Ca2=0.00106d0
	a2Ca2=-0.00664d0
	a3Ca2=0.02266d0

	a1Ca3=-0.01346d0
	a2Ca3=0.04617d0
	a3Ca3=-0.02034d0

	a1alpha=0.382d0
	a2alpha=-1.322d0
	a3alpha=1.466d0

	a1beta=-0.038d0
	a2beta=0.204d0
	a3beta=0.613d0

c Tpion is in GeV
	Tpion=(s-xmpion**2-xM**2)/(2.d0*xM)-xmpion

	if ((Tpion.gt.0.090d0).and.(Tpion.lt.0.315d0)) then
	    xCq=xcoefficient(a1Cq,a2Cq,a3Cq,Tpion)*0.d0
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,Tpion)
	    xCa3=xcoefficient(a1Ca3,a2Ca3,a3Ca3,Tpion)
	    alpha=xcoefficient(a1alpha,a2alpha,a3alpha,Tpion)
	    beta=xcoefficient(a1beta,a2beta,a3beta,Tpion)
	    gamma=2.d0*beta

c Definition of the imaginary part of the Delta Selfenergy
	    xImSelfenergyDelta2=-xCq*(rho/xrho0)**alpha-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**gamma

	 elseif((Tpion.ge.0.315d0)) then
	    xCq=xcoefficient(a1Cq,a2Cq,a3Cq,0.315d0)*0.d0
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,0.315d0)
	    xCa3=xcoefficient(a1Ca3,a2Ca3,a3Ca3,0.315d0)
	    alpha=xcoefficient(a1alpha,a2alpha,a3alpha,0.315d0)
	    beta=xcoefficient(a1beta,a2beta,a3beta,0.315d0)
	    gamma=2.d0*beta

c Definition of the imaginary part of the Delta Selfenergy
	    xImSelfenergyDelta2=-xCq*(rho/xrho0)**alpha-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**gamma

	 elseif((Tpion.gt.0.d0).and.(Tpion.le.0.090d0))then
	    xCq=Tpion/xmpion*(20.2d0-8.58d0*Tpion/xmpion+
     &      0.702d0*(Tpion/xmpion)**2)*1.d-3*0.d0
	    alphaQ=1.d0+Tpion/xmpion*(-0.309d0-0.315d0*Tpion/xmpion+
     &      0.151d0*(Tpion/xmpion)**2)
	    xCa3=Tpion*1.d3/85.d0*3.7d0*1.d-3
	    alphaA3=1.d0+Tpion/xmpion*(0.984d0-0.512d0*Tpion/xmpion+
     &      0.1d0*(Tpion/xmpion)**2)
	    xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,Tpion)
	    beta=xcoefficient(a1beta,a2beta,a3beta,Tpion)

ccccccc 
	    xImSelfenergyDelta2=-xCq*(rho/xrho0)**alphaQ-xCa2*(rho/
     &      xrho0)**beta-xCa3*(rho/xrho0)**alphaA3

	 else
            if(s.gt.xM**2)then
               xCa2=xcoefficient(a1Ca2,a2Ca2,a3Ca2,0.d0)
               beta=xcoefficient(a1beta,a2beta,a3beta,0.d0)
               xImSelfenergyDelta2=-xCa2*(rho/xrho0)**beta
            else
               xImSelfenergyDelta2=0.d0
            endif
         endif

	return
	end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccc PION SELF-ENERGY ccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c xkpion is the four-momentum of the virtual pion in GeV
c rho is the density in fm^(-3)

c cPolarization is the selfenergy of the pion given in GeV^2
      subroutine PolarizationPion(xkpion,rho,cPolarization)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3)

c Constants of interest
      gA=1.26d0
      fpion=0.093d0   !!! in GeV
      gprime=0.63d0
      xlambda=1.d0    !!! in GeV
      hbarc=0.19733d0 !!!factor of conversion

      call LorentzScalarProduct(xkpion,xkpion,xkpion2)
c      formfactor=xlambda**2/(xlambda**2-xkpion2) !!dimensionless

c the energy and momentum of the pion in fm^(-1)
      qzr=xkpion(0)/hbarc
      q=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2)/hbarc
      qgev=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2) !! in GeV
c      formfactor=xlambda**2/(xlambda**2+qgev**2)!!dimensionless
      formfactor=Fpi(xkpion)

      call ULIND2(qzr,q,rho,CUfun) !!where CUfun is in fm^(-2)

      CUfunGeV=CUfun*hbarc**2     !!! in GeV^2

      cPolarization=(gA/(2.d0*fpion))**2*formfactor**2*qgev**2*
     &     CUfunGeV/(1.d0-(gA/(2.d0*fpion))**2*gprime*
     &     CUfunGeV)

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Vt, where xk is the four-momentum in GeV. 
c Vt is returned in GeV^(-2)
      real*8 function Vt(xk)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xk(0:3)
c Constants of interest
      gA=1.26d0
      xcrho=2.d0
      fpion=0.093d0   !!! in GeV
      gprime=0.63d0
      xmrho=0.77d0 !!! in GeV

      xmodk=dsqrt(xk(1)**2+xk(2)**2+xk(3)**2) !!! in GeV
      call LorentzScalarProduct(xk,xk,xk2)
      ff=FormFactorRho(xk)

      Vt=(gA/(2.d0*fpion))**2*(xcrho*ff**2*xmodk**2/(xk2-xmrho**2)+
     &     gprime)
      return
      end
*********************************************************************
*********************************************************************
*********************************************************************
cimaginary part of the polarization of the pion in GeV²
c INPUT VARIABLES
c 'xkpion' is the pion four-momentum in GeV
c 'rho' is the density in fm^(-3)

      real*8 function xImPolarizationPion(xkpion,rho)
      implicit complex*16(c)
      implicit real*8 (a,b,d-h,o-z)
      dimension xkpion(0:3)

c Constants of interest
      gA=1.26d0
      fpion=0.093d0   !!! in GeV
      gprime=0.63d0
      hbarc=0.19733d0 !!!factor of conversion

      xmodkpion=dsqrt(xkpion(1)**2+xkpion(2)**2+xkpion(3)**2) !! in GeV
      xmodkpionfermis=xmodkpion/hbarc   !!! in fm^(-1)
      xkpi0=xkpion(0)/hbarc             !!! in fm^(-1)
      ff=Fpi(xkpion)

      call ULIND2(xkpi0,xmodkpionfermis,rho,CUfuncompleta)
      call ULIND3(xkpi0,xmodkpionfermis,rho,CUnuc)
c Conversion of the Lindhard functions to GeV^2
      CUfuncompletaGeV=CUfuncompleta*hbarc**2  !!! in GeV^2
      CUnucGeV=CUnuc*hbarc**2                  !!! in GeV^2
      

      xImPolarizationPion=(gA/(2.d0*fpion))**2*xmodkpion**2*ff**2*
     &     dimag(CUnucGeV)/(zabs(1.d0-(gA/(2.d0*fpion))**2*gprime*
     &     CUfuncompletaGeV))**2

      return
      end
**********************************************************************
**********************************************************************
**********************************************************************
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Lindhard Function
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        SUBROUTINE ULIND2(QZR,Q,rho,CUfun)
        implicit complex*16 (C)
        implicit real*8 (a,b,d-h,o-z)
        complex*16 lindhard,delta_lind,q_zero,q_zerol

        q_zero=dcmplx(qzr,0.0d0)
c        rho=2.d0*xkf**3/3.d0/3.14159265d0**2
	xkf=(3.d0/2.d0*3.14159265d0**2*rho)**(1.d0/3.d0)
c        xlam=1300.d0/197.33d0
c        ff=(xlam**2/(xlam**2+q**2))**2
        ff=1.d0

        q_zerol = q_zero
        cunuc=lindhard (q_zerol,q,rho,xkf)
        cudel=delta_lind (q_zero,q,rho,xkf)
        cufun=(cunuc+cudel)*ff

        return
        end
****************************************************************
****************************************************************
****************************************************************
        SUBROUTINE ULIND3(QZR,Q,rho,CUfun)
        implicit complex*16 (C)
        implicit real*8 (a,b,d-h,o-z)
        complex*16 lindhard,delta_lind,q_zero,q_zerol

        q_zero=dcmplx(qzr,0.0d0)
c        rho=2.d0*xkf**3/3.d0/3.14159265d0**2
	xkf=(3.d0/2.d0*3.14159265d0**2*rho)**(1.d0/3.d0)
c        xlam=1300.d0/197.33d0
c        ff=(xlam**2/(xlam**2+q**2))**2
        ff=1.d0

        q_zerol = q_zero
        cunuc=lindhard (q_zerol,q,rho,xkf)
c        cudel=delta_lind (q_zero,q,rho,xkf)
        cufun=(cunuc)*ff

        return
        end

        complex*16 function lindhard (q_zero,q_mod,rho,k_fermi)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   complex Lindhard function for symmetric nuclear matter:
c                    from Appendix of
c                    E.Oset et al Phys. Rept. 188:79, 1990
c                    formula A.1 
c
c   input variables: 
c     q_zero [fm^-1] : Energy
c     q_mod  [fm^-1] : Momentum
c     rho    [fm^-3]  : Nuclear density
c     k_fermi[fm^-1] : Fermi momentum
c
c   All variables are real*8
c
c   output variable: 
c     lindhard [fm^-2]
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        implicit none
        real*8 q_mod,rho,k_fermi,m,epsilon,pi,ep
        complex*16 q_zero,z,zp,i,pzeta,pzetap
c
c m = 939/197.3, epsilon is a small number
c
        data m,epsilon,i/4.7592d0,1.d-36,(0.d0,1.d0)/
        data pi/3.14159265358979323846264338d0/
        
        if((imag(q_zero))**2.lt.1.d-36)then
         q_zero = q_zero + i*sign(epsilon,dreal(q_zero))
        endif
c       the following line is useless.
        ep=epsilon
        
        z  = m /(q_mod*k_fermi) *( q_zero - q_mod**2/(2.d0*m))
        zp = m /(q_mod*k_fermi) *(-q_zero - q_mod**2/(2.d0*m))
        
c
c care with limit cases
c
        if(abs(z).gt.100.d0)then
            pzeta =  2.d0/(3.d0*z) +2.d0/(15.d0*z**3)
        else if(abs(z).lt.1.d-2)then
            pzeta =  2.d0*z -2.d0/3.d0*z**3 -i*pi/2.d0*(1.d0 - z**2) 
        else
            pzeta =  z +  (1.d0-z**2) * log((z+1.d0)/(z-1.d0))/2.d0
        endif
        
        if(abs(zp).gt.100.d0)then
            pzetap =  2.d0/(3.d0*zp) +2.d0/(15.d0*zp**3)
        else if(abs(zp).lt.1.d-2)then
        pzetap =  2.d0*zp -2.d0/3.d0*zp**3 -i*pi/2.d0*(1.d0 - zp**2) 
        else
        pzetap =  zp + (1.d0-zp**2) * log((zp+1.d0)/(zp-1.d0))/2.d0
        endif
                  
        lindhard=3.d0/2.d0*rho*m/(q_mod*k_fermi)*( pzeta + pzetap)
     
        return
        end
        
        
        complex*16 function delta_lind (q_zero,q_mod,rho,k_fermi)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   complex Lindhard function for symmetric nuclear matter:
c                    from Appendix of
c                    E.Oset et al Phys. Rept. 188:79, 1990
c                    formula A.4 
c
c   input variables: 
c     q_zero [fm^-1] : Energy
c     q_mod  [fm^-1] : Momentum
c     rho    [fm^-3]  : Nuclear density
c     k_fermi[fm^-1] : Fermi momentum
c
c   All variables are real*8
c
c   output variable: 
c     delta_lind [fm^-2]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            ATTENTION!!!
c Only works properly for real q_zero,
c if q_zero has an imaginary part calculates the L. function
c assuming Gamma= 0.
c Therefore this subroutine provides two different functions
c depending on whether q_zero is real or not!!!!!!!!!!!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        implicit none
        real*8 q_mod,rho,k_fermi,m,rq_zero,gamma,s,srot,wr
        real*8 fdel_f,mpi,gammap,sp,srotp,qcm,qcmp,md,pi
        complex*16 q_zero,z,zp,i,pzeta,pzetap
c
c m = 939/197.3,md = 1232./179.3, mpi = 139./197.3 
c
        data md,m,i/6.2433d0,4.7592d0,(0.d0,1.d0)/
        data mpi,pi/0.7045d0,3.14159265358979323846264338d0/
c
c  f*/f = 2.13 --> f*^2/4pi = .36
c  NOTE: fdel = fdel_f because f =1.0
c

        fdel_f = 2.13d0
        
        wr = md-m
        
        
        gamma = 0.d0 
        gammap = 0.d0
        if(imag(q_zero)**2.lt.1.d-36)then
          rq_zero = dreal(q_zero)
          s = m**2+rq_zero**2-q_mod**2+
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(s.gt.(m+mpi)**2)then
           srot = sqrt(s)
           qcm = sqrt(s**2+mpi**4+m**4-2.d0*(s*mpi**2+s*m**2+
     &    (mpi*m)**2)) /(2.d0*srot)
           gamma = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcm**3/srot*(m+sqrt(m**2+qcm**2))/mpi**2  
         endif          
          sp = m**2+rq_zero**2-q_mod**2-
     &        2.d0*rq_zero *sqrt(m**2+3.d0/5.d0*k_fermi**2)
          if(sp.gt.(m+mpi)**2)then
           srotp = sqrt(sp)
           qcmp=sqrt(sp**2+mpi**4+m**4-2.d0*(sp*mpi**2+sp*m**2+
     &                    (mpi*m)**2))/(2.d0*srotp)
           gammap = 1.d0/3.d0 * 1.d0/(4.d0*pi) * fdel_f**2*
     &             qcmp**3/srotp*(m+sqrt(m**2+qcmp**2))/mpi**2  
          endif          
        endif
              
        z=md/(q_mod*k_fermi)*(q_zero-q_mod**2/(2.d0*md)
     &                         -wr +i*gamma/2.d0)

        zp=md/(q_mod*k_fermi)*(-q_zero-q_mod**2/(2.d0*md)
     &                          -wr +i*gammap/2.d0)
c
c care with limit cases
c
        if(abs(z).gt.50.d0)then
            pzeta =  2.d0/(3.d0*z) +2.d0/(15.d0*z**3)
        else if(abs(z).lt.1.d-2)then
            pzeta =  2.d0*z -2.d0/3.d0*z**3 -i*pi/2.d0*(1.d0 - z**2) 
        else
            pzeta =  z +  (1.d0-z**2) * log((z+1.d0)/(z-1.d0))/2.d0
        endif
        
        if(abs(zp).gt.50.d0)then
            pzetap =  2.d0/(3.d0*zp) +2.d0/(15.d0*zp**3)
        else if(abs(zp).lt.1.d-2)then
            pzetap =  2.d0*zp -2.d0/3.d0*zp**3 -i*pi/2.d0*(1.d0 - zp**2) 
        else
            pzetap =  zp + (1.d0-zp**2) * log((zp+1.d0)/(zp-1.d0))/2.d0
        endif


       delta_lind = 2.d0/3.d0 * rho * md/(q_mod*k_fermi) * ( 
     &        pzeta +pzetap) * fdel_f **2    
     
        return
        end
******************************************************************
******************************************************************
******************************************************************
ccurrents for rho production. They depend now on two
c Lorentz indices and on the momenta of incoming nucleon (p) and
c W^+ (q) and outgoing rho momentum (xkrho).

c INPUT VARIABLES

c 'mu' is the index of the weak current
c 'ialpha' is the index of the rho meson
c 'p' is the momentum of the incoming nucleon
c 'q' is the momentum of the incoming W+
c 'xkrho' is the momentum of the outgoing rho meson
c

c OUTPUT VARIABLES

c cjota is the matrix for the current

      subroutine CurrentCrossedDeltaPoleRho(mu,ialpha,p,q,xkrho,
     &     cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),xkrho(0:3),pdelta(0:3),pprime(0:3)
      dimension qprime(0:3)
      dimension caux(4,4),cjota(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cpdMd(4,4),cGamma(4,4),cProjector(4,4)
      dimension cGammadaga(4,4),caux1(4,4),caux2(4,4),caux3(4,4)
      dimension caux4(4,4),ckrhoslash(4,4),caux5(4,4),caux6(4,4)
      dimension caux7(4,4),caux8(4,4),caux9(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
c new
      common/deltapolex/cdprojx(0:3,0:3,4,4)
!$OMP THREADPRIVATE(/deltapolex/)
      dimension cqprimeslash(4,4)

c     CONSTANTS OF INTEREST
      xcrho=2.d0
      fstar=2.14d0
      xmpion=0.13957d0
      xCabibbo=0.974d0
      xMdelta=1.232d0
      gprime=0.63d0
      xmrho=0.77d0
c first we must define pdelta=p-xkrho
      call RestaFourVectors(p,xkrho,pdelta)
      call LorentzScalarProduct(pdelta,pdelta,pdelta2)
      s=pdelta2
      call LorentzScalarProduct(xkrho,xkrho,xk2)
      xmodk=dsqrt(xkrho(1)**2+xkrho(2)**2+xkrho(3)**2)
      ffrho=FormFactorRho(xkrho)

cconstant that goes in front of the matrix
      constant=dsqrt(xcrho)*fstar/xmpion*1.d0/dsqrt(3.d0)*xCabibbo/
     &     (pdelta2-xMdelta**2+ci*xMdelta*DeltaWidth(s))*(1.d0+
     &     gprime*(xk2-xmrho**2)/(ffrho**2*xcrho*xmodk**2))

c must define pprime=pdelta+q and qprime=-q
      call SumFourVectors(pdelta,q,pprime)

      do i=0,3
         qprime(i)=-q(i)
      enddo

c have to set a matrix equal to zero
      do i=1,4
         do j=1,4
            caux(i,j)=0.d0
         enddo
      enddo
c We will need later ckrhoslash
      call slash(xkrho,ckrhoslash)
c
      call slash(qprime,cqprimeslash)
cmatrix
      do nu=0,3
         call Vertex(nu,mu,pprime,qprime,cqprimeslash,cGamma)
         call Dagger(cGamma,cGammadaga)
         call MatrixProductg0Mg0(cGammadaga,caux2)

      do i=1,4
         do j=1,4
            caux1(i,j)=0.d0
         enddo
      enddo


         do ibeta=0,3
            nuprime=nu
            ibetaprime=ibeta
                  
c new
                  do ii=1,4
                  do jj=1,4
                  cprojector(ii,jj)=cdprojx(nu,ibeta,ii,jj)
                  enddo
                  enddo
c new

c            call MatrixProduct(caux2,cProjector,caux3)
c in caux3 is stored the product Vertex*Projector, which must be
c multiplied by cga5
            call MatrixProductg5(cProjector,caux4)
            xconstant=g(nu,nuprime)*g(ibeta,ibetaprime)
            call ScalarProductMatrix(xconstant,caux4,caux5)
c must store kslash*g(beta,alpha)-k(beta)*gamma(alpha)
            xconstant1=g(ibeta,ialpha)
            call ScalarProductMatrix(xconstant1,ckrhoslash,caux6)
            xconstant2=xkrho(ibeta)
            call ScalarProductMatrix(xconstant2,cga(1,1,ialpha),
     &           caux7)
            call MatrixResta(caux6,caux7,caux8)
c must multiply caux5 by caux8
            call MatrixProduct(caux5,caux8,caux9)

            call MatrixSum(caux1,caux9,caux1)
         enddo
      call MatrixProduct(caux2,caux1,caux3)
      call MatrixSum(caux,caux3,caux)
      enddo
           
      call ComplexScalarProductMatrix(constant,caux,cjota)

      return
      end
******************************************************************
******************************************************************
******************************************************************
c currents for rho production. They depend now on two
c Lorentz indices and on the momenta of incoming nucleon (p) and
c W^+ (q) and outgoing rho momentum (xkrho).

c INPUT VARIABLES

c 'mu' is the index of the weak current
c 'ialpha' is the index of the rho meson
c 'p' is the momentum of the incoming nucleon
c 'q' is the momentum of the incoming W+
c 'xkrho' is the momentum of the outgoing rho meson
c

c OUTPUT VARIABLES

c cjota is the matrix for the current
      subroutine CurrentDeltaPoleRho(mu,ialpha,p,q,xkrho,cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),xkrho(0:3),pdelta(0:3)
      dimension caux(4,4),cjota(4,4)
      dimension cga(4,4,0:3),cga5(4,4),cIdent(4,4),g(0:3,0:3)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),ckrhoslash(4,4),cpdMd(4,4),caux1(4,4)
      dimension caux2(4,4),caux3(4,4),caux4(4,4),cProjector(4,4)
      dimension caux5(4,4),cGamma(4,4),caux6(4,4),caux7(4,4)
      dimension cGamma0(4,4),cGamma1(4,4),cGamma2(4,4),cGamma3(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity

c new
      dimension cdproj(0:3,0:3,4,4)
      common/deltapole/cdproj
!$OMP THREADPRIVATE(/deltapole/)
c new



c     CONSTANTS OF INTEREST
      xcrho=2.d0
      fstar=2.14d0
      xmpion=0.13957d0
      xCabibbo=0.974d0
      xMdelta=1.232d0
      gprime=0.63d0
      xmrho=0.77d0
c first the constant that multiplies the whole matrix
      call SumFourVectors(p,q,pdelta)
      call LorentzScalarProduct(pdelta,pdelta,pdelta2)
      s=pdelta2
      call LorentzScalarProduct(xkrho,xkrho,xk2)
      xmodk=dsqrt(xkrho(1)**2+xkrho(2)**2+xkrho(3)**2)
      ffrho=FormFactorRho(xkrho) 

      constant=(-1.d0)*dsqrt(xcrho)*fstar/xmpion*dsqrt(3.d0)*
     &     xCabibbo/(pdelta2-xMdelta**2+ci*xMdelta*DeltaWidth(s))*
     &     (1.d0+gprime*(xk2-xmrho**2)/(ffrho**2*xcrho*xmodk**2))

c matrix equal to zero
      do i=1,4
         do j=1,4
            caux(i,j)=0.d0
         enddo
      enddo

c we call slash to obtain ckrhoslash
      call slash(xkrho,ckrhoslash)
c (pdeltaslash+Mdelta) is stored in cpdMd
c      call slash(pdelta,cpdeltaslash)
c      call MatrixSum(cpdeltaslash,cMdeltaMatrix,cpdMd)


         call slash(q,cqslash)
         call Vertex(0,mu,p,q,cqslash,cGamma0)
         call Vertex(1,mu,p,q,cqslash,cGamma1)
         call Vertex(2,mu,p,q,cqslash,cGamma2)
         call Vertex(3,mu,p,q,cqslash,cGamma3)

c the matrix
      do nu=0,3
         xcons1=g(nu,ialpha)
         xcons2=xkrho(nu)
         call ScalarProductMatrix(xcons1,ckrhoslash,caux1)
         call ScalarProductMatrix(xcons2,cga(1,1,ialpha),caux2)
         call MatrixResta(caux1,caux2,caux3)
         call MatrixProductg5x(caux3,caux4)
         nuprime=nu


            do ibetaprime=0,3
               ibeta=ibetaprime
c                  call DeltaProjector(nuprime,ibetaprime,pdelta,cpdMd,
c     &                 cProjector)
c new
                  do ii=1,4
                  do jj=1,4
                  cprojector(ii,jj)=cdproj(nu,ibeta,ii,jj)
                   enddo
                  enddo
c new

                  call MatrixProduct(caux4,cProjector,caux5)

                  if(ibeta.eq.0)then
                  call MatrixProduct(caux5,cGamma0,caux6)
                  elseif(ibeta.eq.1)then
                  call MatrixProduct(caux5,cGamma1,caux6)
                  elseif(ibeta.eq.2)then
                  call MatrixProduct(caux5,cGamma2,caux6)
                  elseif(ibeta.eq.3)then
                  call MatrixProduct(caux5,cGamma3,caux6)
                  endif

                  xcons3=g(nu,nuprime)*g(ibetaprime,ibeta)
                  call ScalarProductMatrix(xcons3,caux6,caux7)

                  call MatrixSum(caux,caux7,caux)
         enddo
      enddo

      call ComplexScalarProductMatrix(constant,caux,cjota)

      return
      end
******************************************************************
******************************************************************
******************************************************************
c currents for rho production. They depend now on two
c Lorentz indices and on the momenta of incoming nucleon (p) and
c W^+ (q) and outgoing rho momentum (xkrho).

c INPUT VARIABLES

c 'mu' is the index of the weak current
c 'ialpha' is the index of the rho meson
c 'p' is the momentum of the incoming nucleon
c 'q' is the momentum of the incoming W+
c 'xkrho' is the momentum of the outgoing rho meson
c

c OUTPUT VARIABLES

c cjota is the matrix for the current

      subroutine CurrentNucleonPoleRho(mu,ialpha,p,q,xkrho,cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),xkrho(0:3),pplusq(0:3)
      dimension caux(4,4),cjota(4,4),cga(4,4,0:3),cga5(4,4),cIdent(4,4)
      dimension g(0:3,0:3),cqslash(4,4),cpdeltaslash(4,4)
      dimension cMdeltaMatrix(4,4),cMidentity(4,4),ckrhoslash(4,4)
      dimension caux1(4,4),caux2(4,4),caux3(4,4),caux4(4,4)
      dimension cpplusqslash(4,4),caux5(4,4),cVector(4,4),cAxial(4,4)
      dimension caux6(4,4),caux7(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
c CONSTANTS OF INTEREST
      xcrho=2.d0
      gA=1.26d0
      fpion=0.093d0
      xCabibbo=0.974d0
      xM=0.94d0
      gprime=0.63d0
      xmrho=0.77d0
      width=0.020d0 !!! in GeV ~ 20 MeV
c Firstly we need the four vector p+q and its scalar product
      call SumFourVectors(p,q,pplusq)
      call LorentzScalarProduct(pplusq,pplusq,pplusq2)
      call LorentzScalarProduct(xkrho,xkrho,xk2)
      xmodk=dsqrt(xkrho(1)**2+xkrho(2)**2+xkrho(3)**2)
      ffrho=FormFactorRho(xkrho)
cconstant that goes in front of the current
      constant=-ci*dsqrt(xcrho)*gA*xCabibbo/(dsqrt(2.d0)*fpion)*
     &     1.d0/(pplusq2-xM**2+ci*xM*width)*
     &     (1.d0+gprime*(xk2-xmrho**2)/(ffrho**2*xcrho*xmodk**2))

cmatrices. We will need ckrhoslash
      constant1=ci/2.d0
      call slash(xkrho,ckrhoslash)
      call MatrixProduct(ckrhoslash,cga(1,1,ialpha),caux1)
      call MatrixProduct(cga(1,1,ialpha),ckrhoslash,caux2)
      call MatrixResta(caux1,caux2,caux3)
      call ComplexScalarProductMatrix(constant1,caux3,caux4)
c In caux4 we have stored sigma[nu,alpha]*k[nu]

      call slash(pplusq,cpplusqslash)
      call MatrixSum(cpplusqslash,cMidentity,caux5)
c In caux5 we have stored (pslash+qslash+M)
      call slash(q,cqslash)
      call VectorN(mu,q,cqslash,cVector)
      call AxialN(mu,q,cqslash,cAxial)
      call MatrixResta(cVector,cAxial,caux6)
c In caux6 we have stored (V[mu]-A[mu])

c And finally, we must multiply the stored matrices
      call MatrixProduct(caux4,caux5,caux7)
      call MatrixProduct(caux7,caux6,caux)

      call ComplexScalarProductMatrix(constant,caux,cjota)

      return
      end
******************************************************************
******************************************************************
******************************************************************
c currents for rho production. They depend now on two
c Lorentz indices and on the momenta of incoming nucleon (p) and
c W^+ (q) and outgoing rho momentum (xkrho).

c INPUT VARIABLES

c 'mu' is the index of the weak current
c 'ialpha' is the index of the rho meson
c 'p' is the momentum of the incoming nucleon
c 'q' is the momentum of the incoming W+
c 'xkrho' is the momentum of the outgoing rho meson
c

c OUTPUT VARIABLES

c cjota is the matrix for the current
      subroutine CurrentCrossedNucleonPoleRho(mu,ialpha,p,q,xkrho,
     &     cjota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),xkrho(0:3),pminusk(0:3)
      dimension caux(4,4),cjota(4,4),cga(4,4,0:3),cga5(4,4),cIdent(4,4)
      dimension g(0:3,0:3),cqslash(4,4),cpdeltaslash(4,4)
      dimension cMdeltaMatrix(4,4),cMidentity(4,4),ckrhoslash(4,4)
      dimension caux1(4,4),caux2(4,4),caux3(4,4),caux4(4,4)
      dimension cpminuskslash(4,4),caux5(4,4),cVector(4,4),cAxial(4,4)
      dimension caux6(4,4),caux7(4,4)
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g
      common /complexi/ ci
      common /matrizqslash/cMdeltaMatrix,cMidentity
c CONSTANTS OF INTEREST
      xcrho=2.d0
      gA=1.26d0
      fpion=0.093d0
      xCabibbo=0.974d0
      xM=0.94d0
      gprime=0.63d0
      xmrho=0.77d0
c Firstly we need the four vector p+q and its scalar product
      call RestaFourVectors(p,xkrho,pminusk)
      call LorentzScalarProduct(pminusk,pminusk,pminusk2)
      call LorentzScalarProduct(xkrho,xkrho,xk2)
      xmodk=dsqrt(xkrho(1)**2+xkrho(2)**2+xkrho(3)**2)
      ffrho=FormFactorRho(xkrho)
cconstant that goes in front of the current
      constant=-ci*dsqrt(xcrho)*gA*xCabibbo/(dsqrt(2.d0)*fpion)*
     &     1.d0/(pminusk2-xM**2)*(1.d0+gprime*(xk2-xmrho**2)/(ffrho**2*
     &     xcrho*xmodk**2))

cmatrices. We will need ckrhoslash
      constant1=ci/2.d0
      call slash(xkrho,ckrhoslash)
      call MatrixProduct(ckrhoslash,cga(1,1,ialpha),caux1)
      call MatrixProduct(cga(1,1,ialpha),ckrhoslash,caux2)
      call MatrixResta(caux1,caux2,caux3)
      call ComplexScalarProductMatrix(constant1,caux3,caux4)
c In caux4 we have stored sigma[nu,alpha]*k[nu]

      call slash(pminusk,cpminuskslash)
      call MatrixSum(cpminuskslash,cMidentity,caux5)
c In caux5 we have stored (pslash-kslash+M)
      call slash(q,cqslash)
      call VectorN(mu,q,cqslash,cVector)
      call AxialN(mu,q,cqslash,cAxial)
      call MatrixResta(cVector,cAxial,caux6)
c In caux6 we have stored (V[mu]-A[mu])

c And finally, we must multiply the stored matrices
      call MatrixProduct(caux6,caux5,caux7)
      call MatrixProduct(caux7,caux4,caux)

      call ComplexScalarProductMatrix(constant,caux,cjota)

      return
      end
*************************************************************
*************************************************************
*************************************************************
      subroutine Sum4Currents(xc1,xc2,xc3,xc4,cj1,cj2,cj3,cj4,cJota)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cj1(4,4),cj2(4,4),cj3(4,4),cj4(4,4),cJota(4,4)

      do i=1,4
         do j=1,4
            cJota(i,j)=xc1*cj1(i,j)+xc2*cj2(i,j)+xc3*cj3(i,j)+
     &           xc4*cj4(i,j)
         enddo
      enddo

      return
      end
****************************************************************
****************************************************************
****************************************************************
c subroutine tensorArho, which
c calculates the tensor A(mu,nu) for a given momenta configuration
c (p,q,xkrho)
c
c INPUT VARIABLES
c 'mu' is the Lorentz index of the first weak vertex
c 'nu' is the Lorentz index of the second weak vertex
c 'p' is the momentum of the incoming nucleon
c 'q' is the momentum of the incoming W⁺
c 'xkrho' is the momentum of the outgoing rho meson 

c OUTPUT VARIABLES
c 'cAproton' is the tensor component (for the given Lorentz indices)
c for the reaction p + W⁺ --> p + rho⁺

c 'cAneutronrhoplus' is the tensor component for the reaction
c n + W⁺ --> n + rho⁺

c 'cAneutronrho0' is the tensor component for the reaction
c n + W⁺ --> p + rho⁰

      subroutine TensorArho(mu,nu,p,q,xkrho,cAproton,cAneutronrhoplus,
     &     cAneutronrho0)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension p(0:3),q(0:3),xkrho(0:3),pplusq(0:3),pprime(0:3)
      dimension cqslash(4,4),cpdeltaslash(4,4),cMdeltaMatrix(4,4)
      dimension cMidentity(4,4),cga(4,4,0:3),cga5(4,4),cIdent(4,4)
      dimension g(0:3,0:3),cpslash(4,4),cpprimeslash(4,4),cpM(4,4)
      dimension cpprimeM(4,4),cj1(4,4),cj2(4,4),cj3(4,4),cj4(4,4)
      dimension cj1bis(4,4),cj2bis(4,4),cj3bis(4,4),cj4bis(4,4)
      dimension cJotamuproton(4,4),cJotamuneutron(4,4)
      dimension cJotanuproton(4,4),cJotanuneutron(4,4),cauxpi0(4,4)
      dimension cJotamurho0(4,4),cJotanurho0(4,4),caux1(4,4),caux2(4,4)
      dimension cJtildenuproton(4,4),cJtildenuneutron(4,4),caux3(4,4)
      dimension cJtildenurho0(4,4),cj1tildenuproton(4,4)
      dimension cj1muneutron(4,4),cj1nuneutron(4,4),cj1tildenurho0(4,4)
      dimension cj1tildenuneutron(4,4),cj1murho0(4,4),cj3murho0(4,4)
      dimension cj3tildenuneutron(4,4),cj3muproton(4,4)
      dimension cj3tildenurho0(4,4)
      common /matrizqslash/cMdeltaMatrix,cMidentity
      common /matrices/ cga,cga5,cIdent
      common /metrictensor/ g

c  pprime = p + q - k
      call SumFourVectors(p,q,pplusq)
      call RestaFourVectors(pplusq,xkrho,pprime)
c pslash and pprimeslash
      call slash(p,cpslash)
      call slash(pprime,cpprimeslash)
c (pslash+M) and (pprimeslash+M)
      call MatrixSum(cpslash,cMidentity,cpM)
      call MatrixSum(cpprimeslash,cMidentity,cpprimeM)
c we have stored (pslash+M) in cpM
c and (pprimeslash+M) in cpprimeM. We cannot overwrite these two 
c matrices.

c isospin factors
c First the isospin factors for p --> p + rho⁺
      xCdeltaproton=1.d0
      xCCrossedDeltaproton=1.d0
      xCNucleonPoleproton=0.d0
      xCCrossedNucleonPoleproton=1.d0
c Secondly, the isospin factors for n --> n + rho⁺
      xCdeltaneutron=1.d0/3.d0
      xCCrossedDeltaneutron=3.d0
      xCNucleonPoleneutron=1.d0
      xCCrossedNucleonPoleneutron=0.d0

      csumaproton=0.d0
      csumaneutron=0.d0
      csumarho0=0.d0

      csumaprotonDeltaDelta=0.d0
      csumaneutronDeltaDelta=0.d0
      csumarho0DeltaDelta=0.d0

      csumaneutronNPtimesNP=0.d0
      csumarho0NPtimesNP=0.d0

      do ialpha=0,3
c         do ibeta=0,3
         ibeta=ialpha
c currents in the vertex 'mu'

            call CurrentDeltaPoleRho(mu,ialpha,p,q,xkrho,cj1)
            call CurrentCrossedDeltaPoleRho(mu,ialpha,p,q,xkrho,cj2)
            call CurrentNucleonPoleRho(mu,ialpha,p,q,xkrho,cj3)
            call CurrentCrossedNucleonPoleRho(mu,ialpha,p,q,xkrho,cj4)

c And we also call the currents in the vertex 'nu'
              if(nu.ne.mu)then
            call CurrentDeltaPoleRho(nu,ibeta,p,q,xkrho,cj1bis)
            call CurrentCrossedDeltaPoleRho(nu,ibeta,p,q,xkrho,cj2bis)
            call CurrentNucleonPoleRho(nu,ibeta,p,q,xkrho,cj3bis)
            call CurrentCrossedNucleonPoleRho(nu,ibeta,p,q,xkrho,cj4bis)
              else
              do iii=1,4
              do jjj=1,4
              cj1bis(iii,jjj)=cj1(iii,jjj)
              cj2bis(iii,jjj)=cj2(iii,jjj)
              cj3bis(iii,jjj)=cj3(iii,jjj)
              cj4bis(iii,jjj)=cj4(iii,jjj)
              enddo
              enddo
              endif     
c construct the sum (with the proper isospin factors) of the four
c currents to give cJotamuproton for the vertex 'mu'
            call Sum4Currents(xCdeltaproton,xCCrossedDeltaproton,
     &           xCNucleonPoleProton,xCCrossedNucleonPoleproton,
     &           cj1,cj2,cj3,cj4,cJotamuproton)
c And the same (with the isospin factors for neutron) to give cJotamuneutron
            call Sum4Currents(xCdeltaneutron,xCCrossedDeltaneutron,
     &           xCNucleonPoleNeutron,xCCrossedNucleonPoleneutron,
     &           cj1,cj2,cj3,cj4,cJotamuneutron)

c  construct the sum for the vertex 'nu'. The results will be stored
c in cJotanuproton and cJotanuneutron
            call Sum4Currents(xCdeltaproton,xCCrossedDeltaproton,
     &           xCNucleonPoleProton,xCCrossedNucleonPoleproton,
     &           cj1bis,cj2bis,cj3bis,cj4bis,cJotanuproton)
            call Sum4Currents(xCdeltaneutron,xCCrossedDeltaneutron,
     &           xCNucleonPoleNeutron,xCCrossedNucleonPoleneutron,
     &           cj1bis,cj2bis,cj3bis,cj4bis,cJotanuneutron)

c  construct the currents for the process n --> p + rho⁰
c Firstly, the vertex 'mu'
            xcons=-1.d0/dsqrt(2.d0)
            call MatrixResta(cJotamuproton,cJotamuneutron,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cJotamurho0)
c  the vertex 'nu'
            call MatrixResta(cJotanuproton,cJotanuneutron,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cJotanurho0)

c  matrices Jtilde=gamma0*J^\dagger*gamma0
            call Dagger(cJotanuproton,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cJtildenuproton)
            call MatrixProductg0Mg0(caux1,cJtildenuproton)


            call Dagger(cJotanuneutron,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cJtildenuneutron)
            call MatrixProductg0Mg0(caux1,cJtildenuneutron)

            call Dagger(cJotanurho0,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cJtildenurho0)
            call MatrixProductg0Mg0(caux1,cJtildenurho0)

c matrices jtilde for the substraction of DeltaPoleXDeltaPole
c and NucleonPoleXNucleonPole
c Firstly, the DeltaPoleXDeltaPole
            call Dagger(cj1bis,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cj1tildenuproton)
            call MatrixProductg0Mg0(caux1,cj1tildenuproton)

            call ScalarProductMatrix(xCdeltaneutron,cj1,cj1muneutron)
            call ScalarProductMatrix(xCdeltaneutron,cj1bis,cj1nuneutron)
            call Dagger(cj1nuneutron,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cj1tildenuneutron)
            call MatrixProductg0Mg0(caux1,cj1tildenuneutron)

            call MatrixResta(cj1,cj1muneutron,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cj1murho0)
            call MatrixResta(cj1tildenuproton,cj1tildenuneutron,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cj1tildenurho0)

c  NucleonPoleXNucleonPole
            call Dagger(cj3bis,caux1)
c            call MatrixProduct(cga(1,1,0),caux1,caux2)
c            call MatrixProduct(caux2,cga(1,1,0),cj3tildenuneutron)
            call MatrixProductg0Mg0(caux1,cj3tildenuneutron)

            call ScalarProductMatrix(xCNucleonPoleproton,cj3,
     &           cj3muproton)
            call MatrixResta(cj3muproton,cj3,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cj3murho0)
            call MatrixResta(cj3muproton,cj3tildenuneutron,cauxpi0)
            call ScalarProductMatrix(xcons,cauxpi0,cj3tildenurho0)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c matrices 
c J^{\mu\alpha}*(pslash+M)*Jtilde^{\nu\beta}*(pprimeslash+M)
c in order to calculate the traces
            call MatrixProduct(cj1,cpM,caux1)
            call MatrixProduct(caux1,cj1tildenuproton,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazaproton11)
           call MatrixPt(caux2,cpprimeM,ctrazaproton11)

            csumaprotonDeltaDelta=csumaprotonDeltaDelta+g(ialpha,ibeta)*
     &           ctrazaproton11

            call MatrixProduct(cJotamuproton,cpM,caux1)
            call MatrixProduct(caux1,cJtildenuproton,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazaproton)
           call MatrixPt(caux2,cpprimeM,ctrazaproton)

            csumaproton=csumaproton+g(ialpha,ibeta)*ctrazaproton
c There is no need to calculate the NPtimesNP for protons because the
c currents go multiplied by zero (see xCNucleonPoleproton). 
c But for 
c neutrons it is necessary because j3 goes multiplied by 1.d0
c (see xCNucleonPoleneutron)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Now for the process n + W⁺ --> n + rho⁺
            call MatrixProduct(cj1muneutron,cpM,caux1)
            call MatrixProduct(caux1,cj1tildenuneutron,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazaneutron11)
            call MatrixPt(caux2,cpprimeM,ctrazaneutron11)

            csumaneutronDeltaDelta=csumaneutronDeltaDelta+
     &           g(ialpha,ibeta)*ctrazaneutron11

            call MatrixProduct(cj3,cpM,caux1)
            call MatrixProduct(caux1,cj3tildenuneutron,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazaneutron33)
          call MatrixPt(caux2,cpprimeM,ctrazaneutron33)

            csumaneutronNPtimesNP=csumaneutronNPtimesNP+g(ialpha,ibeta)*
     &           ctrazaneutron33

            call MatrixProduct(cJotamuneutron,cpM,caux1)
            call MatrixProduct(caux1,cJtildenuneutron,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazaneutron)
           call MatrixPt(caux2,cpprimeM,ctrazaneutron)

            csumaneutron=csumaneutron+g(ialpha,ibeta)*ctrazaneutron
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Now for the process n + W⁺ --> p + rho⁰
            call MatrixProduct(cj1murho0,cpM,caux1)
            call MatrixProduct(caux1,cj1tildenurho0,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazarho011)
            call MatrixPt(caux2,cpprimeM,ctrazarho011)

            csumarho0DeltaDelta=csumarho0DeltaDelta+g(ialpha,ibeta)*
     &           ctrazarho011

            call MatrixProduct(cj3murho0,cpM,caux1)
            call MatrixProduct(caux1,cj3tildenurho0,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazarho033)
           call MatrixPt(caux2,cpprimeM,ctrazarho033)


            csumarho0NPtimesNP=csumarho0NPtimesNP+g(ialpha,ibeta)*
     &           ctrazarho033

            call MatrixProduct(cJotamurho0,cpM,caux1)
            call MatrixProduct(caux1,cJtildenurho0,caux2)
c            call MatrixProduct(caux2,cpprimeM,caux3)
c            call trace(caux3,ctrazarho0)
           call MatrixPt(caux2,cpprimeM,ctrazarho0)

            csumarho0=csumarho0+g(ialpha,ibeta)*ctrazarho0

c         enddo
      enddo
            
      cAproton=csumaproton-csumaprotonDeltaDelta
      cAneutronrhoplus=csumaneutron-csumaneutronDeltaDelta-
     &     csumaneutronNPtimesNP
      cAneutronrho0=csumarho0-csumarho0DeltaDelta-csumarho0NPtimesNP

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
c subrutinas integracion
c==============================================
c
C 8PUNTOS
        SUBROUTINE STGS8(XINF,XSUP,Z) 
	 implicit real*8 (a-h,o-z)  
         DIMENSION X(4),Z(2) 
         DATA X/-.9602898564d0,-.7966664774d0,-.5255324099d0,
     x-.1834346424d0/  
         DO 10 I=1,4  
         Z(I)=(XSUP-XINF)*X(I)/2.d0+XSUP/2.+XINF/2.d0 
         L=9-I 
         Z(L)=(XINF-XSUP)/2.d0*X(I)+(XSUP+XINF)/2.d0  
  10        CONTINUE   
         RETURN  
         END     
c 8 PUNTOS funciones complejas
	SUBROUTINE REGS8C(X1,X2,cF,cRES) 
	 implicit real*8 (a,b,d-h,o-z)
	 implicit complex*16 (c)
        DIMENSION cF(8),W(4) 
        DATA W/0.1012285362d0,0.2223810344d0,0.3137066458d0,
     x0.3626837833d0/
        cRES=0.d0
        DO 1 I=1,4
        L=9-I
 1	cRES=cRES+W(I)*(cF(I)+cF(L))
        cRES=cRES*(X2-X1)*.5d0
        RETURN 
        END


cccccccccccccccccccccccccccccccccccccccccccccccccc
c The matrices must be dimensioned in the main program
c This is for calculating the matricial product of matrices cA,gamma5
      subroutine MatrixProductg0Mg0(cA,cC)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)
      dimension cA(4,4),cC(4,4)

      do i=1,4
         do j=1,4
          cc(i,j)=cA(i,j)    
        enddo
      enddo
        cc(1,3)=-cc(1,3)
        cc(2,3)=-cc(2,3)
        cc(3,1)=-cc(3,1)
        cc(3,2)=-cc(3,2)
        cc(1,4)=-cc(1,4)
        cc(4,1)=-cc(4,1)
        cc(4,2)=-cc(4,2)
        cc(2,4)=-cc(2,4)
      return
      end
cccccccccccccccccccccccccccccccccccccc
****************************************************************
****************************************************************
c The matrices must be dimensioned in the main program
c This is for calculating the matricial product of two matrices cA,cB
      subroutine MatrixPT(cA,cB,ctraza)
      implicit complex*16 (c)
      implicit real*8 (a,b,d-h,o-z)

      dimension cA(4,4),cB(4,4)

            ctraza=0.d0
      do i=1,4
            do k=1,4
               ctraza=ctraza+cA(i,k)*cB(k,i)
            enddo
      enddo
      return
      end




