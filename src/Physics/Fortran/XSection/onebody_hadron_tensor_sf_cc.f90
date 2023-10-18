  module onebody_hadron_tensor_sf
    implicit none
  end module onebody_hadron_tensor_sf    

! Computes the hadronic response tensor given form factors
! initial and final hadron momentum, energy transfer
! Output is complex 4x4 resp array
! This function expects that q is along Z!!!
  subroutine compute_hadron_tensor(xmn_in, w, wt, xk_x, xk_y, xk_z, &
        &       q_x, q_y, q_z, ff1v, ff2v, ffa, ffp, resp) &
        &       bind(C,name="compute_hadron_tensor_SF_CC")
    use onebody_hadron_tensor_sf
    use onebody_currents_sf
    implicit none
    integer*4 :: i 
    real*8 :: w, wt, xk, xp, q_x, q_y, q_z, xmn_in
    real*8 :: ek,epf
    real*8 :: xk_x,xk_y,xk_z
    real*8 :: p_4(4),pp_4(4),qt_4(4)
    real*8 :: ff1v, ff2v, ffa, ffp
    complex*16 :: resp(4,4)
    logical :: conserve_current

      conserve_current = .TRUE.

      !Set up all dirac matrices
      call dirac_matrices_in(xmn_in, conserve_current)

      ! Define initial and final nucleon 4 momentum
      ! Remember that q points along z
      xk = sqrt(xk_x**2 + xk_y**2 + xk_z**2)
      xp = sqrt((xk_x+q_x)**2 + (xk_y+q_y)**2 + (xk_z + q_z)**2)

      ek = sqrt(xmn**2 + xk**2)
      epf = sqrt(xmn**2 + xp**2)

      qt_4(1) = wt
      qt_4(2) = q_x
      qt_4(3) = q_y
      qt_4(4) = q_z

      p_4(1)=ek
      p_4(2)=xk_x
      p_4(3)=xk_y
      p_4(4)=xk_z

      pp_4(1) = epf
      pp_4(2) = p_4(2) + qt_4(2)
      pp_4(3) = p_4(3) + qt_4(3)
      pp_4(4) = p_4(4) + qt_4(4)

      call current_init(p_4,pp_4,qt_4,w)
      call define_spinors()
      call sigccc2(resp,ff1v,ff2v,ffa,ffp)
      return
    end subroutine compute_hadron_tensor
    
    subroutine sigccc(resp,ff1v,ff2v,ffa,ffp) bind(C,name="sigccc2")
      use onebody_currents_sf
      implicit none
      real*8 :: ff1v,ff2v,ffa,ffp
      complex*16 :: resp(4,4)
      complex*16 :: inter(4,4)
      integer*4 :: i,j 
      
      call det_Ja(ff1v,ff2v,ffa,ffp)
      call det_res1b(resp)

      !...Spin average gives a factor of 1/2 to the nuclear response tensor 
      resp = resp*0.5d0
      
      call shift2(resp)

      inter = TRANSPOSE(resp)
      resp = inter

      return
    end subroutine sigccc

    !This module shifts the response tensor so that it is in (x,y,z,t), used by
    !TLorentzVector and GENIE convention instead of (t,x,y,z)
    subroutine shift(resp) bind(C,name="shift2")
      use onebody_currents_sf
      implicit none
      complex*16 :: resp(4,4)

      resp = CSHIFT(resp, 1, DIM=1)
      resp = CSHIFT(resp, 1, DIM=2)

      return
    end subroutine shift

