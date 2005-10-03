c______________________________________________________________________
c \lib      [dummy] libneugen3.a
c
c \brief    Dummy NeuGEN's for getting GENIE's NeuGEN Facade compiling
c           even in systems with no working NeuGEN
c
c \author   Costas Andreopoulos (RAL)  <C.V.Andreopoulos@rl.ac.uk>
c
c \created  October 03, 2005
c______________________________________________________________________


c______________________________________________________________________
        SUBROUTINE print_configuration
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE initialize_configuration(s,i,j,b)
        IMPLICIT NONE
        CHARACTER *(10) s
        INTEGER  i,j
        LOGICAL b
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE set_parameters(s,i,f)
        IMPLICIT NONE
        CHARACTER *(10) s
        INTEGER i
        REAL f
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE set_pdfset(i,j,f)
        IMPLICIT NONE
        INTEGER i,j
        REAL f
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE makestate(i1,i2,i3,i4,i5,i6,i7,i8,i9)
        IMPLICIT NONE
        INTEGER i1,i2,i3,i4,i5,i6,i7,i8,i9
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE writestate(i)
        IMPLICIT NONE
        INTEGER i
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE set_default_parameters
        IMPLICIT NONE
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE set_kvcuts(i,f1,f2)
        IMPLICIT NONE
        INTEGER i
        REAL f1,f2
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE set_masks(i,b1,b2,b3)
        IMPLICIT NONE
        INTEGER i
        LOGICAL b1,b2,b3
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE sig_value(f1,i1,i2,f2)
        IMPLICIT NONE
        INTEGER i1,i2
        REAL f1,f2
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE dsig_value(f1,i1,f2,i2,i3,f3)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        REAL f1,f2,f3
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE nu_structurefunctions(
     +                    i1,i2,i3,f1,f2,i4,f3,f4,f5,f6,f7,f8)
        IMPLICIT NONE
        INTEGER i1,i2,i3,i4
        REAL f1,f2,f3,f4,f5,f6,f7,f8
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE e_structurefunctions(i1,i2,f1,f2,i3,f3,f4)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        REAL f1,f2,f3,f4
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE ddsig_e_value(f1,i1,f2,i2,f3,i3,i4,f4)
        IMPLICIT NONE
        INTEGER i1,i2,i3,i4
        REAL f1,f2,f3,f4
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE gen_control(s,i,b)
        IMPLICIT NONE
        CHARACTER *(10) s
        INTEGER i
        LOGICAL b
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE heplst(i)
        IMPLICIT NONE
        INTEGER i
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE generate_nu_event(i1,f,i2,i3)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        REAL f
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_n_stdhep_entries(i)
        IMPLICIT NONE
        INTEGER i
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_stdhep_particle_info(i1,i2,i3)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_stdhep_particle_mothers(i1,i2,i3)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_stdhep_particle_daughters(i1,i2,i3)
        IMPLICIT NONE
        INTEGER i1,i2,i3
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_stdhep_particle_p4(i,d1,d2,d3,d4)
        IMPLICIT NONE
        INTEGER i
        DOUBLE PRECISION d1,d2,d3,d4
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE get_stdhep_particle_v4(i,d1,d2,d3,d4)
        IMPLICIT NONE
        INTEGER i
        DOUBLE PRECISION d1,d2,d3,d4
        CALL dummy_neugen_warning
        RETURN
        END
c______________________________________________________________________
        SUBROUTINE dummy_neugen_warning

        WRITE(*,*) '************************************************'
        WRITE(*,*) '*             <<<< WARNING >>>>                *'
        WRITE(*,*) '* The NeuGenWrapper was asked to drive NeuGEN  *'
        WRITE(*,*) '* but the dummy NeuGEN library is linked in!!  *'
        WRITE(*,*) '*  ** Your installation must be messed-up. **  *'
        WRITE(*,*) '* Are you sure that you installed NeuGEN and   *'
        WRITE(*,*) '* that you enabled at your GENIE installation? *'
        WRITE(*,*) '************************************************'

        RETURN
        END
c______________________________________________________________________

