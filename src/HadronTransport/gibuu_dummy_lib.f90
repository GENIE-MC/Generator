!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! GENIE/GiBUU INTERFACE LIB
! This builds a dummy version of the library, to be used at GENIE
! builds where GiBUU is disabled.
!
! Tina Leitner (Giessen Univ.) & Costas Andreopoulos (Rutherford Lab)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE GiBUU_GENIE
  IMPLICIT NONE
CONTAINS

!  
! Set the nucleus pdg code (GENIE convention 1aaazzz000)
! 
  SUBROUTINE SetNucleus(pdg)
	INTEGER, INTENT(IN) :: pdg
        CALL DummyLibWarning()
  END SUBROUTINE SetNucleus

!  
! Add a GENIE hadron shower particle to the GiBUU particle array 
!
  SUBROUTINE AddHadron(in_had, pdg, px, py, pz, E)
	INTEGER, INTENT(IN) :: in_had, pdg
	REAL,    INTENT(IN) :: px, py, pz, E
        CALL DummyLibWarning()
  END SUBROUTINE AddHadron

!
! The GiBUU hadron transport model run on the specified GENIE event
! 
  SUBROUTINE HadronTransportForGENIE()
        CALL DummyLibWarning()
  END SUBROUTINE HadronTransportForGENIE

!  
! Get the number of daughter (final state) hadrons generated from the input hadron in_had
! 
  FUNCTION NumOfFinalStateHadrons(in_had)
	INTEGER, INTENT(IN) :: in_had 
	INTEGER :: NumOfFinalStateHadrons
	NumOfFinalStateHadrons = 0
        CALL DummyLibWarning()
  END FUNCTION NumOfFinalStateHadrons

!  
! Get the PDG code for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronPDG(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	INTEGER :: FinalStateHadronPDG
	FinalStateHadronPDG = 0
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronPDG

!  
! Get the x momentum component for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronPX(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronPX
	FinalStateHadronPX = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronPX

!  
! Get the y momentum component for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronPY(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronPY
	FinalStateHadronPY = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronPY

!  
! Get the z momentum component for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronPZ(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronPZ
	FinalStateHadronPZ = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronPZ

!  
! Get the energy for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronE(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronE
	FinalStateHadronE = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronE

!  
! Get the x-position for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronX(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronX
	FinalStateHadronX = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronX

!  
! Get the y-position for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronY(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronY
	FinalStateHadronY = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronY

!  
! Get the z-position for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronZ(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had 
	REAL :: FinalStateHadronZ
	FinalStateHadronZ = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronZ

!  
! Get the temporal 4-position component for the out_had daughter of the in_had input hadron
! 
  FUNCTION FinalStateHadronT(in_had,out_had)
	INTEGER, INTENT(IN) :: in_had,out_had
	REAL :: FinalStateHadronT
	FinalStateHadronT = 0.
        CALL DummyLibWarning()
  END FUNCTION FinalStateHadronT

!  
! Print a warning to let the user know that he is not using the actual GiBUU code but
! just a dummy library
! 
  SUBROUTINE DummyLibWarning()
        PRINT *, '****************************************************'
	PRINT *, '*** WARNING:                                     ***'
        PRINT *, '*** YOU ARE USING A DUMMY GIBUU LIBRARY!         ***'
	PRINT *, '*** Please obtain the actual GiBUU code from:    ***'
	PRINT *, '*** http://tp8.physik.uni-giessen.de:8080/GiBUU/ ***'
        PRINT *, '****************************************************'
        CALL abort()
  END SUBROUTINE DummyLibWarning

END MODULE GiBUU_GENIE
  
