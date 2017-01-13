IMPLICIT NONE
REAL     WWST, WWST2, G, VKRM, EXCM, BETA, BTG, ELFC, WOLD, WNEW
REAL     PIHF, EPSU2, EPSUST, EPSIT, EPSA, ZTMIN, ZTMAX, HPBL,     &
& SQVISC
REAL     RIC, RRIC, FHNEU, RFC, RFAC, ZZ, PSLMU, PSLMS, PSLHU,     &
& PSLHS
REAL     XX, PSPMU, YY, PSPMS, PSPHU, PSPHS, ZLM, Z0, THZ0, THLM
REAL     SFCSPD, CZIL, AKMS, AKHS, ZILFC, ZU, ZT, RDZ, CXCH
REAL     DTHV, DU2, BTGH, WSTAR2, USTAR, ZSLU, ZSLT, RLOGU, RLOGT
REAL     RLMO, ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4
!CC   ......REAL ZTFC

REAL     XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN,  &
&         RLMA

INTEGER  ITRMX, ILECH, ITR
PARAMETER                                                         &
&        (WWST = 1.2,WWST2 = WWST * WWST,G = 9.8,VKRM = 0.40,      &
&         EXCM = 0.001                                             &
&        ,BETA = 1./270.,BTG = BETA * G,ELFC = VKRM * BTG          &
&                  ,WOLD =.15,WNEW = 1. - WOLD,ITRMX = 05,         &
&                   PIHF = 3.14159265/2.)
PARAMETER                                                         &
&         (EPSU2 = 1.E-4,EPSUST = 0.07,EPSIT = 1.E-4,EPSA = 1.E-8  &
&         ,ZTMIN = -5.,ZTMAX = 1.,HPBL = 1000.0                    &
&          ,SQVISC = 258.2)
PARAMETER                                                         &
&       (RIC = 0.183,RRIC = 1.0/ RIC,FHNEU = 0.8,RFC = 0.191       &
&        ,RFAC = RIC / (FHNEU * RFC * RFC))

! ----------------------------------------------------------------------
! NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
! ----------------------------------------------------------------------
! LECH'S SURFACE FUNCTIONS
! ----------------------------------------------------------------------
PSLMU (ZZ)= -0.96* log (1.0-4.5* ZZ)
PSLMS (ZZ)= ZZ * RRIC -2.076* (1. -1./ (ZZ +1.))
PSLHU (ZZ)= -0.96* log (1.0-4.5* ZZ)

! ----------------------------------------------------------------------
! PAULSON'S SURFACE FUNCTIONS
! ----------------------------------------------------------------------
PSLHS (ZZ)= ZZ * RFAC -2.076* (1. -1./ (ZZ +1.))
PSPMU (XX)= -2.* log ( (XX +1.)*0.5) - log ( (XX * XX +1.)*0.5)   &
&        +2.* ATAN (XX)                                            &
&- PIHF
PSPMS (YY)= 5.* YY
PSPHU (XX)= -2.* log ( (XX * XX +1.)*0.5)

! ----------------------------------------------------------------------
! THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
! OVER SOLID SURFACE (LAND, SEA-ICE).
! ----------------------------------------------------------------------
PSPHS (YY)= 5.* YY

! ----------------------------------------------------------------------
!     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
!     C......ZTFC=0.1
!     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
! ----------------------------------------------------------------------
ILECH = 0

! ----------------------------------------------------------------------
ZILFC = - CZIL * VKRM * SQVISC
!     C.......ZT=Z0*ZTFC
ZU = Z0
RDZ = 1./ ZLM
CXCH = EXCM * RDZ
DTHV = THLM - THZ0

! ----------------------------------------------------------------------
! BELJARS CORRECTION OF USTAR
! ----------------------------------------------------------------------
DU2 = MAX (SFCSPD * SFCSPD,EPSU2)
!cc   If statements to avoid TANGENT LINEAR problems near zero
BTGH = BTG * HPBL
IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
   WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
ELSE
   WSTAR2 = 0.0
END IF

! ----------------------------------------------------------------------
! ZILITINKEVITCH APPROACH FOR ZT
! ----------------------------------------------------------------------
USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)

! ----------------------------------------------------------------------
ZT = EXP (ZILFC * SQRT (USTAR * Z0))* Z0
ZSLU = ZLM + ZU
!     PRINT*,'ZSLT=',ZSLT
!     PRINT*,'ZLM=',ZLM
!     PRINT*,'ZT=',ZT

ZSLT = ZLM + ZT
RLOGU = log (ZSLU / ZU)

RLOGT = log (ZSLT / ZT)
!     PRINT*,'RLMO=',RLMO
!     PRINT*,'ELFC=',ELFC
!     PRINT*,'AKHS=',AKHS
!     PRINT*,'DTHV=',DTHV
!     PRINT*,'USTAR=',USTAR

RLMO = ELFC * AKHS * DTHV / USTAR **3
! ----------------------------------------------------------------------
! 1./MONIN-OBUKKHOV LENGTH-SCALE
! ----------------------------------------------------------------------
DO ITR = 1,ITRMX
   ZETALT = MAX (ZSLT * RLMO,ZTMIN)
   RLMO = ZETALT / ZSLT
   ZETALU = ZSLU * RLMO
   ZETAU = ZU * RLMO

   ZETAT = ZT * RLMO
   IF (ILECH .eq. 0) THEN
      IF (RLMO .lt. 0.)THEN
         XLU4 = 1. -16.* ZETALU
         XLT4 = 1. -16.* ZETALT
         XU4 = 1. -16.* ZETAU

         XT4 = 1. -16.* ZETAT
         XLU = SQRT (SQRT (XLU4))
         XLT = SQRT (SQRT (XLT4))
         XU = SQRT (SQRT (XU4))

         XT = SQRT (SQRT (XT4))
!     PRINT*,'-----------1------------'
!     PRINT*,'PSMZ=',PSMZ
!     PRINT*,'PSPMU(ZETAU)=',PSPMU(ZETAU)
!     PRINT*,'XU=',XU
!     PRINT*,'------------------------'
         PSMZ = PSPMU (XU)
         SIMM = PSPMU (XLU) - PSMZ + RLOGU
         PSHZ = PSPHU (XT)
         SIMH = PSPHU (XLT) - PSHZ + RLOGT
      ELSE
         ZETALU = MIN (ZETALU,ZTMAX)
         ZETALT = MIN (ZETALT,ZTMAX)
!     PRINT*,'-----------2------------'
!     PRINT*,'PSMZ=',PSMZ
!     PRINT*,'PSPMS(ZETAU)=',PSPMS(ZETAU)
!     PRINT*,'ZETAU=',ZETAU
!     PRINT*,'------------------------'
         PSMZ = PSPMS (ZETAU)
         SIMM = PSPMS (ZETALU) - PSMZ + RLOGU
         PSHZ = PSPHS (ZETAT)
         SIMH = PSPHS (ZETALT) - PSHZ + RLOGT
      END IF
! ----------------------------------------------------------------------
! LECH'S FUNCTIONS
! ----------------------------------------------------------------------
   ELSE
      IF (RLMO .lt. 0.)THEN
!     PRINT*,'-----------3------------'
!     PRINT*,'PSMZ=',PSMZ
!     PRINT*,'PSLMU(ZETAU)=',PSLMU(ZETAU)
!     PRINT*,'ZETAU=',ZETAU
!     PRINT*,'------------------------'
         PSMZ = PSLMU (ZETAU)
         SIMM = PSLMU (ZETALU) - PSMZ + RLOGU
         PSHZ = PSLHU (ZETAT)
         SIMH = PSLHU (ZETALT) - PSHZ + RLOGT
      ELSE
         ZETALU = MIN (ZETALU,ZTMAX)

         ZETALT = MIN (ZETALT,ZTMAX)
!     PRINT*,'-----------4------------'
!     PRINT*,'PSMZ=',PSMZ
!     PRINT*,'PSLMS(ZETAU)=',PSLMS(ZETAU)
!     PRINT*,'ZETAU=',ZETAU
!     PRINT*,'------------------------'
         PSMZ = PSLMS (ZETAU)
         SIMM = PSLMS (ZETALU) - PSMZ + RLOGU
         PSHZ = PSLHS (ZETAT)
         SIMH = PSLHS (ZETALT) - PSHZ + RLOGT
      END IF
! ----------------------------------------------------------------------
! BELJAARS CORRECTION FOR USTAR
! ----------------------------------------------------------------------
   END IF

! ----------------------------------------------------------------------
! ZILITINKEVITCH FIX FOR ZT
! ----------------------------------------------------------------------
   USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)

   ZT = EXP (ZILFC * SQRT (USTAR * Z0))* Z0
   ZSLT = ZLM + ZT
!-----------------------------------------------------------------------
   RLOGT = log (ZSLT / ZT)
   USTARK = USTAR * VKRM
   AKMS = MAX (USTARK / SIMM,CXCH)
!-----------------------------------------------------------------------
! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
!-----------------------------------------------------------------------
   AKHS = MAX (USTARK / SIMH,CXCH)
   IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
      WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
   ELSE
      WSTAR2 = 0.0
   END IF
!-----------------------------------------------------------------------
   RLMN = ELFC * AKHS * DTHV / USTAR **3
!-----------------------------------------------------------------------
!     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
!-----------------------------------------------------------------------
   RLMA = RLMO * WOLD+ RLMN * WNEW
!-----------------------------------------------------------------------
   RLMO = RLMA
!     PRINT*,'----------------------------'
!     PRINT*,'SFCDIF OUTPUT !  ! ! ! ! ! ! ! !  !   !    !'

!     PRINT*,'ZLM=',ZLM
!     PRINT*,'Z0=',Z0
!     PRINT*,'THZ0=',THZ0
!     PRINT*,'THLM=',THLM
!     PRINT*,'SFCSPD=',SFCSPD
!     PRINT*,'CZIL=',CZIL
!     PRINT*,'AKMS=',AKMS
!     PRINT*,'AKHS=',AKHS
!     PRINT*,'----------------------------'

END DO
! ----------------------------------------------------------------------
END SUBROUTINE SFCDIF_off
! ----------------------------------------------------------------------
