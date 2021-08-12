!---------------------------------------------------------------------------------
!
! This program is part of ECLIPSE.
!
!
! ECLIPSE: a fast Quadratic Maximum Likelihood estimator for CMB
!          intensity and polarization power spectra
!
! Copyright (C) 2021      Juan Daniel Bilbao Ahedo
!
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
! For more information about ECLIPSE see https://github.com/CosmoTool/ECLIPSE
!
!---------------------------------------------------------------------------------




!---------------------------------------------------------------------------------
!
!     This program makes use of the HEALPix subroutines (or modified versions):
!     
!     getsize_fits
!     printerror
!     read_bintab
!     
!     HEALPix Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!     
!     For more information about HEALPix see http://healpix.sourceforge.net
!
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
!
!     This program makes use of FITSIO subroutines.
!     
!     FITSIO Copyright (Unpublished--all rights reserved under the copyright laws of
!     the United States), U.S. Government as represented by the Administrator
!     of the National Aeronautics and Space Administration.  No copyright is
!     claimed in the United States under Title 17, U.S. Code.
!     
!     
!     For more information about FITSIO see https://heasarc.gsfc.nasa.gov/fitsio/
!
!---------------------------------------------------------------------------------



MODULE DatosProblema

  IMPLICIT NONE

  INTEGER         :: NSide, Lmax, Dlmax
  INTEGER         :: TipoAnalisis, QuitarSesgo
  REAL(kind=8)    :: RuidoTT, RuidoQQ
  REAL(kind=8)    :: FWHM_Beam
  INTEGER         :: PixelWindowEnMapas
  INTEGER         :: KindOfGrouping
  INTEGER         :: KindOfBinCenter
  INTEGER         :: ComputeFisherMatrix
  INTEGER         :: Binned

  CHARACTER(len=100) :: Lugar
  CHARACTER(len=100) :: ArchivoProblema
  CHARACTER(len=100) :: ArchivoFiducial
  CHARACTER(len=100) :: FileMaskIntensity
  CHARACTER(len=100) :: FileMaskPolarization
  CHARACTER(len=100) :: FileMaps
  CHARACTER(len=100) :: FileMapsCross
  CHARACTER(len=100) :: FileBinLimits
  CHARACTER(len=100) :: DirHealpixData

  INTEGER :: NMapas
  INTEGER :: NumMaskMapIntensity          
  INTEGER :: NumMaskMapPolarization       
  INTEGER :: TipoFileMaps                 !TQU file or only QU file

  INTEGER :: ControlInversaMC             

  INTEGER :: TipoRuido
  INTEGER :: TipoDatosMapaRuido
  CHARACTER(len=100) :: FileMapRuido

  INTEGER :: NB

END MODULE DatosProblema

PROGRAM QMLTEB

  USE ifport
  USE DatosProblema

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE CargaProblema(MuestraMemoria, Problema)
       INTEGER, INTENT(out) :: MuestraMemoria
       CHARACTER(LEN=100), INTENT(IN) :: Problema
     END SUBROUTINE CargaProblema

     SUBROUTINE CargaMascara(nside, Lugar, File, NumMapaCargar, Mascara)
       INTEGER, INTENT(in) :: Nside, NumMapaCargar
       CHARACTER(len=100), INTENT(in) :: Lugar, File
       INTEGER, DIMENSION(12*nside**2), INTENT(out) :: Mascara
     END SUBROUTINE CargaMascara

  END INTERFACE


  INTEGER :: NPixT, NPixP, MuestraMemoria, TotalArmonicos
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MascaraT, MascaraP

  REAL(kind=8) t,p,L, A, B

  !*************************************************************************

  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
     WRITE(*,*)'ERROR, ONE COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
     STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,ArchivoProblema)  

  !Carga problema
  !write(*,110)
  WRITE(*,*) "Loading configuration from: ", TRIM(ArchivoProblema)
  WRITE(*,*)
  CALL CargaProblema( MuestraMemoria, ArchivoProblema)

  WRITE(*,*)
  WRITE(*,*) "**************************************************"
  WRITE(*,*) "Loading masks"
  WRITE(*,*)

  ALLOCATE(MascaraT(12*NSide**2))
  CALL CargaMascara(nside, Lugar, FileMaskIntensity, NumMaskMapIntensity, MascaraT)
  NPixT = SUM(MascaraT)
  WRITE(*,'(A,T60,I8)') "   Number of observed pixels in temperature: ", NPIXT

  ALLOCATE(MascaraP(12*NSide**2))
  CALL CargaMascara(nside, Lugar, FileMaskPolarization, NumMaskMapPolarization, MascaraP)
  NPixP = SUM(MascaraP)
  WRITE(*,'(A,T60,I8)') "   Number of observed pixels in polarization: ", NPIXP

  ! npixt = 29009
  ! npixp = 29009
  ! DlMax = 192

  lmax = DlMax
  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  WRITE(*,'(A,T60,I8)') "   Lmax power spectrum:                       ", lmax
  WRITE(*,'(A,T60,I8)') "   Number of spherical harmonics:             ", TotalArmonicos

  !************************************

  WRITE(*,*)
  WRITE(*,*) "**************************************************"


  t = npixt
  p = npixp
  L = TotalArmonicos

  !TEB
  A = t*t + 2*t*p + 3*p*p + t*L + 4*p*L + 2*(t+2*p)*(3*L)
  B = 2*( 6*L*L + 4*p*L +  4*p*L)

  A = 8D0*A/1024D0**3
  B = 8D0*B/1024D0**3

  WRITE(*,*)
  WRITE(*,*) "**************************************************"
  WRITE(*,*) "ECLIPSE_TEB"

  WRITE(*,*) "Memory equation [Gb]: ", A
  WRITE(*,*) "Memory equation [Gb]: ", B

  !EB
  A = (2*p)**2 + 2*(2*p) * (2*L) + (2*p)*(2*L)
  B = 2*(2*(2*p)*(2*L) + 3*(L*L))

  A = 8D0*A/1024D0**3
  B = 8D0*B/1024D0**3

  WRITE(*,*)
  WRITE(*,*) "**************************************************"
  WRITE(*,*) "ECLIPSE_EB"

  WRITE(*,*) "Memory equation [Gb]: ", A
  WRITE(*,*) "Memory equation [Gb]: ", B



  !T
  L = L + 3
  A = t*t + 2*t*L + t*L
  B = 4*t*L + 2*L*L

  A = 8D0*A/1024D0**3
  B = 8D0*B/1024D0**3

  WRITE(*,*)
  WRITE(*,*) "**************************************************"
  WRITE(*,*) "ECLIPSE_T"

  WRITE(*,'(A,T60,I8)') "   Number of spherical harmonics:             ", TotalArmonicos+3

  WRITE(*,*) "Memory equation [Gb]: ", A
  WRITE(*,*) "Memory equation [Gb]: ", B

  WRITE(*,*)
  WRITE(*,*) "**************************************************"


END PROGRAM QMLTEB


FUNCTION getsize_fits(filename, extno_in, nside, nmaps)


  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)             :: filename
  INTEGER, INTENT(in) :: extno_in
  INTEGER, INTENT(out) :: nside, nmaps

  INTEGER, PARAMETER :: MAXDIM = 200

  INTEGER :: status,unit,readwrite
  CHARACTER(LEN=80) :: comment
  LOGICAL     ::  extend

  INTEGER :: hdunum, naxis, hdutype
  INTEGER :: nmove
  INTEGER :: nrows, tfields, varidat, rowlen

  CHARACTER(LEN=20), DIMENSION(1:MAXDIM) :: ttype, tform, tunit
  INTEGER,      DIMENSION(1:MAXDIM) :: tbcol

  CHARACTER(LEN=20)                      :: extname

  INTEGER :: repeat1, repeat2, datacode, width
  INTEGER :: nsize


  INTEGER           :: nmaps_in, ordering_in, nside_in
  INTEGER           :: obs_npix_in, ftype_in

  CHARACTER(LEN=20)      :: order_val !, object_val, ttype_val, polcconv_val

  INTEGER :: getsize_fits

  EXTERNAL :: ftgiou, ftnopn, ftclos, ftghdn, ftgkyj, ftgkyl, ftghdt, ftmrhd, printerror
  EXTERNAL :: ftghbn, ftghtb, ftbnfm, ftgnrw, ftgkys, ftgerr


  !-----------------------------------------------------------------------
  ttype = ""
  tform = ""
  tunit = ""

  status=0


  readwrite=0

  CALL ftgiou(unit,status)

  CALL ftnopn(unit,filename,readwrite,status)
  IF (status > 0) THEN
     ftype_in = -1
     getsize_fits = -1
     CALL printerror(status)
     CALL ftclos(unit, status)
     RETURN
  ENDIF
  !     -----------------------------------------
  CALL ftghdn(unit, hdunum)
  IF (hdunum == 1) THEN  ! in primary HDU: move to next HDU
     !     determines the presence of image
     CALL ftgkyj(unit,'NAXIS', naxis, comment, status)

     !     determines the presence of an extension
     CALL ftgkyl(unit,'EXTEND', extend, comment, status)
     IF (status > 0) THEN
        ftype_in = 0
        status = 0 ! no extension : 
        !     to be compatible with first version of the code
     ENDIF

  ELSE ! already in non primary HDU
     extend = .TRUE.
     CALL ftghdt(unit, hdutype, status)
  ENDIF


  !*******************
  IF (extend) THEN 

     nmove =  extno_in
     IF (hdunum == 1) nmove = extno_in + 1
     CALL ftmrhd(unit, nmove, hdutype, status)
     IF (status > 0) THEN ! extension not found
        !print*,'Extension #',extno_in,' not found in '//trim(filename)
        WRITE(*,*) "Error, no hay esa extension"
        CALL printerror(status)
        CALL ftclos(unit, status)
     ENDIF

     !        reads all the keywords
     IF (hdutype == 2) THEN ! binary table
        CALL ftghbn(unit, MAXDIM, nrows, tfields, ttype, tform, tunit, extname, varidat, status)
     ELSE ! ASCII table (hdutype = 1)
        ftype_in = 1
        CALL ftghtb(unit, MAXDIM, rowlen, nrows, tfields, ttype, tbcol, tform, tunit, extname, status)
     ENDIF

     WRITE(*,*) "****************"
     WRITE(*,*) "HDUTYPE: ", hdutype
     WRITE(*,*) "MAXDIM: ", MAXDIM
     WRITE(*,*) "nrows: ", nrows
     WRITE(*,*) "tfields: ", tfields 
     WRITE(*,*) "ttype: ", ttype(1:tfields)
     WRITE(*,*) "tform: ", tform(1:tfields)
     WRITE(*,*) "tunit: ", tunit(1:tfields)
     WRITE(*,*) "****************"
     WRITE(*,*) 

     !Revisar
     CALL ftgnrw(unit, nrows, status)

     !  parse TFORM keyword to find out the length of the column vector
     repeat1 = 1
     repeat2 = 1
     CALL ftbnfm(tform(1), datacode, repeat1, width, status)

     IF (tfields > 1) CALL ftbnfm(tform(2), datacode, repeat2, width, status)

     nsize = nrows * MAX(repeat1,repeat2)

     nmaps_in = tfields

     CALL ftgkys(unit,'ORDERING',order_val,comment,status)
     IF (status == 202) THEN ! Ordering not found
        ordering_in = 0
        order_val = ''
        status = 0
     ENDIF

     CALL ftgkyj(unit,'NSIDE',nside_in,comment,status)
     IF (status == 202) THEN ! Nside not found
        nside_in = -1
        status = 0
     ENDIF

     CALL ftgkyj(unit,'OBS_NPIX',obs_npix_in,comment,status)
     IF (status == 202) THEN ! obs_npix not found
        obs_npix_in = -1
        status = 0
     ENDIF

     nside = nside_in
     nmaps = nmaps_in
     getsize_fits = NSIZE

  END IF

  CALL ftclos(unit, status)

END FUNCTION getsize_fits


SUBROUTINE printerror(status)

  IMPLICIT NONE

  !  This subroutine prints out the descriptive text corresponding to the
  !  error status value and prints out the contents of the internal
  !  error message stack generated by FITSIO whenever an error occurs.

  INTEGER status
  CHARACTER errtext*30,errmessage*80

  EXTERNAL :: ftgerr, ftgmsg

  !  Check if status is OK (no error); if so, simply return
  IF (status .LE. 0)RETURN

  !  The FTGERR subroutine returns a descriptive 30-character text string that
  !  corresponds to the integer error status number.  A complete list of all
  !  the error numbers can be found in the back of the FITSIO User's Guide.
  CALL ftgerr(status,errtext)
  PRINT *,'FITSIO Error Status =',status,': ',errtext

  !  FITSIO usually generates an internal stack of error messages whenever
  !  an error occurs.  These messages provide much more information on the
  !  cause of the problem than can be provided by the single integer error
  !  status value.  The FTGMSG subroutine retrieves the oldest message from
  !  the stack and shifts any remaining messages on the stack down one
  !  position.  FTGMSG is called repeatedly until a blank message is
  !  returned, which indicates that the stack is empty.  Each error message
  !  may be up to 80 characters in length.  Another subroutine, called
  !  FTCMSG, is available to simply clear the whole error message stack in
  !  cases where one is not interested in the contents.
  CALL ftgmsg(errmessage)
  DO WHILE (errmessage .NE. ' ')
     PRINT *,errmessage
     CALL ftgmsg(errmessage)
  END DO
END SUBROUTINE printerror


SUBROUTINE CargaMascara(nside, Lugar, File, NumMapaCargar, Mascara)

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: imap
       REAL(kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
       !real(kind=8),                                intent(OUT) :: nullval
       INTEGER                   , INTENT(IN) :: extno
     END SUBROUTINE read_bintab_UnMapa
  END INTERFACE


  INTEGER, INTENT(in) :: Nside, NumMapaCargar
  CHARACTER(len=100), INTENT(in) :: Lugar, File
  INTEGER, DIMENSION(12*nside**2), INTENT(out) :: Mascara

  CHARACTER(len=200) :: ArchivoMascara
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: MascaraReal
  INTEGER :: nmaps, extno, npix_full

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, igebs2d, igebr2d

  !*********************************************************************

  !Codigo
  npix_full = 12*Nside**2

  !if(iam==0) then

  ArchivoMascara = TRIM(Lugar)//"/"//TRIM(File)

  ALLOCATE(MascaraReal(0:npix_full-1))

  nmaps = NumMapaCargar
  extno = 0
  CALL read_bintab_UnMapa(ArchivoMascara, MascaraReal, npix_full, NumMapaCargar, extno)

  WHERE(MascaraReal < -1.637E+030) MascaraReal = 0 

  Mascara(1:npix_full) = INT(MascaraReal(0:npix_full-1))

  DEALLOCATE(MascaraReal)

  !open(unit=20, file=Trim(lugar)//"/Mascara.dat", action="write")
  !write(20,*) Mascara
  !close(20)

  !end if

END SUBROUTINE CargaMascara


SUBROUTINE read_bintab(filename, map, npixtot, nmaps, extno)

  IMPLICIT NONE

  !=======================================================================
  CHARACTER(len=*),                          INTENT(IN)  :: filename
  INTEGER,                              INTENT(IN)  :: npixtot
  INTEGER,                              INTENT(IN)  :: nmaps
  REAL(kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
  INTEGER, INTENT(IN) :: extno

  INTEGER :: status,unit,readwrite,naxes(2),nfound, naxis
  INTEGER :: group, firstpix, i, npix32
  REAL(kind=8)   :: blank, testval
  REAL(kind=8)     :: bscale,bzero
  CHARACTER(len=80) :: comment
  LOGICAL(kind=8) :: extend
  INTEGER :: nmove, hdutype, hdunum
  INTEGER :: frow, imap
  INTEGER :: datacode, width
  LOGICAL(kind=8) ::  anynull_i

  INTEGER,     PARAMETER            :: MAXDIM = 400 !MAXDIM_TOP !number of columns in the extension
  INTEGER                      :: npix_old
  INTEGER, DIMENSION(1:MAXDIM) :: npix
  INTEGER, DIMENSION(1:MAXDIM) :: i0, i1
  INTEGER, DIMENSION(1:MAXDIM) :: repeat
  INTEGER                      :: nrow2read, nelem

  INTEGER                           :: nrows, tfields, varidat
  CHARACTER(len=20), DIMENSION(1:MAXDIM) :: ttype, tform, tunit
  CHARACTER(len=20)                      :: extname
  CHARACTER(len=*), PARAMETER            :: code='read_bintab'

  REAL(kind=8) :: ABS

  LOGICAL :: anynull
  REAL(kind=8) ::  nullval

  EXTERNAL :: ftnopn, ftghdn, ftgkyj, ftmrhd, ftgkyl, ftghdt, ftghbn, ftgkyd
  EXTERNAL :: ftgrsz, ftbnfm, ftgcvd, ftclos, ftgknj, ftgpvd, printerror

  nullval =-1233.33D0

  !-----------------------------------------------------------------------


  status=0

  unit = 146
  naxes(1) = 1
  naxes(2) = 1
  nfound = -1
  anynull = .FALSE.
  bscale = 1.0d0
  bzero = 0.0d0
  blank = -2.e25
  !nullval = bscale*blank + bzero
  comment=''
  ttype=''
  tform=''
  tunit=''
  extname=''


  readwrite=0
  CALL ftnopn(unit,filename,readwrite, status) 
  IF (status > 0) CALL printerror(status)
  !     -----------------------------------------
  CALL ftghdn(unit, hdunum)

  IF (hdunum == 1) THEN  ! in primary HDU
     !     determines the presence of image
     CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
     IF (status > 0) CALL printerror(status)

     !     determines the presence of an extension
     CALL ftgkyl(unit,'EXTEND', extend, comment, status)
     IF (status > 0) THEN
        extend = .FALSE.
        status = 0 ! no extension : 
        !     to be compatible with first version of the code
     ENDIF
  ENDIF

  IF (naxis > 0 .AND. .NOT.extend .AND. hdunum==1) THEN ! there is an image
     !        determine the size of the image (look naxis1 and naxis2)
     CALL ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

     !        check that it found only NAXIS1
     IF (nfound == 2 .AND. naxes(2) > 1) THEN
        PRINT *,'multi-dimensional image'
        PRINT *,'expected 1-D data.'
        !call fatal_error
     END IF

     IF (nfound < 1) THEN
        CALL printerror(status)
        PRINT *,'can not find NAXIS1.'
        !call fatal_error
     ENDIF

     npix(1)=naxes(1)
     IF (npix(1) /= npixtot) THEN
        PRINT *,'WARNING: found ',npix(1),' pixels in '//TRIM(filename)
        PRINT *,'         expected ',npixtot
        npix(1) = MIN(npix(1), npixtot)
        PRINT *,'         only ',npix(1),' will be read'
     ENDIF

     CALL ftgkyd(unit,'BSCALE',bscale,comment,status)
     IF (status == 202) THEN ! BSCALE not found
        bscale = 1.0d0
        status = 0
     ENDIF
     CALL ftgkyd(unit,'BZERO', bzero, comment,status)
     IF (status == 202) THEN ! BZERO not found
        bzero = 0.0d0
        status = 0
     ENDIF
     CALL ftgkyd(unit, 'BLANK', blank, comment, status)
     IF (status == 202) THEN ! BLANK not found 
        ! (according to fitsio BLANK is integer)
        blank = -2.e25
        status = 0
     ENDIF
     !nullval = bscale*blank + bzero

     !        -----------------------------------------

     group=1
     firstpix = 1
     npix32 = npix(1)
     CALL ftgpvd(unit, group, firstpix, npix32, nullval, map(0:npix(1)-1,1), anynull, status)
     ! if there are any NaN pixels, (real data)
     ! or BLANK pixels (integer data) they will take nullval value
     ! and anynull will switch to .true.
     ! otherwise, switch it by hand if necessary
     testval = 1.e-6 * ABS(nullval)
     DO i=0, npix(1)-1
        IF (ABS(map(i,1)-nullval) < testval) THEN
           anynull = .TRUE.
           GOTO 111
        ENDIF
     ENDDO
111  CONTINUE

  ELSE IF (extend .OR. hdunum>1) THEN

     IF (hdunum == 1) THEN
        !nmove = +1
        nmove = +1 + extno
        CALL ftmrhd(unit, nmove, hdutype, status)
     ELSE
        CALL ftghdt(unit, hdutype, status)
     ENDIF

     !call assert (hdutype==2, 'this is not a binary table')

     !        reads all the keywords
     CALL ftghbn(unit, MAXDIM, &
          &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
          &        status)


     IF (tfields < nmaps) THEN
        PRINT *,'found ',tfields,' maps in file '//TRIM(filename)
        PRINT *,'expected ',nmaps
        !call fatal_error
     ENDIF

     !        finds the bad data value
     CALL ftgkyd(unit, 'BAD_DATA', nullval, comment, status)
     IF (status == 202) THEN ! bad_data not found
        !if (KMAP == SP) nullval = s_bad_value ! default value
        !if (KMAP == DP) nullval = d_bad_value ! default value
        status = 0
     ENDIF


     npix_old = npixtot
     DO imap = 1, nmaps
        !parse TFORM keyword to find out the length of the column vector
        CALL ftbnfm(tform(imap), datacode, REPEAT(imap), width, status)
        npix(imap) = nrows * REPEAT(imap)
        IF (npix(imap) /= npixtot .AND. npix_old /= npix(imap)) THEN
           PRINT *,'WARNING: found ',npix(imap),' pixels in '//TRIM(filename)//', column ',imap
           PRINT *,'         expected ',npixtot,' or ',npix_old
           npix_old = npix(imap)
           npix(imap) = MIN(npix(imap), npixtot)
           PRINT *,'         only  ',npix(imap),' will be read'
        ENDIF
     ENDDO


     CALL ftgrsz(unit, nrow2read, status)
     nrow2read = MAX(nrow2read, 1)
     firstpix  = 1  ! starting position in FITS within row, 1 based
     i0(:) = 0  ! starting element in array, 0 based
     DO frow = 1, nrows, nrow2read

        !write(*,*) "Filas: ", frow, nrows, nrow2read

        DO imap = 1, nmaps

           !write(*,*) "Mapa: ", imap

           i1(imap) = MIN(i0(imap) + nrow2read * REPEAT(imap), npix(imap)) - 1
           nelem = i1(imap) - i0(imap) + 1
           !write(*,*) imap, frow, firstpix,nelem,i0(imap), i1(imap)
           CALL ftgcvd(unit, imap, frow, firstpix, nelem, &
                & nullval, map(i0(imap):i1(imap),imap), anynull_i, status)

           anynull = anynull .OR. anynull_i
           i0(imap) = i1(imap) + 1
        ENDDO
     ENDDO

     ! sanity check
     DO imap = 1, nmaps
        IF (i0(imap) /= npix(imap)) THEN
           !call fatal_error('something wrong during piece wise reading')
        ENDIF
     ENDDO

  ELSE ! no image no extension, you are dead, man
     !call fatal_error(' No image, no extension in '//trim(filename))
  ENDIF
  !     close the file
  CALL ftclos(unit, status)

  !     check for any error, and if so print out error messages
  IF (status > 0) CALL printerror(status)

  RETURN
END SUBROUTINE read_bintab


SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)

  IMPLICIT NONE

  !=======================================================================
  CHARACTER(len=*),                          INTENT(IN)  :: filename
  INTEGER,                              INTENT(IN)  :: npixtot
  INTEGER,                              INTENT(IN)  :: imap
  REAL(kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
  INTEGER                   , INTENT(IN) :: extno

  INTEGER :: status,unit,readwrite,naxes(2),nfound, naxis
  INTEGER ::  firstpix
  REAL(kind=8)   :: blank
  REAL(kind=8)     :: bscale,bzero
  CHARACTER(len=80) :: comment
  LOGICAL(kind=8) :: extend
  INTEGER :: nmove, hdutype, hdunum
  INTEGER :: frow
  INTEGER :: datacode, width
  LOGICAL(kind=8) ::  anynull_i

  INTEGER,     PARAMETER            :: MAXDIM = 400 !MAXDIM_TOP !number of columns in the extension
  INTEGER                      :: npix_old
  INTEGER, DIMENSION(1:MAXDIM) :: npix
  INTEGER, DIMENSION(1:MAXDIM) :: i0, i1
  INTEGER, DIMENSION(1:MAXDIM) :: repeat
  INTEGER                      :: nrow2read, nelem

  INTEGER                           :: nrows, tfields, varidat
  CHARACTER(len=20), DIMENSION(1:MAXDIM) :: ttype, tform, tunit
  CHARACTER(len=20)                      :: extname
  CHARACTER(len=*), PARAMETER            :: code='read_bintab'



  EXTERNAL :: ftnopn, ftghdn, ftgkyj, ftmrhd, ftgkyl, ftghdt, ftghbn, ftgkyd
  EXTERNAL :: ftgrsz, ftbnfm, ftgcvd, ftclos,  printerror, ftgiou, ftfiou

  LOGICAL :: anynull
  REAL(kind=8) ::  nullval
  nullval = -1.6375D30


  !-----------------------------------------------------------------------


  status=0


  naxes(1) = 1
  naxes(2) = 1
  nfound = -1
  anynull = .FALSE.
  bscale = 1.0d0
  bzero = 0.0d0
  blank = -2.e25

  comment=''
  ttype=''
  tform=''
  tunit=''
  extname=''


  CALL ftgiou(unit,status)

  readwrite=0
  CALL ftnopn(unit,filename,readwrite, status) 
  IF (status > 0) CALL printerror(status)
  !     -----------------------------------------
  CALL ftghdn(unit, hdunum)

  IF (hdunum == 1) THEN  ! in primary HDU
     !     determines the presence of image
     CALL ftgkyj(unit,'NAXIS', naxis, comment, status)
     IF (status > 0) CALL printerror(status)

     !     determines the presence of an extension
     CALL ftgkyl(unit,'EXTEND', extend, comment, status)
     IF (status > 0) THEN
        extend = .FALSE.
        status = 0 ! no extension : 
     ENDIF
  ENDIF


  IF (extend .OR. hdunum>1) THEN

     IF (hdunum == 1) THEN
        nmove = +1
        nmove = +1 + extno
        CALL ftmrhd(unit, nmove, hdutype, status)
     ELSE
        CALL ftghdt(unit, hdutype, status)
     ENDIF

     !call assert (hdutype==2, 'this is not a binary table')

     !        reads all the keywords
     CALL ftghbn(unit, MAXDIM, &
          &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
          &        status)


     !        finds the bad data value
     CALL ftgkyd(unit, 'BAD_DATA', nullval, comment, status)

     IF (status == 202) THEN ! bad_data not found
        !if (KMAP == SP) nullval = s_bad_value ! default value
        !if (KMAP == DP) nullval = d_bad_value ! default value
        nullval = -1.6375D30
        status = 0
     ENDIF


     npix_old = npixtot
     !parse TFORM keyword to find out the length of the column vector
     CALL ftbnfm(tform(imap), datacode, REPEAT(imap), width, status)
     npix(imap) = nrows * REPEAT(imap)
     IF (npix(imap) /= npixtot .AND. npix_old /= npix(imap)) THEN
        PRINT *,'WARNING: found ',npix(imap),' pixels in '//TRIM(filename)//', column ',imap
        PRINT *,'         expected ',npixtot,' or ',npix_old
        npix_old = npix(imap)
        npix(imap) = MIN(npix(imap), npixtot)
        PRINT *,'         only  ',npix(imap),' will be read'
     ENDIF


     CALL ftgrsz(unit, nrow2read, status)
     nrow2read = MAX(nrow2read, 1)
     firstpix  = 1  ! starting position in FITS within row, 1 based
     i0(:) = 0  ! starting element in array, 0 based
     DO frow = 1, nrows, nrow2read

        i1(imap) = MIN(i0(imap) + nrow2read * REPEAT(imap), npix(imap)) - 1
        nelem = i1(imap) - i0(imap) + 1
        CALL ftgcvd(unit, imap, frow, firstpix, nelem, &
             & nullval, MAPA(i0(imap):i1(imap)), anynull_i, status)

        anynull = anynull .OR. anynull_i
        i0(imap) = i1(imap) + 1

     ENDDO

  ENDIF
  !     close the file
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  !     check for any error, and if so print out error messages
  IF (status > 0) CALL printerror(status)

  RETURN
END SUBROUTINE read_bintab_UnMapa


SUBROUTINE CargaProblema(MuestraMemoria, Problema)

  USE DatosProblema

  IMPLICIT NONE

  INTEGER, INTENT(out) :: MuestraMemoria
  CHARACTER(LEN=100), INTENT(IN) :: Problema

  INTEGER, DIMENSION(20) :: Paso

  CHARACTER(len=200) :: cin
  INTEGER :: Leer

  CHARACTER(LEN=100) :: Texto

  EXTERNAL :: blacs_barrier, igebs2d, igebr2d, dgebs2d, dgebr2d

  !if(iam==0) then

  TipoFileMaps = -1000
  QuitarSesgo = -1000
  RuidoTT = -1000
  RuidoQQ = -1000

  OPEN(unit=10, file=TRIM(Problema), action="read")

  Leer = 1
  DO WHILE(Leer==1)

     READ(10,'(A)', END=100) cin

     IF(cin(1:1)/="#") THEN 

        !Entero
        CALL CargaEntero(cin, "NSide = ", NSide)
        CALL CargaEntero(cin, "Lmax_Covariance_Matrix = ", Lmax)
        CALL CargaEntero(cin, "Lmax_Power_Spectrum = ", DlMax)
        CALL CargaEntero(cin, "Type_of_Analysis = ", TipoAnalisis)
        CALL CargaEntero(cin, "Number_of_Maps = ", NMapas)
        CALL CargaEntero(cin, "Kind_of_Data = ", TipoFileMaps)
        CALL CargaEntero(cin, "Pixel_Window = ", PixelWindowEnMapas)
        CALL CargaEntero(cin, "Intensity_Mask_NumMap = ", NumMaskMapIntensity)
        CALL CargaEntero(cin, "Polarization_Mask_NumMap = ", NumMaskMapPolarization)
        CALL CargaEntero(cin, "Kind_of_Noise = ", TipoRuido)
        CALL CargaEntero(cin, "Kind_of_Noise_Data = ", TipoDatosMapaRuido)
        CALL CargaEntero(cin, "Remove_Noise_Biass = ", QuitarSesgo)
        CALL CargaEntero(cin, "Matrices_Cyclic_Block_Size = ", NB)
        CALL CargaEntero(cin, "Inverse_Covariance_Matrix_Control = ", ControlInversaMC)
        CALL CargaEntero(cin, "Show_Memory_Allocated = ", MuestraMemoria)
        CALL CargaEntero(cin, "Kind_of_Bin_Center = ", KindOfBinCenter)
        CALL CargaEntero(cin, "Kind_of_Grouping = ", KindOfGrouping)
        CALL CargaEntero(cin, "Compute_Fisher_Matrix = ", ComputeFisherMatrix)
        CALL CargaEntero(cin, "Binned = ", Binned)


        !Real
        CALL CargaReal(cin, "Intensity_Noise = ", RuidoTT)
        CALL CargaReal(cin, "Polarization_Noise = ", RuidoQQ)
        CALL CargaReal(cin, "Beam_FWHM = ", FWHM_Beam)

        !Cadenas
        CALL CargaCadena(cin, "Data_Folder = ", Lugar)
        CALL CargaCadena(cin, "Maps_FileName = ", FileMaps)
        CALL CargaCadena(cin, "Maps2Cross_FileName = ", FileMapsCross)
        CALL CargaCadena(cin, "Fiducial_FileName = ", ArchivoFiducial)
        CALL CargaCadena(cin, "Intensity_Mask_FileName = ", FileMaskIntensity)
        CALL CargaCadena(cin, "Polarization_Mask_FileName = ", FileMaskPolarization)
        CALL CargaCadena(cin, "Noise_Map_FileName = ", FileMapRuido)
        CALL CargaCadena(cin, "Healpix_Data_Folfer = ", DirHealpixData)
        CALL CargaCadena(cin, "Binnig_Limits_FileName = ", FileBinLimits)


     END IF

     GOTO 200
100  Leer = 0
200  CONTINUE

  END DO

  CLOSE(10)


  Paso(1) = NSide
  Paso(2) = LMax
  Paso(3) = DlMax
  Paso(4) = TipoAnalisis
  Paso(5) = NMapas
  Paso(6) = QuitarSesgo
  Paso(7) = NB
  Paso(8) = ControlInversaMC
  Paso(9) = TipoRuido
  Paso(10) = TipoDatosMapaRuido
  Paso(11) = MuestraMemoria
  Paso(12) = KindOfBinCenter
  Paso(13) = KindOfGrouping
  Paso(14) = ComputeFisherMatrix
  Paso(15) = Binned

  !end if



  !*********************
300 FORMAT (1X,A35,TR2,I)
310 FORMAT (1X,A35,TR2,F)
320 FORMAT (1X,A35,TR2,A)

  !if(iam==0) then

  WRITE(*,320) "Data_Folder: ",TRIM(Lugar)
  WRITE(TEXTO,*) NSide
  WRITE(*,320) "NSide: ", ADJUSTL(TRIM(Texto))
  !write(*,320) "Fiducial_FileName: ", Trim(ArchivoFiducial)
  ! write(TEXTO,*) Lmax
  ! write(*,320) "Lmax_Covariance_Matrix: ", Adjustl(Trim(Texto))
  WRITE(TEXTO,*) Dlmax
  WRITE(*,320) "Lmax_Power_Spectrum: ", ADJUSTL(TRIM(Texto))
  ! write(TEXTO,*) TipoAnalisis
  ! write(*,320) "Type_of_Analysis: ", Adjustl(Trim(Texto))
  WRITE(*,*)

  ! write(*,320) "Maps_FileName: ", Trim(FileMaps)

  ! if(TipoAnalisis==1) then
  ! write(*,320) "Maps2Cross_FileName: ", Trim(FileMapsCross)
  ! end if

  ! write(TEXTO,*) NMapas
  ! write(*,320) "Number_of_Maps: ", Adjustl(Trim(Texto))
  !write(*,300) "Kind_of_Data: ", TipoFileMaps

  ! write(TEXTO,*) PixelWindowEnMapas
  ! write(*,320) "Pixel_Window: ", Adjustl(Trim(Texto))
  ! write(TEXTO,*) FWHM_Beam
  ! write(*,320) "Beam_FWHM: ", Adjustl(Trim(Texto))
  ! write(*,*)

  WRITE(*,320) "Intensity_Mask_FileName: ", TRIM(FileMaskIntensity)
  WRITE(TEXTO,*) NumMaskMapIntensity
  WRITE(*,320) "Intensity_Mask_NumMap: ", ADJUSTL(TRIM(Texto))
  WRITE(*,320) "Polarization_Mask_FileName: ", TRIM(FileMaskPolarization)
  WRITE(TEXTO,*) NumMaskMapPolarization
  WRITE(*,320) "Polarization_Mask_NumMap: ", ADJUSTL(TRIM(Texto))
  WRITE(*,*)


  ! write(TEXTO,*) TipoRuido
  ! write(*,320) "Kind_of_Noise: ", Adjustl(Trim(Texto))
  ! if(TipoRuido==0) then
  ! write(TEXTO,*) RuidoTT
  ! write(*,320) "Intensity_Noise:", Adjustl(Trim(Texto))
  ! write(TEXTO,*) RuidoQQ
  ! write(*,320) "Polarization_Noise: ", Adjustl(Trim(Texto))
  ! else
  ! write(*,320) "Noise_Map_FileName: ", Trim(FileMapRuido)
  !write(*,300) "Kind_of_Noise_Data: ", TipoDatosMapaRuido
  ! end if
  ! write(TEXTO,*) QuitarSesgo
  ! write(*,320) "Remove_Noise_Biass: ", Adjustl(Trim(Texto))
  ! write(*,*)


  ! write(TEXTO,*) NB
  ! write(*,320) "Matrices_Cyclic_Block_Size: ", Adjustl(Trim(Texto))
  ! write(TEXTO,*) ControlInversaMC
  ! write(*,320) "Inverse_Covariance_Matrix_Control: ", Adjustl(Trim(Texto))
  ! write(TEXTO,*) MuestraMemoria
  ! write(*,320) "Show_Memory_Allocated: ", Adjustl(Trim(Texto))
  ! write(*,320) "Healpix_Data_Folfer: ", Trim(DirHealpixData)

  !end if

  !*********************

  RETURN

CONTAINS
SUBROUTINE CargaEntero(cin, dato, valor)
            
   IMPLICIT NONE
   
   CHARACTER(len=100), INTENT(in) :: cin
   CHARACTER(len=*), INTENT(in) :: dato
   INTEGER, INTENT(inout) :: valor
   
   INTEGER :: p
   CHARACTER(len=100) :: cad
   
   p =INDEX(cin, dato)
   
   IF(p/=0) THEN
       
       p =  LEN(dato)
       cad = cin(p+1:LEN(TRIM(cin)))
     
       if(Len(Trim(cad))>0) then
         READ(cad,*)  valor
       else
         valor = -123456789
       endif
     
   END IF
         
END SUBROUTINE CargaEntero


SUBROUTINE CargaReal(cin, dato, valor)
   
   IMPLICIT NONE
   
   CHARACTER(len=100), INTENT(in) :: cin
   CHARACTER(len=*), INTENT(in) :: dato
   REAL(kind=8), INTENT(inout) :: valor
   
   INTEGER :: p
   CHARACTER(len=100) :: cad
   
   p =INDEX(cin, dato)
   
   IF(p/=0) THEN
       
       p =  LEN(dato)
       cad = cin(p+1:LEN(TRIM(cin)))
     
       if(Len(Trim(cad))>0) then
         READ(cad,*)  valor
       else
         valor = -123456789
       endif
     
   END IF
     
END SUBROUTINE 


SUBROUTINE CargaCadena(cin, dato, valor)
   
   IMPLICIT NONE
   
   CHARACTER(len=100), INTENT(in) :: cin
   CHARACTER(len=*), INTENT(in) :: dato
   CHARACTER(len=100), INTENT(inout) :: valor
   
   INTEGER :: p
   CHARACTER(len=100) :: cad
   
   p =INDEX(cin, dato)
   
   IF(p/=0) THEN
       
       p =  LEN(dato)
       cad = cin(p+1:LEN(TRIM(cin)))

       if(Len(Trim(cad))>0) then
         READ(cad,*)  valor
       else
         valor = "!!!!!!!!!!!!!! Bad data !!!!!!!!!!!!!"
       endif
       
   END IF
   
END SUBROUTINE CargaCadena

END SUBROUTINE CargaProblema
