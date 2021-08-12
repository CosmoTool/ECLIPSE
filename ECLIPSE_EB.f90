!---------------------------------------------------------------------------------
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
!     gaussbeam
!     getsize_fits
!     pix2_ang_ring
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
!     This program makes use of SCALAPack subroutines.
!     
!     ScaLAPACK is a software package provided by Univ. of Tennessee;
!     Univ. of California, Berkeley; Univ. of Colorado Denver; and NAG Ltd..
!     
!     
!     For more information about SCALAPack see http://www.netlib.org/scalapack/
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




MODULE ControlMemoria

  IMPLICIT NONE

  INTEGER :: MuestraMemoria
  REAL    :: TotalProcesador, TotalPrograma

CONTAINS
  SUBROUTINE SumaMemoria(ictxt, iam, Muestra, F, C, Unidades)

    !F: Numero de filas reservadas
    !C: Numero de columnas reservadas
    !Unidades: numero de bloques de esas dimensiones reservados

    IMPLICIT NONE

    INTEGER, INTENT(in) :: ictxt, iam, Muestra, F, C, Unidades

    REAL            :: Sumar, FReal, CReal, Mostrar
    REAL, PARAMETER :: Factor = 1.0/1024D0
    INTEGER         :: Caso

    EXTERNAL :: blacs_barrier, sgsum2d

    Caso = 0
    FReal = F
    CReal = C

    Sumar = 8.0*FReal*CReal*Unidades

    !Lo pasa a kilobytes
    TotalProcesador = TotalProcesador + Sumar * Factor
    TotalPrograma   = TotalProcesador


    CALL blacs_barrier(ictxt, 'All')
    CALL SGSUM2D( ictxt, 'All', '1-tree', 1, 1, TotalPrograma, 1, -1, -1)
    CALL blacs_barrier(ictxt, 'All')

    Mostrar = TotalPrograma
    IF(Mostrar>1048576) THEN
       Mostrar = Mostrar*Factor**2
       Caso = 2
    ELSEIF(Mostrar>1024) THEN
       Mostrar = Mostrar*Factor
       Caso = 1
    END IF

    !Elimina problemas de redondeo
    IF(Mostrar<0.000001) Mostrar = 0 


    IF((iam == 0).AND.(Muestra==1)) THEN

       IF(caso==0) THEN    
          WRITE(*,'(T78, F9.2, " KB")') Mostrar
       ELSEIF(caso==1) THEN
          WRITE(*,'(T78, F9.2, " MB")') Mostrar
       ELSE
          WRITE(*,'(T78, F9.2, " GB")') Mostrar
       ENDIF

    END IF

  END SUBROUTINE SumaMemoria

END MODULE ControlMemoria


MODULE DatosProblema

  IMPLICIT NONE

  INTEGER         :: NSide, Lmax, Dlmax
  INTEGER         :: TipoAnalisis, QuitarSesgo
  REAL(Kind=8)    :: RuidoTT, RuidoQQ
  REAL(Kind=8)    :: FWHM_Beam
  INTEGER         :: PixelWindowEnMapas
  INTEGER         :: TypeOfGrouping
  INTEGER         :: TypeOfBinCenter
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
  INTEGER :: NumMaskMapIntensity          !Numero de mapa en el archivo de las mascaras
  INTEGER :: NumMaskMapPolarization       !Numero de mapa en el archivo de las mascaras
  INTEGER :: TipoFileMaps                 !TQU file or only QU file

  INTEGER :: ControlInversaMC             !Para controlar diagonal producto IMC.MC

  INTEGER :: TipoRuido
  INTEGER :: TipoDatosMapaRuido
   CHARACTER(len=100) :: FileMapRuido

  INTEGER :: NB

END MODULE DatosProblema

PROGRAM ECLIPSE_EB

  USE ifport
  USE ControlMemoria
  USE DatosProblema

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE CargaProblema(ICTXT,iam, MuestraMemoria, Problema)
       INTEGER, INTENT(in) :: ICTXT, iam
       INTEGER, INTENT(out) :: MuestraMemoria
       CHARACTER(LEN=100), INTENT(IN) :: Problema
     END SUBROUTINE CargaProblema

     SUBROUTINE CargaMascara(ictxt, nside, Lugar, File, NumMapaCargar, Mascara)
       INTEGER, INTENT(in) :: ictxt, Nside, NumMapaCargar
       CHARACTER(len=100), INTENT(in) :: Lugar, File
       INTEGER, DIMENSION(12*nside**2), INTENT(out) :: Mascara
     END SUBROUTINE CargaMascara

     SUBROUTINE PuntosPixelesObservados(ictxt, nside, npix, Mascara, NB, NPixelesProcesador, z, a)
       INTEGER, INTENT(in) :: ictxt, nside, npix, NB, NPixelesProcesador
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
       REAL(Kind=8), DIMENSION(NPixelesProcesador) :: z
       REAL(Kind=8), DIMENSION(NPixelesProcesador) :: a
     END SUBROUTINE PuntosPixelesObservados

     SUBROUTINE BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)
       INTEGER, INTENT(in) :: ICTXT, PixelWindowEnMapas, nside, lmax
       REAL(Kind=8) :: FWHM_Beam
       CHARACTER(len=*), INTENT(in) :: DirHealpixData       
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(out) :: BeamPWDl
     END SUBROUTINE BeamPixelWindow

     SUBROUTINE CargaMultipolos(ICTXT, iam, Lugar, FileData, LMax, DlEE, DlBB, DlEB)
       INTEGER, INTENT(in) :: ICTXT, IAM, lmax
       CHARACTER(len=100), INTENT(in) :: Lugar, FileData
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(OUT) :: DlEE, DlBB, DlEB
     END SUBROUTINE CargaMultipolos

     SUBROUTINE CalculaBloquesMatrizArmonicos(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, &
          &YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)
       INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
       REAL(Kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
       REAL(Kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
     END SUBROUTINE CalculaBloquesMatrizArmonicos

     SUBROUTINE CargaMapaRuido(ICTXT, TipoRuido, RuidoQQ, nside, NFMC, NFBloque, NB, Mascara, Lugar, FileData, TipoDatos,&
          &  MapaRuido2, MapaRuido2Bloque)
       INTEGER, INTENT(in) :: ICTXT, nside, NFMC, NFBloque, NB, TipoRuido, TipoDatos
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
       CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
       REAL(Kind=8), INTENT(inout) :: RuidoQQ
       REAL(Kind=8), DIMENSION(NFMC), INTENT(out) :: MapaRuido2
       REAL(Kind=8), DIMENSION(NFBloque), INTENT(out) :: MapaRuido2Bloque
     END SUBROUTINE CargaMapaRuido

     SUBROUTINE CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NFBA, NCBA, DimMC, lmax, MapaRuido2, DlEE, DlBB, DlEB, &
          &YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, MC, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, NB, NFMC, NCMC, NFBA, NCBA, DimMC, lmax
       INTEGER(Kind=8) :: TiempoInicio
       REAL(Kind=8), DIMENSION(NFBA), INTENT(in) :: MapaRuido2
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: DlEE, DlBB, DlEB
       REAL(Kind=8), DIMENSION(NFBA,NCBA), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
       REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(out)  :: MC
     END SUBROUTINE CalculaMatrizCovarianza

     SUBROUTINE InvierteMC(ICTXT, DimMC, NB, NF, NC, MC, ControlInversa, Fallo, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NB, NF, NC, ControlInversa
       REAL(Kind=8), DIMENSION(NF,NC), INTENT(inout) :: MC
       INTEGER, INTENT(out) :: Fallo
       INTEGER(Kind=8), INTENT(in) ::  TiempoInicio
     END SUBROUTINE InvierteMC

     SUBROUTINE CargaMapas(ICTXT, DimMC, NMapas, NB, NFMC, NCMapas, Mapas, Lugar, Caso)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMapas, Caso
       REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CargaMapas

     SUBROUTINE CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileData, TipoDatos, Mapas)
       INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas, TipoDatos
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
       CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
       REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas
     END SUBROUTINE CargaMapasFits

     SUBROUTINE CargaMapasFits2(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileData, TipoDatos, Mapas)
       INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas, TipoDatos
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
       CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
       REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas
     END SUBROUTINE CargaMapasFits2

     SUBROUTINE CalculaIMCMapas(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
       REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
       REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaIMCMapas

     SUBROUTINE CalculaYl(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas, MXR, MXI, Mapas)
       INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas
       CHARACTER(LEN=100), INTENT(in) :: Lugar
       REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: MXR, MXI
       REAL(Kind=8), DIMENSION(NFMX, NCMapas), INTENT(in) :: Mapas
     END SUBROUTINE CalculaYl

     SUBROUTINE CalculaBl(ICTXT, Lugar, DimMC, lmax, NB, NFMX, NCMX, MapaRuido2, IMCMXR, IMCMXI)
       INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NB, NFMX, NCMX
       CHARACTER(LEN=100), INTENT(in) :: Lugar
       REAL(Kind=8), DIMENSION(NFMX), INTENT(in) :: MapaRuido2
       REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: IMCMXR, IMCMXI
     END SUBROUTINE CalculaBl

     SUBROUTINE CalculaMatrizArmonicosCompleja(ICTXT, NB, NF, NFMX, NCMX, npix, lmax, z, Phi, BeamPWDl, MXC)
       INTEGER, INTENT(IN) :: ICTXT, NB, NF, NFMX, NCMX, npix, lmax
       REAL(Kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
       COMPLEX(Kind=8), DIMENSION(NFMX,NCMX), INTENT(OUT) :: MXC
     END SUBROUTINE CalculaMatrizArmonicosCompleja

     SUBROUTINE CalculaMatrizFisher(ICTXT, lugar, lmax, NB, NFBloque, NCBloque, BloqueEE, BloqueBB, BloqueEB)
       INTEGER, INTENT(in) :: ICTXT, lmax, NB, NFBloque, NCBloque
       CHARACTER(LEN=100), INTENT(in) :: Lugar
       COMPLEX(Kind=8), DIMENSION(NFBloque, NCBloque), INTENT(in) :: BloqueEE, BloqueBB, BloqueEB
     END SUBROUTINE CalculaMatrizFisher

     SUBROUTINE CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)
       INTEGER(Kind=4), INTENT(in) :: IAM, lmax, NMapas, QuitarSesgo
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaEspectroPotencia

     SUBROUTINE CalculaIMCMapasMapasCross(ICTXT,Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas, MapasCross)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
       REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
       REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas, MapasCross
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaIMCMapasMapasCross

     SUBROUTINE CalculaYlParesMapas(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas, MXR, MXI, Mapas, MapasCross)
       INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas
       CHARACTER(LEN=100), INTENT(in) :: Lugar
       REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: MXR, MXI
       REAL(Kind=8), DIMENSION(NFMX, NCMapas), INTENT(in) :: Mapas, MapasCross
     END SUBROUTINE CalculaYlParesMapas

     SUBROUTINE CalculaBineado(NSide, npix, lmax, RuidoDiagonal, DlEE, DlBB, DlEB, BeamPWDl,&
          & TipoCompactado, TipoCentrado, QuitarSesgo, NMapas, Lugar, Filename)
       INTEGER, INTENT(in) :: lmax, NSide, npix, QuitarSesgo, NMapas
       INTEGER, INTENT(in) :: TipoCompactado
       INTEGER, INTENT(in) :: TipoCentrado
       REAL(Kind=8), INTENT(in) :: RuidoDiagonal
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: DlEE, DlBB, DlEB
       REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
       CHARACTER(len=100), INTENT(in) :: Lugar, Filename
     END SUBROUTINE CalculaBineado

  END INTERFACE

  !Grid
  INTEGER :: IAM, NPROCS, NPROW, NPCOL, ICTXT, MYROW, MYCOL

  !Tiempos
  INTEGER(Kind=8) :: TiempoInicio, TiempoFinal

  INTEGER, EXTERNAL :: NUMROC
  EXTERNAL :: BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO
  EXTERNAL :: blacs_barrier, BLACS_GRIDEXIT, BLACS_EXIT, pdgemr2d, descinit
  EXTERNAL :: pdsymm, pzgemm

  !Mascara
  INTEGER :: NPix
  INTEGER, ALLOCATABLE, DIMENSION(:) :: Mascara

  !Pixeles observados
  INTEGER :: NPixelesProcesador
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: z, a

  !Matrices de armonicos esfericos
  INTEGER :: TotalArmonicos, NColumnas, NFMX, NCMX, info
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MXReal, MXImaginaria
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: IMCMXR, IMCMXI
  INTEGER, DIMENSION(9) :: DescBloque, DescMX, DescMC


  !Multipolos
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: DlEE, DlBB, DlEB, BeamPWDl

  !Matriz de covarianza
  INTEGER :: DimMC, NFMC, NCMC, Fallo
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MC


  !Mapas
  INTEGER :: NCMapas
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Mapas, MapasCross
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: MapaRuido2                  !Mapa ruido dimensiones MC -> 2Npix
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: MapaRuido2Bloque            !Mapa ruido diemnsiones bloque de MC -> Npix


  !Para la matriz de Fisher
  !IMC.MX en forma compleja
  COMPLEX(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: IMCMX, MXC, BloqueEE, BloqueBB, BloqueEB
  COMPLEX(Kind=8) :: Alpha, Beta
  INTEGER :: NFB, NCB

  !Cargar numero de filas
  CHARACTER(100) :: num1char


  !*************************************************************************

  !Inicia el grid de procesos
  CALL BLACS_PINFO(IAM, NPROCS)

  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT().NE.2)THEN
     WRITE(*,*)'ERROR, TWO COMMAND-LINE ARGUMENTS REQUIRED, STOPPING'
     CALL BLACS_EXIT(0)
     STOP
  ENDIF

  CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the two values
  READ(num1char,*) NPROW
  NPCOL = NPROCS/NPROW

  CALL GET_COMMAND_ARGUMENT(2,ArchivoProblema)  


  TotalProcesador = 0
  TotalPrograma   = 0
  TiempoInicio    = Time()

100 FORMAT(2X, A, T79, I10, "s")
110 FORMAT(/,'*****************************************************************************************')
115 FORMAT(/)

  IF(iam==0) THEN
     WRITE(*,*) 
     WRITE(*,110) 
     WRITE(*,*) "Start time: ", CTime(Time())
     WRITE(*,*)
     WRITE(*,*) "  Number of cores:   ", NPROCS
     WRITE(*,*) "  Number of rows:    ", NPROW 
     WRITE(*,*) "  Number of columns: ", NPCOL
     WRITE(*,*) 
  END IF

  CALL BLACS_GET( -1, 0, ICTXT )
  CALL BLACS_GRIDINIT( ICTXT, 'R', NPROW, NPCOL)
  CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)

  !Carga problema
  IF(iam==0) THEN
     WRITE(*,110)
     WRITE(*,*) "Loading configuration from: ", TRIM(ArchivoProblema)
     WRITE(*,*)
  END IF
  CALL CargaProblema(ictxt, iam, MuestraMemoria, ArchivoProblema)

  IF(iam==0) THEN
     WRITE(*,110)
     WRITE(*,100) "Loading mask"
     WRITE(*,*)
  END IF

  ALLOCATE(Mascara(12*NSide**2))
  CALL CargaMascara(ictxt, nside, Lugar, FileMaskPolarization , NumMaskMapPolarization, Mascara)

  npix = SUM(Mascara)

  IF(iam==0) THEN
     WRITE(*,'(A,T82,I8)') "   Number of observed pixels in polarization: ", NPIX
  END IF

  !Posiciones pixeles observados
  NPixelesProcesador = MAX(1,NUMROC(npix, NB, MYROW, 0, NPROW))
  ALLOCATE(z(NPixelesProcesador))
  ALLOCATE(a(NPixelesProcesador))
  a = 0
  z = 0

  CALL PuntosPixelesObservados(ictxt, nside, npix, Mascara, NB, NPixelesProcesador, z, a)


  !Carga multipolos
  ALLOCATE(DlEE(0:lmax))
  ALLOCATE(DlBB(0:lmax))
  ALLOCATE(DlEB(0:lmax))
  ALLOCATE(BeamPWDl(0:lmax))

  DlEE = 1
  DLBB = 1
  DLEB = 0
  BeamPWDl = 1

  IF((iam==0).AND.(PixelWindowEnMapas==1)) THEN
     WRITE(*,*)
     WRITE(*,100) "Loading HEALPix Pixel Window"
  END IF

  CALL BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)

  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,100) "Loading Fiducial Power Spectrum "
  END IF

  CALL CargaMultipolos(ICTXT, iam, Lugar, ArchivoFiducial, LMax, DlEE, DlBB, DlEB)

  IF(iam==0) THEN
     WRITE(*,110)
     WRITE(*,115)
     WRITE(*,100) "Step 1. Computing the covariance matrix "
     WRITE(*,*)
     WRITE(*,100) "  Computing the spherical harmonics matrix ", Time()-TiempoInicio
  END IF


  !***************************************************************************************
  !Calcula la matriz de los armonicos esfericos para calcular matriz de covarianza
  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  !Los bloques tienen TotalArmonicos columnas, no 2*TotalArmonicos
  NColumnas = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))

  ALLOCATE(YQER(NPixelesProcesador, NColumnas))
  ALLOCATE(YQEI(NPixelesProcesador, NColumnas))
  ALLOCATE(YQBR(NPixelesProcesador, NColumnas))
  ALLOCATE(YQBI(NPixelesProcesador, NColumnas))
  ALLOCATE(YUER(NPixelesProcesador, NColumnas))
  ALLOCATE(YUEI(NPixelesProcesador, NColumnas))
  ALLOCATE(YUBR(NPixelesProcesador, NColumnas))
  ALLOCATE(YUBI(NPixelesProcesador, NColumnas))

  YQER=0
  YQEI=0
  YQBR=0
  YQBI=0
  YUER=0
  YUEI=0
  YUBR=0
  YUBI=0

  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NPixelesProcesador, NColumnas, 8)


  CALL CalculaBloquesMatrizArmonicos(ICTXT, NB, NPixelesProcesador, NColumnas, npix, lmax, z, a, BeamPWDl, &
       & YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)       

  !Matriz de covarianza
  DimMC = 2*NPix
  NFMC = MAX(1,1,NUMROC(DimMC, NB, MYROW, 0, NPROW))
  NCMC = MAX(1,1,NUMROC(DimMC, NB, MYCOL, 0, NPCOL))
  ALLOCATE(MC(NFMC, NCMC))
  MC = 1D0


  !Mapa de ruido en dos formatos
  ALLOCATE(MapaRuido2Bloque(NPixelesProcesador))
  ALLOCATE(MapaRuido2(NFMC))

  CALL blacs_barrier(ICTXT, "ALL")

  IF((iam==0).AND.(TipoRuido==1)) THEN
     WRITE(*,*)
     WRITE(*,100) "  Loading noise maps"
  END IF

  CALL CargaMapaRuido(ICTXT, TipoRuido, RuidoQQ, nside, NFMC, NPixelesProcesador, NB, Mascara, Lugar,&
       & FileMapRuido, TipoDatosMapaRuido, MapaRuido2, MapaRuido2Bloque)


  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,100) "  Computing blocks of the covariance matrix ", Time()-TiempoInicio
  END IF
  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMC, 1)


  !Calcula la matriz de covarianza a partir de los bloques de la matriz de armonicos
  CALL CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NPixelesProcesador, NColumnas, DimMC, lmax, MapaRuido2Bloque,&
       & DlEE, DlBB, DlEB, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, MC, TiempoInicio)

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,*) "   Diagonal element of the covariance matrix: ", MC(1,1)
     WRITE(*,*)
     WRITE(*,100) "Covariance matrix already computed"
  END IF


  !Libera memoria
  DEALLOCATE(YQER)
  DEALLOCATE(YQEI)
  DEALLOCATE(YQBR)
  DEALLOCATE(YQBI)
  DEALLOCATE(YUER)
  DEALLOCATE(YUEI)
  DEALLOCATE(YUBR)
  DEALLOCATE(YUBI)

  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NPixelesProcesador, NColumnas, -8)

  IF(iam==0) THEN
     WRITE(*,110) 
     WRITE(*,115) 
     WRITE(*,100) "Step 2. Inverting the covariance matrix ", Time()-TiempoInicio
  END IF


  Fallo = 0
  CALL InvierteMC(ICTXT, DimMC, NB, NFMC, NCMC, MC, ControlInversaMC, Fallo, TiempoInicio)

  IF(Fallo==1) THEN

     !Matriz de covarianza no regular

     IF(iam==0) THEN
        WRITE(*,*) 
        WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        WRITE(*,*) "The covariance matrix is singular"

        TiempoFinal = Time()

        WRITE(*,*)
        WRITE(*,*) "Elapsed time  (s): ", TiempoFinal - TiempoInicio
        WRITE(*,*) "Ended - Covariance matrix failure: ", CTime(Time())
        WRITE(*,*) "***************************************************************************"
     END IF

  ELSE

     !Matriz de covarianza regular         

     IF(iam==0) THEN
        WRITE(*,*)
        WRITE(*,100) "Covariance matrix inverted ", Time()-TiempoInicio
     END IF

     ! IF(iam==0) THEN
     !    WRITE(*,*) "IMC: ", MC(1,1), MC(1,2)
     ! END IF

     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, 0)

     !Fin de la parte en la que lmax esta determinado por el limite al calcular la matriz de covarianza

     IF(iam==0) THEN
        WRITE(*,110)
        WRITE(*,115)
        WRITE(*,100) "Step 3. Computing coupled power in harmonics space"
        WRITE(*,*)
        WRITE(*,100) "  Computing the blocks of the spherical harmonics matrix "
     END IF


    lmax = Dlmax
    DEALLOCATE(BeamPWDl)
    ALLOCATE(BeamPWDl(0:lmax))
    BeamPWDl = 1

    CALL blacs_barrier(ICTXT, "ALL") 
    CALL BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)

     TotalArmonicos = -3 + 2 * lmax + lmax*lmax

     !Los bloques tienen TotalArmonicos columnas, no 2*TotalArmonicos
     NColumnas = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))


     ALLOCATE(YQER(NPixelesProcesador, NColumnas))
     ALLOCATE(YQEI(NPixelesProcesador, NColumnas))
     ALLOCATE(YQBR(NPixelesProcesador, NColumnas))
     ALLOCATE(YQBI(NPixelesProcesador, NColumnas))
     ALLOCATE(YUER(NPixelesProcesador, NColumnas))
     ALLOCATE(YUEI(NPixelesProcesador, NColumnas))
     ALLOCATE(YUBR(NPixelesProcesador, NColumnas))
     ALLOCATE(YUBI(NPixelesProcesador, NColumnas))

     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NPixelesProcesador, NColumnas, 8)

     YQER=0
     YQEI=0
     YQBR=0
     YQBI=0
     YUER=0
     YUEI=0
     YUBR=0
     YUBI=0



     CALL blacs_barrier(ICTXT, "ALL")

     CALL CalculaBloquesMatrizArmonicos(ICTXT, NB, NPixelesProcesador, NColumnas, npix, lmax, z, a, BeamPWDl,&
          & YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)       

     CALL blacs_barrier(ICTXT, "ALL")

     IF(iam==0) THEN
        WRITE(*,100) "  Blocks of the spherical harmonics matrix done", Time()-TiempoInicio
     END IF


     !Pasa los bloques de la matriz de armonicos a formato de una sola matriz completa

     !Dimesiones por procesador
     NFMX = NFMC
     NCMX = MAX(1,NUMROC(2*TotalArmonicos, NB, MYCOL, 0, NPCOL))

     CALL descinit(DescMX, DimMC, 2*TotalArmonicos, NB, NB, 0, 0, ICTXT, NFMX, INFO)
     CALL descinit(DescBloque, npix, TotalArmonicos, NB, NB, 0, 0, ICTXT, NPixelesProcesador, INFO)

     !Parte real
     ALLOCATE(MXReal(NFMX,NCMX))
     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, 1)

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YQER, 1, 1, DescBloque, MXReal, 1, 1, DescMX, ictxt)
     DEALLOCATE(YQER)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YQBR, 1, 1, DescBloque, MXReal, 1, 1+TotalArmonicos, DescMX, ictxt)
     DEALLOCATE(YQBR)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YUER, 1, 1, DescBloque, MXReal, 1+npix, 1, DescMX, ictxt)
     DEALLOCATE(YUER)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YUBR, 1, 1, DescBloque, MXReal, 1+npix, 1+TotalArmonicos, DescMX, ictxt)
     DEALLOCATE(YUBR)
     CALL SumaMemoria(ictxt, iam, 0, NPixelesProcesador, NColumnas, -4)

     !Parte imaginaria
     ALLOCATE(MXImaginaria(NFMX,NCMX))
     CALL SumaMemoria(ictxt, iam, 0, NFMX, NCMX, 1)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YQEI, 1, 1, DescBloque, MXImaginaria, 1, 1, DescMX, ictxt)
     DEALLOCATE(YQEI)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YQBI, 1, 1, DescBloque, MXImaginaria, 1, 1+TotalArmonicos, DescMX, ictxt)
     DEALLOCATE(YQBI)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YUEI, 1, 1, DescBloque, MXImaginaria, 1+npix, 1, DescMX, ictxt)
     DEALLOCATE(YUEI)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(npix, TotalArmonicos, YUBI, 1, 1, DescBloque, MXImaginaria, 1+npix, 1+TotalArmonicos, DescMX, ictxt)
     DEALLOCATE(YUBI)

     IF(iam==0) THEN
        WRITE(*,100) "  Spherical harmonics matrix done", Time()-TiempoInicio
     END IF

     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NPixelesProcesador, NColumnas, -4)

     !************************************************************************************
     !Carga mapas y calcula Yl
     !Diferencia caso auto de cross-correlacion

     NCMapas = MAX(1,NUMROC(NMapas, NB, MYCOL, 0, NPCOL))

     IF(TipoAnalisis==0) THEN

        IF(iam==0) THEN
           WRITE(*,*) 
           WRITE(*,100) "  Loading maps "
        END IF

        ALLOCATE(Mapas(NFMC, NCMapas))
        Mapas = 0
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        CALL  CargaMapasFits2(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileMaps, TipoFileMaps, Mapas)

        IF(iam==0) THEN
           WRITE(*,100) "  Computing C^-1 Maps ", Time()-TiempoInicio
        END IF

        CALL blacs_barrier(ICTXT, 'All')
        CALL CalculaIMCMapas(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, MC, Mapas)

        CALL blacs_barrier(ICTXT, 'ALL')
        CALL CalculaYl(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas, MXReal,&
             &MXImaginaria, Mapas)
        DEALLOCATE(Mapas)

        IF(iam==0) WRITE(*,100) " Coupled power of auto-correlation already computed ", Time()-TiempoInicio

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, -1)

     ELSE

        IF(iam==0) THEN
           WRITE(*,*) 
           WRITE(*,100) "Loading maps from two files "
        END IF

        ALLOCATE(Mapas(NFMC, NCMapas))
        Mapas = 0
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        !call CargaMapas(ICTXT, DimMC, NMapas, NB, NFMC, NCMapas, Mapas, Lugar, 0)
        CALL  CargaMapasFits2(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileMaps, TipoFileMaps, Mapas)


        ALLOCATE(MapasCross(NFMC, NCMapas))
        MapasCross = 0
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        !call CargaMapas(ICTXT, DimMC, NMapas, NB, NFMC, NCMapas, MapasCross, Lugar, 1)
        CALL  CargaMapasFits2(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileMapsCross, TipoFileMaps, MapasCross)

        CALL blacs_barrier(ICTXT, 'All')


        IF(iam==0) THEN
           !write(*,*) 
           WRITE(*,100) "Computing C^-1 Maps ", Time()-TiempoInicio
        END IF

        CALL CalculaIMCMapasMapasCross(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, MC, Mapas, MapasCross)

        CALL blacs_barrier(ICTXT, 'All')

!!$        IF(iam==0) THEN
!!$           WRITE(*,*) 
!!$           WRITE(*,100) "Computing Yl-Cross ", Time()-TiempoInicio
!!$        END IF

        CALL CalculaYlParesMapas(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas,&
             & MXReal, MXImaginaria, Mapas, MapasCross)

        CALL blacs_barrier(ICTXT, 'All')

        DEALLOCATE(Mapas)
        DEALLOCATE(MapasCross)

        IF(iam==0) WRITE(*,100) " Coupled power cross-correlation already computed ", Time()-TiempoInicio

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, -2)


     END IF !Auto o cross-crorrelation

     !Fin carga mapas y calcula Yl
     !************************************************************************************

     IF(ComputeFisherMatrix==1) THEN

        IF(iam==0) THEN
           WRITE(*,110)
           WRITE(*,115)
           WRITE(*,100) "Step 4. Computing C^-1 Y"
           WRITE(*,*)  
           WRITE(*,100) "  Real part ", Time()-TiempoInicio
        END IF

        !Calcula C^{-1} MX
        !Dimesiones: Mismas que MX
        ALLOCATE(IMCMXR(NFMX, NCMX))
        IMCMXR = 0
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, 1)

        CALL descinit(DescMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)

        !Producto de la parte real
        CALL blacs_barrier(ICTXT, 'ALL')
        CALL pdsymm("L", "U", DimMC, 2*TotalArmonicos, 1.0D0, MC, 1, 1, DescMC, MXReal, 1, 1, DescMX, 0.0D0, IMCMXR, 1, 1, DescMX)

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "  Imaginary part ", Time()-TiempoInicio
        END IF

        DEALLOCATE(MXReal)
        ALLOCATE(IMCMXI(NFMX, NCMX))
        IMCMXI = 0

        !Producto de la parte imaginaria
        CALL blacs_barrier(ICTXT, 'ALL')
        CALL pdsymm("L", "U", DimMC, 2*TotalArmonicos, 1.0D0, MC, 1, 1, DescMC, MXImaginaria, 1, 1, DescMX, 0.0D0, IMCMXI, 1, 1, DescMX)

        DEALLOCATE(MC)
        CALL SumaMemoria(ictxt, iam, 0, NFMC, NCMC, -1)            
        DEALLOCATE(MXImaginaria)



        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "Blocks of C^-1 Y already computed", Time()-TiempoInicio
        END IF

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, -1)

        IF(iam==0) THEN
           WRITE(*,110) 
           WRITE(*,115) 
           WRITE(*,100) "Step 5. Computing noise bias ", Time()-TiempoInicio
        END IF

        !Calcula Bl
        CALL blacs_barrier(ICTXT, 'ALL')
        CALL CalculaBl(ICTXT, Lugar, DimMC, lmax, NB, NFMX, NCMX, MapaRuido2, IMCMXR, IMCMXI)
        CALL blacs_barrier(ICTXT, 'ALL')

        !*********************************************************************************
        !Calcula la Matriz de Fisher

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "Noise bias already computed"
           WRITE(*,*)
           WRITE(*,110) 
           WRITE(*,115) 
           WRITE(*,100) "Step 6. Computing the Fisher matrix"
           WRITE(*,*)
           WRITE(*,100) " Moving blocks of C^-1 Y to complex form ", Time()-TiempoInicio
        END IF

        !Prepara la multiplicacion MX^{\dag} IMC MX

        !Pasa a forma compleja IMC.MX
        ALLOCATE(IMCMX(NFMX, NCMX))
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, 2)
        IMCMX = dcmplx(IMCMXR, IMCMXI)
        DEALLOCATE(IMCMXR)
        DEALLOCATE(IMCMXI)
        CALL SumaMemoria(ictxt, iam, 0, NFMX, NCMX, -2)

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) " Blocks of C^-1 Y in complex form ", Time()-TiempoInicio
           WRITE(*,*)
           WRITE(*,100) " Computing blocks of the harmonic matrix in complex form"
           WRITE(*,*)
        END IF


        !Calcula los bloques de la matriz de armonicos
        ALLOCATE(MXC(NFMX, NCMX))
        CALL SumaMemoria(ictxt, iam, 0, NFMX, NCMX, 2)
        MXC = 0

        CALL blacs_barrier(ICTXT, 'ALL')
        CALL CalculaMatrizArmonicosCompleja(ICTXT, NB, NPixelesProcesador, NFMX, NCMX, npix, lmax, z, a, BeamPWDl,MXC) 
        CALL blacs_barrier(ICTXT, 'ALL')

        DEALLOCATE(a)
        DEALLOCATE(z)

        IF(iam==0) THEN
           WRITE(*,*) "  Multiplications Y^H C^-1 Y"
           WRITE(*,*)
           WRITE(*,*) "   EE ", Time()-TiempoInicio, "s"
        END IF

        NFB = MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))
        NCB = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))
        Alpha = dcmplx(1.0, 0.0)
        Beta  = dcmplx(0.0,0.0)

        CALL DESCINIT(DescBloque, TotalArmonicos, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFB, INFO)

        ALLOCATE(BloqueEE(NFB,NCB))
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFB, NCB, 2)
        BloqueEE = DCMPLX(0D0,0D0)
        CALL blacs_barrier( ICTXT, "A" )
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, DimMC, Alpha, MXC, 1, 1,&
             & DescMX, IMCMX, 1, 1, DescMX, Beta, BloqueEE, 1, 1, DescBloque)
        CALL blacs_barrier( ICTXT, "A" )

        IF (IAM==0) THEN
           WRITE(*,*) "   BB ", Time()-TiempoInicio, "s"
        END IF

        ALLOCATE(BloqueBB(NFB,NCB))
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFB, NCB, 2)
        BloqueBB = DCMPLX(0D0,0D0)

        CALL blacs_barrier( ICTXT, "A" )
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, DimMC, Alpha, MXC, 1, 1+TotalArmonicos,&
             & DescMX, IMCMX, 1, 1+TotalArmonicos, DescMX, Beta, BloqueBB, 1, 1, DescBloque)
        CALL blacs_barrier( ICTXT, "A" )

        IF (IAM==0) THEN
           WRITE(*,*) "   EB ", Time()-TiempoInicio, "s"
        END IF

        !EB
        ALLOCATE(BloqueEB(NFB, NCB))
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFB, NCB, 2)
        BloqueEB = DCMPLX(0D0,0D0)
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, DimMC, Alpha, MXC, 1, 1,&
             & DescMX, IMCMX, 1, 1+TotalArmonicos, DescMX, Beta, BloqueEB, 1, 1, DescBloque)
        CALL blacs_barrier( ICTXT, "A" )

        IF (IAM==0) THEN
           WRITE(*,*) 
           WRITE(*,*) " Computing blocks of the Fisher matrix"
        END IF

        DEALLOCATE(MXC)
        DEALLOCATE(IMCMX)
        !Son matrices con números complejos
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, -4)


        !Va a necesitar otro bloque. Apunto aqui la reserva de esa memoria
        CALL SumaMemoria(ictxt, iam, 0, NFB, NCB, 1)
        CALL CalculaMatrizFisher(ICTXT, lugar, lmax, NB, NFB, NCB, BloqueEE, BloqueBB, BloqueEB)

        DEALLOCATE(BloqueEE)
        DEALLOCATE(BloqueBB)
        DEALLOCATE(BloqueEB)
        !Apunta la memoria que ha liberado dentro de la llamada a la funcion
        CALL SumaMemoria(ictxt, iam, 0, NFB, NCB, -4)


        !Solo le quedan calculos al procesador 0
        !Los demas salen

        CALL blacs_barrier(ICTXT,'All')
        IF(iam==0) THEN


           IF(Binned==0) THEN

              WRITE(*,110)
              WRITE(*,*) "Computing the power spectrum"

              CALL CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)
           ELSE   
              WRITE(*,110)
              WRITE(*,*) "Computing the binned power spectrum"

              CALL CalculaBineado(NSide, npix, lmax, RuidoQQ, DlEE, DlBB, DlEB, BeamPWDl,&
                   & TypeOfGrouping, TypeOfBinCenter, QuitarSesgo, NMapas, Lugar, FileBinLimits)

           END IF !Binear o no binear


           TiempoFinal = Time()

           WRITE(*,*)
           WRITE(*,*) "Elapsed time: ", TiempoFinal - TiempoInicio
           WRITE(*,*) "End: ", CTime(Time())
           WRITE(*,110) 

        END IF !Estimar espectro. Solo procesador 0

     ELSE
        !No calcula la matriz de Fisher
        IF(iam==0) THEN
           WRITE(*,110)
           WRITE(*,*) "The program does not compute the Fisher matrix"
           WRITE(*,*) "Elapsed time: ", Time() - TiempoInicio
           WRITE(*,*) "End: ", CTime(Time())
           WRITE(*,110) 
        END IF

     END IF !Calcular la matriz de Fihser y todo lo demas

  END IF !Matriz de covarianza regular o singular

  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)

END PROGRAM ECLIPSE_EB



!******************************************************************************************************


SUBROUTINE CargaMultipolos(ICTXT, iam, Lugar, FileData, LMax, DlEE, DlBB, DlEB)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, IAM, lmax
  CHARACTER(len=100), INTENT(in) :: Lugar, FileData
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(OUT) :: DlEE, DlBB, DlEB

  INTEGER :: l

  REAL(Kind=8), DIMENSION(7) :: Temp

  EXTERNAL :: blacs_barrier, dgebs2d, dgebr2d

  IF(IAM==0) THEN

     OPEN(Unit=20, File=TRIM(Lugar)//"/"//TRIM(FileData), Action="read")
     DO l=0,lmax
        READ (20,*) temp
        DlEE(l) = temp(3)
        DlBB(l) = temp(4)
        DlEB(l) = temp(7)
     END DO
     CLOSE(20)

  END IF

  CALL blacs_barrier(ICTXT, "ALL")


  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlEE, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlEE, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlBB, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlBB, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlEB, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlEB, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE CargaMultipolos


SUBROUTINE CalculaBloquesMatrizArmonicos(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
  REAL(Kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
  REAL(Kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI

  !Polinomios y bloques
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: X1,X2
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: Alpha1, Alpha2, Beta1, Beta2, Senos, Cosenos
  REAL(Kind=8) :: Factor1, Factor2

  !Descriptores
  INTEGER, DIMENSION(9) :: DescYl, DescY
  !Tamaños
  INTEGER :: TotalArmonicos


  !Bloques intermedios
  REAL(Kind=8), ALLOCATABLE,  DIMENSION(:,:) :: XQER, XQEI, XUER, XUEI, XQBR, XQBI, XUBR, XUBI

  !Contadores
  INTEGER :: l, m, indice, Columna, base

  !Grid
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pdgemr2d

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  !Reserva memoria para los polinomios de Legendre
  ALLOCATE( P1Base1(NF, 0:lmax))
  ALLOCATE( P1Base2(NF, 0:lmax))
  ALLOCATE( Pl1(NF, 0:lmax))
  P1Base1 = 0D0
  P1Base2 = 0D0
  Pl1 = 0D0


  !Reserva memoria para X1, X2
  ALLOCATE( X1(NF, 0:lmax))
  ALLOCATE( X2(NF, 0:lmax))
  ALLOCATE( Alpha1(NF))
  ALLOCATE( Alpha2(NF))
  ALLOCATE( Beta1(NF))
  ALLOCATE( Beta2(NF))
  X1 = 0D0
  X2 = 0D0
  Alpha1 = 0D0
  Alpha2 = 0D0
  Beta1 = 0D0
  Beta2 = 0D0


  !Memoria para funciones trigonometricas
  ALLOCATE( Cosenos(NF))
  ALLOCATE( Senos(NF))
  Cosenos = 0
  Senos = 0

  !Reserva memoria para los armonicos esfericos

  !Para almacenar y pasar la parte de los armonicos de la polarizacion
  ALLOCATE(XQER(NF, 2*lmax+1))
  ALLOCATE(XQBR(NF, 2*lmax+1))
  ALLOCATE(XUER(NF, 2*lmax+1))
  ALLOCATE(XUBR(NF, 2*lmax+1))
  ALLOCATE(XQEI(NF, 2*lmax+1))
  ALLOCATE(XQBI(NF, 2*lmax+1))
  ALLOCATE(XUEI(NF, 2*lmax+1))
  ALLOCATE(XUBI(NF, 2*lmax+1))

  XQER = 0.0D0
  XQBR = 0.0D0
  XUER = 0.0D0
  XUBR = 0.0D0
  XQEI = 0.0D0
  XQBI = 0.0D0
  XUEI = 0.0D0
  XUBI = 0.0D0



  !************************************************

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  CALL DESCINIT(DescYl, npix, NPCol*(2*lmax+1), NB, 2*lmax+1, 0, 0, ICTXT, NF, INFO)
  CALL DESCINIT(DescY,  npix, TotalArmonicos, NB, NB, 0, 0, ICTXT, NF, INFO)

  !Prepara los valores iniciales de los elementos de la iteracion
  IF (MYCOL==0) THEN
     P1Base2(:,0) = 1.D0/SQRT(4D0 * PI)
     P1Base1(:,0) = SQRT(3D0/(4D0 * PI)) * z
     P1Base1(:,1) = -(1D0/2D0) * SQRT(3D0/(2*PI)) * SQRT(1-z**2)
  END IF

  !Calcula armonicos

  !IMPORTANTE: Las variables PlX NO contienen polinomios de Legendre.
  !            Son los valores absolutos de los armonicos esfericos.

  !Nota:       P1Base2 corresponde a l-2
  !            P1Base1 corresponde a l-1

  YQER = 0.0D0
  YQEI = 0.0D0
  YQBR = 0.0D0
  YQBI = 0.0D0
  YUER = 0.0D0
  YUEI = 0.0D0
  YUBR = 0.0D0
  YUBI = 0.0D0

  DO l=2,lMax

     IF(MYCOL==0) THEN
        !Avanza al siguiente polinomio
        Bucle_m: DO m=0,l-2
           Pl1(:,m) = z(:) * SQRT((4D0 * l**2-1D0)/(l**2-m**2)) * P1Base1(:,m) - &
                &SQRT(((1.0D0+2.0D0*l)*(l-m-1.0D0)*(l+m-1.0D0) * 1D0)/((2.0D0*l-3.0D0)*(l-m)*(l+m)*1D0)) * P1Base2(:,m)

           !Calula X1 y X2
           Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
           Beta1  = 2.0*((l+m)*z/(1.0-z**2))
           Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z
           Beta2  = 2.0*m/(1.0-z**2)*(l+m)

           Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

           !Introduce Beam, Pixel y calculo en Dl
           Factor1 = BeamPWDl(l) * Factor1 

           Factor2 = Factor1 * SQRT(1.0*(l-m)*(2.0D0*l+1.0D0)/((l+m)*(2.0D0*l-1.0D0)))

           X1(:,m) = Alpha1 * Factor1 * Pl1(:,m) + Beta1 * Factor2 * P1Base1(:,m)
           X2(:,m) = Alpha2 * Factor1 * Pl1(:,m) + Beta2 * Factor2 * P1Base1(:,m)

        END DO Bucle_m

        !m = l-1
        Pl1(:,l-1) = z(:) * SQRT(1D0*(1+2*l)) * P1Base1(:,l-1)

        !Calula X1 y X2
        m=l-1
        Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
        Beta1  = 2.0*((l+m)*z/(1.0-z**2))
        Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z
        Beta2  = 2.0*m/(1.0-z**2)*(l+m)

        Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

        !Introduce Beam, Pixel y calculo en Dl
        Factor1 = BeamPWDl(l) * Factor1 

        Factor2 = Factor1 * SQRT(1.0*(l-m)*(2.0D0*l+1.0D0)/((l+m)*(2.0D0*l-1.0D0)))

        X1(:,m) = Alpha1 * Factor1 * Pl1(:,m) + Beta1 * Factor2 * P1Base1(:,m)
        X2(:,m) = Alpha2 * Factor1 * Pl1(:,m) + Beta2 * Factor2 * P1Base1(:,m)

        !m = l
        Pl1(:,l) = - SQRT(((1 + 2*l)*(1 - z**2))/(2D0*l)) * P1Base1(:,l-1)

        !Calula X1 y X2
        m = l
        Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
        Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z

        Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

        !Introduce Beam, Pixel y calculo en Dl
        Factor1 = BeamPWDl(l) * Factor1 

        X1(:,m) = Alpha1 * Factor1 * Pl1(:,m)
        X2(:,m) = Alpha2 * Factor1 * Pl1(:,m)
        !Ya tiene el siguiente polinomio

        base = l+1
        indice = base
        Bucle_m2: DO m=0,l

           Cosenos = COS(1.0D0 * m * Phi)
           Senos   = SIN(1.0D0 * m * Phi)

           XQER(:,indice) = -X1(:,m) * Cosenos
           XQER(:,base-m) = (-1.0)**m * XQER(:,indice)

           XQEI(:,indice) = -X1(:,m) * Senos
           XQEI(:,base-m) = -(-1.0)**m * XQEI(:,indice)

           XQBR(:,indice) = X2(:,m) * Senos
           XQBR(:,base-m) = (-1.0)**m * XQBR(:,indice)

           XQBI(:,indice) = -X2(:,m) * Cosenos
           XQBI(:,base-m) = -(-1.0)**m * XQBI(:,indice)

           XUER(:,indice) = -X2(:,m) * Senos
           XUER(:,base-m) = (-1.0)**m * XUER(:,indice)

           XUEI(:,indice) = X2(:,m) * Cosenos
           XUEI(:,base-m) = -(-1.0)**m * XUEI(:,indice)

           XUBR(:,indice) = -X1(:,m) * Cosenos
           XUBR(:,base-m) = (-1.0)**m * XUBR(:,indice)

           XUBI(:,indice) = -X1(:,m) * Senos
           XUBI(:,base-m) = -(-1.0)**m * XUBI(:,indice)

           indice = indice + 1

        END DO Bucle_m2

     END IF

     CALL blacs_barrier( ICTXT, "A" )

     Columna = -3 + l*l

     CALL pdgemr2d(NPix, 2*l+1, XQER, 1, 1, DescYl, YQER, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPix, 2*l+1, XQEI, 1, 1, DescYl, YQEI, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pdgemr2d(NPix, 2*l+1, XQBR, 1, 1, DescYl, YQBR, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPix, 2*l+1, XQBI, 1, 1, DescYl, YQBI, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pdgemr2d(NPix, 2*l+1, XUER, 1, 1, DescYl, YUER, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPix, 2*l+1, XUEI, 1, 1, DescYl, YUEI, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pdgemr2d(NPix, 2*l+1, XUBR, 1, 1, DescYl, YUBR, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPix, 2*l+1, XUBI, 1, 1, DescYl, YUBI, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     IF(MYCOL==0) THEN
        !Adelanta un paso los polinomios base de la iteracion
        P1Base2(:,0:l-1) = P1Base1(:,0:l-1)
        P1Base1(:,0:l) = Pl1(:,0:l)
     END IF

  END DO

END SUBROUTINE CalculaBloquesMatrizArmonicos


SUBROUTINE CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NFBA, NCBA, DimMC, lmax, MapaRuido2Bloque,&
     & DlEE, DlBB, DlEB, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, MC, TiempoInicio)

  USE ifport

  IMPLICIT NONE

  !NFMC : Numero de filas de la matriz de covarianza en cada procesador
  !NCMC : Numero de columnas de la matriz de covarianza en cada procesador
  !NFBA : Numero de filas de cada bloque de armonicos en cada procesador
  !NCBA : Numero de columnas de cada bloque de armonicos en cada procesador
  !DimMC : Dimension global de la matriz de covarianza

  INTEGER, INTENT(in) :: ICTXT, NB, NFMC, NCMC, NFBA, NCBA, DimMC, lmax
  INTEGER(Kind=8) :: TiempoInicio
  REAL(Kind=8), DIMENSION(NFBA), INTENT(in) :: MapaRuido2Bloque
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: DlEE, DlBB, DlEB
  REAL(Kind=8), DIMENSION(NFBA,NCBA), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
  REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(out)  :: MC

  !************************+

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  INTEGER :: iam, TotalArmonicos
  INTEGER :: i, j, indice, l

  !Datos sobre los bloques QQ, UU, QU, UQ de la matriz de covarianza
  INTEGER :: NFGlobalBloqueQQ, NCGlobalBloqueQQ, NFLocalBloqueQQ, NCLocalBloqueQQ

  !Descriptores
  INTEGER, DIMENSION(9) :: DescBloqueArmonicos, DescBloqueProducto, DescMC

  !Tamaños
  INTEGER :: NUMROC, INDXL2G

  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesFBMC, IndicesCBMC, IndicesL

  INTEGER :: NPROW, NPCOL, MYROW, MYCOL, INFO

  !Para productos intermedios
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MTemp, MTemp2, MTemp3, MTemp4
  REAL(Kind=8) :: Factor

  !Componente EB
  INTEGER :: ContribucionEB

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pdgemr2d, pdgemm, pdtran

  !**********************************************************************************

  ContribucionEB = 0
  l = 2
  DO WHILE ((ContribucionEB == 0).AND.(l<=lmax))       
     IF(DlEB(l)/=0) ContribucionEB = 1
     l = l + 1      
  END DO

  !Codigo
  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  NFGlobalBloqueQQ = DimMC/2
  NCGlobalBloqueQQ = DimMC/2

  NFLocalBloqueQQ = MAX(1,NUMROC(NFGlobalBloqueQQ, NB, MYROW, 0, NPROW))
  NCLocalBloqueQQ = MAX(1,NUMROC(NCGlobalBloqueQQ, NB, MYCOL, 0, NPCOL))

  !Indices de las filas y columnas de los bloques de la matriz de covarianza
  ALLOCATE(IndicesFBMC(NFLocalBloqueQQ))
  ALLOCATE(IndicesCBMC(NCLocalBloqueQQ))

  !Valores de l que corresponden a las columnas de los bloques de la matriz de armónicos
  ALLOCATE(IndicesL(NCBA))

  DO i=1,NFLocalBloqueQQ
     IndicesFBMC(i) = INDXL2G(i, NB, MYROW, 0, NPROW)
  END DO

  DO i=1,NCLocalBloqueQQ
     IndicesCBMC(i) = INDXL2G(i, NB, MYCOL, 0, NPCOL)
  END DO

  !Recorre las columnas de los bloques de armonicos
  DO i=1,NCBA
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     Factor = 3.0D0 + 1.0D0 * j
     IndicesL(i) = FLOOR(SQRT(Factor))
  END DO

  ALLOCATE(MTemp(NFBA, NCBA))

  !Dimensiones del bloque de la matriz de covarianza
  ALLOCATE(MTemp2(NFLocalBloqueQQ, NCLocalBloqueQQ))
  ALLOCATE(MTemp3(NFLocalBloqueQQ, NCLocalBloqueQQ))
  ALLOCATE(MTemp4(NFLocalBloqueQQ, NCLocalBloqueQQ))
  MTemp = 0
  MTemp2 = 0
  MTemp3 = 0
  MTemp4 = 0


  !Descriptor del bloque de armonicos: mismas filas que bloque tipo QQ
  CALL DESCINIT(DescBloqueArmonicos, NFGlobalBloqueQQ, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFLocalBloqueQQ, INFO)

  !Descriptor del bloque producto: bloque tipo QQ
  CALL DESCINIT(DescBloqueProducto, NFGlobalBloqueQQ, NCGlobalBloqueQQ, NB, NB, 0, 0, ICTXT, NFLocalBloqueQQ, INFO)

  !Descriptor de la matriz de covarianza completa
  CALL DESCINIT(DescMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  MC = 0.0D0

  !Caso QQ ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YQER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBR(:,i)
        END IF
     END IF
  END DO


  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YQEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma las componentes real e imaginaria
  MTemp3 = MTemp2 + MTemp3

  !Segundo sumando
  !Real
  MTemp = YQBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YQBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp4, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma las componentes real e imaginaria
  MTemp4 = MTemp2 + MTemp4
  MTemp3 = MTemp3 + Mtemp4

  !Suma ruido
  DO i=1,NFLocalBloqueQQ
     indice = IndicesFBMC(i)
     DO j=1,NCLocalBloqueQQ
        IF(indice==IndicesCBMC(j)) THEN
           MTemp3(i,j) = MTemp3(i,j) + MapaRuido2Bloque(i)
        END IF
     END DO
  END DO

  !Lo envia a la matriz de covarianza
  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemr2d(NFGlobalBloqueQQ, NCGlobalBloqueQQ, MTemp3, 1, 1, DescBloqueProducto, MC, 1, 1 , DescMC, ictxt)

  IF(iam==0) WRITE(*,*) "    Block QQ done ", Time()-TiempoInicio, "s"

  CALL blacs_barrier( ICTXT, "A" )


  !Caso QU ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YUER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBR(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma las componentes real e imaginaria
  MTemp3 = MTemp2 + MTemp3

  !Segundo sumnado
  !Real
  MTemp = YUBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YQBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp4, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma las componentes real e imaginaria
  MTemp4 = MTemp2 + MTemp4
  MTemp3 = MTemp3 + Mtemp4

  !Lo envia a la matriz de covarianza
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NFGlobalBloqueQQ, NCGlobalBloqueQQ, MTemp3, 1, 1, DescBloqueProducto, MC, 1, NFGlobalBloqueQQ+1 , DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  !Traspuesta
  CALL pdtran(NFGlobalBloqueQQ, NCGlobalBloqueQQ, 1.0D0, MTemp3, 1, 1, DescBloqueProducto, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)
  !Envia al bloque inferior
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NFGlobalBloqueQQ, NCGlobalBloqueQQ, MTemp2, 1, 1, DescBloqueProducto, MC, NFGlobalBloqueQQ+1, 1 , DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF(iam==0) WRITE(*,*) "    Block QU done ", Time()-TiempoInicio, "s"

  !Caso UU ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YUER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBR(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YUER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YUEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma las componentes real e imaginaria
  MTemp3 = MTemp2 + MTemp3

  !Segundo sumnado
  !Real
  MTemp = YUBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YUBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescBloqueProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBA
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        !write(*,*) i, l, Factor
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NFGlobalBloqueQQ, NCGlobalBloqueQQ, TotalArmonicos, 1.0D0, YUBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp4, 1, 1, DescBloqueProducto)

  !Suma las componentes real e imaginaria
  MTemp4 = MTemp2 + MTemp4
  MTemp3 = MTemp3 + Mtemp4

  !Suma ruido
  DO i=1,NFLocalBloqueQQ
     indice = IndicesFBMC(i)
     DO j=1,NCLocalBloqueQQ
        IF(indice==IndicesCBMC(j)) THEN
           MTemp3(i,j) = MTemp3(i,j) + MapaRuido2Bloque(i)
        END IF
     END DO
  END DO

  !Lo envia a la matriz de covarianza
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NFGlobalBloqueQQ, NCGlobalBloqueQQ, MTemp3, 1, 1, DescBloqueProducto, MC,&
       & NFGlobalBloqueQQ+1, NCGlobalBloqueQQ+1, DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF(iam==0) WRITE(*,*) "    Block UU done ", Time()-TiempoInicio, "s"

  DEALLOCATE(IndicesL)
  DEALLOCATE(IndicesFBMC)
  DEALLOCATE(IndicesCBMC)

  DEALLOCATE(MTemp)
  DEALLOCATE(MTemp2)
  DEALLOCATE(MTemp3)
  DEALLOCATE(MTemp4)

END SUBROUTINE CalculaMatrizCovarianza


SUBROUTINE InvierteMC(ICTXT, DimMC, NB, NF, NC, MC, ControlInversa, Fallo, TiempoInicio)

  !Invierte la matriz de covarianza y comprueba la diagonal del producto IMC.MC
  !Tal como esta escrita, necesita recservar memomria para otra matriz de covarianza.
  !Si no se quiere hacer esa prueba (la diagonal del producto IMC.MC) no es necesario reservar una
  !matriz de covarianza extra.

  USE ifport
  USE ControlMemoria

  IMPLICIT NONE
  INTEGER, INTENT(in) :: ICTXT, DimMC, NB, NF, NC, ControlInversa
  REAL(Kind=8), DIMENSION(NF,NC), INTENT(inout) :: MC
  INTEGER, INTENT(out) :: Fallo
  INTEGER(Kind=8), INTENT(in) ::  TiempoInicio

  !********************************************************

  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: IMC
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesC
  INTEGER, DIMENSION(9) :: DESC, DescFilaRepartida, DescFila, DescColumnaRepartida, DescColumna
  INTEGER :: info
  REAL(Kind=8) :: ProductoDiagonal, SumaDiagonal
  INTEGER ::  IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i,j, ii, jj, Seguir
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: SumaFila, SumaColumna, TodasColumnas, TodasFilas

  EXTERNAL :: BLACS_GRIDINFO, descinit, pdpotrf, pdpotri, pdtran, blacs_barrier, DGSUM2D, pdgemr2d
  INTEGER, EXTERNAL :: INDXL2G

  Fallo = 0

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  CALL descinit(desc, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NF, INFO)

  IF(ControlInversa==0) THEN

     !Invierte matriz
     CALL pdpotrf("U", DimMC, MC, 1, 1, desc, info)
     IF(IAM==0) THEN
        WRITE(*,*)
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Cholesky factorization: ", info, Time()-TiempoInicio
     END IF

     CALL blacs_barrier(ICTXT, "All")

     IF(info/=0) THEN
        Fallo = 1
        RETURN
     END IF

     CALL blacs_barrier(ICTXT, "All")

     CALL pdpotri("U", DimMC, MC, 1, 1, desc, info)
     IF(IAM==0) THEN
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Inversion result:        ", info, Time()-TiempoInicio
     END IF

     CALL blacs_barrier(ICTXT, "All")

     IF(info/=0) THEN
        Fallo = 1
        RETURN
     END IF

  ELSE

     ALLOCATE (IMC(NF,NC))
     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NF, NC, 1)
     IMC = MC


     !Invierte matriz
     CALL pdpotrf("U", DimMC, IMC, 1, 1, desc, info)
     IF(IAM==0) THEN
        WRITE(*,*)
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Cholesky factorization: ", info, Time()-TiempoInicio
     END IF

     CALL blacs_barrier(ICTXT, "All")

     IF(info/=0) THEN
        Fallo = 1
        DEALLOCATE(IMC)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NF, NC, -2)
        RETURN
     END IF

     CALL blacs_barrier(ICTXT, "All")

     CALL pdpotri("U", DimMC, IMC, 1, 1, desc, info)
     IF(IAM==0) THEN
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Inversion result:        ", info, Time()-TiempoInicio
     END IF

     CALL blacs_barrier(ICTXT, "All")

     IF(info/=0) THEN
        Fallo = 1
        DEALLOCATE(IMC)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NF, NC, -2)
        RETURN
     END IF

     ALLOCATE(IndicesC(NC))
     DO i=1,NC
        IndicesC(i) = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     END DO

     DO i=1,NF
        seguir = 1
        ii = INDXL2G(i, NB, MYROW, 0, NPROW)
        j = 1
        DO WHILE (seguir == 1)
           jj = IndicesC(j)
           IF (jj < ii) THEN
              IMC(i,j) = 0.0D0
           ELSE IF (jj > ii) THEN
              seguir = 0
           ELSE IF (jj == ii) THEN
              IMC(i,j) = IMC(i,j)/2.0D0
              seguir = 0
           END IF
           j = j + 1
           IF(j>NC) seguir = 0
        END DO
     END DO


     CALL blacs_barrier(ICTXT, "All")

     MC = IMC*MC

     ALLOCATE(SumaColumna(NC))
     SumaColumna = 0
     SumaColumna= SUM(MC,DIM=1)

     ALLOCATE(SumaFila(NF))
     SumaFila = 0
     SumaFila= SUM(MC,DIM=2)

     CALL DGSUM2D(ICTXT,  'C', '1-tree', 1, NC, SumaColumna, 1, 0, MYCOL)
     CALL blacs_barrier(ICTXT, "All")
     CALL DGSUM2D(ICTXT,  'R', '1-tree', 1, NF, SumaFila, 1, MYROW, 0)
     CALL blacs_barrier(ICTXT, "All")

     !Arregla la diagonal de IMC
     DO i=1,NF
        seguir = 1
        ii = INDXL2G(i, NB, MYROW, 0, NPROW)
        j = 1
        DO WHILE (seguir == 1)
           jj = IndicesC(j)
           IF (jj == ii) THEN
              IMC(i,j) = 2* IMC(i,j)
              seguir = 0
           END IF
           j = j + 1
           IF(j>NC) seguir = 0
        END DO
     END DO

     !Pasa la inversa a MC para devolver la inversa
     MC = IMC
     DEALLOCATE(IMC)

     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NF, NC, -1)

     ALLOCATE(TodasColumnas(DimMC))
     ALLOCATE(TodasFilas(DimMC))

     CALL DESCINIT(DescFilaRepartida, 1, DimMC, NB, NB, 0, 0, ICTXT, 1, INFO)
     CALL DESCINIT(DescFila, 1, DimMC, NB, DimMC, 0, 0, ICTXT, 1, INFO)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(1, DimMC, SumaColumna, 1, 1, DescFilaRepartida, TodasColumnas, 1, 1, DescFila, ictxt)

     CALL DESCINIT(DescColumnaRepartida, DimMC, 1, NB, NB, 0, 0, ICTXT, NF, INFO)
     CALL DESCINIT(DescColumna, DimMC, 1, DimMC, NB, 0, 0, ICTXT, DimMC, INFO)
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(DimMC, 1, SumaFila, 1, 1, DescColumnaRepartida, TodasFilas, 1, 1, DescColumna, ictxt)

     CALL blacs_barrier(ICTXT, "ALL")

     IF(iam==0) THEN

        TodasFilas = TodasFilas + TodasColumnas

        SumaDiagonal = SUM(TodasFilas)
        ProductoDiagonal = PRODUCT(TodasFilas)

        WRITE (*,*) "   Diagonal product: ", ProductoDiagonal
        WRITE (*,*) "   Diagonal sum:     ", SumaDiagonal
        WRITE (*,*) "   Matrix size:      ", DimMC

     END IF


     DEALLOCATE(TodasFilas)
     DEALLOCATE(TodasColumnas)
     DEALLOCATE(SumaColumna)
     DEALLOCATE(SumaFila)
     DEALLOCATE(IndicesC)

  END IF

END SUBROUTINE InvierteMC



SUBROUTINE CalculaIMCMapas(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
  REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
  REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  INTEGER, DIMENSION(9) :: DESCIMC, DESCMapas, DescTodos, DescDistribuidos
  INTEGER :: info
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Producto
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: SumaParte, RecibeSumas, ListaChi2
  REAL(Kind=8) :: SumaChi2
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER, EXTERNAL :: NUMROC

  EXTERNAL :: BLACS_GRIDINFO, descinit, pdsymm, DGSUM2D, dgerv2d, dgesd2d, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  CALL descinit(DESCIMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  ALLOCATE(Producto(NFMC, NCMapas))
  Producto = 0

  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, Mapas, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)

  Mapas = Mapas * Producto
  ALLOCATE(SumaParte(1,NCMapas))
  ALLOCATE(RecibeSumas(1,NCMapas))
  RecibeSumas = 0
  SumaParte   = 0

  SumaParte(1,:) = SUM(Mapas, DIM=1)

  !Suma los datos en cada columna
  CALL blacs_barrier(ICTXT, 'ALL')
  CALL DGSUM2D(ICTXT, 'C', '1-tree', 1, NCMapas, SumaParte, 1, -1, MYCOL)

  ALLOCATE(ListaChi2(1,NMapas))
  ListaChi2 = 0 

  CALL descinit(DescTodos, 1, NMapas, 1, NMapas, 0, 0, ICTXT, 1, INFO)
  CALL descinit(DescDistribuidos, 1, NMAPAS, 1, NB, 0, 0, ICTXT, 1, INFO)

  !Pasa todos los datos a la ListaChi del primer procesador
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(1, NMapas, SumaParte, 1, 1, DescDistribuidos, ListaChi2, 1, 1, DescTodos, ictxt)
  CALL blacs_barrier(ICTXT, "ALL")

  IF((MYCOl==0).AND.(MYROW==0)) THEN
     OPEN(unit=20, file=TRIM(Lugar)//"/ChiSquare.dat", action="write")
     WRITE (20,*) ListaChi2
     CLOSE(20)
     SumaChi2 = SUM(ListaChi2)/NMapas
     WRITE(*,*) "    <x^t C^-1 x>: ", SumaChi2
     WRITE(*,*) "    Matrix size:  ", DimMC
  END IF

  !Devuelve IMC.Mapas
  Mapas = Producto

  DEALLOCATE(Producto)
  DEALLOCATE(SumaParte)
  DEALLOCATE(ListaChi2)

  CALL blacs_barrier(ICTXT, 'All')

END SUBROUTINE CalculaIMCMapas


SUBROUTINE CalculaYl(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas, MXR, MXI, Mapas)

  USE ifport
  USE ControlMemoria

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas
  CHARACTER(LEN=100), INTENT(in) :: Lugar
  REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: MXR, MXI
  REAL(Kind=8), DIMENSION(NFMX, NCMapas), INTENT(in) :: Mapas

  !************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  !Descriptores
  INTEGER, DIMENSION(9) :: DescMX, DescMapas, DescPMXTMapas, DescTemp, DescUnMapa, DescTodosMapas

  !Tamaños
  INTEGER :: TotalArmonicos, NPixelesProcesador, NFPMXTMapas, NFTemp

  !Contadores
  INTEGER :: l, i,j, salto

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  !Para calcular bl de ruido diagonal en pixeles
  INTEGER :: Caso
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2, Temp3

  !Para calcular Yl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MPYMR, MPYMI
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: ListaYl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:)   :: YlReestructurado, RecibeMapa

  INTEGER, EXTERNAL :: INDXL2G, NUMROC
  EXTERNAL :: BLACS_GRIDINFO, descinit, pdgemm, pdgemr2d, dgesd2d, dgerv2d, DGSUM2D, blacs_barrier


  !Codigo
  NPixelesProcesador = DimMC/2

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  !write(*,*) "Yl. ", iam, NPixelesProcesador, NB, TotalArmonicos

  NFPMXTMapas  = MAX(1,NUMROC(2*TotalArmonicos, NB, MYROW, 0, NPROW))

  !Multiplica Y^t por IMC.Mapas
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMX, INFO)
  CALL descinit(DescPMXTMapas, 2*TotalArmonicos, NMapas, NB, NB, 0, 0, ICTXT, NFPMXTMapas, INFO)
  CALL DESCINIT(DescMX, DimMC, 2*TotalArmonicos, NB, NB, 0, 0, ICTXT, NFMX, INFO)

  ALLOCATE(MPYMR(NFPMXTMapas, NCMapas))
  ALLOCATE(MPYMI(NFPMXTMapas, NCMapas))
  MPYMR = 0
  MPYMI = 0


  !Temp1 tiene filas como MX traspuesta y columnas como NMapas
  ALLOCATE(Temp1(NFPMXTMapas, NCMapas))
  Temp1 = 0
  !call SumaMemoria(ictxt, iam, MuestraMemoria, NFPMXTMapas, NCMapas, 3)


  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXR, 1, 1, DescMX, Mapas, 1, 1, DescMapas,&
       & 0.0D0, MPYMR, 1, 1, DescPMXTMapas)

  !*************

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXI, 1, 1, DescMX, Mapas, 1, 1, DescMapas,&
       & 0.0D0, MPYMI, 1, 1, DescPMXTMapas)

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,*) "   Maps transformed to harmonis space"
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  !Yl EE y BB
  Temp1 = MPYMR**2+MPYMI**2

  ALLOCATE(ListaYl(2:lMax,3,NCMapas))
  ListaYl = 0.0D0

  DO i=1,NFPMXTMapas
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     IF(j<= TotalArmonicos) THEN
        Caso = 1
     ELSE
        Caso = 2
        j = j-TotalArmonicos
     ENDIF
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) ListaYl(l,Caso,:) = ListaYl(l,Caso,:) + Temp1(i,:)
  END DO


  !Yl EB

  DEALLOCATE(Temp1)
  !call SumaMemoria(ictxt, iam, MuestraMemoria, NFPMXTMapas, NCMapas, -1)

  !Ahora los bloques tienen numero de filas: TotalArmonicos
  NFTemp =MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))

  ALLOCATE(Temp1(NFTemp, NCMapas))
  ALLOCATE(Temp2(NFTemp, NCMapas))
  ALLOCATE(Temp3(NFTemp, NCMapas))
  !call SumaMemoria(ictxt, iam, MuestraMemoria, NFTemp, NCMapas, 3)


  Temp1 = 0
  Temp2 = 0
  Temp3 = 0

  CALL DESCINIT(DescTemp, TotalArmonicos, NMapas, NB, NB, 0, 0, ICTXT, NFTemp, INFO)

  !Pasa la parte real de Q y de U
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMR, 1, 1, DescPMXTMapas, Temp1, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMR, 1+TotalArmonicos, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)

  Temp1 = Temp1*Temp2

  !Pasa la parte imaginaria de Q de U
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMI, 1, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMI, 1+TotalArmonicos, 1, DescPMXTMapas, Temp3, 1, 1, DescTemp, ictxt)

  Temp1 = 2.0*(Temp1 + Temp2*Temp3)

  DO i=1,NFTemp
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) ListaYl(l,3,:) = ListaYl(l,3,:) + Temp1(i,:)
  END DO


  ListaYl = ListaYl/2.0D0

  DEALLOCATE(Temp1)
  DEALLOCATE(Temp2)
  DEALLOCATE(Temp3)
  !call SumaMemoria(ictxt, iam, MuestraMemoria, NFTemp, NCMapas, -3)

  !Cada procesador tiene una parte de Yl

  !Suma Yl en las columnas
  !La suma esta lanzada en una forma que no tiene en cuenta la forma de la variable ListaYl
  CALL blacs_barrier( ICTXT, "A" )
  CALL DGSUM2D( ictxt, 'C', '1-tree', 3*(lMax-1)*NCMapas, 1, ListaYl, 3*(lMax-1)*NCMapas, -1, MYCOL)

  !Todos los procesadores de la misma columna ya tienen completos los Yl de su parte de los mapas
  !Los procesadores de la primera fila le envían los datos al procesador 0, y los guarda

  salto = lmax-1
  ALLOCATE(YlReestructurado(3*salto, NCMapas))
  DO i=1,NCMapas
     YlReestructurado(1:salto,i) = ListaYl(:,1,i)
     YlReestructurado(salto+1:2*salto,i) = ListaYl(:,2,i)
     YlReestructurado(2*salto+1:3*salto,i) = ListaYl(:,3,i)
  END DO

  ALLOCATE(RecibeMapa(3*salto,1))

  CALL descinit(DescTodosMapas, 3*salto, NMapas, 3*salto, NB, 0, 0, ICTXT, 3*salto, INFO)
  CALL descinit(DescUnMapa, 3*salto, 1, 3*salto, NB, 0, 0, ICTXT, 3*salto, INFO)


  IF (IAM==0) THEN
     OPEN(Unit=20, File=TRIM(Lugar)//"/CoupledPower.dat", Action="write")
  END IF

  DO i=1,NMapas

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(3*salto, 1, YlReestructurado, 1, i, DescTodosMapas, RecibeMapa, 1, 1, DescUnMapa, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")
     IF(iam==0) WRITE(20,*) RecibeMapa(:,1)

  END DO

  IF (IAM==0) THEN
     CLOSE(20)
  END IF

  DEALLOCATE(ListaYl)
  DEALLOCATE(YlReestructurado)
  DEALLOCATE(RecibeMapa)    

  !Fin calcula Yl
  CALL blacs_barrier( ICTXT, "A" )

END SUBROUTINE CalculaYl


SUBROUTINE CalculaBl(ICTXT, Lugar, DimMC, lmax, NB, NFMX, NCMX, MapaRuido2, IMCMXR, IMCMXI)

  USE ifport
  USE ControlMemoria

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NB, NFMX, NCMX
  CHARACTER(LEN=100), INTENT(in) :: Lugar
  REAL(Kind=8), DIMENSION(NFMX), INTENT(in) :: MapaRuido2
  REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: IMCMXR, IMCMXI

  INTEGER   :: TotalArmonicos, NFB, NCB

  !Contadores
  INTEGER :: l,i,j

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO
  INTEGER, EXTERNAL :: INDXL2G, NUMROC

  !Para calcular bl de ruido diagonal en pixeles
  INTEGER :: Caso
  INTEGER, DIMENSION(9) :: DescTemp, DescMX
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MAE, Bl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2, Temp3
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: SumasBl

  EXTERNAL :: DESCINIT, pdgemr2d, DGSUM2D, BLACS_GRIDINFO, blacs_barrier

  !*********************************************************************

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  CALL blacs_barrier( ICTXT, "A" )
  ALLOCATE(MAE(NFMX,NCMX))
  MAE = 0
  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMX, NCMX, 1)

  !Calcula (IMC.YP)^h.(IMC.YP)
  MAE = IMCMXR**2+IMCMXI**2

  !Multiplica las filas por el ruido en los píxeles
  DO i=1,NCMX
     MAE(:,i) = MAE(:,i)*MapaRuido2(:)
  END DO

  ALLOCATE (SumasBl(NCMX))

  !Suma las filas
  SumasBl = SUM(MAE,Dim=1)

  ALLOCATE(Bl(3,2:lmax))
  Bl = 0

  !Coloca las sumas repartidas por los procesadores en
  !los Bl que les corresponden
  DO i=1,NCMX
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     IF(j<= TotalArmonicos) THEN
        Caso = 1
     ELSE
        Caso = 2
        j = j-TotalArmonicos
     ENDIF
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) Bl(Caso,l) = Bl(Caso, l) + SumasBl(i)
  END DO

  DEALLOCATE(MAE)
  CALL SumaMemoria(ictxt, iam, 0, NFMX, NCMX, -1)

  IF(iam==0) THEN
     WRITE(*,*) "    EE & BB"
  END IF

  !Terminos BlEB
  !Los  bloques ahora tienen dimensiones DimMC x TotalArmonicos
  NFB = NFMX
  NCB = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))

  ALLOCATE(Temp1(NFB, NCB))
  ALLOCATE(Temp2(NFB, NCB))
  ALLOCATE(Temp3(NFB, NCB))
  Temp1 = 0
  Temp2 = 0
  Temp3 = 0
  CALL SumaMemoria(ictxt, iam, 0, NFB, NCB, 3)

  CALL DESCINIT(DescTemp, DimMC, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFB, INFO)
  CALL DESCINIT(DescMX, DimMC, 2*TotalArmonicos, NB, NB, 0, 0, ICTXT, NFMX, INFO)

  !Pasa la parte real de la primera mitad de los armonicos a Temp1
  CALL pdgemr2d(DimMC, TotalArmonicos, IMCMXR, 1, 1, DescMX, Temp1, 1, 1 , DescTemp, ictxt)
  !Pasa la parte real de la segunda mitad de los armonicos a Temp2
  CALL pdgemr2d(DimMC, TotalArmonicos, IMCMXR, 1, 1+TotalArmonicos, DescMX, Temp2, 1, 1 , DescTemp, ictxt)

  !Multiplica las partes reales
  Temp3 = Temp1*Temp2

  !Pasa la parte imaginaria de la primera mitad de los armonicos a Temp1
  CALL pdgemr2d(DimMC, TotalArmonicos, IMCMXI, 1, 1, DescMX, Temp1, 1, 1 , DescTemp, ictxt)
  !Pasa la parte imaginaria de la segunda mitad de los armonicos a Temp2
  CALL pdgemr2d(DimMC, TotalArmonicos, IMCMXI, 1, 1+TotalArmonicos, DescMX, Temp2, 1, 1 , DescTemp, ictxt)

  !Multiplica las partes imaginarias
  Temp1 = Temp1*Temp2

  !Suma las multiplicaciones
  Temp3 = Temp3+Temp1

  !Multiplica por dos
  Temp3 = 2.0D0*Temp3

  DO i=1,NCB
     Temp3(:,i) = Temp3(:,i)*MapaRuido2(:)
  END DO

  DEALLOCATE(Temp1)
  DEALLOCATE(Temp2)

  !Elimna el vector anterior para sumas 2*TotalArmonicos resultados
  DEALLOCATE(SumasBl)

  !Reserva memoria para vector de TotalArmonicos resultados
  ALLOCATE (SumasBl(NCB))

  !Suma las filas
  SumasBl = SUM(Temp3,Dim=1)


  !Multiplica las sumas por el valor del ruido diagonal
  !Esto solo vale para ruido isotropo
  !SumasBl = SumasBl*RuidoQQ**2

  !Suma las columnas
  DO i=1,NCB
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) Bl(3,l) = Bl(3, l) + SumasBl(i)
  END DO

  IF(iam==0) THEN
     WRITE(*,*) "    EB"
  END IF

  DEALLOCATE(Temp3)
  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFB, NCB, -3)

  !Suma todos los Bl repartidos por los procesadores
  CALL DGSUM2D( ictxt, 'All', '1-tree', 3, lMax-1, Bl, 3, -1, -1)

  Bl = Bl/2.0D0

  !Muestra resultado
  IF (iam==0) THEN        
     OPEN(Unit=20, File=TRIM(Lugar)//"/NoiseBias.dat", Action="write")
     DO Caso=1,3
        DO l=2,lMax
           WRITE (20,*) Bl(Caso,l)
        END DO
     END DO
     CLOSE(20)
  END IF


END SUBROUTINE CalculaBl

SUBROUTINE CalculaMatrizArmonicosCompleja(ICTXT, NB, NF, NFMX, NCMX, npix, lmax, z, Phi, BeamPWDl, MXC)

  !NF : Numero de pixeles por procesador

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NF, NFMX, NCMX, npix, lmax
  REAL(Kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
  COMPLEX(Kind=8), DIMENSION(NFMX,NCMX), INTENT(OUT) :: MXC

  !Polinomios y bloques
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: X1,X2
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: Alpha1, Alpha2, Beta1, Beta2, Senos, Cosenos
  REAL(Kind=8) :: Factor1, Factor2

  !Descriptores
  INTEGER, DIMENSION(9) :: DescYl, DescMX
  !Tamaños
  INTEGER :: TotalArmonicos


  !Bloques intermedios
  COMPLEX(Kind=8), ALLOCATABLE,  DIMENSION(:,:) :: XQE, XQB, XUE, XUB

  !Contadores
  INTEGER :: l, m, indice, Columna, base

  !Grid
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pzgemr2d

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  !Reserva memoria para los polinomios de Legendre
  ALLOCATE( P1Base1(NF, 0:lmax))
  ALLOCATE( P1Base2(NF, 0:lmax))
  ALLOCATE( Pl1(NF, 0:lmax))
  P1Base1 = 0D0
  P1Base2 = 0D0
  Pl1 = 0D0


  !Reserva memoria para X1, X2
  ALLOCATE( X1(NF, 0:lmax))
  ALLOCATE( X2(NF, 0:lmax))
  ALLOCATE( Alpha1(NF))
  ALLOCATE( Alpha2(NF))
  ALLOCATE( Beta1(NF))
  ALLOCATE( Beta2(NF))
  X1 = 0D0
  X2 = 0D0
  Alpha1 = 0D0
  Alpha2 = 0D0
  Beta1 = 0D0
  Beta2 = 0D0


  !Memoria para funciones trigonometricas
  ALLOCATE(Cosenos(NF))
  ALLOCATE(Senos(NF))
  Cosenos = 0
  Senos = 0

  !Reserva memoria para los armonicos esfericos

  !Para almacenar y pasar la parte de los armonicos de la polarizacion
  ALLOCATE(XQE(NF, 2*lmax+1))
  ALLOCATE(XUE(NF, 2*lmax+1))
  ALLOCATE(XQB(NF, 2*lmax+1))
  ALLOCATE(XUB(NF, 2*lmax+1))

  XQE = 0.0D0
  XUE = 0.0D0
  XQB = 0.0D0
  XUB = 0.0D0

  !************************************************

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  CALL DESCINIT(DescYl, npix, NPCol*(2*lmax+1), NB, 2*lmax+1, 0, 0, ICTXT, NF, INFO)
  CALL DESCINIT(DescMX, 2*npix, 2*TotalArmonicos, NB, NB, 0, 0, ICTXT, NFMX, INFO)

  !Prepara los valores iniciales de los elementos de la iteracion
  IF (MYCOL==0) THEN
     P1Base2(:,0) = 1.D0/SQRT(4D0 * PI)
     P1Base1(:,0) = SQRT(3D0/(4D0 * PI)) * z
     P1Base1(:,1) = -(1D0/2D0) * SQRT(3D0/(2*PI)) * SQRT(1-z**2)
  END IF

  !Calcula armonicos

  !IMPORTANTE: Las variables PlX NO contienen polinomios de Legendre.
  !            Son los valores absolutos de los armonicos esfericos.

  !Nota:       P1Base2 corresponde a l-2
  !            P1Base1 corresponde a l-1

  MXC = 0

  DO l=2,lMax

     IF(MYCOL==0) THEN

        !Avanza al siguiente polinomio
        Bucle_m: DO m=0,l-2
           Pl1(:,m) = z(:) * SQRT((4D0 * l**2-1D0)/(l**2-m**2)) * P1Base1(:,m) - &
                &SQRT(((1.0D0+2.0D0*l)*(l-m-1.0D0)*(l+m-1.0D0) * 1D0)/((2.0D0*l-3.0D0)*(l-m)*(l+m)*1D0)) * P1Base2(:,m)

           !Calula X1 y X2
           Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
           Beta1  = 2.0*((l+m)*z/(1.0-z**2))
           Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z
           Beta2  = 2.0*m/(1.0-z**2)*(l+m)

           Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

           !Introduce Beam, Pixel Window y calculo en Dl
           Factor1 = BeamPWDl(l) * Factor1

           Factor2 = Factor1 * SQRT(1.0*(l-m)*(2.0D0*l+1.0D0)/((l+m)*(2.0D0*l-1.0D0)))

           X1(:,m) = Alpha1 * Factor1 * Pl1(:,m) + Beta1 * Factor2 * P1Base1(:,m)
           X2(:,m) = Alpha2 * Factor1 * Pl1(:,m) + Beta2 * Factor2 * P1Base1(:,m)                
        END DO Bucle_m

        !m = l-1
        Pl1(:,l-1) = z(:) * SQRT(1D0*(1+2*l)) * P1Base1(:,l-1)

        !Calula X1 y X2
        m=l-1
        Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
        Beta1  = 2.0*((l+m)*z/(1.0-z**2))
        Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z
        Beta2  = 2.0*m/(1.0-z**2)*(l+m)

        Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

        !Introduce Beam, Pixel Window y calculo en Dl
        Factor1 = BeamPWDl(l) * Factor1

        Factor2 = Factor1 * SQRT(1.0*(l-m)*(2.0D0*l+1.0D0)/((l+m)*(2.0D0*l-1.0D0)))

        X1(:,m) = Alpha1 * Factor1 * Pl1(:,m) + Beta1 * Factor2 * P1Base1(:,m)
        X2(:,m) = Alpha2 * Factor1 * Pl1(:,m) + Beta2 * Factor2 * P1Base1(:,m)

        !m = l
        Pl1(:,l) = - SQRT(((1 + 2*l)*(1 - z**2))/(2D0*l)) * P1Base1(:,l-1)

        !Calula X1 y X2
        m = l
        Alpha1 = -2.0*((l-m*m)/(1.0-z**2)+l*(l-1)/2.0)
        Alpha2 = -2.0*m/(1.0-z**2)*(l-1)*z

        Factor1 = SQRT(1.0/((l+2.0D0)*(l+1.0D0)*l*(l-1.0D0)))

        !Introduce Beam, Pixel Window y calculo en Dl
        Factor1 = BeamPWDl(l) * Factor1

        X1(:,m) = Alpha1 * Factor1 * Pl1(:,m)
        X2(:,m) = Alpha2 * Factor1 * Pl1(:,m)
        !Ya tiene el siguiente polinomio

        base = l+1
        indice = base

        Bucle_m2: DO m=0,l

           Cosenos = COS(1.0D0 * m * Phi)
           Senos   = SIN(1.0D0 * m * Phi)

           XQE(:,indice) = dcmplx(-X1(:,m) * Cosenos, -X1(:,m) * Senos)
           XQE(:,base-m) = (-1)**m * dconjg(XQE(:,indice))

           XUE(:,indice) = -dcmplx(X2(:,m) * Senos, -X2(:,m) * Cosenos)
           XUE(:,base-m) = (-1)**m * dconjg(XUE(:,indice))

           XQB(:,indice) = -dcmplx(-X2(:,m) * Senos, X2(:,m) * Cosenos)
           XQB(:,base-m) =  (-1)**m * dconjg(XQB(:,indice))

           XUB(:,indice) = dcmplx(-X1(:,m) * Cosenos, -X1(:,m) * Senos)
           XUB(:,base-m) =  (-1)**m * dconjg(XUB(:,indice))

           indice = indice + 1

        END DO Bucle_m2

     END IF

     CALL blacs_barrier( ICTXT, "A" )

     Columna = -3 + l*l

     CALL pzgemr2d(NPix, 2*l+1, XQE, 1, 1, DescYl, MXC, 1, Columna , DescMX, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XUE, 1, 1, DescYl, MXC, 1+Npix, Columna, DescMX, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XQB, 1, 1, DescYl, MXC, 1, TotalArmonicos + Columna , DescMX, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XUB, 1, 1, DescYl, MXC, 1+NPix, TotalArmonicos + Columna, DescMX, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     IF(MYCOL==0) THEN
        !Adelanta un paso los polinomios base de la iteracion
        P1Base2(:,0:l-1) = P1Base1(:,0:l-1)
        P1Base1(:,0:l) = Pl1(:,0:l)
     END IF

  END DO

END SUBROUTINE CalculaMatrizArmonicosCompleja


SUBROUTINE CalculaMatrizFisher(ICTXT, lugar, lmax, NB, NFBloque, NCBloque, BloqueEE, BloqueBB, BloqueEB)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE CalculaBloqueMatrizFisher(ICTXT, Bloque, NFilas, NColumnas, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
       INTEGER, INTENT(in) :: ICTXT, NFilas, NColumnas, Lmax
       REAL(Kind=8), DIMENSION(NFilas, NColumnas), INTENT(in) :: Bloque
       INTEGER, DIMENSION(Nfilas), INTENT(in) :: IndicesFila
       INTEGER, DIMENSION(NColumnas), INTENT(in) :: IndicesColumna
       REAL(Kind=8), DIMENSION(2:lmax, 2:lmax), INTENT(out) :: BloqueFisher
     END SUBROUTINE CalculaBloqueMatrizFisher

     SUBROUTINE MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
       INTEGER, INTENT(in) :: Lmax, Caso1, Caso2
       REAL(Kind=8), INTENT(in) :: Factor
       REAL(Kind=8), DIMENSION(2:Lmax, 2:Lmax), INTENT(in) :: BloqueFisher
       REAL(Kind=8), DIMENSION(3, 2:Lmax, 3, 2:Lmax), INTENT(inout) :: MatrizFisher
     END SUBROUTINE MatrizFisherGlobal
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, lmax, NB, NFBloque, NCBloque
  CHARACTER(LEN=100), INTENT(in) :: Lugar
  COMPLEX(Kind=8), DIMENSION(NFBloque, NCBloque), INTENT(in) :: BloqueEE, BloqueBB, BloqueEB

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  !Tiempos
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  !Bloques
  INTEGER :: TotalArmonicos

  !Contadores
  INTEGER :: l1, l2, i,j

  !Operaciones
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Producto

  !Bloque traspuesto
  COMPLEX(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: BloqueBE
  INTEGER, DIMENSION(9) :: DescBloque
  COMPLEX(Kind=8) :: Alpha, Beta

  !Para calculo de la matriz de Fisher
  INTEGER :: Caso1, Caso2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesFila, IndicesColumna
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: BloqueFisher
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: MatrizFisher
  REAL(Kind=8) :: Factor

  INTEGER, EXTERNAL :: INDXL2G
  EXTERNAL :: BLACS_GRIDINFO, dgsum2d, descinit, pztranc, blacs_barrier

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  ALLOCATE(BloqueFisher(2:lmax, 2:lmax))
  IF(IAM==0) ALLOCATE(MatrizFisher(3, 2:lmax, 3, 2:lmax))

  ALLOCATE(IndicesFila(NFBloque))
  ALLOCATE(IndicesColumna(NCBloque))

  DO i=1,NFBloque
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     IndicesFila(i) = FLOOR(SQRT(3.0 + j))
  END DO
  DO i=1,NCBloque
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     IndicesColumna(i) = FLOOR(SQRT(3.0 + j))
  END DO

  !Calcula productos y los bloques de la matriz de Fisher
  ALLOCATE(Producto(NFBloque, NCBloque))

  !EEEE
  CALL blacs_barrier( ICTXT, "A" )
  Caso1 = 1
  Caso2 = 1
  Factor = 1.0D0
  Producto = DREAL(BloqueEE)**2 + DIMAG(BloqueEE)**2
  CALL blacs_barrier( ICTXT, "A" )
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     WRITE (*,*) "     EEEE"
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  !BBBB
  Caso1 = 2
  Caso2 = 2
  Factor = 1.0D0
  Producto = DREAL(BloqueBB)**2 + DIMAG(BloqueBB)**2
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     WRITE (*,*) "     BBBB"
  END IF

  !EEBB
  Caso1 = 1
  Caso2 = 2
  Factor = 1.0D0
  Producto =  Dreal(BloqueEB)**2 + DIMAG(BloqueEB)**2
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "     EEBB"
  END IF

  !EEEB
  Caso1 = 1
  Caso2 = 3
  Factor = 2.0D0
  Producto =  Dreal(BloqueEB) * Dreal(BloqueEE) + DIMAG(BloqueEB) * DIMAG(BloqueEE)
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "     EEEB"
  END IF

  !EBBB
  !Lo que sigue calcula el bloque 3,2 de la matriz de Fisher!!
  Caso1 = 3
  Caso2 = 2
  Factor = 2.0D0
  Producto =  Dreal(BloqueEB) * Dreal(BloqueBB) + DIMAG(BloqueEB) * DIMAG(BloqueBB)
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "     BBEB"
  END IF

  !EBEB
  Caso1 = 3
  Caso2 = 3
  Factor = 2.0D0
  !Primera parte
  Producto =  Dreal(BloqueBB) * Dreal(BloqueEE) + Dimag(BloqueBB) * Dimag(BloqueEE) 

  !Segunda parte
  !Se encesita el traspuesto de EB
  ALLOCATE(BloqueBE(NFBloque, NCBloque))
  BloqueBE = DCMPLX(0D0,0D0)
  CALL blacs_barrier( ICTXT, "A" )

  IF(iam == 0) THEN
     WRITE(*,*) "     Building the block BE"
  END IF

  Alpha = dcmplx(1.0, 0.0)
  Beta = dcmplx(0.0, 0.0)
  CALL descinit(DescBloque, TotalArmonicos, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBloque, INFO)

  !Crea el bloque conjugado traspuesto!!
  CALL pztranc(TotalArmonicos, TotalArmonicos, Alpha, BloqueEB, 1, 1, DescBloque, Beta, BloqueBE, 1, 1, DescBloque)
  Producto = Producto +  Dreal(BloqueEB) * Dreal(BloqueBE) + Dimag(BloqueEB) * Dimag(BloqueBE)


  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBloque, NCBloque, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "     EBEB"
  END IF

  DEALLOCATE(BloqueBE)

  IF(IAM==0) THEN
     OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrix.dat", Action="write")

     DO Caso1 = 1, 3
        DO l1 = 2, lmax

           DO Caso2 = 1, 3
              DO l2 = 2, lmax

                 !  MatrizFisher(Caso1,l1,Caso2,l2) = 2*PI/(l1*(l1+1))*BeamPWDl(l1)*2*PI/(l2*(l2+1))*&
                 !  &BeamPWDl(l2)*MatrizFisher(Caso1,l1,Caso2,l2)

                 WRITE (25, *)  MatrizFisher(Caso1,l1,Caso2,l2)

              END DO
           END DO

        END DO
     END DO

     CLOSE(25)
  END IF

  CALL blacs_barrier( ICTXT, "A" )

END SUBROUTINE CalculaMatrizFisher

SUBROUTINE CalculaBloqueMatrizFisher(ICTXT, Bloque, NFilas, NColumnas, IndicesFila, IndicesColumna, BloqueFisher, Lmax)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFilas, NColumnas, Lmax
  REAL(Kind=8), DIMENSION(NFilas, NColumnas), INTENT(in) :: Bloque
  INTEGER, DIMENSION(Nfilas), INTENT(in) :: IndicesFila
  INTEGER, DIMENSION(NColumnas), INTENT(in) :: IndicesColumna
  REAL(Kind=8), DIMENSION(2:lmax, 2:lmax), INTENT(out) :: BloqueFisher

  !***********

  INTEGER :: i,j,l1,l2

  EXTERNAL :: blacs_barrier,  dgsum2d

  !Suma los elementos repartidos por el bloque
  BloqueFisher = 0.0D0
  DO i=1,NFilas
     l1 = IndicesFila(i)
     IF(l1<=lmax) THEN
        DO j=1,NColumnas
           l2 = IndicesColumna(j)
           IF(l2<=lmax) BloqueFisher(l1,l2) = BloqueFisher(l1,l2) + Bloque(i,j)
        END DO
     END IF
  END DO

  BloqueFisher = BloqueFisher/2.0

  CALL blacs_barrier(ICTXT, "ALL")

  !Suma los datos repartidos por los procesadores
  CALL dgsum2d( ictxt, "ALL", "I", Lmax-1, Lmax-1, BloqueFisher, Lmax-1, 0, 0)

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE CalculaBloqueMatrizFisher



SUBROUTINE MatrizFisherGlobal(MatrizFisher, BloqueFisher, Lmax, Caso1, Caso2, Factor)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: Lmax, Caso1, Caso2
  REAL(Kind=8), INTENT(in) :: Factor
  REAL(Kind=8), DIMENSION(2:Lmax, 2:Lmax), INTENT(in) :: BloqueFisher
  REAL(Kind=8), DIMENSION(3, 2:Lmax, 3, 2:Lmax), INTENT(inout) :: MatrizFisher

  INTEGER :: l1, l2

  DO l1=2,Lmax
     DO l2=2,Lmax
        MatrizFisher(Caso1, l1, Caso2, l2) = Factor * BloqueFisher(l1,l2)
        MatrizFisher(Caso2, l2, Caso1, l1) = MatrizFisher(Caso1, l1, Caso2, l2)
     END DO
  END DO

END SUBROUTINE MatrizFisherGlobal


SUBROUTINE CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)

  IMPLICIT NONE

  INTEGER(Kind=4), INTENT(in) :: IAM, lmax, NMapas, QuitarSesgo
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  !*****************************

  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MF, IMF, Yl, Dl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: Bl, DlMedio, DlSigma, SumaDl, SumaCuardadoDl, ErrorTeorico
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: DlEscalar, DlTensorial, P1, r
  REAL(Kind=8) :: divisor, RMedio, SigmaR, ErrorR
  INTEGER :: info, dim, i, AllocateStatus

  EXTERNAL :: dpotrf, dpotri, dsymm

  dim = 3*(lmax-1)

  IF(iam==0) THEN

     !Carga la matriz de Fisher

     ALLOCATE(MF(dim, dim), STAT = AllocateStatus)

     WRITE(*,*)
     WRITE(*,*) " Loading Fisher matrix"

     OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrix.dat", Action="read")
     READ (25, *)  MF
     CLOSE(25)

     ALLOCATE(IMF(dim, dim), STAT = AllocateStatus)
     IMF = MF

     !Invierte
     WRITE(*,*) " Inverting Fisher matrix"

     info = 0
     CALL dpotrf("U", dim, IMF, dim, info)
     IF(IAM==0) THEN
        WRITE(*,*) " Cholesky factorization: ", info
     END IF

     IF(info/=0) THEN
        WRITE(*,*)
        WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        WRITE(*,*) "Fisher matrix is singular"
        !            return
     END IF

     CALL dpotri("U", dim, IMF, dim, info)
     IF(IAM==0) THEN
        WRITE(*,*) " Inversion result:       ", info
     END IF

     IF(info/=0) THEN
        WRITE(*,*)
        WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        WRITE(*,*) "Fisher matrix is singular"
        !           return
     END IF

     !Carga Bl
     WRITE(*,*) " Loading Noise Bias"
     ALLOCATE(Bl(dim), STAT = AllocateStatus)
     OPEN(Unit=25, File=TRIM(Lugar)//"/NoiseBias.dat", Action="read")
     READ (25, *)  Bl
     CLOSE(25)

     !Carga Yl y resta Bl
     WRITE(*,*) " Loading Coupled Power"
     ALLOCATE(Yl(dim,NMapas), STAT = AllocateStatus)
     OPEN(Unit=25, File=TRIM(Lugar)//"/CoupledPower.dat", Action="read")
     IF(QuitarSesgo==1) THEN
        DO i=1, NMapas
           READ (25, *)  Yl(:,i)
           Yl(:,i) = Yl(:,i) - Bl
        END DO
     ELSE
        DO i=1, NMapas
           READ (25, *)  Yl(:,i)
        END DO
     END IF
     CLOSE(25)

     !Multiplica por la inversa de MF
     WRITE(*,*) " Computing Dl"
     ALLOCATE(Dl(dim, NMapas), STAT = AllocateStatus)
     CALL dsymm("L", "U", dim, NMapas, 1.0D0, IMF, dim, Yl, dim, 0.0D0, Dl, dim)

     IF(NMapas>1) THEN

        !Medias y desviacion estandar Dl
        ALLOCATE(DlMedio(dim))
        ALLOCATE(DlSigma(dim))

        ALLOCATE(SumaDl(dim))
        ALLOCATE(SumaCuardadoDl(dim))

        DlMedio = 0.0D0
        DlSigma = 0.0D0
        SumaDl = 0.0D0
        SumaCuardadoDl = 0.0D0


        DO i=1,NMapas
           SumaDl = SumaDl + Dl(:,i)
           SumaCuardadoDl = SumaCuardadoDl + Dl(:,i)**2
        END DO

        DlMedio = SumaDl/NMapas
        DlSigma = SQRT(SumaCuardadoDl/NMapas-DlMedio**2)
        DlSigma = SQRT((1.0D0 * NMapas)/(1.0D0 * (NMapas-1))) * DlSigma

        !Guarda
        OPEN(Unit=25, File=TRIM(Lugar)//"/MeanDl.dat", Action="write")
        WRITE (25, *)  DlMedio
        CLOSE(25)

        OPEN(Unit=25, File=TRIM(Lugar)//"/SigmaDl.dat", Action="write")
        WRITE (25, *)  DlSigma
        CLOSE(25)

        DEALLOCATE(DlMedio)
        DEALLOCATE(DlSigma)
        DEALLOCATE(SumaDl)
        DEALLOCATE(SumaCuardadoDl)

     END IF

     !Guarda Dl
     WRITE(*,*) " Saving Dl"
     OPEN(Unit=25, File=TRIM(Lugar)//"/Dl.dat", Action="write")
     WRITE (25, *)  Dl
     CLOSE(25)

     !Guarda los errores teóricos
     ALLOCATE(ErrorTeorico(dim))
     DO i=1,dim
        ErrorTeorico(i) = SQRT(IMF(i,i))
     END DO
     OPEN(Unit=25, File=TRIM(Lugar)//"/FisherErrorDl.dat", Action="write")
     WRITE(25,*) ErrorTeorico
     CLOSE(25)
     DEALLOCATE(ErrorTeorico)

  END IF

END SUBROUTINE CalculaEspectroPotencia


SUBROUTINE CalculaIMCMapasMapasCross(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas, MapasCross)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
  REAL(Kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
  REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas, MapasCross
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  INTEGER, DIMENSION(9) :: DESCIMC, DESCMapas, DescTodos, DescDistribuidos
  INTEGER :: info
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Producto
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: SumaParte, RecibeSumas, ListaChi2
  REAL(Kind=8) :: SumaChi2
  INTEGER :: NPROW, NPCOL, MYROW, MYCOL

  EXTERNAL :: BLACS_GRIDINFO, descinit, pdsymm, DGSUM2D, dgerv2d, dgesd2d, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)

  CALL descinit(DESCIMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  ALLOCATE(Producto(NFMC, NCMapas))
  Producto = 0

  CALL blacs_barrier(ICTXT, 'ALL')

  !Calcula IMC.Mapas
  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, Mapas, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)
  !Devuelve IMC.Mapas
  Mapas = Producto

  !Calcula Chi2
  Producto = Mapas * MapasCross
  ALLOCATE(SumaParte(1,NCMapas))
  ALLOCATE(RecibeSumas(1,NCMapas))
  RecibeSumas = 0
  SumaParte   = 0

  SumaParte(1,:) = SUM(Producto, DIM=1)

  !Suma los datos en cada columna
  CALL blacs_barrier(ICTXT, 'ALL')
  CALL DGSUM2D(ICTXT, 'C', '1-tree', 1, NCMapas, SumaParte, 1, -1, MYCOL)

  ALLOCATE(ListaChi2(1,NMapas))
  ListaChi2 = 0 

  CALL descinit(DescTodos, 1, NMapas, 1, NMapas, 0, 0, ICTXT, 1, INFO)
  CALL descinit(DescDistribuidos, 1, NMAPAS, 1, NB, 0, 0, ICTXT, 1, INFO)

  !Pasa todos los datos a la ListaChi del primer procesador
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(1, NMapas, SumaParte, 1, 1, DescDistribuidos, ListaChi2, 1, 1, DescTodos, ictxt)
  CALL blacs_barrier(ICTXT, "ALL")

  IF((MYCOl==0).AND.(MYROW==0)) THEN
     OPEN(unit=20, file=TRIM(Lugar)//"/ChiSquare.dat", action="write")
     WRITE (20,*) ListaChi2
     CLOSE(20)
     SumaChi2 = SUM(ListaChi2)/NMapas
     WRITE(*,*) "    <x^t C^-1 z>: ", SumaChi2
     WRITE(*,*) "    Matrix size:  ", DimMC
  END IF

  !Calcula IMC.MapasCross
  Producto = 0
  CALL blacs_barrier(ICTXT, 'ALL')
  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, MapasCross, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)
  !Devuelve IMC.Mapas
  MapasCross = Producto

  DEALLOCATE(Producto)

END SUBROUTINE CalculaIMCMapasMapasCross


SUBROUTINE CalculaYlParesMapas(ICTXT, Lugar, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas, MXR, MXI, Mapas, MapasCross)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, lmax, NMapas, NB, NFMX, NCMX, NCMapas
  CHARACTER(LEN=100), INTENT(in) :: Lugar
  REAL(Kind=8), DIMENSION(NFMX, NCMX), INTENT(in) :: MXR, MXI
  REAL(Kind=8), DIMENSION(NFMX, NCMapas), INTENT(in) :: Mapas, MapasCross

  !************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  !Descriptores
  INTEGER, DIMENSION(9) :: DescMX, DescMapas, DescPMXTMapas, DescTemp, DescTodosMapas, DescUnMapa

  !Tamaños
  INTEGER :: TotalArmonicos, NPixelesProcesador, NFPMXTMapas, NFTemp

  !Contadores
  INTEGER :: l, i,j, salto

  REAL(Kind=8), PARAMETER :: Pi = 3.141592653589793D0

  !Para calcular bl de ruido diagonal en pixeles
  INTEGER :: Caso
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2, Temp3

  !Para calcular Yl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MPYMR, MPYMI, MPYMRCross, MPYMICross
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: ListaYl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:)   :: YlReestructurado, RecibeMapa

  !CHARACTER(len=24) :: systime

  INTEGER, EXTERNAL :: INDXL2G, NUMROC
  EXTERNAL :: BLACS_GRIDINFO, descinit, pdgemm, pdgemr2d, dgesd2d, dgerv2d, DGSUM2D, blacs_barrier


  !Codigo
  NPixelesProcesador = DimMC/2

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  NFPMXTMapas  = MAX(1,NUMROC(2*TotalArmonicos, NB, MYROW, 0, NPROW))

  !Multiplica Y^t por IMC.Mapas
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMX, INFO)
  CALL descinit(DescPMXTMapas, 2*TotalArmonicos, NMapas, NB, NB, 0, 0, ICTXT, NFPMXTMapas, INFO)
  CALL DESCINIT(DescMX, DimMC, 2*TotalArmonicos, NB, NB, 0, 0, ICTXT, NFMX, INFO)

  ALLOCATE(MPYMR(NFPMXTMapas, NCMapas))
  ALLOCATE(MPYMI(NFPMXTMapas, NCMapas))
  ALLOCATE(MPYMRCross(NFPMXTMapas, NCMapas))
  ALLOCATE(MPYMICross(NFPMXTMapas, NCMapas))
  MPYMR = 0
  MPYMI = 0
  MPYMRCross = 0
  MPYMICross = 0


  !Temp1 tiene filas como MX traspuesta y columnas como NMapas
  ALLOCATE(Temp1(NFPMXTMapas, NCMapas))
  Temp1 = 0


  IF(IAM==0) THEN
     WRITE (*,*) " Multiplicacion Parte real Y por mapas " !, CTIME(TIME())
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXR, 1, 1, DescMX, Mapas, 1, 1, DescMapas,&
       & 0.0D0, MPYMR, 1, 1, DescPMXTMapas)

  IF(iam==0) WRITE(*,*) " ", ctime(time())

  IF(IAM==0) THEN
     WRITE (*,*) " Multiplicacion Parte real Y por mapas cross " !, CTIME(TIME())
  END IF


  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXR, 1, 1, DescMX, MapasCross, 1, 1, DescMapas,&
       & 0.0D0, MPYMRCross, 1, 1, DescPMXTMapas)

  IF(iam==0) WRITE(*,*) " ", ctime(time())

  !*************

  IF(IAM==0) THEN
     WRITE (*,*) " Multiplicacion Parte imaginaria Y por mapas " !, CTIME(TIME())
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXI, 1, 1, DescMX, Mapas, 1, 1, DescMapas,&
       & 0.0D0, MPYMI, 1, 1, DescPMXTMapas)

  IF(iam==0) WRITE(*,*) " ", ctime(time())

  IF(IAM==0) THEN
     WRITE (*,*) " Multiplicacion Parte imaginaria Y por mapas cross " !, CTIME(TIME())
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("T", "N", 2*TotalArmonicos, NMapas, DimMC, 1.0D0, MXI, 1, 1, DescMX, MapasCross, 1, 1, DescMapas,&
       & 0.0D0, MPYMICross, 1, 1, DescPMXTMapas)

  IF(iam==0) WRITE(*,*) " ", ctime(time())

  IF(IAM==0) THEN
     WRITE (*,*) " Calcula Yl cross sumando partes de bloques"
  END IF

  CALL blacs_barrier( ICTXT, "A" )

  !Yl EE y BB
  Temp1 = MPYMR*MPYMRCross+MPYMI*MPYMICross

  ALLOCATE(ListaYl(2:lMax,3,NCMapas))
  ListaYl = 0.0D0

  DO i=1,NFPMXTMapas
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     IF(j<= TotalArmonicos) THEN
        Caso = 1
     ELSE
        Caso = 2
        j = j-TotalArmonicos
     ENDIF
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) ListaYl(l,Caso,:) = ListaYl(l,Caso,:) + Temp1(i,:)
  END DO


  !Yl EB

  DEALLOCATE(Temp1)

  !Ahora los bloques tienen numero de filas: TotalArmonicos
  NFTemp =MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))

  ALLOCATE(Temp1(NFTemp, NCMapas))
  ALLOCATE(Temp2(NFTemp, NCMapas))
  ALLOCATE(Temp3(NFTemp, NCMapas))

  Temp1 = 0
  Temp2 = 0
  Temp3 = 0

  CALL DESCINIT(DescTemp, TotalArmonicos, NMapas, NB, NB, 0, 0, ICTXT, NFTemp, INFO)

  !Pasa la parte real de Q y de U
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMR, 1, 1, DescPMXTMapas, Temp1, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMRCross, 1+TotalArmonicos, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')

  Temp1 = Temp1*Temp2

  !Pasa la parte imaginaria de Q de U
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMI, 1, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMICross, 1+TotalArmonicos, 1, DescPMXTMapas, Temp3, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')

  Temp1 = Temp1 + Temp2*Temp3

  DO i=1,NFTemp
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) ListaYl(l,3,:) = ListaYl(l,3,:) + Temp1(i,:)
  END DO

  !Pasa la parte real de Q y de U
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMRCross, 1, 1, DescPMXTMapas, Temp1, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMR, 1+TotalArmonicos, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')

  Temp1 = Temp1*Temp2

  !Pasa la parte imaginaria de Q de U
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMICross, 1, 1, DescPMXTMapas, Temp2, 1, 1, DescTemp, ictxt)
  CALL blacs_barrier(ICTXT,'All')
  CALL pdgemr2d(TotalArmonicos, NMapas, MPYMI, 1+TotalArmonicos, 1, DescPMXTMapas, Temp3, 1, 1, DescTemp, ictxt)

  Temp1 = Temp1 + Temp2*Temp3

  DO i=1,NFTemp
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) ListaYl(l,3,:) = ListaYl(l,3,:) + Temp1(i,:)
  END DO


  ! !Aplica fator m_u_l
  ! do l=2,lMax
  !     ListaYl(l,:,:) = ListaYl(l,:,:) * BeamPWDl(l)*2*Pi/(l*(l+1))/2.0D0
  ! end do

  ListaYl = ListaYl/2.0D0

  DEALLOCATE(Temp1)
  DEALLOCATE(Temp2)
  DEALLOCATE(Temp3)

  !Cada procesador tiene una parte de Yl

  !Suma Yl en las columnas
  !La suma esta lanzada en una forma que no tiene en cuenta la forma de la variable ListaYl
  CALL blacs_barrier( ICTXT, "A" )
  CALL DGSUM2D( ictxt, 'C', '1-tree', 3*(lMax-1)*NCMapas, 1, ListaYl, 3*(lMax-1)*NCMapas, -1, MYCOL)

  !Todos los procesadores de la misma columna ya tienen completos los Yl de su parte de los mapas
  !Los procesadores de la primera fila le envían los datos al procesador 0, y los guarda

  salto = lmax-1
  ALLOCATE(YlReestructurado(3*salto, NCMapas))
  DO i=1,NCMapas
     YlReestructurado(1:salto,i) = ListaYl(:,1,i)
     YlReestructurado(salto+1:2*salto,i) = ListaYl(:,2,i)
     YlReestructurado(2*salto+1:3*salto,i) = ListaYl(:,3,i)
  END DO

  ALLOCATE(RecibeMapa(3*salto,1))

  CALL descinit(DescTodosMapas, 3*salto, NMapas, 3*salto, NB, 0, 0, ICTXT, 3*salto, INFO)
  CALL descinit(DescUnMapa, 3*salto, 1, 3*salto, NB, 0, 0, ICTXT, 3*salto, INFO)


  IF (IAM==0) THEN
     OPEN(Unit=20, File=TRIM(Lugar)//"/CoupledPower.dat", Action="write")
  END IF

  DO i=1,NMapas

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(3*salto, 1, YlReestructurado, 1, i, DescTodosMapas, RecibeMapa, 1, 1, DescUnMapa, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")
     IF(iam==0) WRITE(20,*) RecibeMapa(:,1)

  END DO

  IF (IAM==0) THEN
     CLOSE(20)
  END IF

  DEALLOCATE(ListaYl)
  DEALLOCATE(YlReestructurado)
  DEALLOCATE(RecibeMapa)   

  !Fin calcula Yl
  CALL blacs_barrier( ICTXT, "A" )

END SUBROUTINE CalculaYlParesMapas


!*************************************************************************


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

     !write(*,*) "NAXIS: ", naxis

     !     determines the presence of an extension
     CALL ftgkyl(unit,'EXTEND', extend, comment, status)
     IF (status > 0) THEN
        ftype_in = 0
        status = 0 ! no extension : 
        !     to be compatible with first version of the code
     ENDIF

     !write(*,*) "EXTEND: ", extend

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

     !write(*,*) "NSIZE: ", NSIZE
     !write(*,*) "NMAPS: ", nmaps_in

     CALL ftgkys(unit,'ORDERING',order_val,comment,status)
     IF (status == 202) THEN ! Ordering not found
        ordering_in = 0
        order_val = ''
        status = 0
     ENDIF

     !write(*,*) "ORDERING: ", order_val

     CALL ftgkyj(unit,'NSIDE',nside_in,comment,status)
     IF (status == 202) THEN ! Nside not found
        nside_in = -1
        status = 0
     ENDIF

     !write(*,*) "NSIDE: ", nside_in

     CALL ftgkyj(unit,'OBS_NPIX',obs_npix_in,comment,status)
     IF (status == 202) THEN ! obs_npix not found
        obs_npix_in = -1
        status = 0
     ENDIF

     !write(*,*) "OBS_NPIX: ", obs_npix_in

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


SUBROUTINE CargaMascara(ictxt, nside, Lugar, File, NumMapaCargar, Mascara)

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: imap
       REAL(Kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
       !real(Kind=8),                                intent(OUT) :: nullval
       INTEGER                   , INTENT(IN) :: extno
     END SUBROUTINE read_bintab_UnMapa
  END INTERFACE


  INTEGER, INTENT(in) :: ictxt, Nside, NumMapaCargar
  CHARACTER(len=100), INTENT(in) :: Lugar, File
  INTEGER, DIMENSION(12*nside**2), INTENT(out) :: Mascara

  CHARACTER(len=200) :: ArchivoMascara
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: MascaraReal
  INTEGER :: nmaps, extno, npix_full

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, igebs2d, igebr2d

  !*********************************************************************

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  npix_full = 12*Nside**2

  IF(iam==0) THEN

     ArchivoMascara = TRIM(Lugar)//"/"//TRIM(File)

     ALLOCATE(MascaraReal(0:npix_full-1))

     nmaps = NumMapaCargar
     extno = 0
     CALL read_bintab_UnMapa(ArchivoMascara, MascaraReal, npix_full, NumMapaCargar, extno)

     WHERE(MascaraReal < -1.637E+030) MascaraReal = 0 

     Mascara(1:npix_full) = INT(MascaraReal(0:npix_full-1))

     DEALLOCATE(MascaraReal)

  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  !Pasa la máscara

  IF(iam==0) THEN
     CALL igebs2d(ICTXT, "All","I" , npix_full, 1, Mascara, npix_full)
  ELSE
     CALL igebr2d(ICTXT, "All","I" , npix_full, 1, Mascara, npix_full, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")


END SUBROUTINE CargaMascara


SUBROUTINE read_bintab(filename, map, npixtot, nmaps, extno)

  IMPLICIT NONE

  !=======================================================================
  CHARACTER(len=*),                          INTENT(IN)  :: filename
  INTEGER,                              INTENT(IN)  :: npixtot
  INTEGER,                              INTENT(IN)  :: nmaps
  REAL(Kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
  INTEGER, INTENT(IN) :: extno

  INTEGER :: status,unit,readwrite,naxes(2),nfound, naxis
  INTEGER :: group, firstpix, i, npix32
  REAL(Kind=8)   :: blank, testval
  REAL(Kind=8)     :: bscale,bzero
  CHARACTER(len=80) :: comment
  LOGICAL(Kind=8) :: extend
  INTEGER :: nmove, hdutype, hdunum
  INTEGER :: frow, imap
  INTEGER :: datacode, width
  LOGICAL(Kind=8) ::  anynull_i

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

  REAL(Kind=8) :: ABS

  LOGICAL :: anynull
  REAL(Kind=8) ::  nullval

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


     !  if (nl_header > 0) then
     !     do i=1,nl_header
     !        header(i)(1:len_header) = ""
     !     enddo
     !     call get_clean_header(unit, header, filename, status)
     !     status = 0
     !  endif

     !  if (nl_units > 0) then
     !     do i=1,nl_units
     !        units(i)(1:len_units) = 'unknown' ! default
     !     enddo
     !     do imap = 1, min(nmaps, nl_units)
     !        units(imap)(1:len_units) = adjustl(tunit(imap))
     !     enddo
     !  endif


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
  REAL(Kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
  INTEGER                   , INTENT(IN) :: extno

  INTEGER :: status,unit,readwrite,naxes(2),nfound, naxis
  INTEGER ::  firstpix
  REAL(Kind=8)   :: blank
  REAL(Kind=8)     :: bscale,bzero
  CHARACTER(len=80) :: comment
  LOGICAL(Kind=8) :: extend
  INTEGER :: nmove, hdutype, hdunum
  INTEGER :: frow
  INTEGER :: datacode, width
  LOGICAL(Kind=8) ::  anynull_i

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
  EXTERNAL :: ftgrsz, ftbnfm, ftgcvd, ftclos,  printerror

  LOGICAL :: anynull
  REAL(Kind=8) ::  nullval
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

     ! sanity check
     !do imap = 1, nmaps
     IF (i0(imap) /= npix(imap)) THEN
        !call fatal_error('something wrong during piece wise reading')
     ENDIF
     !enddo

  ELSE ! no image no extension, you are dead, man
     !call fatal_error(' No image, no extension in '//trim(filename))
  ENDIF
  !     close the file
  CALL ftclos(unit, status)

  !     check for any error, and if so print out error messages
  IF (status > 0) CALL printerror(status)

  RETURN
END SUBROUTINE read_bintab_UnMapa


SUBROUTINE pix2ang_ring(nside, ipix, z, phi)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nside
  INTEGER, INTENT(IN)  :: ipix
  REAL(Kind=8),     INTENT(OUT) :: z, phi

  INTEGER ::  nl2, nl4, iring, iphi
  INTEGER ::  npix, ncap, ip
  REAL(Kind=8) ::  fodd, dnside
  REAL(Kind=8), PARAMETER :: half = 0.500000000000000D0
  !REAL(Kind=8), parameter :: one  = 1.000000000000000D0
  !REAL(Kind=8), parameter :: three = 3.00000000000000D0
  REAL(Kind=8), PARAMETER :: threehalf = 1.50000000000000D0

  REAL(Kind=8), PARAMETER :: PI = 3.1415926535897932384626433832795D0
  REAL(Kind=8), PARAMETER :: HALFPI = 1.5707963267948966192313216916398D0

  CHARACTER(len=*), PARAMETER :: code = "pix2ang_ring"
  !-----------------------------------------------------------------------

  npix = 12*nside**2       ! total number of points
  !if (ipix <0 .or. ipix>npix-1) call fatal_error (code//"> ipix out of range")

  nl2  = 2*nside
  ncap = nl2*(nside-1 ) ! points in each polar cap, =0 for nside =1
  dnside = REAL(nside, Kind=8)

  IF (ipix < ncap) THEN ! North Polar cap -------------

     iring = (SQRT(2.0D0*ipix+2.0D0) + 1D0)/2D0
     iphi  = ipix - 2*iring*(iring - 1 )

     z = COS(2.0D0 * ASIN(iring / (SQRT(6.0D0)*dnside)))
     phi   = (REAL(iphi,Kind=8) + half) * HALFPI/iring

  ELSEIF (ipix < npix-ncap) THEN ! Equatorial region ------

     ip    = ipix - ncap
     nl4   = 4*nside
     iring = INT( ip / nl4 ) + nside ! counted from North pole
     iphi  = IAND(ip, nl4-1 )

     fodd  = half * ( IAND(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
     z =  (nl2 - iring) / (threehalf*dnside)
     phi   = (REAL(iphi,Kind=8) + fodd) * HALFPI / dnside

  ELSE ! South Polar cap -----------------------------------

     ip    = npix - ipix
     iring = (SQRT(2D0*ip) + 1) / 2
     iphi  = 2*iring*(iring + 1 ) - ip


     z = COS(PI - 2.d0 * ASIN(iring / (SQRT(6.0D0)*dnside)))
     phi   = (REAL(iphi,Kind=8) + half) * HALFPI/iring

  ENDIF

  RETURN

END SUBROUTINE pix2ang_ring



SUBROUTINE PuntosPixelesObservados(ictxt, nside, npix, Mascara, NB, NPixelesProcesador, z, a)

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE pix2ang_ring(nside, ipix, z, phi)
       INTEGER, INTENT(IN)  :: nside
       INTEGER, INTENT(IN)  :: ipix
       REAL(Kind=8),     INTENT(OUT) :: z, phi
     END SUBROUTINE pix2ang_ring

  END INTERFACE

  INTEGER, INTENT(in) :: ictxt, nside, npix, NB, NPixelesProcesador
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
  REAL(Kind=8), DIMENSION(NPixelesProcesador) :: z
  REAL(Kind=8), DIMENSION(NPixelesProcesador) :: a

  !***************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO
  INTEGER :: npix_full, ipix, j
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: TodosZ, TodosA
  REAL(Kind=8) :: TempZ, TempA

  INTEGER, DIMENSION(9) :: DescCarga, DescReparto

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, descinit, pdgems2d, pdgemr2d

  !*********************************************************************

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  npix_full = 12*nside**2

  ALLOCATE(TodosZ(npix))
  ALLOCATE(TodosA(npix))

  TodosZ = -1
  TOdosA = -2
  z = -1000
  a = -2000

  IF(iam==0) THEN

     j = 1
     DO ipix=1,npix_full
        IF(Mascara(ipix)==1) THEN
           CALL pix2ang_ring(nside, ipix-1, TempZ, TempA)
           ! write(*,*) "pix: ", nside, ipix-1, j, TempZ, TempA
           TodosZ(j) = TempZ
           TodosA(j) = TempA
           j = j + 1
        END IF
     END DO

  END IF



  CALL descinit(DescCarga, npix, 1, npix, 1, 0, 0, ICTXT, npix, INFO)
  CALL descinit(DescReparto, npix, 1, NB, 1, 0, 0, ICTXT, NPixelesProcesador, INFO)

  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(npix, 1, TodosZ, 1, 1, DescCarga, z, 1, 1, DescReparto, ictxt)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(npix, 1, TodosA, 1, 1, DescCarga, a, 1, 1, DescReparto, ictxt)

  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN

     DEALLOCATE(TodosZ)
     DEALLOCATE(TodosA)

  END IF


END SUBROUTINE PuntosPixelesObservados


SUBROUTINE BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE gaussbeam(fwhm_arcmin, lmax, gb)
       REAL(Kind=8),                   INTENT(in)  :: fwhm_arcmin
       INTEGER,               INTENT(in)  :: lmax
       REAL(Kind=8), DIMENSION(0:lmax,1:3), INTENT(out) :: gb
     END SUBROUTINE gaussbeam

     SUBROUTINE get_pixel_window(dir, nside, pixel)
       CHARACTER(len=*), INTENT(in) :: dir
       INTEGER, INTENT(in) :: nside
       REAL(Kind=8), DIMENSION(0:4*NSide, 3), INTENT(out) :: pixel
     END SUBROUTINE get_pixel_window

  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, PixelWindowEnMapas, Nside, lmax
  REAL(Kind=8) :: FWHM_Beam
  CHARACTER(len=*), INTENT(in) :: DirHealpixData                          !Lugar donde estan los datos Healpix
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(out) :: BeamPWDl

  !********

  INTEGER :: LmaxTemp, l
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Beam, Pixel

  REAL(Kind=8), PARAMETER :: DosPi = 6.2831853071795864769D0

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL

  EXTERNAL ::  BLACS_GRIDINFO, blacs_barrier, dgebs2d, dgebr2d

  !***************************

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  LmaxTemp = 4*nside

  BeamPWDl = 1

  IF(iam==0) THEN
     !Beam
     ALLOCATE(beam(0:LmaxTemp, 3))
     CALL gaussbeam(FWHM_Beam, LmaxTemp, beam)

     !Pixel_Window
     ALLOCATE(Pixel(0:LmaxTemp,3))

     IF(PixelWindowEnMapas==1) THEN
        CALL get_pixel_window(TRIM(DirHealpixData), nside, Pixel)
        !CALL get_pixel_window(TRIM(DirHealpixData), 64, Pixel)
     ELSE
        Pixel = 1D0
     ENDIF

     !Prepara la salida de los datos
     DO l=1,lmax
        BeamPWDl(l) = SQRT(DosPi/(REAL(l)*(REAL(l)+1))) * Beam(l,2) * pixel(l,2)
     END DO

  END IF

  CALL blacs_barrier(ictxt, "all")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , Lmax+1, 1, BeamPWDl, Lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , Lmax+1, 1, BeamPWDl, Lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ictxt, "all")


END SUBROUTINE BeamPixelWindow


SUBROUTINE gaussbeam(fwhm_arcmin, lmax, gb)

  IMPLICIT NONE


  !===========================================================
  ! gaussbeam(fwhm_arcmin, gb)
  !   returns in gb the window function on [0,lmax] corresponding
  !   to the gaussian beam of FWHM = fwhm_arcmin
  ! The returned beam function has up to 3 components,
  ! gb(*,1) = bt                  : temperature
  ! gb(*,2) = bt * exp(2*sigma^2) : grad
  ! gb(*,3) = bt * exp(2*sigma^2) : curl
  ! with sigma = gaussian rms in radian
  ! and bt = exp( l(l+1) * sigma^2 / 2)
  !===========================================================
  REAL(Kind=8),                   INTENT(in)  :: fwhm_arcmin
  REAL(Kind=8), DIMENSION(0:lmax,1:3), INTENT(out) :: gb
  INTEGER,               INTENT(in)  :: lmax

  INTEGER :: l, nd
  REAL(Kind=8)     :: sigma, arcmin2rad, sigma2fwhm, fact_pol
  !===========================================================

  REAL(Kind=8), PARAMETER :: PI = 3.1415926535897932384626433832795D0


  nd   = SIZE(gb,2)

  arcmin2rad = PI / (180.0D0 * 60.0D0)
  sigma2fwhm = SQRT(8.0D0 * LOG(2.0D0))

  sigma    = fwhm_arcmin * arcmin2rad / sigma2fwhm ! in radians

  fact_pol = EXP(2.0D0*sigma**2) ! correction for polarised fields

  ! temperature
  DO l=0,lmax
     gb(l,1) = EXP(-.5D0 * l*(l+1.0D0) * sigma**2)
  ENDDO
  ! electric or gradient
  IF (nd > 1) gb(0:lmax,2) = gb(0:lmax,1) * fact_pol
  ! magnetic or curl
  IF (nd > 2) gb(0:lmax,3) = gb(0:lmax,1) * fact_pol

  RETURN
END SUBROUTINE gaussbeam



SUBROUTINE get_pixel_window(dir, nside, pixel)

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE read_bintab(filename, map, npixtot, nmaps, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: nmaps
       REAL(Kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
       INTEGER,                   INTENT(IN) :: extno
     END SUBROUTINE read_bintab

  END INTERFACE

  CHARACTER(len=*), INTENT(in) :: dir
  INTEGER, INTENT(in) :: nside
  REAL(Kind=8), DIMENSION(0:4*NSide, 3), INTENT(out) :: pixel

  !******************************

  CHARACTER(100) :: file
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: PixelTemp


  SELECT CASE(NSide)
  CASE(1)
     Pixel = 1.0D0
     RETURN
  CASE(2)
     file = "/pixel_window_n0002.fits"
  CASE(4)
     file = "/pixel_window_n0004.fits"
  CASE(8)
     file = "/pixel_window_n0008.fits"
  CASE(16)
     file = "/pixel_window_n0016.fits"
  CASE(32)
     file = "/pixel_window_n0032.fits"
  CASE(64)
     file = "/pixel_window_n0064.fits"
  CASE(128)
     file = "/pixel_window_n0128.fits"
  CASE(256)
     file = "/pixel_window_n0256.fits"
  CASE(512)
     file = "/pixel_window_n0512.fits"
  END SELECT

  file = TRIM(dir)//TRIM(file)

  ALLOCATE(PixelTemp(0:4*NSide,2))

  CALL read_bintab(file, PixelTemp, 4*NSide+1, 2, 0)

  Pixel(:,1) = PixelTemp(:,1)
  Pixel(:,2) = PixelTemp(:,2)
  Pixel(:,3) = PixelTemp(:,2)

  DEALLOCATE(PixelTemp)

  RETURN

END SUBROUTINE get_pixel_window


SUBROUTINE CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileData, TipoDatos, Mapas)

  USE ifport


  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: imap
       REAL(Kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
       !real(Kind=8),                                intent(OUT) :: nullval
       INTEGER                   , INTENT(IN) :: extno
     END SUBROUTINE read_bintab_UnMapa
  END INTERFACE


  INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas, TipoDatos
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
  CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
  REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas

  !******************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i, imap, INFO, NpixFull, Salto, extno, nmaps
  INTEGER :: npix, ipixQ, ipixU
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MapasFullSkyQ, MapasFullSkyU, MapasQUMascara
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: Mapa
  CHARACTER(len=200)  :: FileName
  INTEGER, DIMENSION(9) :: DescCarga, DescMapas

  EXTERNAL :: BLACS_GRIDINFO, descinit, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  FileName = TRIM(Lugar)//"/"//TRIM(FileData)

  NpixFull = 12*nside**2

  ALLOCATE(Mapa(0:NpixFull-1))

  IF(MYROW==0) THEN
     ALLOCATE(MapasFullSkyQ(NpixFull, NCMapas))
     ALLOCATE(MapasFullSkyU(NpixFull, NCMapas))
  ELSE    
     ALLOCATE(MapasFullSkyQ(1, NCMapas))
     ALLOCATE(MapasFullSkyU(1, NCMapas))
  END IF


  IF(IAM==0) THEN
     WRITE (*,*) "Cargando Mapas: "
     WRITE(*,*) " ", ctime(time())
  END IF

  !Los datos pueden ser TQU o solo QU

  !Caso TQU
  IF(TipoDatos==0) THEN
     nmaps = 3*NMapas
     salto = 1
     !Caso QU    
  ELSE
     nmaps = 2*NMapas
     salto = 0
  END IF

  !Descriptores para trabajar solo en la primera fila
  CALL descinit(DescCarga, NpixFull, 1, NpixFull, 1, 0, 0, ICTXT, NpixFull, INFO)
  CALL descinit(DescMapas, NpixFull, NMAPAS, NpixFull, NB, 0, 0, ICTXT, NpixFull, INFO)

  extno = 0
  imap = 1 + salto
  DO i=1,NMapas

     !Carga Q
     IF(iam==0) THEN
        !write(*,*) "Carga Q: ", i, imap, extno
        CALL read_bintab_UnMapa(FileName, Mapa, NpixFull, imap, extno)
     END IF

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, 1, Mapa, 1, 1, DescCarga, MapasFullSkyQ, 1, i, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")

     !Carga U
     imap = imap + 1
     IF(iam==0) THEN
        !write(*,*) "Carga U: ", i, imap, extno
        CALL read_bintab_UnMapa(FileName, Mapa, NpixFull, imap, extno)
     END IF

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, 1, Mapa, 1, 1, DescCarga, MapasFullSkyU, 1, i, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")

     !Prepara indice para el siguiente mapa    
     imap = imap + 1 + salto

     !Controla cambio de extension
     IF(MOD(i,50)==0) THEN
        imap = 1 + salto
        extno = extno + 1
        IF(iam==0) WRITE(*,*) "Salto de extension: ", extno 
     END IF

  END DO

  DEALLOCATE(Mapa)

  !Aplica mascara
  npix = SUM(Mascara)
  IF(MYROW==0) THEN
     ALLOCATE(MapasQUMascara(2*npix, NCMapas))
  ELSE
     ALLOCATE(MapasQUMascara(1, NCMapas))
  END IF
  MapasQUMascara = 0

  ipixQ = 1
  ipixU = npix+1
  IF(MYROW==0) THEN
     DO i=1,NpixFull
        IF(Mascara(i)==1) THEN
           MapasQUMascara(ipixQ,:) = MapasFullSkyQ(i,:)
           ipixQ = ipixQ + 1
           MapasQUMascara(ipixU,:) = MapasFullSkyU(i,:)
           ipixU = ipixU + 1 
        END IF
     END DO
  END IF

  DEALLOCATE(MapasFullSkyQ)
  DEALLOCATE(MapasFullSkyU)

  CALL blacs_barrier(ICTXT, "ALL")

  !Reparte los mapas entre los procesadores de la columna
  Mapas = 0
  !Lo pasa todo a la vez
  IF(IAM==0) WRITE(*,*) "Envia"
  CALL descinit(DescCarga, 2*NPix, NMapas, 2*NPix, NB, 0, 0, ICTXT, 2*NPix, INFO)
  CALL descinit(DescMapas, 2*NPix, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL pdgemr2d(2*NPix, NMapas, MapasQUMascara, 1, 1, DescCarga, Mapas, 1, 1, DescMapas, ictxt)

  DEALLOCATE(MapasQUMascara)

END SUBROUTINE CargaMapasFits


SUBROUTINE CargaMapaRuido(ICTXT, TipoRuido, RuidoQQ, nside, NFMC, NFBloque, NB, Mascara, Lugar, FileData, TipoDatos,&
     &  MapaRuido2, MapaRuido2Bloque)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: imap
       REAL(Kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
       INTEGER                   , INTENT(IN) :: extno
     END SUBROUTINE read_bintab_UnMapa
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, nside, NFMC, NFBloque, NB, TipoRuido, TipoDatos
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
  CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
  REAL(Kind=8), INTENT(inout) :: RuidoQQ
  REAL(Kind=8), DIMENSION(NFMC), INTENT(out) :: MapaRuido2
  REAL(Kind=8), DIMENSION(NFBloque), INTENT(out) :: MapaRuido2Bloque


  !******************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i, k, imap, NpixFull, Salto, extno
  INTEGER :: npix
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: MapaDoble, MapaRuidoMascara
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Mapa
  CHARACTER(len=200)  :: FileName


  EXTERNAL :: BLACS_GRIDINFO, descinit, blacs_barrier, pdgemrsd, dgebr2d, dgebs2d
  INTEGER, EXTERNAL :: INDXL2G

  IF(TipoRuido==0) THEN
     MapaRuido2 = RuidoQQ**2
     MapaRuido2Bloque = RuidoQQ**2
     RETURN
  END IF


  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  FileName = TRIM(Lugar)//"/"//TRIM(FileData)

  NpixFull = 12*nside**2

  ALLOCATE(Mapa(NpixFull,1))
  Mapa = 0

  !Los datos pueden ser TQU o solo QU

  !Caso Temperatura y Polarizacion
  IF(TipoDatos==0) THEN
     salto = 1
     !Caso Polarizacion   
  ELSE
     salto = 0
  END IF

  extno = 0
  imap = 1 + salto

  !Carga Mapa
  IF(iam==0) THEN
     CALL read_bintab_UnMapa(FileName, Mapa, NpixFull, imap, extno)
  END IF

  CALL blacs_barrier(ICTXT,"all")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , NpixFull, 1, Mapa, NpixFull)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , NpixFull, 1, Mapa, NpixFull, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT,"all")

  npix = SUM(Mascara)
  ALLOCATE(MapaRuidoMascara(npix))

  k=1
  DO i=1,NpixFull
     IF(Mascara(i)==1) THEN
        MapaRuidoMascara(k) = Mapa(i,1)
        k = k + 1
     END IF
  END DO

  !Calcula la media armonica, para tener el ruido por pixel equivalente
  !Lo usa para binear
  RuidoQQ = REAL(npix)/SUM(1.0D0/MapaRuidoMascara)

  !Lo pasa a cuadrados
  MapaRuidoMascara = MapaRuidoMascara**2

  !Toma elementos en MapaRuido2Bloque
  DO i=1,NFBloque
     MapaRuido2Bloque(i) = MapaRuidoMascara(INDXL2G(i, NB, MYROW, 0, NPROW))
  END DO

  ALLOCATE(MapaDoble(2*npix))
  MapaDoble(1:npix)        = MapaRuidoMascara(1:npix)
  MapaDoble(npix+1:2*npix) = MapaRuidoMascara(1:npix)

  !Toma elementos en MapaRuido2
  DO i=1,NFMC
     MapaRuido2(i) = MapaDoble(INDXL2G(i, NB, MYROW, 0, NPROW))
  END DO

  DEALLOCATE(Mapa)
  DEALLOCATE(Mapadoble)
  DEALLOCATE(MapaRuidoMascara)

  CALL blacs_barrier(ICTXT,"all")

END SUBROUTINE CargaMapaRuido


SUBROUTINE CargaMapasFits2(ICTXT, nside, NMapas, NB, NFMC, NCMapas, Mascara, Lugar, FileData, TipoDatos, Mapas)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab(filename, map, npixtot, nmaps, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: nmaps
       REAL(Kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
       INTEGER, INTENT(IN) :: extno
     END SUBROUTINE read_bintab
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas, TipoDatos
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
  CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
  REAL(Kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas

  !******************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i,j,k, imap, INFO, NpixFull, extno, numext
  INTEGER :: npix, ipixQ, ipixU, NumMapasCargar, NumMapasPasar, Faltan, Pasados
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: MapasFullSkyQ, MapasFullSkyU, MapasQUMascara
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Mapa, MapasPasar
  CHARACTER(len=200)  :: FileName
  INTEGER, DIMENSION(9) :: DescCarga, DescMapas

  EXTERNAL :: BLACS_GRIDINFO, descinit, blacs_barrier, pdgemr2d

  Mapas = -1234

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  FileName = TRIM(Lugar)//"/"//TRIM(FileData)

  NpixFull = 12*nside**2

  NumExt = INT(NMapas/50)
  IF(NumExt>0) THEN
     NumMapasCargar = 50
  ELSE
     NumMapasCargar = NMapas
  END IF

  !Pasa Q y U por separado
  NumMapasPasar = NumMapasCargar

  !Pasa a dimensiones TQU o QU
  IF(TipoDatos==0) THEN
     NumMapasCargar = 3*NumMapasCargar
  ELSE
     NumMapasCargar = 2*NumMapasCargar
  END IF

  ALLOCATE(Mapa(0:NpixFull-1, NumMapasCargar))
  ALLOCATE(MapasPasar(0:NpixFull-1, NumMapasPasar))

  IF(MYROW==0) THEN
     ALLOCATE(MapasFullSkyQ(NpixFull, NCMapas))
     ALLOCATE(MapasFullSkyU(NpixFull, NCMapas))
  ELSE    
     ALLOCATE(MapasFullSkyQ(1, NCMapas))
     ALLOCATE(MapasFullSkyU(1, NCMapas))
  END IF




  !Descriptores para trabajar solo en la primera fila
  CALL descinit(DescCarga, NpixFull, NumMapasPasar, NpixFull, NumMapasPasar, 0, 0, ICTXT, NpixFull, INFO)
  CALL descinit(DescMapas, NpixFull, NMAPAS, NpixFull, NB, 0, 0, ICTXT, NpixFull, INFO)

  extno = 0
  imap = 1
  i=0
  Pasados = 0
  DO i=0,NumExt-1

     !if(iam==0) write(*,*) "Cargando mapas hasta: ", Pasados+50 

     !Carga
     IF(iam==0) THEN
        CALL read_bintab(filename, Mapa, NPixFull, NumMapasCargar, extno)

        !Prepara el paso Q
        IF(TipoDatos==0) THEN
           !Deja fuera los mapas T
           k = 1
           DO j=1,150,3
              MapasPasar(:,k) = Mapa(:,j+1)
              !write(*,*) "exto, Q: ", i, j, k  
              k = k + 1

           END DO
        ELSE
           k = 1
           DO j=1,100,2
              MapasPasar(:,k) = Mapa(:,j)
              !write(*,*) "exto, Q: ", i, j, k 
              k = k + 1

           END DO
        END IF

     END IF

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasPasar, MapasPasar, 1, 1, DescCarga, MapasFullSkyQ, 1,&
          & 1+NumMapasPasar*i, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")

     !Prepara el paso U
     IF(iam==0) THEN
        IF(TipoDatos==0) THEN
           !Deja fuera los mapas T y Q
           k = 1
           DO j=1,150,3
              MapasPasar(:,k) = Mapa(:,j+2)
              !write(*,*) "exto, U: ", i, j, k   
              k = k + 1

           END DO
        ELSE
           k = 1
           DO j=1,100,2
              MapasPasar(:,k) = Mapa(:,j+1)
              !write(*,*) "exto, U: ", i, j, k   
              k = k + 1

           END DO
        END IF
     END IF


     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasPasar, MapasPasar, 1, 1, DescCarga, MapasFullSkyU, 1,&
          & 1+NumMapasPasar*i, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")

     !Prepara indice para el siguiente mapa    
     extno = extno + 1

     Pasados = Pasados + 50

  END DO

  !Carga los mapas restantes
  IF(Pasados < NMapas) THEN

     Faltan = NMapas - (NumExt*50)

     !if(iam==0) write(*,*) "Cargando de la ultima extension ", Faltan, " mapas"

     !Pasa a dimensiones TQU o QU
     IF(TipoDatos==0) THEN
        NumMapasCargar = 3*Faltan
     ELSE
        NumMapasCargar = 2*Faltan
     END IF

     NumMapasPasar  = Faltan

     DEALLOCATE(Mapa)
     DEALLOCATE(MapasPasar)

     ALLOCATE(Mapa(0:NpixFull-1, NumMapasCargar))
     ALLOCATE(MapasPasar(0:NpixFull-1, NumMapasPasar))

     CALL descinit(DescCarga, NpixFull, NumMapasPasar, NpixFull, NumMapasPasar, 0, 0, ICTXT, NpixFull, INFO)

     !Carga
     IF(iam==0) THEN
        CALL read_bintab(filename, Mapa, NPixFull, NumMapasCargar, extno)

        !Prepara el paso Q
        IF(TipoDatos==0) THEN
           !Deja fuera los mapas T
           k = 1
           DO j=1,NumMapasCargar,3
              MapasPasar(:,k) = Mapa(:,j+1)
              !write(*,*) "exto, Q: ", i, j, k  
              k = k + 1

           END DO
        ELSE
           k = 1
           DO j=1,NumMapasCargar,2
              MapasPasar(:,k) = Mapa(:,j)
              !write(*,*) "exto, Q: ", i, j, k 
              k = k + 1

           END DO
        END IF

     END IF

     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasPasar, MapasPasar, 1, 1, DescCarga, MapasFullSkyQ, 1,&
          & 1+Pasados, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")

     !Prepara el paso U
     IF(iam==0) THEN
        IF(TipoDatos==0) THEN
           !Deja fuera los mapas T y Q
           k = 1
           DO j=1,NumMapasCargar,3
              MapasPasar(:,k) = Mapa(:,j+2)
              !write(*,*) "exto, U: ", i, j, k   
              k = k + 1

           END DO
        ELSE
           k = 1
           DO j=1,NumMapasCargar,2
              MapasPasar(:,k) = Mapa(:,j+1)
              !write(*,*) "exto, U: ", i, j, k   
              k = k + 1

           END DO
        END IF
     END IF


     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasPasar, MapasPasar, 1, 1, DescCarga, MapasFullSkyU, 1,&
          & 1+Pasados, DescMapas, ictxt)
     CALL blacs_barrier(ICTXT, "ALL")


  END IF
  !Fin de carga los mapas rtestantes

  DEALLOCATE(Mapa)
  DEALLOCATE(MapasPasar)

  !Aplica mascara
  npix = SUM(Mascara)
  IF(MYROW==0) THEN
     ALLOCATE(MapasQUMascara(2*npix, NCMapas))
  ELSE
     ALLOCATE(MapasQUMascara(1, NCMapas))
  END IF
  MapasQUMascara = 0

  ipixQ = 1
  ipixU = npix+1
  IF(MYROW==0) THEN
     DO i=1,NpixFull
        IF(Mascara(i)==1) THEN
           MapasQUMascara(ipixQ,:) = MapasFullSkyQ(i,:)
           ipixQ = ipixQ + 1
           MapasQUMascara(ipixU,:) = MapasFullSkyU(i,:)
           ipixU = ipixU + 1 
        END IF
     END DO
  END IF

  DEALLOCATE(MapasFullSkyQ)
  DEALLOCATE(MapasFullSkyU)

  !Reparte los mapas entre los procesadores de la columna
  Mapas = -23445
  !Lo pasa todo a la vez
  CALL descinit(DescCarga, 2*NPix, NMapas, 2*NPix, NB, 0, 0, ICTXT, 2*NPix, INFO)
  CALL descinit(DescMapas, 2*NPix, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(2*NPix, NMapas, MapasQUMascara, 1, 1, DescCarga, Mapas, 1, 1, DescMapas, ictxt)
  CALL blacs_barrier(ICTXT, "ALL")

  DEALLOCATE(MapasQUMascara)

END SUBROUTINE CargaMapasFits2


SUBROUTINE CargaProblema(ICTXT, iam, MuestraMemoria, Problema)

  USE DatosProblema

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, iam
  INTEGER, INTENT(out) :: MuestraMemoria
  CHARACTER(LEN=100), INTENT(IN) :: Problema

  INTEGER, DIMENSION(20) :: Paso

  CHARACTER(len=200) :: cin
  INTEGER :: Leer

  CHARACTER(LEN=100) :: Texto

  EXTERNAL :: blacs_barrier, igebs2d, igebr2d, dgebs2d, dgebr2d

  IF(iam==0) THEN

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
           CALL CargaEntero(cin, "Type_of_Data = ", TipoFileMaps)
           CALL CargaEntero(cin, "Pixel_Window = ", PixelWindowEnMapas)
           !call CargaEntero(cin, "Intensity_Mask_NumMap = ", NumMaskMapIntensity)
           CALL CargaEntero(cin, "Polarization_Mask_NumMap = ", NumMaskMapPolarization)
           CALL CargaEntero(cin, "Type_of_Noise = ", TipoRuido)
           CALL CargaEntero(cin, "Type_of_Noise_Data = ", TipoDatosMapaRuido)
           CALL CargaEntero(cin, "Remove_Noise_Bias = ", QuitarSesgo)
           CALL CargaEntero(cin, "Matrices_Cyclic_Block_Size = ", NB)
           CALL CargaEntero(cin, "Inverse_Covariance_Matrix_Control = ", ControlInversaMC)
           CALL CargaEntero(cin, "Show_Memory_Allocated = ", MuestraMemoria)
           CALL CargaEntero(cin, "Type_of_Bin_Center = ", TypeOfBinCenter)
           CALL CargaEntero(cin, "Type_of_Grouping = ", TypeOfGrouping)
           CALL CargaEntero(cin, "Compute_Fisher_Matrix = ", ComputeFisherMatrix)
           CALL CargaEntero(cin, "Binned = ", Binned)


           !Real
           !call CargaReal(cin, "Intensity_Noise = ", RuidoTT)
           CALL CargaReal(cin, "Polarization_Noise = ", RuidoQQ)
           CALL CargaReal(cin, "Beam_FWHM = ", FWHM_Beam)

           !Cadenas
           CALL CargaCadena(cin, "Data_Folder = ", Lugar)
           CALL CargaCadena(cin, "Maps_FileName = ", FileMaps)
           CALL CargaCadena(cin, "Maps2Cross_FileName = ", FileMapsCross)
           CALL CargaCadena(cin, "Fiducial_FileName = ", ArchivoFiducial)
           !call CargaCadena(cin, "Intensity_Mask_FileName = ", FileMaskIntensity)
           CALL CargaCadena(cin, "Polarization_Mask_FileName = ", FileMaskPolarization)
           CALL CargaCadena(cin, "Noise_Map_FileName = ", FileMapRuido)
           CALL CargaCadena(cin, "Healpix_Data_Folfer = ", DirHealpixData)
           CALL CargaCadena(cin, "Binnig_Limits_FileName = ", FileBinLimits)


        END IF

        GOTO 200
100     Leer = 0
200     CONTINUE

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
     Paso(12) = TypeOfBinCenter
     Paso(13) = TypeOfGrouping
     Paso(14) = ComputeFisherMatrix
     Paso(15) = Binned

  END IF


  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN
     CALL igebs2d(ICTXT, "All","I" , 20, 1, Paso, 20)
  ELSE
     CALL igebr2d(ICTXT, "All","I" , 20, 1, Paso, 20, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , 1, 1, RuidoTT, 1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , 1, 1, RuidoTT, 1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , 1, 1, RuidoQQ, 1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , 1, 1, RuidoQQ, 1, 0, 0)

     NSide =  Paso(1)
     LMax = Paso(2)
     DlMax =  Paso(3) 
     TipoAnalisis = Paso(4)
     NMapas = Paso(5)
     QuitarSesgo = Paso(6)
     NB = Paso(7)
     ControlInversaMC =   Paso(8)
     TipoRuido = Paso(9)
     TipoDatosMapaRuido = Paso(10)
     MuestraMemoria = Paso(11)
     TypeOfBinCenter = Paso(12)
     TypeOfGrouping = Paso(13)
     ComputeFisherMatrix = Paso(14)
     Binned = Paso(15)

  END IF

  CALL blacs_barrier(ICTXT, "ALL")


  !*********************
  !    300 Format (1X, A, T82, I8)
  !    310 Format (1X, A, T75 F15.6)
  !    320 Format (1X, A, T30, A60)

300 FORMAT (1X,A35,TR2,I)
310 FORMAT (1X,A35,TR2,F)
320 FORMAT (1X,A35,TR2,A)

  IF(iam==0) THEN

     WRITE(*,320) "Data_Folder: ",TRIM(Lugar)
     WRITE(TEXTO,*) NSide
     WRITE(*,320) "NSide: ", ADJUSTL(TRIM(Texto))
     WRITE(*,320) "Fiducial_FileName: ", TRIM(ArchivoFiducial)
     WRITE(TEXTO,*) Lmax
     WRITE(*,320) "Lmax_Covariance_Matrix: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) Dlmax
     WRITE(*,320) "Lmax_Power_Spectrum: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) TipoAnalisis
     WRITE(*,320) "Type_of_Analysis: ", ADJUSTL(TRIM(Texto))
     WRITE(*,*)

     WRITE(*,320) "Maps_FileName: ", TRIM(FileMaps)

     IF(TipoAnalisis==1) THEN
        WRITE(*,320) "Maps2Cross_FileName: ", TRIM(FileMapsCross)
     END IF

     WRITE(TEXTO,*) NMapas
     WRITE(*,320) "Number_of_Maps: ", ADJUSTL(TRIM(Texto))
     !write(*,300) "Type_of_Data: ", TipoFileMaps

     WRITE(TEXTO,*) PixelWindowEnMapas
     WRITE(*,320) "Pixel_Window: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) FWHM_Beam
     WRITE(*,320) "Beam_FWHM: ", ADJUSTL(TRIM(Texto))
     WRITE(*,*)

     !write(*,320) "Intensity_Mask_FileName: ", Trim(FileMaskIntensity)
     !write(TEXTO,*) NumMaskMapIntensity
     !write(*,320) "Intensity_Mask_NumMap: ", Adjustl(Trim(Texto))
     WRITE(*,320) "Polarization_Mask_FileName: ", TRIM(FileMaskPolarization)
     WRITE(TEXTO,*) NumMaskMapPolarization
     WRITE(*,320) "Polarization_Mask_NumMap: ", ADJUSTL(TRIM(Texto))
     WRITE(*,*)


     WRITE(TEXTO,*) TipoRuido
     WRITE(*,320) "Type_of_Noise: ", ADJUSTL(TRIM(Texto))
     IF(TipoRuido==0) THEN
        !write(TEXTO,*) RuidoTT
        !write(*,320) "Intensity_Noise:", Adjustl(Trim(Texto))
        WRITE(TEXTO,*) RuidoQQ
        WRITE(*,320) "Polarization_Noise: ", ADJUSTL(TRIM(Texto))
     ELSE
        WRITE(*,320) "Noise_Map_FileName: ", TRIM(FileMapRuido)
        !write(*,300) "Type_of_Noise_Data: ", TipoDatosMapaRuido
     END IF
     WRITE(TEXTO,*) QuitarSesgo
     WRITE(*,320) "Remove_Noise_Bias: ", ADJUSTL(TRIM(Texto))
     WRITE(*,*)


     WRITE(TEXTO,*) NB
     WRITE(*,320) "Matrices_Cyclic_Block_Size: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) ControlInversaMC
     WRITE(*,320) "Inverse_Covariance_Matrix_Control: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) MuestraMemoria
     WRITE(*,320) "Show_Memory_Allocated: ", ADJUSTL(TRIM(Texto))
     WRITE(*,320) "Healpix_Data_Folfer: ", TRIM(DirHealpixData)

  END IF

  !*********************

  CALL blacs_barrier(ICTXT, "ALL")

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


SUBROUTINE CalculaBineado(NSide, npix, lmax, RuidoDiagonal, DlEE, DlBB, DlEB, BeamPWDl,&
     & TipoCompactado, TipoCentrado, QuitarSesgo, NMapas, Lugar, Filename)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: lmax, NSide, npix, QuitarSesgo, NMapas
  INTEGER, INTENT(in) :: TipoCompactado                                    !Compacta Fisher pesando o no por fiducial
  INTEGER, INTENT(in) :: TipoCentrado                                      !\ell del bin pesado por errores o \ell medio   
  REAL(Kind=8), INTENT(in) :: RuidoDiagonal
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: DlEE, DlBB, DlEB
  REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
  CHARACTER(len=100), INTENT(in) :: Lugar, Filename


  !Internas
  CHARACTER(len=100) :: File


  INTEGER :: l, i, li, ls, Caso, f, c, dim
  REAL(Kind=8) :: fsky, temp

  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: ErrorDlEE, ErrorDlBB, ErrorDlEB

  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: DLBin, LBin, Factores, MatrizR
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: SumaPesos, SumasDl, SumasL
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: PesosCompactado

  REAL(Kind=8), ALLOCATABLE, DIMENSION(:,:) :: Yl, MatrizProducto, MFBines, YlBines, Dl
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: Bl, MF
  REAL(Kind=8), ALLOCATABLE, DIMENSION(:) :: DlMedio, DlSigma, SumaDl, SumaCuardadoDl, ErrorTeoricoDl

  REAL(Kind=8) :: npixfull

  INTEGER :: Leer, Limite, NBines, info
  INTEGER, DIMENSION(1000) :: LimiteSuperiorBines

  EXTERNAL :: dsymm, dgemm, dpotrf, dpotri


  File = TRIM(Lugar)//"/"//TRIM(Filename)
  OPEN(unit=10, file=TRIM(File), action="read")
  NBines = 0

  Leer = 1
  DO WHILE(Leer==1)

     READ(10, *, END=100) Limite
     IF(Limite<lmax) THEN
        NBines = NBines + 1
        LimiteSuperiorBines(NBines) = Limite
     ELSE
        NBines = NBines + 1
        LimiteSuperiorBines(NBines) = Lmax
        Leer = 0
     END IF

     GOTO 200
100  Leer = 0
200  CONTINUE
  END DO

  CLOSE(10)


  ALLOCATE(ErrorDlEE(2:lmax))
  ALLOCATE(ErrorDlBB(2:lmax))
  ALLOCATE(ErrorDlEB(2:lmax))

  fsky = (1.0D0*npix)/(12D0*NSide**2)
  npixfull = 12D0*nside**2

  IF(TipoCentrado==1) THEN

     WRITE(*,*) 
     WRITE(*,*) " Noise per pixel: ", RuidoDiagonal

     DO l=2,lMax
        ErrorDlEE(l) =  ErrorTeorico(1, l, fsky, npixfull, DlEE, DlBB, DlEB, BeamPWDl, RuidoDiagonal)
        ErrorDlBB(l) =  ErrorTeorico(2, l, fsky, npixfull, DlEE, DlBB, DlEB, BeamPWDl, RuidoDiagonal)
        ErrorDlEB(l) =  ErrorTeorico(3, l, fsky, npixfull, DlEE, DlBB, DlEB, BeamPWDl, RuidoDiagonal)
     END DO
  ELSE
     ErrorDlEE(:) = 1D0
     ErrorDlBB(:) = 1D0
     ErrorDlEB(:) = 1D0
  END IF

  ALLOCATE(DLBin(3, NBines))
  ALLOCATE(LBin(3, NBines))
  ALLOCATE(SumaPesos(3))
  ALLOCATE(SumasDl(3))
  ALLOCATE(SumasL(3))
  ALLOCATE(PesosCompactado(3, 2:lmax))

  li = 2
  DO i=1,NBines

     ls = LimiteSuperiorBines(i)

     SumaPesos(:) = 0d0
     SumasDl(:) = 0d0
     SumasL(:) = 0d0

     DO l=li,ls

        temp = 1.0d0/ErrorDlEE(l)**2
        SumaPesos(1) = SumaPesos(1) + temp
        SumasDl(1) = SumasDl(1) + DlEE(l) * temp
        SumasL(1) = SumasL(1) + l * temp
        PesosCompactado(1,l) = temp

        temp = 1.0d0/ErrorDlBB(l)**2
        SumaPesos(2) = SumaPesos(2) + temp
        SumasDl(2) = SumasDl(2) + DlBB(l) * temp
        SumasL(2) = SumasL(2) + l * temp
        PesosCompactado(2,l) = temp

        temp = 1.0d0/ErrorDlEB(l)**2
        SumaPesos(3) = SumaPesos(3) + temp
        SumasDl(3) = SumasDl(3) + DlEB(l) * temp
        SumasL(3) = SumasL(3) + l * temp
        PesosCompactado(3,l) = temp

     END DO

     SumasDl = SumasDl/SumaPesos
     SumasL = SumasL/SumaPesos

     DLBin(:,i) = SumasDl(:)
     LBin(:,i)  = SumasL(:)

     PesosCompactado(1,li:ls) = PesosCompactado(1,li:ls)/SumaPesos(1)
     PesosCompactado(2,li:ls) = PesosCompactado(2,li:ls)/SumaPesos(2)
     PesosCompactado(3,li:ls) = PesosCompactado(3,li:ls)/SumaPesos(3)

     !-----------------
     li = ls + 1

  END DO



  DEALLOCATE(ErrorDlEE)
  DEALLOCATE(ErrorDlBB)
  DEALLOCATE(ErrorDlEB)

  OPEN(Unit=25, File=TRIM(Lugar)//"/Positions.dat", Action="write")
  DO Caso=1,3
     DO i=1,NBines
        WRITE(25,*) LBin(Caso, i)
     END DO
  END DO
  CLOSE(25)

  OPEN(Unit=25, File=TRIM(Lugar)//"/BinnedFiducial.dat", Action="write")
  DO Caso=1,3
     DO i=1,NBines
        WRITE(25,*) DLBin(Caso, i)
     END DO
  END DO
  CLOSE(25)

  !Calcula los factores para pasar de Dl a DLBin
  ALLOCATE(Factores(3,2:lmax))
  Factores = 0

  IF(TipoCompactado==1) THEN
     li = 2
     DO i=1,NBines

        ls = LimiteSuperiorBines(i)
        DO l=li,ls

          IF(DlBin(1,i)/=0) THEN
           Factores(1,l) = DlEE(l)/DlBin(1,i)
          ELSE
           Factores(1,l) = 1.0D0
          END IF

          IF(DlBin(2,i)/=0) THEN
           Factores(2,l) = DlBB(l)/DlBin(2,i)
          ELSE
           Factores(2,l) = 1.0D0
          END IF

          IF(DlBin(3,i)/=0) THEN
              Factores(3,l) = DlEB(l)/DlBin(3,i)
          ELSE
              Factores(3,l) = 1D0
          END IF

        END DO
        li = ls + 1
     END DO
  ELSE
     li = 2
     DO i=1,NBines
        ls = LimiteSuperiorBines(i)
        DO l=li,ls
           Factores(1,l) = 1D0
           Factores(2,l) = 1D0
           Factores(3,l) = 1D0
        END DO
        li = ls + 1
     END DO

  END IF

  !Construye la matriz de los Factores
  ALLOCATE(MatrizR(3*(lmax-1),3*NBines))

  MatrizR = 0.0D0
  f = 1
  c = 1
  DO Caso = 1,3
     li = 2
     DO i=1,NBines
        ls = LimiteSuperiorBines(i)
        DO l=li,ls
           MatrizR(c,f) =  Factores(Caso,l)
           c = c + 1
        END DO
        li = ls+1
        f = f + 1
     END DO
  END DO

  DEALLOCATE(Factores)

  !Carga la matriz de Fisher
  dim = 3*(lmax-1)
  ALLOCATE(MF(dim*dim))
  WRITE(*,*) " Loading Fisher matrix"
  OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrix.dat", Action="read")
  READ (25, *)  MF
  CLOSE(25)

  !Compacta la matriz de Fisher
  ALLOCATE(MatrizProducto(dim,3*NBines))  
  CALL dsymm("L", "U", dim, 3*NBines, 1.0D0, MF, dim, MatrizR, dim, 0.0D0, MatrizProducto, dim)

  ALLOCATE(MFBines(3*Nbines, 3*Nbines))
  CALL dgemm("T","N", 3*NBines, 3*NBines, dim, 1.0D0, MatrizR, dim, MatrizProducto, dim, 0.0D0, MFBines, 3*NBines)

  DEALLOCATE(MatrizProducto)

  OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrixBinnedDl.dat", Action="write")
  WRITE(25,*)  MFBines
  CLOSE(25)

  !Carga Bl
  WRITE(*,*) " Loading Noise Bias"
  ALLOCATE(Bl(dim))
  OPEN(Unit=25, File=TRIM(Lugar)//"/NoiseBias.dat", Action="read")
  READ (25, *)  Bl
  CLOSE(25)

  !Carga Yl y resta Bl
  WRITE(*,*) " Loading Coupled Power"
  ALLOCATE(Yl(dim,NMapas))
  OPEN(Unit=25, File=TRIM(Lugar)//"/CoupledPower.dat", Action="read")
  READ (25, *)  Yl
  CLOSE(25)
   
  IF(QuitarSesgo==1) THEN
     DO i=1, NMapas     
        Yl(:,i) = Yl(:,i) - Bl
     END DO
  END IF

  !Invierte la matriz de Fisher compactada
  WRITE(*,*) " Inverting the compacted Fisher matrix"

  info = 0
  dim = 3*NBines
  CALL dpotrf("U", dim, MFBines, dim, info)
  !IF(IAM==0) THEN
  WRITE(*,*) " Cholesky factorization: ", info
  !END IF

  IF(info/=0) THEN
     WRITE(*,*)
     WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     WRITE(*,*) "Fisher matrix is singular"
     RETURN
  END IF

  CALL dpotri("U", dim, MFBines, dim, info)
  !IF(IAM==0) THEN
  WRITE(*,*) " Inversion result:       ", info
  !END IF

  IF(info/=0) THEN
     WRITE(*,*)
     WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
     WRITE(*,*) "Fisher matrix is singular"
     RETURN
  END IF

  !Compacta el vertor Yl
  ALLOCATE(YlBines(dim, NMapas))
  CALL dgemm("T","N", dim, NMapas, 3*(lmax-1), 1.0D0, MatrizR, 3*(lmax-1), Yl, 3*(lmax-1), 0.0D0, YlBines, dim)

  !Multiplica por la inversa de MF
  WRITE(*,*) " Computing Dl"
  ALLOCATE(Dl(dim, NMapas))
  CALL dsymm("L", "U", dim, NMapas, 1.0D0, MFBines, dim, YlBines, dim, 0.0D0, Dl, dim)

  !Guarda
  WRITE(*,*) " Saving Dl"
  OPEN(Unit=25, File=TRIM(Lugar)//"/BinnedDl.dat", Action="write")
  WRITE (25, *)  Dl
  CLOSE(25)

  IF(NMapas>1) THEN

     !Medias y desviacion estandar Dl
     ALLOCATE(DlMedio(dim))
     ALLOCATE(DlSigma(dim))

     ALLOCATE(SumaDl(dim))
     ALLOCATE(SumaCuardadoDl(dim))

     DlMedio = 0.0D0
     DlSigma = 0.0D0
     SumaDl = 0.0D0
     SumaCuardadoDl = 0.0D0


     DO i=1,NMapas
        SumaDl = SumaDl + Dl(:,i)
        SumaCuardadoDl = SumaCuardadoDl + Dl(:,i)**2
     END DO

     DlMedio = SumaDl/NMapas
     DlSigma = SQRT(SumaCuardadoDl/NMapas-DlMedio**2)
     DlSigma = SQRT((1.0D0 * NMapas)/(1.0D0 * (NMapas-1))) * DlSigma

     !Guarda
     OPEN(Unit=25, File=TRIM(Lugar)//"/MeanBinnedDl.dat", Action="write")
     WRITE (25, *)  DlMedio
     CLOSE(25)

     OPEN(Unit=25, File=TRIM(Lugar)//"/SigmaBinnedDl.dat", Action="write")
     WRITE (25, *)  DlSigma
     CLOSE(25)

     DEALLOCATE(DlMedio)
     DEALLOCATE(DlSigma)
     DEALLOCATE(SumaDl)
     DEALLOCATE(SumaCuardadoDl)

  END IF


  !Guarda los errores teóricos
  ALLOCATE(ErrorTeoricoDl(dim))
  DO i=1,dim
     ErrorTeoricoDl(i) = SQRT(MFBines(i,i))
  END DO
  OPEN(Unit=25, File=TRIM(Lugar)//"/FisherErrorBinnedDl.dat", Action="write")
  WRITE(25,*) ErrorTeoricoDl
  CLOSE(25)
  DEALLOCATE(ErrorTeoricoDl)        

  RETURN

CONTAINS

  REAL(Kind=8) FUNCTION ErrorTeorico(Caso, l, fsky, npixfull, DlEE, DlBB, DlEB, BeamPWDl, RuidoDiagonal)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: Caso, l
    REAL(Kind=8), INTENT(in) :: RuidoDiagonal, fsky, npixfull
    REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: DlEE, DlBB, DlEB
    REAL(Kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl

    !Internas
    REAL(Kind=8) error
    REAL(Kind=8), PARAMETER :: DosPi = 6.2831853071796D0


    SELECT CASE (Caso)
    CASE (1)
       error = DlEE(l) + (2D0*DosPi/npixfull)*RuidoDiagonal**2/BeamPWDl(l)**2
       error = error**2
       error = 2.0D0/(2.0*l+1.0) * 1.0D0/fsky * error
       ErrorTeorico = SQRT(error)
    CASE (2)
       error = DlBB(l) + (2D0*DosPi/npixfull) * RuidoDiagonal**2/BeamPWDl(l)**2
       error = error**2
       error = 2.0D0/(2.0*l+1.0) * 1.0D0/fsky * error
       ErrorTeorico = SQRT(error)
    CASE (3)
       error = DlEB(l)**2 + (DlEE(l) + (2D0*DosPi/npixfull) * RuidoDiagonal**2/BeamPWDl(l)**2) *&
            & (DlBB(l) + (2D0*DosPi/npixfull) * RuidoDiagonal**2/BeamPWDl(l)**2)
       error = 1.0D0/(2.0*l+1.0) * 1.0D0/fsky * error
       ErrorTeorico = SQRT(error)
    END SELECT

    RETURN

  END FUNCTION ErrorTeorico

END SUBROUTINE CalculaBineado
