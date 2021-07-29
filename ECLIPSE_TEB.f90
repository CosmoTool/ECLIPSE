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
    IF(Mostrar<0.1) Mostrar = 0 


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
  REAL(kind=8)    :: RuidoTT, RuidoQQ
  REAL(kind=8)    :: FWHM_Beam
  INTEGER         :: PixelWindowEnMapas
  INTEGER         :: KindOfGrouping
  INTEGER         :: KindOfBinCenter
  INTEGER         :: ComputeFisherMatrix
  INTEGER         :: Binned
  INTEGER         :: ComputeSpectrum

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

PROGRAM ECLIPSE_TEB

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

     SUBROUTINE PuntosPixelesObservados(ictxt, nside, npix, Mascara, NB, NPixProc, z, a)
       INTEGER, INTENT(in) :: ictxt, nside, npix, NB, NPixProc
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
       REAL(kind=8), DIMENSION(NPixProc) :: z, a
     END SUBROUTINE PuntosPixelesObservados

     SUBROUTINE BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)
       INTEGER, INTENT(in) :: ICTXT, PixelWindowEnMapas, nside, lmax
       REAL(kind=8) :: FWHM_Beam
       CHARACTER(len=*), INTENT(in) :: DirHealpixData       
       REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(out) :: BeamPWDl
     END SUBROUTINE BeamPixelWindow

     SUBROUTINE CargaMultipolos(ICTXT, iam, Lugar, FileData, LMax, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB)
       INTEGER, INTENT(in) :: ICTXT, IAM, lmax
       CHARACTER(len=100), INTENT(in) :: Lugar, FileData
       REAL(kind=8), DIMENSION(0:lmax), INTENT(OUT) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
     END SUBROUTINE CargaMultipolos

     SUBROUTINE CalculaBloquesMatrizArmonicosY(ICTXT, NB, NFP, NC, npix, lmax, z, Phi, BeamPWDl, YTTR, YTTI)
       INTEGER, INTENT(IN) :: ICTXT, NB, NFP, NC, npix, lmax
       REAL(kind=8), DIMENSION(NFP), INTENT(IN)   :: z, phi
       REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(in) :: BeamPWDl
       REAL(kind=8), DIMENSION(NFP,NC), INTENT(OUT) :: YTTR, YTTI
     END SUBROUTINE CalculaBloquesMatrizArmonicosY

     SUBROUTINE CalculaBloquesMatrizArmonicosYFC(ICTXT, NB, NFP, NC, npix, lmax, z, Phi, BeamPWDl, YTTC)
       INTEGER, INTENT(IN) :: ICTXT, NB, NFP, NC, npix, lmax
       REAL(kind=8), DIMENSION(NFP), INTENT(IN)   :: z, phi
       REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(in) :: BeamPWDl
       COMPLEX(kind=8), DIMENSION(NFP,NC), INTENT(OUT) :: YTTC
     END SUBROUTINE CalculaBloquesMatrizArmonicosYFC

     SUBROUTINE CalculaBloquesMatrizArmonicosX(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, &
          &YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)
       INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
       REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
       REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
       REAL(kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
     END SUBROUTINE CalculaBloquesMatrizArmonicosX

     SUBROUTINE CalculaBloquesMatrizArmonicosXFC(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, &
          &YQEC, YQBC, YUEC, YUBC)
       INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
       REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
       REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
       COMPLEX(kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQEC, YQBC, YUEC, YUBC
     END SUBROUTINE CalculaBloquesMatrizArmonicosXFC

     SUBROUTINE CargaMapasRuido(ICTXT, TipoRuido, RuidoTT, RuidoQQ, nside, NFBY, NFBX, NFQUMC, NB,&
          & MascaraT, MascaraP, Lugar, FileMapRuido, MapaRuidoT2Bloque, MapaRuidoP2Bloque, MapaRuidoP2)
       INTEGER, INTENT(in) :: ICTXT, TipoRuido, nside,  NFBY, NFBX, NFQUMC, NB
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: MascaraT, MascaraP
       CHARACTER(LEN=100), INTENT(in) :: Lugar, FileMapRuido
       REAL(kind=8), INTENT(inout) :: RuidoTT, RuidoQQ
       REAL(kind=8), DIMENSION(NFBY), INTENT(out) :: MapaRuidoT2Bloque
       REAL(kind=8), DIMENSION(NFBX), INTENT(out) :: MapaRuidoP2Bloque
       REAL(kind=8), DIMENSION(NFQUMC), INTENT(out) :: MapaRuidoP2
     END SUBROUTINE CargaMapasRuido

     SUBROUTINE CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NFBY, NCBY, NFBX, lmax, NPixT, NPixP,&
          & MapaRuidoT2Bloque, MapaRuidoP2Bloque, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB,&
          & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, ControlInversaMC, MC, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, NB, NFMC, NCMC, NFBY, NCBY, NFBX, lmax, NPixT, NPixP
       INTEGER, INTENT(in) :: ControlInversaMC
       INTEGER(kind=8), INTENT(in) ::  TiempoInicio
       REAL(kind=8), DIMENSION(NFBY), INTENT(in) :: MapaRuidoT2Bloque
       REAL(kind=8), DIMENSION(NFBX), INTENT(in) :: MapaRuidoP2Bloque
       REAL(kind=8), DIMENSION(0:lmax), INTENT(in) ::  DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
       REAL(KIND=8), DIMENSION(NFBY,NCBY), INTENT(in)  :: YTTR, YTTI
       REAL(KIND=8), DIMENSION(NFBX,NCBY), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
       REAL(KIND=8), DIMENSION(NFMC,NCMC), INTENT(out)  :: MC
     END SUBROUTINE CalculaMatrizCovarianza

     SUBROUTINE InvierteMC(ICTXT, DimMC, NB, NF, NC, MC, ControlInversa, Fallo, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NB, NF, NC, ControlInversa
       REAL(kind=8), DIMENSION(NF,NC), INTENT(inout) :: MC
       INTEGER, INTENT(out) :: Fallo
       INTEGER(kind=8), INTENT(in) ::  TiempoInicio
     END SUBROUTINE InvierteMC

     SUBROUTINE CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, MascaraT, MascaraP, Lugar, FileData,  Mapas)
       INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas
       INTEGER, DIMENSION(12*nside**2), INTENT(in) :: MascaraT, MascaraP
       CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
       REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas
     END SUBROUTINE CargaMapasFits

     SUBROUTINE CalculaIMCMapas(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
       REAL(kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
       REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaIMCMapas

     SUBROUTINE CalculaYDagIMCMapas(ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, lmax, NMapas, NCMapas, NB,&
          & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, IMCMapas, PTR, PER, PBR, PTI, PEI, PBI)
       INTEGER, INTENT(in) :: ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, Lmax, NMapas, NCMapas, NB
       REAL(KIND=8), DIMENSION(NFBY,NCBA), INTENT(in)  :: YTTR, YTTI
       REAL(KIND=8), DIMENSION(NFBX,NCBA), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
       REAL(KIND=8), DIMENSION(NFMC,NCMapas), INTENT(in)  :: IMCMapas
       REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(out) :: PTR, PER, PBR, PTI, PEI, PBI
     END SUBROUTINE CalculaYDagIMCMapas

     SUBROUTINE CalculaYlAuto(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI, Lugar)
       INTEGER, INTENT(in) :: ICTXT, NFMapasA, Lmax, NCMapas, NMapas, NB
       REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR, MATI, MAER, MAEI, MABR, MABI
       CHARACTER(len=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaYlAuto

     SUBROUTINE CalculaYlCross(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI,&
          &MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2, Lugar)
       INTEGER, INTENT(in) :: ICTXT, NFMapasA, Lmax, NCMapas, NMapas, NB
       REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR, MATI, MAER, MAEI, MABR, MABI
       REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2
       CHARACTER(len=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaYlCross

     SUBROUTINE BloquesMatrizCovarianza(ICTXT, DimMC, NFMC, NCMC, NpixT, NpixP, NFBMCT, NCBMCT, NFBMCP, NCBMCP,&
          & NB, MC, BTT, BTQ, BTU, BQQ, BQU, BUU)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NFMC, NCMC, NpixT, NpixP, NFBMCT, NCBMCT, NFBMCP, NCBMCP, NB
       REAL(kind=8), DIMENSION(NFMC, NCMC), INTENT(in) :: MC
       REAL(kind=8), DIMENSION(NFBMCT, NCBMCT), INTENT(out) :: BTT
       REAL(kind=8), DIMENSION(NFBMCT, NCBMCP), INTENT(out) :: BTQ, BTU
       REAL(kind=8), DIMENSION(NFBMCP, NCBMCP), INTENT(out) :: BQQ, BQU, BUU
     END SUBROUTINE BloquesMatrizCovarianza

     SUBROUTINE  ProductosBloquesIMCArmonicos(ICTXT, NCBMCT, NCBMCP, BTT, BTQ, BTU, BQQ, BQU, BUU,&
          & lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP, YTT, YQE, YQB, YUE, YUB,&
          & BPTT, BPTE, BPTB, BPQT, BPQE, BPQB, BPUT, BPUE, BPUB, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, NCBMCT, NCBMCP, lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP
       !Bloques de IMC
       REAL(kind=8), DIMENSION(NFBY, NCBMCT), INTENT(in) :: BTT
       REAL(kind=8), DIMENSION(NFBY, NCBMCP), INTENT(in ) :: BTQ, BTU
       REAL(kind=8), DIMENSION(NFBX, NCBMCP), INTENT(in) :: BQQ, BQU, BUU
       !Bloques de la matriz de armonicos
       REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(in) :: YTT
       REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(in) :: YQE, YQB, YUE, YUB

       INTEGER(kind=8), INTENT(in) ::  TiempoInicio

       !Bloques producto IMC.Y
       REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(out) :: BPTT, BPTE, BPTB
       REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(out) :: BPQT, BPQE, BPQB, BPUT, BPUE, BPUB
     END SUBROUTINE  ProductosBloquesIMCArmonicos

     SUBROUTINE  CalculaBl(ICTXT, Lugar, NFBY, NFBX, NCBA, Lmax, NB, BTTR, BTER, BTBR, BQTR, BQER, BQBR, BUTR, BUER, BUBR,&
          & BTTI, BTEI, BTBI, BQTI, BQEI, BQBI, BUTI, BUEI, BUBI, RuidoTT, RuidoQQ)
       INTEGER, INTENT(in) :: ICTXT, NFBY, NFBX, NCBA, Lmax, NB
       CHARACTER(len=100), INTENT(in) :: Lugar
       REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(IN) :: BTTR, BTER, BTBR, BTTI, BTEI, BTBI
       REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(IN) :: BQTR, BQER, BQBR, BQTI, BQEI, BQBI
       REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(IN) :: BUTR, BUER, BUBR, BUTI, BUEI, BUBI
       REAL(kind=8), DIMENSION(NFBY), INTENT(IN) :: RuidoTT
       REAL(kind=8), DIMENSION(NFBX), INTENT(IN) :: RuidoQQ 
     END SUBROUTINE CalculaBl

     SUBROUTINE CalculaMatrizArmonicosCompleja(ICTXT, NB, NF, NFMX, NCMX, npix, lmax, z, Phi, BeamPWDl, MXC)
       INTEGER, INTENT(IN) :: ICTXT, NB, NF, NFMX, NCMX, npix, lmax
       REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
       REAL(kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
       COMPLEX(kind=8), DIMENSION(NFMX,NCMX), INTENT(OUT) :: MXC
     END SUBROUTINE CalculaMatrizArmonicosCompleja

     SUBROUTINE CalculaMatrizFisher(ICTXT, NFBA, NCBA, lmax, NB, TTR, TTI, TER, TEI, TBR,&
          & TBI, EER, EEI, BBR, BBI, EBR, EBI, lugar, TiempoInicio)
       INTEGER, INTENT(in) :: ICTXT, NFBA, NCBA, lmax, NB
       INTEGER(kind=8), INTENT(in) ::  TiempoInicio
       REAL(kind=8), DIMENSION(NFBA, NCBA), INTENT(in) ::  TTR, TTI, TER, TEI, TBR, TBI, EER, EEI, BBR, BBI, EBR, EBI
       CHARACTER(len=100), INTENT(in) :: lugar
     END SUBROUTINE CalculaMatrizFisher

     SUBROUTINE CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)
       INTEGER(kind=4), INTENT(in) :: IAM, lmax, NMapas, QuitarSesgo
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaEspectroPotencia

     SUBROUTINE CalculaIMCMapasMapasCross(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, IMC, Mapas, MapasCross)
       INTEGER, INTENT(in) :: ICTXT, DimMC, NMapas, NB, NFMC, NCMC, NCMapas
       REAL(kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
       REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas, MapasCross
       CHARACTER(LEN=100), INTENT(in) :: Lugar
     END SUBROUTINE CalculaIMCMapasMapasCross

     SUBROUTINE CalculaBineado(NSide, npixT, npixP, lmax, RuidoTT, RuidoQQ, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB, BeamPWDl,&
          & TipoCompactado, TipoCentrado, QuitarSesgo, NMapas, Lugar, FileBinLimits)
       INTEGER, INTENT(in) :: lmax, NSide, npixT, npixP, QuitarSesgo, NMapas
       INTEGER, INTENT(in) :: TipoCompactado                                    
       INTEGER, INTENT(in) :: TipoCentrado            
       REAL(kind=8), INTENT(in) :: RuidoTT, RuidoQQ
       REAL(kind=8), DIMENSION(0:lmax), INTENT(in) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
       REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
       CHARACTER(len=100), INTENT(in) :: Lugar, FileBinLimits
     END SUBROUTINE CalculaBineado

  END INTERFACE

  !Grid
  INTEGER :: IAM, NPROCS, NPROW, NPCOL, ICTXT, MYROW, MYCOL

  !Tiempos
  INTEGER(kind=8) :: TiempoInicio, TiempoFinal

  INTEGER, EXTERNAL :: NUMROC, INDXG2L, INDXG2P
  EXTERNAL :: BLACS_PINFO, BLACS_GET, BLACS_GRIDINIT, BLACS_GRIDINFO
  EXTERNAL :: blacs_barrier, BLACS_GRIDEXIT, BLACS_EXIT, pdgemr2d, descinit
  EXTERNAL :: pdsymm, pzgemm, pdgemm

  !Mascara
  INTEGER :: NPixT, NPixP
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MascaraT, MascaraP

  !Pixeles observados
  INTEGER :: NFBY, NFBX
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: zt, at, zp, ap

  !Matrices de armonicos esfericos
  INTEGER :: TotalArmonicos, NCBA
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI

  !Multipolos
  REAL(kind = 8), ALLOCATABLE, DIMENSION(:) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
  REAL(kind = 8), ALLOCATABLE, DIMENSION(:,:) :: BeamPWDl

  !Matriz de covarianza
  INTEGER :: DimMC, NFMC, NCMC, NFQUMC, Fallo
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MC

  !Mapas
  INTEGER :: NCMapas
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Mapas, MapasCross
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: MapaRuidoT2Bloque, MapaRuidoP2Bloque
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: MapaRuidoP2

  !Mapas en armonicos
  INTEGER :: NFMapasA
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MATR, MATI, MAER, MAEI, MABR, MABI
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2

  !Bloques de la matriz de covarianza
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BTT, BTQ, BTU, BQQ, BQU, BUU
  INTEGER :: NFBMCT, NFBMCP, NCBMCT, NCBMCP

  !Bloques Producto IMC.Armonicos
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPTTR, BPTER, BPTBR
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPQTR, BPQER, BPQBR
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPUTR, BPUER, BPUBR
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPTTI, BPTEI, BPTBI
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPQTI, BPQEI, BPQBI
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPUTI, BPUEI, BPUBI

  !Bloques Producto IMC.Armonicos en forma compleja
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPTTC, BPTEC, BPTBC
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPQEC, BPQBC
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BPUEC, BPUBC

  !Para la matriz de Fisher
  !Bloques de armonicos en forma compleja
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:,:) :: YTTC, YQEC, YQBC, YUEC, YUBC
  !Bloques de Y^dag.IMC.Y
  COMPLEX(kind=8), ALLOCATABLE, DIMENSION(:,:) :: TT, TE, TB, EE, BB, EB
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: TTR, TER, TBR, EER, BBR, EBR
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: TTI, TEI, TBI, EEi, BBI, EBI
  COMPLEX(kind=8) :: Alpha, Beta
  INTEGER :: NFBA, INFO
  INTEGER, DIMENSION(9) :: DescBYIMCY, DescBIMCY, DescBY

  INTEGER :: fp, cp, Temp_f, Temp_c

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
     WRITE(*,100) "Loading masks"
     WRITE(*,*)
  END IF
  ALLOCATE(MascaraT(12*NSide**2))
  CALL CargaMascara(ictxt, nside, Lugar, FileMaskIntensity, NumMaskMapIntensity, MascaraT)
  NPixT = SUM(MascaraT)

  IF(iam==0) THEN
     WRITE(*,'(A,T82,I8)') "   Number of observed pixels in temperature: ", NPIXT
  END IF

  ALLOCATE(MascaraP(12*NSide**2))
  CALL CargaMascara(ictxt, nside, Lugar, FileMaskPolarization, NumMaskMapPolarization, MascaraP)
  NPixP = SUM(MascaraP)

  IF(iam==0) THEN
     WRITE(*,'(A,T82,I8)') "   Number of observed pixels in polarization: ", NPIXP
  END IF

  !Posiciones pixeles observados en temperatura
  NFBY = MAX(1,NUMROC(NPixT, NB, MYROW, 0, NPROW))
  ALLOCATE(zt(NFBY))
  ALLOCATE(at(NFBY))
  at = 0
  zt = 0  
  CALL PuntosPixelesObservados(ictxt, nside, NPixT, MascaraT, NB, NFBY, zt, at)

  !Posiciones pixeles observados en polarizacion
  NFBX = MAX(1,NUMROC(NPixP, NB, MYROW, 0, NPROW))
  ALLOCATE(zp(NFBX))
  ALLOCATE(ap(NFBX))
  ap = 0
  zp = 0
  CALL PuntosPixelesObservados(ictxt, nside, NPixP, MascaraP, NB, NFBX, zp, ap)


  !Carga multipolos
  ALLOCATE(DlTT(0:lmax))
  ALLOCATE(DlEE(0:lmax))
  ALLOCATE(DlBB(0:lmax))
  ALLOCATE(DlTE(0:lmax))
  ALLOCATE(DlTB(0:lmax))
  ALLOCATE(DlEB(0:lmax))
  ALLOCATE(BeamPWDl(3,0:lmax))

  IF((iam==0).AND.(PixelWindowEnMapas==1)) THEN
     WRITE(*,*)
     WRITE(*,100) "Loading HEALPix Pixel Window"
  END IF
  CALL BeamPixelWindow(ICTXT, PixelWindowEnMapas, nside, lmax, FWHM_Beam, DirHealpixData, BeamPWDl)

  CALL blacs_barrier(ICTXT, "ALL")

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,100) "Loading Fiducial Power Spectrum"
  END IF

  CALL CargaMultipolos(ICTXT, iam, Lugar, ArchivoFiducial, LMax, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB)

  !***************************************************************************************
  !Calcula los bloques de la matriz de armonicos esfericos

  IF(iam==0) THEN
     WRITE(*,110)
     WRITE(*,115)
     WRITE(*,100) "Step 1. Computing covariance matrix"
     WRITE(*,*)
     WRITE(*,100) "  Computing blocks of the spherical harmonics matrix ", Time()-TiempoInicio
  END IF

  TotalArmonicos = -3 + 2 * lmax + lmax*lmax

  !Los bloques tienen TotalArmonicos columnas, no 2*TotalArmonicos
  NCBA = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))

  ALLOCATE(YTTR(NFBY, NCBA))
  ALLOCATE(YTTI(NFBY, NCBA))

  CALL SumaMemoria(ictxt, iam, 0, NFBY, NCBA, 2)
  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBX, NCBA, 8)

  CALL CalculaBloquesMatrizArmonicosY(ICTXT, NB, NFBY, NCBA, NPixT, lmax, zt, at, BeamPWDl, YTTR, YTTI)

  IF(iam==0) THEN
     WRITE(*,100) "   Block YTT done ", Time()-TiempoInicio
  END IF

  ALLOCATE(YQER(NFBX, NCBA))
  ALLOCATE(YQEI(NFBX, NCBA))
  ALLOCATE(YQBR(NFBX, NCBA))
  ALLOCATE(YQBI(NFBX, NCBA))
  ALLOCATE(YUER(NFBX, NCBA))
  ALLOCATE(YUEI(NFBX, NCBA))
  ALLOCATE(YUBR(NFBX, NCBA))
  ALLOCATE(YUBI(NFBX, NCBA))


  CALL CalculaBloquesMatrizArmonicosX(ICTXT, NB, NFBX, NCBA, NPixP, lmax, zp, ap, BeamPWDl, &
       & YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)       

  IF(iam==0) THEN
     WRITE(*,100) "   Block YPP done ", Time()-TiempoInicio
  END IF


  IF((iam==0).AND.(TipoRuido==1)) THEN
     WRITE(*,*)
     WRITE(*,100) "  Loading noise maps"
  END IF

  !Mapas de ruido en formatos bloque armonicos y bloque matriz covarianza
  ALLOCATE(MapaRuidoT2Bloque(NFBY))
  ALLOCATE(MapaRuidoP2Bloque(NFBX))
  NFQUMC = MAX(1,1,NUMROC(2*NPixP, NB, MYROW, 0, NPROW))
  ALLOCATE(MapaRuidoP2(NFQUMC))
  MapaRuidoT2Bloque = -1000
  MapaRuidoP2Bloque = -1000
  MapaRuidoP2 = -1000

  CALL blacs_barrier(ICTXT, "ALL")

  CALL CargaMapasRuido(ICTXT, TipoRuido, RuidoTT, RuidoQQ, nside, NFBY, NFBX, NFQUMC, NB,&
       & MascaraT, MascaraP, Lugar, FileMapRuido, MapaRuidoT2Bloque, MapaRuidoP2Bloque, MapaRuidoP2)


  !Matriz de covarianza
  DimMC = NPixT + 2*NPixP
  NFMC = MAX(1,1,NUMROC(DimMC, NB, MYROW, 0, NPROW))
  NCMC = MAX(1,1,NUMROC(DimMC, NB, MYCOL, 0, NPCOL))
  ALLOCATE(MC(NFMC, NCMC))
  MC = 1D0


  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,100) "  Computing blocks of the covariance matrix ", Time()-TiempoInicio
  END IF
  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMC, 1)

  !Calcula la matriz de covarianza a partir de los bloques de la matriz de armonicos

  CALL CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NFBY, NCBA, NFBX, lmax,&
       & NPixT, NPixP,&
       & MapaRuidoT2Bloque, MapaRuidoP2Bloque, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB,&
       & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, ControlInversaMC, MC,&
       & TiempoInicio)

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,*) "  Intensity diagonal element:   ", MC(1,1)
  END IF

  CALL blacs_barrier(ICTXT,"ALL")

  fp = INDXG2P(NPixT+1, NB, MYROW, 0, NPROW)
  cp = INDXG2P(NPixT+1, NB, MYCOL, 0, NPCOL)

  IF((MYROW==fp).AND.(MYCOL==cp)) THEN
     temp_f = INDXG2L(NPixT+1, NB, MYROW, 0, NPROW)
     temp_c = INDXG2L(NPixT+1, NB, MYCOL, 0, NPCOL)
     WRITE(*,*) "  Polarization diagonal element: ", MC(temp_f, temp_c)
  END IF

  CALL blacs_barrier(ICTXT,"ALL")

  IF(iam==0) THEN
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

  CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, -8)

  DEALLOCATE(YTTR)
  DEALLOCATE(YTTI)

  CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBY, NCBA, -2)

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
        !write(*,*)
        !write(*,*) "  Two elements of C^-1: ", MC(1,1), MC(1,2)
        WRITE(*,*)
        WRITE(*,100) "Covariance matrix inverted ", Time()-TiempoInicio
     END IF

     CALL SumaMemoria(ictxt, iam, MuestraMemoria, 0,0,1)

     !Fin de la parte en la que lmax esta determinado por el limite al calcular la matriz de covarianza

     !Comienza los calculos en los que lmax esta determinado por el limite al estimar el espectro

     lmax = Dlmax

     TotalArmonicos = -3 + 2 * lmax + lmax*lmax

     !Los bloques tienen TotalArmonicos columnas, no 2*TotalArmonicos
     NCBA = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))

     ALLOCATE(YQER(NFBX, NCBA))
     ALLOCATE(YQEI(NFBX, NCBA))
     ALLOCATE(YQBR(NFBX, NCBA))
     ALLOCATE(YQBI(NFBX, NCBA))
     ALLOCATE(YUER(NFBX, NCBA))
     ALLOCATE(YUEI(NFBX, NCBA))
     ALLOCATE(YUBR(NFBX, NCBA))
     ALLOCATE(YUBI(NFBX, NCBA))

     ALLOCATE(YTTR(NFBY, NCBA))
     ALLOCATE(YTTI(NFBY, NCBA))


     IF(iam==0) THEN
        WRITE(*,110)
        WRITE(*,115)
        WRITE(*,100) "Step 3. Computing coupled power in the harmonics space"
        WRITE(*,*)
        WRITE(*,100) "  Computing blocks of the spherical harmonics matrix "
     END IF

     CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, 8)
     CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBY, NCBA, 2)

     CALL blacs_barrier(ICTXT, "ALL")

     CALL CalculaBloquesMatrizArmonicosY(ICTXT, NB, NFBY, NCBA, NPixT, lmax, zt, at, BeamPWDl, YTTR, YTTI)

     IF(iam==0) THEN
        WRITE(*,100) "   Block YTT done ", Time()-TiempoInicio
     END IF

     CALL CalculaBloquesMatrizArmonicosX(ICTXT, NB, NFBX, NCBA, NPixP, lmax, zp, ap, BeamPWDl, &
          & YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)  

     CALL blacs_barrier(ICTXT, "ALL")

     IF(iam==0) THEN
        WRITE(*,100) "   Block YPP done", Time()-TiempoInicio
     END IF


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
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        CALL CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, MascaraT, MascaraP, Lugar, FileMaps, Mapas)

        IF(iam==0) THEN
           WRITE(*,100) "  Computing C^-1 Maps ", Time()-TiempoInicio
        END IF

        CALL blacs_barrier(ICTXT, 'All')

        CALL CalculaIMCMapas(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, MC, Mapas)

        !Filas por procesador de los mapas en armonicos
        NFMapasA = MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))
        ALLOCATE(MATR(NFMapasA, NCMAPAS))
        ALLOCATE(MATI(NFMapasA, NCMAPAS))
        ALLOCATE(MAER(NFMapasA, NCMAPAS))
        ALLOCATE(MAEI(NFMapasA, NCMAPAS))
        ALLOCATE(MABR(NFMapasA, NCMAPAS))
        ALLOCATE(MABI(NFMapasA, NCMAPAS))
        MATR = 0
        MATI = 0
        MAER = 0
        MAEI = 0
        MABR = 0
        MABI = 0

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMapasA, NCMapas, 6)

        CALL CalculaYDagIMCMapas(ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, lmax, NMapas, NCMapas, NB,&
             & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, Mapas, MATR, MAER, MABR, MATI, MAEI, MABI)

        DEALLOCATE(Mapas)
        CALL SumaMemoria(ictxt, iam, 0, NFMC, NCMapas, -1)

        IF(iam==0) WRITE(*,100) "  Maps transformed to harmonis space"

        CALL CalculaYlAuto(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI, Lugar)

        DEALLOCATE(MATR)
        DEALLOCATE(MATI)
        DEALLOCATE(MAER)
        DEALLOCATE(MAEI)
        DEALLOCATE(MABR)
        DEALLOCATE(MABI)

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "Auto-corrrelation power already computed", Time()-TiempoInicio
        END IF

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMapasA, NCMapas, -6)

     ELSE

        IF(iam==0) THEN
           WRITE(*,*) 
           WRITE(*,100) "Loading maps from two files"
        END IF

        ALLOCATE(Mapas(NFMC, NCMapas))
        Mapas = 0
        CALL SumaMemoria(ictxt, iam, 0, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        CALL  CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, MascaraT, MascaraP, Lugar, FileMaps, Mapas)


        ALLOCATE(MapasCross(NFMC, NCMapas))
        MapasCross = 0
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, 1)
        CALL blacs_barrier(ICTXT, 'All')
        CALL  CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, MascaraT, MascaraP, Lugar, FileMapsCross, MapasCross)

        CALL blacs_barrier(ICTXT, 'All')


        IF(iam==0) THEN
           !write(*,*) 
           WRITE(*,100) "Computing C^-1 Maps ", Time()-TiempoInicio
        END IF

        CALL CalculaIMCMapasMapasCross(ICTXT, Lugar, DimMC, NMapas, NB, NFMC, NCMC, NCMapas, MC, Mapas, MapasCross)

        CALL blacs_barrier(ICTXT, 'All')


        !Filas por procesador de los mapas en armonicos
        NFMapasA = MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))
        ALLOCATE(MATR(NFMapasA, NCMAPAS))
        ALLOCATE(MATI(NFMapasA, NCMAPAS))
        ALLOCATE(MAER(NFMapasA, NCMAPAS))
        ALLOCATE(MAEI(NFMapasA, NCMAPAS))
        ALLOCATE(MABR(NFMapasA, NCMAPAS))
        ALLOCATE(MABI(NFMapasA, NCMAPAS))

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMapasA, NCMapas, 6)

        CALL CalculaYDagIMCMapas(ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, lmax, NMapas, NCMapas, NB,&
             & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, Mapas, MATR, MAER, MABR, MATI, MAEI, MABI)

        DEALLOCATE(Mapas)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, -1)

        IF(iam==0) WRITE(*,*) " Maps transformed to harmonic space"

        ALLOCATE(MATR2(NFMapasA, NCMAPAS))
        ALLOCATE(MATI2(NFMapasA, NCMAPAS))
        ALLOCATE(MAER2(NFMapasA, NCMAPAS))
        ALLOCATE(MAEI2(NFMapasA, NCMAPAS))
        ALLOCATE(MABR2(NFMapasA, NCMAPAS))
        ALLOCATE(MABI2(NFMapasA, NCMAPAS))

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMapasA, NCMapas, 6)

        CALL CalculaYDagIMCMapas(ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, lmax, NMapas, NCMapas, NB,&
             & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, MapasCross, MATR2, MAER2, MABR2, MATI2, MAEI2, MABI2)

        DEALLOCATE(MapasCross)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMapas, -1)

        IF(iam==0) WRITE(*,*) " Maps to cross transformed to harmonic space"

        CALL CalculaYlCross(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI,&
             & MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2, Lugar)

        DEALLOCATE(MATR)
        DEALLOCATE(MATI)
        DEALLOCATE(MAER)
        DEALLOCATE(MAEI)
        DEALLOCATE(MABR)
        DEALLOCATE(MABI)
        DEALLOCATE(MATR2)
        DEALLOCATE(MATI2)
        DEALLOCATE(MAER2)
        DEALLOCATE(MAEI2)
        DEALLOCATE(MABR2)
        DEALLOCATE(MABI2)

        IF(iam==0) WRITE(*,100) " Cross-corrrelation power already computed ", Time()-TiempoInicio

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMapasA, NCMapas, -12)

     END IF !Auto o cross-crorrelation

     !Fin carga mapas y calcula Yl
     !************************************************************************************

     IF(ComputeFisherMatrix==1) THEN

        !Aqui puede liberar memoria de los bloques de la matriz de armonicos


        IF(iam==0) THEN
           WRITE(*,110)
           WRITE(*,115)
           WRITE(*,100) "Step 4. Computing product C^-1 Y"
           WRITE(*,*)  
           WRITE(*,100) " Moving from matrix C^-1 to blocks of the matrix C^-1 ", Time()-TiempoInicio
        END IF

        NFBMCT = MAX(1,NUMROC(NpixT, NB, MYROW, 0, NPROW))
        NCBMCT = MAX(1,NUMROC(NpixT, NB, MYCOL, 0, NPCOL))
        NFBMCP = MAX(1,NUMROC(NpixP, NB, MYROW, 0, NPROW))
        NCBMCP = MAX(1,NUMROC(NpixP, NB, MYCOL, 0, NPCOL))

        ALLOCATE(BTT(NFBMCT, NCBMCT))
        ALLOCATE(BTQ(NFBMCT, NCBMCP))
        ALLOCATE(BTU(NFBMCT, NCBMCP))
        ALLOCATE(BQQ(NFBMCP, NCBMCP))
        ALLOCATE(BQU(NFBMCP, NCBMCP))
        ALLOCATE(BUU(NFBMCP, NCBMCP))

        CALL SumaMemoria(ictxt, iam, 0, NFBMCT, NCBMCT, 1)
        CALL SumaMemoria(ictxt, iam, 0, NFBMCT, NCBMCP, 2)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBMCP, NCBMCP, 3)

        CALL BloquesMatrizCovarianza(ICTXT, DimMC, NFMC, NCMC, NpixT, NpixP, NFBMCT, NCBMCT, NFBMCP, NCBMCP,&
             & NB, MC, BTT, BTQ, BTU, BQQ, BQU, BUU)

        IF(iam==0) THEN
           WRITE(*,*) 
           WRITE(*,100) " Matrix C^-1 moved to blocks ", Time()-TiempoInicio
        END IF

        DEALLOCATE(MC)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFMC, NCMC, -1)


        IF(iam==0) THEN
           !write(*,110) 
           WRITE(*,100) " Computing blocks of the product C^-1 Y"
           WRITE(*,*)
           WRITE(*,100) "  Real part ", Time()-TiempoInicio
        END IF

        !Bloques de armonicos - Parte Real

        ALLOCATE(BPTTR(NFBY, NCBA))
        ALLOCATE(BPTER(NFBY, NCBA))
        ALLOCATE(BPTBR(NFBY, NCBA))
        ALLOCATE(BPQTR(NFBX, NCBA))
        ALLOCATE(BPQER(NFBX, NCBA))
        ALLOCATE(BPQBR(NFBX, NCBA))
        ALLOCATE(BPUTR(NFBX, NCBA))
        ALLOCATE(BPUER(NFBX, NCBA))
        ALLOCATE(BPUBR(NFBX, NCBA))

        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, 6)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBY, NCBA, 3)

        CALL ProductosBloquesIMCArmonicos(ICTXT, NCBMCT, NCBMCP, BTT, BTQ, BTU, BQQ, BQU, BUU,&
             & lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP, YTTR, YQER, YQBR, YUER, YUBR,&
             & BPTTR, BPTER, BPTBR, BPQTR, BPQER, BPQBR, BPUTR, BPUER, BPUBR, TiempoInicio)


        DEALLOCATE(YQER)
        DEALLOCATE(YQBR)
        DEALLOCATE(YUER)
        DEALLOCATE(YUBR)
        DEALLOCATE(YTTR)

        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, -4)
        CALL SumaMemoria(ictxt, iam, 0, NFBY, NCBA, -1)

        !Bloques de armonicos - Parte Imaginaria

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "  Imaginary part ", Time()-TiempoInicio
        END IF

        ALLOCATE(BPTTI(NFBY, NCBA))
        ALLOCATE(BPTEI(NFBY, NCBA))
        ALLOCATE(BPTBI(NFBY, NCBA))
        ALLOCATE(BPQTI(NFBX, NCBA))
        ALLOCATE(BPQEI(NFBX, NCBA))
        ALLOCATE(BPQBI(NFBX, NCBA))
        ALLOCATE(BPUTI(NFBX, NCBA))
        ALLOCATE(BPUEI(NFBX, NCBA))
        ALLOCATE(BPUBI(NFBX, NCBA))

        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, 6)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBY, NCBA, 3)

        CALL ProductosBloquesIMCArmonicos(ICTXT, NCBMCT, NCBMCP, BTT, BTQ, BTU, BQQ, BQU, BUU,&
             & lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP, YTTI, YQEI, YQBI, YUEI, YUBI,&
             & BPTTI, BPTEI, BPTBI, BPQTI, BPQEI, BPQBI, BPUTI, BPUEI, BPUBI, TiempoInicio)

        DEALLOCATE(YQEI)
        DEALLOCATE(YQBI)
        DEALLOCATE(YUEI)
        DEALLOCATE(YUBI)

        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, -4)

        DEALLOCATE(YTTI)

        CALL SumaMemoria(ictxt, iam, 0, NFBY, NCBA, -1)

        DEALLOCATE(BTT)
        DEALLOCATE(BTQ)
        DEALLOCATE(BTU)
        DEALLOCATE(BQQ)
        DEALLOCATE(BQU)
        DEALLOCATE(BUU)

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "Blocks of product C^-1 Y already computed", Time()-TiempoInicio
        END IF

        CALL SumaMemoria(ictxt, iam, 0, NFBMCT, NCBMCT, -1)
        CALL SumaMemoria(ictxt, iam, 0, NFBMCT, NCBMCP, -2)
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBMCP, NCBMCP, -3)

        !En este punto la unica matriz que tiene en memoria es C^-1.Matriz armonicos, por bloques 

        IF(iam==0) THEN
           WRITE(*,110) 
           WRITE(*,115) 
           WRITE(*,100) "Step 5. Computing noise bias ", Time()-TiempoInicio
        END IF

        !Calcula Bl

        CALL blacs_barrier(ICTXT, 'ALL')

        CALL CalculaBl(ICTXT, Lugar, NFBY, NFBX, NCBA, Lmax, NB,&
             & BPTTR, BPTER, BPTBR, BPQTR, BPQER, BPQBR, BPUTR, BPUER, BPUBR,&
             & BPTTI, BPTEI, BPTBI, BPQTI, BPQEI, BPQBI, BPUTI, BPUEI, BPUBI, &
             & MapaRuidoT2Bloque, MapaRuidoP2Bloque)

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

        !Solo necesita siete bloques
        !Elimina las partes reales e imaginarias de los que no hacen falta
        DEALLOCATE(BPQTR)
        DEALLOCATE(BPQTI)
        DEALLOCATE(BPUTR)
        DEALLOCATE(BPUTI)

        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, -4)

        ALLOCATE(BPTTC(NFBY, NCBA))
        BPTTC = dcmplx(BPTTR, BPTTI)
        DEALLOCATE(BPTTR)
        DEALLOCATE(BPTTI)

        ALLOCATE(BPTEC(NFBY, NCBA))
        BPTEC = dcmplx(BPTER, BPTEI)
        DEALLOCATE(BPTER)
        DEALLOCATE(BPTEI)

        ALLOCATE(BPTBC(NFBY, NCBA))
        BPTBC = dcmplx(BPTBR, BPTBI)
        DEALLOCATE(BPTBR)
        DEALLOCATE(BPTBI)

        ALLOCATE(BPQEC(NFBX, NCBA))
        BPQEC = dcmplx(BPQER, BPQEI)
        DEALLOCATE(BPQER)
        DEALLOCATE(BPQEI)

        ALLOCATE(BPQBC(NFBX, NCBA))
        BPQBC = dcmplx(BPQBR, BPQBI)
        DEALLOCATE(BPQBR)
        DEALLOCATE(BPQBI)

        ALLOCATE(BPUEC(NFBX, NCBA))
        BPUEC = dcmplx(BPUER, BPUEI)
        DEALLOCATE(BPUER)
        DEALLOCATE(BPUEI)

        ALLOCATE(BPUBC(NFBX, NCBA))
        BPUBC = dcmplx(BPUBR, BPUBI)
        DEALLOCATE(BPUBR)
        DEALLOCATE(BPUBI)

        IF(iam==0) THEN
           !write(*,*)
           !write(*,100) " Blocks of C^-1 Y already in complex form ", Time()-TiempoInicio
           WRITE(*,*)
           WRITE(*,100) " Computing blocks of the harmonic matrix in complex form"
           WRITE(*,100) " and multiplications Y^H (C^-1 Y)"
           WRITE(*,*)
        END IF

        !Calcula los bloques de la matriz de armonicos en forma compleja
        !Va haciendo las multiplicaciones Y^dag.(IMC.Y)

        ALLOCATE(YTTC(NFBY, NCBA))

        CALL SumaMemoria(ictxt, iam, 0, NFBY, NCBA, 2)

        IF(iam==0) THEN
           WRITE(*,100) "  Computing block TT of Y in complex form ", Time()-TiempoInicio
        END IF

        CALL blacs_barrier(ICTXT, "ALL")

        CALL CalculaBloquesMatrizArmonicosYFC(ICTXT, NB, NFBY, NCBA, NPixT, lmax, zt, at, BeamPWDl, YTTC)


        IF(iam==0) THEN
           WRITE(*,100) "  Computing blocks TT, TE and TB of Y^H C^-1 Y ", Time()-TiempoInicio
        END IF

        Alpha = dcmplx(1.0, 0.0)
        Beta  = dcmplx(0.0,0.0)

        NFBA = MAX(1,NUMROC(TotalArmonicos, NB, MYROW, 0, NPROW))

        ALLOCATE(TT(NFBA, NCBA))
        ALLOCATE(TE(NFBA, NCBA))
        ALLOCATE(TB(NFBA, NCBA))

        !Lo bloques son de valores complejos -> por 6
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBA, NCBA, 6)

        !Bloque TT
        CALL DESCINIT(DescBYIMCY, TotalArmonicos, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBA, INFO)
        CALL DESCINIT(DescBIMCY, NpixT, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBY, INFO)
        CALL DESCINIT(DescBY, NpixT, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBY, INFO)

        TT = 0
        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixT, Alpha, YTTC, 1, 1,&
             &DescBY, BPTTC, 1, 1, DescBIMCY, Beta, TT, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   TT ", Time()-TiempoInicio, "s"
        END IF

        !Bloque TE
        TE = 0
        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixT, Alpha, YTTC, 1, 1,&
             &DescBY, BPTEC, 1, 1, DescBIMCY, Beta, TE, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   TE ", Time()-TiempoInicio, "s"
        END IF

        !Bloque TB
        TB = 0
        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixT, Alpha, YTTC, 1, 1,&
             &DescBY, BPTBC, 1, 1, DescBIMCY, Beta, TB, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   TB ", Time()-TiempoInicio, "s"
        END IF

        !Librera espacio de los bloques que ya no necesita
        DEALLOCATE(YTTC)
        DEALLOCATE(BPTTC)
        DEALLOCATE(BPTEC)
        DEALLOCATE(BPTBC)

        !Los cuatro bloques son complejos
        CALL SumaMemoria(ictxt, iam, 0, NFBY, NCBA, -8)

        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) "  Computing blocks QE, QB, UE y UB of Y in complex form ", Time()-TiempoInicio
        END IF

        ALLOCATE(YQEC(NFBX, NCBA))
        ALLOCATE(YQBC(NFBX, NCBA))
        ALLOCATE(YUEC(NFBX, NCBA))
        ALLOCATE(YUBC(NFBX, NCBA))

        !Numero de datos reservados, junto con lo anterior
        CALL SumaMemoria(ictxt, iam, 0, NFBX, NCBA, 8)

        CALL CalculaBloquesMatrizArmonicosXFC(ICTXT, NB, NFBX, NCBA, NPixP, lmax, zp, ap, BeamPWDl, &
             & YQEC, YQBC, YUEC, YUBC)  

        CALL blacs_barrier(ICTXT, "ALL")

        IF(iam==0) THEN
           WRITE(*,100) "  Computing blocks EE, BB and EB of Y^H C^-1 Y  ", Time()-TiempoInicio
        END IF

        !Reserva memoria para los productos EE, BB y EB

        ALLOCATE(EE(NFBA, NCBA))
        ALLOCATE(BB(NFBA, NCBA))
        ALLOCATE(EB(NFBA, NCBA))

        !Los bloques son complejos
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBA, NCBA, 6)

        !Multiplicaciones

        !Bloque EE
        EE = 0
        CALL DESCINIT(DescBYIMCY, TotalArmonicos, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBA, INFO)
        CALL DESCINIT(DescBIMCY, NpixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBX, INFO)
        CALL DESCINIT(DescBY, NpixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBX, INFO)

        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YQEC, 1, 1,&
             &DescBY, BPQEC, 1, 1, DescBIMCY, Beta, EE, 1, 1, DescBYIMCY)

        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YUEC, 1, 1,&
             &DescBY, BPUEC, 1, 1, DescBIMCY, Alpha, EE, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   EE ", Time()-TiempoInicio, "s"
        END IF

        !Bloque BB
        BB = 0
        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YQBC, 1, 1,&
             &DescBY, BPQBC, 1, 1, DescBIMCY, Beta, BB, 1, 1, DescBYIMCY)

        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YUBC, 1, 1,&
             &DescBY, BPUBC, 1, 1, DescBIMCY, Alpha, BB, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   BB ", Time()-TiempoInicio, "s"
        END IF

        !Bloque EB
        EB = 0
        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YQEC, 1, 1,&
             &DescBY, BPQBC, 1, 1, DescBIMCY, Beta, EB, 1, 1, DescBYIMCY)

        CALL blacs_barrier(ICTXT, "ALL")
        CALL pzgemm("C", "N", TotalArmonicos, TotalArmonicos, NPixP, Alpha, YUEC, 1, 1,&
             &DescBY, BPUBC, 1, 1, DescBIMCY, Alpha, EB, 1, 1, DescBYIMCY)

        IF(iam==0) THEN
           WRITE(*,*) "   EB ", Time()-TiempoInicio, "s"
        END IF

        !Librera espacio de los bloques que ya no necesita
        DEALLOCATE(YQEC)
        DEALLOCATE(YQBC)
        DEALLOCATE(YUEC)
        DEALLOCATE(YUBC)
        DEALLOCATE(BPQEC)
        DEALLOCATE(BPQBC)
        DEALLOCATE(BPUEC)
        DEALLOCATE(BPUBC)

        !Los bloues son complejos
        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBX, NCBA, -16)


        IF(iam==0) THEN
           WRITE(*,*)
           WRITE(*,100) " Moving blocks to real and imaginary part ", Time()-TiempoInicio
        END IF

        !Pasa a bloques de Y.IMC.Y en forma real e imaginaria

        ALLOCATE(TTR(NFBA, NCBA))
        ALLOCATE(TTI(NFBA, NCBA))
        TTR = DREAL(TT)
        TTI = DIMAG(TT)
        DEALLOCATE(TT)

        ALLOCATE(TER(NFBA, NCBA))
        ALLOCATE(TEI(NFBA, NCBA))
        TER = DREAL(TE)
        TEI = DIMAG(TE)
        DEALLOCATE(TE)

        ALLOCATE(TBR(NFBA, NCBA))
        ALLOCATE(TBI(NFBA, NCBA))
        TBR = DREAL(TB)
        TBI = DIMAG(TB)
        DEALLOCATE(TB)

        ALLOCATE(EER(NFBA, NCBA))
        ALLOCATE(EEI(NFBA, NCBA))
        EER = DREAL(EE)
        EEI = DIMAG(EE)
        DEALLOCATE(EE)

        ALLOCATE(BBR(NFBA, NCBA))
        ALLOCATE(BBI(NFBA, NCBA))
        BBR = DREAL(BB)
        BBI = DIMAG(BB)
        DEALLOCATE(BB)

        ALLOCATE(EBR(NFBA, NCBA))
        ALLOCATE(EBI(NFBA, NCBA))
        EBR = DREAL(EB)
        EBI = DIMAG(EB)
        DEALLOCATE(EB)

        IF(iam==0) THEN
           WRITE(*,*) 
           WRITE(*,100) " Building transposed blocks ", Time()-TiempoInicio
        END IF

        CALL CalculaMatrizFisher(ICTXT, NFBA, NCBA, lmax, NB, TTR, TTI, TER, TEI, TBR,&
             & TBI, EER, EEI, BBR, BBI, EBR, EBI, Lugar, TiempoInicio)


        DEALLOCATE(TTR)
        DEALLOCATE(TTI)
        DEALLOCATE(TER)
        DEALLOCATE(TEI)
        DEALLOCATE(TBR)
        DEALLOCATE(TBI)
        DEALLOCATE(EER)
        DEALLOCATE(EEI)
        DEALLOCATE(BBR)
        DEALLOCATE(BBI)
        DEALLOCATE(EBR)
        DEALLOCATE(EBI)

        IF(iam==0) THEN
           WRITE(*,100) "Fisher matrix already computed ", Time()-TiempoInicio
        END IF

        CALL SumaMemoria(ictxt, iam, MuestraMemoria, NFBA, NCBA, -12)

     ELSE

        IF(iam==0) THEN
           WRITE(*,110)
           WRITE(*,*) "The program does not compute the Fisher matrix"
           WRITE(*,*) "Elapsed time: ", Time() - TiempoInicio, "s"
           WRITE(*,110) 
        END IF

     END IF

     CALL blacs_barrier(ICTXT,'All')

     IF(iam==0) THEN

        IF(ComputeSpectrum==1) THEN

           IF(Binned==0) THEN

              WRITE(*,110)
              WRITE(*,*) "Computing power spectrum"

              CALL CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)
           ELSE   
              WRITE(*,110)
              WRITE(*,*) "Computing binned power espectrum"

              CALL CalculaBineado(NSide, npixT, npixP, lmax, RuidoTT, RuidoQQ, DlTT, DlEE, DlBB, DlTE,&
                   & DlTB, DlEB, BeamPWDl, KindOfGrouping, KindOfBinCenter, QuitarSesgo, NMapas, Lugar, FileBinLimits)

           END IF !Binear o no binear

        ELSE

           WRITE(*,*) "The program does not compute the power spectrum"

        END IF !Calcular espectro

     END IF !Soy cero y calculo espectro

     IF(iam==0) THEN

        TiempoFinal = Time()
        WRITE(*,*)
        WRITE(*,*) "Elapsed time: ", TiempoFinal - TiempoInicio, "s"
        WRITE(*,*) "End: ", CTime(Time())
        WRITE(*,110) 
     END IF

  END IF !Matriz de covarianza regular o singular

  CALL BLACS_GRIDEXIT(ICTXT)
  CALL BLACS_EXIT(0)

END PROGRAM ECLIPSE_TEB


!******************************************************************************************************

SUBROUTINE CalculaMatrizFisher(ICTXT, NFBA, NCBA, lmax, NB, TTR, TTI, TER, TEI, TBR,&
     & TBI, EER, EEI, BBR, BBI, EBR, EBI, lugar, TiempoInicio)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE CalculaBloqueMatrizFisher(ICTXT, Bloque, NFilas, NColumnas, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
       INTEGER, INTENT(in) :: ICTXT, NFilas, NColumnas, Lmax
       REAL(kind=8), DIMENSION(NFilas, NColumnas), INTENT(in) :: Bloque
       INTEGER, DIMENSION(Nfilas), INTENT(in) :: IndicesFila
       INTEGER, DIMENSION(NColumnas), INTENT(in) :: IndicesColumna
       REAL(kind=8), DIMENSION(2:lmax, 2:lmax), INTENT(out) :: BloqueFisher
     END SUBROUTINE CalculaBloqueMatrizFisher

     SUBROUTINE MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
       INTEGER, INTENT(in) :: Lmax, Caso1, Caso2
       REAL(kind=8), INTENT(in) :: Factor
       REAL(kind=8), DIMENSION(2:Lmax, 2:Lmax), INTENT(in) :: BloqueFisher
       REAL(kind=8), DIMENSION(6, 2:Lmax, 6, 2:Lmax), INTENT(inout) :: MF
     END SUBROUTINE MatrizFisherGlobal
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, NFBA, NCBA, lmax, NB
  INTEGER(kind=8), INTENT(in) ::  TiempoInicio
  REAL(kind=8), DIMENSION(NFBA, NCBA), INTENT(in) ::  TTR, TTI, TER, TEI, TBR, TBI, EER, EEI, BBR, BBI, EBR, EBI
  CHARACTER(len=100), INTENT(in) :: lugar

  !*********************************

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:,:) :: MF
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: Salida
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: TERT, TEIT, TBRT, TBIT, EBRT, EBIT
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: BloqueFisher, Producto
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesFila, IndicesColumna
  INTEGER :: TotalArmonicos, i,j, Caso1, caso2, l1
  REAL(kind=8) :: Factor
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO
  INTEGER, DIMENSION(9) :: Desc

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, descinit, pdtran
  INTEGER, EXTERNAL :: INDXL2G

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2*lmax + lmax**2

  CALL DESCINIT(Desc, TotalArmonicos, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBA, INFO)

  !Calcula las trazas
  ALLOCATE(MF(6, 2:lmax, 6, 2:lmax))
  MF = 0D0

  ALLOCATE(Producto(NFBA, NCBA))
  ALLOCATE(BloqueFisher(2:lmax,2:lmax))
  ALLOCATE(IndicesFila(NFBA))
  ALLOCATE(IndicesColumna(NCBA))
  Producto = 0
  BloqueFisher = 0

  !Localica las posiciones globales de los elementos de la submatriz de cada procesador
  DO i=1,NFBA
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     IndicesFila(i) = FLOOR(SQRT(3.0 + j))
  END DO

  DO i=1,NCBA
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     IndicesColumna(i) = FLOOR(SQRT(3.0 + j))
  END DO

  !Necesita bloque traspuesto

  ALLOCATE(TERT(NFBA, NCBA))
  ALLOCATE(TEIT(NFBA, NCBA))
  TERT = 0
  TEIT = 0

  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, TER, 1, 1, Desc, 0.0D0, TERT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, TEI, 1, 1, Desc, 0.0D0, TEIT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")

  ALLOCATE(TBRT(NFBA, NCBA))
  ALLOCATE(TBIT(NFBA, NCBA))
  TBRT = 0
  TBIT = 0

  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, TBR, 1, 1, Desc, 0.0D0, TBRT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, TBI, 1, 1, Desc, 0.0D0, TBIT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")

  !Necesita traspuestos
  ALLOCATE(EBRT(NFBA, NCBA))
  ALLOCATE(EBIT(NFBA, NCBA))
  EBRT = 0
  EBIT = 0

  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, EBR, 1, 1, Desc, 0.0D0, EBRT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")
  CALL pdtran(TotalArmonicos, TotalArmonicos, 1.0D0, EBI, 1, 1, Desc, 0.0D0, EBIT, 1, 1, Desc)
  CALL blacs_barrier(ICTXT, "All")

  IF(iam==0) THEN
     WRITE(*,*)
     WRITE(*,*) "  Computing blocks of the Fisher matrix"
     WRITE(*,*)
  END IF

  !TTTT
  Caso1 = 1
  Caso2 = 1
  Factor = 1.0D0
  Producto = TTR * TTR + TTI * TTI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TTTT "
  END IF


  !EEEE
  Caso1 = 2
  Caso2 = 2
  Factor = 1.0D0
  Producto =  EER * EER + EEI * EEI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EEEE"
  END IF

  !BBBB
  Caso1 = 3
  Caso2 = 3
  Factor = 1.0D0
  Producto =  BBR * BBR + BBI * BBI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   BBBB"
  END IF

  !TTEE
  Caso1 = 1
  Caso2 = 2
  Factor = 1.0D0
  Producto =  TER * TER + TEI * TEI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TTEE"
  END IF

  !TTBB
  Caso1 = 1
  Caso2 = 3
  Factor = 1.0D0
  Producto =  TBR * TBR + TBI * TBI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TTBB"
  END IF

  !EEBB
  Caso1 = 2
  Caso2 = 3
  Factor = 1.0D0
  Producto =  EBR * EBR + EBI * EBI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EEBB"
  END IF

  !TTTE
  Caso1 = 1
  Caso2 = 4
  Factor = 2.0D0
  Producto =  TTR * TER + TTI * TEI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN

     !systime = CTIME(TIME())
     WRITE (*,*) "   TTTE"
  END IF

  !TTTB
  Caso1 = 1
  Caso2 = 5
  Factor = 2.0D0
  Producto =  TTR * TBR + TTI * TBI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TTTB"
  END IF

  !TTEB
  Caso1 = 1
  Caso2 = 6
  Factor = 2.0D0
  Producto =  TER * TBR + TEI * TBI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TTEB"
  END IF

  !EETE
  Caso1 = 2
  Caso2 = 4
  Factor = 2.0D0
  Producto =  EER * TERT - EEI * TEIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EETE"
  END IF

  !EETB
  Caso1 = 2
  Caso2 = 5
  Factor = 2.0D0
  Producto =  EBR * TERT - EBI * TEIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EETB"
  END IF

  !EEEB
  Caso1 = 2
  Caso2 = 6
  Factor = 2.0D0
  Producto =  EBR * EER + EBI * EEI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EEEB"
  END IF


  !BBTE
  Caso1 = 3
  Caso2 = 4
  Factor = 2.0D0
  Producto =  EBRT * TBRT + EBIT * TBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   BBTE"
  END IF

  !BBTB
  Caso1 = 3
  Caso2 = 5
  Factor = 2.0D0
  Producto =  BBR * TBRT - BBI * TBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   BBTB"
  END IF


  !BBEB
  Caso1 = 3
  Caso2 = 6
  Factor = 2.0D0
  Producto =  BBR * EBRT - BBI * EBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   BBEB"
  END IF

  !TETE
  Caso1 = 4
  Caso2 = 4
  Factor = 2.0D0
  Producto =  EER * TTR + EEI * TTI + TER * TERT - TEI * TEIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TETE"
  END IF

  !TETB
  Caso1 = 4
  Caso2 = 5
  Factor = 2.0D0
  Producto =  EBR * TTR + EBI * TTI + TBR * TERT - TBI * TEIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TETB"
  END IF

  !TEEB
  Caso1 = 4
  Caso2 = 6
  Factor = 2.0D0
  Producto =  EER * TBR + EEI * TBI + EBR * TER + EBI * TEI
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TEEB"
  END IF


  !TBTB
  Caso1 = 5
  Caso2 = 5
  Factor = 2.0D0
  Producto =  TTR * BBR + TTI * BBI + TBR * TBRT - TBI * TBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TBTB"
  END IF


  !TBEB
  Caso1 = 5
  Caso2 = 6
  Factor = 2.0D0
  Producto =  TER * BBR + TEI * BBI + TBR * EBRT - TBI * EBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   TBEB"
  END IF

  !EBEB
  Caso1 = 6
  Caso2 = 6
  Factor = 2.0D0
  Producto =  BBR * EER + BBI * EEI + EBR * EBRT - EBI * EBIT
  CALL CalculaBloqueMatrizFisher(ICTXT, Producto, NFBA, NCBA, IndicesFila, IndicesColumna, BloqueFisher, Lmax)
  IF (iam==0) CALL MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)
  IF (IAM==0) THEN
     !systime = CTIME(TIME())
     WRITE (*,*) "   EBEB"
  END IF



  !**************************************************************************************

  IF(IAM==0) THEN

     WRITE(*,*)
     WRITE(*,*) "   Blocks of the Fisher matrix already computed ", Time()-TiempoInicio, "s"  
     WRITE(*,*)
     WRITE(*,*) "   Saving the Fisher matrix"
     WRITE(*,*)

     ALLOCATE(Salida(2:lmax))

     OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrix.dat", Action="write")

     DO Caso1 = 1,6
        DO l1 = 2, lmax
           DO Caso2 = 1,6
              Salida(2:lmax) = MF(Caso1,l1,Caso2,2:lmax)
              WRITE (25, *) Salida(2:lmax)
           END DO
        END DO
     END DO

     DEALLOCATE(Salida)
     CLOSE(25)
  END IF

  CALL blacs_barrier( ICTXT, "A" )


END SUBROUTINE CalculaMatrizFisher

SUBROUTINE CalculaBloqueMatrizFisher(ICTXT, Bloque, NFilas, NColumnas, IndicesFila, IndicesColumna, BloqueFisher, Lmax)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFilas, NColumnas, Lmax
  REAL(kind=8), DIMENSION(NFilas, NColumnas), INTENT(in) :: Bloque
  INTEGER, DIMENSION(Nfilas), INTENT(in) :: IndicesFila
  INTEGER, DIMENSION(NColumnas), INTENT(in) :: IndicesColumna
  REAL(kind=8), DIMENSION(2:lmax, 2:lmax), INTENT(out) :: BloqueFisher

  !***********

  INTEGER :: i,j,l1,l2

  EXTERNAL :: blacs_barrier, dgsum2d

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



SUBROUTINE MatrizFisherGlobal(MF, BloqueFisher, Lmax, Caso1, Caso2, Factor)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: Lmax, Caso1, Caso2
  REAL(kind=8), INTENT(in) :: Factor
  REAL(kind=8), DIMENSION(2:Lmax, 2:Lmax), INTENT(in) :: BloqueFisher
  REAL(kind=8), DIMENSION(6, 2:Lmax, 6, 2:Lmax), INTENT(inout) :: MF

  INTEGER :: l1, l2

  DO l1=2,Lmax
     DO l2=2,Lmax
        MF(Caso1, l1, Caso2, l2) = Factor * BloqueFisher(l1,l2)
        MF(Caso2, l2, Caso1, l1) = MF(Caso1, l1, Caso2, l2)
     END DO
  END DO

END SUBROUTINE MatrizFisherGlobal



SUBROUTINE BloquesMatrizCovarianza(ICTXT, DimMC, NFMC, NCMC, NpixT, NpixP, NFBMCT, NCBMCT, NFBMCP, NCBMCP,&
     & NB, MC, BTT, BTQ, BTU, BQQ, BQU, BUU)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, DimMC, NFMC, NCMC, NpixT, NpixP, NFBMCT, NCBMCT, NFBMCP, NCBMCP, NB
  REAL(kind=8), DIMENSION(NFMC, NCMC), INTENT(in) :: MC
  REAL(kind=8), DIMENSION(NFBMCT, NCBMCT), INTENT(out) :: BTT
  REAL(kind=8), DIMENSION(NFBMCT, NCBMCP), INTENT(out) :: BTQ, BTU
  REAL(kind=8), DIMENSION(NFBMCP, NCBMCP), INTENT(out) :: BQQ, BQU, BUU

  !**************************************************************
  INTEGER, DIMENSION(9) :: DescMC, DescBloque

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, Fila, Columna, INFO

  EXTERNAL :: BLACS_GRIDINFO, DESCINIT, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  CALL DESCINIT(DescMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)


  IF(iam==0) WRITE(*,*) "   TT"
  CALL DESCINIT(DescBloque, NpixT, NpixT, NB, NB, 0, 0, ICTXT, NFBMCT, INFO)
  Fila = 1
  Columna = 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixT, NpixT, MC, Fila, Columna, DescMC, BTT, 1, 1, DescBloque, ictxt)

  IF(iam==0) WRITE(*,*) "   TQ"
  CALL DESCINIT(DescBloque, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBMCT, INFO)
  Fila = 1
  Columna = NPixT + 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixT, NpixP, MC, Fila, Columna, DescMC, BTQ, 1, 1, DescBloque, ictxt)

  IF(iam==0) WRITE(*,*) "   TU"
  CALL DESCINIT(DescBloque, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBMCT, INFO)
  Fila = 1
  Columna = NPixT + NPixP + 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixT, NpixP, MC, Fila, Columna, DescMC, BTU, 1, 1, DescBloque, ictxt)

  IF(iam==0) WRITE(*,*) "   QQ"
  CALL DESCINIT(DescBloque, NpixP, NpixP, NB, NB, 0, 0, ICTXT, NFBMCP, INFO)
  Fila = NPixT + 1
  Columna = NPixT + 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixP, NpixP, MC, Fila, Columna, DescMC, BQQ, 1, 1, DescBloque, ictxt)

  IF(iam==0) WRITE(*,*) "   QU"
  Fila = NPixT + 1
  Columna = NPixT + NPixP + 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixP, NpixP, MC, Fila, Columna, DescMC, BQU, 1, 1, DescBloque, ictxt)

  IF(iam==0) WRITE(*,*) "   UU"
  Fila = NPixT + NPixP + 1
  Columna = NPixT + NPixP + 1
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemr2d(NPixP, NpixP, MC, Fila, Columna, DescMC, BUU, 1, 1, DescBloque, ictxt)

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE BloquesMatrizCovarianza


SUBROUTINE  CalculaBl(ICTXT, Lugar, NFBY, NFBX, NCBA, Lmax, NB, BTTR, BTER, BTBR, BQTR, BQER, BQBR, BUTR, BUER, BUBR,&
     & BTTI, BTEI, BTBI, BQTI, BQEI, BQBI, BUTI, BUEI, BUBI, RuidoTT, RuidoQQ)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFBY, NFBX, NCBA, Lmax, NB
  CHARACTER(len=100), INTENT(in) :: Lugar
  REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(IN) :: BTTR, BTER, BTBR, BTTI, BTEI, BTBI
  REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(IN) :: BQTR, BQER, BQBR, BQTI, BQEI, BQBI
  REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(IN) :: BUTR, BUER, BUBR, BUTI, BUEI, BUBI
  REAL(kind=8), DIMENSION(NFBY), INTENT(IN) :: RuidoTT
  REAL(kind=8), DIMENSION(NFBX), INTENT(IN) :: RuidoQQ 

  !****************************

  INTEGER :: TotalArmonicos
  INTEGER :: i, j, l, Caso
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER, EXTERNAL :: NUMROC

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: TempT, TempP, Bl
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: SumaT, SumaQ, SumaU, SumaTT, SumaEE, SumaBB, SumaTE, SumaTB, SumaEB

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, DGSUM2D
  INTEGER, EXTERNAL :: INDXL2G 

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2*lmax + lmax**2

  !***********************************************

  ALLOCATE(SumaTT(NCBA))
  ALLOCATE(SumaEE(NCBA))
  ALLOCATE(SumaBB(NCBA))
  ALLOCATE(SumaTE(NCBA))
  ALLOCATE(SumaTB(NCBA))
  ALLOCATE(SumaEB(NCBA))

  !TT

  !Parte T
  ALLOCATE(SumaT(NCBA))
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTTR**2+BTTI**2
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(SumaQ(NCBA))
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQTR**2+BQTI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  ALLOCATE(SumaU(NCBA))
  TempP = BUTR**2+BUTI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaTT = SumaT + SumaQ + SumaU

  !***********************************************
  !EE

  !Parte T
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTER**2+BTEI**2
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQER**2+BQEI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  TempP = BUER**2+BUEI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaEE = SumaT + SumaQ + SumaU

  !***********************************************
  !BB

  !Parte T
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTBR**2+BTBI**2
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQBR**2+BQBI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  TempP = BUBR**2+BUBI**2
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaBB = SumaT + SumaQ + SumaU

  !***********************************************
  !TE

  !Parte T
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTTR*BTER+BTTI*BTEI
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQTR*BQER+BQTI*BQEI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  TempP = BUTR*BUER+BUTI*BUEI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaTE = SumaT + SumaQ + SumaU
  SumaTE = 2 * SumaTE

  !***********************************************
  !TB

  !Parte T
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTTR*BTBR+BTTI*BTBI
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQTR*BQBR+BQTI*BQBI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  TempP = BUTR*BUBR+BUTI*BUBI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaTB = SumaT + SumaQ + SumaU
  SumaTB = 2*SumaTB

  !***********************************************
  !EB

  !Parte T
  ALLOCATE(TempT(NFBY, NCBA))
  TempT = BTER*BTBR+BTEI*BTBI
  DO i=1,NCBA
     TempT(:,i) = RuidoTT(:) * TempT(:,i)
  END DO
  SumaT = SUM(TempT, Dim=1)
  DEALLOCATE(TempT)

  !Parte Q
  ALLOCATE(TempP(NFBX, NCBA))
  TempP = BQER*BQBR+BQEI*BQBI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaQ = SUM(TempP, Dim=1)

  !Parte U
  TempP = BUER*BUBR+BUEI*BUBI
  DO i=1,NCBA
     TempP(:,i) = RuidoQQ(:) * TempP(:,i)
  END DO
  SumaU = SUM(TempP, Dim=1)
  DEALLOCATE(TempP)

  SumaEB = SumaT + SumaQ + SumaU
  SumaEB = 2*SumaEB

  DEALLOCATE(SumaT)
  DEALLOCATE(SumaQ)
  DEALLOCATE(SumaU)

  !Recoge

  ALLOCATE(Bl(6,2:lmax))
  Bl = 0

  !Coloca las sumas repartidas por los procesadores en
  !los Bl que les corresponden
  DO i=1,NCBA
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) THEN
        Bl(1,l) = Bl(1,l) + SumaTT(i)
        Bl(2,l) = Bl(2,l) + SumaEE(i)
        Bl(3,l) = Bl(3,l) + SumaBB(i)
        Bl(4,l) = Bl(4,l) + SumaTE(i)
        Bl(5,l) = Bl(5,l) + SumaTB(i)
        Bl(6,l) = Bl(6,l) + SumaEB(i)
     END IF
  END DO

  DEALLOCATE(SumaTT)
  DEALLOCATE(SumaEE)
  DEALLOCATE(SumaBB)
  DEALLOCATE(SumaTE)
  DEALLOCATE(SumaTB)
  DEALLOCATE(SumaEB)

  !Suma todos los Bl repartidos por los procesadores
  CALL blacs_barrier(ICTXT, "ALL")
  CALL DGSUM2D( ictxt, 'All', '1-tree', 6, lMax-1, Bl, 6, -1, -1)

  Bl = Bl/2.0D0

  IF (iam==0) THEN        
     OPEN(Unit=20, File=TRIM(Lugar)//"/NoiseBias.dat", Action="write")
     DO Caso=1,6
        DO l=2,lMax
           WRITE (20,*) Bl(Caso,l)
        END DO
     END DO
     CLOSE(20)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")


END SUBROUTINE CalculaBl


SUBROUTINE  ProductosBloquesIMCArmonicos(ICTXT, NCBMCT, NCBMCP, BTT, BTQ, BTU, BQQ, BQU, BUU,&
     & lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP, YTT, YQE, YQB, YUE, YUB,&
     & BPTT, BPTE, BPTB, BPQT, BPQE, BPQB, BPUT, BPUE, BPUB, TiempoInicio)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NCBMCT, NCBMCP, lmax, NFBY, NFBX, NCBA, NB, NpixT, NPixP
  REAL(kind=8), DIMENSION(NFBY, NCBMCT), INTENT(in) :: BTT
  REAL(kind=8), DIMENSION(NFBY, NCBMCP), INTENT(in ) :: BTQ, BTU
  REAL(kind=8), DIMENSION(NFBX, NCBMCP), INTENT(in) :: BQQ, BQU, BUU

  REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(in) :: YTT
  REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(in) :: YQE, YQB, YUE, YUB

  INTEGER(kind=8), INTENT(in) ::  TiempoInicio

  REAL(kind=8), DIMENSION(NFBY, NCBA), INTENT(out) :: BPTT, BPTE, BPTB
  REAL(kind=8), DIMENSION(NFBX, NCBA), INTENT(out) :: BPQT, BPQE, BPQB, BPUT, BPUE, BPUB

  !*****************************************

  INTEGER :: TotalArmonicos
  INTEGER, DIMENSION(9) :: DescBIMC, DescFactor, DescProducto
  INTEGER :: INFO

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER, EXTERNAL :: NUMROC

  EXTERNAL :: BLACS_GRIDINFO, blacs_barrier, descinit, pdsymm, pdgemm, pdgemr2d

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  TotalArmonicos = -3 + 2*lmax + lmax**2

  !Bloques de filas -> pixel, columnas -> armonicos

  !*****************************************************
  !Fila 1

  !Bloque TT
  !Calcula TT.YTT
  BPTT = 0
  CALL descinit(DescBIMC, NpixT, NpixT, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescFactor, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescProducto, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdsymm("L", "U", NPixT, TotalArmonicos, 1.0D0, BTT, 1, 1, DescBIMC,&
       &YTT, 1, 1, DescFactor, 0.0D0, BPTT, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    TT ", Time()-TiempoInicio, "s"

  !Bloque TE - Parte 1
  BPTE = 0
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixT, TotalArmonicos, NPixP, 1.0D0, BTQ, 1, 1, DescBIMC,&
       & YQE, 1, 1, DescFactor, 0.0D0, BPTE, 1, 1, DescProducto)

  !Bloque TE - Parte 2
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixT, TotalArmonicos, NPixP, 1.0D0, BTU, 1, 1, DescBIMC,&
       & YUE, 1, 1, DescFactor, 1.0D0, BPTE, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    TE ", Time()-TiempoInicio, "s"

  !Bloque TB - Parte 1
  BPTB = 0
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixT, TotalArmonicos, NPixP, 1.0D0, BTQ, 1, 1, DescBIMC,&
       & YQB, 1, 1, DescFactor, 0.0D0, BPTB, 1, 1, DescProducto)

  !Bloque TE - Parte 2
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixT, TotalArmonicos, NPixP, 1.0D0, BTU, 1, 1, DescBIMC,&
       & YUB, 1, 1, DescFactor, 1.0D0, BPTB, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    TB ", Time()-TiempoInicio, "s"

  !*****************************************************
  !Fila 2

  !Bloque QT
  !Calcula: QT.YTT
  BPQT = 0
  !El bloque natrual de IMC es QT. No lo tengo -> Multpiplica por trastuesto de TQ
  !Descriptor de TQ!!!!
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)

  CALL descinit(DescFactor, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("T", "N", NPixP, TotalArmonicos, NPixT, 1.0D0, BTQ, 1, 1, DescBIMC,&
       & YTT, 1, 1, DescFactor, 0.0D0, BPQT, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    QT ", Time()-TiempoInicio, "s"


  !Bloque QE
  !Calcula: QQ.YQE + QU.YUE
  !Como multipilica por QQ -> multiplicacion matriz simetrica
  BPQE = 0
  CALL descinit(DescBIMC, NpixP, NpixP, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdsymm("L", "U", NPixP, TotalArmonicos, 1.0D0, BQQ, 1, 1, DescBIMC,&
       &YQE, 1, 1, DescFactor, 0.0D0, BPQE, 1, 1, DescProducto)

  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixP, TotalArmonicos, NPixP, 1.0D0, BQU, 1, 1, DescBIMC,&
       & YUE, 1, 1, DescFactor, 1.0D0, BPQE, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    QE ", Time()-TiempoInicio, "s"

  !Bloque QB
  !Calcula: QQ.YQB + QU.YUB
  !Como multpilica por QQ -> multiplicacion matriz simetrica
  BPQB = 0
  CALL descinit(DescBIMC, NpixP, NpixP, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdsymm("L", "U", NPixP, TotalArmonicos, 1.0D0, BQQ, 1, 1, DescBIMC,&
       &YQB, 1, 1, DescFactor, 0.0D0, BPQB, 1, 1, DescProducto)

  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("N", "N", NPixP, TotalArmonicos, NPixP, 1.0D0, BQU, 1, 1, DescBIMC,&
       & YUB, 1, 1, DescFactor, 1.0D0, BPQB, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    QB ", Time()-TiempoInicio, "s"


  !*****************************************************
  !Fila 3

  !Bloque UT
  !Calcula: UT.YTT
  BPUT = 0
  !El bloque natrual de IMC es UT. No lo tengo -> Multpiplica por trastuesto de TU
  !Descriptor de TU!!!!
  CALL descinit(DescBIMC, NpixT, NpixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)

  CALL descinit(DescFactor, NPixT, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBY, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("T", "N", NPixP, TotalArmonicos, NPixT, 1.0D0, BTU, 1, 1, DescBIMC,&
       & YTT, 1, 1, DescFactor, 0.0D0, BPUT, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    UT ", Time()-TiempoInicio, "s"


  !Bloque UE
  !Calcula: UQ.YQE + UU.YUE
  !Como multipilica por UU -> multiplicacion matriz simetrica
  !Primer producto: UU.YUE
  BPUE = 0
  CALL descinit(DescBIMC, NpixP, NpixP, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdsymm("L", "U", NPixP, TotalArmonicos, 1.0D0, BUU, 1, 1, DescBIMC,&
       &YUE, 1, 1, DescFactor, 0.0D0, BPUE, 1, 1, DescProducto)

  !Segundo producto: UQ.YQE
  !Usa la traspuesta de QU
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("T", "N", NPixP, TotalArmonicos, NPixP, 1.0D0, BQU, 1, 1, DescBIMC,&
       & YQE, 1, 1, DescFactor, 1.0D0, BPUE, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    UE ", Time()-TiempoInicio, "s"

  !Bloque UB
  !Calcula: QQ.YQB + QU.YUB
  !Como multpilica por UU -> multiplicacion matriz simetrica
  !Primer producto: UU.YUB
  BPUB = 0
  CALL descinit(DescBIMC, NpixP, NpixP, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescFactor, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL descinit(DescProducto, NPixP, TotalArmonicos, NB ,NB, 0, 0, ICTXT, NFBX, INFO)
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdsymm("L", "U", NPixP, TotalArmonicos, 1.0D0, BUU, 1, 1, DescBIMC,&
       &YUB, 1, 1, DescFactor, 0.0D0, BPUB, 1, 1, DescProducto)

  !Segundo producto: UQ.YQB
  !Usa la traspuesta de QU
  CALL blacs_barrier(ICTXT, "ALL")
  CALL pdgemm("T", "N", NPixP, TotalArmonicos, NPixP, 1.0D0, BQU, 1, 1, DescBIMC,&
       & YQB, 1, 1, DescFactor, 1.0D0, BPUB, 1, 1, DescProducto)

  IF(iam==0) WRITE(*,*) "    UB ", Time()-TiempoInicio, "s"

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE ProductosBloquesIMCArmonicos


SUBROUTINE CalculaYDagIMCMapas(ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, lmax, NMapas, NCMapas, NB,&
     & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, IMCMapas, PTR, PER, PBR, PTI, PEI, PBI)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFBY, NFBX, NFMC, NFMapasA, NCBA, NpixT, NpixP, Lmax, NMapas, NCMapas, NB
  REAL(KIND=8), DIMENSION(NFBY,NCBA), INTENT(in)  :: YTTR, YTTI
  REAL(KIND=8), DIMENSION(NFBX,NCBA), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
  REAL(KIND=8), DIMENSION(NFMC,NCMapas), INTENT(in)  :: IMCMapas

  REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(out) :: PTR, PER, PBR, PTI, PEI, PBI 

  !************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: TotalArmonicos, DimMC, Fila, INFO
  INTEGER(kind=4), DIMENSION(9) :: DescBloqueArmonicos, DescIMCMapas, DescProducto

  REAL (kind=8), ALLOCATABLE, DIMENSION(:,:) ::  P1

  EXTERNAL :: BLACS_GRIDINFO, DESCINIT, blacs_barrier, pdgemm

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  PTR = 0
  PTI = 0
  PER = 0
  PEI = 0
  PBR = 0
  PBI = 0

  TotalArmonicos = -3 + 2*lmax + lmax**2
  DimMC = NpixT + 2*NPixP


  !*******************************************************************
  !Parte T
  CALL DESCINIT(DescBloqueArmonicos, NpixT, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL DESCINIT(DescIMCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL DESCINIT(DescProducto, TotalArmonicos, NMapas, NB, NB, 0, 0, ICTXT, NFMapasA, INFO)

  !Producto YTTR^t.IMCMapasT
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixT, 1.0D0, YTTR, 1, 1, DescBloqueArmonicos, IMCMapas, 1, 1,&
       & DescIMCMapas, 0.0D0, PTR, 1, 1, DescProducto)

  !Producto YTTI^t.IMCMapasT
  CALL blacs_barrier(ICTXT, "All")
  !Intriduce YTTI conjugada (-1 delante)
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixT, -1.0D0, YTTI, 1, 1, DescBloqueArmonicos, IMCMapas, 1, 1,&
       & DescIMCMapas, 0.0D0, PTI, 1, 1, DescProducto)

  !*******************************************************************
  !Parte E

  !Real
  !En dos etapas
  ALLOCATE(P1(NFMapasA, NCMAPAS))
  P1 = 0

  !Cambia el numero de filas del bloque
  CALL DESCINIT(DescBloqueArmonicos, NpixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  !Multiplica sobre la parte Q del IMC.Mapas
  Fila = NpixT+1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, 1.0D0, YQER, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, P1, 1, 1, DescProducto)

  !Multiplica sobre la parte U del IMC.Mapas
  Fila = NpixT + NpixP + 1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, 1.0D0, YUER, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, PER, 1, 1, DescProducto)

  PER = PER + P1

  !Imaginaria
  Fila = NpixT+1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, -1.0D0, YQEI, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, P1, 1, 1, DescProducto)

  Fila = NpixT + NpixP + 1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, -1.0D0, YUEI, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, PEI, 1, 1, DescProducto)

  PEI = PEI + P1

  !*******************************************************************
  !Parte B

  Fila = NpixT+1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, 1.0D0, YQBR, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, P1, 1, 1, DescProducto)

  Fila = NpixT + NpixP + 1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, 1.0D0, YUBR, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, PBR, 1, 1, DescProducto)

  PBR = PBR + P1


  Fila = NpixT+1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, -1.0D0, YQBI, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, P1, 1, 1, DescProducto)

  Fila = NpixT + NpixP + 1
  CALL blacs_barrier(ICTXT, "All")
  CALL pdgemm("T", "N", TotalArmonicos, NMapas, NPixP, -1.0D0, YUBI, 1, 1, DescBloqueArmonicos, IMCMapas, Fila, 1,&
       & DescIMCMapas, 0.0D0, PBI, 1, 1, DescProducto)

  PBI = PBI + P1

  CALL blacs_barrier(ICTXT, "All")

  DEALLOCATE(P1)


END SUBROUTINE CalculaYDagIMCMapas



!*******************

SUBROUTINE CalculaYlAuto(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI, Lugar)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFMapasA, Lmax, NCMapas, NMapas, NB
  REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR, MATI, MAER, MAEI, MABR, MABI
  CHARACTER(len=100), INTENT(in) :: Lugar

  !*************

  INTEGER :: iam, NPROW, NPCOL, MYROW, MYCOL, INFO
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: ListaYl
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2, Temp3, YlFinal

  INTEGER(kind=4), DIMENSION(9) :: DescTodos, DescRepartidos

  INTEGER :: i,j,l

  INTEGER, EXTERNAL :: INDXL2G, NUMROC
  EXTERNAL :: BLACS_GRIDINFO, dgesd2d, dgerv2d, DGSUM2D, blacs_barrier, descinit, pdgemr2d


  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  ALLOCATE(ListaYl(2:lMax,6,NCMapas))
  ListaYl = 0.0D0

  !TT EE BB
  ALLOCATE(Temp1(NFMapasA, NCMapas))
  ALLOCATE(Temp2(NFMapasA, NCMapas))
  ALLOCATE(Temp3(NFMapasA, NCMapas))
  Temp1 = 0
  Temp2 = 0
  Temp3 = 0

  Temp1 = MATR**2 + MATI**2
  Temp2 = MAER**2 + MAEI**2
  Temp3 = MABR**2 + MABI**2

  DO i=1,NFMapasA
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) THEN
        ListaYl(l,1,:) = ListaYl(l,1,:) + Temp1(i,:)
        ListaYl(l,2,:) = ListaYl(l,2,:) + Temp2(i,:)
        ListaYl(l,3,:) = ListaYl(l,3,:) + Temp3(i,:)
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  !TE TB EB
  Temp1 = 2*(MATR*MAER + MATI*MAEI)
  Temp2 = 2*(MATR*MABR + MATI*MABI)
  Temp3 = 2*(MAER*MABR + MAEI*MABI)

  DO i=1,NFMapasA
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) THEN
        ListaYl(l,4,:) = ListaYl(l,4,:) + Temp1(i,:)
        ListaYl(l,5,:) = ListaYl(l,5,:) + Temp2(i,:)
        ListaYl(l,6,:) = ListaYl(l,6,:) + Temp3(i,:)
     END IF
  END DO


  DEALLOCATE(Temp1)
  DEALLOCATE(Temp2)
  DEALLOCATE(Temp3)

  ListaYl = ListaYl/2.0D0

  CALL blacs_barrier( ICTXT, "A" )

  !Suma Dl en las columnas
  !La suma esta lanzada en una forma que no tiene en cuenta la forma de la variable ListaYl
  CALL DGSUM2D( ictxt, 'C', '1-tree', 6*(lMax-1)*NCMapas, 1, ListaYl, 6*(lMax-1)*NCMapas, -1, MYCOL)

  !Todos los procesadores de la misma columna ya tienen completos los Yl de su parte de los mapas
  !Los procesadores de la primera fila le envían los datos al procesador 0, y los guarda

  IF(iam==0) THEN
     ALLOCATE(YlFinal(6*(lmax-1), NMapas))
  ELSE
     ALLOCATE(YlFinal(1, NMapas))
  END IF

  CALL DESCINIT(DescTodos, 6*(lmax-1), NMapas, 6*(lmax-1), NMapas, 0, 0, ICTXT, 6*(lmax-1), INFO)
  CALL DESCINIT(DescRepartidos, 6*(lmax-1), NMapas, 6*(lmax-1), NB, 0, 0, ICTXT, 6*(lmax-1), INFO)

  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(6*(lmax-1), NMapas, ListaYl, 1, 1, DescRepartidos, YlFinal, 1, 1, DescTodos, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF (IAM==0) THEN
     OPEN(Unit=20, File=TRIM(Lugar)//"/CoupledPower.dat", Action="write")
     DO i=1,NMapas
        WRITE(20,*) YlFinal(:,i)
     END DO
     CLOSE(20)
  END IF

  DEALLOCATE(YlFinal)
  DEALLOCATE(ListaYl)

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE CalculaYlAuto


SUBROUTINE CalculaYlCross(ICTXT, NFMapasA, lmax, NCMapas, NMapas, NB, MATR, MATI, MAER, MAEI, MABR, MABI,&
     &MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2, Lugar)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NFMapasA, Lmax, NCMapas, NMapas, NB
  REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR, MATI, MAER, MAEI, MABR, MABI
  REAL(kind=8), DIMENSION(NFMapasA, NCMapas), INTENT(in) ::  MATR2, MATI2, MAER2, MAEI2, MABR2, MABI2
  CHARACTER(len=100), INTENT(in) :: Lugar

  !*************

  INTEGER :: iam, NPROW, NPCOL, MYROW, MYCOL, INFO
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:,:) :: ListaYl
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Temp1, Temp2, Temp3, YlFinal

  INTEGER(kind=4), DIMENSION(9) :: DescTodos, DescRepartidos

  INTEGER :: i,j,l

  INTEGER, EXTERNAL :: INDXL2G, NUMROC
  EXTERNAL :: BLACS_GRIDINFO, dgesd2d, dgerv2d, DGSUM2D, blacs_barrier, descinit, pdgemr2d


  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  ALLOCATE(ListaYl(2:lMax,6,NCMapas))
  ListaYl = 0.0D0

  !TT EE BB
  ALLOCATE(Temp1(NFMapasA, NCMapas))
  ALLOCATE(Temp2(NFMapasA, NCMapas))
  ALLOCATE(Temp3(NFMapasA, NCMapas))

  Temp1 = 0
  Temp2 = 0
  Temp3 = 0

  Temp1 = MATR*MATR2 + MATI*MATI2
  Temp2 = MAER*MAER2 + MAEI*MAEI2
  Temp3 = MABR*MABR2 + MABI*MABI2



  DO i=1,NFMapasA
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) THEN
        ListaYl(l,1,:) = ListaYl(l,1,:) + Temp1(i,:)
        ListaYl(l,2,:) = ListaYl(l,2,:) + Temp2(i,:)
        ListaYl(l,3,:) = ListaYl(l,3,:) + Temp3(i,:)
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  !TE TB EB
  Temp1 = MATR*MAER2 + MATI2*MAEI + MATR2*MAER + MATI*MAEI2
  Temp2 = MATR*MABR2 + MATI2*MABI + MATR2*MABR + MATI*MABI2
  Temp3 = MAER*MABR2 + MAEI2*MABI + MAER2*MABR + MAEI*MABI2

  DO i=1,NFMapasA
     j = INDXL2G(i, NB, MYROW, 0, NPROW)
     l = FLOOR(SQRT(3.0D0+j))
     IF(l<=lmax) THEN
        ListaYl(l,4,:) = ListaYl(l,4,:) + Temp1(i,:)
        ListaYl(l,5,:) = ListaYl(l,5,:) + Temp2(i,:)
        ListaYl(l,6,:) = ListaYl(l,6,:) + Temp3(i,:)
     END IF
  END DO



  DEALLOCATE(Temp1)
  DEALLOCATE(Temp2)
  DEALLOCATE(Temp3)


  ListaYl = ListaYl/2.0D0

  CALL blacs_barrier( ICTXT, "A" )

  !Suma Dl en las columnas
  !La suma esta lanzada en una forma que no tiene en cuenta la forma de la variable ListaYl
  CALL DGSUM2D( ictxt, 'C', '1-tree', 6*(lMax-1)*NCMapas, 1, ListaYl, 6*(lMax-1)*NCMapas, -1, MYCOL)

  !Todos los procesadores de la misma columna ya tienen completos los Yl de su parte de los mapas
  !Los procesadores de la primera fila le envían los datos al procesador 0, y los guarda

  CALL blacs_barrier( ICTXT, "A" )

  IF(iam==0) THEN
     ALLOCATE(YlFinal(6*(lmax-1), NMapas))
  ELSE
     ALLOCATE(YlFinal(1, NMapas))
  END IF

  CALL DESCINIT(DescTodos, 6*(lmax-1), NMapas, 6*(lmax-1), NMapas, 0, 0, ICTXT, 6*(lmax-1), INFO)
  CALL DESCINIT(DescRepartidos, 6*(lmax-1), NMapas, 6*(lmax-1), NB, 0, 0, ICTXT, 6*(lmax-1), INFO)

  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(6*(lmax-1), NMapas, ListaYl, 1, 1, DescRepartidos, YlFinal, 1, 1, DescTodos, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF (IAM==0) THEN
     OPEN(Unit=20, File=TRIM(Lugar)//"/CoupledPower.dat", Action="write")
     DO i=1,NMapas
        WRITE(20,*) YlFinal(:,i)
     END DO
     CLOSE(20)
  END IF

  DEALLOCATE(YlFinal)
  DEALLOCATE(ListaYl)

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE CalculaYlCross


SUBROUTINE CargaMultipolos(ICTXT, iam, Lugar, FileData, LMax, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, IAM, lmax
  CHARACTER(len=100), INTENT(in) :: Lugar, FileData
  REAL(kind=8), DIMENSION(0:lmax), INTENT(OUT) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB

  INTEGER :: l

  REAL(kind=8), DIMENSION(7) :: Temp

  EXTERNAL :: blacs_barrier, dgebs2d, dgebr2d

  DlTT = 1
  DlEE = 1
  DLBB = 1
  DlTE = 1
  DlTB = 0
  DLEB = 0

  IF(IAM==0) THEN

     OPEN(Unit=20, File=TRIM(Lugar)//"/"//TRIM(FileData), Action="read")
     DO l=0,lmax
        READ (20,*) temp
        DlTT(l) = temp(2)
        DlEE(l) = temp(3)
        DlBB(l) = temp(4)
        DlTE(l) = temp(5)
        DlTB(l) = temp(6)
        DlEB(l) = temp(7)
     END DO
     CLOSE(20)

  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlTT, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlTT, lmax+1, 0, 0)
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
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlTE, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlTE, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlTB, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlTB, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , lmax+1, 1, DlEB, lmax+1)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , lmax+1, 1, DlEB, lmax+1, 0, 0)
  END IF

  CALL blacs_barrier(ICTXT, "ALL")

END SUBROUTINE CargaMultipolos


SUBROUTINE CalculaBloquesMatrizArmonicosY(ICTXT, NB, NFP, NC, npix, lmax, z, Phi, BeamPWDl, YTTR, YTTI)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NFP, NC, npix, lmax
  REAL(kind=8), DIMENSION(NFP), INTENT(IN)   :: z, phi
  REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(in) :: BeamPWDl
  REAL(kind=8), DIMENSION(NFP,NC), INTENT(OUT) :: YTTR, YTTI

  !Polinomios y bloques
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Senos, Cosenos
  !   REAL(KIND=8) :: Factor1, Factor2

  !Armonicos esfericos
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: YlR, YlI

  !Descriptores
  INTEGER(kind=4), DIMENSION(9) :: DescYl, DescY
  !Tamaños
  INTEGER(kind=4) :: NUMROC, TotalArmonicos, NYColumna

  !Contadores
  INTEGER(kind=4) :: l,  m, indice, Columna, base

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

  INTEGER(kind=4) :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO


  EXTERNAL :: BLACS_GRIDINFO, DESCINIT, blacs_barrier, pdgemr2d

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  YTTR = 0
  YTTI = 0

  !Reserva memoria para los polinomios de Legendre
  ALLOCATE( P1Base1(NFP, 0:lmax))
  ALLOCATE( P1Base2(NFP, 0:lmax))
  ALLOCATE( Pl1(NFP, 0:lmax))
  P1Base1 = 0D0
  P1Base2 = 0D0
  Pl1 = 0D0

  !Memoria para funciones trigonometricas
  ALLOCATE(Cosenos(NFP))
  ALLOCATE(Senos(NFP))

  Cosenos = 0
  Senos = 0

  !Para almacenar temporalmente los armonicos y pasarlos a los procesadores
  ALLOCATE(YlR(NFP, 2*lmax+1))
  ALLOCATE(YlI(NFP, 2*lmax+1))

  YlR = 0.0D0
  YlI = 0.0D0

  !Para los bloques con el total de los armonicos
  TotalArmonicos = -3 + 2*lmax + lmax**2
  NYColumna = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))

  !************************************************

  YTTR = 0.0D0
  YTTI = 0.0D0

  CALL DESCINIT(DescYl, npix, NPROW*(2*lmax+1), NB, 2*lmax+1, 0, 0, ICTXT, NFP, INFO)
  CALL DESCINIT(DescY,  npix, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFP, INFO)

  !Prepara los valores iniciales de los elementos de la iteracion
  IF (MYCOL==0) THEN
     P1Base2(:,0) = 1.D0/SQRT(4D0 * PI)
     P1Base1(:,0) = SQRT(3D0/(4D0 * PI)) * z
     P1Base1(:,1) = -(1D0/2D0) * SQRT(3D0/(2*PI)) * SQRT(1-z**2)
  END IF

  !Guarda l=1

  Cosenos = COS(Phi)
  Senos   = SIN(Phi)

  CALL blacs_barrier( ICTXT, "A" )

  DO l=2,lMax

     IF(MYCOL==0) THEN
        !Avanza al siguiente polinomio
        Bucle_m: DO m=0,l-2
           Pl1(:,m) = z(:) * SQRT((4D0 * l**2-1D0)/(l**2-m**2)) * P1Base1(:,m) - &
                &SQRT(((1.0D0+2.0D0*l)*(l-m-1.0D0)*(l+m-1.0D0) * 1D0)/((2.0D0*l-3.0D0)*(l-m)*(l+m)*1D0)) * P1Base2(:,m)
        END DO Bucle_m

        !m = l-1
        Pl1(:,l-1) = z(:) * SQRT(1D0*(1+2*l)) * P1Base1(:,l-1) 

        !m = l
        Pl1(:,l) = - SQRT(((1 + 2*l)*(1 - z**2))/(2D0*l)) * P1Base1(:,l-1)


        base = l+1
        indice = base
        DO m=0,l

           Cosenos = COS(1.0D0 * m * Phi)
           Senos   = SIN(1.0D0 * m * Phi)

           YlR(:,indice)  = BeamPWDl(1,l) * Pl1(:,m) * Cosenos
           YlR(:, base-m) = (-1.0)**m * YlR(:,indice)

           YlI(:,indice)  = BeamPWDl(1,l) * Pl1(:,m) * Senos
           YlI(:, base-m) = -(-1.0)**m * YlI(:,indice)

           indice = indice + 1

        END DO

     END IF

     CALL blacs_barrier( ICTXT, "A" )

     Columna = -3 + l**2
     CALL pdgemr2d(npix, 2*l+1, YlR, 1, 1, DescYl, YTTR, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(npix, 2*l+1, YlI, 1, 1, DescYl, YTTI, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     IF(MYCOL==0) THEN
        !Adelanta un paso los polinomios base de la iteracion
        P1Base2(:,0:l-1) = P1Base1(:,0:l-1)
        P1Base1(:,0:l) = Pl1(:,0:l)
     END IF

  END DO

END SUBROUTINE CalculaBloquesMatrizArmonicosY

SUBROUTINE CalculaBloquesMatrizArmonicosYFC(ICTXT, NB, NFP, NC, npix, lmax, z, Phi, BeamPWDl, YTTC)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NFP, NC, npix, lmax
  REAL(kind=8), DIMENSION(NFP), INTENT(IN)   :: z, phi
  REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(in) :: BeamPWDl
  COMPLEX(kind=8), DIMENSION(NFP,NC), INTENT(OUT) :: YTTC

  !Polinomios y bloques
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Senos, Cosenos
  !   REAL(KIND=8) :: Factor1, Factor2

  !Armonicos esfericos
  COMPLEX(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: YlC

  !Descriptores
  INTEGER(kind=4), DIMENSION(9) :: DescYl, DescY
  !Tamaños
  INTEGER(kind=4) :: NUMROC, TotalArmonicos, NYColumna

  !Contadores
  INTEGER(kind=4) :: l,  m, indice, Columna, base

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

  INTEGER(kind=4) :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO


  EXTERNAL :: BLACS_GRIDINFO, DESCINIT, blacs_barrier, pzgemr2d

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  YTTC = 0

  !Reserva memoria para los polinomios de Legendre
  ALLOCATE( P1Base1(NFP, 0:lmax))
  ALLOCATE( P1Base2(NFP, 0:lmax))
  ALLOCATE( Pl1(NFP, 0:lmax))
  P1Base1 = 0D0
  P1Base2 = 0D0
  Pl1 = 0D0

  !Memoria para funciones trigonometricas
  ALLOCATE(Cosenos(NFP))
  ALLOCATE(Senos(NFP))

  Cosenos = 0
  Senos = 0

  !Para almacenar temporalmente los armonicos y pasarlos a los procesadores
  ALLOCATE(YlC(NFP, 2*lmax+1))

  YlC = 0.0D0

  !Para los bloques con el total de los armonicos
  TotalArmonicos = -3 + 2*lmax + lmax**2
  NYColumna = MAX(1,NUMROC(TotalArmonicos, NB, MYCOL, 0, NPCOL))


  !************************************************

  CALL DESCINIT(DescYl, npix, NPROW*(2*lmax+1), NB, 2*lmax+1, 0, 0, ICTXT, NFP, INFO)
  CALL DESCINIT(DescY,  npix, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFP, INFO)

  !Prepara los valores iniciales de los elementos de la iteracion
  IF (MYCOL==0) THEN
     P1Base2(:,0) = 1.D0/SQRT(4D0 * PI)
     P1Base1(:,0) = SQRT(3D0/(4D0 * PI)) * z
     P1Base1(:,1) = -(1D0/2D0) * SQRT(3D0/(2*PI)) * SQRT(1-z**2)
  END IF

  !Guarda l=1

  Cosenos = COS(Phi)
  Senos   = SIN(Phi)

  CALL blacs_barrier( ICTXT, "A" )

  DO l=2,lMax

     IF(MYCOL==0) THEN
        !Avanza al siguiente polinomio
        DO m=0,l-2
           Pl1(:,m) = z(:) * SQRT((4D0 * l**2-1D0)/(l**2-m**2)) * P1Base1(:,m) - &
                &SQRT(((1.0D0+2.0D0*l)*(l-m-1.0D0)*(l+m-1.0D0) * 1D0)/((2.0D0*l-3.0D0)*(l-m)*(l+m)*1D0)) * P1Base2(:,m)
        END DO

        !m = l-1
        Pl1(:,l-1) = z(:) * SQRT(1D0*(1+2*l)) * P1Base1(:,l-1) 

        !m = l
        Pl1(:,l) = - SQRT(((1 + 2*l)*(1 - z**2))/(2D0*l)) * P1Base1(:,l-1)


        base = l+1
        indice = base
        DO m=0,l

           Cosenos = COS(1.0D0 * m * Phi)
           Senos   = SIN(1.0D0 * m * Phi)

           YlC(:,indice)  = dcmplx(BeamPWDl(1,l)*Pl1(:,m) * Cosenos, BeamPWDl(1,l)*Pl1(:,m) * Senos)
           YlC(:, base-m) = (-1.0)**m * CONJG(YlC(:,indice))

           indice = indice + 1

        END DO

     END IF

     CALL blacs_barrier( ICTXT, "A" )

     Columna = -3 + l**2
     CALL pzgemr2d(npix, 2*l+1, YlC, 1, 1, DescYl, YTTC, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )


     IF(MYCOL==0) THEN
        !Adelanta un paso los polinomios base de la iteracion
        P1Base2(:,0:l-1) = P1Base1(:,0:l-1)
        P1Base1(:,0:l) = Pl1(:,0:l)
     END IF

  END DO


END SUBROUTINE CalculaBloquesMatrizArmonicosYFC


SUBROUTINE CalculaBloquesMatrizArmonicosX(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
  REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
  REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
  REAL(kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI

  !Polinomios y bloques
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: X1,X2
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Alpha1, Alpha2, Beta1, Beta2, Senos, Cosenos
  REAL(KIND=8) :: Factor1, Factor2

  !Descriptores
  INTEGER, DIMENSION(9) :: DescYl, DescY
  !Tamaños
  INTEGER :: TotalArmonicos


  !Bloques intermedios
  REAL(kind=8), ALLOCATABLE,  DIMENSION(:,:) :: XQER, XQEI, XUER, XUEI, XQBR, XQBI, XUBR, XUBI

  !Contadores
  INTEGER :: l, m, indice, Columna, base

  !Grid
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pdgemr2d

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  YQER=0
  YQEI=0
  YQBR=0
  YQBI=0
  YUER=0
  YUEI=0
  YUBR=0
  YUBI=0

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
           Factor1 = BeamPWDl(2,l) * Factor1 

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
        Factor1 = BeamPWDl(2,l) * Factor1 

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
        Factor1 = BeamPWDl(2,l) * Factor1 

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

END SUBROUTINE CalculaBloquesMatrizArmonicosX

SUBROUTINE CalculaBloquesMatrizArmonicosXFC(ICTXT, NB, NF, NC, npix, lmax, z, Phi, BeamPWDl, YQEC, YQBC, YUEC, YUBC)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NF, NC, npix, lmax
  REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
  REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
  COMPLEX(kind=8), DIMENSION(NF,NC), INTENT(OUT) :: YQEC, YQBC, YUEC, YUBC

  !Polinomios y bloques
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: X1,X2
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Alpha1, Alpha2, Beta1, Beta2, Senos, Cosenos
  REAL(KIND=8) :: Factor1, Factor2

  !Descriptores
  INTEGER, DIMENSION(9) :: DescYl, DescY
  !Tamaños
  INTEGER :: TotalArmonicos

  !Bloques intermedios
  COMPLEX(kind=8), ALLOCATABLE,  DIMENSION(:,:) :: XQEC, XUEC, XQBC, XUBC

  !Contadores
  INTEGER :: l, m, indice, Columna, base

  !Grid
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pzgemr2d

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

  !Codigo
  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  YQEC=0
  YQBC=0
  YUEC=0
  YUBC=0

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
  ALLOCATE(XQEC(NF, 2*lmax+1))
  ALLOCATE(XUEC(NF, 2*lmax+1))
  ALLOCATE(XQBC(NF, 2*lmax+1))
  ALLOCATE(XUBC(NF, 2*lmax+1))

  XQEC = 0.0D0
  XQBC = 0.0D0
  XUEC = 0.0D0
  XUBC = 0.0D0


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
           Factor1 = BeamPWDl(2,l) * Factor1 

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
        Factor1 = BeamPWDl(2,l) * Factor1 

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
        Factor1 = BeamPWDl(2,l) * Factor1 

        X1(:,m) = Alpha1 * Factor1 * Pl1(:,m)
        X2(:,m) = Alpha2 * Factor1 * Pl1(:,m)
        !Ya tiene el siguiente polinomio

        base = l+1
        indice = base
        Bucle_m2: DO m=0,l

           Cosenos = COS(1.0D0 * m * Phi)
           Senos   = SIN(1.0D0 * m * Phi)

           XQEC(:,indice) = dcmplx( -X1(:,m) * Cosenos, -X1(:,m) * Senos)
           XQEC(:,base-m) = (-1.0)**m * CONJG(XQEC(:,indice))

           XQBC(:,indice) = dcmplx( X2(:,m) * Senos, -X2(:,m) * Cosenos)
           XQBC(:,base-m) = (-1.0)**m * CONJG(XQBC(:,indice))

           XUEC(:,indice) = dcmplx( -X2(:,m) * Senos,  X2(:,m) * Cosenos)
           XUEC(:,base-m) = (-1.0)**m * CONJG(XUEC(:,indice))

           XUBC(:,indice) = dcmplx(-X1(:,m) * Cosenos, -X1(:,m) * Senos)
           XUBC(:,base-m) = (-1.0)**m * CONJG(XUBC(:,indice))

           indice = indice + 1

        END DO Bucle_m2

     END IF

     CALL blacs_barrier( ICTXT, "A" )

     Columna = -3 + l*l

     CALL pzgemr2d(NPix, 2*l+1, XQEC, 1, 1, DescYl, YQEC, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XQBC, 1, 1, DescYl, YQBC, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XUEC, 1, 1, DescYl, YUEC, 1, Columna , DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     CALL pzgemr2d(NPix, 2*l+1, XUBC, 1, 1, DescYl, YUBC, 1, Columna, DescY, ictxt)
     CALL blacs_barrier( ICTXT, "A" )

     IF(MYCOL==0) THEN
        !Adelanta un paso los polinomios base de la iteracion
        P1Base2(:,0:l-1) = P1Base1(:,0:l-1)
        P1Base1(:,0:l) = Pl1(:,0:l)
     END IF

  END DO

END SUBROUTINE CalculaBloquesMatrizArmonicosXFC


SUBROUTINE CalculaMatrizCovarianza(ICTXT, NB, NFMC, NCMC, NFBY, NCBY, NFBX, lmax, NPixT, NPixP,&
     & MapaRuidoT2Bloque, MapaRuidoP2Bloque, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB,&
     & YTTR, YTTI, YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI, ControlInversaMC, MC, TiempoInicio)

  USE ifport

  IMPLICIT NONE

  INTEGER, INTENT(in) :: ICTXT, NB, NFMC, NCMC, NFBY, NCBY, NFBX, lmax, NPixT, NPixP
  INTEGER, INTENT(in) :: ControlInversaMC
  INTEGER(kind=8), INTENT(in) ::  TiempoInicio
  REAL(kind=8), DIMENSION(NFBY), INTENT(in) :: MapaRuidoT2Bloque
  REAL(kind=8), DIMENSION(NFBX), INTENT(in) :: MapaRuidoP2Bloque
  REAL(kind=8), DIMENSION(0:lmax), INTENT(in) ::  DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
  REAL(KIND=8), DIMENSION(NFBY,NCBY), INTENT(in)  :: YTTR, YTTI
  REAL(KIND=8), DIMENSION(NFBX,NCBY), INTENT(in)  :: YQER, YQEI, YQBR, YQBI, YUER, YUEI, YUBR, YUBI
  REAL(KIND=8), DIMENSION(NFMC,NCMC), INTENT(out)  :: MC

  !************************+

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

  INTEGER :: iam, TotalArmonicos
  INTEGER :: i, j, indice, l

  !Datos sobre los bloques QQ, UU, QU, UQ de la matriz de covarianza
  INTEGER :: NCBTraspuesto, NFLocalBloqueQQ, NCLocalBloqueQQ

  !Descriptores
  INTEGER, DIMENSION(9) :: DescBloqueArmonicos, DescBTraspuesto, DescBA2, DescProducto, DescMC

  !Tamaños
  INTEGER :: NUMROC, INDXL2G

  INTEGER :: NCBMC, DimMC
  INTEGER :: PuntoF, PuntoC
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesFMC, IndicesCMC, IndicesL

  INTEGER :: NPROW, NPCOL, MYROW, MYCOL, INFO

  !Para productos intermedios
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MTemp, MTemp2, MTemp3
  REAL(kind=8) :: Factor

  !Componentes TB y EB
  INTEGER :: ContribucionTB, ContribucionEB

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pdgemr2d, pdgemm, pdtran

  !**********************************************************************************

  MC = 0.0D0

  ContribucionTB = 0
  l = 2
  DO WHILE ((ContribucionTB == 0).AND.(l<=lmax))       
     IF(DlTB(l)/=0) ContribucionTB = 1
     l = l + 1      
  END DO

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

  !Bloque TT

  NCBMC = MAX(1,NUMROC(NPixT, NB, MYCOL, 0, NPCOL))

  !Indices de las filas y columnas del bloque TT de la matriz de covarianza
  ALLOCATE(IndicesFMC(NFBY))
  ALLOCATE(IndicesCMC(NCBMC))

  !Valores de l que corresponden a las columnas de los bloques de la matriz de armónicos
  ALLOCATE(IndicesL(NCBY))

  DO i=1,NFBY
     IndicesFMC(i) = INDXL2G(i, NB, MYROW, 0, NPROW)
  END DO

  DO i=1,NCBMC
     IndicesCMC(i) = INDXL2G(i, NB, MYCOL, 0, NPCOL)
  END DO

  !Recorre las columnas de los bloques de armonicos
  DO i=1,NCBY
     !Encuentra la columna de la matriz global a la que pertence esta columna local
     j = INDXL2G(i, NB, MYCOL, 0, NPCOL)
     !Encuentra el l del armonico esferico de esa columna
     Factor = 3.0D0 + 1.0D0 * j
     IndicesL(i) = FLOOR(SQRT(Factor))
  END DO

  !Dimensiones de bloque de matriz de armonicos
  ALLOCATE(MTemp(NFBY, NCBY))
  MTemp = 0

  !Dimensiones del bloque de la matriz de covarianza
  ALLOCATE(MTemp2(NFBY, NCBMC))
  MTemp2 = 0

  !Descriptor del bloque TTTT
  CALL DESCINIT(DescBloqueArmonicos, NPixT, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBY, INFO)

  !Descriptor del bloque de la matriz de covarianza resultado del producto
  CALL DESCINIT(DescProducto, NPixT, NPixT, NB, NB, 0, 0, ICTXT, NFBY, INFO)

  !Primer sumando
  MTemp = YTTR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTT(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixT, NPixT, TotalArmonicos, 1.0D0, YTTR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YTTI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTT(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixT, NPixT, TotalArmonicos, 1.0D0, YTTI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Suma ruido
  DO i=1,NFBY
     indice = IndicesFMC(i)
     DO j=1,NCBMC
        IF(indice==IndicesCMC(j)) THEN
           MTemp2(i,j) = MTemp2(i,j) + MapaRuidoT2Bloque(i)
        END IF
     END DO
  END DO

  DEALLOCATE(IndicesFMC)
  DEALLOCATE(IndicesCMC)
  DEALLOCATE(MTemp)

  !Pasa el bloque a la matriz de covarianza

  DimMC = NPixT + 2*NPixP
  CALL DESCINIT(DescMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL pdgemr2d(NPixT, NPixT, MTemp2, 1, 1, DescProducto, MC, 1, 1 , DescMC, ictxt)
  DEALLOCATE(MTemp2)

  IF(iam==0) WRITE(*,*) "    Block TT done ", Time()-TiempoInicio, "s"

  !Fin TT

  !******
  !Caso TQ
  !Real
  !Dimensiones de bloque de matriz de armonicos
  ALLOCATE(MTemp(NFBY, NCBY))
  Mtemp = 0

  !Dimensiones del bloque de la matriz de covarianza de tipo TQ
  NCBMC = MAX(1,NUMROC(NPixP, NB, MYCOL, 0, NPCOL))
  ALLOCATE(MTemp2(NFBY, NCBMC))
  MTemp2 = 0

  MTemp = YQER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionTB==1) THEN
           Factor = DlTB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBR(:,i)
        END IF
     END IF
  END DO

  CALL DESCINIT(DescBA2, NPixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL DESCINIT(DescProducto, NPixT, NPixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemm("N", "T", NPixT, NPixP, TotalArmonicos, 1.0D0, YTTR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBA2, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YQEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionTB==1) THEN
           Factor = DlTB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixT, NPixP, TotalArmonicos, 1.0D0, YTTI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBA2, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Lo envia a la matriz de covarianza
  PuntoF = 1
  PuntoC = NPixT + 1 
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NPixT, NPixP, MTemp2, 1, 1, DescProducto, MC, PuntoF, PuntoC, DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  DEALLOCATE(MTemp)

  IF(iam==0) WRITE(*,*)   "    Block TQ done ", Time()-TiempoInicio, "s"

  !Traspuesta 
  IF(ControlInversaMC==1) THEN

     NCBTraspuesto = MAX(1,NUMROC(NPixT, NB, MYCOL, 0, NPCOL))
     ALLOCATE(MTemp(NFBX, NCBTraspuesto))
     CALL DESCINIT(DescBTraspuesto, NPixP, NPixT, NB, NB, 0, 0, ICTXT, NFBX, INFO)

     CALL blacs_barrier( ICTXT, "A" )
     CALL pdtran(NPixP, NPixT, 1.0D0, MTemp2, 1, 1, DescProducto, 0.0D0, MTemp, 1, 1, DescBTraspuesto)

     !Envia al bloque inferior
     PuntoF = NPixT + 1
     PuntoC = 1 
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPixP, NPixT, MTemp, 1, 1, DescBTraspuesto, MC, PuntoF, PuntoC , DescMC, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     DEALLOCATE(MTemp)
     IF(iam==0) WRITE(*,*) "    Block QT done ", Time()-TiempoInicio, "s"                          
  END IF

  DEALLOCATE(MTemp2)

  !Fin TQ


  !******
  !Caso TU
  !Real
  !Dimensiones de bloque de matriz de armonicos
  ALLOCATE(MTemp(NFBY, NCBY))

  !Dimensiones del bloque de la matriz de covarianza de tipo TQ
  NCBMC = MAX(1,NUMROC(NPixP, NB, MYCOL, 0, NPCOL))
  ALLOCATE(MTemp2(NFBY, NCBMC))

  MTemp = YUER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionTB==1) THEN
           Factor = DlTB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBR(:,i)
        END IF
     END IF
  END DO

  CALL DESCINIT(DescBA2, NPixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFBX, INFO)
  CALL DESCINIT(DescProducto, NPixT, NPixP, NB, NB, 0, 0, ICTXT, NFBY, INFO)
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemm("N", "T", NPixT, NPixP, TotalArmonicos, 1.0D0, YTTR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBA2, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlTE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionTB==1) THEN
           Factor = DlTB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixT, NPixP, TotalArmonicos, 1.0D0, YTTI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBA2, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Lo envia a la matriz de covarianza
  PuntoF = 1
  PuntoC = NPixT + NPixP + 1 
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NPixT, NPixP, MTemp2, 1, 1, DescProducto, MC, PuntoF, PuntoC, DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  DEALLOCATE(MTemp)

  IF(iam==0) WRITE(*,*) "    Block TU done ", Time()-TiempoInicio, "s"

  !Traspuesta 
  IF(ControlInversaMC==1) THEN

     NCBTraspuesto = MAX(1,NUMROC(NPixT, NB, MYCOL, 0, NPCOL))
     ALLOCATE(MTemp(NFBX, NCBTraspuesto))
     MTemp = 0
     CALL DESCINIT(DescBTraspuesto, NPixP, NPixT, NB, NB, 0, 0, ICTXT, NFBX, INFO)
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdtran(NPixP, NPixT, 1.0D0, MTemp2, 1, 1, DescProducto, 0.0D0, MTemp, 1, 1, DescBTraspuesto)

     !Envia al bloque inferior
     PuntoF = NPixT + NPixP + 1
     PuntoC = 1 
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPixP, NPixT, MTemp, 1, 1, DescBTraspuesto, MC, PuntoF, PuntoC , DescMC, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     DEALLOCATE(MTemp)
     IF(iam==0) WRITE(*,*) "    Block UT done ", Time()-TiempoInicio, "s"
  END IF

  DEALLOCATE(MTemp2)

  !Fin TU


  !******************

  NFLocalBloqueQQ = MAX(1,NUMROC(NPixP, NB, MYROW, 0, NPROW))
  NCLocalBloqueQQ = MAX(1,NUMROC(NPixP, NB, MYCOL, 0, NPCOL))

  !Indices de las filas y columnas de los bloques de la matriz de covarianza
  ALLOCATE(IndicesFMC(NFLocalBloqueQQ))
  ALLOCATE(IndicesCMC(NCLocalBloqueQQ))

  DO i=1,NFLocalBloqueQQ
     IndicesFMC(i) = INDXL2G(i, NB, MYROW, 0, NPROW)
  END DO

  DO i=1,NCLocalBloqueQQ
     IndicesCMC(i) = INDXL2G(i, NB, MYCOL, 0, NPCOL)
  END DO

  ALLOCATE(MTemp(NFBX, NCBY))
  MTemp = 0

  !Dimensiones del bloque de la matriz de covarianza
  ALLOCATE(MTemp2(NFLocalBloqueQQ, NCLocalBloqueQQ))
  ALLOCATE(MTemp3(NFLocalBloqueQQ, NCLocalBloqueQQ))
  !allocate(MTemp4(NFLocalBloqueQQ, NCLocalBloqueQQ))
  MTemp2 = 0
  MTemp3 = 0

  !Descriptor del bloque de armonicos: mismas filas que bloque tipo QQ
  CALL DESCINIT(DescBloqueArmonicos, NPixP, TotalArmonicos, NB, NB, 0, 0, ICTXT, NFLocalBloqueQQ, INFO)

  !Descriptor del bloque producto: bloque tipo QQ
  CALL DESCINIT(DescProducto, NPixP, NPixP, NB, NB, 0, 0, ICTXT, NFLocalBloqueQQ, INFO)

  !Descriptor de la matriz de covarianza completa
  CALL DESCINIT(DescMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)


  !Caso QQ ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YQER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBR(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YQEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Segundo sumando
  !Real
  MTemp = YQBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YQBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YQEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp3, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  MTemp3 = MTemp3 + Mtemp2

  !Suma ruido
  DO i=1,NFLocalBloqueQQ
     indice = IndicesFMC(i)
     DO j=1,NCLocalBloqueQQ
        IF(indice==IndicesCMC(j)) THEN
           MTemp3(i,j) = MTemp3(i,j) + MapaRuidoP2Bloque(i)
        END IF
     END DO
  END DO

  !Lo envia a la matriz de covarianza
  CALL blacs_barrier( ICTXT, "A" )
  PuntoF = NPixT+1
  PuntoC = NPixT+1
  CALL pdgemr2d(NPixP, NPixP, MTemp3, 1, 1, DescProducto, MC, PuntoF, PuntoC, DescMC, ictxt)

  IF(iam==0) WRITE(*,*) "    Block QQ done ", Time()-TiempoInicio, "s"

  CALL blacs_barrier( ICTXT, "A" )


  !Caso QU ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YUER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
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

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Segundo sumnado
  !Real
  MTemp = YUBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YQBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp3, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  MTemp3 = MTemp3 + Mtemp2

  !Lo envia a la matriz de covarianza
  PuntoF = NPixT + 1
  PuntoC = NPixT + NPixP + 1 
  CALL blacs_barrier( ICTXT, "A" )
  CALL pdgemr2d(NPixP, NPixP, MTemp3, 1, 1, DescProducto, MC, PuntoF, PuntoC, DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF(iam==0) WRITE(*,*) "    Block QU done ", Time()-TiempoInicio, "s"

  !Traspuesta
  IF(ControlInversaMC==1) THEN
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdtran(NPixP, NPixP, 1.0D0, MTemp3, 1, 1, DescProducto, 0.0D0, MTemp2, 1, 1, DescProducto)
     !Envia al bloque inferior
     PuntoC = NPixT + 1
     PuntoF = NPixT + NPixP + 1 
     CALL blacs_barrier( ICTXT, "A" )
     CALL pdgemr2d(NPixP, NPixP, MTemp2, 1, 1, DescProducto, MC, PuntoF, PuntoC , DescMC, ictxt)
     CALL blacs_barrier( ICTXT, "A" )
     IF(iam==0) WRITE(*,*) "    Block UQ done ", Time()-TiempoInicio, "s"
  END IF



  !Caso UU ********************************************************************************

  !Primer sumando
  !Real
  MTemp = YUER
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBR(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YUER, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUEI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlEE(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUBI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YUEI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp2, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Segundo sumnado
  !Real
  MTemp = YUBR
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUER(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YUBR, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 0.0D0, MTemp3, 1, 1, DescProducto)

  CALL blacs_barrier( ICTXT, "A" )

  !Imaginaria
  MTemp = YUBI
  !Multiplica las columnas por el espectro de potencia
  DO i=1,NCBY
     l = IndicesL(i)
     IF(l<=lmax) THEN
        Factor = DlBB(IndicesL(i))
        MTemp(:,i) = Factor*MTemp(:,i)
        IF(ContribucionEB==1) THEN
           Factor = DlEB(IndicesL(i))
           MTemp(:,i) = MTemp(:,i) + Factor*YUEI(:,i)
        END IF
     END IF
  END DO

  CALL blacs_barrier( ICTXT, "A" )

  CALL pdgemm("N", "T", NPixP, NPixP, TotalArmonicos, 1.0D0, YUBI, 1, 1, DescBloqueArmonicos, MTemp, 1, 1,&
       & DescBloqueArmonicos, 1.0D0, MTemp3, 1, 1, DescProducto)

  MTemp3 = Mtemp3 + Mtemp2

  !Suma ruido
  DO i=1,NFLocalBloqueQQ
     indice = IndicesFMC(i)
     DO j=1,NCLocalBloqueQQ
        IF(indice==IndicesCMC(j)) THEN
           MTemp3(i,j) = MTemp3(i,j) + MapaRuidoP2Bloque(i)
        END IF
     END DO
  END DO

  !Lo envia a la matriz de covarianza
  CALL blacs_barrier( ICTXT, "A" )
  PuntoF = NPixT + NPixP + 1
  PuntoC = NPixT + NPixP + 1 
  CALL pdgemr2d(NPixP, NPixP, MTemp3, 1, 1, DescProducto, MC,&
       & PuntoF, PuntoC, DescMC, ictxt)
  CALL blacs_barrier( ICTXT, "A" )

  IF(iam==0) WRITE(*,*) "    Block UU done ", Time()-TiempoInicio, "s"

  DEALLOCATE(IndicesFMC)
  DEALLOCATE(IndicesCMC)

  DEALLOCATE(MTemp)
  DEALLOCATE(MTemp2)
  DEALLOCATE(MTemp3)

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
  REAL(kind=8), DIMENSION(NF,NC), INTENT(inout) :: MC
  INTEGER, INTENT(out) :: Fallo
  INTEGER(kind=8), INTENT(in) ::  TiempoInicio

  !********************************************************

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: IMC
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IndicesC
  INTEGER, DIMENSION(9) :: DESC, DescFilaRepartida, DescFila, DescColumnaRepartida, DescColumna
  INTEGER :: info
  REAL(kind=8) :: ProductoDiagonal, SumaDiagonal
  INTEGER ::  IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i,j, ii, jj, Seguir
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: SumaFila, SumaColumna, TodasColumnas, TodasFilas

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
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Inversion result:       ", info, Time()-TiempoInicio
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
        WRITE(*,'(4X, A, I6,  T79, I10, "s")') "Inversion result:       ", info, Time()-TiempoInicio
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
  REAL(kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
  REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  INTEGER, DIMENSION(9) :: DESCIMC, DESCMapas, DescTodos, DescDistribuidos
  INTEGER :: info
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Producto
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: SumaParte, RecibeSumas, ListaChi2
  REAL(kind=8) :: SumaChi2
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER, EXTERNAL :: NUMROC

  EXTERNAL :: BLACS_GRIDINFO, descinit, pdsymm, DGSUM2D, dgerv2d, dgesd2d, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  CALL descinit(DESCIMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  ALLOCATE(Producto(NFMC, NCMapas))
  Producto = 0

  CALL blacs_barrier(ICTXT, "All")
  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, Mapas, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)
  CALL blacs_barrier(ICTXT, "All")

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
     WRITE(*,*)
     WRITE(*,*) "    <m^t C^-1 m>: ", SumaChi2
     WRITE(*,*) "    Matrix size:  ", DimMC
  END IF

  !Devuelve IMC.Mapas
  Mapas = Producto

  DEALLOCATE(Producto)
  DEALLOCATE(SumaParte)
  DEALLOCATE(ListaChi2)

  CALL blacs_barrier(ICTXT, 'All')

END SUBROUTINE CalculaIMCMapas


SUBROUTINE CalculaMatrizArmonicosCompleja(ICTXT, NB, NF, NFMX, NCMX, npix, lmax, z, Phi, BeamPWDl, MXC)

  !NF : Numero de pixeles por procesador

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ICTXT, NB, NF, NFMX, NCMX, npix, lmax
  REAL(kind=8), DIMENSION(NF), INTENT(IN)   :: z, phi
  REAL(kind=8), DIMENSION(0:lmax), INTENT(in) :: BeamPWDl
  COMPLEX(kind=8), DIMENSION(NFMX,NCMX), INTENT(OUT) :: MXC

  !Polinomios y bloques
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: P1Base1, P1Base2, Pl1
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: X1,X2
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Alpha1, Alpha2, Beta1, Beta2, Senos, Cosenos
  REAL(KIND=8) :: Factor1, Factor2

  !Descriptores
  INTEGER, DIMENSION(9) :: DescYl, DescMX
  !Tamaños
  INTEGER :: TotalArmonicos


  !Bloques intermedios
  COMPLEX(kind=8), ALLOCATABLE,  DIMENSION(:,:) :: XQE, XQB, XUE, XUB

  !Contadores
  INTEGER :: l, m, indice, Columna, base

  !Grid
  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO

  EXTERNAL :: BLACS_GRIDINFO,  DESCINIT,  blacs_barrier, pzgemr2d

  REAL(KIND=8), PARAMETER :: Pi = 3.141592653589793D0

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


SUBROUTINE CalculaEspectroPotencia(IAM, Lugar, lmax, NMapas, QuitarSesgo)

  IMPLICIT NONE

  INTEGER(kind=4), INTENT(in) :: IAM, lmax, NMapas, QuitarSesgo
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  !*****************************

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: MF, IMF, Yl, Dl
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: Bl, DlMedio, DlSigma, SumaDl, SumaCuardadoDl, ErrorTeorico

  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: DlEscalar, DlTensorial, P1, r
  REAL(kind=8) :: divisor, RMedio, SigmaR, ErrorR

  INTEGER :: info, dim, i, AllocateStatus

  EXTERNAL :: dpotrf, dpotri, dsymm

  dim = 6*(lmax-1)

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
        RETURN
     END IF

     CALL dpotri("U", dim, IMF, dim, info)
     IF(IAM==0) THEN
        WRITE(*,*) " Inversion result:       ", info
     END IF

     IF(info/=0) THEN
        WRITE(*,*)
        WRITE(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        WRITE(*,*) "Fisher matrix is singular"
        RETURN
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
  REAL(kind=8), DIMENSION(NFMC,NCMC), INTENT(IN) :: IMC
  REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(inout) :: Mapas, MapasCross
  CHARACTER(LEN=100), INTENT(in) :: Lugar

  INTEGER, DIMENSION(9) :: DESCIMC, DESCMapas, DescTodos, DescDistribuidos
  INTEGER :: info
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Producto
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: SumaParte, RecibeSumas, ListaChi2
  REAL(kind=8) :: SumaChi2
  INTEGER :: NPROW, NPCOL, MYROW, MYCOL

  EXTERNAL :: BLACS_GRIDINFO, descinit, pdsymm, DGSUM2D, dgerv2d, dgesd2d, blacs_barrier, pdgemr2d

  CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL)

  CALL descinit(DESCIMC, DimMC, DimMC, NB, NB, 0, 0, ICTXT, NFMC, INFO)
  CALL descinit(DESCMapas, DimMC, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  ALLOCATE(Producto(NFMC, NCMapas))
  Producto = 0

  !Calcula IMC.Mapas
  CALL blacs_barrier(ICTXT, 'ALL')
  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, Mapas, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)
  !Devuelve IMC.Mapas
  Mapas = Producto


  !***************************************
  !Calcula Chi2

  !Multiplica MapasCross.(IMC.Mapas)
  Producto = MapasCross * Mapas

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
     WRITE(*,*) "    <m^t C^-1 m_cross>: ", SumaChi2
     WRITE(*,*) "    Matrix size:        ", DimMC
  END IF

  !Calcula IMC.MapasCross
  Producto = 0
  CALL blacs_barrier(ICTXT, 'ALL')
  CALL pdsymm("L", "U", DimMC, NMapas, 1.0D0, IMC, 1, 1, DESCIMC, MapasCross, 1, 1, DESCMapas, 0.0D0, Producto, 1, 1, DESCMapas)
  !Devuelve IMC.Mapas
  MapasCross = Producto

  DEALLOCATE(Producto)

  CALL blacs_barrier(ICTXT, 'ALL')

END SUBROUTINE CalculaIMCMapasMapasCross


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


SUBROUTINE CargaMascara(ictxt, nside, Lugar, File, NumMapaCargar, Mascara)

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


  INTEGER, INTENT(in) :: ictxt, Nside, NumMapaCargar
  CHARACTER(len=100), INTENT(in) :: Lugar, File
  INTEGER, DIMENSION(12*nside**2), INTENT(out) :: Mascara

  CHARACTER(len=200) :: ArchivoMascara
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: MascaraReal
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

     !open(unit=20, file=Trim(lugar)//"/Mascara.dat", action="write")
     !write(20,*) Mascara
     !close(20)

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

        DO imap = 1, nmaps

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


SUBROUTINE pix2ang_ring(nside, ipix, z, phi)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: nside
  INTEGER, INTENT(IN)  :: ipix
  REAL(KIND=8),     INTENT(OUT) :: z, phi

  INTEGER ::  nl2, nl4, iring, iphi
  INTEGER ::  npix, ncap, ip
  REAL(KIND=8) ::  fodd, dnside
  REAL(KIND=8), PARAMETER :: half = 0.500000000000000D0
  REAL(KIND=8), PARAMETER :: threehalf = 1.50000000000000D0

  REAL(kind=8), PARAMETER :: PI = 3.1415926535897932384626433832795D0
  REAL(kind=8), PARAMETER :: HALFPI = 1.5707963267948966192313216916398D0

  CHARACTER(len=*), PARAMETER :: code = "pix2ang_ring"
  !-----------------------------------------------------------------------

  npix = 12*nside**2       ! total number of points
  !if (ipix <0 .or. ipix>npix-1) call fatal_error (code//"> ipix out of range")

  nl2  = 2*nside
  ncap = nl2*(nside-1 ) ! points in each polar cap, =0 for nside =1
  dnside = REAL(nside, kind=8)

  IF (ipix < ncap) THEN ! North Polar cap -------------

     iring = (SQRT(2.0D0*ipix+2.0D0) + 1D0)/2D0
     iphi  = ipix - 2*iring*(iring - 1 )

     z = COS(2.0D0 * ASIN(iring / (SQRT(6.0D0)*dnside)))
     phi   = (REAL(iphi,kind=8) + half) * HALFPI/iring

  ELSEIF (ipix < npix-ncap) THEN ! Equatorial region ------

     ip    = ipix - ncap
     nl4   = 4*nside
     iring = INT( ip / nl4 ) + nside ! counted from North pole
     iphi  = IAND(ip, nl4-1 )

     fodd  = half * ( IAND(iring+nside+1,1) )  ! 0 if iring+nside is odd, 1/2 otherwise
     z =  (nl2 - iring) / (threehalf*dnside)
     phi   = (REAL(iphi,kind=8) + fodd) * HALFPI / dnside

  ELSE ! South Polar cap -----------------------------------

     ip    = npix - ipix
     iring = (SQRT(2D0*ip) + 1) / 2
     iphi  = 2*iring*(iring + 1 ) - ip


     z = COS(PI - 2.d0 * ASIN(iring / (SQRT(6.0D0)*dnside)))
     phi   = (REAL(iphi,kind=8) + half) * HALFPI/iring

  ENDIF

  RETURN

END SUBROUTINE pix2ang_ring


SUBROUTINE PuntosPixelesObservados(ictxt, nside, npix, Mascara, NB, NFBX, z, a)

  IMPLICIT NONE

  INTERFACE

     SUBROUTINE pix2ang_ring(nside, ipix, z, phi)
       INTEGER, INTENT(IN)  :: nside
       INTEGER, INTENT(IN)  :: ipix
       REAL(KIND=8),     INTENT(OUT) :: z, phi
     END SUBROUTINE pix2ang_ring

  END INTERFACE

  INTEGER, INTENT(in) :: ictxt, nside, npix, NB, NFBX
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: Mascara
  REAL(kind=8), DIMENSION(NFBX) :: z
  REAL(kind=8), DIMENSION(NFBX) :: a

  !***************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL, INFO
  INTEGER :: npix_full, ipix, j
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: TodosZ, TodosA
  REAL(kind=8) :: TempZ, TempA

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
  CALL descinit(DescReparto, npix, 1, NB, 1, 0, 0, ICTXT, NFBX, INFO)

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
       REAL(kind=8),                   INTENT(in)  :: fwhm_arcmin
       INTEGER,               INTENT(in)  :: lmax
       REAL(kind=8), DIMENSION(0:lmax,1:3), INTENT(out) :: gb
     END SUBROUTINE gaussbeam

     SUBROUTINE get_pixel_window(dir, nside, pixel)
       CHARACTER(len=*), INTENT(in) :: dir
       INTEGER, INTENT(in) :: nside
       REAL(kind=8), DIMENSION(0:4*NSide, 3), INTENT(out) :: pixel
     END SUBROUTINE get_pixel_window

  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, PixelWindowEnMapas, Nside, lmax
  REAL(kind=8) :: FWHM_Beam
  CHARACTER(len=*), INTENT(in) :: DirHealpixData                          !Lugar donde estan los datos Healpix
  REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(out) :: BeamPWDl

  !********

  INTEGER :: LmaxTemp, l
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Beam, Pixel
  REAL(kind=8), PARAMETER :: DosPi = 6.2831853071795864769D0
  REAL(kind=8) :: f

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL

  EXTERNAL ::  BLACS_GRIDINFO, blacs_barrier, dgebs2d, dgebr2d

  !***************************

  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  LmaxTemp = 4*nside

  BeamPWDl = 0

  IF(iam==0) THEN
     !Beam
     ALLOCATE(beam(0:LmaxTemp, 3))
     CALL gaussbeam(FWHM_Beam, LmaxTemp, beam)

     !Pixel_Window
     ALLOCATE(Pixel(0:LmaxTemp,3))

     IF(PixelWindowEnMapas==1) THEN
        CALL get_pixel_window(TRIM(DirHealpixData), nside, Pixel)
     ELSE
        Pixel = 1D0
     END IF

     !Prepara la salida de los datos
     DO l=1,lmax
        f = SQRT(DosPi/(REAL(l)*(REAL(l)+1)))
        BeamPWDl(1,l) = f * Beam(l,1) * pixel(l,1)
        BeamPWDl(2,l) = f * Beam(l,2) * pixel(l,2)
        BeamPWDl(3,l) = f * Beam(l,3) * pixel(l,3)
     END DO

  END IF

  CALL blacs_barrier(ictxt, "all")

  IF(IAM==0) THEN
     CALL dgebs2d(ICTXT, "All","I" , 3, Lmax+1, BeamPWDl, 3)
  ELSE
     CALL dgebr2d(ICTXT, "All","I" , 3, Lmax+1, BeamPWDl, 3, 0, 0)
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
  REAL(kind=8),                   INTENT(in)  :: fwhm_arcmin
  REAL(kind=8), DIMENSION(0:lmax,1:3), INTENT(out) :: gb
  INTEGER,               INTENT(in)  :: lmax

  INTEGER :: l, nd
  REAL(kind=8)     :: sigma, arcmin2rad, sigma2fwhm, fact_pol
  !===========================================================

  REAL(kind=8), PARAMETER :: PI = 3.1415926535897932384626433832795D0


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
       REAL(kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
       INTEGER,                   INTENT(IN) :: extno
     END SUBROUTINE read_bintab

  END INTERFACE

  CHARACTER(len=*), INTENT(in) :: dir
  INTEGER, INTENT(in) :: nside
  REAL(kind=8), DIMENSION(0:4*NSide, 3), INTENT(out) :: pixel

  !******************************

  CHARACTER(100) :: file
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: PixelTemp


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
  CASE(1024)
     file = "/pixel_window_n1024.fits"
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


SUBROUTINE CargaMapasRuido(ICTXT, TipoRuido, RuidoTT, RuidoQQ, nside, NFBY, NFBX, NFQUMC, NB,&
     & MascaraT, MascaraP, Lugar, FileData, MapaRuidoT2Bloque, MapaRuidoP2Bloque, MapaRuidoP2)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab_UnMapa(filename, MAPA, npixtot, imap, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: imap
       REAL(kind=8),      DIMENSION(0:npixtot-1),         INTENT(OUT) :: MAPA
       INTEGER                   , INTENT(IN) :: extno
     END SUBROUTINE read_bintab_UnMapa
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, TipoRuido, nside,  NFBY, NFBX, NFQUMC, NB
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: MascaraT, MascaraP
  CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
  REAL(kind=8), INTENT(inout) :: RuidoTT, RuidoQQ
  REAL(kind=8), DIMENSION(NFBY), INTENT(out) :: MapaRuidoT2Bloque
  REAL(kind=8), DIMENSION(NFBX), INTENT(out) :: MapaRuidoP2Bloque
  REAL(kind=8), DIMENSION(NFQUMC), INTENT(out) :: MapaRuidoP2


  !******************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i, k, imap, NpixFull, extno
  INTEGER :: npix
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: MapaDoble, MapaRuidoMascara
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: Mapa
  CHARACTER(len=200)  :: FileName


  EXTERNAL :: BLACS_GRIDINFO, descinit, blacs_barrier, pdgemrsd, dgebr2d, dgebs2d
  INTEGER, EXTERNAL :: INDXL2G

  IF(TipoRuido==0) THEN
     MapaRuidoT2Bloque = RuidoTT**2
     MapaRuidoP2Bloque = RuidoQQ**2
     MapaRuidoP2 = RuidoQQ**2
     RETURN
  END IF


  CALL BLACS_GRIDINFO(ICTXT, NPROW, NPCOL, MYROW, MYCOL)
  IAM = MYROW*NPCOL+MYCOL

  FileName = TRIM(Lugar)//"/"//TRIM(FileData)

  NpixFull = 12*nside**2

  ALLOCATE(Mapa(NpixFull,1))
  Mapa = 0


  extno = 0
  imap = 1

  !Carga mapa temperatura

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

  npix = SUM(MascaraT)
  ALLOCATE(MapaRuidoMascara(npix))

  k=1
  DO i=1,NpixFull
     IF(MascaraT(i)==1) THEN
        MapaRuidoMascara(k) = Mapa(i,1)
        k = k + 1
     END IF
  END DO

  !Calcula la media armonica, para tener el ruido por pixel equivalente
  !Lo usa para binear
  RuidoTT = REAL(npix)/SUM(1.0D0/MapaRuidoMascara)

  !Lo pasa a cuadrados
  MapaRuidoMascara = MapaRuidoMascara**2

  !Toma elementos
  DO i=1,NFBY
     MapaRuidoT2Bloque(i) = MapaRuidoMascara(INDXL2G(i, NB, MYROW, 0, NPROW))
  END DO

  DEALLOCATE(MapaRuidoMascara)

  !Carga mapa polarizacion             
  imap = 2

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


  npix = SUM(MascaraP)
  ALLOCATE(MapaRuidoMascara(npix))

  k=1
  DO i=1,NpixFull
     IF(MascaraP(i)==1) THEN
        MapaRuidoMascara(k) = Mapa(i,1)
        k = k + 1
     END IF
  END DO

  !Calcula la media armonica, para tener el ruido por pixel equivalente
  !Lo usa para binear
  RuidoQQ = REAL(npix)/SUM(1.0D0/MapaRuidoMascara)

  !Lo pasa a cuadrados
  MapaRuidoMascara = MapaRuidoMascara**2

  !Toma elementos
  DO i=1,NFBX
     MapaRuidoP2Bloque(i) = MapaRuidoMascara(INDXL2G(i, NB, MYROW, 0, NPROW))
  END DO

  ALLOCATE(MapaDoble(2*npix))
  MapaDoble(1:npix)        = MapaRuidoMascara(1:npix)
  MapaDoble(npix+1:2*npix) = MapaRuidoMascara(1:npix)

  !Toma elementos en MapaRuidoP2
  DO i=1,NFQUMC
     MapaRuidoP2(i) = MapaDoble(INDXL2G(i, NB, MYROW, 0, NPROW))
  END DO

  DEALLOCATE(Mapa)
  DEALLOCATE(Mapadoble)
  DEALLOCATE(MapaRuidoMascara)

  CALL blacs_barrier(ICTXT,"all")

END SUBROUTINE CargaMapasRuido


SUBROUTINE CargaMapasFits(ICTXT, nside, NMapas, NB, NFMC, NCMapas, MascaraT, MascaraP, Lugar, FileData, Mapas)

  USE ifport

  IMPLICIT NONE

  INTERFACE
     SUBROUTINE read_bintab(filename, map, npixtot, nmaps, extno)
       CHARACTER(len=*),                          INTENT(IN)  :: filename
       INTEGER,                              INTENT(IN)  :: npixtot
       INTEGER,                              INTENT(IN)  :: nmaps
       REAL(kind=8),      DIMENSION(0:npixtot-1,1:nmaps),         INTENT(OUT) :: map
       INTEGER, INTENT(IN) :: extno
     END SUBROUTINE read_bintab
  END INTERFACE

  INTEGER, INTENT(in) :: ICTXT, nside, NMapas, NB, NFMC, NCMapas
  INTEGER, DIMENSION(12*nside**2), INTENT(in) :: MascaraT, MascaraP
  CHARACTER(LEN=100), INTENT(in) :: Lugar, FileData
  REAL(kind=8), DIMENSION(NFMC,NCMapas), INTENT(out) :: Mapas

  !******************************************

  INTEGER :: IAM, NPROW, NPCOL, MYROW, MYCOL
  INTEGER :: i,j,k, imap, INFO, NpixFull, extno, numext
  INTEGER :: npix, npixt, npixp, ipixT, ipixQ, ipixU
  INTEGER :: NumMapasCargar, Faltan, Pasados
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MapasFullSkyT, MapasFullSkyQ, MapasFullSkyU, MapasMascara
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: MapasCargar, MapasPasarT, MapasPasarQ, MapasPasarU
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

  IF(iam==0) THEN
     ALLOCATE(MapasCargar(0:NpixFull-1, 3*NumMapasCargar))
     ALLOCATE(MapasPasarT(0:NpixFull-1, NumMapasCargar))
     ALLOCATE(MapaspasarQ(0:NpixFull-1, NumMapasCargar))
     ALLOCATE(MapasPasarU(0:NpixFull-1, NumMapasCargar))
  ELSE
     ALLOCATE(MapasCargar(1, 3*NumMapasCargar))
     ALLOCATE(MapasPasarT(1, NumMapasCargar))
     ALLOCATE(MapaspasarQ(1, NumMapasCargar))
     ALLOCATE(MapasPasarU(1, NumMapasCargar))
  END IF

  IF(MYROW==0) THEN
     ALLOCATE(MapasFullSkyT(NpixFull, NCMapas))
     ALLOCATE(MapasFullSkyQ(NpixFull, NCMapas))
     ALLOCATE(MapasFullSkyU(NpixFull, NCMapas))
  ELSE    
     ALLOCATE(MapasFullSkyT(1, NCMapas))
     ALLOCATE(MapasFullSkyQ(1, NCMapas))
     ALLOCATE(MapasFullSkyU(1, NCMapas))
  END IF

  !Descriptores para trabajar solo en la primera fila
  CALL descinit(DescCarga, NpixFull, NumMapasCargar, NpixFull, NumMapasCargar, 0, 0, ICTXT, NpixFull, INFO)
  CALL descinit(DescMapas, NpixFull, NMAPAS, NpixFull, NB, 0, 0, ICTXT, NpixFull, INFO)

  extno = 0
  imap = 1
  i=0
  Pasados = 0
  DO i=0,NumExt-1

     !if(iam==0) write(*,*) "Cargando mapas hasta: ", Pasados+50 

     !Carga
     IF(iam==0) THEN
        CALL read_bintab(filename, MapasCargar, NPixFull, 3*NumMapasCargar, extno) 

        k = 1
        DO j=1,50
           MapasPasarT(:,j) = MapasCargar(:,k)
           k = k + 1
           MapasPasarQ(:,j) = MapasCargar(:,k)
           k = k + 1
           MapasPasarU(:,j) = MapasCargar(:,k)
           k = k + 1
        ENDDO

     END IF


     !Pasa T
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarT, 1, 1, DescCarga, MapasFullSkyT, 1, Pasados+1, DescMapas, ictxt)

     !Pasa Q
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarQ, 1, 1, DescCarga, MapasFullSkyQ, 1, Pasados+1, DescMapas, ictxt)

     !Pasa U
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarU, 1, 1, DescCarga, MapasFullSkyU, 1, Pasados+1, DescMapas, ictxt)

     !Prepara indice para el siguiente mapa    
     extno = extno + 1

     Pasados = Pasados + 50

  END DO


  !Carga los mapas restantes
  IF(Pasados < NMapas) THEN

     Faltan = NMapas - (NumExt*50)

     ! if(iam==0) write(*,*) "Cargando de la ultima extension ", Faltan, " mapas"

     NumMapasCargar = Faltan

     CALL descinit(DescCarga, NpixFull, NumMapasCargar, NpixFull, NumMapasCargar, 0, 0, ICTXT, NpixFull, INFO)

     !Carga
     IF(iam==0) THEN
        CALL read_bintab(filename, MapasCargar, NPixFull, 3*NumMapasCargar, extno) 

        k = 1
        DO j=1,NumMapasCargar
           MapasPasarT(:,j) = MapasCargar(:,k)
           k = k + 1
           MapasPasarQ(:,j) = MapasCargar(:,k)
           k = k + 1
           MapasPasarU(:,j) = MapasCargar(:,k)
           k = k + 1
        ENDDO

     END IF


     !Pasa T
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarT, 1, 1, DescCarga, MapasFullSkyT, 1, Pasados+1, DescMapas, ictxt)

     !Pasa Q
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarQ, 1, 1, DescCarga, MapasFullSkyQ, 1, Pasados+1, DescMapas, ictxt)

     !Pasa U
     CALL blacs_barrier(ICTXT, "ALL")
     CALL pdgemr2d(NpixFull, NumMapasCargar, MapasPasarU, 1, 1, DescCarga, MapasFullSkyU, 1, Pasados+1, DescMapas, ictxt)

     Pasados = Pasados + Faltan

  END IF
  !Fin de carga los mapas rtestantes

  DEALLOCATE(MapasCargar)

  !Aplica mascara
  npixt = SUM(MascaraT)
  npixp = SUM(MascaraP)
  npix = npixt + 2*npixp
  IF(MYROW==0) THEN
     ALLOCATE(MapasMascara(npix, NCMapas))
  ELSE
     ALLOCATE(MapasMascara(1, NCMapas))
  END IF
  MapasMascara = 0



  ipixT = 1
  ipixQ = npixt + 1
  ipixU = npixt + npixp + 1
  IF(MYROW==0) THEN
     DO i=1,NpixFull
        IF(MascaraT(i)==1) THEN
           MapasMascara(ipixT,:) = MapasFullSkyT(i,:)
           ipixT = ipixT + 1
        END IF
        IF(MascaraP(i)==1) THEN
           MapasMascara(ipixQ,:) = MapasFullSkyQ(i,:)
           ipixQ = ipixQ + 1
           MapasMascara(ipixU,:) = MapasFullSkyU(i,:)
           ipixU = ipixU + 1 
        END IF
     END DO
  END IF

  CALL blacs_barrier(ICTXT, "ALL")


  DEALLOCATE(MapasFullSkyT)
  DEALLOCATE(MapasFullSkyQ)
  DEALLOCATE(MapasFullSkyU)


  !Reparte los mapas entre los procesadores de la columna
  Mapas = 0
  !Lo pasa todo a la vez
  CALL descinit(DescCarga, npix, NMapas, npix, NB, 0, 0, ICTXT, npix, INFO)
  CALL descinit(DescMapas, npix, NMapas, NB, NB, 0, 0, ICTXT, NFMC, INFO)

  CALL blacs_barrier(ICTXT, "ALL")
  !if(iam==0) write(*,*) "XXX: ", npix, nmapas, nb, nfmc
  CALL pdgemr2d(npix, NMapas, MapasMascara, 1, 1, DescCarga, Mapas, 1, 1, DescMapas, ictxt)
  CALL blacs_barrier(ICTXT, "ALL")

  DEALLOCATE(MapasMascara)


END SUBROUTINE CargaMapasFits


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
           CALL CargaEntero(cin, "Intensity_Mask_NumMap = ", NumMaskMapIntensity)
           CALL CargaEntero(cin, "Polarization_Mask_NumMap = ", NumMaskMapPolarization)
           CALL CargaEntero(cin, "Type_of_Noise = ", TipoRuido)
           CALL CargaEntero(cin, "Type_of_Noise_Data = ", TipoDatosMapaRuido)
           CALL CargaEntero(cin, "Remove_Noise_Bias = ", QuitarSesgo)
           CALL CargaEntero(cin, "Matrices_Cyclic_Block_Size = ", NB)
           CALL CargaEntero(cin, "Inverse_Covariance_Matrix_Control = ", ControlInversaMC)
           CALL CargaEntero(cin, "Show_Memory_Allocated = ", MuestraMemoria)
           CALL CargaEntero(cin, "Type_of_Bin_Center = ", KindOfBinCenter)
           CALL CargaEntero(cin, "Type_of_Grouping = ", KindOfGrouping)
           CALL CargaEntero(cin, "Compute_Fisher_Matrix = ", ComputeFisherMatrix)
           CALL CargaEntero(cin, "Binned = ", Binned)
           CALL CargaEntero(cin, "Compute_Spectrum = ", ComputeSpectrum)


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
     Paso(12) = KindOfBinCenter
     Paso(13) = KindOfGrouping
     Paso(14) = ComputeFisherMatrix
     Paso(15) = Binned
     Paso(16) = ComputeSpectrum

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
     KindOfBinCenter = Paso(12)
     KindOfGrouping = Paso(13)
     ComputeFisherMatrix = Paso(14)
     Binned = Paso(15)
     ComputeSpectrum = Paso(16)

  END IF

  CALL blacs_barrier(ICTXT, "ALL")


  !*********************
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
     WRITE(TEXTO,*) ComputeFisherMatrix
     WRITE(*,320) "Compute_Fisher_Matrix: ", ADJUSTL(TRIM(Texto))
     WRITE(TEXTO,*) ComputeSpectrum
     WRITE(*,320) "Compute_Spectrum: ", ADJUSTL(TRIM(Texto))
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

     WRITE(*,320) "Intensity_Mask_FileName: ", TRIM(FileMaskIntensity)
     WRITE(TEXTO,*) NumMaskMapIntensity
     WRITE(*,320) "Intensity_Mask_NumMap: ", ADJUSTL(TRIM(Texto))
     WRITE(*,320) "Polarization_Mask_FileName: ", TRIM(FileMaskPolarization)
     WRITE(TEXTO,*) NumMaskMapPolarization
     WRITE(*,320) "Polarization_Mask_NumMap: ", ADJUSTL(TRIM(Texto))
     WRITE(*,*)


     WRITE(TEXTO,*) TipoRuido
     WRITE(*,320) "Type_of_Noise: ", ADJUSTL(TRIM(Texto))
     IF(TipoRuido==0) THEN
        WRITE(TEXTO,*) RuidoTT
        WRITE(*,320) "Intensity_Noise:", ADJUSTL(TRIM(Texto))
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
       !write(*,*) "Toma entero de: ", cad
       READ(cad,*)  valor

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
       !write(*,*) "Toma real de: ", cad
       READ(cad,*)  valor

    END IF


  END SUBROUTINE CargaReal

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
       !write(*,*) "Toma real de: ", cad
       READ(cad,*)  valor

    END IF

  END SUBROUTINE CargaCadena

END SUBROUTINE CargaProblema


SUBROUTINE CalculaBineado(NSide, npixT, npixP, lmax, RuidoTT, RuidoQQ, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB, BeamPWDl,&
     & TipoCompactado, TipoCentrado, QuitarSesgo, NMapas, Lugar, FileBinLimits)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: lmax, NSide, npixT, npixP, QuitarSesgo, NMapas
  INTEGER, INTENT(in) :: TipoCompactado                                    !Compacta Fisher pesando o no por fiducial
  INTEGER, INTENT(in) :: TipoCentrado                                      !\ell del bin pesado por errores o \ell medio   
  REAL(kind=8), INTENT(in) :: RuidoTT, RuidoQQ
  REAL(kind=8), DIMENSION(0:lmax), INTENT(in) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
  REAL(kind=8), DIMENSION(3,0:lmax), INTENT(in) :: BeamPWDl
  CHARACTER(len=100), INTENT(in) :: Lugar, FileBinLimits


  !Internas
  CHARACTER(len=100) :: File


  INTEGER :: l, i, li, ls, Caso, f, c, dim
  REAL(kind=8) :: fskyT, fskyP, temp

  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: ErrorDlTT, ErrorDlEE, ErrorDlBB
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: ErrorDlTE, ErrorDlTB, ErrorDlEB

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: DLBin, LBin, Factores, MatrizR
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: SumaPesos, SumasDl, SumasL
  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: PesosCompactado

  REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) ::  Yl, MatrizProducto, MFBines, YlBines, Dl
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: Bl , MF
  REAL(kind=8), ALLOCATABLE, DIMENSION(:) :: DlMedio, DlSigma, SumaDl, SumaCuardadoDl, ErrorTeoricoDl

  REAL(kind=8) :: npixfull

  INTEGER :: Leer, Limite, NBines, info
  INTEGER, DIMENSION(1000) :: LimiteSuperiorBines

  EXTERNAL :: dsymm, dgemm, dpotrf, dpotri

  File = TRIM(Lugar)//"/"//TRIM(FileBinLimits)
  OPEN(unit=10, file=TRIM(File), action="read")
  NBines = 0

  Leer = 1
  DO WHILE(Leer==1)

     READ(10, *, END=100) Limite
     IF(Limite<lmax) THEN
        NBines = NBines + 1
        LimiteSuperiorBines(NBines) = Limite
        !WRITE(*,*) "Bin: ", NBines, Limite
     ELSE
        NBines = NBines + 1
        LimiteSuperiorBines(NBines) = Lmax
        !WRITE(*,*) "Bin: ", NBines, Limite
        Leer = 0
     END IF

     GOTO 200
100  Leer = 0
200  CONTINUE
  END DO

  CLOSE(10)

  ALLOCATE(ErrorDlTT(2:lmax))
  ALLOCATE(ErrorDlEE(2:lmax))
  ALLOCATE(ErrorDlBB(2:lmax))
  ALLOCATE(ErrorDlTE(2:lmax))
  ALLOCATE(ErrorDlTB(2:lmax))
  ALLOCATE(ErrorDlEB(2:lmax))

  npixfull = 12D0*nside**2
  fskyT = (1.0D0*npixT)/npixfull
  fskyP = (1.0D0*npixP)/npixfull

  IF(TipoCentrado==1) THEN

     WRITE(*,*)
     WRITE(*,*) " Intensity noise per pixel:    ", RuidoTT
     WRITE(*,*) " Polarization noise per pixel: ", RuidoQQ

     DO l=2,lMax
        ErrorDlTT(l) =  ErrorTeorico(1, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             & DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
        ErrorDlEE(l) =  ErrorTeorico(2, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             & DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
        ErrorDlBB(l) =  ErrorTeorico(3, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             & DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
        ErrorDlTE(l) =  ErrorTeorico(4, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             & DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
        ErrorDlTB(l) =  ErrorTeorico(5, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             & DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
        ErrorDlEB(l) =  ErrorTeorico(6, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB,&
             DlTE, DlTB, DlEB, BeamPWDl, RuidoTT, RuidoQQ)
     END DO
  ELSE
     ErrorDlTT(:) = 1D0
     ErrorDlEE(:) = 1D0
     ErrorDlBB(:) = 1D0
     ErrorDlTE(:) = 1D0
     ErrorDlTB(:) = 1D0
     ErrorDlEB(:) = 1D0
  END IF

  ALLOCATE(DLBin(6, NBines))
  ALLOCATE(LBin(6, NBines))
  ALLOCATE(SumaPesos(6))
  ALLOCATE(SumasDl(6))
  ALLOCATE(SumasL(6))
  ALLOCATE(PesosCompactado(6, 2:lmax))

  li = 2
  DO i=1,NBines

     ls = LimiteSuperiorBines(i)

     SumaPesos(:) = 0d0
     SumasDl(:) = 0d0
     SumasL(:) = 0d0

     DO l=li,ls

        temp = 1.0d0/ErrorDlTT(l)**2
        SumaPesos(1) = SumaPesos(1) + temp
        SumasDl(1) = SumasDl(1) + DlTT(l) * temp
        SumasL(1) = SumasL(1) + l * temp
        PesosCompactado(1,l) = temp

        temp = 1.0d0/ErrorDlEE(l)**2
        SumaPesos(2) = SumaPesos(2) + temp
        SumasDl(2) = SumasDl(2) + DlEE(l) * temp
        SumasL(2) = SumasL(2) + l * temp
        PesosCompactado(2,l) = temp

        temp = 1.0d0/ErrorDlBB(l)**2
        SumaPesos(3) = SumaPesos(3) + temp
        SumasDl(3) = SumasDl(3) + DlBB(l) * temp
        SumasL(3) = SumasL(3) + l * temp
        PesosCompactado(3,l) = temp

        temp = 1.0d0/ErrorDlTE(l)**2
        SumaPesos(4) = SumaPesos(4) + temp
        SumasDl(4) = SumasDl(4) + DlTE(l) * temp
        SumasL(4) = SumasL(4) + l * temp
        PesosCompactado(4,l) = temp

        temp = 1.0d0/ErrorDlTB(l)**2
        SumaPesos(5) = SumaPesos(5) + temp
        SumasDl(5) = SumasDl(5) + DlTB(l) * temp
        SumasL(5) = SumasL(5) + l * temp
        PesosCompactado(5,l) = temp

        temp = 1.0d0/ErrorDlEB(l)**2
        SumaPesos(6) = SumaPesos(6) + temp
        SumasDl(6) = SumasDl(6) + DlEB(l) * temp
        SumasL(6) = SumasL(6) + l * temp
        PesosCompactado(6,l) = temp

     END DO

     SumasDl = SumasDl/SumaPesos
     SumasL = SumasL/SumaPesos

     DLBin(:,i) = SumasDl(:)
     LBin(:,i)  = SumasL(:)

     PesosCompactado(1,li:ls) = PesosCompactado(1,li:ls)/SumaPesos(1)
     PesosCompactado(2,li:ls) = PesosCompactado(2,li:ls)/SumaPesos(2)
     PesosCompactado(3,li:ls) = PesosCompactado(3,li:ls)/SumaPesos(3)
     PesosCompactado(4,li:ls) = PesosCompactado(4,li:ls)/SumaPesos(4)
     PesosCompactado(5,li:ls) = PesosCompactado(5,li:ls)/SumaPesos(5)
     PesosCompactado(6,li:ls) = PesosCompactado(6,li:ls)/SumaPesos(6)

     !-----------------
     li = ls + 1

  END DO


  DEALLOCATE(ErrorDlTT)
  DEALLOCATE(ErrorDlEE)
  DEALLOCATE(ErrorDlBB)
  DEALLOCATE(ErrorDlTE)
  DEALLOCATE(ErrorDlTB)
  DEALLOCATE(ErrorDlEB)

  OPEN(Unit=25, File=TRIM(Lugar)//"/Positions.dat", Action="write")
  DO Caso=1,6
     DO i=1,NBines
        WRITE(25,*) LBin(Caso, i)
     END DO
  END DO
  CLOSE(25)

  OPEN(Unit=25, File=TRIM(Lugar)//"/BinnedFiducial.dat", Action="write")
  DO Caso=1,6
     DO i=1,NBines
        WRITE(25,*) DLBin(Caso, i)
     END DO
  END DO
  CLOSE(25)

  !Calcula los factores para pasar de Dl a DLBin
  ALLOCATE(Factores(6,2:lmax))
  Factores = 0

  IF(TipoCompactado==1) THEN
     li = 2
     DO i=1,NBines

        ls = LimiteSuperiorBines(i)
        DO l=li,ls

           IF(DlBin(1,i)/=0) THEN
              Factores(1,l) = DlTT(l)/DlBin(1,i)
           ELSE
              Factores(1,l) = 1D0
           END IF

           IF(DlBin(2,i)/=0) THEN
              Factores(2,l) = DlEE(l)/DlBin(2,i)
           ELSE
              Factores(2,l) = 1D0
           END IF

           IF(DlBin(3,i)/=0) THEN
              Factores(3,l) = DlBB(l)/DlBin(3,i)
           ELSE
              Factores(3,l) = 1D0
           END IF

           IF(DlBin(4,i)/=0) THEN
              Factores(4,l) = DlTE(l)/DlBin(4,i)
           ELSE
              Factores(4,l) = 1D0
           END IF

           IF(DlBin(5,i)/=0) THEN
              Factores(5,l) = DlTB(l)/DlBin(5,i)
           ELSE
              Factores(5,l) = 1D0
           END IF

           IF(DlBin(6,i)/=0) THEN
              Factores(6,l) = DlEB(l)/DlBin(6,i)
           ELSE
              Factores(6,l) = 1D0
           END IF
        END DO
        li = ls + 1
     END DO
  ELSE
     li = 2
     DO i=1,NBines
        ls = LimiteSuperiorBines(i)
        DO l=li,ls
           Factores(:,l) = 1D0
        END DO
        li = ls + 1
     END DO

  END IF

  !Construye la matriz de los Factores
  ALLOCATE(MatrizR(6*(lmax-1),6*NBines))

  MatrizR = 0.0D0
  f = 1
  c = 1
  DO Caso = 1,6
     li = 2
     DO i=1,NBines
        ls = LimiteSuperiorBines(i)
        DO l=li,ls
           MatrizR(c,f) =  Factores(Caso,l)
           !write(*,*) "R: ", c, f, MatrizR(c,f)
           c = c + 1
        END DO
        li = ls+1
        f = f + 1
     END DO
  END DO

  DEALLOCATE(Factores)

  !Carga la matriz de Fisher
  dim = 6*(lmax-1)
  ALLOCATE(MF(dim*dim))
  WRITE(*,*)
  WRITE(*,*) " Loading Fisher matrix"
  OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrix.dat", Action="read")
  READ (25, *)  MF
  CLOSE(25)

  !Compacta la matriz de Fisher
  ALLOCATE(MatrizProducto(dim,6*NBines))  

  CALL dsymm("L", "U", dim, 6*NBines, 1.0D0, MF, dim, MatrizR, dim, 0.0D0, MatrizProducto, dim)

  ALLOCATE(MFBines(6*Nbines, 6*Nbines))
  CALL dgemm("T","N", 6*NBines, 6*NBines, dim, 1.0D0, MatrizR, dim, MatrizProducto, dim, 0.0D0, MFBines, 6*NBines)

  DEALLOCATE(MatrizProducto)

  OPEN(Unit=25, File=TRIM(Lugar)//"/FisherMatrixBinnedDl.dat", Action="write")
  !DO i=1,6*NBines
  !WRITE(25,*)  MFBines(i,:)
  !END DO
  WRITE(25,*) MFBines
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
  READ(25,*) Yl
  CLOSE(25)

  IF(QuitarSesgo==1) THEN
     DO i=1, NMapas
        !READ (25, *)  Yl(:,i)
        Yl(:,i) = Yl(:,i) - Bl
     END DO
     !ELSE
     !   DO i=1, NMapas
     !      READ (25, *)  Yl(:,i)
     !   END DO
  END IF



  !Invierte la matriz de Fisher compactada
  WRITE(*,*) " Inverting the compacted Fisher matrix"

  info = 0
  dim = 6*NBines
  CALL dpotrf("U", dim, MFBines, dim, info)
  WRITE(*,*) " Cholesky factorization: ", info

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
  CALL dgemm("T","N", dim, NMapas, 6*(lmax-1), 1.0D0, MatrizR, 6*(lmax-1), Yl, 6*(lmax-1), 0.0D0, YlBines, dim)

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

  REAL(kind=8) FUNCTION ErrorTeorico(Caso, l, fskyT, fskyP, npixfull, DlTT, DlEE, DlBB, DlTE, DlTB, DlEB,&
       & BeamPWDl, RuidoTT, RuidoQQ)

    IMPLICIT NONE

    INTEGER, INTENT(in) :: Caso, l
    REAL(kind=8), INTENT(in) :: RuidoTT, RuidoQQ, fskyT, fskyP, npixfull
    REAL(kind=8), DIMENSION(0:lmax), INTENT(in) :: DlTT, DlEE, DlBB, DlTE, DlTB, DlEB
    REAL(kind=8), DIMENSION(3, 0:lmax), INTENT(in) :: BeamPWDl

    !Internas
    REAL(kind=8) error
    REAL(kind=8), PARAMETER :: DosPi = 6.2831853071796D0


    SELECT CASE (Caso)
    CASE (1)
       error = DlTT(l) + (2D0*DosPi/npixfull)*RuidoTT**2/BeamPWDl(1,l)**2
       error = error**2
       error = 2.0D0/(2.0*l+1.0) * 1.0D0/fskyT * error
       ErrorTeorico = SQRT(error)
    CASE (2)
       error = DlEE(l) + (2D0*DosPi/npixfull)*RuidoQQ**2/BeamPWDl(2,l)**2
       error = error**2
       error = 2.0D0/(2.0*l+1.0) * 1.0D0/fskyP * error
       ErrorTeorico = SQRT(error)
    CASE (3)
       error = DlBB(l) + (2D0*DosPi/npixfull) * RuidoQQ**2/BeamPWDl(2,l)**2
       error = error**2
       error = 2.0D0/(2.0*l+1.0) * 1.0D0/fskyP * error
       ErrorTeorico = SQRT(error)
    CASE (4)
       error = DlTE(l)**2 + (DlTT(l) + (2D0*DosPi/npixfull) * RuidoTT**2/BeamPWDl(1,l)**2) *&
            & (DlEE(l) + (2D0*DosPi/npixfull) * RuidoQQ**2/BeamPWDl(2,l)**2)
       error = 1.0D0/(2.0*l+1.0) * 1.0D0/((fskyT+fskyP)/2.0D0) * error
       ErrorTeorico = SQRT(error)
    CASE (5)
       error = DlTB(l)**2 + (DlTT(l) + (2D0*DosPi/npixfull) * RuidoTT**2/BeamPWDl(1,l)**2) *&
            & (DlBB(l) + (2D0*DosPi/npixfull) * RuidoQQ**2/BeamPWDl(2,l)**2)
       error = 1.0D0/(2.0*l+1.0) * 1.0D0/((fskyT+fskyP)/2.0D0) * error
       ErrorTeorico = SQRT(error)
    CASE (6)
       error = DlEB(l)**2 + (DlEE(l) + (2D0*DosPi/npixfull) * RuidoQQ**2/BeamPWDl(2,l)**2) *&
            & (DlBB(l) + (2D0*DosPi/npixfull) * RuidoQQ**2/BeamPWDl(2,l)**2)
       error = 1.0D0/(2.0*l+1.0) * 1.0D0/fskyP * error
       ErrorTeorico = SQRT(error)
    END SELECT

    RETURN

  END FUNCTION ErrorTeorico

END SUBROUTINE CalculaBineado
