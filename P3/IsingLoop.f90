!23456789012345678901234567890123456789012345678901234567890123456789012
!
!     SIMULACIO MODEL D'ISING 2D
!     ALGORISME DE METROPOLIS
!     EV 05-10-2019
!     EV 03-12-2020 (Revisat)
!     EV 03-11-2022 (Revisat)
!     JB 20-05-2024 (nova versió) 
!     JB 17-06-2024 (Revisat)
!     JB 1-03-2025  (optimitzat i parallel OMP)
!    
!***********************************************************************
!
!     Compilar amb
!
!     > gfortran -O3 -o MC3.exe MC3.f
program IsingLoop
!
!     Declaracio de variables 
!
      IMPLICIT NONE
!
!     Dades d'entrada
!
      INTEGER(kind=4) :: L

      DOUBLE PRECISION :: TMin, Tmax, DT, TEMP
      INTEGER (kind=4) :: NCONF
      INTEGER (kind=4) :: MCTOT,MCINI, MCD
      DOUBLE PRECISION :: genrand_real2
!
!     Spin matrix
!
      INTEGER (kind=2), dimension (:,:), allocatable :: S
!
!     Vector for periodic boundary conditions
!
      INTEGER (kind=4), dimension (:), allocatable :: PBC
!
!

 
!     Variables calculades 
!
      INTEGER (kind=2) :: Sold
      INTEGER (kind=4) :: NTEMP
      INTEGER (kind=4) :: ICONF
      INTEGER (kind=4) :: N, IJ
      INTEGER (kind=8) :: ENE, MAG
      INTEGER (kind=4) :: SUMA
      INTEGER (kind=4) :: DE
      REAL  (kind=4) :: TIME1,TIME2
      CHARACTER*30 DATE
!
!     Indexos interns 
!
      INTEGER :: I,J,IPAS,IMC, ITEMP
!
!     Matriu probabilitats de transicio
!
      DOUBLE PRECISION  :: w4, w8
      logical :: accept
!
!     Variables per fer els promitjos d' interes
!  
! --- Estadístiques Welford globals per temperatura
      INTEGER (kind=8) ::  COUNT_AM, COUNT_E
      Double Precision ::  MEAN_E,  M2_E,  M2_AM,  MEAN_AM

! --- (si emprem OpenMP) estadístiques locals per fil (es faran private)
      INTEGER (kind=8) ::  TCOUNT
      Double Precision ::  TMEAN_E, TM2_E, TM2_AM, TMEAN_AM

! --- auxiliar per a updates
      Double Precision :: DELTA
      Double Precision :: VARE, VARM

!
!     Funcions externes
!
 
       INTEGER (kind=8) :: MAGNE, ENERG
!
!     Name of the output file
!
      CHARACTER*128 NOM
      CHARACTER*32 sL
      

!    llegeix input parameters de IsingLoop.nml
     namelist /inParams/ L, TMin, Tmax, DT, NCONF, MCINI, MCD, MCTOT
     open(10, file="IsingLoop.nml", status="old")
     read(10, nml=inParams)
     close(10)


      print*,"use: $ IsingLoop  L "


      print*, "# ---------------------------------------------------"
!     si hi ha arguments, afaga el primer com a L
      if (iargc()>0) then
         print*, "# L OVERRIDE BY INLINE ARGUMENT"
         CALL GETARG(1 , sL)
         READ (sL,*) L
      end if
      
     allocate(PBC(0:L+1))



      
      Write(NOM,'("resultsLoop/SIM-L",I0,"-MCTOT",I0)') L,MCTOT
      OPEN(UNIT=13,FILE=trim(NOM)//".res")



!
      CALL init_genrand(123456)

!
!     Control de temps de CPU inicial
!
      CALL CPU_TIME(TIME1)
      
! 
!     Nombre total de particules
!

       N=L*L


      print*, "# ---------------------------------------------------"
      print '("#  L = ", I0 ," ( N = LxL = ",I0," )")', L , N 
      print '("#  Temp: ", F5.3, " ... ", F5.3, " ...  ", F5.3)',  Tmin, DT, TMax
      print '("#  Nconf = ", I0 )', Nconf
      print '("#  Measure starts at MC: ", I0, " , every: ", I0, " , to: ", I0)',   MCINI, MCD, MCTOT
      print*, "# ---------------------------------------------------"

      
!
!     Condicions periodiques
!
      PBC(0)=L
      PBC(L+1)=1
      do i=1,L
         PBC(i)=i
      enddo

!
!     Comencem el bucle de temperatures
!

      Ntemp = ceiling((TMax-TMin)/dT)
      DO itemp=1,Ntemp,1

            TEMP=TMin+(itemp-1)*Dt
            PRINT '("TEMP = ", F5.3, " ( ",F5.1 ," % )")', TEMP,  itemp*100.0d0/ntemp
      !
      !     vector with the transitions probabilities
      !

            w4 = exp(-4.0d0 / TEMP)
            w8 = exp(-8.0d0 / TEMP)

      !
      !     Posem a zero els promitjos d'interes
      !


            COUNT_AM    = 0_8
            COUNT_E    = 0_8
            MEAN_E  = 0.0D0
            M2_E    = 0.0D0
            M2_AM    = 0.0D0
            MEAN_AM = 0.0D0






            ! --- Bucle de configuracions (paral·lel si -fopenmp; serial si no)
            !$omp parallel default(none) &
            !$omp& shared(L,PBC,TEMP,w4,w8,N,MCINI,MCD,MCTOT,NCONF,ITEMP, &
            !$omp&        MEAN_E,M2_E,COUNT_E,COUNT_AM,M2_AM, MEAN_AM) &
            !$omp& private(ICONF,I,J,IPAS,IMC,S,ENE,MAG,SUMA,DE,accept, &
            !$omp&         IJ,Sold,TCOUNT,TMEAN_E,TM2_E,TM2_AM,TMEAN_AM,DELTA)

            ! -- Inicialitza estadístiques locals de fil
                  TCOUNT  = 0_8
                  TMEAN_E = 0.0D0
                  TM2_E   = 0.0D0
                  TM2_AM   = 0.0D0
                  TMEAN_AM= 0.0D0

            !$omp do schedule(static)
                  DO ICONF=1,NCONF,1

            !        Llavor independent per simulació/temperatura (determinista)
                     call init_genrand(12345 + 1000*ITEMP + ICONF*2)

            !        ---- array de spins privat per fil ----
                     allocate(S(1:L,1:L))

            !        Inicialitza spins a l'atzar per fil
                     DO J=1,L
                       DO I=1,L
                         IF (genrand_real2().GT.0.5D0) THEN
                            S(I,J) = 1
                         ELSE
                            S(I,J) = -1
                         ENDIF
                       ENDDO
                     ENDDO

            !        Energia i magnetització inicials
                     ENE = ENERG(S,L,PBC)
                     MAG = MAGNE(S,L)

            !        Bucle de passes de MC
                     DO IMC=1,MCTOT

            !           N intents
                        DO IPAS = 1,N

            !             Selecció aleatoria d'un lloc (1 RNG)
                           IJ = int(genrand_real2()*N)
                           I  = IJ / L + 1
                           J  = IJ - (I-1)*L + 1

                           Sold = S(I,J)

            !              Suma veïns (PBC per vector)
                           SUMA = S(PBC(I+1),J) + S(I,PBC(J+1)) + S(I,PBC(J-1)) + S(PBC(I-1),J)

            !              Canvi d'energia
                           DE = 2*Sold*SUMA

            !              Decisió Metropolis (branching compacte)
                           select case (DE)
                             case (:0)                   ! DE <= 0
                               accept = .true.
                             case (4)
                               accept = (genrand_real2() .LT. w4)
                             case (8)
                               accept = (genrand_real2() .LT. w8)
                             case default
                               accept = .false.          ! DE {+4,+8} únics positius possibles
                           end select

                           if (accept) then
                              S(I,J) = -Sold
                              ENE    = ENE + DE
                              MAG    = MAG - 2*Sold
                           end if

                        ENDDO  ! IPAS

            !           Mostreig (cada MCD a partir de MCINI)
                        IF ( (IMC.GT.MCINI) .AND. (MOD(IMC,MCD).EQ.1) ) THEN

                          TCOUNT  = TCOUNT + 1_8

                          ! --- Welford per ENERGIA ---
                          DELTA   = DBLE(ENE) - TMEAN_E
                          TMEAN_E = TMEAN_E + DELTA / DBLE(TCOUNT)
                          TM2_E   = TM2_E   + DELTA * (DBLE(ENE) - TMEAN_E)

                          ! --- Welford per AM = |MAG| ---
                          DELTA    = DABS(DBLE(MAG)) - TMEAN_AM
                          TMEAN_AM = TMEAN_AM + DELTA / DBLE(TCOUNT)
                          TM2_AM   = TM2_AM   + DELTA * (DABS(DBLE(MAG)) - TMEAN_AM)

                        END IF
                        

                     ENDDO  ! IMC

                     deallocate(S)

                  ENDDO    ! ICONF
            !$omp end do

            ! --- Fusió Welford (per fil -> globals). Per seguretat, si es creuen 2 arribares alhora.
            !$omp critical
                  call WELFORD_MERGE(COUNT_E, MEAN_E, M2_E, TCOUNT, TMEAN_E, TM2_E)
                  call WELFORD_MERGE(COUNT_AM, MEAN_AM, M2_AM, TCOUNT, TMEAN_AM, TM2_AM)
            !$omp end critical

            !$omp end parallel


                                    
      ! --- Variàncies (per espín; evitem mesures buides (per si de cas))
            if (  COUNT_E    .GT. 0_8) then
              VARE = M2_E / DBLE(  COUNT_E   )
            else
              VARE = 0.0D0
            end if
            
            if (  COUNT_AM    .GT. 0_8) then
              VARM = M2_AM / DBLE(  COUNT_AM   )
            else
              VARM = 0.0D0
            end if



      ! --- Força no-negativitat per rodoneig
      !     VARE = MAX(0.0D0, VARE)
      !     VARM = MAX(0.0D0, VARM)

      
      ! --- comprova consistencia contadors
            if (COUNT_E /= COUNT_AM) then
              write(*,*) "WARNING: COUNT_E != COUNT_AM at TEMP=", TEMP, COUNT_E, COUNT_AM
            end if
      ! --- Sortida d'arxiu

            WRITE(13,*) N, TEMP,   COUNT_AM   , MEAN_E, VARE, MEAN_AM, VARM
                  
      !
      ENDDO
!
!     Tanquem arxiu
!
      CLOSE(13)
!
!     Control de temps de CPU final
!
      CALL CPU_TIME(TIME2)
!
!     Data i hora final
!
      CALL FDATE(DATE)

      WRITE (*,*) DATE
      WRITE (*,*) 'CPUTIME = ', TIME2-TIME1
      WRITE (*,*) MCINI, MCTOT

      STOP
       
      END program
!
!
!     ******************************************************************
!     *                    FUNCTION MAGNE                              *
!     ******************************************************************
!
      integer (kind=8)  FUNCTION MAGNE(S,L)
      INTEGER*2 S(1:L,1:L)
      INTEGER*4 I,J,L
      integer (kind=8) :: MAG
      MAG=0
      DO J=1,L
         DO I=1,L
            MAG=MAG+S(I,J)
         ENDDO
      ENDDO
     
      MAGNE=MAG

      RETURN
      END
!
!     ******************************************************************
!     *                    FUNCTION ENERG                              *
!     ******************************************************************
!
      integer (kind=8) FUNCTION ENERG(S,L,PBC)
      INTEGER*2 S(1:L,1:L)
      INTEGER*4 I,J,L
      INTEGER*4 PBC(0:L+1)
      integer (kind=8) ENE
      ENE=0
      DO J=1,L
         DO I=1,L
            ENE=ENE-S(I,J)*S(PBC(I+1),J)-S(I,J)*S(I,PBC(J+1))
         ENDDO
      ENDDO
     
      ENERG=ENE

      RETURN
      END


! *********************************************************************
! *                  SUBROUTINE WELFORD_MERGE                         *
! *  Fusiona (n_a,mean_a,M2_a) amb (n_b,mean_b,M2_b) -> in-place      *
! *********************************************************************
      SUBROUTINE WELFORD_MERGE(N_A, MEAN_A, M2_A, N_B, MEAN_B, M2_B)
      INTEGER*8 N_A, N_B, N
      REAL*8    MEAN_A, M2_A, MEAN_B, M2_B, DELTA

      IF (N_B .LE. 0_8) RETURN

      IF (N_A .EQ. 0_8) THEN
         MEAN_A = MEAN_B
         M2_A   = M2_B
         N_A    = N_B
         RETURN
      ENDIF

      DELTA = MEAN_B - MEAN_A
      N     = N_A + N_B

      MEAN_A = MEAN_A + DELTA * DBLE(N_B) / DBLE(N)
      M2_A   = M2_A + M2_B + (DELTA*DELTA) * DBLE(N_A) * DBLE(N_B) / DBLE(N)
      N_A    = N

      RETURN
      END

!
!*****************************************************************************
!
!  A C-program for MT19937, with initialization improved 2002/1/26.
!  Coded by Takuji Nishimura and Makoto Matsumoto.
!
!  Before using, initialize the state by using init_genrand(seed)  
!  or init_by_array(init_key, key_length).
!
!  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
!  All rights reserved.                          
!  Copyright (C) 2005, Mutsuo Saito,
!  All rights reserved.                          
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions
!  are met:
!
!    1. Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer.
!
!    2. Redistributions in binary form must reproduce the above copyright
!       notice, this list of conditions and the following disclaimer in the
!       documentation and/or other materials provided with the distribution.
!
!    3. The names of its contributors may not be used to endorse or promote 
!       products derived from this software without specific prior written 
!       permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
!  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!  Any feedback is very welcome.
!  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!  email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
!
!-----------------------------------------------------------------------
!  FORTRAN77 translation by Tsuyoshi TADA. (2005/12/19)
!
!     ---------- initialize routines ----------
!  subroutine init_genrand(seed): initialize with a seed
!  subroutine init_by_array(init_key,key_length): initialize by an array
!
!     ---------- generate functions ----------
!  integer function genrand_int32(): signed 32-bit integer
!  integer function genrand_int31(): unsigned 31-bit integer
!  double precision function genrand_real1(): [0,1] with 32-bit resolution
!  double precision function genrand_real2(): [0,1) with 32-bit resolution
!  double precision function genrand_real3(): (0,1) with 32-bit resolution
!  double precision function genrand_res53(): (0,1) with 53-bit resolution
!
!  This program uses the following non-standard intrinsics.
!    ishft(i,n): If n>0, shifts bits in i by n positions to left.
!                If n<0, shifts bits in i by n positions to right.
!    iand (i,j): Performs logical AND on corresponding bits of i and j.
!    ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!    ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!-----------------------------------------------------------------------
!     initialize mt(0:N-1) with a seed
!-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
!
!$omp threadprivate(/mt_state1/,/mt_state2/,/mt_mask1/)

      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
!
      return
      end
!-----------------------------------------------------------------------
!     initialize by an array with array-length
!     init_key is the array for initializing keys
!     key_length is its length
!-----------------------------------------------------------------------
      subroutine init_by_array(init_key,key_length)
      integer init_key(0:*)
      integer key_length
      integer N
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      parameter (N=624)
      integer i,j,k
      integer mt(0:N-1)
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
!
!$omp threadprivate(/mt_state2/,/mt_mask1/,/mt_mask2/)

      call init_genrand(19650218)
      i=1
      j=0
      do 100 k=max(N,key_length),1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1664525)+init_key(j)+j
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        j=j+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
        if(j.ge.key_length)then
          j=0
        endif
  100 continue
      do 200 k=N-1,1,-1
        mt(i)=ieor(mt(i),ieor(mt(i-1),ishft(mt(i-1),-30))*1566083941)-i
        mt(i)=iand(mt(i),ALLBIT_MASK)
        i=i+1
        if(i.ge.N)then
          mt(0)=mt(N-1)
          i=1
        endif
  200 continue
      mt(0)=TOPBIT_MASK
!
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0xffffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
!
!$omp threadprivate(/mt_state1/,/mt_state2/,/mt_mask3/,/mt_mag01/)

      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
!
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
!
      y=mt(mti)
      mti=mti+1
!
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
!
      genrand_int32=y
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,0x7fffffff]-interval
!-----------------------------------------------------------------------
      function genrand_int31()
      integer genrand_int31
      integer genrand_int32
      genrand_int31=int(ishft(genrand_int32(),-1))
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1]-real-interval
!-----------------------------------------------------------------------
      function genrand_real1()
      double precision genrand_real1,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real1=r/4294967295.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on (0,1)-real-interval
!-----------------------------------------------------------------------
      function genrand_real3()
      double precision genrand_real3,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real3=(r+0.5d0)/4294967296.d0
      return
      end
!-----------------------------------------------------------------------
!     generates a random number on [0,1) with 53-bit resolution
!-----------------------------------------------------------------------
      function genrand_res53()
      double precision genrand_res53
      integer genrand_int32
      double precision a,b
      a=dble(ishft(genrand_int32(),-5))
      b=dble(ishft(genrand_int32(),-6))
      if(a.lt.0.d0)a=a+2.d0**32
      if(b.lt.0.d0)b=b+2.d0**32
      genrand_res53=(a*67108864.d0+b)/9007199254740992.d0
      return
      end
!-----------------------------------------------------------------------
!     initialize large number (over 32-bit constant number)
!-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
      !$omp threadprivate(/mt_mask1/,/mt_mask2/,/mt_mask3/,/mt_mag01/)
!C    TOPBIT_MASK = Z'80000000'
!C    ALLBIT_MASK = Z'ffffffff'
!C    UPPER_MASK  = Z'80000000'
!C    LOWER_MASK  = Z'7fffffff'
!C    MATRIX_A    = Z'9908b0df'
!C    T1_MASK     = Z'9d2c5680'
!C    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
