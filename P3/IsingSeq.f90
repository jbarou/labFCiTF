!C23456789012345678901234567890123456789012345678901234567890123456789012
!     ******************************************************************
!     *                                                                *
!     *     MC1-2.f                                                    *
!     *     IMPROVED VERSION                                           *
!     *     MONTE CARLOS SIMULATION OF THE 2D ISING MODEL              *
!     *     METROPOLIS ALGORITHM                                       *
!     *     EV 05-10-2019                                              *
!     *     EV 03-12-2010 (Revisat)                                    *
!     *     EV 03-11-2022 (Revisat)                                    *
!     *     JB 27-11-2022 (Revisat)                                    *
!     *     JB 27-02-2026 (Revisat i extracte de conf. )               *
!     *                                                                *
!     ******************************************************************
!
!     compile with
!
!     > gfortran -O3 -o IsingSeq IsingSeq.f90
!
!     * Variable declaration ********************************************
!
      IMPLICIT NONE
!
!     Input data
!
      integer (kind=4) :: L
!!      PARAMETER (L=48)
      double precision TEMP, t4, t8
      integer (kind=4) :: SEED
      integer (kind=4) :: MCTOT
      double precision genrand_real2
      logical accept
!
!     Spin matrix
!
      INTEGER (kind=2), dimension (:,:), allocatable :: S
      integer (kind=2) :: Sold
!
!     Vector per les condicions periodiques de contorn
!
      INTEGER (kind=4), dimension (:), allocatable :: PBC
!
!     Variables per calculs
!
      integer (kind=4) :: N, MAG
      double precision ENE
      INTEGER DE
      integer (kind=4) :: SUMA
!
!     Lectura de temps
!      
      REAL*4 TIME1,TIME2
      CHARACTER*64 DATE, cmd
!
!     Index interns
!
      integer (kind=4) :: I,J,IJ, IMC, IPAS
!
!     Funcions externes d'ininicialitzacio
!
      double precision ENERG
      INTEGER (kind=4) :: MAGNE
!
!     arxiu output
!
      CHARACTER*128 NOM
      CHARACTER*32 sL, sTEMP, sMCTOT, sSEED
!
!     dades per linia de comandes
!     > IsingSeq  L TEMP MCTOT

      print*,"use: $ IsingSeq  L TEMP MCTOT SEED"

      CALL GETARG(1 , sL)
      CALL GETARG(2 , sTEMP)
      CALL GETARG(3 , sMCTOT)
      CALL GETARG(4 , sSEED)
      
      READ (sL,*) L
      READ (sTEMP,*) TEMP
      READ (sMCTOT,*) MCTOT
      READ (sSEED,*) SEED
      
      allocate(S(1:L,1:L))
      allocate(PBC(0:L+1))
      
      Write(NOM,'("resultsSeq/SIM-L",I0,"-TEMP",F5.3)') L,TEMP
!
!     Variables calculades
! 
!     Nombre total de particules
! 
      N=L*L
!
!     Cond. periodiques de contorn
!
      PBC(0)=L
      PBC(L+1)=1
      do i=1,L
         PBC(i)=i
      enddo
!
!     Probabilitats de transicio no trivial
!
      t4 = exp(-4.0d0 / TEMP)
      t8 = exp(-8.0d0 / TEMP)
!
!     Obre output files
!
      open(UNIT=12,FILE=trim(nom)//"_EM.seq", status="unknown", action="write")
      open(UNIT=13,FILE=trim(nom)//"_map.conf", status="unknown", action="write")
      open(UNIT=14,FILE="logIsingSeq.out",  status="unknown", position='append', action="write")
 
      CALL CPU_TIME(TIME1)      
      CALL FDATE(DATE) 
      call get_command(cmd)
      

      WRITE (14,*) ""
      WRITE (14,*) "========================="
      WRITE (14,*) trim(DATE)
      WRITE (14,*) trim(cmd)
      WRITE (14,*) ""
      WRITE (14,*) "L = ", L, "TEMP = ", TEMP, "MCTOT= " ,MCTOT , "SEED = ", SEED
      WRITE (14,*) "========================="      
      WRITE (14,*) ""


      WRITE (*,*) ""
      WRITE (*,*) "========================="
      WRITE (*,*) trim(DATE)
      WRITE (*,*) trim(cmd)
      WRITE (*,*) ""
      WRITE (*,*) "L = ", L, "TEMP = ", TEMP, "MCTOT= " ,MCTOT , "SEED = ", SEED
      WRITE (*,*) "========================="      
      WRITE (*,*) ""

!
!     Initializatio del generador amb valor SEED
!
      CALL init_genrand(1235+SEED*2)
!
!     generacio de matriu inicial
!

      DO J=1,L
         DO I=1,L
            IF (genrand_real2().GT.0.5D0) THEN
                S(I,J) = 1
            ELSE
                S(I,J)=-1
            ENDIF
         ENDDO
         write(13,*) S(:,J)
      ENDDO
      write(13,*) ''
      write(13,*) ''
!
!     Initial energy
!
      ENE=ENERG(S,L,PBC)
      MAG=MAGNE(S,L)
      IMC=0

      WRITE(*,*) 'MC=',IMC, ' ENERGIA =', ENE, 'MAGNE =',MAG      
      WRITE(14,*) 'MC=',IMC, ' ENERGIA =', ENE, 'MAGNE =',MAG
      WRITE(12,*) IMC, ENE, MAG, N
!
!     Main Monte Carlo loop
!
      DO IMC=1,MCTOT
!
!        Loop de passes MC (N intents de canvi)
!
         DO IPAS = 1,N
!
!           Escull spin aleatori
            IJ=int(genrand_real2()*N)
            I  = IJ / L + 1
            J  = IJ - (I-1)*L + 1  

            Sold = S(i,j)
            !
            !  Suma dels veins de l'spi S(I,J)
            !      
            SUMA=S(PBC(I+1),J)+S(I,PBC(J+1))+S(I,PBC(J-1))+S(PBC(I-1),J) 
            !
            !  Canvi d'energia associat al gir de l'spi S(I,J)
            !
            DE=2*Sold*SUMA
            !
            !  Decisio Metropolis
            !
      
            select case (de)
                    case (: 0)                ! <= 0
                       accept = .true.
                    case (4)
                       accept = (genrand_real2() < t4)
                    case (8)
                       accept = (genrand_real2() < t8)
                    case default
                       accept = .false.        ! de = {+4,+8} únics positius possibles
            end select
            if (accept) then
               S(i,j) = -Sold
               ENE    = ENE + de
               MAG    = MAG - 2*sold
            end if
            !

!
!
!          Final de passa MC
!
         ENDDO
!
!            
         WRITE(12,*) IMC, ENE, MAG, N
!
!        Bolca la configuracio cada 1000 pasos
!
      if(mod(IMC,1000).eq.0) then
            WRITE(*,*) 'MC=',IMC, ' ENERGIA =', ENE, 'MAGNE =',MAG
            WRITE(14,*) 'MC=',IMC, ' ENERGIA =', ENE, 'MAGNE =',MAG
            do J=1,L
                  write(13,*) S(:,J)
            enddo
            write(13,*) ''
            write(13,*) ''
      endif
      ENDDO





      CLOSE(12)
!
!     Escriu conf. final       
!
      do J=1,L
            write(13,*) S(:,J)
      enddo
      
      CLOSE(13)

      CALL CPU_TIME(TIME2)

      CALL FDATE(DATE)
      
      WRITE (14,*) ""
      WRITE (14,*) DATE
      WRITE (14,*) ""
      WRITE (14,*) 'CPUTIME = ', TIME2-TIME1
      WRITE (14,*) ""
      WRITE (14,*) "========================="      
      WRITE (14,*) ""


      WRITE (*,*) ""
      WRITE (*,*) DATE
      WRITE (*,*) ""
      WRITE (*,*) 'CPUTIME = ', TIME2-TIME1
      WRITE (*,*) ""
      WRITE (*,*) "========================="      
      WRITE (*,*) ""
      
      STOP
       
      END
!
!
!     ******************************************************************
!     *                    FUNCTION MAGNE                              *
!     ******************************************************************
!
      integer (kind=4) FUNCTION MAGNE(S,L)
      integer (kind=2) :: S(1:L,1:L)
      integer (kind=4) :: I,J,L
      integer (kind=4) :: MAG
      MAG=0
      DO I=1,L
         DO J=1,L
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
      double precision FUNCTION ENERG(S,L,CPC)
      integer (kind=2) :: S(1:L,1:L)
      integer (kind=4) :: I,J,L
      integer (kind=4) :: CPC(0:L+1)
      double precision ENE
      ENE=0.0D0
      DO I=1,L
         DO J=1,L
            ENE=ENE-S(I,J)*S(CPC(I+1),J)-S(I,J)*S(I,CPC(J+1))
         ENDDO
      ENDDO
     
      ENERG=ENE

      RETURN
      END
      
      
!
!     ******************************************************************
!     *               SUBROUTINE WRITE CONF. MATRIX                    *
!     ******************************************************************
!
!      subroutine PrintConf(S)      
      

      
      
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




     
