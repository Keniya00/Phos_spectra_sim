        PROGRAM FRANCK_CONDON
                IMPLICIT NONE
                INTEGER,PARAMETER::N=96   !N= NO. OF NORMAL MODES 
                INTEGER,PARAMETER::N1=192  !N1=2*N
                INTEGER::I,J,II,NN,INFO,KK,L
                INTEGER,DIMENSION(N1)::IPIV
                REAL*8::EV_AU,CM_AU,SEC_AU,AU_HZ,J_AU,KB,BETA,SOCME, &
                DEL_E,DIPOLE,SL
                REAL*8::TI,TF,T,H,TEMP,ETA,Z,RATE,T1
                REAL*8,DIMENSION(N,N)::JMAT,WSI,WSF,ASF,BSF,P
                REAL*8,DIMENSION(N,1)::DMAT,WS1,WT1,WS,WT,DP_DER,DP_DER1
                REAL*8,DIMENSION(1,N)::DP_DER_TR,DP_DER_TR1
                COMPLEX*16::DET,SQDET,FUNC,NFUNC,EXT,EXT1,EXT2,DAMP, &
                        X,TR_KINV,HT,HT_FUNC,CO_FUNC,XT,XT1,XT2
                COMPLEX*16,DIMENSION(N1)::WORK
                COMPLEX*16,DIMENSION(N,N)::JTC,JC,ASFC,BSFC,ASID,BSID, &
                        BSIN,ASI,BSI,ASI_J,A2,A,BSI_J,B2,B,C
                COMPLEX*16,DIMENSION(N,N)::BSQ,B_INV,BINVA,A_BINVA, &
                ASQ,KMAT,KMAT_INV,ASFASI,ASIFK,PC,PPC,ASIFKZ
                COMPLEX*16,DIMENSION(N,1)::DC,CD,JTCD
                COMPLEX*16,DIMENSION(1,N)::DTC,DTCC
                COMPLEX*16,DIMENSION(N1,N1)::K
                COMPLEX*16,DIMENSION(N1,1)::F,KINV_F
                COMPLEX*16,DIMENSION(1,N1)::FT,FT_KINV,KINV_FTRANS
                COMPLEX*16,DIMENSION(1,1)::DTCD,FT_KINV_F,KFKF

                EV_AU=0.036749844D0
                CM_AU=4.55634D0*(10.0D0**(-6.0D0))
                SEC_AU=0.4134D0*(10.0D0**17.0D0)
                AU_HZ=6.579D0*(10.0D0**15.0D0)
                J_AU=2.2937126583579D0*(10.0D0**17.0D0)
                DEL_E=2.285714286d0 !ENERGY GAP BETWEEN 2 STATES
                DEL_E=DEL_E*EV_AU
                SL=137.03599D0
                DIPOLE=0.00003d0  !TRANSITION MOMENT
                ETA=100.0D0 !DAMPING PARAMETER
                ETA=ETA*CM_AU
                KB=1.380649D0*(10.0D0**(-23.0D0))
                TEMP=78.0D0 !TEMPERATURE
                BETA=1.0D0/(KB*TEMP*J_AU)

                OPEN(1,FILE='WS0_DPPZ_THF.DAT') !FINAL FREQUENCY FILE
                OPEN(2,FILE='WT1_DPPZ_THF.DAT') !INITIAL FREQUENCY FILE 
                OPEN(3,FILE='displacement_t1_s0.dat') !DISPLACEMENT
                  !VECTOR
                OPEN(4,FILE='duschinsky_t1_s0.dat') !DUSCHINSKY ROT
                   !MATRIX
                OPEN(90,FILE='REALFUNC_EMS_DPPZ_THF.DAT')
                OPEN(91,FILE='IMAG_FUNC_EMS_DPPZ_THF.DAT')
                OPEN(92,FILE='REAL_FUNC_EMS_DPPZ_THF.DAT')
                OPEN(93,FILE='TIME.DAT')
                READ(1,*)(WS1(I,1),I=1,N)
                READ(2,*)(WT1(I,1),I=1,N)
                READ(3,*)(DMAT(I,1),I=1,N)
                READ(4,*)((JMAT(I,J),J=1,N),I=1,N)

                NN=100000  !NO. OF POINTS
                TI=-10E-12  !INITIAL TIME
                TF=10E-12   !FINAL TIME
                TI=TI*SEC_AU
                TF=TF*SEC_AU
                H=(TF-TI)/FLOAT(NN)

                DO I=1,N
                WT(I,1)=WT1(I,1)*CM_AU
                WS(I,1)=WS1(I,1)*CM_AU
                ENDDO
                DO I=1,N
                DO J=1,N
                 IF(I.EQ.J)THEN
                         WSI(I,J)=WT(I,1)
                         WSF(I,J)=WS(I,1)
                 ELSE
                      WSI(I,J)=0.0D0
                      WSF(I,J)=0.0D0
              ENDIF
              ENDDO
              ENDDO   
                
                DO I=1,N
                DC(I,1)=CMPLX((DMAT(I,1)),0.0D0)
               ! DC(I,1)=CMPLX(0.0D0,0.0D0)
                ENDDO
                DTC=TRANSPOSE(DC)
                DO I=1,N
                DO J=1,N
                JC(I,J)=CMPLX((JMAT(I,J)),0.0D0)
                ENDDO
                ENDDO
                JTC=TRANSPOSE(JC)

                NFUNC=(0.0D0,0.0D0)
              
                DO I=1,N
                DO J=1,N
                IF(I.EQ.J)THEN
                P(I,J)=2.0d0*SINH((WSI(I,J)*BETA)/2.0D0)
                PC(I,J)=CMPLX((P(I,I)),0.0D0)
                   ELSE
                           P(I,J)=0.0D0
                           PC(I,J)=(0.0D0,0.0D0)
                   ENDIF
                           ENDDO
                ENDDO
                DP_DER_TR=TRANSPOSE(DP_DER)
                !!TIME LOOP STARTS HERE=============================================================================
                DO II=1,NN+1
                T=TI+(II-1)*H
                 IF(T.EQ.0.0D0)THEN
                        T=(10.0D0**(-30.0D0))*SEC_AU
                 ENDIF
                DO I=1,N
                 DO J=1,N
                   IF(I.EQ.J)THEN
                        ASF(I,J)=WSF(I,J)/SIN(WSF(I,J)*T)
                        BSF(I,J)=WSF(I,J)/TAN(WSF(I,J)*T)
                        ASID(I,J)=CMPLX((SIN(-WSI(I,J)*T)* &
         COSH(-WSI(I,J)*BETA)),(COS(-WSI(I,J)*T)*SINH(-WSI(I,J)*BETA)))
                        ASI(I,J)=WSI(I,J)/ASID(I,J)
         BSID(I,J)=CMPLX(TAN(-WSI(I,J)*T),(TANH(-WSI(I,J)*BETA)))/ &
          (1.0D0-((0.0D0,1.0D0)*TAN(-WSI(I,J)*T)*TANH(-WSI(I,J)*BETA)))
                         BSI(I,J)=WSI(I,J)/BSID(I,J)
                   ELSE
                        ASF(I,J)=0.0D0
                        BSF(I,J)=0.0D0
                        ASID(I,J)=(0.0D0,0.0D0)
                        ASI(I,J)=(0.0D0,0.0D0)
                        BSID(I,J)=(0.0D0,0.0D0)
                        BSI(I,J)=(0.0D0,0.0D0)
                   ENDIF
                  ENDDO
                ENDDO
                DO I=1,N
                  DO J=1,N
                   ASFC(I,J)=CMPLX((ASF(I,J)),0.0D0)
                   BSFC(I,J)=CMPLX((BSF(I,J)),0.0D0)
                  ENDDO
                ENDDO
              !!GENERATION OF A,B,C=======================================================================================
             
               CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),ASI, &
                N,JC,N,(0.0D0,0.0D0),ASI_J,N)
               CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),JTC,N,ASI_J,N,(0.0D0,0.0D0),A2,N)
               DO I=1,N
                 DO J=1,N
                   A(I,J)=ASFC(I,J)+A2(I,J)     !GENERATION OF A
                 ENDDO
               ENDDO
               CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),BSI,N,JC, &
                N,(0.0D0,0.0D0),BSI_J,N)
               CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),JTC,N,BSI_J,N,(0.0D0,0.0D0),B2,N)
               DO I=1,N
                  DO J=1,N
                  B(I,J)=BSFC(I,J)+B2(I,J)         !GENERATION OF B
                  B_INV(I,J)=B(I,J)
                ENDDO
               ENDDO
              !!REPRESENTATION OF 2N*2N MATRIX IN N*N FORM===============================================
                CALL ZGETRF(N,N,B_INV,N,IPIV,INFO)
                CALL ZGETRI(N,B_INV,N,IPIV,WORK,N,INFO)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),B,N,B,N,(0.0D0,0.0D0),BSQ,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),B_INV,N,A,N,(0.0D0,0.0D0),BINVA,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),A,N,BINVA,N,(0.0D0,0.0D0),A_BINVA,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),B,N,A_BINVA,N,(0.0D0,0.0D0),ASQ,N)
                DO I=1,N
                 DO J=1,N
                 KMAT(I,J)=BSQ(I,J)-ASQ(I,J)
                 KMAT_INV(I,J)=KMAT(I,J)
                 ENDDO
                ENDDO
                CALL ZGETRF(N,N,KMAT_INV,N,IPIV,INFO)
                CALL ZGETRI(N,KMAT_INV,N,IPIV,WORK,N,INFO)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),ASFC,N,ASI,N,(0.0D0,0.0D0),ASFASI,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),ASFASI,N, &
                KMAT_INV,N,(0.0D0,0.0D0),ASIFK,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0), &
                PC,N,PC,N,(0.0D0,0.0D0),PPC,N)
                CALL ZGEMM('N','N',N,N,N,(1.0D0,0.0D0),ASIFK,N,PPC,N,(0.0D0,0.0D0),ASIFKZ,N)
                !!CALCULATION OF DETERMINANT==================================================================
                CALL ZGETRF(N,N,ASIFKZ,N,IPIV,INFO)
                DET=(1.0D0,0.0D0)
                DO I=1,N
                DET=DET*ASIFKZ(I,I)
                ENDDO
                SQDET=SQRT(ABS(DET))*EXP((0.0D0,1.0D0)*((ATAN &
                (AIMAG(DET)/REAL(DET)))/2.0D0))      !SQUAREROOT OF DETERMINANT
               ! WRITE(*,*)T,(SQDET)
                DO I=1,N
                 DO J=1,N
                 C(I,J)=BSI(I,J)-ASI(I,J)               !GENERATION OF C 
                 ENDDO
                ENDDO
            
                !!CALCULATION OF EXPONENTIAL TERMS==============================================================
                !!GENERATION OF EXP(iDTCD)

                CALL ZGEMM('N','N',1,N,N,(1.0D0,0.0D0),DTC,1,C,N, &
                (0.0D0,0.0D0),DTCC,1)
                CALL ZGEMM('N','N',1,1,N,(1.0D0,0.0D0),DTCC,1,DC,N,(0.0D0,0.0D0),DTCD,1)
                EXT=EXP((1.0D0)*((0.0D0,1.0D0)*DTCD(1,1)))
                XT=(1.0D0)*((0.0D0,1.0D0)*DTCD(1,1))
                !!GENERATION OF EXP(-(i/2)*FT_KINV_F)

                CALL ZGEMM('N','N',N,1,N,(1.0D0,0.0D0),C,N,DC,N,(0.0D0,0.0D0),CD,N)
                CALL ZGEMM('N','N',N,1,N,(1.0D0,0.0D0),JTC,N,CD, &
                N,(0.0D0,0.0D0),JTCD,N)
                DO I=1,N
                F(I,1)=JTCD(I,1)
                F(I+N,1)=JTCD(I,1)
                ENDDO
                FT=TRANSPOSE(F)
                DO I=1,N
                 DO J=1,N
                 K(I,J)=B(I,J)
                 K(I,J+N)=-A(I,J)
                 K(I+N,J)=-A(I,J)
                 K(I+N,J+N)=B(I,J)
                 ENDDO
                ENDDO
                CALL ZGETRF(N1,N1,K,N1,IPIV,INFO)
                CALL ZGETRI(N1,K,N1,IPIV,WORK,N1,INFO)
                CALL ZGEMM('N1','N1',1,N1,N1,(1.0D0,0.0D0),FT,1,K,N1,(0.0D0,0.0D0),FT_KINV,1)
                CALL ZGEMM('N1','N1',1,1,N1,(1.0D0,0.0D0),FT_KINV, &
                1,F,N1,(0.0D0,0.0D0),FT_KINV_F,1)
                X=(FT_KINV_F(1,1)/2.0D0)
               EXT1=EXP((-1.0D0)*((0.0D0,1.0D0)*(FT_KINV_F(1,1)/2.0D0)))
               XT1=(-1.0D0)*((0.0D0,1.0D0)*(FT_KINV_F(1,1)/2.0D0))
                
                !!GENERATION OF EXP(-iET) & DAMPING FACTOR

                EXT2=EXP((1.0D0)*((0.0D0,1.0D0)*DEL_E*T))
                XT2=(1.0D0)*((0.0D0,1.0D0)*DEL_E*T)
                DAMP=EXP((-1.0D0)*((1.0D0,0.0D0)*ETA*ABS(T)))
                T1=ABS(T)
                FUNC=SQDET*EXP(XT+XT1+XT2)*EXP(-T1*ETA)*(DIPOLE**2.0D0)
                FUNC=FUNC*((2.0d0*(DEL_E**3.0d0))/(3.0d0*3.14d0*(SL**3.0D0)))
                WRITE(90,*)REAL(FUNC)
                WRITE(91,*)T,AIMAG(FUNC)
                WRITE(92,*)T,REAL(FUNC)
                WRITE(93,*)T

               ENDDO    !TIME LOOP ENDS HERE
                END PROGRAM

                
