      subroutine setPDF2gnnlo
      implicit none 
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
!
      double precision alphasPDF,zmass,leanalphamz
      logical validPDF
      character*30 oldPDFname,PDFname
      integer i,PDFmember
      logical lhapdfs,PDFerrors
      common/lhapdfs/lhapdfs

      common/leanalphamz/leanalphamz
      integer maxPDFsets,currentPDF
c--- 40 is my choice for the maximum number of PDF error sets
c--- 50 is my choice for the maximum number of dipoles
c      double precision PDFxsec(0:40),PDFxsec_nd(0:40,0:50),
c     . PDFwgt(0:40)
c      common/PDFerrors/PDFerrors,maxPDFsets,PDFxsec,PDFwgt  
      character*7 pdlabel
      
      
      iset2g = 0
      ih12g  = 1 ! 1 FOR PROTON , -1 FOR ANTIPROTON 
      ih22g  = 1 ! 1 FOR PROTON , -1 FOR ANTIPROTON

      zmass=91.1876d0

c      PDFname='MSTW2008lo68cl.LHgrid'
c       PDFname='cteq6m.LHpdf'
      PDFname='NNPDF31_nnlo_as_0118'
      PDFmember=0



!      include 'masses.f'
!      include 'lhapdf.f'
!      include 'PDFerrors.f'
!      include 'pdlabel.f'
!      double precision amz,alphasPDF

      lhapdfs=.true.
      
c      if (newinput .eqv. .false.) then
c        open(unit=21,file='lhapdf.DAT',status='old',err=999)
c        call checkversion(21,'lhapdf.DAT')
c        read(21,*) PDFname
c        read(21,*) PDFmember            
c        close(21)
c      endif
      
      oldPDFname=PDFname
      validPDF=.false.
      i=0
   20 continue
      i=i+1    
      if ((oldPDFname(i:i) .eq. '.') .or.
     .    (oldPDFname(i:i) .eq. ' ') .or.
     .    (oldPDFname(i:i) .eq. '[')) then
        validPDF=.true.
        if (oldPDFname(i:i+6) .eq. '.LHgrid') then        
          PDFname=oldPDFname(1:i-1)//'.LHgrid'
        else
          PDFname=oldPDFname(1:i-1)//'.LHpdf'
        endif
      endif  
      if ((i .lt. 20) .and. (validPDF .eqv. .false.)) goto 20
      
c$$$      if (validPDF .eqv. .false.) then
c$$$        write(6,*) 'Problem with PDFname'
c$$$        write(6,*)
c$$$        stop
c$$$      endif
      

      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*) 'C                                                  C'
      write(6,*) 'C            2gNNLO now calling LHAPDF             C'
      write(6,*) 'C                                                  C'
      write(6,98) 'PDFname',PDFname(1:20)
      write(6,99) 'PDFmember',PDFmember
      write(6,*) 'C                                                  C'
      write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write(6,*)

      call InitPDFset('/home/lcieri/LHC/share/LHAPDF/'//PDFname)
      
      if (PDFmember .lt. 0) then
        PDFerrors=.true.
        call numberPDF(maxPDFsets)
        if (maxPDFsets .gt. 50) then
          write(6,*) 'ERROR: Max. number of error sets is 50!'
          stop
        endif
        write(6,*)
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'        
        write(6,*) 'C        Calculating errors using      C'
        write(6,*) 'C        ',maxPDFsets,' sets of error PDFs        C'
        write(6,*) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
        call InitPDF(0)
        leanalphamz=alphasPDF(zmass)
        currentPDF=0
      else  
        call InitPDF(PDFmember)
        leanalphamz=alphasPDF(zmass)
      endif

c--- rename pdlabel to get sensible output name
      pdlabel=PDFname(1:7)

      return
 
   98 format(' C            ',a7,' ',a20,'          C')
   99 format(' C                ',a10,i3,'                     C')

!  999 write(6,*) 'Error reading lhapdf.DAT'
!      call flush(6)
!      stop

      end
 

c          84: MSTW 2008       LO
C          82: MSTW 2008      NLO
C          83: MSTW 2008      NNLO
C          50: MRST 2002      NLO
C          40: MRST 2002       LO 
C          43: CTEQ 6         NLO
c          46: CTEQ 6          LO 





C  SI SE USA CTEQ SOLO DISPONIBLE VERSIONES 6

C IF A NEW SET IS INCLUDED IN "MLMPDF", ADD HERE THE CORRESPONDING ALPHA_S CALL
C ACTUALLY ALPHA_S/4/PI GOES HERE
      double precision function alphass(q2,isett)
      implicit real * 8 (a-h,o-z)
      INTEGER ISETP
      double precision leanalphamz,order
      data pi/3.14159d0/
      data dpi/6.28318530717958647692d0/
c      common/xlam/xlam
c      common/isetp/isetp
      common/pdf/ih1,ih2,ndns1,ndns2
      
      integer iset2g,ih12g,ih22g,iifborn,iinnlo,iifnnlo
      common/parpdfs2gamma/iset2g,ih12g,ih22g
     
      common/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      common/leanalphamz/leanalphamz
      common/calcula/iifborn,iinnlo,iifnnlo
      isetp=iset2g
c       isetp=ndns1
      qq2=q2
      qq=sqrt(q2)
c      write(*,*) 'q',qq 
       
      if(iifnnlo.eq.1) then
       order=2
      elseif(iifnlo.eq.1) then      
       order=1
      else
       order=0
      endif

      alphaSMZLN = leanalphamz

         CALL INITALPHAS(order,1.D0,91.1876d0, 
     &  alphaSMZLN,1.4d0,4.75d0,1.d10)

c        alphass=ALPHAS(qq)/2.d0/dpi !MSTW2008
        alphass=alphasPDF(qq)/2.d0/dpi ! Alphas LHA
       
      return
      end






C--------------------------------------------------
C- STRUCTURE FUNCTION MAIN PROGRAM
C--------------------------------------------------
*****************
* LHAPDF version*
*****************
      subroutine MLMPDF(NDNS,IH,Q2,X,FX,NF) 
! fdist(ih,x,xmu,fx)
      implicit none
      double precision fx(-NF:NF),x,xmu,Q2,fPDF(-6:6)
      integer Iprtn,ih,nf,ndns
c---  ih1=+1 proton 
c---  ih1=-1 pbar 
      xmu=dsqrt(Q2)

C---set to zero if x out of range
      if (x .ge. 1d0) then
          do Iprtn=-NF,NF
             fx(Iprtn)=0d0
          enddo
          return
      endif
 
      call evolvePDF(x,xmu,fPDF)
c Checkear que de lo mismo numericamente
c Checkeado formalmente
c      if (ih.eq.1) then
        do Iprtn=-NF,NF
          fx(+Iprtn)=fPDF(+Iprtn)/x
        enddo
c      elseif(ih.eq.-1) then
c        do Iprtn=-NF,NF
c          fx(+Iprtn)=fPDF(-Iprtn)/x
c        enddo
c      endif
      
      return
      end

c-------------------------------------------------------------------------
      function zgmu2(xmur2,isetp)
c sets the mass scales and returns the strong coupling squared
c-------------------------------------------------------------------------
      implicit real * 8 (a-h,o-z)
      parameter (pi=3.14159265358979312D0)

        as = 4d0*pi*alphass(xmur2,isetp)
      zgmu2 = 4.d0*pi*as
      zg = sqrt(zgmu2)
      return
      end

 
