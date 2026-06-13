      subroutine setPDF2gnnlo
      implicit none 
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
      
      
      iset2g = 82
      ih12g  = 1 ! 1 FOR PROTON , -1 FOR ANTIPROTON 
      ih22g  = 1 ! 1 FOR PROTON , -1 FOR ANTIPROTON

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
      INTEGER ISETP,IORD,leanasorder
      double precision leanalphamz,alphal
      data pi/3.14159d0/
      data dpi/6.28318530717958647692d0/
c      common/xlam/xlam
c      common/isetp/isetp
      common/pdf/ih1,ih2,ndns1,ndns2
      
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
     
      common/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      common/leanasorder/leanasorder
      common/leanalphamz/leanalphamz
      isetp=iset2g
c       isetp=ndns1
      qq2=q2
      qq=sqrt(q2)
c      write(*,*) 'q',qq 
      alphass=0d0
      if(isetp.eq.40) then
        iord=0
        alambda=0.220
        alphass=alphal(qq2,isett)/4/PI !mrst lo
       elseif(isetp.eq.43) then !alpha NLO Cteq6
        iord=1
c        qqm=67.5*67.5
        alphass=alphal(qq2,iord)/4.0d0/PI
c        alphass=0.118d0/2.0d0/dpi 
c        write (*,*) alphal(qqm,iord)
        elseif(isetp.eq.46) then !alpha LO Cteq6
         iord=4
       alphass=alphal(qq2,iord)/4.0d0/PI
       elseif(isetp.eq.50) then !mrst2002 nlo
        iord=1
        alambda=0.3342
        alphass=alpha2002(qq2,alambda,iord)/4/PI !mrst2002 nlo


       elseif((isetp.eq.82).or.(isetp.eq.83).or.(isetp.eq.84))then  !MSTW2008
 
         if((isetp.eq.84))then
         alphaSMZLN = 0.13938999999999999d0
         leanasorder=0
         elseif((isetp.eq.82))then
          alphaSMZLN = 0.12018000000000000d0
          leanasorder=1
         else
          alphaSMZLN = 0.11706999999999999d0
          leanasorder=2
          endif
          
          
         CALL INITALPHAS(leanasorder,1.D0,91.1876d0, 
     &  alphaSMZLN,1.4d0,4.75d0,1.d10)

        alphass=ALPHAS(qq)/2.d0/dpi !MSTW2008

       

        else 
           write(6,*) "wrong alpha_s"
           stop
       endif
      return
      end






C--------------------------------------------------
C- STRUCTURE FUNCTION MAIN PROGRAM
C--------------------------------------------------
      SUBROUTINE MLMPDF(NDNS,IH,Q2,X,FX,NF)
C USE ONE OF THE FREE NDNS NUMBERS TO INCLUDE A NEW SET      
       implicit double precision (b-h,o-z)
       implicit real*8 (a) 
c       double precision afext(11)   
       double precision FX(-NF:NF)
        INTEGER iset,
     &     alphaSorder,alphaSnfmax
      
      
      integer iset2g,ih12g,ih22g
      common/parpdfs2gamma/iset2g,ih12g,ih22g
      
      CHARACTER prefix*50,flavours(-5:5)*10

      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
      DATA flavours /"bbar","cbar","sbar","ubar","dbar","glu",
     &     "dn","up","str","chm","bot"/

c      INTEGER iset
c      real * 8 dx,dq,duv,ddv,ddel,dudb,dsb,dgl,chr,bot,ub,db,dq1
c      real*8 xmin,xmax,qsqmin,qsqmax,qsq,ctq4pdf,u,d
c      real*8 ixmin,ixmax,iqsqmin,iqsqmax,x,q2
       data xmin,xmax,qsqmin,qsqmax/1.d-5,1.d0,0.4d0,1.d6/
       data ini/0/
       common/pdfsext/afext(11)
       save/pdfsext/
       iset=iset2g
c      dx=dble(x)
c      dq=dble(q2)
      
      dx=x
      dq=q2
      
      dq1=dq**0.5      
C.....Fix to prevent undefined math operations for x=1.
C.....Assumes that all structure functions vanish for x=1.
      if(1-x.eq.0) then
         do j=-nf,nf
            fx(j) = 0
         enddo
         return
      endif
C.....Initialization, all to 0
      do J=-nf,nf
         fx(j) = 0
      enddo
      chr=0d0
      bot=0d0
      
      if ((iset.gt.42).and.(iset.lt.47)) then
C.....new CTEQ6 PDFs
         mode=iset-42
         Call SetCtq6(mode)
         dq1=dq**0.5
         dgl=Ctq6Pdf (0, dx, dq1)*dx
         ub=Ctq6Pdf (-1, dx, dq1)*dx
         db=Ctq6Pdf (-2, dx, dq1)*dx
         DSB=Ctq6Pdf (-3, dx, dq1)*dx
         str=Ctq6Pdf (3, dx, dq1)*dx
         chr=Ctq6Pdf (4, dx, dq1)*dx
         chrb=Ctq6Pdf (-4, dx, dq1)*dx
         botb=Ctq6Pdf (-5, dx, dq1)*dx  
         bot=Ctq6Pdf (5, dx, dq1)*dx
c        duv=Ctq6Pdf (1, dx, dq1)*dx 
c        ddv=Ctq6Pdf (2, dx, dq1)*dx 
         duv=Ctq6Pdf (1, dx, dq1)*dx - ub
         ddv=Ctq6Pdf (2, dx, dq1)*dx - db
     
           elseif (iset.eq.50) then 
C NEW MRST2002 NLO
         mode=1
         dq1=dq**0.5 
         call mrst2002(dx,dq1,mode,DUV,DDV,UB,DB,DSB,CHR,BOT,DGL)



      elseif (iset.eq.82) then 
      prefix = "Grids/mstw2008nlo" ! prefix for the grid files
C NEW MSTW2008 NLO
      isettt = 0
C--   First the traditional MRST-like interface
C--   (but note the "sbar", "cbar", "bbar" and "phot").
      CALL GetAllPDFs(prefix,isettt,dx,dq1,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot)
      dgl=glu
      ub=usea
      duv=upv
      ddv=dnv
      db=dsea
      chr=chm  
       
      elseif (iset.eq.83) then 
      prefix = "Grids/mstw2008nnlo" ! prefix for the grid files
      
      
C NEW MSTW2008 NNLO
      isettt = 0
      CALL GetAllPDFs(prefix,isettt,dx,dq1,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot)

      dgl=glu
      ub=usea
      duv=upv
      ddv=dnv
      db=dsea
      chr=chm    
       
      elseif (iset.eq.84) then 
      prefix = "Grids/mstw2008lo" ! prefix for the grid files
      
      
C NEW MSTW2008 LO
      isettt = 0
      CALL GetAllPDFs(prefix,isettt,dx,dq1,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot)

      dgl=glu
      ub=usea
      duv=upv
      ddv=dnv
      db=dsea
      chr=chm    
       else
         write(*,*)'No such distribution'
         stop
      endif
    
c          write(*,*) nf = Esta en 5
c.....for protons
       if (ih.eq.1) then
         fx(0)=(dgl)
         fx(1)=(duv+ub)
         fx(2)=(ddv+db)
       if (iset.eq.43) then 
            fx(3)=(str)
         elseif (iset.eq.82) then  
             fx(3)=(str)
         elseif (iset.eq.83) then  
             fx(3)=(str)
         elseif (iset.eq.84) then  
             fx(3)=(str)
         else 
         fx(3)=(str)
       endif
         if(nf.ge.4) fx(4)=(chr)
         if(nf.ge.5) fx(5)=(bot)
         if(nf.eq.6) fx(6)=0
         fx(-1)=(ub)
         fx(-2)=(db)
         if (iset.eq.43) then 
            fx(-3)=(dsb)
         elseif (iset.eq.82) then  
             fx(-3)=(sbar)
         elseif (iset.eq.83) then  
             fx(-3)=(sbar)
         elseif (iset.eq.84) then  
             fx(-3)=(sbar)
         else 
         fx(-3)=(dsb)
         endif
         if(nf.ge.4) then
               if (iset.eq.43) then 
            fx(-4)=(chrb)
         elseif (iset.eq.82) then  
             fx(-4)=(cbar)
         elseif (iset.eq.83) then  
            fx(-4)=(cbar)
         elseif (iset.eq.84) then  
            fx(-4)=(cbar)  
            else
            fx(-4)=(chrb)
            endif
            endif
         if(nf.ge.5) then
               if (iset.eq.43) then 
            fx(-5)=(botb)
         elseif (iset.eq.82) then  
             fx(-5)=(bbar)
         elseif (iset.eq.83) then  
            fx(-5)=(bbar)
         elseif (iset.eq.84) then  
            fx(-5)=(bbar) 
            else
            fx(-5)=(botb)
            endif
            endif 
       if(nf.eq.6) fx(-6)=0
       elseif (ih.eq.-1) then
         fx(0)=(dgl)
         fx(-1)=(duv+ub)
         fx(-2)=(ddv+db)
       if (iset.eq.43) then 
            fx(-3)=(str)
         elseif (iset.eq.82) then  
             fx(-3)=(str)
         elseif (iset.eq.83) then  
             fx(-3)=(str)
         elseif (iset.eq.84) then  
             fx(-3)=(str)
         else 
         fx(-3)=(str)
       endif
         if(nf.ge.4) fx(-4)=(chr)
         if(nf.ge.5) fx(-5)=(bot)
         if(nf.eq.6) fx(-6)=0
         fx(1)=(ub)
         fx(2)=(db)
         if (iset.eq.43) then 
            fx(3)=(dsb)
         elseif (iset.eq.82) then  
             fx(3)=(sbar)
         elseif (iset.eq.83) then  
             fx(3)=(sbar)
         elseif (iset.eq.84) then  
             fx(3)=(sbar)
         else 
         fx(3)=(dsb)
         endif
         if(nf.ge.4) then
               if (iset.eq.43) then 
            fx(4)=(chrb)
         elseif (iset.eq.82) then  
             fx(4)=(cbar)
         elseif (iset.eq.83) then  
            fx(4)=(cbar)
         elseif (iset.eq.84) then  
            fx(4)=(cbar)  
            else
            fx(4)=(chrb)
            endif
            endif
         if(nf.ge.5) then
               if (iset.eq.43) then 
            fx(5)=(botb)
         elseif (iset.eq.82) then  
             fx(5)=(bbar)
         elseif (iset.eq.83) then  
            fx(5)=(bbar)
         elseif (iset.eq.84) then  
            fx(5)=(bbar) 
            else
            fx(5)=(botb)
            endif
            endif 
         if(nf.eq.6) fx(6)=0
      endif 
      do i=-nf,nf
         afext(i+6)=fx(i)/dx
         fx(i)=fx(i)/X
      enddo
      
      
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

 
