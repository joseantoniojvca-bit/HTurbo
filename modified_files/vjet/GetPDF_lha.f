      subroutine getapdfs(x,q,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot,order,ih)

      IMPLICIT NONE
      INTEGER iset2g,ih12g,ih22g,
     &     order,nf,ih
      
      DOUBLE PRECISION x,q,upv,dnv,usea,dsea,str,sbar,
     &     chm,cbar,bot,bbar,glu,phot,
     &      fx(-5:5),q2
c      double precision ALPHAS
     
       CHARACTER flavours(-5:5)*10
c       character prefix*50 !,prefix1*55
       DATA flavours /"bbar","cbar","sbar","ubar","dbar","glu",
     &     "dn","up","str","chm","bot"/
c      common/parpdfs2gamma/iset2g,ih12g,ih22g

       
       NF=5
       Q2=q*q
c       ih=1

      call MLMPDF(1,IH,Q2,X,FX,NF) 

       if (ih.eq.1) then
         glu=fx(0)*x
         str=fx(3)*x
         chm=fx(4)*x
         bot=fx(5)*x
         usea=fx(-2)*x
         dsea=fx(-1)*x
         upv=-usea+fx(2)*x
         dnv=-dsea+fx(1)*x
         sbar=fx(-3)*x
         cbar=fx(-4)*x
         bbar=fx(-5)*x
         phot=0d0
       elseif (ih.eq.-1) then
         glu=fx(0)*x
         str=fx(-3)*x
         chm=fx(-4)*x
         bot=fx(-5)*x
         usea=fx(2)*x
         dsea=fx(1)*x
         upv=-usea+fx(-2)*x
         dnv=-dsea+fx(-1)*x
         sbar=fx(3)*x
         cbar=fx(4)*x
         bbar=fx(5)*x
         phot=0d0
       endif


!      prefix = "Grids/mstw2008nnlo" ! prefix for the grid files
!      write(*,*) "----------------"
!      write(*,*)"M",x,q,glu,str,chm,bot,cbar,bbar,usea,dsea,upv,dnv,phot
!
!      CALL GetAllPDFs(prefix,0,x,q,upv,dnv,usea,dsea,str,sbar,
!     &        chm,cbar,bot,bbar,glu,phot)
!
!      write(*,*)"G",x,q,glu,str,chm,bot,cbar,bbar,usea,dsea,upv,dnv,phot
!      write(*,*) "----------------"
     
     
      END
C----------------------------------------------------------------------




