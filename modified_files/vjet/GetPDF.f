      subroutine getapdfs(x,q,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot,order,ih)

      IMPLICIT NONE
      INTEGER iset,
     &     alphaSorder,alphaSnfmax,order,ih
      
      DOUBLE PRECISION x,q,upv,dnv,usea,dsea,str,sbar,
     &     chm,cbar,bot,bbar,glu,phot,
     &     
     &     
     &     distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ
c      double precision ALPHAS
     
       CHARACTER prefix*50,flavours(-5:5)*10
c       character prefix1*55
      COMMON/mstwCommon/distance,tolerance,
     &     mCharm,mBottom,alphaSQ0,alphaSMZ,alphaSorder,alphaSnfmax
       DATA flavours /"bbar","cbar","sbar","ubar","dbar","glu",
     &     "dn","up","str","chm","bot"/

      if(order.eq.83)then
      prefix = "Grids/mstw2008nnlo" ! prefix for the grid files
      elseif(order.eq.82)then
      prefix = "Grids/mstw2008nlo" ! prefix for the grid files
      elseif(order.eq.84)then
      prefix = "Grids/mstw2008lo" ! prefix for the grid files 
      else
      write(*,*) 'Mala eleccion orden PDFs!!!!!'  
      endif
      iset = 0

      
      CALL GetAllPDFs(prefix,iset,x,q,upv,dnv,usea,dsea,str,sbar,
     &        chm,cbar,bot,bbar,glu,phot)

      !    write(*,*) 'desde fortran', upv
     
     
      END
C----------------------------------------------------------------------




