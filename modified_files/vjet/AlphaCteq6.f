      double precision function alphal(qmu2,iset)
      implicit real * 8 (a-h,o-z)
      INTEGER ISETP,IORD,leanasorder
      double precision lam4,lam,b0,b1,mu
      double precision lam5nlo,lam4nlo
      double precision aaa,nfl,ll,aln,bb0,bb
      data pi/3.14159265358979323846d0/
      data dpi/6.28318530717958647692d0/
      data lam4nlo/0.326d0/
      data lam5nlo/0.226d0/
      data lam4lo/0.192d0/
      data lam5lo/0.146d0/
      data lam3ffn/0.395d0/
      data lam4ffn/0.309d0/
c Cteq6 LO lam's differ from CTEQ5:
      data lamC64lo/0.215d0/
      data lamC65lo/0.165d0/
c I've shifted mrst lam4's so that they match at 4.5 GeV
c with the lam5's, giving the correct alpha_s(MZ).
      data lam4mrst/0.318/
      data lam4mrst4/0.244/
      data lam4mrst5/0.405/
      data lam4mrst6/0.3227/
      data lam4mrst7/0.3115/
      data lam5mrst/0.220/
      data lam5mrst4/0.164/
      data lam5mrst5/0.289/
      data lam5mrst6/0.2237/
      data lam5mrst7/0.215/
      data mb/4.5d0/
        
       mu=sqrt(qmu2) 
       nfl=5.0d0
       if(iset.eq.4) then
       
       lam=lamC65lo
       lam4=lamC64lo
       else 
       
       lam=lam5nlo
       lam4=lam4nlo
       endif
        b0=(11.0d0-2.0d0/3.0d0*nfl)
        b1=102.0d0-38.0d0/3.0d0*nfl 
        ll=lam
        bb0=b0
        bb1=b1
  
       aln=log(qmu2/ll/ll)
c       write(*,*) aln, b1, b0, ll, lam, mu
       if(aln.le.0.0d0) then   
       write(*,*) ' Value of Q given was  mu '
       write(*,*) ' Value of iset was  iset '
       write(*,*) ' You have reached the Landau singularity'
       write(*,*) ' All hope is lost!'
       endif
       
       aaa=4.0d0*pi/bb0/aln
       if(iset.eq.4) then
c       write(*,*)  'Usando alpha_s a Iset=4  ', aaa
       alphal=aaa
       return 
       else 
       alphal=aaa*(1.0d0-bb1*log(aln)/aln/bb0/bb0)
c       write(*,*)  'Usando alpha_s a un loop, a Iset=1 ', alphal
       return
       endif
       
       end 
     
      
