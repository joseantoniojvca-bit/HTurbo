             subroutine Setname
             
             character * 25 prefn
             common/vname/prefn
             
             prefn = 'Histo' ! Name of the top file (output file)
               
                  
             
             return
             end


C INIJET INICIALIZA LOS HISTOGRAMAS
C Y DEFINE NOMBRE Y CARACTERISTICAS DE LOS MISMOS
  
             subroutine inijet
             implicit real * 8 (a-h,o-z)
c Reset histograms            
             call inihist
c Here book the needed histograms, using:
c call bookplot(hist. number, title, num. of bins, xmin, xmax, ....
c              ... xtitle, ytitle, yscale, nrep)
c  where yscale can be either   'LOG'
c                         or    'LIN'
c and for nrep=1  a new figure is created
c whereas for nrep=2 the hist. is included in the figure defined before. 
c The maximum number of allowed histos is 50 
c The maximum number of allowed bins is 100 
c 
c
c Notice that some plots only differ by the bin size

            call bookmbobjects
          
c GM: it would be nice to define the histograms to be filled
c     with the variables and the titles, ranges, bins, ... 
c     in an external card file
!        call bookplot(1,'ds/dMgg',20,100.d0,160.d0,'Mgg'
!     &          ,'fb/','LIN',1)!
!
!        call bookplot(2,'ds/deta',20,-2.5d0,2.5d0,'Mgg'
!     &          ,'fb/','LIN',1)
!     
!        call bookplot(3,'ds/Mgg',50,100.d0,160.d0,'Mgg'
!     &          ,'fb/','LIN',1)
!     
!        call bookplot(4,'ds/dPHI',20,0.d0,3.16d0,'Mgg'
!     &          ,'fb/','LIN',1)
!   
!        call bookplot(5,'ds/dpT12',40,0.d0,2.d0,'pT','fb/',
!     &          'LIN',1)
!     
!        call bookplot(6,'ds/dpT',40,0.d0,2.d0,'pT','fb/',
!     &          'LIN',1)
!     
!        call bookplot(7,'ds/dpT12',30,0.d0,500.d0,'pT','fb/',
!     &          'LIN',1)
!     
!        call bookplot(8,'ds/dpT12',60,0.d0,500.d0,'pT','fb/',
!     &          'LIN',1)
!       
!        call bookplot(9,'ds/dPHIGG',12,0.d0,3.1416d0,
!     &         'PHIGG' ,'fb/','LIN',1)
!    
!        call bookplot(10,'ds/dPHIGG',20,0.d0,3.1416d0,
!     &         'PHIGG' ,'fb/','LIN',1)
!     
!        call bookplot(11,'ds/dPHIGG',30,0.d0,3.1416d0,
!     &         'PHIGG' ,'fb/','LIN',1)
    ! 
!        call bookplot(12,'ds/dPHIGG',60,0.d0,3.1416d0,!
!     &         'PHIGG' ,'fb/','LIN',1)   
!
!        call bookplot(13,'ds/dCosP',20,-1.d0,1.d0,'Cos*'
!     &          ,'fb/','LIN',1)
!     
!        call bookplot(14,'ds/dCosP',25,-1.d0,1.d0,'Cos*'
!     &          ,'fb/','LIN',1)
!     
!        call bookplot(15,'ds/dCosP',40,-1.d0,1.d0,'Cos*'
!     &          ,'fb/','LIN',1)
!    
!        call bookplot(16,'ds/dCosP',50,-1.d0,1.d0,'Cos*'
!     &          ,'fb/','LIN',1)
!       
!        call bookplot(17,'ds/dCosA',25,-1.d0,1.d0,'Cos*'
!     &          ,'fb/','LIN',1)!

              return
        end



