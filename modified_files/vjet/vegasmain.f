       subroutine vegasmain()

             implicit real*8 (a-h,o-z)
             character * 25 prefn
             character * 25 fname
             character * 25 cppname
             character * 7 newver

             integer n0,i0,ncl0,dimenb,iifborn
             integer nn0,inn,nncl,dimenn,iinnlo
             integer fnn0,ifnn,fnncl,dimefnn,iifnnlo,iihyylo,iirhyy

             integer dimct,dimfini

             integer iset2g,ih12g,ih22g

             logical isave

             external inijet,wab

             external wnn,wfnn,wct,wctnn,wfini,wboxh,wh2
             external wggffg,wggHg,wggINTg,whint,wctgg,whintqq

             common/vname/prefn

             common/newver/newver

             common/dimenb/dimenb
             common/dimenn/dimenn
             common/dimefnn/dimefnn
             common/dimct/dimct
             common/dimfini/dimfini

             common/cppname/cppname
             common/restartcg/i0,inn,ifnn
             common/niter/n0,nn0,fnn0
             common/ncalls/ncl0,nncl,fnncl
             common/calcula/iifborn,iifnnlo,iinnlo,iihyylo,iirhyy 

             common/parpdfs2gamma/iset2g,ih12g,ih22g

             call setname
             call setPDF2gnnlo


c DIMENSION DEL ESPACIO DE INTEGRACION
             dimfini=16
             dimenb=7
             dimenn=16
             dimefnn=8 ! 7 del REAL y 7 del BORN + 1 CT
             dimct  =8
C  SI SAVE = TRUE GUARDA LA INFO
             isave = .true.
C INICIALIZA LAS VARIABLES "DE DESPUES" DE INTEGRAR

             ifiles=0

             avtot=0.d0
             dtot=0.d0

             avtot0 = 0.d0
             dtot0 =  0.d0
             av0 =0.d0
             dd0 =0.d0
             chi2a=0.d0


             avtotn = 0.d0
             dtotn =  0.d0
             avn =0.d0
             dn =0.d0
             chi2n=0.d0

             avtotfn = 0.d0
             dtotfn =  0.d0
             avfn =0.d0
             dfn =0.d0
             chi2fn=0.d0

             avct =0.d0
             dct =0.d0
             chi2ct=0.d0

             avctnn =0.d0
             dctnn =0.d0
             chi2ctnn=0.d0

             avfini =0.d0
             dfini =0.d0
             chi2fini=0.d0

             avbox =0.d0
             dbox =0.d0
             chi2box=0.d0

             avh2 =0.d0
             dh2 =0.d0
             chi2h2=0.d0

             avggffg =0.d0
             dggffg =0.d0
             chi2ggffg=0.d0

             avggHg =0.d0
             dggHg =0.d0
             chi2ggHg=0.d0

             avggINTg =0.d0
             dggINTg =0.d0
             chi2ggINTg=0.d0
             
             avhint =0.d0
             dhint =0.d0
             chi2hint=0.d0
             
             avctgg =0.d0
             dctgg =0.d0
             chi2ctgg=0.d0
             
             avhintqq =0.d0
             dhintqq =0.d0
             chi2hintqq=0.d0

C FIN DE LAS VARIABLES "DE DESPUES" DE INTEGRAR
C LLAMA A PARÁMETROS DEPENDIENTES DEL SISTEMA
             call sysdep

C SI IIBORN=1 CALCULA EL BORN Y EL H1^F
C SI IIH1=1 CALCULA EL H1 SIN EL H1^F
C SI IINLO=1 CALCULA EL NLO+CT

c            iifborn=0
c            iinnlo=1
c            iifnnlo=0

          !  write(*,*) 'iinnlo',iinnlo
          !  write(*,*) 'iihyylo',iihyylo
          !  write(*,*) 'iirhyy',iirhyy

           if(iifborn.eq.0) goto 100
           write(*,*) ' Start BORN  '
C INTEGRANDO PARA EL Full NNLO Cinematica BORN 2->2
          call strcat(prefn,'_lo',fname)
C COMIENZA EL PROCESO DE INTEGRACION
          call integrate
     #  (inijet,wab,fname, n0,i0,dimenb,
     #    ncl0,av0 ,dd0,chi2a,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'BORN  :  ',av0,' +- ',dd0

 100       if(iifnnlo.eq.0) goto 200


           write(*,*) ' Start CT NLO + H1'
C INTEGRANDO PARA EL Full NNLO Cinematca BORN 2->2
          call strcat(prefn,'_ctnlo',fname)
C COMIENZA EL PROCESO DE INTEGRACION
          call integrate
     #  (inijet,wct,fname,fnn0,ifnn,dimct,
     #    fnncl,avct,dct,chi2ct,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'CT NLO + H1:  ',avct,' +- ',dct


          write(*,*) 'Start Real NLO '
          call strcat(prefn,'_Rnlo',fname)
C COMIENZA EL PROCESO DE INTEGRACION
          call integrate
     #  (inijet,wfnn,fname, fnn0,ifnn,dimefnn,
     #    fnncl,avfn ,dfn,chi2fn,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'Real NLO:  ',avfn,' +- ',dfn


C INTEGRANDO PARA EL Full NNLO Cinematica 2->3 y 2->4 : 2Real + Dip +
c                                                       fini + RealNLO
 200      if(iinnlo.eq.0) goto 300



C INTEGRANDO PARA EL NNLO Cinematica 2->4 y 2->3: 2Real + Dip
          write(*,*) 'Start Real NNLO + Dip'
          call strcat(prefn,'_nnlo',fname)
C COMIENZA EL PROCESO DE INTEGRACION
         call integrate
     #  (inijet,wnn,fname, nn0,inn,dimenn,
     #    nncl,avn ,dn,chi2n,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'Real NNLO + Dip:  ',avn,' +- ',dn


C INTEGRANDO PARA EL NNLO Cinematica 2->3:fini
          write(*,*) 'Start Real NNLO (fini)'
          call strcat(prefn,'_fini',fname)
C COMIENZA EL PROCESO DE INTEGRACION
          call integrate
     #  (inijet,wfini,fname,fnn0,inn,dimfini,
     #    fnncl,avfini ,dfini,chi2fini,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'Real NNLO (fini) :  ',avfini,' +- ',dfini

C INTEGRANDO PARA EL NNLO Cinematica 2->2
          write(*,*) 'Start CTNNLO '
          call strcat(prefn,'_ctnn',fname)
C COMIENZA EL PROCESO DE INTEGRACION
              call integrate
     #  (inijet,wctnn,fname, fnn0,ifnn,dimct,
     #    fnncl,avctnn ,dctnn,chi2ctnn,isave)
            ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'CTNNLO :  ',avctnn,' +- ',dctnn


C INTEGRANDO PARA EL NNLO Cinematica 2->2
          write(*,*) 'Start  H2f + BOX'
          call strcat(prefn,'_boxh',fname)
C COMIENZA EL PROCESO DE INTEGRACION
                     call integrate
     #  (inijet,wboxh,fname, n0,i0,dimenb,
     #    ncl0,avbox ,dbox,chi2box,isave)
            ifiles=ifiles+1
        call suma(fname,ifiles)
         write(*,*) ' H2f + BOX :  ',avbox,' +- ',dbox

C INTEGRANDO PARA EL NNLO Cinematica 2->2
          write(*,*) 'Start  H2'
          call strcat(prefn,'_h2',fname)
c COMIENZA EL PROCESO DE INTEGRACION
              call integrate
     #  (inijet,wh2,fname, fnn0,ifnn,dimct,
     #    fnncl,avh2 ,dh2,chi2h2,isave)
            ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'H2 :  ',avh2,' +- ',dh2

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C INTEGRANDO PARA H int a LO
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           
 300      if(iihyylo.eq.0) goto 400
       
          write(*,*) 'Signal H and interference at LO (gg channel)'
          call strcat(prefn,'_HintLO',fname)
C COMIENZA EL PROCESO DE INTEGRACION whint
          call integrate
     #  (inijet,whint,fname,n0,i0,dimenb,
     #    ncl0,avhint,dHint,chi2Hint,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'Hint LO:  ',avHint,' +- ',dHint
         
  
c$$$         write(*,*)'Interference at LO (qqbar and qg channel)'
c$$$          call strcat(prefn,'_HintLOqq',fname)
c$$$C COMIENZA EL PROCESO DE INTEGRACION whint
c$$$          call integrate
c$$$     #  (inijet,whintqq,fname,fnn0,ifnn,dimefnn,
c$$$     #    ncl0,avhintqq,dHintqq,chi2Hintqq,isave)
c$$$          ifiles=ifiles+1
c$$$         call suma(fname,ifiles)
c$$$         write(*,*) 'Hint LO qqbar and qg :  ',avHintqq,' +- ',dHintqq
c$$$         
c$$$         
c$$$          write(*,*) ' Start CTgg NLO'
c$$$C INTEGRANDO PARA EL Full NNLO Cinematca BORN 2->2
c$$$          call strcat(prefn,'_ctgg',fname)
c$$$C COMIENZA EL PROCESO DE INTEGRACION
c$$$          call integrate
c$$$     #  (inijet,wctgg,fname,fnn0,ifnn,dimct,
c$$$     #    fnncl,avctgg,dctgg,chi2ctgg,isave)
c$$$          ifiles=ifiles+1
c$$$         call suma(fname,ifiles)
c$$$         write(*,*) 'CT NLO gg :  ',avctgg,' +- ',dctgg





Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C INTEGRANDO PARA H int REALS
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
 400      if(iirhyy.eq.0) goto 9100

c INTEGRANDOS PARA EL NNNLO (partial) Cinematica 2->3
c          write(*,*) 'Start Real gg -> gam gam g'
c          call strcat(prefn,'_ggffg',fname)
C COMIENZA EL PROCESO DE INTEGRACION
c          call integrate
c     #  (inijet,wggffg,fname, fnn0,inn,dimefnn,
c     #    fnncl,avggffg ,dggffg,chi2ggffg,isave)
c          ifiles=ifiles+1
c         call suma(fname,ifiles)
c         write(*,*) 'Real ggffg:  ',avggffg,' +- ',dggffg

c          write(*,*) 'Start Real gg -> (H->gam gam) g '
c          call strcat(prefn,'_ggHg',fname)
C COMIENZA EL PROCESO DE INTEGRACION
c          call integrate
c     #  (inijet,wgghg,fname, fnn0,inn,dimefnn,
c     #    fnncl,avgghg ,dggHg,chi2ggHg,isave)
c          ifiles=ifiles+1
c         call suma(fname,ifiles)
c         write(*,*) 'Real ggHg:  ',avggHg,' +- ',dggHg

          write(*,*) 'Start Real NLO interference'
          call strcat(prefn,'_ggINTg',fname)
C COMIENZA EL PROCESO DE INTEGRACION
          call integrate
     #  (inijet,wggintg,fname, fnn0,inn,dimefnn,
     #    fnncl,avggintg ,dggINTg,chi2ggINTg,isave)
          ifiles=ifiles+1
         call suma(fname,ifiles)
         write(*,*) 'Real ggINTg:  ',avggINTg,' +- ',dggINTg

Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



 9100         avtot=avctnn+avn+avfn+avct+av0+avfini+avbox+avh2
     #         +avggffg+avHint+avctgg+avHintqq+avgghg+avggintg
         dtot=sqrt(dctnn**2+dn**2+dfn**2+dct**2+dd0**2
     #    +dfini**2+dh2**2+dbox**2+dggffg**2+dHint**2+dctgg**2
     #    +dHintqq**2+dggHg**2+dggINTg**2)

           

c Write the .top file if the option was selected
       if(isave) then
          call inijet
c          call mclear
          call mclear2
          call strcat(prefn,'.top',fname)
c          fname='si.top'
          open(unit=99,file=fname,status=newver)
          call topout

          close(99)
        endif
c Write the total integral on the default output device
        write(*,*) 'Sigma(Tot) =  ', avtot,' +- ',dtot

         return

         end















C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando : Full NNLO cinematica Born: 2-> 2
           function wab(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision winteb
          integer dimenb,i
          common/dimenb/dimenb
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimenb
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandob(winteb,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=winteb*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wab=sig_evb


           return
            end


C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando : Full NNLO cinematica Born: 2-> 2
           function wboxh(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision winteb
          integer dimenb,i
          common/dimenb/dimenb
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimenb
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandoboxh(winteb,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=winteb*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wboxh=sig_evb


           return
            end

C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C    Integrando NNLO : 2Real + Dip + fin
           function wnn(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision winten
          integer dimenn,i
          common/dimenn/dimenn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimenn
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandon(winten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=winten*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wnn=sig_evb


           return
            end


C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C    Integrando NNLO : 2Real + Dip + fin
           function wfini(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision winten
          integer dimfini,i
          common/dimfini/dimfini
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimfini
             xxextern(i)=xx(i)
       !      write(*,*) 'xx(i)',xx(i),i
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandofini(winten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=winten*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wfini=sig_evb


           return
            end

C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando Full NNLO 2Real + dip + fini + RealNLO
           function wfnn(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision wfinten
          integer dimefnn,i
          common/dimefnn/dimefnn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimefnn
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandofullnlo(wfinten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=wfinten*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wfnn=sig_evb


           return
            end


C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando Full NNLO 2Real + dip + fini + RealNLO
           function wct(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision wwct
          integer dimct,i
          common/dimct/dimct
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimct
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandoctnlo(wwct,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=wwct*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wct=sig_evb


           return
            end



C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando Full NNLO 2Real + dip + fini + RealNLO
           function wctnn(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision wwctnn
          integer dimct,i
          common/dimct/dimct
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimct
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandoctnnlo(wwctnn,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=wwctnn*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wctnn=sig_evb


           return
            end




C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
C Integrando Full NNLO 2Real + dip + fini + RealNLO
           function wh2(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision wwh2
          integer dimct,i
          common/dimct/dimct
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimct
             xxextern(i)=xx(i)
            enddo
             

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandoh2(wwh2,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=wwh2*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wh2=sig_evb


           return
            end
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            

C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
       function wggffg(xx,wght)
         implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
         double precision xx(25),wght
         double precision wggffginten
         integer dimefnn,i
         common/dimefnn/dimefnn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
         common/xxrands/xxextern(25)
         save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

         do i=1,dimefnn
           xxextern(i)=xx(i)
         enddo
          
c         write(*,*) 'desde   function wggffg(xx,wght) ' 

C LLAMA AL INTEGRANDO DEFINIDO EN C++
         call integgffg(wggffginten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
         vv_evb=wggffginten*wght
         if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
            sig_evb=tot/wght
         endif


         wggffg=sig_evb

         return
       end
       
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C
C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
       function wgghg(xx,wght)
         implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
         double precision xx(25),wght
         double precision wgghginten
         integer dimefnn,i
         common/dimefnn/dimefnn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
         common/xxrands/xxextern(25)
         save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

         do i=1,dimefnn
           xxextern(i)=xx(i)
         enddo
         !     write(*,*) 'desde  function wgghg(xx,wght) ' 

C LLAMA AL INTEGRANDO DEFINIDO EN C++
         call integghg(wgghginten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
         vv_evb=wgghginten*wght
         if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
            sig_evb=tot/wght
         endif


         wgghg=sig_evb

         return
       end
       
       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       
       
C
C
C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
       function wggintg(xx,wght)
         implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
         double precision xx(25),wght
         double precision wggintginten
         integer dimefnn,i
         common/dimefnn/dimefnn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
         common/xxrands/xxextern(25)
         save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

         do i=1,dimefnn
           xxextern(i)=xx(i)
         enddo
c          write(*,*) 'desde function wggintg(xx,wght) ' 
C LLAMA AL INTEGRANDO DEFINIDO EN C++
         call integgintg(wggintginten,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
         vv_evb=wggintginten*wght
         if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
            sig_evb=tot/wght
         endif


         wggintg=sig_evb

         return
       end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc H production and signal LO  ------LEAN
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
C
C
C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
       function whint(xx,wght)
         implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
         double precision xx(25),wght
         double precision wgghyy
         integer dimenb,i
         common/dimenb/dimenb
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
         common/xxrands/xxextern(25)
         save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

         do i=1,dimenb
           xxextern(i)=xx(i)
         enddo
       !  write(*,*) 'sono qua'
C LLAMA AL INTEGRANDO DEFINIDO EN C++
         call integghyy(wgghyy,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
         vv_evb=wgghyy*wght
         if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
            sig_evb=tot/wght
         endif


         whint=sig_evb

         return
       end 

C
C
C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
C Y UN PESO QUE USAN DE MULETO LOS HISTOGRAMAS WGHT
       function whintqq(xx,wght)
         implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
         double precision xx(25),wght
         double precision wqqbarhyy
         integer dimefnn,i
         common/dimefnn/dimefnn
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
         common/xxrands/xxextern(25)
         save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

         do i=1,dimefnn
           xxextern(i)=xx(i)
         enddo
       !  write(*,*) 'sono qua'
C LLAMA AL INTEGRANDO DEFINIDO EN C++
         call inteqqbarhyy(wqqbarhyy,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
         vv_evb=wqqbarhyy*wght
         if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
            sig_evb=tot/wght
         endif


         whintqq=sig_evb

         return
       end 


C ESTA FUNCION ES LA QUE COLECTA TODOS LOS PESOS
C PARA INTEGRAR EL BORN
C POR PARÁMETRO Sscott E PASAN LOS RANDOMS NUMBERS XX[I]
           function wctgg(xx,wght)
           implicit real*8 (a-e,f-z)
c La dimension en los rand's debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL
          double precision xx(25),wght
          double precision wwctgg
          integer dimct,i
          common/dimct/dimct
c La dimension en los rand's externos debe ser la misma que se
c eligió para la dimensión de la llamada a INTEGRAL.
c Estos rand's se pasan por common al archivo de C++
c que los recupera como "struct".
          common/xxrands/xxextern(25)
          save/xxrands/
c Copia los rand's en los del common para pasarlos a C++

            do i=1,dimct
             xxextern(i)=xx(i)
            enddo

C LLAMA AL INTEGRANDO DEFINIDO EN C++
          call integrandoctggnlo(wwctgg,wght)

C SI ESTÁ BIEN LLAMA A LAS RUTINAS DE HISTOS
             vv_evb=wwctgg*wght
               if (vv_evb.ne.0.d0) then
            call outall(tot,vv_evb)
             sig_evb=tot/wght
               endif


              wctgg=sig_evb


           return
            end
       
C
C
C
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C SUBROUTINE AL PEDO SI LAS HAY SOLO CHEQUEA CUTS
C CONVENDRIA REDEFINIR
C LOS CUTS SE IMPLEMENTASN EN C++
C LLAMA A OUTFUN PARA LLENAR LOS HISTOS
       subroutine outall(tot,res)
        implicit real * 8 (a-h,o-z)
c Total for the events, to be returned to the calling function
          tot = 0.d0
          if( (abs(res).lt.10000d6)) then
            if(res.ne.0) then
             tot = res

             call outfun(tot,iout)
             if (iout.eq.0) tot=0d0
            endif
          else
            write(6,*)res,' OJO en outall, res excede 1000d6 '
          endif

       end

C ESTA LLENA LOS HISTOS CON LOS PESOS Y LAS VARIABLES
C QUE CORRESPONDEN

       subroutine outfun(www,iout)
         implicit real * 8 (a-h,o-z)


              iout=1


              return
       end

