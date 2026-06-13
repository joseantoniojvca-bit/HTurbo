      integer:: jets
      character(len=2):: jetlabel(mcfm_mxpart)
      common/parts_int/jets
      common/parts_char/jetlabel
!$omp threadprivate(/parts_int/,/parts_char/)
