#ifndef SetVEGAS_H
#define SetVEGAS_H 1



extern "C" {

extern struct {
	
	int i0,inn,ifnn;
	
	
} restartcg_;

extern struct {
	
	int n0,nn0,fnn0;
	
	
} niter_;

extern struct {
	
	int ncl0,nncl,fnncl;
	
	
} ncalls_;


extern struct {
	
	int iifborn,iifnnlo,iinnlo,iihyylo,iirhyy;
	
	
} calcula_;

extern struct {
	
	char *cppname[25];
	
	
} cppname_;


}

 void charge_VEGAS(_Evento event){
	 
//	 int sizename = sizeof(event.cppname[25]);
	// event.cppname[sizename--] = '\0';
	 
	 
	 
//	 cppname_.cppname[sizename] = event.cppname[sizename];
	 
	
	 
	 restartcg_.i0   = event.i0;
	 restartcg_.inn  = event.inn;
	 restartcg_.ifnn = event.ifnn;
	 
	 niter_.n0   = event.n0;
	 niter_.nn0  = event.nn0;
	 niter_.fnn0 = event.fnn0;
	 
	 ncalls_.ncl0   = event.ncl0;
	 ncalls_.nncl   = event.nncl;
	 ncalls_.fnncl  = event.fnncl;
	 
	 calcula_.iifborn   = event.iifborn;
	 calcula_.iifnnlo   = event.iifnnlo;
	 calcula_.iinnlo    = event.iinnlo;
	 //calcula_.iilean    = event.iilean;
	calcula_.iihyylo   = event.iihyylo;
	calcula_.iirhyy    = event.iirhyy;
	 
	 
 }

#endif


