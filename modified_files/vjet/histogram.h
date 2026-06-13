#ifndef HISTOHDR
#define HISTOHDR


#include <math.h>
#include "tdio.h"
#include "average.h"

/*
   Defines a histogram class. 
*/

class histogram : public output {
 private:
    int nBins;       // Number of Bins.

    double binSize;  // Size of each Bin.

    double minBin;   // Lower edge of first Bin.

    double maxBin;   // Upper edge of last Bin.
                     //     maxBin=minBin + nBins*binSize. 

    int* hcall;      // Array holding the number in 
                     //     each bin.

    double* hwgt;    // Array holding the sum of the weights
                     //     of the points within each bin.

    double* hwsq;    // Array holding the sum of the square
                     //     of the weights.

    int totIn;       // Total number of points falling in the
                     // histogram.

    int totOut;      // Total number of points falling outside the
                     // histogram.

    int totCalls;    // Total number of points thrown, whether or not
                     // they land in range of the histogram.

 public: 
    inline histogram(const double &min, const double &size, 
                     const int &numBins);
    inline histogram(ofstream &os, const double &min, const double &size, 
                     const int &numBins);
    ~histogram(void) {delete [] hcall; delete [] hwgt; delete [] hwsq;}

    inline void bin(double x, double wgt=1.0);  
                     // Bin the point x. The default weight is 1.0.
                     // To average a function f(x) within each bin,
                     // then use wgt = weight(x) * f(x). 

    inline void bin(int n, double *x, double *wgt);  
                     // Bin n points given in the array x, each
                     // with weight given in the array wgt.
                     // If 2 or more points fall in the same bin,
                     // their weight is treated as the sum of
                     // the weights.  This gives a more reasonable
                     // error estimate when using the subtraction
                     // method.  Note that each time binVegas is
                     // called, totalCalls increases by 1, not n.

    inline void print(double norm=1.0);
                     // Print in TopDraw format the histogram
                     // with some normalization (norm). 

    inline void printVegas(double norm=1.0);
                     // Print in TopDraw format the histogram
                     // with some normalization (norm).
                     // The points were assumed to be binned 
                     // using binVegas. 

    inline void printAvg(double &totAvg, double &totSd);
                     // Print in TopDraw format the average and
                     // standard deviation of the function f(x) within
                     // each bin.  It also sets totAvg and totSd
                     // to the average and standard deviation of
                     // the function over the total range of the
                     // histogram.

    inline int totalIn(void) {return totIn;}

    inline int totalOut(void) {return totOut;}

    inline int totalCalls(void) {return totCalls;}
};


inline histogram::histogram(const double &min, 
                     const double &size, const int &numBins)
     : output()
{
   nBins=numBins;
   binSize=size;
   minBin=min;
   maxBin=minBin + binSize * (double) nBins;
   hcall = new int[nBins];
   hwgt = new double[nBins];
   hwsq= new double[nBins];
   totIn = totOut = totCalls = 0;
   for(int i=0;i<nBins;i++) {
      hcall[i] = 0; 
      hwgt[i] = hwsq[i] = 0.0;
   } 
}

inline histogram::histogram(ofstream &os, const double &min, 
                     const double &size, const int &numBins)
     : output(os)
{
   nBins=numBins;
   binSize=size;
   minBin=min;
   maxBin=minBin + binSize * (double) nBins;
   hcall = new int[nBins];
   hwgt = new double[nBins];
   hwsq= new double[nBins];
   totIn = totOut = totCalls = 0;
   for(int i=0;i<nBins;i++) {
      hcall[i] = 0; 
      hwgt[i] = hwsq[i] = 0.0;
   } 
}


inline void histogram::bin(double x, double wgt)
{
   totCalls +=1;
   if (x>minBin && x<maxBin) 
   {
      int i = (int) ((x-minBin)/binSize);
      hcall[i] += 1;
      hwgt[i] += wgt;
      hwsq[i] += wgt*wgt;
      totIn += 1;
   }
   else totOut += 1;
}   

inline void histogram::bin(int n, double *x, double *wgt)
{
   totCalls +=1;

   double *wgtSum = new double[n];
   int *i = new int[n];

   for(int j=0;j<n;j++) {
     if (x[j]>minBin && x[j]<maxBin) 
     {
        i[j] = (int) ((x[j]-minBin)/binSize);
	wgtSum[j] = wgt[j];
     }
     for(int k=0;k<j;k++) {
       if (x[k]>minBin && x[k]<maxBin) 
	 {
           if (i[j]==i[k]) {
             wgtSum[k]+=wgtSum[j];
             wgtSum[j]=0.0;
           }
         }
      }
   }
   for(int j=0;j<n;j++) {
     if (x[j]>minBin && x[j]<maxBin) 
     {
        hwgt[i[j]] +=wgtSum[j];
        hwsq[i[j]] +=wgtSum[j]*wgtSum[j];
     }
   }
   delete [] wgtSum;
   delete [] i;
}   

inline void histogram::print(double norm)
{
   double binCenter = minBin-binSize/2.0;
   for(int i=0;i<nBins;i++) {
      binCenter += binSize;
      pData(binCenter, hwgt[i]/norm, sqrt(hwsq[i])/norm);
   }
}

inline void histogram::printVegas(double norm)
{
   double binCenter = minBin-binSize/2.0;
   double dcalls = (double) totCalls;
   for(int i=0;i<nBins;i++) {
      binCenter += binSize;
      double mean = hwgt[i]/dcalls;
      double sd = sqrt(fabs(hwsq[i]/dcalls-mean*mean)/(dcalls-1.0));
      pData(binCenter, mean/norm, sd/norm);
   }
}

inline void histogram::printAvg(double &totAvg, double &totSd)
{
   double binCenter = minBin-binSize/2.0;
   double avg, sd;
   totAvg=0.0;
   double totSq=0.0;
   for(int i=0;i<nBins;i++) {
      totAvg += hwgt[i];
      totSq +=  hwsq[i];
      binCenter += binSize;
      if (hcall[i] == 0) {
           avg = 0.0;
           sd = 0.0;
      } else if (hcall[i] == 1) {
           avg = hwgt[i];
           sd = 0.0;
      } else {
           double calls = (double) hcall[i];
           avg = hwgt[i]/calls;
           sd = sqrt(fabs(hwsq[i]/calls-avg*avg)/(calls-1.0));
      }
      pData(binCenter, avg, sd);
   }
   if (totIn > 1) {
        double sum = (double) totIn;
        totAvg /= sum;
        totSd = sqrt(fabs(totSq/sum-totAvg*totAvg)/(sum-1.0));
   } 
}


#endif
