/* -------------------------------------
   
  Various math routines in template format

  Coded by J. de la Cruz Rodriguez (ISP-SU 2016)

  -------------------------------------- */

#ifndef MATH_H
#define MATH_H

#include <cmath>
#include <algorithm>
#include <cstring>
#include <cstdio>

namespace mth{

  inline double sqr(double a){return a*a;};
  inline float sqr(float a){return a*a;};
  inline int sqr(int a){return a*a;};

  
  /* --- sqr --- */
  
  template <class T> T sqr(T a){return a*a;};
  /* ------------------------------------------------------------------------------- */

  /* --- Sum, mean, stdev --- */
  
  template <class T> T sum(size_t n, T* a){
    T tot = 0;
    for(size_t ii=0;ii<n;ii++) tot += a[ii];
    return tot;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> T mean(size_t n, T* a){
    T tot = sum<T>(n, a);
    return tot/n;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> double stdev(size_t n, T* a){
    double res = 0.0;
    double a_mean = sum<T>(n, a) / (double)n;
    //
    for(size_t ii=0; ii<n; ii++) res += sqr<double>((double)a[ii]-a_mean);
    //
    return sqrt(res/(n-1.0));
  }
  /* ------------------------------------------------------------------------------- */
  
  /* --- more accurate Kahan addition --- */
  
  template <class T> double ksum(size_t n, T* arr){
    
    long double sum = 0.0L, c = 0.0L;
    
    for(size_t kk = 0; kk<n; kk++){
      long double y = arr[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */

  template <class T> double ksum2(size_t n, T* arr){
    
    long double sum = 0.0L, c = 0.0L;
    
    for(size_t kk = 0; kk<n; kk++){
      long double y = arr[kk]*arr[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */

  /* --- Dot product of vectors --- */

  template <class T> double dot(size_t n, T *a, T *b){
    double res = 0.0;
    for(size_t ii=0; ii<n; ii++) res += a[ii]*b[ii];
    return res;
  }
  /* ------------------------------------------------------------------------------- */

  template <class T> double kdot(size_t n, T *a, T* b){
    long double sum = 0.0L, c = 0.0L;
    for(size_t kk = 0; kk<n; kk++){
      long double y = a[kk]*b[kk] - c;
      long double t = sum + y;
      c = (t - sum) - y;
      sum = t;
    }
    return (double)sum;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> size_t argmax(size_t n, T *a){
    T mm = 0;
    size_t idx = 0;
    for(size_t ii=1; ii<n; ii++) if(a[ii] > mm){
	mm = a[ii];
	idx = ii;
      }
    return idx;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> size_t argmin(size_t n, T *a){
    T mm = 0;
    size_t idx = 0;
    for(size_t ii=1; ii<n; ii++) if(a[ii] < mm){
	mm = a[ii];
	idx = ii;
      }
    return idx;
  }
  /* ------------------------------------------------------------------------------- */
  
  template <class T> void convolve1D(size_t n, T *sp, size_t npsf, T *psf){
    
    size_t npsf2 = npsf / 2, ii = 0;
    double *tmp = new double [n];
    
    for(size_t ww = 0; ww < n; ww++){
      
      /* --- Define limits to integrate with the PSF --- */
      
      int w0 = std::max((size_t)(ww-npsf2), (size_t)0), w1 = std::min((size_t)(ww+npsf2), (size_t)(n-1)), ii = w0 - (ww - npsf2);
      double sum = 0.0;

      /* --- Accumulate both the area of the PSF and the integrated data --- */
      
      for(int kk = w0; kk <= w1; kk++){
	tmp[ww] += psf[ii+kk-w0]*sp[kk];
	sum += psf[ii+kk-w0];
      }

      /* --- nomalize by the area of the PSF --- */
      
      tmp[ww] /= sum; 
    }

    
    /* --- Copy back array, cannot assume same type to use memcpy --- */
    
    for(size_t ii=0; ii<n; ii++) sp[ii] = (T)tmp[ii];

    
    /* --- Clean-up --- */
    
    delete [] tmp;
  }
  /* ------------------------------------------------------------------------------- */

  template <class T> void cent_der(size_t n, T *x, T *y, T *yp)
    {

      if(n == 1){
	yp[0] = 0;
	return;
      }
      
      double dx = x[1] - x[0], oder = 0, odx = 0;
      double der = (y[1] - y[0]) / dx;

      /* --- Fill in first point with regular finite difference --- */

      yp[0] = der;
      
      for(size_t k = 1; k<(n-1); k++){

	/* --- copy derivative from upwind interval --- */
	
	odx = dx, oder = der;

	/* --- Compute downwind derivative --- */
	
	dx = x[k+1] - x[k];
	der = (y[k+1] - y[k]) / dx;

	
	/* --- If not max/min then compute harmonic centered derivative --- */
	
	if(der*oder > 0.0){
	  double lambda = (1.0 + dx / (dx + odx)) / 3.0;
	  yp[k] = (der*oder) / ((1.0-lambda)*oder + lambda * der);
	}else{
	  yp[k] = 0;
	}
      } //k

      
      /* --- Fill in last point with regular finite difference --- */
      
      yp[n-1] = der; 
    }
  

  /* ------------------------------------------------------------------------------- */

  template <class T> void cmul(size_t n, T *arr, T c)
    {for(size_t ii=0;ii<n;ii++) arr[ii] *= c;};
  
  
  /* ------------------------------------------------------------------------------- */
  
  template <class T> void csum(size_t n, T *arr, T c)
    {for(size_t ii=0;ii<n;ii++) arr[ii] += c;};
  
  /* ------------------------------------------------------------------------------- */
  
  template <class T> void cdiv(size_t n, T *arr, T c)
    {for(size_t ii=0;ii<n;ii++) arr[ii] /= c;};
  
  /* ------------------------------------------------------------------------------- */
  template <class T> void csub(size_t n, T *arr, T c)
    {for(size_t ii=0;ii<n;ii++) arr[ii] -= c;};
  /* ------------------------------------------------------------------------------- */

  template <class T> void smooth(size_t n, T *arr, int win)
    {

      /* --- get init point --- */

      bool even = ((((win/2)*2) == win)?true:false);
      
      
      
      if(!even){
	T *tmp = new T [n];
	memcpy(tmp, arr, n*sizeof(T));
	
	int k0 = win/2, k1 = n-k0-1;
	double save0 = 0.0, save1 = 0.0, iwin = 1.0 / double(win);
	int w2 = win/2;
	
	for(int ii=0;ii<win;ii++){
	  save0 += arr[ii] * iwin;
	  save1 += arr[n-1-win+ii] * iwin;
	}

	double sum = save0;
	arr[k0] = sum;
	
	for(int ii=k0+1;ii<=k1;ii++){
	  sum += (tmp[ii+w2] - tmp[ii-w2-1]) * iwin;
	  arr[ii] = sum;
	}

	delete [] tmp;
      } // even
      
      
    }
  /* ------------------------------------------------------------------------------- */

  template <class T> inline int sign(const T val)
    {
      // Not exactly the same as sign in fortran, as it can return zero if val is zero.
      int res = (T(0) < val) - (val < T(0));
      
      if(res == 0) res = 1;
      return res;
    }

  /* ------------------------------------------------------------------------------- */
  
} // namespace


#endif
