/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTWCommon.h,v $
  Language:  C++
  Date:      $Date: 2008-12-21 19:13:11 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTWCommon_h
#define __itkFFTWCommon_h

#if defined(USE_FFTWF) || defined(USE_FFTWD)
#include "fftw3.h"
#endif
#include "itkFFTWLock.h"

namespace itk
{
namespace fftw
{
/**
 * \class Interface
 * \brief Wrapper for FFTW API
 */
template <typename TPixel>
class  Proxy
{
  // empty -- only double and float specializations work
protected:
  Proxy() {}; 
  ~Proxy() {};
};

#if defined(USE_FFTWF)

template <>
class Proxy<float>
{
public:
  typedef float         PixelType;
  typedef fftwf_complex ComplexType;
  typedef fftwf_plan    PlanType;
  typedef Proxy<float>  Self;

  static PlanType Plan_dft_c2r_1d(int n,
                                  ComplexType *in,
                                  PixelType *out,
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_c2r_1d(n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r_2d(int nx, 
                                  int ny,
                                  ComplexType *in,
                                  PixelType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_c2r_2d(nx,ny,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r_3d(int nx, 
                                  int ny,
                                  int nz,
                                  ComplexType *in, 
                                  PixelType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_c2r_3d(nx,ny,nz,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r(int rank, 
                               const int *n,
                               ComplexType *in, 
                               PixelType *out, 
                               unsigned flags,
                               int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_c2r(rank,n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }

  static PlanType Plan_dft_r2c_1d(int n,
                                  PixelType *in,
                                  ComplexType *out,
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_r2c_1d(n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c_2d(int nx, 
                                  int ny,
                                  PixelType *in,
                                  ComplexType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_r2c_2d(nx,ny,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c_3d(int nx, 
                                  int ny,
                                  int nz,
                                  PixelType *in, 
                                  ComplexType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_r2c_3d(nx,ny,nz,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c(int rank, 
                               const int *n,
                               PixelType *in, 
                               ComplexType *out, 
                               unsigned flags,
                               int threads=1)
    {
    FFTWLock::Lock();
    fftwf_plan_with_nthreads(threads);
    PlanType plan = fftwf_plan_dft_r2c(rank,n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }

  static void Execute(PlanType p)
    {
    fftwf_execute(p);
    }
  static void DestroyPlan(PlanType p)
    {
    fftwf_destroy_plan(p);
    }
};

#endif // USE_FFTWF


#if defined(USE_FFTWD)
template <> 
class Proxy<double>
{
public:
  typedef double        PixelType;
  typedef fftw_complex  ComplexType;
  typedef fftw_plan     PlanType;
  typedef Proxy<double> Self;

  static PlanType Plan_dft_c2r_1d(int n,
                                  ComplexType *in,
                                  PixelType *out,
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_c2r_1d(n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r_2d(int nx, 
                                  int ny,
                                  ComplexType *in,
                                  PixelType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_c2r_2d(nx,ny,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r_3d(int nx, 
                                  int ny,
                                  int nz,
                                  ComplexType *in, 
                                  PixelType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_c2r_3d(nx,ny,nz,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_c2r(int rank, 
                               const int *n,
                               ComplexType *in, 
                               PixelType *out, 
                               unsigned flags,
                               int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_c2r(rank,n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }

  static PlanType Plan_dft_r2c_1d(int n,
                                  PixelType *in,
                                  ComplexType *out,
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_r2c_1d(n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c_2d(int nx, 
                                  int ny,
                                  PixelType *in,
                                  ComplexType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_r2c_2d(nx,ny,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c_3d(int nx, 
                                  int ny,
                                  int nz,
                                  PixelType *in, 
                                  ComplexType *out, 
                                  unsigned flags,
                                  int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_r2c_3d(nx,ny,nz,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static PlanType Plan_dft_r2c(int rank, 
                               const int *n,
                               PixelType *in, 
                               ComplexType *out, 
                               unsigned flags,
                               int threads=1)
    {
    FFTWLock::Lock();
    fftw_plan_with_nthreads(threads);
    PlanType plan = fftw_plan_dft_r2c(rank,n,in,out,flags);
    FFTWLock::Unlock();
    return plan;
    }
  static void Execute(PlanType p)
    {
    fftw_execute(p);
    }
  static void DestroyPlan(PlanType p)
    {
    fftw_destroy_plan(p);
    }
};

#endif
}
}
#endif
