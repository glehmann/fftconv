/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFastMutexLock.cxx,v $
  Language:  C++
  Date:      $Date: 2003-09-10 14:29:07 $
  Version:   $Revision: 1.10 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

  Portions of this code are covered under the VTK copyright.
  See VTKCopyright.txt or http://www.kitware.com/VTKCopyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkFFTWLock.h"
#if defined(USE_FFTWF) || defined(USE_FFTWD)
#include "fftw3.h"
#endif

namespace itk
{
SimpleFastMutexLock FFTWLock::m_Lock;
FFTWLock            FFTWLock::m_Singleton;

FFTWLock
::FFTWLock()
{ 
  // std::cout << "======== init fftw stuff =========" << std::endl;
#if defined(USE_FFTWF)
  fftwf_init_threads();
#endif
#if defined(USE_FFTWD)
  fftw_init_threads();
#endif
}
  
FFTWLock
::~FFTWLock()
{
  // std::cout << "======== cleanup fftw stuff =========" << std::endl;
#if defined(USE_FFTWF)
  fftwf_cleanup_threads();
  fftwf_cleanup();
#endif
#if defined(USE_FFTWD)
  fftw_cleanup_threads();
  fftw_cleanup();
#endif
}
  
void
FFTWLock
::Lock()
{
  FFTWLock::m_Lock.Lock();
}
  
void
FFTWLock
::Unlock()
{
  FFTWLock::m_Lock.Unlock();
}
  

}//end namespace itk
