/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTWLock.h,v $
  Language:  C++
  Date:      $Date: 2008-12-21 19:13:11 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTWLock_h
#define __itkFFTWLock_h

#include "itkSimpleFastMutexLock.h"

namespace itk
{
/**
 * A simple global lock which must be called before calling FFTW unsafe functions.
 * It also handle cleanly the initialization and cleanup of FFTW.
 */
class  FFTWLock
{

  public:
  
  static void Lock();
  static void Unlock();
  
  protected:
  
  FFTWLock();
  ~FFTWLock();
  
  static FFTWLock            m_Singleton;
  static SimpleFastMutexLock m_Lock;

};
}
#endif
