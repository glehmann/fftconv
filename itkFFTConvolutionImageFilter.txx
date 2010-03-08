/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilter_txx
#define __itkFFTConvolutionImageFilter_txx

#include "itkFFTConvolutionImageFilter.h"
#include "itkProgressAccumulator.h"
#include "itkFlipImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"

namespace itk {

template <class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::FFTConvolutionImageFilter()
{
  m_Normalize = true;
  m_GreatestPrimeFactor = 13;
  m_PadMethod = ZERO_FLUX_NEUMANN;
  this->SetNumberOfRequiredInputs(2);
}

template <class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
void 
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  InputImageType * input0 = const_cast<InputImageType *>(this->GetInput());
  InputImageType * input1 = const_cast<KernelImageType *>(this->GetKernelImage());
  if ( !input0 || !input1 )
    { 
    return;
    }
  input0->SetRequestedRegion( input0->GetLargestPossibleRegion() );
  input1->SetRequestedRegion( input1->GetLargestPossibleRegion() );
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
void
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::GenerateData()
{
  this->AllocateOutputs();
  const InputImageType * input = this->GetInput();
  const KernelImageType * kernel = this->GetKernelImage();
  OutputImageType * output = this->GetOutput();

  typedef typename itk::Image< FFTPrecisionType, ImageDimension > InternalImageType;
  
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);

  typedef itk::FlipImageFilter< KernelImageType > FlipType;
  typename FlipType::Pointer flip = FlipType::New();
  flip->SetInput( kernel );
  typename FlipType::FlipAxesArrayType axes;
  axes.Fill( true ); // we must flip all the axes
  flip->SetFlipAxes( axes );
  flip->SetNumberOfThreads( this->GetNumberOfThreads() );
  flip->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( flip, 0.005f );

  typedef itk::ImageToImageFilter< KernelImageType, InternalImageType > NormType;
  typedef itk::NormalizeToConstantImageFilter< KernelImageType, InternalImageType > NormConstType;
  typedef itk::CastImageFilter< KernelImageType, InternalImageType > CastType;
  typename NormType::Pointer norm;
  if( m_Normalize )
    {
    norm = NormConstType::New();
    }
  else
    {
    norm = CastType::New();
    }
  norm->SetInput( flip->GetOutput() );
  norm->SetNumberOfThreads( this->GetNumberOfThreads() );
  norm->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( norm, 0.005f );

  typedef itk::FFTPadImageFilter< InputImageType, InternalImageType, InternalImageType, InternalImageType > PadType;
  typename PadType::Pointer pad = PadType::New();
  pad->SetInput( input );
  pad->SetInputKernel( norm->GetOutput() );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  pad->SetReleaseDataFlag( true );
  pad->SetPadMethod( m_PadMethod );
  progress->RegisterInternalFilter( pad, 0.05f );

  typedef itk::FFTRealToComplexConjugateImageFilter< FFTPrecisionType, ImageDimension > FFTType;
  typename FFTType::Pointer fft = FFTType::New();
  fft->SetInput( pad->GetOutput() );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( fft, 0.25f );

  // vnl filters need a size which is a power of 2
  if( std::string(fft->GetNameOfClass()).find("Vnl") == 0 )
    {
    pad->SetPadToPowerOfTwo( true );
    // fake the cyclic behavior
    if( m_PadMethod == NO_PADDING )
      {
      pad->SetPadMethod( WRAP );
      }
    }
  else
    {
    pad->SetGreatestPrimeFactor( m_GreatestPrimeFactor );
    }

  typedef itk::FFTShiftImageFilter< InternalImageType, InternalImageType > ShiftType;
  typename ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pad->GetOutputKernel() );
  shift->SetNumberOfThreads( this->GetNumberOfThreads() );
  shift->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( shift, 0.04f );
  
  typename FFTType::Pointer fftk = FFTType::New();
  fftk->SetInput( shift->GetOutput() );
  fftk->SetNumberOfThreads( this->GetNumberOfThreads() );
  fftk->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( fftk, 0.25f );
  
  typedef itk::MultiplyImageFilter< typename FFTType::OutputImageType,
                                    typename FFTType::OutputImageType,
                                    typename FFTType::OutputImageType > MultType;
  typename MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, fftk->GetOutput() );
  mult->SetNumberOfThreads( this->GetNumberOfThreads() );
  mult->SetReleaseDataFlag( true );
  mult->SetInPlace( true );
  progress->RegisterInternalFilter( mult, 0.1f );
  
  typedef itk::FFTComplexConjugateToRealImageFilter< FFTPrecisionType, ImageDimension > IFFTType;
  typename IFFTType::Pointer ifft = IFFTType::New();
  ifft->SetInput( mult->GetOutput() );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );
  progress->RegisterInternalFilter( ifft, 0.25f );
  
  typedef itk::RegionFromReferenceImageFilter< InternalImageType, OutputImageType > CropType;
  typename CropType::Pointer crop = CropType::New();
  crop->SetInput( ifft->GetOutput() );
  crop->SetReferenceImage( input );
  crop->SetNumberOfThreads( this->GetNumberOfThreads() );
  progress->RegisterInternalFilter( crop, 0.05f );
  
  crop->GraftOutput( output );
  crop->Update();
  this->GraftOutput( crop->GetOutput() );
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TFFTPrecision>
void
FFTConvolutionImageFilter<TInputImage, TKernelImage, TOutputImage, TFFTPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}

}// end namespace itk
#endif
