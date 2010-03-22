/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilterBase.txx,v $
  Language:  C++
  Date:      $Date: 2007-01-20 20:05:58 $
  Version:   $Revision: 1.3 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilterBase_txx
#define __itkFFTConvolutionImageFilterBase_txx

#include "itkFFTConvolutionImageFilterBase.h"
#include "itkFlipImageFilter.h"
#include "itkFFTPadImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"
#include "itkIntensityWindowingImageFilter.h"

namespace itk {

template <class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::FFTConvolutionImageFilterBase()
{
  m_Normalize = true;
  m_GreatestPrimeFactor = 13;
  m_PadMethod = ZERO_FLUX_NEUMANN;
  this->SetNumberOfRequiredInputs(2);
}

template <class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void 
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
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


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::Init( InternalImagePointerType & paddedInput, InternalImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight )
{
  const InputImageType * input = this->GetInput();
  const KernelImageType * kernel = this->GetKernelImage();

  // Create a process accumulator for tracking the progress of this minipipeline
  m_ProgressAccumulator = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typedef itk::FlipImageFilter< KernelImageType > FlipType;
  typename FlipType::Pointer flip = FlipType::New();
  flip->SetInput( kernel );
  typename FlipType::FlipAxesArrayType axes;
  axes.Fill( true ); // we must flip all the axes
  flip->SetFlipAxes( axes );
  flip->SetNumberOfThreads( this->GetNumberOfThreads() );
  flip->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( flip, 0.33 * progressWeight );
    }

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
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( norm,  0.33 * progressWeight );
    }

  typedef itk::FFTPadImageFilter< InputImageType, InternalImageType, InternalImageType, InternalImageType > PadType;
  typename PadType::Pointer pad = PadType::New();
  pad->SetInput( input );
  pad->SetInputKernel( norm->GetOutput() );
  pad->SetNumberOfThreads( this->GetNumberOfThreads() );
  pad->SetReleaseDataFlag( true );
  pad->SetPadMethod( m_PadMethod );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( pad,  0.34 * progressWeight );
    }
  // instantiate a fft filter to know if it is a vnl implementation
  typename FFTFilterType::Pointer fft = FFTFilterType::New();
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
  
  pad->Update();
  
  paddedInput = pad->GetOutput();
  paddedInput->DisconnectPipeline();
  paddedKernel = pad->GetOutputKernel();
  paddedKernel->DisconnectPipeline();
  xIsOdd = paddedInput->GetLargestPossibleRegion().GetSize()[0] % 2;
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::Init( InternalImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight )
{
  InternalImagePointerType pk;
  
  this->Init( paddedInput, pk, xIsOdd, 0.15 * progressWeight );

  typedef itk::FFTShiftImageFilter< InternalImageType, InternalImageType > ShiftType;
  typename ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pk );
  shift->SetInverse( true );
  shift->SetNumberOfThreads( this->GetNumberOfThreads() );
  shift->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( shift, 0.05 * progressWeight );
    }
  
  typename FFTFilterType::Pointer fftk = FFTFilterType::New();
  fftk->SetInput( shift->GetOutput() );
  fftk->SetNumberOfThreads( this->GetNumberOfThreads() );
  fftk->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( fftk, 0.8 * progressWeight );
    }
  
  fftk->Update();
  paddedKernel = fftk->GetOutput();
  paddedKernel->DisconnectPipeline();
  
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::Init( ComplexImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight )
{
  // Create a process accumulator to track the progress of this minipipeline
  m_ProgressAccumulator = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  InternalImagePointerType pi;
  
  this->Init( pi, paddedKernel, xIsOdd, progressWeight * 0.6 );

  typename FFTFilterType::Pointer fft = FFTFilterType::New();
  fft->SetInput( pi );
  fft->SetNumberOfThreads( this->GetNumberOfThreads() );
  fft->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( fft, progressWeight * 0.4 );
    }

  fft->Update();
  paddedInput = fft->GetOutput();
  paddedInput->DisconnectPipeline();

}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::RegisterInternalFilter( ProcessObject * filter, float progressWeight )
{
  m_ProgressAccumulator->RegisterInternalFilter( filter, progressWeight );
}


template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::End( ComplexImageType * paddedOutput, bool xIsOdd, float progressWeight )
{
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typename IFFTFilterType::Pointer ifft = IFFTFilterType::New();
  ifft->SetInput( paddedOutput );
  ifft->SetActualXDimensionIsOdd( xIsOdd );
  ifft->SetNumberOfThreads( this->GetNumberOfThreads() );
  ifft->SetReleaseDataFlag( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( ifft, progressWeight * 0.8 );
    }
  
  this->End( ifft->GetOutput(), progressWeight * 0.2 );
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::End( InternalImageType * paddedOutput, float progressWeight )
{
  // Create a process accumulator for tracking the progress of this minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  m_ProgressAccumulator->SetMiniPipelineFilter(this);

  typedef itk::IntensityWindowingImageFilter< InternalImageType, OutputImageType > WindowType;
  typename WindowType::Pointer window = WindowType::New();
  window->SetInput( paddedOutput );
  window->SetWindowMinimum( NumericTraits< OutputImagePixelType >::Zero );
  window->SetWindowMaximum( NumericTraits< OutputImagePixelType >::max() );
  window->SetNumberOfThreads( this->GetNumberOfThreads() );
  window->SetReleaseDataFlag( true );
  window->SetInPlace( true );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( window, 0.1 * progressWeight );
    }

  typedef itk::RegionFromReferenceImageFilter< OutputImageType, OutputImageType > CropType;
  typename CropType::Pointer crop = CropType::New();
  crop->SetInput( window->GetOutput() );
  crop->SetReferenceImage( this->GetInput() );
  crop->SetNumberOfThreads( this->GetNumberOfThreads() );
  if( progressWeight != 0 )
    {
    m_ProgressAccumulator->RegisterInternalFilter( crop, 0.1 * progressWeight );
    }
  
  this->AllocateOutputs();
  crop->GraftOutput( this->GetOutput() );
  crop->Update();
  this->GraftOutput( crop->GetOutput() );
  m_ProgressAccumulator = NULL;
}

template<class TInputImage, class TKernelImage, class TOutputImage, class TInternalPrecision>
void
FFTConvolutionImageFilterBase<TInputImage, TKernelImage, TOutputImage, TInternalPrecision>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Normalize: "  << m_Normalize << std::endl;
  os << indent << "GreatestPrimeFactor: "  << m_GreatestPrimeFactor << std::endl;
  os << indent << "PadMethod: "  << m_PadMethod << std::endl;
}

}// end namespace itk
#endif
