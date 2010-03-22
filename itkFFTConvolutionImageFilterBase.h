/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkFFTConvolutionImageFilterBase.h,v $
  Language:  C++
  Date:      $Date: 2009-01-28 18:14:36 $
  Version:   $Revision: 1.7 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkFFTConvolutionImageFilterBase_h
#define __itkFFTConvolutionImageFilterBase_h

#include "itkImageToImageFilter.h"
#include "itkConceptChecking.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkProgressAccumulator.h"

namespace itk {

/** \class FFTConvolutionImageFilterBase
 * \brief Pad two images with zeros to make them suitable for a convolution in the
 * frequency domain.
 *
 * The requires two input images. The first one is supposed to be the image to convolve
 * and the second one the kernel.
 *
 * The filter pad with zeros the first image by the size of the second image. The second image
 * is centered and padded with zeros to have the same size than the first output image.
 * After the transform, both images have the same LargestPossibleRegion.
 *
 * The option PadToPowerOfTwo can be set to true to force the size of the
 * images be a power of two - if the size of the padded image is 274 without this option, it
 * would be increased to 512 when PadToPowerOfTwo is true.
 * This option is makes the images usable with vnl's implementation of FFT.
 * PadToPowerOfTwo is false by default.
 * 
 * \author Gaetan Lehmann
 *
 * \sa FFTShiftImageFilter NormalizeToConstantImageFilter FFTRealToComplexConjugateImageFilter
 */
template<class TInputImage, class TKernelImage=TInputImage, class TOutputImage=TInputImage, class TInternalPrecision=float>
class ITK_EXPORT FFTConvolutionImageFilterBase : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef FFTConvolutionImageFilterBase Self;

  typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;

  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage                              InputImageType;
  typedef TKernelImage                             KernelImageType;
  typedef KernelImageType                          PointSpreadFunctionType;
  typedef TOutputImage                             OutputImageType;
  typedef TInternalPrecision                       InternalPrecisionType;
  typedef typename InputImageType::Pointer         InputImagePointer;
  typedef typename InputImageType::ConstPointer    InputImageConstPointer;
  typedef typename InputImageType::PixelType       InputImagePixelType;
  typedef typename KernelImageType::Pointer        KernelImagePointer;
  typedef typename KernelImageType::ConstPointer   KernelImageConstPointer;
  typedef typename KernelImageType::PixelType      KernelImagePixelType;
  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  typedef typename InputImageType::RegionType      RegionType;
  typedef typename InputImageType::IndexType       IndexType;
  typedef typename InputImageType::SizeType        SizeType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  typedef itk::FFTRealToComplexConjugateImageFilter< InternalPrecisionType, ImageDimension > FFTFilterType;
  typedef itk::FFTComplexConjugateToRealImageFilter< InternalPrecisionType, ImageDimension > IFFTFilterType;
  typedef typename FFTFilterType::OutputImageType                                            ComplexImageType;
  typedef typename ComplexImageType::PixelType                                               ComplexType;
  typedef typename ComplexImageType::Pointer                                                 ComplexImagePointerType;

  typedef typename itk::Image< InternalPrecisionType, ImageDimension > InternalImageType;
  typedef typename InternalImageType::Pointer                          InternalImagePointerType;

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(FFTConvolutionImageFilterBase, ImageToImageFilter);

   /** Set the kernel image */
  void SetKernelImage(const TKernelImage *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<TKernelImage *>(input) );
    }

  /** Get the kernel image */
  const KernelImageType * GetKernelImage() const
    {
    return static_cast<KernelImageType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

  /** Set the input image */
  void SetInput1(const TInputImage *input)
    {
    this->SetInput( input );
    }

  /** Set the kernel image */
  void SetInput2(const TKernelImage *input)
    {
    this->SetKernelImage( input );
    }

   /** Set the kernel image */
  void SetPointSpreadFunction(const PointSpreadFunctionType *input)
    {
    // Process object is not const-correct so the const casting is required.
    this->SetNthInput( 1, const_cast<PointSpreadFunctionType *>(input) );
    }

  /** Get the kernel image */
  const PointSpreadFunctionType * GetPointSpreadFunction() const
    {
    return static_cast<PointSpreadFunctionType*>(
      const_cast<DataObject *>(this->ProcessObject::GetInput(1)));
    }

  /**
   * Set/Get the greatest prime factor allowed on the size of the padded image.
   * The filter increase the size of the image to reach a size with the greatest
   * prime factor smaller or equal to the specified value. The default value is
   * 13, which is the greatest prime number for which the FFT are precomputed
   * in FFTW, and thus gives very good performance.
   * A greatest prime factor of 2 produce a size which is a power of 2, and thus
   * is suitable for vnl base fft filters.
   * A greatest prime factor of 1 or less - typically 0 - disable the extra padding.
   *
   * Warning: this parameter is not used (and useful) only when ITK is built with
   * FFTW support.
   */
  itkGetConstMacro(GreatestPrimeFactor, int);
  itkSetMacro(GreatestPrimeFactor, int);
  
  /**
   * Set/Get the padding method.
   */
  typedef enum { NO_PADDING=0, ZERO_FLUX_NEUMANN=1, ZERO=2, MIRROR=3, WRAP=4 } PadMethod;
  itkGetConstMacro(PadMethod, int);
  itkSetMacro(PadMethod, int);
  
  /**
   * Set/Get whether the kernel should be normalized to one or not.
   * Default is true.
   */
  itkGetConstMacro(Normalize, bool);
  itkSetMacro(Normalize, bool);
  itkBooleanMacro(Normalize);
  
#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasPixelTraitsCheck,
    (Concept::HasPixelTraits<InputImagePixelType>));
  itkConceptMacro(InputHasNumericTraitsCheck,
    (Concept::HasNumericTraits<InputImagePixelType>));
  /** End concept checking */
#endif


protected:
  FFTConvolutionImageFilterBase();
  ~FFTConvolutionImageFilterBase() {};

  void GenerateInputRequestedRegion();
  
  void Init( InternalImagePointerType & paddedInput, InternalImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight=0.66 );
  void Init( ComplexImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight=0.66 );
  void Init( InternalImagePointerType & paddedInput, ComplexImagePointerType & paddedKernel, bool & xIsOdd, float progressWeight=0.66 );
  
  void RegisterInternalFilter( ProcessObject * filter, float progressWeight );
  
  void End( ComplexImageType * paddedOutput, bool xIsOdd, float progressWeight=0.34 );
  void End( InternalImageType * paddedOutput, float progressWeight=0.34 );

  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  FFTConvolutionImageFilterBase(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  int                          m_GreatestPrimeFactor;
  int                          m_PadMethod;
  bool                         m_Normalize;
  ProgressAccumulator::Pointer m_ProgressAccumulator;

}; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkFFTConvolutionImageFilterBase.txx"
#endif

#endif
