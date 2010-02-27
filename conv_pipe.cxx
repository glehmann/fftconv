#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"

#include "itkFFTZeroPaddingImageFilter.h"
#include "itkNormalizeToConstantImageFilter.h"
#include "itkFFTShiftImageFilter.h"
#include "itkFFTRealToComplexConjugateImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkFFTComplexConjugateToRealImageFilter.h"
#include "itkRegionFromReferenceImageFilter.h"


int main(int argc, char * argv[])
{

  if( argc != 4 )
    {
    std::cerr << "usage: " << argv[0] << " intput kernel output" << std::endl;
    std::cerr << " input: the input image" << std::endl;
    std::cerr << " output: the output image" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }

  const int dim = 2;
  
  typedef float PType;
  typedef itk::Image< PType, dim > IType;

  typedef itk::ImageFileReader< IType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  // itk::SimpleFilterWatcher watcher_reader(reader, "reader");

  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( argv[2] );
  // itk::SimpleFilterWatcher watcher_reader2(reader2, "reader2");

  typedef itk::NormalizeToConstantImageFilter< IType, IType > NormType;
  NormType::Pointer norm = NormType::New();
  norm->SetInput( reader2->GetOutput() );
  itk::SimpleFilterWatcher watcher_norm(norm, "norm");

  typedef itk::FFTZeroPaddingImageFilter< IType, IType > PadType;
  PadType::Pointer pad = PadType::New();
  pad->SetInput( 0, reader->GetOutput() );
  pad->SetInput( 1, norm->GetOutput() );
  // TODO: check if this is needed
  pad->SetPadToPowerOfTwo( true );
  itk::SimpleFilterWatcher watcher_pad(pad, "pad");

  typedef itk::FFTShiftImageFilter< IType, IType > ShiftType;
  ShiftType::Pointer shift = ShiftType::New();
  shift->SetInput( pad->GetOutput(1) );
  itk::SimpleFilterWatcher watcher_shift(shift, "shift");
  
  typedef itk::FFTRealToComplexConjugateImageFilter< PType, dim > FFTType;
  FFTType::Pointer fft = FFTType::New();
  fft->SetInput( pad->GetOutput() );
  itk::SimpleFilterWatcher watcher_fft(fft, "fft");

  FFTType::Pointer fftk = FFTType::New();
  fftk->SetInput( shift->GetOutput() );
  itk::SimpleFilterWatcher watcher_fftk(fftk, "fftk");
  
  typedef itk::MultiplyImageFilter< FFTType::OutputImageType, FFTType::OutputImageType, FFTType::OutputImageType > MultType;
  MultType::Pointer mult = MultType::New();
  mult->SetInput( 0, fft->GetOutput() );
  mult->SetInput( 1, fftk->GetOutput() );
  itk::SimpleFilterWatcher watcher_mult(mult, "mult");
  
  typedef itk::FFTComplexConjugateToRealImageFilter< PType, dim > IFFTType;
  IFFTType::Pointer ifft = IFFTType::New();
  ifft->SetInput( mult->GetOutput() );
  itk::SimpleFilterWatcher watcher_ifft(ifft, "ifft");
  
  typedef itk::RegionFromReferenceImageFilter< IType, IType > CropType;
  CropType::Pointer crop = CropType::New();
  crop->SetInput( ifft->GetOutput() );
  crop->SetReferenceImage( reader->GetOutput() );
  itk::SimpleFilterWatcher watcher_crop(crop, "crop");

  typedef itk::ImageFileWriter< IType > WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput( crop->GetOutput() );
  writer->SetFileName( argv[3] );
  // itk::SimpleFilterWatcher watcher_writer(writer, "writer");
  
  writer->Update();

  return 0;
}

