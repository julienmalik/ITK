/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkInverseDeconvolutionImageFilter_h
#define __itkInverseDeconvolutionImageFilter_h

#include "itkFFTConvolutionImageFilter.h"

namespace itk
{
/** \class InverseDeconvolutionImageFilter
 * \brief The direct linear inverse deconvolution filter.
 *
 * The inverse filter is the most straightforward deconvolution
 * method. Considering that convolution of two images in the spatial domain is
 * equivalent to multiplying the Fourier transform of the two images,
 * the inverse filter consists of inverting the multiplication. In
 * other words, this filter computes the following:
 * \f[ hat{F}(\omega) =
 *       \begin{cases}
 *         G(\omega) / H(\omega) & \text{if $|H(\omega)| \geq \epsilon$} \\
 *         0                     & \text{otherwise}
 *       \end{cases}
 * \f]
 * where \f$\hat{F}(\omega)\f$ is the Fourier transform of the
 * estimate produced by this filter, \f$G(\omega)\f$ is the Fourier
 * transform of the input blurred image, \f$H(\omega)\f$ is the
 * Fourier transform of the blurring kernel, and \f$\epsilon\f$ is a
 * constant real non-negative threshold (called
 * KernelZeroMagnitudeThreshold in this filter) that determines when
 * the magnitude of a complex number is considered zero.
 *
 * \author Gaeten Lehmann, Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France
 * \author Cory Quammen, The University of North Carolina at Chapel Hill
 *
 * \ingroup ITKDeconvolution
 *
 */
template< class TInputImage, class TKernelImage = TInputImage, class TOutputImage = TInputImage, class TInternalPrecision=double >
class ITK_EXPORT InverseDeconvolutionImageFilter :
  public FFTConvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
{
public:
  typedef InverseDeconvolutionImageFilter                        Self;
  typedef FFTConvolutionImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                                   Pointer;
  typedef SmartPointer< const Self >                             ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information ( and related methods ) */
  itkTypeMacro(InverseDeconvolutionImageFilter, FFTConvolutionImageFilter);

  /** Dimensionality of input and output data is assumed to be the same. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  typedef TInputImage                           InputImageType;
  typedef TOutputImage                          OutputImageType;
  typedef TKernelImage                          KernelImageType;
  typedef typename Superclass::InputPixelType   InputPixelType;
  typedef typename Superclass::OutputPixelType  OutputPixelType;
  typedef typename Superclass::KernelPixelType  KernelPixelType;
  typedef typename Superclass::InputIndexType   InputIndexType;
  typedef typename Superclass::OutputIndexType  OutputIndexType;
  typedef typename Superclass::KernelIndexType  KernelIndexType;
  typedef typename Superclass::InputSizeType    InputSizeType;
  typedef typename Superclass::OutputSizeType   OutputSizeType;
  typedef typename Superclass::KernelSizeType   KernelSizeType;
  typedef typename Superclass::SizeValueType    SizeValueType;
  typedef typename Superclass::InputRegionType  InputRegionType;
  typedef typename Superclass::OutputRegionType OutputRegionType;
  typedef typename Superclass::KernelRegionType KernelRegionType;

  /** Internal image types. */
  typedef typename Superclass::InternalImageType               InternalImageType;
  typedef typename Superclass::InternalImagePointerType        InternalImagePointerType;
  typedef typename Superclass::InternalComplexType             InternalComplexType;
  typedef typename Superclass::InternalComplexImageType        InternalComplexImageType;
  typedef typename Superclass::InternalComplexImagePointerType InternalComplexImagePointerType;

  /** Set/get the threshold value uused to determine whether a
  * frequency of the Fourier transform of the blurring kernel is
  * considered to be zero. Default value is 1.0e-4. */
  itkSetMacro(KernelZeroMagnitudeThreshold, double);
  itkGetConstMacro(KernelZeroMagnitudeThreshold, double);

protected:
  InverseDeconvolutionImageFilter();
  ~InverseDeconvolutionImageFilter() {}

  /** This filter uses a minipipeline to compute the output. */
  void GenerateData();

  virtual void PrintSelf(std::ostream & os, Indent indent) const;

private:
  InverseDeconvolutionImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                  //purposely not implemented

  double m_KernelZeroMagnitudeThreshold;
};

namespace Functor
{
template< class TInput1, class TInput2, class TOutput >
class InverseDeconvolutionFunctor
{
public:
  InverseDeconvolutionFunctor() { m_KernelZeroMagnitudeThreshold = 0.0; }
  ~InverseDeconvolutionFunctor() {}

  bool operator!=( const InverseDeconvolutionFunctor & ) const
  {
    return false;
  }
  bool operator==( const InverseDeconvolutionFunctor & other) const
  {
    return !(*this != other);
  }
  inline TOutput operator()(const TInput1 & I, const TInput2 & H) const
  {
    double absH = std::abs( H );
    TOutput value = NumericTraits< TOutput >::ZeroValue();
    if ( absH >= m_KernelZeroMagnitudeThreshold )
      {
      value = static_cast< TOutput >( I / H );
      }
    return value;
  }

  /** Set/get the threshold value below which complex magnitudes are considered
   * to be zero. */
  void SetKernelZeroMagnitudeThreshold(double mu)
  {
    m_KernelZeroMagnitudeThreshold = mu;
  }
  double GetKernelZeroMagnitudeThreshold() const
  {
    return m_KernelZeroMagnitudeThreshold;
  }

private:
   double m_KernelZeroMagnitudeThreshold;
};
} //namespace Functor

}

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkInverseDeconvolutionImageFilter.hxx"
#endif

#endif
