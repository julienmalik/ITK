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
#ifndef __itkLandweberDeconvolutionImageFilter_hxx
#define __itkLandweberDeconvolutionImageFilter_hxx

#include "itkLandweberDeconvolutionImageFilter.h"

namespace itk
{

template< class TInputImage, class TKernelImage, class TOutputImage >
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::LandweberDeconvolutionImageFilter()
{
  m_Alpha = 0.1;
  m_TransformedInput = NULL;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::~LandweberDeconvolutionImageFilter()
{
  m_TransformedInput = NULL;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::Initialize(ProgressAccumulator * progress, float progressWeight,
             float iterationProgressWeight)
{
  this->Superclass::Initialize( progress, 0.5f * progressWeight,
                                iterationProgressWeight );

  this->PrepareInput( this->GetInput(), m_TransformedInput, progress,
                      0.5f * progressWeight );

  // Set up minipipeline to compute estimate at each iteration
  m_LandweberFilter = LandweberFilterType::New();
  // Transform of current estimate will be set as input 1 in Iteration()
  m_LandweberFilter->SetInput2( this->m_TransferFunction );
  m_LandweberFilter->SetInput3( m_TransformedInput );
  m_LandweberFilter->GetFunctor().m_Alpha = m_Alpha;
  progress->RegisterInternalFilter( m_LandweberFilter,
                                    0.3f * iterationProgressWeight );

  m_IFFTFilter = IFFTFilterType::New();
  m_IFFTFilter->SetActualXDimensionIsOdd( this->GetXDimensionIsOdd() );
  m_IFFTFilter->SetInput( m_LandweberFilter->GetOutput() );
  progress->RegisterInternalFilter( m_IFFTFilter,
                                    0.7f * iterationProgressWeight );
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::Iteration(ProgressAccumulator * progress, float iterationProgressWeight)
{
  // Set up minipipeline to compute the new estimate
  InternalComplexImagePointerType transformedEstimate;
  this->TransformPaddedInput( this->m_CurrentEstimate,
                              transformedEstimate, progress,
                              0.1f * iterationProgressWeight );

  // Set the inputs
  m_LandweberFilter->SetInput1( transformedEstimate );

  // Trigger the update
  m_IFFTFilter->UpdateLargestPossibleRegion();

  // Store the current estimate
  this->m_CurrentEstimate = m_IFFTFilter->GetOutput();
  this->m_CurrentEstimate->DisconnectPipeline();

  progress->ResetFilterProgressAndKeepAccumulatedProgress();
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::Finish(ProgressAccumulator * progress, float progressWeight)
{
  this->Superclass::Finish( progress, progressWeight );

  m_LandweberFilter = NULL;
  m_IFFTFilter = NULL;
}

template< class TInputImage, class TKernelImage, class TOutputImage >
void
LandweberDeconvolutionImageFilter< TInputImage, TKernelImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  this->Superclass::PrintSelf( os, indent );

  os << indent << "Alpha: " << m_Alpha << std::endl;
}
} // end namespace itk

#endif
