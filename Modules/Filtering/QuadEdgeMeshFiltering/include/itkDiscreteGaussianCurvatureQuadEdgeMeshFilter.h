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
#ifndef __itkDiscreteGaussianCurvatureQuadEdgeMeshFilter_h
#define __itkDiscreteGaussianCurvatureQuadEdgeMeshFilter_h

#include "itkDiscreteCurvatureQuadEdgeMeshFilter.h"
#include "vnl/vnl_math.h"


namespace itk
{
/**
 * \class DiscreteGaussianCurvatureQuadEdgeMeshFilter
 * \brief see the following paper
 * title: Discrete Differential-Geometry Operators for Triangulated 2-Manifolds
 * authors: Mark Meyer, Mathieu Desbrun, Peter Schroder, Alan H. Barr
 * conference: VisMath '02
 * location: Berlin (Germany)
 * \author: Arnaud Gelas, Alexandre Gouaillard
 * \ingroup ITKQuadEdgeMeshFiltering
 */
template< class TInputMesh, class TOutputMesh=TInputMesh >
class ITK_EXPORT DiscreteGaussianCurvatureQuadEdgeMeshFilter:
  public DiscreteCurvatureQuadEdgeMeshFilter< TInputMesh, TOutputMesh >
{
public:
  typedef DiscreteGaussianCurvatureQuadEdgeMeshFilter   Self;
  typedef SmartPointer< Self >                          Pointer;
  typedef SmartPointer< const Self >                    ConstPointer;
  typedef DiscreteCurvatureQuadEdgeMeshFilter<
    TInputMesh, TOutputMesh >                           Superclass;

  typedef typename Superclass::InputMeshType    InputMeshType;
  typedef typename Superclass::InputMeshPointer InputMeshPointer;

  typedef typename Superclass::OutputMeshType                OutputMeshType;
  typedef typename Superclass::OutputMeshPointer             OutputMeshPointer;
  typedef typename Superclass::OutputPointsContainerPointer  OutputPointsContainerPointer;
  typedef typename Superclass::OutputPointsContainerIterator OutputPointsContainerIterator;
  typedef typename Superclass::OutputPointType               OutputPointType;
  typedef typename Superclass::OutputVectorType              OutputVectorType;
  typedef typename Superclass::OutputCoordType               OutputCoordType;
  typedef typename Superclass::OutputPointIdentifier         OutputPointIdentifier;
  typedef typename Superclass::OutputCellIdentifier          OutputCellIdentifier;
  typedef typename Superclass::OutputQEType                  OutputQEType;
  typedef typename Superclass::OutputMeshTraits              OutputMeshTraits;
  typedef typename Superclass::OutputCurvatureType           OutputCurvatureType;
  typedef typename Superclass::TriangleType                  TriangleType;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(DiscreteGaussianCurvatureQuadEdgeMeshFilter, DiscreteCurvatureQuadEdgeMeshFilter);

  /** New macro for creation of through a Smart Pointer   */
  itkNewMacro(Self);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( OutputIsFloatingPointCheck,
                   ( Concept::IsFloatingPoint< OutputCurvatureType > ) );
  /** End concept checking */
#endif

protected:
  DiscreteGaussianCurvatureQuadEdgeMeshFilter() {}
  ~DiscreteGaussianCurvatureQuadEdgeMeshFilter() {}

  virtual OutputCurvatureType EstimateCurvature(const OutputPointType & iP)
  {
    OutputMeshPointer output = this->GetOutput();

    OutputQEType *qe = iP.GetEdge();

    if ( qe != 0 )
      {
      OutputQEType *qe_it = qe;
      OutputQEType *qe_it2;

      OutputPointType q0, q1;

      OutputCurvatureType sum_theta = 0.;
      OutputCurvatureType area = 0.;

      do
        {
        // cell_id = qe_it->GetLeft();
        qe_it2 = qe_it->GetOnext();
        q0 = output->GetPoint( qe_it->GetDestination() );
        q1 = output->GetPoint( qe_it2->GetDestination() );

        // Compute Angle;
        sum_theta += static_cast< OutputCurvatureType >(
          TriangleType::ComputeAngle(q0, iP, q1) );
        area += this->ComputeMixedArea(qe_it, qe_it2);
        qe_it = qe_it2;
        }
      while ( qe_it != qe );

      return ( 2.0 * vnl_math::pi - sum_theta ) / area;
      }

    return 0.;
  }

private:
  DiscreteGaussianCurvatureQuadEdgeMeshFilter(const Self &); // purposely not
                                                                // implemented
  void operator=(const Self &);                                 // purposely not
                                                                // implemented
};
}

#endif
