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
#include "itkGDALImageIOFactory.h"
#include "itkCreateObjectFunction.h"
#include "itkGDALImageIO.h"
#include "itkVersion.h"

namespace itk
{
GDALImageIOFactory::GDALImageIOFactory()
{
  this->RegisterOverride( "itkImageIOBase",
                          "itkGDALImageIO",
                          "GDAL Image IO",
                          1,
                          CreateObjectFunction< GDALImageIO >::New() );
}

GDALImageIOFactory::~GDALImageIOFactory()
{
}

const char *
GDALImageIOFactory::GetITKSourceVersion(void) const
{
  return ITK_SOURCE_VERSION;
}

const char *
GDALImageIOFactory::GetDescription(void) const
{
  return "GDAL ImageIO Factory, allows the loading of GDAL images into insight";
}

// Undocumented API used to register during static initialization.
// DO NOT CALL DIRECTLY.

static bool GDALImageIOFactoryHasBeenRegistered;

void GDALImageIOFactoryRegister__Private(void)
{
  if( !GDALImageIOFactoryHasBeenRegistered )
    {
    GDALImageIOFactoryHasBeenRegistered = true;
    GDALImageIOFactory::RegisterOneFactory();
    }
}

} // end namespace itk
