set(DOCUMENTATION "This module contains ImageIO classes for reading via GDAL (http://www.gdal.org)")

itk_module(ITKIOGDAL
  DEPENDS
    ITKGDAL
    ITKIOImageBase
  TEST_DEPENDS
    ITKTestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
