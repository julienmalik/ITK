itk_module_test()
set(ITKIOGDALTests
itkGDALImageIOTest.cxx
itkIOGDALHeaderTest.cxx
)

CreateTestDriver(ITKIOGDAL  "${ITKIOGDAL-Test_LIBRARIES}" "${ITKIOGDALTests}")

add_test(NAME itkIOGDALHeaderTest
      COMMAND ITKIOGDALTestDriver itkIOGDALHeaderTest)

# TIFF
add_test(NAME itkGDALImageIOTestTIFF1
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/cthead1.tif}
              ${ITK_TEST_OUTPUT_DIR}/cthead1.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/cthead1.tif} ${ITK_TEST_OUTPUT_DIR}/cthead1.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF2
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImage.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImage.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImage.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImage.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF3
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageCCITTFax3.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax3.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageCCITTFax3.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax3.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF4
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageCCITTFax4.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax4.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageCCITTFax4.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax4.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF5
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageGray.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageGray.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageGray.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageGray.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF6
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageHuffmanRLE.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageHuffmanRLE.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageHuffmanRLE.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageHuffmanRLE.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF7
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageJPEG.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageJPEG.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageJPEG.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageJPEG.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF8
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageLZW.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageLZW.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageLZW.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageLZW.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF9
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageNone.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageNone.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageNone.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageNone.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF10
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImagePackBits.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImagePackBits.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImagePackBits.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImagePackBits.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF11
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageZIP.tif}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageZIP.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageZIP.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageZIP.tif 2 1)
add_test(NAME itkGDALImageIOTestTIFF12
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/ShortTestImage.tiff}
              ${ITK_TEST_OUTPUT_DIR}/ShortTestImage.tiff
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/ShortTestImage.tif}f ${ITK_TEST_OUTPUT_DIR}/ShortTestImage.tiff 2 2)
add_test(NAME itkGDALImageIOTestTIFF13
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/48BitTestImage.tif}
              ${ITK_TEST_OUTPUT_DIR}/48BitTestImage.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/48BitTestImage.tif} ${ITK_TEST_OUTPUT_DIR}/48BitTestImage.tif 2 2)
add_test(NAME itkGDALImageIOTestTIFF14
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageCCITTFax3.mha}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax3.mha
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageCCITTFax3.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax3.mha 2 1)
add_test(NAME itkGDALImageIOTestTIFF15
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/RGBTestImageCCITTFax4.mha}
              ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax4.mha
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/RGBTestImageCCITTFax4.tif} ${ITK_TEST_OUTPUT_DIR}/RGBTestImageCCITTFax4.mha 2 1)
add_test(NAME itkGDALImageIOMultiPagesTestTIFF1
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/ramp.tif}
              ${ITK_TEST_OUTPUT_DIR}/ramp.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/ramp.tif} ${ITK_TEST_OUTPUT_DIR}/ramp.tif 3 1)
add_test(NAME itkGDALImageIOMultiPagesTestTIFF2
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/rampShort.tif}
              ${ITK_TEST_OUTPUT_DIR}/rampShort.tif
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/rampShort.tif} ${ITK_TEST_OUTPUT_DIR}/rampShort.tif 3 3)

# BMP
add_test(NAME itkGDALImageIOTestBMP1
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Input/smallRGBA.bmp}
              ${ITK_TEST_OUTPUT_DIR}/smallRGBA.bmp
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/smallRGBA.bmp} ${ITK_TEST_OUTPUT_DIR}/smallRGBA.bmp 2 1)
add_test(NAME itkGDALImageIOTestBMP2
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/cthead1.bmp}
              ${ITK_TEST_OUTPUT_DIR}/cthead1.bmp
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/cthead1.bmp} ${ITK_TEST_OUTPUT_DIR}/cthead1.bmp 2 1)
add_test(NAME itkGDALImageIOTestBMP3
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/cthead254x254.bmp}
              ${ITK_TEST_OUTPUT_DIR}/cthead254x254.bmp
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/cthead254x254.bmp} ${ITK_TEST_OUTPUT_DIR}/cthead254x254.bmp 2 1)
add_test(NAME itkGDALImageIOTestBMP4
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/image_grayscale.bmp}
              ${ITK_TEST_OUTPUT_DIR}/image_grayscale.bmp
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/image_grayscale.bmp} ${ITK_TEST_OUTPUT_DIR}/image_grayscale.bmp 2 1)
add_test(NAME itkGDALImageIOTestBMP5
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/image_color.bmp}
              ${ITK_TEST_OUTPUT_DIR}/image_color.bmp
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/image_color.bmp} ${ITK_TEST_OUTPUT_DIR}/image_color.bmp 2 1)

# JPEG
add_test(NAME itkGDALImageIOTestJPEG1
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/cthead1.png}
              ${ITK_TEST_OUTPUT_DIR}/cthead1.png
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/cthead1.jpg} ${ITK_TEST_OUTPUT_DIR}/cthead1.png 2 1)
add_test(NAME itkGDALImageIOTestJPEG2
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Baseline/IO/cthead1.jpg}
              ${ITK_TEST_OUTPUT_DIR}/cthead1.jpg
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/cthead1.jpg} ${ITK_TEST_OUTPUT_DIR}/cthead1.jpg 2 1)

# PNG
add_test(NAME itkGDALImageIOTestPNG1
      COMMAND ITKIOGDALTestDriver itkGDALImageIOTest
              DATA{${ITK_DATA_ROOT}/Input/cthead1.png} ${ITK_TEST_OUTPUT_DIR}/itkGDALImageIOTest1.png 2 1)
add_test(NAME itkGDALImageIOTestPNG2
      COMMAND ITKIOGDALTestDriver itkGDALImageIOTest
              DATA{${ITK_DATA_ROOT}/Input/VisibleWomanEyeSlice.png} ${ITK_TEST_OUTPUT_DIR}/itkGDALImageIOTest2.png 2 1)

#JPEG2000
add_test(NAME itkGDALImageIOTestJP2K_1
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Input/Bretagne1.j2k}
              ${ITK_TEST_OUTPUT_DIR}/Bretagne1_01.png
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/Bretagne1.j2k} ${ITK_TEST_OUTPUT_DIR}/Bretagne1_01.png 2 1)
add_test(NAME itkGDALImageIOTestJP2K_2
      COMMAND ITKIOGDALTestDriver
    --compare DATA{${ITK_DATA_ROOT}/Input/Cevennes2.jp2}
              ${ITK_TEST_OUTPUT_DIR}/Cevennes2_01.png
    itkGDALImageIOTest DATA{${ITK_DATA_ROOT}/Input/Cevennes2.jp2} ${ITK_TEST_OUTPUT_DIR}/Cevennes2_01.png 2 1)
add_test(NAME itkGDALImageIOTestJP2K_3
      COMMAND ITKIOGDALTestDriver
    --compare ${ITK_TEST_OUTPUT_DIR}/Bretagne1_01.png
              ${ITK_TEST_OUTPUT_DIR}/Bretagne1_02.j2k
    itkGDALImageIOTest ${ITK_TEST_OUTPUT_DIR}/Bretagne1_01.png ${ITK_TEST_OUTPUT_DIR}/Bretagne1_02.j2k 2 1)
set_tests_properties( itkGDALImageIOTestJP2K_3 PROPERTIES DEPENDS itkGDALImageIOTestJP2K_1 )
