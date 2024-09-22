#ifndef OPENEXR_OUTPUT_H
#define OPENEXR_OUTPUT_H

#include <ImfNamespace.h>
#include <ImfOutputFile.h>
#include <ImfInputFile.h>
#include <ImfChannelList.h>
#include <ImfStringAttribute.h>
#include <ImfMatrixAttribute.h>
#include <ImfArray.h>
#include <ImfArray.h>
#include <half.h>
#include <ImfRgbaFile.h>

namespace EXR = OPENEXR_IMF_NAMESPACE;
using namespace IMATH_NAMESPACE;

#if defined(_MSC_VER)
#include <io.h> // For access()
#else
#include <sys/io.h> // For access()
#include <stdlib.h>
#include <unistd.h>
#endif

inline bool dirExists(const char* path)
{
#ifdef _MSC_VER
    if (_access(path, 0) == 0)
    {
        return 1;
    }
    else
        return 0;
#else
    if (access(path, 0) == 0)
    {
        return 1;
    }
    else
        return 0;
#endif
}

void genFullFilePath(char* dest, const char* outputFolder, const char* postfix, const int frame)
{
#ifdef _MSC_VER
    char slash[5];
    strcpy(slash, "\\");
#else
    char slash[5];
    strcpy(slash, "/");
#endif

    char cFrame[10];
    char prefixFrame[10];
    if (frame < 10)
        strcpy(prefixFrame, "000");
    else if (frame < 100)
        strcpy(prefixFrame, "00");
    else if (frame < 1000)
        strcpy(prefixFrame, "0");
    sprintf(cFrame, "%s%d", prefixFrame, frame);

    strcpy(dest, outputFolder);
    strcat(dest, slash);

    if (strcmp(postfix, "") == 0)
        strcat(dest, "aaOceanData");
    else
        strcat(dest, "aaOceanData_");

    strcat(dest, postfix);
    strcat(dest, ".");
    strcat(dest, cFrame);
    strcat(dest, ".exr");
}

// writes floating point data to disk
void writeExr(const char fileName[], const float *rPixels, const float *gPixels,
    const float *bPixels, const float *aPixels, int width, int height)
{
    EXR::Header header(width, height);
    header.compression() = EXR::PIZ_COMPRESSION;

    header.channels().insert("R", EXR::Channel(EXR::FLOAT));
    header.channels().insert("G", EXR::Channel(EXR::FLOAT));
    header.channels().insert("B", EXR::Channel(EXR::FLOAT));
    header.channels().insert("A", EXR::Channel(EXR::FLOAT));

    EXR::OutputFile file(fileName, header);

    EXR::FrameBuffer frameBuffer;

    frameBuffer.insert("R",             // name
        EXR::Slice(EXR::FLOAT,          // type
        (char *)rPixels,                // base
            sizeof(*rPixels) * 1,       // xStride
            sizeof(*rPixels) * width)); // yStride

    frameBuffer.insert("G",             // name
        EXR::Slice(EXR::FLOAT,          // type
        (char *)gPixels,                // base
            sizeof(*gPixels) * 1,       // xStride
            sizeof(*gPixels) * width)); // yStride

    frameBuffer.insert("B",             // name
        EXR::Slice(EXR::FLOAT,          // type
        (char *)bPixels,                // base
            sizeof(*bPixels) * 1,       // xStride
            sizeof(*bPixels) * width)); // yStride

    frameBuffer.insert("A",             // name
        EXR::Slice(EXR::FLOAT,          // type
        (char *)aPixels,                // base
            sizeof(*aPixels) * 1,       // xStride
            sizeof(*aPixels) * width)); // yStride

    file.setFrameBuffer(frameBuffer);
    file.writePixels(height);
}

// writes floating single channel data point data to disk
void writeSingleChannelExr(const char fileName[], const float *rPixels, int width, int height)
{
    EXR::Header header(width, height);
    header.channels().insert("R", EXR::Channel(EXR::FLOAT));
    header.compression() = EXR::PIZ_COMPRESSION;

    EXR::OutputFile file(fileName, header);

    EXR::FrameBuffer frameBuffer;

    frameBuffer.insert("R",                     // name
        EXR::Slice(EXR::FLOAT,      // type
        (char *)rPixels,            // base
            sizeof(*rPixels) * 1,       // xStride
            sizeof(*rPixels) * width)); // yStride

    file.setFrameBuffer(frameBuffer);
    file.writePixels(height);
}

// prepares aaOcean arrays to be written to OpenEXR format
void oceanDataToEXR(aaOcean *&pOcean, const char *outputFolder, const char *postfix, int frame, char *outputFileName)
{
    int dimension = pOcean->getResolution();
    int arraySize = dimension * dimension;

    float *red, *green, *blue, *alpha = 0;
    green = (float*)malloc(arraySize * sizeof(float));
    pOcean->getOceanArray(green, aaOcean::eHEIGHTFIELD);

    if (pOcean->isChoppy())
    {
        red = (float*)malloc(arraySize * sizeof(float));
        blue = (float*)malloc(arraySize * sizeof(float));
        alpha = (float*)malloc(arraySize * sizeof(float));

        pOcean->getOceanArray(red, aaOcean::eCHOPX);
        pOcean->getOceanArray(blue, aaOcean::eCHOPZ);
        pOcean->getOceanArray(alpha, aaOcean::eFOAM);
    }

    EXR::Array2D<float> rPixels(dimension, dimension);
    EXR::Array2D<float> gPixels(dimension, dimension);
    EXR::Array2D<float> bPixels(dimension, dimension);
    EXR::Array2D<float> aPixels(dimension, dimension);

    #pragma omp parallel for
    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < dimension; j++)
        {
            gPixels[j][i] = green[i*dimension + j];
            if (red)
            {
                rPixels[j][i] = red[i*dimension + j];
                bPixels[j][i] = blue[i*dimension + j];
                aPixels[j][i] = alpha[i*dimension + j];
            }
            else
                rPixels[j][i] = bPixels[j][i] = aPixels[j][i] = 0.0f;
        }
    }

    genFullFilePath(&outputFileName[0], &outputFolder[0], &postfix[0], frame);
    writeExr(&outputFileName[0], &rPixels[0][0], &gPixels[0][0], &bPixels[0][0], &aPixels[0][0], dimension, dimension);

    free(green);
    if (red)
        free(red);
    if (blue)
        free(blue);
    if (alpha)
        free(alpha);
}
#endif /* OPENEXR_OUTPUT_H */
