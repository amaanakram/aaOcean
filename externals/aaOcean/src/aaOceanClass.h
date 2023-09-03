// aaOcean
// Author: Amaan Akram 
// www.amaanakram.com
//
// LICENSE: 
// aaOcean is free software and can be redistributed and modified under the terms of the 
// GNU General Public License (Version 3) as provided by the Free Software Foundation.
// GNU General Public License http://www.gnu.org/licenses/gpl.html
//
// A "New BSD" License for aaOcean can be obtained by contacting the author
// For more details on aaOcean and associated 3rd Party licenses, please see
// license.txt file that is part of the aaOcean repository:
// https://github.com/amaanakram/aaOcean

#ifndef AAOCEANCLASS_H
#define AAOCEANCLASS_H

#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fft.c"
#include "kissfft/tools/kiss_fftnd.h"
#include "kissfft/tools/kiss_fftnd.c"

class aaOcean
{
public:
    aaOcean();
    aaOcean(const aaOcean &cpy);
    ~aaOcean();

    // for retrieving array pointers in getOcean() from host app
    // declared here for convenience in identifying array names
    enum arrayType
    {
        eHEIGHTFIELD,
        eCHOPX,
        eCHOPZ,
        eFOAM,
        eEIGENPLUSX,
        eEIGENPLUSZ,
        eEIGENMINUSX,
        eEIGENMINUSZ,
        eSPECTRUM
    };

    // array for holding the current state of aaOcean object
    char m_state[512];

    // cleans up any left over data after ocean arrays are ready
    void clearResidualArrays();

    // main input function
    // should be called from host app
    void input(
        int     resolution,
        unsigned int  spectrum,
        unsigned int  seed,
        float   oceanScale,
        float   oceanDepth,
        float   surfaceTension,
        float   velocity,
        float   cutoff,
        float   windDir,
        int     windAlign,
        float   damp,
        float   waveSpeed,
        float   waveHeight,
        float   chopAmount,
        float   time,
        float   repeatTime,
        bool    doFoam,
        float   randWeight = 0.f,
        float   spectrumMult = 1.f,
        float   pmWaveSize = 1.f,
        float   jswpfetch = 100.f,
        float   swell = 0.0f);

    // main output function
    // should be called from host app
    float getOceanData(float uCoord, float vCoord, aaOcean::arrayType type) const;

    // retrieves input float arrays to write to open-exr format
    void getOceanArray(float *&outArray, aaOcean::arrayType type) const;

    // retrieves eigenvalue (foam) array bounds
    void getFoamBounds(float& outBoundsMin, float& outBoundsMax) const;

    // Is Valid if input is verified and 
    // ocean is correctly initialized
    bool isValid();

    // if choppiness is greater than 1, isChoppy() returns 1
    bool isChoppy();

    // returns resolution in powers of 2
    int getResolution();

    // Returns current ocean state
    // Needs better implementation
    char* getState();

    // initialization functions
private:
    int     m_resolution;     // resolution in powers of 2
    unsigned int m_seed;      // seed for random number generator
    unsigned int m_spectrum;  // philips, JONSWAP, PM , TMA etc.
    float   m_oceanScale;     // size of the ocean patch to generate in meters
    float   m_velocity;       // exposed as 'Wave Size' in some aaOcean plugins
    float   m_windDir;        // wind direction in degrees
    int     m_windAlign;      // stretches waves perpendicular to wind direction
    float   m_cutoff;         // cuts off waves smaller than this wavelength (smoothes ocean surface)
    float   m_damp;           // Damps waves travelling opposite the wind direction (wave reflection)
    float   m_chopAmount;     // squeezes the wave peaks to make them appear choppy
    float   m_waveHeight;     // wave height field multiplier
    float   m_waveSpeed;      // wave movement multiplier
    float   m_time;           // current time in seconds
    float   m_loopTime;       // time in seconds before the ocean shape repeats/loops
    float   m_oceanDepth;     // slows waves down with decreasing depth
    float   m_surfaceTension; // generates fast moving, high frequency waves
    float   m_foamBoundmin;   // stores eigenvalue array bounds
    float   m_foamBoundmax;   // stores eigenvalue array bounds
    float   m_randWeight;     // control blend between rand distributions

    // optional variables
    float   m_spectrumMult;         // multiplier for generated spectrum
    float   m_peakSharpening;   // JONSWAP Peak Sharpening
    float   m_jswpfetch;            // wind region
    float   m_swell;                // swell

    // ocean array pointers
    int     *m_xCoord;  // holds ocean grid coordinates
    int     *m_zCoord;  // holds ocean grid coordinates
    float   *m_hokReal; // real component of HoK (see Tessendorf paper)
    float   *m_hokImag; // imaginary component of HoK (see Tessendorf paper)
    float   *m_hktReal; // real component of HkT (see Tessendorf paper)
    float   *m_hktImag; // real component of HkT (see Tessendorf paper)
    float   *m_kX;      // x-component of wave vector
    float   *m_kZ;      // z-component of wave vector
    float   *m_omega;   // omega (see Tessendorf paper)
    float   *m_rand1;   // random number array 
    float   *m_rand2;   // random number array 
    float   *m_fftSpectrum; // spectrum array 

    // ocean output array pointers
    float *m_out_fft_htField;   // y displacement
    float *m_out_fft_chopX;     // x displacement
    float *m_out_fft_chopZ;     // z displacement
    float *m_out_fft_jxxX;      // eigenvector X component
    float *m_out_fft_jxxZ;      // eigenvector Z component
    float *m_out_fft_jzzX;      // eigenvector X component
    float *m_out_fft_jzzZ;      // eigenvector Z component
    float *m_out_fft_jxz;       // eigenvalue

    // array of pointers pointing to m_out* arrays
    // used with 'enum arrayType' and in getOceanData() and getOceanArray()
    float *m_arrayPointer[9];

    // bool types for various checks during run-time
    bool    m_isAllocated;      // working arrays memory allocation check
    bool    m_isFoamAllocated;  // eigenvalues/vectors memory allocation check
    bool    m_doSetup;          // triggers setupGrid()
    bool    m_doHoK;            // triggers evaluateHokData()
    bool    m_doChop;           // triggers evaluateChopField()
    bool    m_doFoam;           // triggers evaluateJacobians()

    // memory tracking -- needs better implementation
    int     m_memory;

    // kissfft arrays
    kiss_fft_cpx *m_fft_htField;
    kiss_fft_cpx *m_fft_chopX;
    kiss_fft_cpx *m_fft_chopZ;
    kiss_fft_cpx *m_fft_jxx; // eigenvector
    kiss_fft_cpx *m_fft_jzz; // eigenvector
    kiss_fft_cpx *m_fft_jxz; // eigenvalue

    kiss_fftnd_cfg m_planHeightField;
    kiss_fftnd_cfg m_planChopX;
    kiss_fftnd_cfg m_planChopZ;
    kiss_fftnd_cfg m_planJxx;
    kiss_fftnd_cfg m_planJxz;
    kiss_fftnd_cfg m_planJzz;

    void prepareOcean();
    void setupGrid();
    unsigned int generateUID(const float, const float) const;
    float philips(float k_sq);
    float piersonMoskowitz(float omega, float k_sq);
    float tma(float omega, int index = 0);
    float swell(float omega, float theta, float k_mag);

    // tessendorf ocean functions
    void evaluateHokData();
    void evaluateHieghtField();
    void evaluateChopField();
    void evaluateJacobians();
    void evaluateNormal();

    // interpolation functions
    inline float catmullRom(const float t, const float a, const float b, const float c, const float d) const;
    inline int wrap(int x) const;

    // memory management functions
    void allocateBaseArrays();
    void allocateFoamArrays();
    void clearArrays();
};

#endif  /* AAOCEANCLASS_H */

