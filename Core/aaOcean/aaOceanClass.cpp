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

#ifndef AAOCEANCLASS_CPP
#define AAOCEANCLASS_CPP

#include <cmath>
#include <omp.h>
#include <climits>
#include <float.h>
#include <string>
#include <algorithm>
#include "dSFMT/dSFMT.c"
#include "constants.h"
#include "functionLib.h"
#include "timer.h"
#include "aaOceanClass.h"

aaOcean::aaOcean() :
    // input variables
    m_resolution(-1),
    m_spectrum(0),
    m_seed(1),
    m_windAlign(2),
    m_velocity(15.0f),
    m_windDir(45.0f),
    m_cutoff(-1.0f),
    m_damp(0.985f),
    m_oceanScale(100.0f),
    m_oceanDepth(10000.f),// forcing to be power of two, setting minimum resolution of 2^4
    m_surfaceTension(0.0f),
    m_chopAmount(1.0f),
    m_waveHeight(1.0f),
    m_waveSpeed(1.0f),
    m_time(1.0f),
    m_loopTime(10000.0f),
    m_foamBoundmin(FLT_MAX),
    m_foamBoundmax(-FLT_MAX),
    m_randWeight(0.f),
    m_spectrumMult(1.f),
    m_peakSharpening(1.f),
    m_jswpfetch(1.f),
    m_swell(0.0f),

    // working arrays
    m_xCoord(0),
    m_zCoord(0),
    m_hokReal(0),
    m_hokImag(0),
    m_hktReal(0),
    m_hktImag(0),
    m_kX(0),
    m_kZ(0),
    m_omega(0),
    m_rand1(0),
    m_rand2(0),
    m_fftSpectrum(0),

    // output arrays
    m_out_fft_htField(0),
    m_out_fft_chopX(0),
    m_out_fft_chopZ(0),
    m_out_fft_jxxX(0),
    m_out_fft_jxxZ(0),
    m_out_fft_jzzX(0),
    m_out_fft_jzzZ(0),
    m_out_fft_jxz(0),

    // bools to check ocean state
    m_isAllocated(0),
    m_isFoamAllocated(0),
    m_doHoK(0),
    m_doSetup(0),
    m_doChop(0),
    m_doFoam(0),
    m_isShader(0),

    // memory tracking
    m_memory(0),

    // fftw arrays
    m_fft_htField(0),
    m_fft_chopX(0),
    m_fft_chopZ(0),
    m_fft_jxx(0),
    m_fft_jzz(0),
    m_fft_jxz(0)
{
    strcpy(m_state, "");
}

aaOcean::aaOcean(const aaOcean &cpy)
{
    // empty copy constructor
    input(cpy.m_resolution,
        cpy.m_spectrum,
        cpy.m_seed,
        cpy.m_oceanScale,
        cpy.m_oceanDepth,
        cpy.m_surfaceTension,
        cpy.m_velocity,
        cpy.m_cutoff,
        cpy.m_windDir,
        cpy.m_windAlign,
        cpy.m_damp,
        cpy.m_waveSpeed,
        cpy.m_waveHeight,
        cpy.m_chopAmount,
        cpy.m_time,
        cpy.m_loopTime,
        cpy.m_doFoam
    );
}

aaOcean::~aaOcean()
{
    if(m_isShader)
        clearShaderArrays(true);
    else
        clearArrays();
}

bool aaOcean::isChoppy()
{
    return m_doChop;
}

bool aaOcean::isShader()
{
    return m_isShader;
}

void aaOcean::setShaderMode(bool mode)
{
    m_isShader = mode;
}

char* aaOcean::getState()
{
    return &m_state[0];
}

int aaOcean::getResolution()
{
    return m_resolution;
}

void aaOcean::input(
    int resolution,
    unsigned int spectrum,
    unsigned int seed,
    float oceanScale,
    float oceanDepth,
    float surfaceTension,
    float velocity,
    float cutoff,
    float windDir,
    int windAlign,
    float damp,
    float waveSpeed,
    float waveHeight,
    float chopAmount,
    float time,
    float loopTime,
    bool doFoam,
    float randWeight,
    float spectrumMult,
    float peakSharpening,
    float jswpfetch,
    float swell)
{
    Timer timer;

    // forcing to be power of two, setting minimum resolution of 2^4
    resolution = static_cast<int>(pow(2.0f, (4 + abs(resolution))));
    if (m_resolution != resolution)
    {
        m_resolution = resolution;
        m_max_threads = omp_get_max_threads();
        if(m_resolution < 1024)
            m_max_threads = 1;
        
        if(!m_isShader)
            allocateBaseArrays();
    }

    // scaled for better UI control
    if (spectrum < 2 )
    {
        // for Philips and PM spectrums
        waveHeight *= 0.01f;
        chopAmount *= 0.01f;
    }

    m_waveHeight    = waveHeight;
    m_time          = time;
    m_waveSpeed     = waveSpeed;
    m_doFoam        = doFoam;

    if (chopAmount > aa_EPSILON)
    {
        m_chopAmount = chopAmount;
        m_doChop = 1;
    }
    else
    {
        m_doChop = 0;
        m_chopAmount = 0.0f;
    }

    oceanScale  = maximum<float>(oceanScale, aa_EPSILON);
    velocity    = maximum<float>(velocity,   aa_EPSILON);
    oceanDepth  = maximum<float>(oceanDepth, aa_EPSILON);
    cutoff = fabs(cutoff * 0.01f);                      // scaling for better UI control
    windDir = windDir * aa_PIBY180;                     // to radians
    windAlign = std::max<int>((windAlign + 1) & ~1, 0); // forcing to even numbers
    damp = std::clamp(damp, 0.0f, 1.0f);

    if (m_oceanScale != oceanScale          ||
        m_spectrum != spectrum              ||
        m_oceanDepth != oceanDepth          ||
        m_surfaceTension != surfaceTension  ||
        m_windDir != windDir                ||
        m_cutoff != cutoff                  ||
        m_velocity != velocity              ||
        m_windAlign != windAlign            ||
        m_damp != damp                      ||
        m_loopTime != loopTime              ||
        m_spectrumMult != spectrumMult      ||
        m_peakSharpening != peakSharpening  ||
        m_jswpfetch != jswpfetch            ||
        m_swell != swell             )
    {
        m_oceanScale = oceanScale;
        m_spectrum = spectrum;
        m_oceanDepth = oceanDepth;
        m_surfaceTension = surfaceTension;
        m_windDir = windDir;
        m_cutoff = cutoff;
        m_velocity = velocity;
        m_windAlign = windAlign;
        m_damp = damp;
        m_loopTime = loopTime;
        m_spectrumMult = spectrumMult;
        m_peakSharpening = peakSharpening;
        m_jswpfetch = jswpfetch;
        m_swell = swell;

        // re-evaluate HoK if any of these vars change
        m_doHoK = 1;
    }

    if (m_seed != seed || m_randWeight != randWeight)
    {
        m_seed = seed;
        m_randWeight = randWeight;
        m_doSetup = m_doHoK = 1;
    }

    if(m_isShader)
    {
        shaderEvaluate();
        return;
    }

    if (!m_doHoK || !m_doSetup){
        snprintf(m_state, sizeof(m_state), "[aaOcean Core] Ocean base state unchanged. Re-evaluating ocean with cached data");
        prepareOcean();
    }
    else
    {
        prepareOcean();
        snprintf(m_state, sizeof(m_state), "[aaOcean Core] input processed. now proceeding with prepareOcean() at %dx%d resolution", m_resolution, m_resolution);
        timer.printElapsed(m_state, true);
    }
}

void aaOcean::prepareOcean()
{
    if (m_doSetup)
        setupGrid();

    if (m_doHoK)
        evaluateHokData();

    evaluateHieghtField();

    if (m_doChop)
        evaluateChopField();

    if (m_doFoam)
    {
        if (!m_isFoamAllocated)
            allocateFoamArrays();
        evaluateJacobians();
    }
    char temp_buffer[sizeof(m_state)];

    snprintf(temp_buffer, sizeof(temp_buffer), "%s", m_state);
    snprintf(m_state, sizeof(m_state), "%s\n[aaOcean Core] Working memory allocated: %.2f megabytes", temp_buffer, float(m_memory) / 1048576.f);
}

void aaOcean::shaderEvaluate()
{
    Timer timer;
    m_memory = 0;
    const int n = m_resolution;
    const int half_n = (-n / 2) - ((n - 1) / 2);
    const int nn = m_resolution * m_resolution;
    const int dims[2] = { m_resolution, m_resolution };

    m_fft_htField   = (kiss_fft_cpx *)malloc(nn * sizeof(kiss_fft_cpx));
    m_memory += nn * sizeof(kiss_fft_cpx);
    m_out_fft_htField = (float*)malloc(nn * sizeof(float));
    m_out_fft_chopX = (float*)malloc(nn * sizeof(float));
    m_out_fft_chopZ = (float*)malloc(nn * sizeof(float));
    m_memory += nn * sizeof(float) * 3;
    m_arrayPointer[eHEIGHTFIELD] = m_out_fft_htField;
    m_arrayPointer[eCHOPX] = m_out_fft_chopX;
    m_arrayPointer[eCHOPZ] = m_out_fft_chopZ;

    float *hokReal = m_out_fft_chopX;   // temporary storage
    float *hokImag = m_out_fft_chopZ;   // temporary storage

    if(m_doChop)
    {
        m_fft_chopX     = (kiss_fft_cpx *)malloc(nn * sizeof(kiss_fft_cpx));
        m_fft_chopZ     = (kiss_fft_cpx *)malloc(nn * sizeof(kiss_fft_cpx));
        m_memory        += nn * sizeof(kiss_fft_cpx) * 2;
    }

    // begin ocean spectrum build
    const float k_mult = aa_TWOPI / m_oceanScale;
    const float windx = cos(m_windDir);
    const float windz = sin(m_windDir);
    const float omega0 = aa_TWOPI / m_loopTime;
    const float signs[2] = { 1.0f, -1.0f };
    bool bDamp = false;
    if (m_damp > aa_EPSILON)
        bDamp = true;

    dsfmt_t dsfmt;
    #pragma omp parallel for firstprivate(dsfmt) num_threads(m_max_threads)
    for (int index = 0; index < nn; ++index)
    {
        unsigned int uID;
        float rand1, rand2, omega, spectrum;

        int i = index / n; // Calculate row index
        int j = index % n; // Calculate column index

        int xCoord = half_n + i * 2;
        int zCoord = half_n + j * 2;
        float x = static_cast<float>(xCoord);
        float z = static_cast<float>(zCoord);
        uID = static_cast<u_int32_t>(generateUID(x, z));
        dsfmt_init_gen_rand(&dsfmt, uID + m_seed);

        const float g1 = static_cast<float>(gaussian(dsfmt)); // gaussian rand draw 1
        const float g2 = static_cast<float>(gaussian(dsfmt)); // gaussian rand draw 2

        if (m_randWeight > aa_EPSILON)
        {
            const float u1 = static_cast<float>(uniform(dsfmt)); // uniform rand draw 1
            const float u2 = static_cast<float>(uniform(dsfmt)); // uniform rand draw 2

            rand1 = lerp(m_randWeight, g1, u1);
            rand2 = lerp(m_randWeight, g2, u2);
        }
        else
        {
            rand1 = g1;
            rand2 = g2;
        }

        float kX = (float)xCoord * k_mult;
        float kZ = (float)zCoord * k_mult;
        float k_sq = kX * kX + kZ * kZ;
        float k_mag = sqrt(k_sq);

        if (k_mag < aa_EPSILON)
                spectrum = aa_EPSILON;

        // calculate spectrum
        if (m_spectrum == 1) // Pierson-Moskowitz
        {
            // build dispersion relationship with oceanDepth relationship and capillary waves
            omega = sqrt(aa_GRAVITY * k_mag * tanh(k_mag * m_oceanDepth));
            spectrum = piersonMoskowitz(omega, k_sq);

        }
        else if (m_spectrum == 2) //  Texel MARSEN ARSLOE (TMA) SpectruM
        {
            // build dispersion relationship with oceanDepth relationship and capillary waves
            omega = sqrt(aa_GRAVITY * k_mag * tanh(k_mag * m_oceanDepth));
            spectrum = tma(omega, index);
        }
        else if (m_spectrum == 3) //  JONSWAP
        {
            // Constants for the simulation
            const float g = aa_GRAVITY;                     // Gravitational acceleration (m/s^2)
            const float alpha = 0.0081f;                    // Empirical constant for wave growth
            const float gamma = m_peakSharpening * 3.3f;    // Peak enhancement factor
            const float sigma_low = 0.07f;                  // Sigma for f < f_peak
            const float sigma_high = 0.09f;                 // Sigma for f > f_peak
            const float wind_speed = m_velocity * 10.f;     // Wind speed at 10 meters above sea level (m/s)
            const float fetch = m_jswpfetch * 100.f;        // Fetch distance in meters

            const float f_peak = (3.5f * std::pow(g, 0.67f)) / (std::pow(wind_speed, 0.34f) * std::pow(fetch, 0.33f));
            const float k_peak = 4.0f * aa_PI * aa_PI * f_peak * f_peak / g;
            float sigma = (k_mag <= k_peak) ? sigma_low : sigma_high;
            float S = alpha * g * g / std::pow(k_mag, 5.0f) * std::exp(-5.0f / 4.0f * std::pow(k_mag / k_peak, -4.0f));
            float peak_factor = std::pow(gamma, std::exp(-std::pow((k_mag - k_peak), 2) / (2.0f * sigma * sigma * k_peak * k_peak)));

            spectrum = std::sqrt(S * peak_factor / 2.0f) * m_spectrumMult * 0.01f;
        }
        else // Philips
        { 
            omega = aa_GRAVITY * k_mag * (1.0f + k_sq * m_surfaceTension * m_surfaceTension); // modifying dispersion for capillary waves
            omega = sqrt(omega);
            omega = (int(omega / omega0)) * omega0; // add time looping support with OmegaNought

            spectrum = sqrt(philips(k_sq));
        }

        // bias towards wind dir
        float k_dot_w = (kX * (1.0f / k_mag) * windx) + (kZ * (1.0f / k_mag) * windz);
        spectrum *= pow(k_dot_w, m_windAlign);  

        // eliminate wavelengths smaller than cutoff
        spectrum *= exp(-k_sq * m_cutoff);     

        // swell
        if (m_spectrum == 2)
            spectrum = swell(spectrum, k_dot_w, k_mag); 

       // reduce reflected waves
        if (bDamp && (k_dot_w < -aa_EPSILON))
            spectrum *= (1.0f - m_damp);

        m_fft_htField[index].r = omega;     // temporary storage
        m_fft_htField[index].i = kX;        // temporary storage
        m_out_fft_htField[index] = kZ;      // temporary storage

        hokReal[index] = aa_INV_SQRTTWO * rand1 * spectrum; // store hokReal
        hokImag[index] = aa_INV_SQRTTWO * rand2 * spectrum; // store hokImag
    }

    // prepare heightfield and chopfields for FFTs
    if(m_doFoam)
        allocateFoamArrays();
        
    const int n_sq_minus1 = n * n - 1;
    #pragma omp parallel for num_threads(m_max_threads)
    for (size_t index = 0; index < nn; ++index)
    {
        float omega = m_fft_htField[index].r;   // retrieve temp values
        float kX = m_fft_htField[index].i;      // retrieve temp values
        float kZ = m_out_fft_htField[index];    // retrieve temp values

        size_t index_rev = n_sq_minus1 - index; //tail end 
        float hokRealOpp = hokReal[index_rev];
        float hokImagOpp = hokImag[index_rev];

        float coswt = cos(omega * m_waveSpeed * m_time);
        float sinwt = sin(omega * m_waveSpeed * m_time);

        float hktReal =  (hokReal[index]    *  coswt) + (hokImag[index] * sinwt) +
                         (hokRealOpp        *  coswt) - (hokImagOpp *  sinwt);  //complex conjugage

        float hktImag  =  (-hokReal[index]   *  sinwt) + (hokImag[index] * coswt) +
                         (hokRealOpp        *  sinwt) + (hokImagOpp *  coswt);  //complex conjugage

        m_fft_htField[index].r = hktReal;
        m_fft_htField[index].i = hktImag;

        if(m_doChop)
        {
            float kMag = sqrt(kX * kX + kZ * kZ);
            float kXZ = kX * kZ / kMag;
            kX = kX / kMag;
            kZ = kZ / kMag;

            m_fft_chopX[index].r =  hktImag * kX;
            m_fft_chopX[index].i = -hktReal * kX;

            m_fft_chopZ[index].r =  hktImag * kZ;
            m_fft_chopZ[index].i = -hktReal * kZ;

            if(m_doFoam)
            {
                m_fft_jxx[index].r = hktReal * kX;
                m_fft_jxx[index].i = hktImag * kX;

                m_fft_jzz[index].r = hktReal * kZ;
                m_fft_jzz[index].i = hktImag * kZ;

                m_fft_jxz[index].r = hktReal * kXZ;
                m_fft_jxz[index].i = hktImag * kXZ;
            }
        }
    }

    size_t mem_needed = 0;
    kiss_fftnd_alloc(dims, 2, 1, NULL, &mem_needed); // Calculate the memory requirement without allocating
    m_memory += mem_needed;

    // Do the FFTs
    m_fft_plan = kiss_fftnd_alloc(dims, 2, 1, 0, 0);
    kiss_fftnd(m_fft_plan, m_fft_htField, m_fft_htField);
    if(m_doChop)
    {
        kiss_fftnd(m_fft_plan, m_fft_chopX, m_fft_chopX);
        kiss_fftnd(m_fft_plan, m_fft_chopZ, m_fft_chopZ);

        if(m_doFoam)
        {
            kiss_fftnd(m_fft_plan, m_fft_jxx, m_fft_jxx);
            kiss_fftnd(m_fft_plan, m_fft_jzz, m_fft_jzz);
            kiss_fftnd(m_fft_plan, m_fft_jxz, m_fft_jxz);
        }
    }

    // process the FFT results
    #pragma omp parallel for num_threads(m_max_threads)
    for (size_t index = 0; index < nn; ++index)
    {
        size_t i = index / n;  // Row index
        size_t j = index % n;  // Column index

        m_out_fft_htField[index] = m_fft_htField[index].r * signs[(i + j) & 1] * m_waveHeight;

        if(m_doChop)
        {
            float multiplier = -m_chopAmount * signs[(i + j) & 1];
            m_out_fft_chopX[index] = m_fft_chopX[index].r * multiplier;
            m_out_fft_chopZ[index] = m_fft_chopZ[index].r * multiplier;

            if(m_doFoam)
            {
                float Jxx = m_fft_jxx[index].r * multiplier + 1.0f;
                float Jzz = m_fft_jzz[index].r * multiplier + 1.0f;
                float Jxz = m_fft_jxz[index].r * multiplier;

                float temp = (0.5f * sqrt(((Jxx - Jzz) * (Jxx - Jzz)) + 4.0f * (Jxz*Jxz)));
                float jPlus = (0.5f * (Jxx + Jzz)) + temp;
                float jMinus = (0.5f * (Jxx + Jzz)) - temp;
                float qPlus = (jPlus - Jxx) / Jxz;
                float qMinus = (jMinus - Jxx) / Jxz;

                temp = sqrt(1.0f + qPlus * qPlus);
                m_out_fft_jxxX[index] = 1.0f / temp;
                m_out_fft_jxxZ[index] = qPlus / temp;

                temp = sqrt(1.0f + qMinus * qMinus);
                m_out_fft_jzzX[index] = 1.0f / temp;
                m_out_fft_jzzZ[index] = qMinus / temp;

                //store foam back in this array for convenience
                m_out_fft_jxz[index] = ((Jxx + Jzz) - sqrt(((Jxx - Jzz) * (Jxx - Jzz)) + 4.0f * (Jxz*Jxz))) * 0.5f;
            }
        }
    }

    bool clearAllArrays = false;
    clearShaderArrays(clearAllArrays);

    m_isAllocated = 1;
    snprintf(m_state, sizeof(m_state), "[aaOcean Shader] Generated ocean vector displacement at %dx%d resolution, spectrum: %d", m_resolution, m_resolution, m_spectrum);
    timer.printElapsed(m_state);
}

void aaOcean::clearShaderArrays(bool clearAll)
{
    if(m_fft_htField)
	    free(m_fft_htField);
    if(m_fft_chopX)
        free(m_fft_chopX);
    if(m_fft_chopZ)
        free(m_fft_chopZ);
    if(m_fft_plan)
        free(m_fft_plan);
    if(m_fft_jxx)
        free(m_fft_jxx);
    if(m_fft_jzz)
        free(m_fft_jzz);
    if(m_fft_jxz)
        free(m_fft_jxz);

    m_fft_htField = m_fft_chopX = m_fft_chopZ = m_fft_jxz = m_fft_jzz = m_fft_jxx = nullptr;
    m_fft_plan = nullptr;

    if(clearAll)
    {
        if(m_out_fft_htField)
            free(m_out_fft_htField);
        if(m_out_fft_chopX)
            free(m_out_fft_chopX);
        if(m_out_fft_chopZ)
            free(m_out_fft_chopZ);
        if(m_out_fft_jxxX)
            free(m_out_fft_jxxX);
        if(m_out_fft_jxxZ)
            free(m_out_fft_jxxZ);
        if(m_out_fft_jzzX)
            free(m_out_fft_jzzX);
        if(m_out_fft_jzzZ)
            free(m_out_fft_jzzZ);
        if(m_out_fft_jxz)
            free(m_out_fft_jxz);

        m_out_fft_htField = m_out_fft_chopX = m_out_fft_chopZ = m_out_fft_jxxX = m_out_fft_jxxZ = m_out_fft_jzzX = m_out_fft_jzzZ = m_out_fft_jxz = nullptr;
    }
}

void aaOcean::allocateBaseArrays()
{
    if (m_isAllocated)
    {
        snprintf(m_state, sizeof(m_state), "[aaOcean Core] Reallocating memory for ocean data structures for resolution %dx%d", m_resolution, m_resolution);
        clearArrays();
    }
    else
        snprintf(m_state, sizeof(m_state), "[aaOcean Core] Allocating memory for ocean data structures for resolution %dx%d", m_resolution, m_resolution);

    m_memory = 0;
    int size = m_resolution * m_resolution;
    int dims[2] = { m_resolution, m_resolution };

    m_fft_htField = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_fft_chopX = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_fft_chopZ = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_memory += size * sizeof(kiss_fft_cpx) * 3;

    m_xCoord = (int*)malloc(size * sizeof(int));
    m_zCoord = (int*)malloc(size * sizeof(int));
    m_memory += size * sizeof(int) * 2;

    m_hokReal = (float*)malloc(size * sizeof(float));
    m_hokImag = (float*)malloc(size * sizeof(float));
    m_hktReal = (float*)malloc(size * sizeof(float));
    m_hktImag = (float*)malloc(size * sizeof(float));
    m_kX = (float*)malloc(size * sizeof(float));
    m_kZ = (float*)malloc(size * sizeof(float));
    m_omega = (float*)malloc(size * sizeof(float));
    m_rand1 = (float*)malloc(size * sizeof(float));
    m_rand2 = (float*)malloc(size * sizeof(float));
    m_fftSpectrum = (float*)malloc(size * sizeof(float));
    m_memory += size * sizeof(float) * 10;

    m_out_fft_htField = (float*)malloc(size * sizeof(float));
    m_out_fft_chopX = (float*)malloc(size * sizeof(float));
    m_out_fft_chopZ = (float*)malloc(size * sizeof(float));
    m_memory += size * sizeof(float) * 3;

    size_t mem_needed = 0;
    kiss_fftnd_alloc(dims, 2, 1, NULL, &mem_needed); // Calculate the memory requirement without allocating
    m_memory += mem_needed;
    m_fft_plan = kiss_fftnd_alloc(dims, 2, 1, 0, 0);

    m_arrayPointer[eHEIGHTFIELD] = m_out_fft_htField;
    m_arrayPointer[eCHOPX] = m_out_fft_chopX;
    m_arrayPointer[eCHOPZ] = m_out_fft_chopZ;
    m_arrayPointer[eSPECTRUM] = m_fftSpectrum;

    m_doHoK = 1;
    m_doSetup = 1;
    m_isAllocated = 1;
}

void aaOcean::allocateFoamArrays()
{
    int size = m_resolution * m_resolution;

    m_fft_jxx = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_fft_jzz = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_fft_jxz = (kiss_fft_cpx *)malloc(size * sizeof(kiss_fft_cpx));
    m_memory += size * sizeof(kiss_fft_cpx) * 3;

    m_out_fft_jxxX = (float*)malloc(size * sizeof(float));
    m_out_fft_jxxZ = (float*)malloc(size * sizeof(float));
    m_out_fft_jzzX = (float*)malloc(size * sizeof(float));
    m_out_fft_jzzZ = (float*)malloc(size * sizeof(float));
    m_out_fft_jxz = (float*)malloc(size * sizeof(float));
    m_memory += size * sizeof(float) * 5;

    m_arrayPointer[eFOAM] = m_out_fft_jxz;
    m_arrayPointer[eEIGENPLUSX] = m_out_fft_jxxX;
    m_arrayPointer[eEIGENPLUSZ] = m_out_fft_jxxZ;
    m_arrayPointer[eEIGENMINUSX] = m_out_fft_jzzX;
    m_arrayPointer[eEIGENMINUSZ] = m_out_fft_jzzZ;

    m_isFoamAllocated = 1;
}

void aaOcean::clearResidualArrays()
{
    bool isResidualAllocated = 1;

    if (m_fftSpectrum)
    {
        free(m_fftSpectrum);
        m_fftSpectrum = 0;
    }
    else
        isResidualAllocated = 0;

    if (m_rand2)
    {
        free(m_rand2);
        m_rand2 = 0;
    }

    if (m_rand1)
    {
        free(m_rand1);
        m_rand1 = 0;
    }
    if (m_omega)
    {
        free(m_omega);
        m_omega = 0;
    }
    if (m_kZ)
    {
        free(m_kZ);
        m_kZ = 0;
    }
    if (m_kX)
    {
        free(m_kX);
        m_kX = 0;
    }
    if (m_hktImag)
    {
        free(m_hktImag);
        m_hktImag = 0;
    }
    if (m_hktReal)
    {
        free(m_hktReal);
        m_hktReal = 0;
    }
    if (m_hokImag)
    {
        free(m_hokImag);
        m_hokImag = 0;
    }
    if (m_hokReal)
    {
        free(m_hokReal);
        m_hokReal = 0;
    }
    if (m_zCoord)
    {
        free(m_zCoord);
        m_zCoord = 0;
    }
    if (m_xCoord)
    {
        free(m_xCoord);
        m_xCoord = 0;
    }

    if (isResidualAllocated)
    {
        int size = m_resolution * m_resolution;
        float cleared_memory = float(m_memory);
        m_memory = m_memory - (size * sizeof(float) * 9 + size * sizeof(int) * 2);
        cleared_memory = (cleared_memory - float(m_memory)) / 1048576.f;
        float working_memory = (float)m_memory / 1048576.f;

        char temp_buffer[sizeof(m_state)];
        strncpy(temp_buffer, m_state, sizeof(temp_buffer));
        snprintf(m_state, sizeof(m_state), "%s\n[aaOcean Core] Clearing %.2f megabytes of working memory. Current usage %.2f megabytes",
                temp_buffer, cleared_memory, working_memory);
    }
}

void aaOcean::clearArrays()
{
    if (m_isAllocated)
    {
        if (m_isFoamAllocated)
        {
            if (m_fft_jxx)
            {
                free(m_fft_jxx);
                free(m_out_fft_jxxX);
                free(m_out_fft_jxxZ);
                m_out_fft_jxxZ = m_out_fft_jxxX = 0;
                m_fft_jxx = 0;
            }
            if (m_fft_jzz)
            {
                free(m_fft_jzz);
                free(m_out_fft_jzzX);
                free(m_out_fft_jzzZ);
                m_out_fft_jzzX = m_out_fft_jzzZ = 0;
                m_fft_jzz = 0;
            }
            if (m_fft_jxz)
            {
                free(m_fft_jxz);
                free(m_out_fft_jxz);
                m_out_fft_jxz = 0;
                m_fft_jxz = 0;
            }
            m_isFoamAllocated = 0;
        }
        if (m_fft_chopZ)
        {
            free(m_fft_chopZ);
            free(m_out_fft_chopZ);
            m_out_fft_chopZ = 0;
            m_fft_chopZ = 0;
        }
        if (m_fft_chopX)
        {
            free(m_fft_chopX);
            free(m_out_fft_chopX);
            m_out_fft_chopX = 0;
            m_fft_chopX = 0;
        }
        if (m_fft_htField)
        {
            free(m_fft_htField);
            free(m_out_fft_htField);
            free(m_fft_plan);
            m_out_fft_htField = 0;
            m_fft_htField = 0;
        }
        m_isAllocated = 0;
    }
    clearResidualArrays();
}

u_int32_t aaOcean::generateUID(const float xCoord, const float zCoord) const
{
    // a very simple hash function. should probably do a better one at some point
    float angle;
    float length;
    float coordSq;
    float id_out;
    u_int32_t returnVal = 1;

    if (zCoord != 0.0f && xCoord != 0.0f)
    {
        coordSq = zCoord * zCoord + xCoord * xCoord;
        length = sqrt(coordSq);
        angle = xCoord / length;
        angle = acos(angle);
        angle = RadsToDegs(angle);

        if (angle > 180.0f)
            angle = 360.0f - angle;

        id_out = coordSq + (length * angle) + 0.5f;

        if (angle == 0.0f)
            returnVal = (u_int32_t)coordSq;
        else if (zCoord <= 0.0f)
            returnVal = (u_int32_t)floor(id_out);
        else
            returnVal = INT_MAX - (u_int32_t)floor(id_out);
    }
    return returnVal;
}

u_int32_t aaOcean::generateUIDHash(const float xCoord, const float zCoord) const
{
    // Return a default value for zero coordinates
    if (xCoord == 0.0f && zCoord == 0.0f)
        return 1;

    // Use integers to calculate a simple hash
    int ix = static_cast<int>(xCoord * 10000); // scale to avoid precision issues
    int iz = static_cast<int>(zCoord * 10000);

    // Generate a unique identifier using a simple but effective hashing technique
    u_int32_t hash = 2166136261u;
    hash = (hash ^ ix) * 16777619u;
    hash = (hash ^ iz) * 16777619u;

    return hash == 0 ? 1 : hash; // Ensure a non-zero return value
}

void aaOcean::setupGrid()
{
    Timer timer;

    const int n = m_resolution;
    const int half_n = (-n / 2) - ((n - 1) / 2);

    dsfmt_t dsfmt;
    #pragma omp parallel for firstprivate(dsfmt) num_threads(m_max_threads)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            unsigned int index, uID;

            index = (i * n) + j;

            m_xCoord[index] = half_n + i * 2;
            m_zCoord[index] = half_n + j * 2;
            float x = static_cast<float>(m_xCoord[index]);
            float z = static_cast<float>(m_zCoord[index]);
            uID = static_cast<u_int32_t>(generateUID(x, z));
            dsfmt_init_gen_rand(&dsfmt, uID + m_seed);

            float g1 = static_cast<float>(gaussian(dsfmt));
            float g2 = static_cast<float>(gaussian(dsfmt));

            if (m_randWeight > 0.0f)
            {
                float u1 = static_cast<float>(uniform(dsfmt));
                float u2 = static_cast<float>(uniform(dsfmt));

                m_rand1[index] = lerp(m_randWeight, g1, u1);
                m_rand2[index] = lerp(m_randWeight, g2, u2);
            }
            else
            {
                m_rand1[index] = g1;
                m_rand2[index] = g2;
            }
        }
    }
    timer.printElapsed("[aaOcean Core] Setup Grid Done");
    m_doSetup = 0;
}

float aaOcean::philips(float k_sq)
{
    const float L = m_velocity * m_velocity / aa_GRAVITY;
    return ((exp(-1.0f / (L * L * k_sq))) / (k_sq * k_sq)) * m_spectrumMult;
}

float aaOcean::piersonMoskowitz(float omega, float k_sq)
{
    float peakOmegaPM = 0.87f * aa_GRAVITY / (m_velocity / 2.f);
    
    const float alpha = 0.0081f;
    const float beta = 5.29f;

    float spectrum = 100.f * alpha * aa_GRAVITYSQ / pow(omega, 5.f);
    spectrum *= exp(-beta * pow(peakOmegaPM / omega, 4.f));

    return spectrum * m_spectrumMult * 5.f;
}

float aaOcean::tma(float omega, int index)
{
    // This function is based on the following paper
    // Empirical Directional Wave Spectra for Computer Graphics
    // Christopher J. Horvath

    float dimensionlessFetch = abs(aa_GRAVITY * m_jswpfetch / (m_velocity * m_velocity));
    float alpha = 0.076f * std::pow(dimensionlessFetch, -0.22f);
    
    float peakOmega = 2.84f * pow(aa_GRAVITY, 0.7f) * pow(m_jswpfetch, -0.3f) * pow(m_velocity, -0.4f);
    float sigma = (omega <= peakOmega) ? (0.07f) : (0.09f);
    float peakSharpening = pow(m_peakSharpening, exp(-pow((omega - peakOmega) / (sigma * peakOmega), 2.0f) / 2.0f));

    float tma = peakSharpening * (alpha * aa_GRAVITYSQ / std::pow(omega, 5.0f)) * \
                 std::exp(-1.25f * std::pow(peakOmega / omega, 4.0f));

    float m_kdGain = sqrt(m_oceanDepth / aa_GRAVITY);
    float wh = omega * m_kdGain;
    float kitaigorodskiiDepth =  0.5f + (0.5f * tanh(1.8f * (wh - 1.125f)));

    tma *= kitaigorodskiiDepth;
    return  tma * m_spectrumMult;
}

float aaOcean::swell(float spectrumValue, float kdotw, float k_mag)
{
    // TODO: implement frequency-dependent swell generation
    float freq = std::sqrt(k_mag * aa_GRAVITY) / (2.f * aa_PI);

    // Low-frequency threshold for identifying swells
    const float swellThreshold = 0.2f;  // Adjust this threshold as needed

    if (freq < swellThreshold) 
        return spectrumValue * (1.0f + m_swell);

    return spectrumValue;
}

void aaOcean::evaluateHokData()
{
    Timer timer;

    const int      n = m_resolution * m_resolution;
    const float    k_mult = aa_TWOPI / m_oceanScale;
    const float    windx = cos(m_windDir);
    const float    windz = sin(m_windDir);
    const float    omega0 = aa_TWOPI / m_loopTime;

    bool bDamp = 0;
    if (m_damp > 0.0f)
        bDamp = 1;

    #pragma omp parallel for 
    for (int index = 0; index < n; ++index)
    {
        float k_sq, k_mag, k_dot_w;
        m_kX[index] = (float)m_xCoord[index] * k_mult;
        m_kZ[index] = (float)m_zCoord[index] * k_mult;
        k_sq = (m_kX[index] * m_kX[index]) + (m_kZ[index] * m_kZ[index]);
        k_mag = sqrt(k_sq);
        if (k_mag < aa_EPSILON)
                m_fftSpectrum[index] = 0.f;
        
        // calculate spectrum
        if (m_spectrum == 1) // Pierson-Moskowitz
        {
            // build dispersion relationship with oceanDepth relationship and capillary waves
            m_omega[index] = aa_GRAVITY * k_mag * tanh(k_mag * m_oceanDepth);
            m_omega[index] = sqrt(m_omega[index]);
            m_fftSpectrum[index] = piersonMoskowitz(m_omega[index], k_sq);

        }
        else if (m_spectrum == 2) //  Texel MARSEN ARSLOE (TMA) SpectruM
        {
            // build dispersion relationship with oceanDepth relationship and capillary waves
            m_omega[index] = aa_GRAVITY * k_mag * tanh(k_mag * m_oceanDepth);
            m_omega[index] = sqrt(m_omega[index]);
            m_fftSpectrum[index] = tma(m_omega[index], index);
        }
        else if (m_spectrum == 3) //  JONSWAP
        {
            // Constants for the simulation
            const float g = aa_GRAVITY;                     // Gravitational acceleration (m/s^2)
            const float alpha = 0.0081f;                    // Empirical constant for wave growth
            const float gamma = m_peakSharpening * 3.3f;    // Peak enhancement factor
            const float sigma_low = 0.07f;                  // Sigma for f < f_peak
            const float sigma_high = 0.09f;                 // Sigma for f > f_peak
            const float wind_speed = m_velocity * 10.f;     // Wind speed at 10 meters above sea level (m/s)
            const float fetch = m_jswpfetch * 100.f;        // Fetch distance in meters

            const float f_peak = (3.5f * std::pow(g, 0.67f)) / (std::pow(wind_speed, 0.34f) * std::pow(fetch, 0.33f));
            const float k_peak = 4.0f * aa_PI * aa_PI * f_peak * f_peak / g;
            float sigma = (k_mag <= k_peak) ? sigma_low : sigma_high;
            float S = alpha * g * g / std::pow(k_mag, 5.0f) * std::exp(-5.0f / 4.0f * std::pow(k_mag / k_peak, -4.0f));
            float peak_factor = std::pow(gamma, std::exp(-std::pow((k_mag - k_peak), 2) / (2.0f * sigma * sigma * k_peak * k_peak)));

            m_fftSpectrum[index] = std::sqrt(S * peak_factor / 2.0f) * m_spectrumMult * 0.01f;
        }
        else // Philips
        { 
            // modifying dispersion for capillary waves
            m_omega[index] = aa_GRAVITY * k_mag * (1.0f + k_sq * m_surfaceTension * m_surfaceTension);
            m_omega[index] = sqrt(m_omega[index]);
            // add time looping support with OmegaNought
            m_omega[index] = (int(m_omega[index] / omega0)) * omega0;

            m_fftSpectrum[index] = sqrt(philips(k_sq));
        }

        // bias towards wind dir
        float k_mag_inv = 1.0f / k_mag;
        k_dot_w = (m_kX[index] * k_mag_inv * windx) + (m_kZ[index] * k_mag_inv * windz);
        m_fftSpectrum[index] *= pow(k_dot_w, m_windAlign);  

        // eliminate wavelengths smaller than cutoff
        m_fftSpectrum[index] *= exp(-k_sq * m_cutoff);     

        // swell
        if (m_spectrum == 2)
            m_fftSpectrum[index] = swell(m_fftSpectrum[index], k_dot_w, k_mag); 

       // reduce reflected waves
        if (bDamp && (k_dot_w < 0.0f))
            m_fftSpectrum[index] *= (1.0f - m_damp);

        m_hokReal[index] = aa_INV_SQRTTWO * m_rand1[index] * m_fftSpectrum[index];
        m_hokImag[index] = aa_INV_SQRTTWO * m_rand2[index] * m_fftSpectrum[index];
    }

    snprintf(m_state, sizeof(m_state), "\n[aaOcean Core] Finished initializing all ocean data");

    timer.printElapsed("[aaOcean Core] Evaluated HoK Data");
    m_doHoK = 0;
}

void aaOcean::evaluateHieghtField()
{
    Timer timer;

    const float wt = m_waveSpeed * m_time;
    const int n = m_resolution;
    const int nn = n * n;
    const int n_sq = n * n - 1;
    const float signs[2] = { 1.0f, -1.0f };

    #pragma omp parallel for
    for (size_t index = 0; index < nn; ++index)
    {
        float hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt;

        size_t index_rev = n_sq - index; //tail end 
        hokReal = m_hokReal[index];
        hokImag = m_hokImag[index];
        hokRealOpp = m_hokReal[index_rev];
        hokImagOpp = m_hokImag[index_rev];

        coswt = cos(m_omega[index] * wt);
        sinwt = sin(m_omega[index] * wt);

        m_hktReal[index] =  (hokReal    *  coswt) + (hokImag    *  sinwt) +
                            (hokRealOpp *  coswt) - (hokImagOpp *  sinwt);  //complex conjugage

        m_hktImag[index] =  (-hokReal    *  sinwt) + (hokImag    *  coswt) +
                            (hokRealOpp *  sinwt) + (hokImagOpp *  coswt);  //complex conjugage

        m_fft_htField[index].r = m_hktReal[index];
        m_fft_htField[index].i = m_hktImag[index];
    }
    
    timer.printElapsed("[aaOcean Core] Heightfield data prepared");

    kiss_fftnd(m_fft_plan, m_fft_htField, m_fft_htField);

    timer.printElapsed("[aaOcean Core] Heightfield FFT done");

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
            m_out_fft_htField[(i*n) + j] = m_fft_htField[(i*n) + j].r * signs[(i + j) & 1] * m_waveHeight;
    }
    timer.printElapsed("[aaOcean Core] Heightfield FFT copied with signs swapped");
}

void aaOcean::evaluateChopField()
{
    int n = m_resolution * m_resolution;
    const float signs[2] = { 1.0f, -1.0f };

    Timer timer;

    for (size_t index = 0; index < n; ++index)
    {
        float  kX, kZ, kMag;
        kMag = sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
        kX = m_kX[index] / kMag;
        kZ = m_kZ[index] / kMag;

        m_fft_chopX[index].r = m_hktImag[index] * kX;
        m_fft_chopX[index].i = -m_hktReal[index] * kX;

        m_fft_chopZ[index].r = m_hktImag[index] * kZ;
        m_fft_chopZ[index].i = -m_hktReal[index] * kZ;
    }

    timer.printElapsed("[aaOcean Core] chopfield prepared");

    kiss_fftnd(m_fft_plan, m_fft_chopX, m_fft_chopX);
    kiss_fftnd(m_fft_plan, m_fft_chopZ, m_fft_chopZ);

    timer.printElapsed("[aaOcean Core] chopfield fft done");

    n = m_resolution;
    for (size_t i = 0; i < n; ++i)
    {
        int index;
        float multiplier;
        for (size_t j = 0; j < n; ++j)
        {
            index = (i*n) + j;
            multiplier = m_chopAmount * signs[(i + j) & 1] * -1.0f;
            m_out_fft_chopX[index] = m_fft_chopX[index].r * multiplier;
            m_out_fft_chopZ[index] = m_fft_chopZ[index].r * multiplier;
        }
    }
    timer.printElapsed("[aaOcean Core] chopfield copied with signs swapped");
}

void aaOcean::evaluateJacobians()
{
    const float signs[2] = { 1.0f, -1.0f };
    int n = m_resolution * m_resolution;

    Timer timer;

    #pragma omp parallel for
    for (size_t index = 0; index < n; ++index)
    {
        float kX, kZ, kMag, kXZ;
        kMag = 1.0f / sqrt(m_kX[index] * m_kX[index] + m_kZ[index] * m_kZ[index]);
        kX = (m_kX[index] * m_kX[index]) * kMag;
        kZ = (m_kZ[index] * m_kZ[index]) * kMag;
        kXZ = (m_kX[index] * m_kZ[index]) * kMag;

        m_fft_jxx[index].r = m_hktReal[index] * kX;
        m_fft_jxx[index].i = m_hktImag[index] * kX;

        m_fft_jzz[index].r = m_hktReal[index] * kZ;
        m_fft_jzz[index].i = m_hktImag[index] * kZ;

        m_fft_jxz[index].r = m_hktReal[index] * kXZ;
        m_fft_jxz[index].i = m_hktImag[index] * kXZ;
    }
    timer.printElapsed("[aaOcean Core] Jacobians setup complete");

    kiss_fftnd(m_fft_plan, m_fft_jxx, m_fft_jxx);
    kiss_fftnd(m_fft_plan, m_fft_jzz, m_fft_jzz);
    kiss_fftnd(m_fft_plan, m_fft_jxz, m_fft_jxz);

    timer.printElapsed("[aaOcean Core] Jacobians FFT done");

    n = m_resolution;
    #pragma omp parallel for 
    for (size_t i = 0; i < n; ++i)
    {
        float multiplier;
        size_t index;
        for (size_t j = 0; j < n; ++j)
        {
            index = (i*n) + j;
            multiplier = -m_chopAmount * signs[(i + j) & 1];
            m_fft_jxx[index].r = m_fft_jxx[index].r * multiplier + 1.0f;
            m_fft_jzz[index].r = m_fft_jzz[index].r * multiplier + 1.0f;
            m_fft_jxz[index].r = m_fft_jxz[index].r * multiplier;
        }
    }

    timer.printElapsed("[aaOcean Core] Jacobians FFT copied");

    #pragma omp parallel for 
    for (size_t index = 0; index < n*n; ++index)
    {
        float jPlus, jMinus, qPlus, qMinus, Jxx, Jzz, Jxz, temp;

        Jxx = m_fft_jxx[index].r;
        Jzz = m_fft_jzz[index].r;
        Jxz = m_fft_jxz[index].r;

        temp = (0.5f * sqrt(((Jxx - Jzz) * (Jxx - Jzz)) + 4.0f * (Jxz*Jxz)));
        jPlus = (0.5f * (Jxx + Jzz)) + temp;
        jMinus = (0.5f * (Jxx + Jzz)) - temp;
        qPlus = (jPlus - Jxx) / Jxz;
        qMinus = (jMinus - Jxx) / Jxz;

        temp = sqrt(1.0f + qPlus * qPlus);
        m_out_fft_jxxX[index] = 1.0f / temp;
        m_out_fft_jxxZ[index] = qPlus / temp;

        temp = sqrt(1.0f + qMinus * qMinus);
        m_out_fft_jzzX[index] = 1.0f / temp;
        m_out_fft_jzzZ[index] = qMinus / temp;

        //store foam back in this array for convenience
        m_out_fft_jxz[index] = ((Jxx + Jzz) - sqrt(((Jxx - Jzz) * (Jxx - Jzz)) + 4.0f * (Jxz*Jxz))) * 0.5f;
    }

    timer.printElapsed("[aaOcean Core] Jacobians Foam completed");
}

void aaOcean::getFoamBounds(float& outBoundsMin, float& outBoundsMax) const
{
    outBoundsMax = -FLT_MAX;
    outBoundsMin = FLT_MAX;
    Timer timer;
    const size_t n = m_resolution * m_resolution;
    for (size_t index = 0; index < n; index++)
    {
        if (outBoundsMax < m_out_fft_jxz[index])
            outBoundsMax = m_out_fft_jxz[index];

        if (outBoundsMin > m_out_fft_jxz[index])
            outBoundsMin = m_out_fft_jxz[index];
    }

    timer.printElapsed("[aaOcean Core] Foam Bounds calculated");
}

void aaOcean::getOceanArray(float *&outArray, aaOcean::arrayType type) const
{
    // get the pointer to the aaOcean array that we want to pull data from
    Timer timer;

    float *arrayPointer = m_arrayPointer[type];
    const int arraySize = m_resolution * m_resolution;

    for (int i = 0; i < arraySize; ++i)
        outArray[i] = arrayPointer[i];

    timer.printElapsed("[aaOcean Core] got Ocean Arrays");
}

float aaOcean::getOceanData(float uCoord, float vCoord, aaOcean::arrayType type) const
{
    float u, v, du, dv = 0;
    int xMinus1, yMinus1, x, y, xPlus1, yPlus1, xPlus2, yPlus2;

    // get the pointer to the aaOcean array that we want to pull data from
    float *arrayPointer = m_arrayPointer[type];

    // begin UV coordinate wrapping to [0-1] interval
    uCoord = fmod(uCoord, 1.0f);
    vCoord = fmod(vCoord, 1.0f);
    if(uCoord < 0.0f)
        uCoord = (fabs(uCoord) >= 1e-6) ? 1.0f + uCoord : 0.0f;
    if(vCoord < 0.0f)
        vCoord = (fabs(vCoord) >= 1e-6) ? 1.0f + vCoord : 0.0f;

    // use UV coordinates to work out ocean array indeces
    u = uCoord * static_cast<float>(m_resolution);
    v = vCoord * static_cast<float>(m_resolution);
    x = static_cast<int>(floor(u));
    y = static_cast<int>(floor(v));
    du = u - x; 
    dv = v - y; 

    // prepare catmul-rom end points for interpolation
    // wrap any indices that are outside the array boundaries
    xMinus1 = wrap(x - 1) * m_resolution;
    xPlus1  = wrap(x + 1) * m_resolution;
    xPlus2  = wrap(x + 2) * m_resolution;
    x       = wrap(x)   * m_resolution;
    y       = wrap(y);
    yMinus1 = wrap(y - 1);
    yPlus1  = wrap(y + 1);
    yPlus2  = wrap(y + 2);
    
    // prepare for catmul-rom interpolation
    const float a1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yMinus1],
                            arrayPointer[x          + yMinus1],
                            arrayPointer[xPlus1     + yMinus1],
                            arrayPointer[xPlus2     + yMinus1]);

    const float b1 = catmullRom(du, 
                            arrayPointer[xMinus1    +   y],
                            arrayPointer[x          +   y],
                            arrayPointer[xPlus1     +   y],
                            arrayPointer[xPlus2     +   y]);

    const float c1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yPlus1],
                            arrayPointer[x          + yPlus1],
                            arrayPointer[xPlus1     + yPlus1],
                            arrayPointer[xPlus2     + yPlus1]);

    const float d1 = catmullRom(du, 
                            arrayPointer[xMinus1    + yPlus2],
                            arrayPointer[x          + yPlus2],
                            arrayPointer[xPlus1     + yPlus2],
                            arrayPointer[xPlus2     + yPlus2]);

    return catmullRom(dv, a1, b1, c1, d1);
}

inline float aaOcean::catmullRom(const float t, const float a, const float b, const float c, const float d) const
{
    return  0.5f * ( ( 2.0f * b ) + ( -a + c ) * t + 
            ( 2.0f * a - 5.0f * b + 4.0f * c - d ) * t * t + 
            ( -a + 3.0f * b - 3.0f * c + d )* t * t * t );
}

inline int aaOcean::wrap(int x) const
{
    // fast modulo for power of 2 numbers
    x = x & (m_resolution - 1);
    return x;
}

size_t aaOcean::getMemory()
{
    return m_memory;
}

#endif  /* AAOCEANCLASS_CPP */

