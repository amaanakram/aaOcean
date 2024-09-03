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

#ifndef AAOCEANSPECTRUM_H
#define AAOCEANSPECTRUM_H

#include <vector>
#include <cmath>
#include <algorithm> 
#include <limits.h>

#include "dSFMT/dSFMT.c"
#include "constants.h"
#include "functionLib.h"
#include "timer.h"

class aaSpectrum
{
public:
    aaSpectrum();
    ~aaSpectrum();

    // main input function
    // should be called from SpectrumEvaluate
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
        float   loopTime,
        bool    doFoam,
        float   randWeight = 0.f,
        float   spectrumMult = 1.f,
        float   peakSharpening = 1.f,
        float   jswpfetch = 100.f,
        float   swell = 0.0f)
    {
        resolution = static_cast<int>(powf(2.0f, (4 + abs(resolution))));
        if (m_resolution != resolution)
            m_resolution = resolution;

        // scaled for better UI control
        if (spectrum < 2 )
        {
            // for Philips and PM spectrums
            waveHeight *= 0.01f;
            chopAmount *= 0.01f;
        }
        else
        {
            // TMA spectrum
            waveHeight *= 0.1f;
            chopAmount *= 0.1f;
        }

        m_waveHeight    = waveHeight;
        m_time          = time;
        m_waveSpeed     = waveSpeed;

        if (chopAmount > aa_FLT_EPSILON)
            m_chopAmount = chopAmount;
        else
            m_chopAmount = 0.0f;

        // clamping to minimum value
        oceanScale  = std::max<float>(oceanScale, aa_FLT_EPSILON);
        velocity    = std::max<float>(velocity, aa_FLT_EPSILON);
        oceanDepth  = std::max<float>(oceanDepth, aa_FLT_EPSILON);
        // scaling for better UI control
        cutoff = fabsf(cutoff * 0.01f);
        // to radians
        windDir = windDir * aa_PIBY180;
        // forcing to even numbers
        windAlign = std::max<int>(((windAlign + 1) * 2), 2);
        // clamping to a maximum value of 1
        damp = std::min<float>(damp, 1.0f);
        
        m_seed = seed;
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
        m_randWeight = randWeight;
    }

    void SetupGrid()
    {
        Timer timer;
        const int n = m_resolution;
        const int half_n = (-n / 2) - ((n - 1) / 2);
        const int half_n_plus_1 = n/2 + 1;

        m_xCoord.resize(n * half_n_plus_1);
        m_zCoord.resize(n * half_n_plus_1);
        m_rand1.resize(n * half_n_plus_1);
        m_rand2.resize(n * half_n_plus_1);

        dsfmt_t dsfmt;
        #pragma omp parallel for firstprivate(dsfmt)
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
                uID = static_cast<u_int32_t>(GenerateUID(x, z));
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
        timer.printElapsed("[aaOcean Spectrum] Setup Grid Done");
    }

    u_int32_t GenerateUID(const float xCoord, const float zCoord) const
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
            length = sqrtf(coordSq);
            angle = xCoord / length;
            angle = acosf(angle);
            angle = RadsToDegs(angle);

            if (angle > 180.0f)
                angle = 360.0f - angle;

            id_out = coordSq + (length * angle) + 0.5f;

            if (isFlEqual(angle, 0.0f))
                returnVal = (u_int32_t)coordSq;
            else if (zCoord <= 0.0f)
                returnVal = (u_int32_t)floor(id_out);
            else
                returnVal = INT_MAX - (u_int32_t)floor(id_out);
        }
        return returnVal;
    }

    float spPhilips(float k_sq)
    {
        const float L = m_velocity * m_velocity / aa_GRAVITY;
        return ((expf(-1.0f / (L * L * k_sq))) / (k_sq * k_sq)) * m_spectrumMult;
    }

    float spPiersonMoskowitz(float omega, float k_sq)
    {
        // https://wikiwaves.org/Ocean-Wave_Spectra

        float peakOmegaPM = 0.87f * aa_GRAVITY / m_velocity;
        
        const float alpha = 0.0081f;
        const float beta = 1.29f;

        float spectrum = 100.f * alpha * aa_GRAVITYSQ / powf(omega, 5.f);
        spectrum *= expf(-beta * powf(peakOmegaPM / omega, 4.f));

        return spectrum * m_spectrumMult;
    }

    float spTMA(float omega, int index)
    {
        // This function is based on the following paper
        // Empirical Directional Wave Spectra for Computer Graphics
        // Christopher J. Horvath

        float dimensionlessFetch = abs(aa_GRAVITY * m_jswpfetch / (m_velocity * m_velocity));
        float alpha = 0.076f * powf(dimensionlessFetch, -0.22f);
        
        float peakOmega = 2.84f * powf(aa_GRAVITY, 0.7f) * powf(m_jswpfetch, -0.3f) * powf(m_velocity, -0.4f);
        float sigma = (omega <= peakOmega) ? (0.07f) : (0.09f);
        float peakSharpening = powf(m_peakSharpening, expf(-powf((omega - peakOmega) / (sigma * peakOmega), 2.0f) / 2.0f));

        float tma = peakSharpening * (alpha * aa_GRAVITYSQ / std::pow(omega, 5.0f)) * \
                    expf(-1.25f * powf(peakOmega / omega, 4.0f));

        float m_kdGain = sqrt(m_oceanDepth / aa_GRAVITY);
        float wh = omega * m_kdGain;
        float kitaigorodskiiDepth =  0.5f + (0.5f * tanhf(1.8f * (wh - 1.125f)));

        tma *= kitaigorodskiiDepth;
        return  tma * m_spectrumMult;
    }

    float swellTMA(float omega, float kdotw, float k_mag)
    {
        // TODO: implement frequency-dependent swell generation
        return 1.f;
    }

    float spJONSWAP()
    {
        // https://wikiwaves.org/Ocean-Wave_Spectra
    }

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
    float   m_randWeight;     // control blend between rand distributions

    // optional variables
    float   m_spectrumMult;     // multiplier for generated spectrum
    float   m_peakSharpening;   // JONSWAP Peak Sharpening
    float   m_jswpfetch;        // wind region
    float   m_swell;            // swell

    std::vector<int>    m_xCoord;  // holds ocean grid coordinates
    std::vector<int>    m_zCoord;  // holds ocean grid coordinates
    std::vector<float>  m_hokReal; // real component of HoK (see Tessendorf paper)
    std::vector<float>  m_hokImag; // imaginary component of HoK (see Tessendorf paper)
    std::vector<float>  m_hktReal; // real component of HkT (see Tessendorf paper)
    std::vector<float>  m_hktImag; // real component of HkT (see Tessendorf paper)
    std::vector<float>  m_kX;      // x-component of wave vector
    std::vector<float>  m_kZ;      // z-component of wave vector
    std::vector<float>  m_omega;   // omega (see Tessendorf paper)
    std::vector<float>  m_rand1;   // random number array 
    std::vector<float>  m_rand2;   // random number array 
};
#endif