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
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include "dSFMT/dSFMT.c"
#include "constants.h"
#include "functionLib.h"
#include "timer.h"

class aaSpectrum
{
public:

    aaSpectrum()
    {
        m_hokReal = m_hokImag = m_hktReal = m_hktImag = m_omega = nullptr;
    }

    ~aaSpectrum()
    {
        if(m_hokReal)
            free(m_hokReal);
        if(m_hokImag)
            free(m_hokImag);
        if(m_hktReal)
            free(m_hktReal);
        if(m_hktImag)
            free(m_hktImag);
        if(m_omega)
            free(m_omega);
    }

    void allocate(size_t alignment = 16)
    {
        size_t size_flt_array = m_resolution * m_resolution * sizeof(float*);
        m_hokReal = static_cast<float*>(aligned_alloc(alignment, size_flt_array ));
        m_hokImag = static_cast<float*>(aligned_alloc(alignment, size_flt_array ));
        m_hktReal = static_cast<float*>(aligned_alloc(alignment, size_flt_array ));
        m_hktImag = static_cast<float*>(aligned_alloc(alignment, size_flt_array ));
        m_omega   = static_cast<float*>(aligned_alloc(alignment, size_flt_array ));
    }

    // main input function
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
        m_resolution = static_cast<int>(powf(2.0f, (4 + abs(resolution))));

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

        m_chopAmount    = std::max(chopAmount, aa_EPSILON);
        oceanScale      = std::max(oceanScale, aa_EPSILON);
        velocity        = std::max(velocity,   aa_EPSILON);
        oceanDepth      = std::max(oceanDepth, aa_EPSILON);
        cutoff          = std::fabs(cutoff * 0.01f);              // scaling for better UI control
        windDir         = windDir * aa_PIBY180;                   // to radians
        windAlign       = std::max<int>((windAlign + 1) & ~1, 0); // forcing to even numbers
        damp            = std::clamp(damp, 0.0f, 1.0f);
        
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

    void BuildSpectrum()
    {
        Timer timer;
        allocate();

        const int n = m_resolution;
        const int nn = n * n;
        const int n_sq = n * n - 1;

        const int half_n = (-n / 2) - ((n - 1) / 2);
        const int half_n_plus_1 = n/2 + 1;

        const float k_mult = aa_TWOPI / m_oceanScale;
        const float windx = cosf(m_windDir);
        const float windz = sinf(m_windDir);
        const float omega0 = aa_TWOPI / m_loopTime;
        const float wt = m_waveSpeed * m_time;
        
        bool bDamp = 0;
        if (m_damp > aa_EPSILON)
            bDamp = 1;

        dsfmt_t dsfmt;
        #pragma omp parallel for firstprivate(dsfmt, bDamp)
        for (int i = 0; i < n; ++i)
        {
            int xCoord, zCoord;
            unsigned int index, uID;
            float rand1, rand2;
            float kX, kZ, k_sq, k_mag, k_dot_w, omega, spectrum;
            float hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt;

            for (int j = 0; j < n; ++j)
            {
                index = (i * n) + j;

                // build random numbers
                xCoord = half_n + i * 2;
                zCoord = half_n + j * 2;
                float x = static_cast<float>(xCoord);
                float z = static_cast<float>(zCoord);
                uID = static_cast<u_int32_t>(GenerateUID(x, z));
                dsfmt_init_gen_rand(&dsfmt, uID + m_seed);
                rand1 = static_cast<float>(gaussian(dsfmt));
                rand2 = static_cast<float>(gaussian(dsfmt));
                if (m_randWeight > aa_EPSILON)
                {
                    float u1 = static_cast<float>(uniform(dsfmt));
                    float u2 = static_cast<float>(uniform(dsfmt));
                    rand1 = lerp(m_randWeight, rand1, u1);
                    rand2 = lerp(m_randWeight, rand2, u2);
                }

                // build spectrum
                kX = static_cast<float>(xCoord) * k_mult;
                kZ = static_cast<float>(zCoord) * k_mult;
                k_sq = kX *kX + kZ * kZ;
                k_mag = sqrtf(k_sq);
                float k_mag_inv = 1.0f / k_mag;

                // build dispersion relationship with oceanDepth relationship and capillary waves
                omega = aa_GRAVITY * k_mag * tanhf(k_mag * m_oceanDepth);

                // calculate spectrum
                if (m_spectrum == 1) // Pierson-Moskowitz
                {
                    omega = sqrtf(omega);
                    spectrum = spPiersonMoskowitz(omega, k_sq);
                }
                else if (m_spectrum == 2) //  Texel MARSEN ARSLOE (TMA) SpectruM
                {
                    omega = sqrtf(omega);
                    spectrum = spTMA(omega, index);
                }
                else // Philips
                { 
                    // modifying dispersion for capillary waves
                    omega = aa_GRAVITY * k_mag * (1.0f + k_sq * m_surfaceTension * m_surfaceTension);
                    omega = sqrt(omega);
                    omega = (int(omega / omega0)) * omega0;
                    spectrum = sqrtf(spPhilips(k_sq));
                }

                // spectrum indenpendant modifications
                k_dot_w = (kX * k_mag_inv * windx) + (kZ * k_mag_inv * windz);
                spectrum *= powf(k_dot_w, m_windAlign);  // bias towards wind dir
                spectrum *= expf(-k_sq * m_cutoff);      // eliminate wavelengths smaller than cutoff

                // reduce reflected waves
                if (bDamp && (k_dot_w < 0.0f))
                    spectrum *= (1.0f - m_damp);

                m_hokReal[index] = aa_INV_SQRTTWO * rand1 * spectrum;
                m_hokImag[index] = aa_INV_SQRTTWO * rand2 * spectrum;
                m_omega[index] = omega;
            }
        }

        #pragma omp parallel for
        for (size_t index = 0; index < nn; ++index)
        {
            float hokReal, hokImag, hokRealOpp, hokImagOpp, sinwt, coswt;

            size_t index_rev = n_sq - index; //tail end 
            hokReal = m_hokReal[index];
            hokImag = m_hokImag[index];
            hokRealOpp = m_hokReal[index_rev];
            hokImagOpp = m_hokImag[index_rev];

            coswt = cosf(m_omega[index] * wt);
            sinwt = sinf(m_omega[index] * wt);

            m_hktReal[index] =  (hokReal    *  coswt) + (hokImag    *  sinwt) +
                                (hokRealOpp *  coswt) - (hokImagOpp *  sinwt);  //complex conjugage

            m_hktImag[index] =  (-hokReal    *  sinwt) + (hokImag    *  coswt) +
                                (hokRealOpp *  sinwt) + (hokImagOpp *  coswt);  //complex conjugage
        }

        char message[256];
        sprintf(message, "[aaOcean Spectrum] Spectrum done for [%d x %d] ocean", m_resolution, m_resolution);
        timer.printElapsed(message, true);
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
        // The Phillips Spectrum serves as a theoretical baseline for fully developed waves

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
        // for shallow seas

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
        // JONSWAP Spectrum (Joint North Sea Wave Project)
        // JONSWAP is a more realistic model for young, growing seas, capturing the peaked energy distribution observed in nature, 
        // while Phillips is a more simplified, theoretical approach for fully developed seas.
        // TODO: implement JONSWAP
        return 1.f;
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

    float *m_hokReal; // real component of HoK (see Tessendorf paper)
    float *m_hokImag; // imaginary component of HoK (see Tessendorf paper)
    float *m_hktReal; // real component of HkT (see Tessendorf paper)
    float *m_hktImag; // real component of HkT (see Tessendorf paper)
    float *m_omega;   // omega (see Tessendorf paper)
};
#endif