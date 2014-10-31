#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "FilterClass.h"
#include <iostream>

using namespace std;

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/LowPassFilter"
#define N_SAMPLES_DEFAULT 64
enum {IN, OUT_1, FREQ, ORDER, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class LowPassFilter
{
public:
    LowPassFilter(double samplerate, uint32_t n_samples){Construct(samplerate, n_samples);}
    ~LowPassFilter(){Destruct();}
    void Construct(double samplerate, uint32_t n_samples)
    {
        SampleRate = samplerate;
        lpf = new FilterClass(samplerate, n_samples);
    }
    void Destruct()
    {
        delete lpf;
    }
    void Realloc(uint32_t n_samples)
    {
        Destruct();
        Construct(SampleRate, n_samples);
    }

    static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features);
    static void activate(LV2_Handle instance);
    static void deactivate(LV2_Handle instance);
    static void connect_port(LV2_Handle instance, uint32_t port, void *data);
    static void run(LV2_Handle instance, uint32_t n_samples);
    static void cleanup(LV2_Handle instance);
    static const void* extension_data(const char* uri);
    float *in;
    float *out_1;
    float *freq;
    float *order;

    double SampleRate;
    FilterClass *lpf;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    LowPassFilter::instantiate,
    LowPassFilter::connect_port,
    LowPassFilter::activate,
    LowPassFilter::run,
    LowPassFilter::deactivate,
    LowPassFilter::cleanup,
    LowPassFilter::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle LowPassFilter::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    LowPassFilter *plugin = new LowPassFilter(samplerate, N_SAMPLES_DEFAULT);
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void LowPassFilter::activate(LV2_Handle instance){}

/**********************************************************************************************************************************************************/

void LowPassFilter::deactivate(LV2_Handle instance){}

/**********************************************************************************************************************************************************/

void LowPassFilter::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    LowPassFilter *plugin;
    plugin = (LowPassFilter *) instance;

    switch (port)
    {
        case IN:
            plugin->in = (float*) data;
            break;
        case OUT_1:
            plugin->out_1 = (float*) data;
            break;
        case FREQ:
            plugin->freq = (float*) data;
            break;
        case ORDER:
            plugin->order = (float*) data;
            break;
    }
}

/**********************************************************************************************************************************************************/

void LowPassFilter::run(LV2_Handle instance, uint32_t n_samples)
{
    LowPassFilter *plugin;
    plugin = (LowPassFilter *) instance;

    float soma_abs = 0;
    for (uint32_t i = 0; i < n_samples; i++) soma_abs += abs(plugin->in[i]);

    if (soma_abs == 0)
        fill_n(plugin->out_1, n_samples, 0);
    else
    {
        double f = (double)(*(plugin->freq));
        float Order = (float)(*(plugin->order));
        Order = round(Order)+1;

        if ( (plugin->lpf)->N != (int)n_samples )
            plugin->Realloc(n_samples);
        
        for (uint32_t i=0; i < n_samples;i++) (plugin->lpf)->u[i] = plugin->in[i];

            switch ((int)Order)
            {
                case 1:
                    (plugin->lpf)->LPF1_Bilinear(f);
                    break;
                case 2:
                    (plugin->lpf)->LPF2_Bilinear(f);
                    break;
                case 3:
                    (plugin->lpf)->LPF3_Bilinear(f);
                    break;
            }

        for (uint32_t i=0; i < n_samples;i++) plugin->out_1[i] = (plugin->lpf)->y[i];
    }
}

/**********************************************************************************************************************************************************/

void LowPassFilter::cleanup(LV2_Handle instance)
{
    delete (LowPassFilter *) instance;
}

/**********************************************************************************************************************************************************/

const void* LowPassFilter::extension_data(const char* uri)
{
    return NULL;
}