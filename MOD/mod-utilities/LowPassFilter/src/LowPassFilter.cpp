#include <stdlib.h>
#include <cmath>
#include <lv2.h>
#include "FilterClass.h"
#include <algorithm>
#include "SomeUtilities.h"

using namespace std;

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/LowPassFilter"
#define N_SAMPLES_DEFAULT 64
enum {IN, OUT, FREQ, ORDER, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class LowPassFilter
{
public:
    LowPassFilter(double samplerate, uint32_t n_samples){Construct(samplerate, n_samples);}
    ~LowPassFilter(){Destruct();}
    void Construct(double samplerate, uint32_t n_samples)
    {
    	this->n_samples = n_samples;
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
    float *ports[PLUGIN_PORT_COUNT];

    uint32_t n_samples;
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
    plugin->ports[port] = (float*) data;
}

/**********************************************************************************************************************************************************/

void LowPassFilter::run(LV2_Handle instance, uint32_t n_samples)
{
    LowPassFilter *plugin;
    plugin = (LowPassFilter *) instance;

    float *in   = plugin->ports[IN];
    float *out  = plugin->ports[OUT];
    double f    = (double)(*(plugin->ports[FREQ]));
    int Order   = (int)(*(plugin->ports[ORDER]));

    FilterClass *lpf = plugin->lpf;

    if ( plugin->n_samples != n_samples )
    {
        plugin->Realloc(n_samples);
        return;
    }

    if (VectorNorm1(in, n_samples) == 0)
    {
        fill_n(out, n_samples, 0);
        return;
    }

    copy_n(in, n_samples, lpf->u);
    switch (Order+1)
    {
        case 1:
            lpf->LPF1_Bilinear(f);
            break;
        case 2:
            lpf->LPF2_Bilinear(f);
            break;
        case 3:
            lpf->LPF3_Bilinear(f);
            break;
    }
    copy_n(lpf->y, n_samples, out);
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