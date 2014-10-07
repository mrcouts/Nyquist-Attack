#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include "MatchedFilters.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/StopBandFilter"
#define TAMANHO_DO_BUFFER 128
enum {IN, OUT_1, FREQ, Q, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class StopBandFilter
{
public:
    StopBandFilter() {}
    ~StopBandFilter() {}
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
    float *q_factor;
    
    double *y;
    double *u;
    
    double T;
    double f;
    double q;
    double bw;
    double y_1;
    double u_1;
    double y_2;
    double u_2;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    StopBandFilter::instantiate,
    StopBandFilter::connect_port,
    StopBandFilter::activate,
    StopBandFilter::run,
    StopBandFilter::deactivate,
    StopBandFilter::cleanup,
    StopBandFilter::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle StopBandFilter::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    StopBandFilter *plugin = new StopBandFilter();
    
    plugin->y = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->u = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    
    plugin->T = 1/samplerate;
    plugin->f = 1;
    plugin->q = 0.5;
    plugin->bw = 2;
    plugin->y_1 = 0;
    plugin->u_1 = 0;
    plugin->y_2 = 0;
    plugin->u_2 = 0;
    	
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void StopBandFilter::activate(LV2_Handle instance)
{
    // TODO: include the activate function code here
}

/**********************************************************************************************************************************************************/

void StopBandFilter::deactivate(LV2_Handle instance)
{
    // TODO: include the deactivate function code here
}

/**********************************************************************************************************************************************************/

void StopBandFilter::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    StopBandFilter *plugin;
    plugin = (StopBandFilter *) instance;

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
        case Q:
            plugin->q_factor = (float*) data;
            break;
    }
}

/**********************************************************************************************************************************************************/

void StopBandFilter::run(LV2_Handle instance, uint32_t n_samples)
{
    StopBandFilter *plugin;
    plugin = (StopBandFilter *) instance;    
    double f_before = plugin->f;
    double bw_before = plugin->bw;
    
    plugin->f = (float)(*(plugin->freq));
    plugin->q = (float)(*(plugin->q_factor));
    plugin->bw = plugin->f/plugin->q;
    
    for (uint32_t i=0; i<= n_samples-1;i++)
    {
		plugin->u[i] = plugin->in[i];
	}
    
    SBF1(plugin->u, plugin->y, n_samples, f_before, plugin->f, bw_before, plugin->bw, plugin->T, &plugin->u_1, &plugin->u_2, &plugin->y_1, &plugin->y_2);
    
    for (uint32_t i=0; i<= n_samples-1;i++)
    {
		plugin->out_1[i] = plugin->y[i];
	}
}

/**********************************************************************************************************************************************************/

void StopBandFilter::cleanup(LV2_Handle instance)
{
	StopBandFilter *plugin;
	plugin = (StopBandFilter *) instance;
	
	free(plugin->u);
	free(plugin->y);
	
    delete ((StopBandFilter *) instance);
}

/**********************************************************************************************************************************************************/

const void* StopBandFilter::extension_data(const char* uri)
{
    return NULL;
}
