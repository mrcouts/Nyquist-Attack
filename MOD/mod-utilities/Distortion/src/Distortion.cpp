#include <stdlib.h>
#include <cmath>
#include <lv2.h>
#include <algorithm>
#include "SomeUtilities.h"
#include "DistortionClass.h"
#include "GainClass.h"
#include "Oversample8xClass.h"
#include "Downsample8xClass.h"
#include "FilterClass.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/Distortion"
#define N_SAMPLES_DEFAULT 64
enum {IN, OUT, PRE, DIST, POST, FC1, FC2, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class Plugin
{
public:
    Plugin(double samplerate, uint32_t n_samples){Construct(samplerate, n_samples);}
    ~Plugin(){Destruct();}
    void Construct(double samplerate, uint32_t n_samples)
    {
    	this->n_samples = n_samples;
        SampleRate = samplerate;
        u.zeros(n_samples);
        Dist = new DistortionClass(n_samples);
        Pre = new GainClass(n_samples);
        Post = new GainClass(n_samples);
        Over = new Oversample8xClass(n_samples);
        Down = new Downsample8xClass(n_samples);
        Hpf = new FilterClass3(samplerate, n_samples);
        Lpf = new FilterClass3(samplerate, n_samples);
    }
    void Destruct()
    {
        u.clear();
        delete Dist;
        delete Pre;
        delete Post;
        delete Over;
        delete Down;
        delete Hpf;
        delete Lpf;
    }
    void Realloc(uint32_t n_samples)
    {
        Destruct();
        Construct(SampleRate, n_samples);
    }
    void SetInput(float *in)
    {
        for(uint32_t i = 0; i < n_samples; i++) u(i) = in[i];
    }

    static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features);
    static void activate(LV2_Handle instance);
    static void deactivate(LV2_Handle instance);
    static void connect_port(LV2_Handle instance, uint32_t port, void *data);
    static void run(LV2_Handle instance, uint32_t n_samples);
    static void cleanup(LV2_Handle instance);
    static const void* extension_data(const char* uri);
    float *ports[PLUGIN_PORT_COUNT];
    vec u;

    uint32_t n_samples;
    double SampleRate;
    DistortionClass *Dist;
    GainClass *Pre;
    GainClass *Post;
    Oversample8xClass *Over;
    Downsample8xClass *Down;
    FilterClass3 *Hpf;
    FilterClass3 *Lpf;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    Plugin::instantiate,
    Plugin::connect_port,
    Plugin::activate,
    Plugin::run,
    Plugin::deactivate,
    Plugin::cleanup,
    Plugin::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle Plugin::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    Plugin *plugin = new Plugin(samplerate, N_SAMPLES_DEFAULT);
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void Plugin::activate(LV2_Handle instance){}

/**********************************************************************************************************************************************************/

void Plugin::deactivate(LV2_Handle instance){}

/**********************************************************************************************************************************************************/

void Plugin::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    Plugin *plugin;
    plugin = (Plugin *) instance;
    plugin->ports[port] = (float*) data;
}

/**********************************************************************************************************************************************************/

void Plugin::run(LV2_Handle instance, uint32_t n_samples)
{
    Plugin *plugin;
    plugin = (Plugin *) instance;

    float *in   = plugin->ports[IN];
    float *out  = plugin->ports[OUT];
    float pre = *(plugin->ports[PRE]);
    int dist_curve = (int)*(plugin->ports[DIST]);
    float post = *(plugin->ports[POST]);
    float fc1 = *(plugin->ports[FC1]);
    float fc2 = *(plugin->ports[FC2]);

    vec u = plugin->u;
    DistortionClass *Dist = plugin->Dist;
    GainClass *Pre = plugin->Pre;
    GainClass *Post = plugin->Post;
    Oversample8xClass *Over = plugin->Over;
    Downsample8xClass *Down = plugin->Down;
    FilterClass3 *Hpf = plugin->Hpf;
    FilterClass3 *Lpf = plugin->Lpf;

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

    //Algoritmo
    plugin->SetInput(in);
    Hpf->HPF3(fc1,&u);
    Pre->Gain(pre, &Hpf->y);
    Over->Oversample8x(&Pre->y);
    switch(dist_curve)
    {
        case 0:
            Dist->SoftClip(&Over->y);
            break;
        case 1:
            Dist->ArcTg(&Over->y);
            break;
        case 2:
            Dist->TgH(&Over->y);
            break;
        case 3:
            Dist->HardClip(&Over->y);
            break;
        default:
            Dist->SoftClip(&Over->y);
            break;
    }
    Down->Downsample8x(&Dist->y);
    Post->Gain(post, &Down->y);
    Lpf->LPF3(fc2, &Post->y);
    for (uint32_t i = 0; i < n_samples; i++) out[i] = Lpf->y(i);
}

/**********************************************************************************************************************************************************/

void Plugin::cleanup(LV2_Handle instance)
{
    delete (Plugin *) instance;
}

/**********************************************************************************************************************************************************/

const void* Plugin::extension_data(const char* uri)
{
    return NULL;
}