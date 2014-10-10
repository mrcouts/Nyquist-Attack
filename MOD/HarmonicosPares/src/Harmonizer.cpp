#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lv2.h>
#include <complex>
#include "PitchDetection.h"
#include "shift.h"
#include "window.h"
#include "angle.h"
#include <fftw3.h>
#include "MatchedFilters.h"

/**********************************************************************************************************************************************************/

#define PLUGIN_URI "http://portalmod.com/plugins/mod-devel/Harmonizer"
#define TAMANHO_DO_BUFFER 128
enum {IN, OUT_1, OUT_2, TONE, SCALE, INTERVAL, MODE, LOWNOTE, GAIN_1, GAIN_2, PLUGIN_PORT_COUNT};

/**********************************************************************************************************************************************************/

class PitchShifter
{
public:
    PitchShifter() {}
    ~PitchShifter() {}
    static LV2_Handle instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features);
    static void activate(LV2_Handle instance);
    static void deactivate(LV2_Handle instance);
    static void connect_port(LV2_Handle instance, uint32_t port, void *data);
    static void run(LV2_Handle instance, uint32_t n_samples);
    static void cleanup(LV2_Handle instance);
    static const void* extension_data(const char* uri);
    //Ports
    float *in;
    float *out_1;
    float *out_2;
    float *tone;
    float *scale;
    float *interval;
    float *mode;
    float *lownote;
    float *gain_1;
    float *gain_2;
    
    int nBuffers;
    int nBuffers2;
    int Qcolumn;
    int hopa;
    int N;
    int N2;
    int cont;
    
    double s;
    float g1;
    float g2;
    
    int *Hops;
    double *frames;
    double *ysaida;
    double *ysaida2;
    double *yshift;
    double **b;
    
	float *frames2;
	float *q;
    fftwf_complex *fXa;
    fftwf_complex *fXs;
    
    cx_vec Xa;
    cx_vec Xs;
    cx_vec XaPrevious;
    vec Xa_arg;
    vec XaPrevious_arg;
    vec Phi;
    vec PhiPrevious;
    vec d_phi;
    vec d_phi_prime;
    vec d_phi_wrapped;
    vec omega_true_sobre_fs;
    vec AUX;
    vec Xa_abs;
    vec I;
    vec w;
    
    float *frames3;
    float *q2;
    fftwf_complex *fXa2;
    fftwf_complex *fXs2;
    cx_vec Xa2;
    cx_vec Xs2;
        
    vec R;
	vec NORM;
	vec F;
	vec AUTO;
	
	int note;
	int oitava;
	double frequency;
	double frequency_1;
	
	double SampleRate;
        
    fftwf_plan p;
	fftwf_plan p2;
	fftwf_plan p3;
	fftwf_plan p4;

	double *y;
    double *u;

    double y_1;
    double u_1;
    double y_2;
    double u_2;

    double y1_1;
    double u1_1;
    double y1_2;
    double u1_2;

    double y2_1;
    double u2_1;
    double y2_2;
    double u2_2;

    double y3_1;
    double u3_1;
    double y3_2;
    double u3_2;

    double y4_1;
    double u4_1;
    double y4_2;
    double u4_2;

    double y5_1;
    double u5_1;
    double y5_2;
    double u5_2;

    double y6_1;
    double u6_1;
    double y6_2;
    double u6_2;

    double y7_1;
    double u7_1;
    double y7_2;
    double u7_2;

    double y8_1;
    double u8_1;
    double y8_2;
    double u8_2;

    double y9_1;
    double u9_1;
    double y9_2;
    double u9_2;
};

/**********************************************************************************************************************************************************/

static const LV2_Descriptor Descriptor = {
    PLUGIN_URI,
    PitchShifter::instantiate,
    PitchShifter::connect_port,
    PitchShifter::activate,
    PitchShifter::run,
    PitchShifter::deactivate,
    PitchShifter::cleanup,
    PitchShifter::extension_data
};

/**********************************************************************************************************************************************************/

LV2_SYMBOL_EXPORT
const LV2_Descriptor* lv2_descriptor(uint32_t index)
{
    if (index == 0) return &Descriptor;
    else return NULL;
}

/**********************************************************************************************************************************************************/

LV2_Handle PitchShifter::instantiate(const LV2_Descriptor* descriptor, double samplerate, const char* bundle_path, const LV2_Feature* const* features)
{
    PitchShifter *plugin = new PitchShifter();
    
    plugin->SampleRate = samplerate;
        
    plugin->nBuffers = 16;
    plugin->nBuffers2 = 16;
    plugin->Qcolumn = 1*plugin->nBuffers;
    plugin->hopa = TAMANHO_DO_BUFFER;
    plugin->N = plugin->nBuffers*plugin->hopa;
    plugin->N2 = plugin->nBuffers2*plugin->hopa;
    plugin->cont = 0;
    
    plugin->s = 0;
    plugin->g1 = 1;
    plugin->g2 = 1;  
    
    plugin->Hops = (int*)calloc(plugin->Qcolumn,sizeof(int));
    for (int i=1;i<=plugin->Qcolumn;i++) plugin->Hops[i-1] = plugin->hopa;
        
    plugin->frames = (double*)calloc(plugin->N,sizeof(double));
    plugin->ysaida = (double*)calloc(2*(plugin->N + 2*(plugin->Qcolumn-1)*plugin->hopa),sizeof(double));
    plugin->yshift = (double*)calloc(plugin->hopa,sizeof(double));
    plugin->b = (double**)calloc(plugin->hopa,sizeof(double*));
    
    plugin->frames2 = fftwf_alloc_real(plugin->N);
    plugin->q = fftwf_alloc_real(plugin->N);
    plugin->fXa = fftwf_alloc_complex(plugin->N/2 + 1);
	plugin->fXs = fftwf_alloc_complex(plugin->N/2 + 1);
	
    
    plugin->Xa.zeros(plugin->N/2 + 1); 
	plugin->Xs.zeros(plugin->N/2 + 1);
	plugin->Xa2.zeros(plugin->N + 1);
	plugin->Xs2.zeros(plugin->N + 1);
	plugin->XaPrevious.zeros(plugin->N/2 + 1);
	plugin->Xa_arg.zeros(plugin->N/2 + 1);
	plugin->XaPrevious_arg.zeros(plugin->N/2 + 1);
	plugin->Phi.zeros(plugin->N/2 + 1);
	plugin->PhiPrevious.zeros(plugin->N/2 + 1);
    plugin->d_phi.zeros(plugin->N/2 + 1);
	plugin->d_phi_prime.zeros(plugin->N/2 + 1);
	plugin->d_phi_wrapped.zeros(plugin->N/2 + 1);
	plugin->omega_true_sobre_fs.zeros(plugin->N/2 + 1);
	plugin->AUX.zeros(plugin->N/2 + 1);
	plugin->Xa_abs.zeros(plugin->N/2 + 1);
	plugin->w.zeros(plugin->N); hann(plugin->N,&plugin->w);
	plugin->I.zeros(plugin->N/2 + 1); plugin->I = linspace(0, plugin->N/2,plugin->N/2 + 1);
	
	plugin->frames3 = fftwf_alloc_real(2*plugin->N2); memset(plugin->frames3, 0, 2*plugin->N2 );
	plugin->q2 = fftwf_alloc_real(2*plugin->N2);	
	plugin->fXa2 = fftwf_alloc_complex(plugin->N2 + 1);
	plugin->fXs2 = fftwf_alloc_complex(plugin->N2 + 1);
	
	plugin->R.zeros(plugin->N2);
	plugin->NORM.zeros(plugin->N2);
	plugin->F.zeros(plugin->N2);
	plugin->AUTO.zeros(plugin->N2);
    
    for (int i=1 ; i<= (plugin->nBuffers); i++)
    {
		plugin->b[i-1] = &plugin->frames[(i-1)*plugin->hopa];
	}
	
	plugin->p = fftwf_plan_dft_r2c_1d(plugin->N, plugin->frames2, plugin->fXa, FFTW_ESTIMATE);
	plugin->p2 = fftwf_plan_dft_c2r_1d(plugin->N, plugin->fXs, plugin->q, FFTW_ESTIMATE);
	plugin->p3 = fftwf_plan_dft_r2c_1d(2*plugin->N2, plugin->frames3, plugin->fXa2, FFTW_ESTIMATE);
	plugin->p4 = fftwf_plan_dft_c2r_1d(2*plugin->N2, plugin->fXs2, plugin->q2, FFTW_ESTIMATE);

	plugin->y = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    plugin->u = (double*)malloc(TAMANHO_DO_BUFFER*sizeof(double));
    
    plugin->y_1 = 0;
    plugin->u_1 = 0;
    plugin->y_2 = 0;
    plugin->u_2 = 0;

    plugin->y1_1 = 0;
    plugin->u1_1 = 0;
    plugin->y1_2 = 0;
    plugin->u1_2 = 0;

    plugin->y2_1 = 0;
    plugin->u2_1 = 0;
    plugin->y2_2 = 0;
    plugin->u2_2 = 0;

    plugin->y3_1 = 0;
    plugin->u3_1 = 0;
    plugin->y3_2 = 0;
    plugin->u3_2 = 0;

    plugin->y4_1 = 0;
    plugin->u4_1 = 0;
    plugin->y4_2 = 0;
    plugin->u4_2 = 0;

    plugin->y5_1 = 0;
    plugin->u5_1 = 0;
    plugin->y5_2 = 0;
    plugin->u5_2 = 0;

    plugin->y6_1 = 0;
    plugin->u6_1 = 0;
    plugin->y6_2 = 0;
    plugin->u6_2 = 0;

    plugin->y7_1 = 0;
    plugin->u7_1 = 0;
    plugin->y7_2 = 0;
    plugin->u7_2 = 0;

    plugin->y8_1 = 0;
    plugin->u8_1 = 0;
    plugin->y8_2 = 0;
    plugin->u8_2 = 0;

    plugin->y9_1 = 0;
    plugin->u9_1 = 0;
    plugin->y9_2 = 0;
    plugin->u9_2 = 0;

    plugin->frequency = 80;
	plugin->frequency_1 = 80;
	
    return (LV2_Handle)plugin;
}

/**********************************************************************************************************************************************************/

void PitchShifter::activate(LV2_Handle instance)
{
    // TODO: include the activate function code here
}

/**********************************************************************************************************************************************************/

void PitchShifter::deactivate(LV2_Handle instance)
{
    // TODO: include the deactivate function code here
}

/**********************************************************************************************************************************************************/

void PitchShifter::connect_port(LV2_Handle instance, uint32_t port, void *data)
{
    PitchShifter *plugin;
    plugin = (PitchShifter *) instance;

    switch (port)
    {
        case IN:
            plugin->in = (float*) data;
            break;
        case OUT_1:
            plugin->out_1 = (float*) data;
            break;
        case OUT_2:
            plugin->out_2 = (float*) data;
            break;
        case TONE:
            plugin->tone = (float*) data;
            break;
        case SCALE:
            plugin->scale = (float*) data;
            break;
        case INTERVAL:
            plugin->interval = (float*) data;
            break;
        case MODE:
            plugin->mode = (float*) data;
            break;
        case LOWNOTE:
            plugin->lownote = (float*) data;
            break;
        case GAIN_1:
            plugin->gain_1 = (float*) data;
            break;
        case GAIN_2:
            plugin->gain_2 = (float*) data;
            break;
    }
}

/**********************************************************************************************************************************************************/

void PitchShifter::run(LV2_Handle instance, uint32_t n_samples)
{
    PitchShifter *plugin;
    plugin = (PitchShifter *) instance;
    
    if ( ((plugin->hopa) != (int)n_samples) )
    {
		
		switch ((int)n_samples)
		{
			case 64:
				plugin->nBuffers = 16;
				plugin->nBuffers2 = 16;
				break;
			case 128:
				plugin->nBuffers = 8;
				plugin->nBuffers2 = 8;
				break;
			case 256:
				plugin->nBuffers = 4;
				plugin->nBuffers2 = 4;
				break;
			case 512:
				plugin->nBuffers = 2;
				plugin->nBuffers2 = 2;
				break;
		}
		
		plugin->Qcolumn = plugin->nBuffers;
		
		plugin->hopa = n_samples;
		plugin->N = plugin->nBuffers*plugin->hopa;
		plugin->N2 = plugin->nBuffers2*plugin->hopa;
		
		plugin->Hops = (int*)realloc(plugin->Hops,plugin->Qcolumn*sizeof(int));
		for (int i=1;i<=plugin->Qcolumn;i++) plugin->Hops[i-1] = plugin->hopa;
		
		plugin->frames = (double*)realloc(plugin->frames,plugin->N*sizeof(double));                                          memset(plugin->frames, 0, plugin->N*sizeof(double) );
		plugin->ysaida = (double*)realloc(plugin->ysaida,2*(plugin->N + 2*(plugin->Qcolumn-1)*plugin->hopa)*sizeof(double)); memset(plugin->ysaida, 0, 2*(plugin->N + 2*(plugin->Qcolumn-1)*plugin->hopa)*sizeof(double) );
		plugin->yshift = (double*)realloc(plugin->yshift,plugin->hopa*sizeof(double));                                       memset(plugin->yshift, 0, plugin->hopa*sizeof(double) );
		plugin->b = (double**)realloc(plugin->b,plugin->hopa*sizeof(double*));
		
		fftwf_free(plugin->frames2); plugin->frames2 = fftwf_alloc_real(plugin->N);
		fftwf_free(plugin->q);       plugin->q = fftwf_alloc_real(plugin->N);
		fftwf_free(plugin->fXa);     plugin->fXa = fftwf_alloc_complex(plugin->N/2 + 1);
		fftwf_free(plugin->fXs);     plugin->fXs = fftwf_alloc_complex(plugin->N/2 + 1);
		
		plugin->Xa.zeros(plugin->N/2 + 1); 
		plugin->Xs.zeros(plugin->N/2 + 1); 
		plugin->XaPrevious.zeros(plugin->N/2 + 1);
		plugin->Xa_arg.zeros(plugin->N/2 + 1);
		plugin->XaPrevious_arg.zeros(plugin->N/2 + 1);
		plugin->Phi.zeros(plugin->N/2 + 1);
		plugin->PhiPrevious.zeros(plugin->N/2 + 1);
		plugin->d_phi.zeros(plugin->N/2 + 1);
		plugin->d_phi_prime.zeros(plugin->N/2 + 1);
		plugin->d_phi_wrapped.zeros(plugin->N/2 + 1);
		plugin->omega_true_sobre_fs.zeros(plugin->N/2 + 1);
		plugin->AUX.zeros(plugin->N/2 + 1);
		plugin->Xa_abs.zeros(plugin->N/2 + 1);
		plugin->w.zeros(plugin->N); hann(plugin->N,&plugin->w);
		plugin->I.zeros(plugin->N/2 + 1); plugin->I = linspace(0, plugin->N/2,plugin->N/2 + 1);
		
		fftwf_free(plugin->frames3); plugin->frames3 = fftwf_alloc_real(2*plugin->N2); memset(plugin->frames3, 0, 2*plugin->N2 );
		fftwf_free(plugin->q2);      plugin->q2 = fftwf_alloc_real(2*plugin->N2);
		fftwf_free(plugin->fXa2);    plugin->fXa2 = fftwf_alloc_complex(plugin->N2 + 1);
		fftwf_free(plugin->fXs2);    plugin->fXs2 = fftwf_alloc_complex(plugin->N2 + 1);		
		
		plugin->R.zeros(plugin->N2);
		plugin->NORM.zeros(plugin->N2);
		plugin->F.zeros(plugin->N2);
		plugin->AUTO.zeros(plugin->N2);	
		
		for (int i=1 ; i<= plugin->nBuffers; i++)
		{
			plugin->b[i-1] = &plugin->frames[(i-1)*plugin->hopa];
		}
		
		fftwf_destroy_plan(plugin->p);  plugin->p = fftwf_plan_dft_r2c_1d(plugin->N, plugin->frames2, plugin->fXa, FFTW_ESTIMATE);
		fftwf_destroy_plan(plugin->p2); plugin->p2 = fftwf_plan_dft_c2r_1d(plugin->N, plugin->fXs, plugin->q, FFTW_ESTIMATE);
		fftwf_destroy_plan(plugin->p3); plugin->p3 = fftwf_plan_dft_r2c_1d(2*plugin->N2, plugin->frames3, plugin->fXa2, FFTW_ESTIMATE);
		fftwf_destroy_plan(plugin->p4); plugin->p4 = fftwf_plan_dft_c2r_1d(2*plugin->N2, plugin->fXs2, plugin->q2, FFTW_ESTIMATE);
		
		return;
	}

    float media = 0;
    
    for (uint32_t i=1; i<=n_samples; i++)
    {
		media = media + abs(plugin->in[i-1]);
	}
	
	if (media == 0)
	{
		for (uint32_t i=1; i<=n_samples; i++)
		{
			plugin->out_1[i-1] = 0;
			plugin->out_2[i-1] = 0;
		}
	}
	else
	{
		
	int Tone = (int)(*(plugin->tone));
	int Scale = (int)(*(plugin->scale));
	int Interval = (int)(*(plugin->interval));
	int Mode = (int)(*(plugin->mode));
	int LowNote = (int)(*(plugin->lownote));
		    
    double g1_before = plugin->g1;
    double g2_before = plugin->g2;
    plugin->g1 = pow(10, (float)(*(plugin->gain_1))/20.0);
    plugin->g2 = pow(10, (float)(*(plugin->gain_2))/20.0);
    
	for (int k=1; k<= plugin->Qcolumn-1; k++)
    {
		plugin->Hops[k-1] = plugin->Hops[k];
	}
    
    
		for (int i=1; i<=plugin->hopa; i++)
		{
			for (int j=1; j<=(plugin->nBuffers-1); j++)
			{
				plugin->b[j-1][i-1] = plugin->b[j][i-1];
			}
			plugin->b[plugin->nBuffers-1][i-1] = plugin->in[i-1];
		}
		
		if ( plugin->cont < plugin->nBuffers-1)
		{
			plugin->cont = plugin->cont + 1;
		}
		else
		{
			double freq = FindNote(plugin->N2, plugin->frames, plugin->frames3, &plugin->Xa2, &plugin->Xs2, plugin->q2, plugin->Qcolumn, plugin->p3, plugin->p4, plugin->fXa2, plugin->fXs2, &plugin->R, &plugin->NORM, &plugin->F, &plugin->AUTO, plugin->SampleRate, &plugin->note, &plugin->oitava );
			//cout << "freq = " << freq << "\n";
			if (freq != -1) plugin->frequency = freq;
			cout << "freq = " << plugin->frequency << "\n";

			for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->frames[i];
			}
		    
		    SBF1(plugin->u, plugin->y, n_samples, 2*plugin->frequency_1, 2*plugin->frequency, 20, 20, 1/plugin->SampleRate, &plugin->u_1, &plugin->u_2, &plugin->y_1, &plugin->y_2);
		    plugin->frequency_1 = plugin->frequency;

		    

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 4*plugin->frequency_1, 4*plugin->frequency, 40, 40, 1/plugin->SampleRate, &plugin->u1_1, &plugin->u1_2, &plugin->y1_1, &plugin->y1_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 6*plugin->frequency_1, 6*plugin->frequency, 60, 60, 1/plugin->SampleRate, &plugin->u2_1, &plugin->u2_2, &plugin->y2_1, &plugin->y2_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 8*plugin->frequency_1, 8*plugin->frequency, 80, 80, 1/plugin->SampleRate, &plugin->u3_1, &plugin->u3_2, &plugin->y3_1, &plugin->y3_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 10*plugin->frequency_1, 10*plugin->frequency, 100, 100, 1/plugin->SampleRate, &plugin->u4_1, &plugin->u4_2, &plugin->y4_1, &plugin->y4_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 12*plugin->frequency_1, 12*plugin->frequency, 120, 120, 1/plugin->SampleRate, &plugin->u5_1, &plugin->u5_2, &plugin->y5_1, &plugin->y5_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 14*plugin->frequency_1, 14*plugin->frequency, 140, 140, 1/plugin->SampleRate, &plugin->u6_1, &plugin->u6_2, &plugin->y6_1, &plugin->y6_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 16*plugin->frequency_1, 16*plugin->frequency, 160, 160, 1/plugin->SampleRate, &plugin->u7_1, &plugin->u7_2, &plugin->y7_1, &plugin->y7_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 18*plugin->frequency_1, 18*plugin->frequency, 180, 180, 1/plugin->SampleRate, &plugin->u8_1, &plugin->u8_2, &plugin->y8_1, &plugin->y8_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->u[i] = plugin->y[i];
			}

			SBF1(plugin->u, plugin->y, n_samples, 20*plugin->frequency_1, 20*plugin->frequency, 200, 200, 1/plugin->SampleRate, &plugin->u9_1, &plugin->u9_2, &plugin->y9_1, &plugin->y9_2);
		    plugin->frequency_1 = plugin->frequency;

		    for (uint32_t i=0; i < n_samples;i++)
		    {
				plugin->out_1[i] = (float)plugin->y[i];
				plugin->out_2[i] = plugin->frames[i];
			}
			
		}
		
	}

}

/**********************************************************************************************************************************************************/

void PitchShifter::cleanup(LV2_Handle instance)
{
	PitchShifter *plugin;
	plugin = (PitchShifter *) instance;
	
	free(plugin->Hops);
	free(plugin->frames);
	free(plugin->ysaida);
	free(plugin->yshift);
	free(plugin->b);
	
	fftwf_free(plugin->frames2);
	fftwf_free(plugin->q);
	fftwf_free(plugin->fXa);
	fftwf_free(plugin->fXs);
	
	fftwf_free(plugin->frames3);
	fftwf_free(plugin->q2);
	fftwf_free(plugin->fXa2);
	fftwf_free(plugin->fXs2);
	
	plugin->Xa.clear();
	plugin->Xs.clear();
	plugin->Xa2.clear();
	plugin->Xs2.clear();
	plugin->XaPrevious.clear();
	plugin->Xa_arg.clear();
	plugin->XaPrevious_arg.clear();
	plugin->Phi.clear();
	plugin->PhiPrevious.clear();
    plugin->d_phi.clear();
	plugin->d_phi_prime.clear();
	plugin->d_phi_wrapped.clear();
	plugin->omega_true_sobre_fs.clear();
	plugin->AUX.clear();
	plugin->Xa_abs.clear();
	plugin->w.clear();
	plugin->I.clear();
	
	plugin->R.clear();
	plugin->NORM.clear();
	plugin->F.clear();
	plugin->AUTO.clear();
	
	fftwf_destroy_plan(plugin->p);
	fftwf_destroy_plan(plugin->p2);
	fftwf_destroy_plan(plugin->p3);
	fftwf_destroy_plan(plugin->p4);
	
    delete ((PitchShifter *) instance);
}

/**********************************************************************************************************************************************************/

const void* PitchShifter::extension_data(const char* uri)
{
    return NULL;
}
