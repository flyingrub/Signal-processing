#include "Wave.cpp"
#include <math.h>
#include <iostream>
#include <stdlib.h>

using namespace std;

#define SAMPLE_RATE 22050.0
#define PI 3.141592654

unsigned char doubleToUnsignedChar(double d) {
    return (unsigned char) floor((d+1.0)*127.5);
}

double unsignedCharToDouble(unsigned char c) {
    return (double) (c / 128.0) - 1.0;
}

void tabDoubleToUCNorm(unsigned char* res, double* d, int size) {
    double mini = d[0], maxi = d[0] ;
    for (int i = 0; i < size; i++) {
        maxi = maxi<d[i] ? d[i] : maxi ;
        mini = mini>d[i] ? d[i] : mini ;
    }
    maxi = maxi - mini ;
    maxi = maxi<1e-16 ? 1e-16 : maxi ;
    for (int i = 0; i < size; i++) {
        res[i] = doubleToUnsignedChar((d[i]-mini)/maxi);
    }
}

void tabDoubleToUC(unsigned char* res, double* d, int size) {
    for (int i = 0; i < size; i++) {
        res[i] = doubleToUnsignedChar(d[i]);
    }
}

void tabUnsignedCtoDouble(double* res, unsigned char* c, int size) {
    for (int i = 0; i < size; i++) {
        res[i] = unsignedCharToDouble(c[i]);
    }
}

const int getBufferSize(double duration) {
    return SAMPLE_RATE * duration;
}

int freqGeneratorUC(unsigned char* data, int buffer_size, double frequency, double amplitude) {
    double freq_num = frequency / SAMPLE_RATE ;

    for (int i = 0; i < buffer_size; i++) {
        double d = amplitude * sin(2.0 * PI * freq_num * (double)i);
        data[i] = doubleToUnsignedChar(d);
    }
    return buffer_size;
}

int freqGeneratorDouble(double* data, int buffer_size, double frequency, double amplitude) {
    double freq_num = frequency / SAMPLE_RATE ;

    for (int i = 0; i < buffer_size; i++) {
        data[i] = amplitude * sin(2.0 * PI * freq_num * (double)i);
    }
    return buffer_size;
}

void DFT(double *signal_data, double *partie_reelle, double *partie_imaginaire, int N)  {
    for (int k = 0; k < N; k++) {
        partie_reelle[k] = 0;
        partie_imaginaire[k] = 0;
        for (int n = 0; n < N; n++) {
            double angle = 2 * PI * k * n / N;
            partie_reelle[k] += signal_data[n] * cos(angle);
            partie_imaginaire[k] += -signal_data[n] * sin(angle);
        }
    }
}

void IDFT(double *signal_data, double *partie_reelle, double *partie_imaginaire, int N) {
    for (int k = 0; k < N; k++) {
        signal_data[k] = 0;
        for (int n = 0; n < N; n++) {
            double angle = 2 * PI * k * n / N;
            signal_data[k] += partie_reelle[n] * cos(angle) - partie_imaginaire[n] * sin(angle);
        }
        signal_data[k] /= N;
    }
}

double getAlpha(double frequency, double sample_rate) {
    return PI * frequency / sample_rate;
}

void filterDFTsin(double* partie_reelle, double* partie_imaginaire, int buffer_size, double freqCoupure, bool passe_haut) {
    double freqMin = SAMPLE_RATE / buffer_size;
    double freqMax = SAMPLE_RATE / 2;
    double freqBinStep = (freqMax - freqMin) / buffer_size;
    double freq = freqMin;
    int i=0;
    for(; freq <= freqCoupure; i++) {
        freq = freqMin + freqBinStep * i;
        int repli = SAMPLE_RATE-i;
        partie_reelle[i] = 0;
        partie_reelle[repli] = 0;
        partie_imaginaire[i] = 0;
        partie_imaginaire[repli] = 0;
    }
    cout << "FREQ " << freq  << " | i " << i << endl;
}

void filtreButterworth3(double *input, double *output, int N, double fCoupure, double sample_rate, bool passe_haut) {
    double alpha = getAlpha(fCoupure, sample_rate);
    double a0 = 1 + 2 * alpha + 2 * pow(alpha, 2) + pow(alpha, 3);
    double a1 = -3 - 2 * alpha + 2 * pow(alpha, 2) + 3 * pow(alpha, 3);
    double a2 = 3 - 2 * alpha - 2 * pow(alpha, 2) + 3 * pow(alpha, 3);
    double a3 = -1 + 2 * alpha - 2 * pow(alpha, 2) + pow(alpha, 3);
    double b0 = 1 * pow(alpha, 3);
    double b1 = 3 * b0;
    double b2 = 3 * b0;
    double b3 = b0;
    for (int i = 0; i < N; i++) {
        if (i < 3) {
            output[i] = input[i];
        } else {
            double out = (b0*input[i] + b1 * input[i-1] + b2 * input[i-2] + b3 * input[i-3]
                        - a1 * output[i-1] - a2 * output[i-2] - a3 * output[i-3]);
            output[i] = (out / a0);
        }
    }
    if (passe_haut) {
        for (int i = 0; i < N; i++) {
            output[i] = input[i] - output[i];
        }
    }
}

void test_write() {
    const int buffer_size = getBufferSize(6);
    unsigned char s_dataUC[buffer_size];

    freqGeneratorUC(s_dataUC, buffer_size, 440.0, 0.5);

    Wave* wave = new Wave(s_dataUC, buffer_size, 1, SAMPLE_RATE);
    wave->write("test.wav");
}

void test_dft() {
    double *s_data;
    unsigned char *s_data_UC;
    double *s_data_butter;
    unsigned char* s_data_read;
    int buffer_size;
    double sample_rate;

    Wave* read = new Wave();
    read->read("GammePiano.wav");
    read->getData8(&s_data_read, &buffer_size);
    s_data = new double [buffer_size];
    s_data_UC = new unsigned char [buffer_size];
    s_data_butter = new double [buffer_size];
    sample_rate = (double) read->getSampleRate();
    double partie_imaginaire[buffer_size];
    double partie_reelle[buffer_size];

    tabUnsignedCtoDouble(s_data, s_data_read, buffer_size);
    DFT(s_data, partie_reelle, partie_imaginaire, buffer_size);

    filterDFTsin(partie_reelle, partie_imaginaire, buffer_size, 400.0, false);



    IDFT(s_data, partie_reelle, partie_imaginaire, buffer_size);
    tabDoubleToUC(s_data_UC, s_data, buffer_size);

    Wave* waveDFT = new Wave(s_data_UC, buffer_size, 1, SAMPLE_RATE);
    waveDFT->write("test_DFT_IDFT.wav");
}

void test_bw() {
    double *s_data;
    unsigned char *s_data_UC;
    double *s_data_butter;
    unsigned char* s_data_read;
    int buffer_size;
    double sample_rate;

    Wave* read = new Wave();
    read->read("GammePiano.wav");
    read->getData8(&s_data_read, &buffer_size);
    s_data = new double [buffer_size];
    s_data_UC = new unsigned char [buffer_size];
    s_data_butter = new double [buffer_size];
    sample_rate = (double) read->getSampleRate();

    tabUnsignedCtoDouble(s_data, s_data_read, buffer_size);
    filtreButterworth3(s_data, s_data_butter, buffer_size, 10000, sample_rate, true);
    tabDoubleToUC(s_data_UC, s_data_butter, buffer_size);

    Wave* waveButter = new Wave(s_data_UC, buffer_size, 1, sample_rate);
    waveButter->write("test_butter.wav");
}

int main() {
    test_dft();
    return 0;
}


