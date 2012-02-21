#ifndef _WAVE_
#define _WAVE_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define HEADER_LENGTH 44

typedef struct {
    char riff[4];
    int file_length;
    char wave[4];
    char fmt[4] ;
    int fmt_length;
    short format;
    short channels;
    int sample_rate;
    int bytes_per_second;
    short block_align;
    short bits_sample;
    char data_init[4];
    int data_length;
    short* data;
} wave_t;

wave_t* wave_read(const char* file_name);
void wave_write(const char* file_name, wave_t* wave);
void wave_print(wave_t* wave);
void wave_duplicate(wave_t* wave);
void wave_encrypt(wave_t* wave);
void wave_reduce(wave_t* wave);
void wave_enlarge(wave_t* wave);
void wave_trim(wave_t* wave, int start, int end);
void wave_no_silence(wave_t* wave);
void wave_reverse(wave_t* wave);
void wave_echo(wave_t* wave);
int wave_is_equal(wave_t* wave1, wave_t* wave2);
void wave_free(wave_t *wave);
void wave_reload(wave_t *wave, const char* file_name);

void wave_filter_highpass(wave_t *wave, double dt, double rc);
void wave_filter_mean(wave_t* wave, int range);
void wave_filter_median(wave_t* wave, int range);
void wave_filter_gaussian(wave_t* wave, int range);

void quick_sort(short *array, int begin, int end);
void merge_sort(short *array, int size);

void print_chars(const char* title, const char* chars, int length);
#endif

