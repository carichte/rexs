/*
    Filename: rebin.c
    To be used with rebin.py, as an imported library. Use Scons to compile,
    simply type 'scons' in the same directory as this file (see www.scons.org).
*/

#define max(a,b) (((a) > (b)) ? (a) : (b))
#define min(a,b) (((a) < (b)) ? (a) : (b))
#include <stdlib.h>

void PyInit_librebin() {} //Python 3.5
void initlibrebin() {} //Python 2.7


double* get_borders(double *x, int length) {
    int i;
    double *borders;
    double right, left;
    borders = malloc(sizeof(double) * (length+1));
    for (i=0; i < length; i++) {
        if (i==0) {
            right = (x[i+1] + x[i])/2;
            left  = 2*x[i] - right;
            borders[i] = left;
        }
        else if (i==(length-1)) {
            left  = (x[i-1] + x[i])/2;
            right = 2*x[i] - left;
        }
        else {
            /*left  = (x[i-1] + x[i])/2;*/
            right = (x[i+1] + x[i])/2;
        }
        borders[i+1] = right;
    }
    return borders;
}


void rebin(double *borders, double *y, double *newborders, double *newy, int length, int newlength) {
    int i,j;
    float right, left, newright, newleft, width, overlap;
    if ((length < 2) || (newlength < 2)) return;
    for (j=0; j < newlength; j++) {
        newleft  = newborders[j];
        newright = newborders[j+1];
        newy[j] = 0;
        for (i=0; i < length; i++) {
            left  = borders[i];
            right = borders[i+1];
            if ((newright < left) || (newleft > right)) continue;
            width = right - left;
            if (width <= 0) continue;
            overlap = min(right, newright) - max(left, newleft);
            overlap = max(0, overlap);
            
            newy[j] = newy[j] + overlap/width * y[i];
        }
    }
}

void rebin_from_centers(double *x, double *y, double *newx, double *newy, int length, int newlength) {
    double *borders;
    double *newborders;
    borders = get_borders(x, length);
    newborders = get_borders(newx, newlength);
    rebin(borders, y, newborders, newy, length, newlength);
    free(borders);
    free(newborders);
    
}


