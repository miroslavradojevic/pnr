#include "frangi.h"
#include <string>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cfloat>
#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a)>(b)?(a):(b))
#define n 3
//#ifdef PI
//#undef PI
//#endif
#define PI 3.1415926535897932
using namespace std;

template<typename T>
T clamp(T x, T x1, T x2) {
    T xC = (x<x1)?x1:x;
    return (xC>x2)?x2:xC;
}

Frangi::Frangi(float _sig_min, float _sig_stp, float _sig_max, float _zdist, float _alpha, float _beta, float _C, float _beta_one, float _beta_two) {

    alpha   = _alpha;
    beta    = _beta;
    C       = _C;
    BetaOne = _beta_one;
    BetaTwo = _beta_two;
    zdist   = _zdist;

    for (float s = _sig_min; s <= _sig_max+FLT_MIN; s+=_sig_stp) {
        sig.push_back(s);
    }

    blackwhite = false;

}

Frangi::~Frangi(){}

void Frangi::filter3d(unsigned char* I, int w, int h, int l, unsigned char* J8, unsigned char* Sc, float* Vx, float* Vy, float* Vz) {

    cout << "Frangi 3d... I["<<w<<","<<h<<","<<l<<"] sig = "<< flush;

    float* F     = new float[w*h*l];
    float* J     = new float[w*h*l];
    float* Dzz   = new float[w*h*l];
    float* Dyy   = new float[w*h*l];
    float* Dyz   = new float[w*h*l];
    float* Dxx   = new float[w*h*l];
    float* Dxy   = new float[w*h*l];
    float* Dxz   = new float[w*h*l];

    double  Ma[3][3];
    double  Davec[3][3];
    double  Daeig[3];

    double Lambda1, Lambda2, Lambda3;
    double Vecx, Vecy, Vecz, Vecn;
    double LambdaAbs1, LambdaAbs2, LambdaAbs3;
    double Ra, Rb, S, expRa, expRb, expS;
    double Voxel_data;

    for (int si = 0; si < sig.size(); ++si) {

        cout << sig[si] << ((si<sig.size()-1)?", ":"") << flush;

        hessian3d(I, w, h, l, sig[si], zdist, Dzz, Dyy, Dyz, Dxx, Dxy, Dxz, F);

        for (long i = 0; i < (w*h*l); ++i) {

            Ma[0][0]=Dxx[i]; Ma[0][1]=Dxy[i]; Ma[0][2]=Dxz[i];
            Ma[1][0]=Dxy[i]; Ma[1][1]=Dyy[i]; Ma[1][2]=Dyz[i];
            Ma[2][0]=Dxz[i]; Ma[2][1]=Dyz[i]; Ma[2][2]=Dzz[i];

            eigen_decomposition(Ma, Davec, Daeig);

            Lambda1 = Daeig[0];
            Lambda2 = Daeig[1];
            Lambda3 = Daeig[2];
            Vecx = Davec[0][0];
            Vecy = Davec[1][0];
            Vecz = Davec[2][0];

            LambdaAbs1 = abs(Lambda1);
            LambdaAbs2 = abs(Lambda2);
            LambdaAbs3 = abs(Lambda3);

            Ra = LambdaAbs2/LambdaAbs3;
            Rb = LambdaAbs1/sqrt(LambdaAbs2*LambdaAbs3);

            S = sqrt(LambdaAbs1*LambdaAbs1+LambdaAbs2*LambdaAbs2+LambdaAbs3*LambdaAbs3);

            expRa = (1 - exp(-((Ra*Ra)/(2*alpha*alpha))));
            expRb =      exp(-((Rb*Rb)/(2*beta*beta)));
            expS  = (1 - exp(-(S*S)/(2*C*C)));

            Voxel_data = expRa * expRb * expS;

            if (blackwhite) {
                Voxel_data = (Lambda2<0)? 0 : Voxel_data;
                Voxel_data = (Lambda3<0)? 0 : Voxel_data;
            }
            else {
                Voxel_data = (Lambda2>0)? 0 : Voxel_data;
                Voxel_data = (Lambda3>0)? 0 : Voxel_data;
            }

            // remove NaN
            Voxel_data = isnan(Voxel_data)? 0 : Voxel_data ;

            // add result of this scale to output
            if (si==0) {
                J[i] = Voxel_data;
                Sc[i] = (unsigned char) si;
                // normalize before adding
                Vecn = sqrt(Vecx*Vecx+Vecy*Vecy+Vecz*Vecz);
                Vx[i] = Vecx/Vecn;
                Vy[i] = Vecy/Vecn;
                Vz[i] = Vecz/Vecn;
            }
            else {
                if (Voxel_data>J[i]) {
                    J[i] = Voxel_data; // keep maximum filter response
                    Sc[i] = (unsigned char) si;
                    Vecn = sqrt(Vecx*Vecx+Vecy*Vecy+Vecz*Vecz);
                    Vx[i] = Vecx/Vecn;
                    Vy[i] = Vecy/Vecn;
                    Vz[i] = Vecz/Vecn;
                }
            }
        } // pix
    }

    // convert min-max normalized J to J8 (byte8 image)
    float Jmin = J[0];
    float Jmax = J[0];

    for (long i = 1; i < (w*h*l); ++i) {
        if (J[i]<Jmin) Jmin = J[i];
        if (J[i]>Jmax) Jmax = J[i];
    }

    if (abs(Jmax-Jmin)<=FLT_MIN) J8 = 0;

    // export into byte8
    unsigned char J8min = 255;
    unsigned char J8max = 0;
    for (long i = 0; i < (w*h*l); ++i) {
        int val = round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
        val = (val<0)?0:val;
        val = (val>255)?255:val;
        J8[i] = (unsigned char) val;
        if (J8[i]<J8min) J8min = J8[i];
        if (J8[i]>J8max) J8max = J8[i];
    }

    delete [] F;   F = 0;
    delete [] J;   J = 0;
    delete [] Dzz; Dzz = 0;
    delete [] Dyy; Dyy = 0;
    delete [] Dyz; Dyz = 0;
    delete [] Dxx; Dxx = 0;
    delete [] Dxy; Dxy = 0;
    delete [] Dxz; Dxz = 0;

}

void Frangi::hessian3d(unsigned char* I, int w, int h, int l, float sig, float zdist, float* Dzz, float* Dyy, float* Dyz, float* Dxx, float* Dxy, float* Dxz, float* F) {

    imgaussian(I, w, h, l, sig, zdist, F);

    int x,y,z;

    // Dz
    float* Dz = new float[w*h*l];
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (z==0)           Dz[i] =     F[(z+1)*w*h+y*w+x]  - F[i];
        else if (z<(l-1))   Dz[i] = .5*(F[(z+1)*w*h+y*w+x]  - F[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dz[i] =     F[i]                - F[(z-1)*w*h+y*w+x];
    }
    // Dzz
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (z==0)           Dzz[i] =     Dz[(z+1)*w*h+y*w+x]  - Dz[i];
        else if (z<(l-1))   Dzz[i] = .5*(Dz[(z+1)*w*h+y*w+x]  - Dz[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dzz[i] =     Dz[i]                - Dz[(z-1)*w*h+y*w+x];

        Dzz[i] *= (sig*sig); // correct for scaling
    }

    delete [] Dz; Dz = 0;

    // Dy
    float* Dy = new float[w*h*l];
    for (long i = 0; i < (w*h*l); ++i) {
       x = i%w; z = i/(w*h); y = i/w-z*h;
       if (y==0)            Dy[i] =     F[z*w*h+(y+1)*w+x] - F[i];
       else if (y<(h-1))    Dy[i] = .5*(F[z*w*h+(y+1)*w+x]  - F[z*w*h+(y-1)*w+x]);
       else if (y==(h-1))   Dy[i] =     F[i]                - F[z*w*h+(y-1)*w+x];
    }

    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        // Dyy
        if (y==0)           Dyy[i] =     Dy[z*w*h+(y+1)*w+x] - Dy[i];
        else if (y<(h-1))   Dyy[i] = .5*(Dy[z*w*h+(y+1)*w+x] - Dy[z*w*h+(y-1)*w+x]);
        else if (y==(h-1))  Dyy[i] =     Dy[i]               - Dy[z*w*h+(y-1)*w+x];

        Dyy[i] *= (sig*sig); // correct for scaling
        // Dyz
        if (z==0)           Dyz[i] =     Dy[(z+1)*w*h+y*w+x] - Dy[i];
        else if (z<(l-1))   Dyz[i] = .5*(Dy[(z+1)*w*h+y*w+x] - Dy[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dyz[i] =     Dy[i]               - Dy[(z-1)*w*h+y*w+x];

        Dyz[i] *= (sig*sig); // correct for scaling
    }

    delete [] Dy; Dy = 0;

    // Dx
    float* Dx = new float[w*h*l];
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (x==0)           Dx[i] =      F[z*w*h+y*w+(x+1)]  - F[i];
        else if (x<(w-1))   Dx[i] =  .5*(F[z*w*h+y*w+(x+1)]  - F[z*w*h+y*w+(x-1)]);
        else if (x==(w-1))  Dx[i] =      F[i]                - F[z*w*h+y*w+(x-1)];
    }

    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        // Dxx
        if (x==0)           Dxx[i] =      Dx[z*w*h+y*w+(x+1)] - Dx[i];
        else if (x<(w-1))   Dxx[i] =  .5*(Dx[z*w*h+y*w+(x+1)] - Dx[z*w*h+y*w+(x-1)]);
        else if (x==(w-1))  Dxx[i] =      Dx[i]               - Dx[z*w*h+y*w+(x-1)];

        Dxx[i] *= (sig*sig); // correct for scaling
        // Dxy
        if (y==0)           Dxy[i] =     Dx[z*w*h+(y+1)*w+x] - Dx[i];
        else if (y<(h-1))   Dxy[i] = .5*(Dx[z*w*h+(y+1)*w+x] - Dx[z*w*h+(y-1)*w+x]);
        else if (y==(h-1))  Dxy[i] =     Dx[i]               - Dx[z*w*h+(y-1)*w+x];

        Dxy[i] *= (sig*sig); // correct for scaling
        // Dxz
        if (z==0)           Dxz[i] =     Dx[(z+1)*w*h+y*w+x] - Dx[i];
        else if (z<(l-1))   Dxz[i] = .5*(Dx[(z+1)*w*h+y*w+x] - Dx[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dxz[i] =     Dx[i]               - Dx[(z-1)*w*h+y*w+x];

        Dxz[i] *= (sig*sig); // correct for scaling
    }

    delete [] Dx;   Dx = 0;

}

void Frangi::filter2d(unsigned char* I, int w, int h, int l, unsigned char* J8, unsigned char* Sc, float* Vx, float* Vy, float* Vz) {
    // unsigned char* Dxx8, unsigned char* Dxy8, unsigned char* Dyy8, unsigned char* L18, unsigned char* L28
    cout << "Frangi 2d... I["<<w<<","<<h<<","<<l<<"] sig = "<< flush;

    float* F     = new float[w*h*l];
    float* J     = new float[w*h*l];
    float* Dxx   = new float[w*h*l];
    float* Dxy   = new float[w*h*l];
    float* Dyy   = new float[w*h*l];
//    float* L1    = new float[w*h*l];
//    float* L2    = new float[w*h*l];

    float tmp, mag;
    float v2x, v2y, v1x, v1y;
    float mu1, mu2;
    float Lambda1, Lambda2;
    float Vecx, Vecy, Vecn;
    float Rb, S2, Ifiltered;

    float beta  = 2*pow(BetaOne,2);
    float c     = 2*pow(BetaTwo,2);

    for (int si = 0; si < sig.size(); ++si) {

        cout << sig[si] << ((si<sig.size()-1)?", ":"") << flush;

        hessian2d(I, w, h, sig[si], Dyy, Dxy, Dxx, F);

        for (long i = 0; i < (w*h*l); ++i) {

            // compute the eigenvectors of J, v1 and v2
            tmp = sqrt(pow(Dxx[i]-Dyy[i],2) + 4*pow(Dxy[i],2));
            v2x = 2*Dxy[i];
            v2y = Dyy[i] - Dxx[i] + tmp;

            // normalize
            mag = sqrt(pow(v2x,2) + pow(v2y,2));
            if (mag>0) {
                v2x /= mag;
                v2y /= mag;
            }

            v1x = -v2y; // eigenvectors
            v1y = v2x;

            // eigenvalues
            mu1 = 0.5 * (Dxx[i] + Dyy[i] + tmp);
            mu2 = 0.5 * (Dxx[i] + Dyy[i] - tmp);

            //sort eigenvalues abs(Lambda1)>abs(Lambda2)
            bool check = abs(mu1)<abs(mu2); // abs(mu1)>abs(mu2) was the original, it was swithched as the x-y axes were switched
            Lambda1 = (check)?mu2:mu1;
            Lambda2 = (check)?mu1:mu2;
            Vecx    = (check)?v2x:v1x;
            Vecy    = (check)?v2y:v1y;

//            L1[i] = Lambda1;
//            L2[i] = Lambda2;

            // extract vesselness measure using Lambda1 and Lambda2
            Lambda1 = (Lambda1==0)?FLT_MIN:Lambda1;
            Rb = pow(Lambda2/Lambda1,2);
            S2 = pow(Lambda1,2) + pow(Lambda2,2);
            Ifiltered = exp(-Rb/beta) * (1-exp(-S2/c));

            if (blackwhite)
                Ifiltered = (Lambda1<0)?0:Ifiltered;
            else
                Ifiltered = (Lambda1>0)?0:Ifiltered;

            // add result of this scale to output
            if (si==0) {
                J[i] = Ifiltered;
                Sc[i] = si;
                // normalize before adding
                Vecn = sqrt(Vecx*Vecx+Vecy*Vecy);
                Vx[i] = Vecx/Vecn;
                Vy[i] = Vecy/Vecn;
                Vz[i] = 0;
            }
            else {
                if (Ifiltered>J[i]) {
                    J[i] = Ifiltered; // keep maximum filter response
                    Sc[i] = si;
                    Vecn = sqrt(Vecx*Vecx+Vecy*Vecy);
                    Vx[i] = Vecx/Vecn;
                    Vy[i] = Vecy/Vecn;
                    Vz[i] = 0;
                }
            }
        }
    } // sigma loop

    // convert min-max normalized J to J8 (byte8 image)
    float Jmin = J[0];
    float Jmax = J[0];

//    float Dxxmin = Dxx[0];
//    float Dxxmax = Dxx[0];

//    float Dyymin = Dyy[0];
//    float Dyymax = Dyy[0];

//    float Dxymin = Dxy[0];
//    float Dxymax = Dxy[0];

//    float L1min = L1[0];
//    float L1max = L1[0];

//    float L2min = L2[0];
//    float L2max = L2[0];

    for (long i = 1; i < (w*h*l); ++i) {

        Jmin = (J[i]<Jmin)?J[i]:Jmin;
        Jmax = (J[i]>Jmax)?J[i]:Jmax;

//        Dxxmin = (Dxx[i]<Dxxmin)?Dxx[i]:Dxxmin;
//        Dxxmax = (Dxx[i]>Dxxmax)?Dxx[i]:Dxxmax;

//        Dyymin = (Dyy[i]<Dyymin)?Dyy[i]:Dyymin;
//        Dyymax = (Dyy[i]>Dyymax)?Dyy[i]:Dyymax;

//        Dxymin = (Dxy[i]<Dxymin)?Dxy[i]:Dxymin;
//        Dxymax = (Dxy[i]>Dxymax)?Dxy[i]:Dxymax;

//        L1min = (L1[i]<L1min)?L1[i]:L1min;
//        L1max = (L1[i]>L1max)?L1[i]:L1max;

//        L2min = (L2[i]<L2min)?L2[i]:L2min;
//        L2max = (L2[i]>L2max)?L2[i]:L2max;

    }

    if (abs(Jmax-Jmin)<=FLT_MIN) J8 = NULL;

    // export into byte8
    for (long i = 0; i < (w*h*l); ++i) {
        int val = round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
        val = (val<0)?0: (val>255)?255:val  ;
        J8[i] = (unsigned char) val;

//        val = round(((Dxx[i]-Dxxmin)/(Dxxmax-Dxxmin)) * 255);
//        val = (val<0)?0: (val>255)?255:val  ;
//        Dxx8[i] = (unsigned char) val;

//        val = round(((Dyy[i]-Dyymin)/(Dyymax-Dyymin)) * 255);
//        val = (val<0)?0: (val>255)?255:val  ;
//        Dyy8[i] = (unsigned char) val;

//        val = round(((Dxy[i]-Dxymin)/(Dxymax-Dxymin)) * 255);
//        val = (val<0)?0: (val>255)?255:val  ;
//        Dxy8[i] = (unsigned char) val;

//        val = round(((L1[i]-L1min)/(L1max-L1min)) * 255);
//        val = (val<0)?0: (val>255)?255:val  ;
//        L18[i] = (unsigned char) val;

//        val = round(((L2[i]-L2min)/(L2max-L2min)) * 255);
//        val = (val<0)?0: (val>255)?255:val  ;
//        L28[i] = (unsigned char) val;

    }

    delete [] F;   F = 0;
    delete [] J;   J = 0;
    delete [] Dyy; Dyy = 0;
    delete [] Dxx; Dxx = 0;
    delete [] Dxy; Dxy = 0;
//    delete [] L1;   L1 = 0;
//    delete [] L2;   L2 = 0;
}

void Frangi::hessian2d(unsigned char* I, int w, int h, float sig, float* Dyy, float* Dxy, float* Dxx, float* F) {

    imgaussian(I, w, h, sig, F);

    int x,y;

    // Dy
    float* Dy = new float[w*h];
    for (long i = 0; i < (w*h); ++i) {
       x = i%w; y = i/w;
       if (y==0)            Dy[i] =     F[(y+1)*w+x] - F[i];
       else if (y<(h-1))    Dy[i] = .5*(F[(y+1)*w+x] - F[(y-1)*w+x]);
       else if (y==(h-1))   Dy[i] =     F[i]         - F[(y-1)*w+x];
    }

    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        // Dyy
        if (y==0)           Dyy[i] =     Dy[(y+1)*w+x] - Dy[i];
        else if (y<(h-1))   Dyy[i] = .5*(Dy[(y+1)*w+x] - Dy[(y-1)*w+x]);
        else if (y==(h-1))  Dyy[i] =     Dy[i]         - Dy[(y-1)*w+x];

        Dyy[i] *= (sig*sig); // correct for scaling
    }

    delete [] Dy; Dy = 0;

    // Dx
    float* Dx = new float[w*h];
    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        if (x==0)           Dx[i] =      F[y*w+(x+1)]  - F[i];
        else if (x<(w-1))   Dx[i] =  .5*(F[y*w+(x+1)]  - F[y*w+(x-1)]);
        else if (x==(w-1))  Dx[i] =      F[i]          - F[y*w+(x-1)];
    }

    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        // Dxx
        if (x==0)           Dxx[i] =      Dx[y*w+(x+1)] - Dx[i];
        else if (x<(w-1))   Dxx[i] =  .5*(Dx[y*w+(x+1)] - Dx[y*w+(x-1)]);
        else if (x==(w-1))  Dxx[i] =      Dx[i]         - Dx[y*w+(x-1)];

        Dxx[i] *= (sig*sig); // correct for scaling
        // Dxy
        if (y==0)           Dxy[i] =     Dx[(y+1)*w+x] - Dx[i];
        else if (y<(h-1))   Dxy[i] = .5*(Dx[(y+1)*w+x] - Dx[(y-1)*w+x]);
        else if (y==(h-1))  Dxy[i] =     Dx[i]         - Dx[(y-1)*w+x];

        Dxy[i] *= (sig*sig); // correct for scaling
    }

    delete [] Dx;   Dx = 0;

    // kernel
//    vector<int> X,Y;  // offsets
//    vector<float> dGxx, dGxy, dGyy; // gaussian derivative kernel weights (at flipped x,y offsets for convolution)
//    int L = ceil(3*sig);
//    for (int x = -L; x <= L; ++x) {
//        for (int y = -L; y <= L; ++y) {
//            X.push_back(x);
//            Y.push_back(y);
//            // use -x,-y for convolution, x,y for correlation
//            dGxx.push_back( (1/(2*PI*pow(sig,4))) * ((x*x)/pow(sig,2) - 1) * exp(-((x*x) + (y*y))/(2*pow(sig,2))));
//            dGxy.push_back( (1/(2*PI*pow(sig,6))) * ( x*y )                * exp(-((x*x) + (y*y))/(2*pow(sig,2))));
//            dGyy.push_back( (1/(2*PI*pow(sig,4))) * ((y*y)/pow(sig,2) - 1) * exp(-((x*x) + (y*y))/(2*pow(sig,2))));
//        }
//    }

    // convolve I with dGxx,dGxy,dGyy to obtain Dxx, Dxy, Dyy
//    long i0, i1;
//    for (int x = 0; x < w; ++x) {
//        for (int y = 0; y < h; ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig);
//            Dxy[i0] *= (sig*sig);
//            Dyy[i0] *= (sig*sig);
//        }
//    }

//    long i0, i1;
//    for (int x = 0; x < min(L,w); ++x) {
//        for (int y = 0; y < min(L,h); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//        for (int y = L; y < (h-L); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = (y+Y[k])*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }

//        for (int y = max(h-L,L); y < h; ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//    }
//    //----------------------------
//    for (int x = L; x < (w-L); ++x) {
//        for (int y = 0; y < min(L,h); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+(x+X[k]);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//        for (int y = L; y < (h-L); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = (y+Y[k])*w+(x+X[k]);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }

//        for (int y = max(h-L,L); y < h; ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+(x+X[k]);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//    }
//    //----------------------------
//    for (int x = max(w-L,L); x < w; ++x) {
//        for (int y = 0; y < min(L,h); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//        for (int y = L; y < (h-L); ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = (y+Y[k])*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//        for (int y = max(h-L,L); y < h; ++y) {
//            i0 = y*w+x;
//            Dxx[i0] = Dxy[i0] = Dyy[i0] = 0;
//            for (int k = 0; k < X.size(); ++k) {
//                i1 = clamp(y+Y[k],0,h-1)*w+clamp(x+X[k],0,w-1);
//                Dxx[i0] += I[i1] * dGxx[k];
//                Dxy[i0] += I[i1] * dGxy[k];
//                Dyy[i0] += I[i1] * dGyy[k];
//            }
//            Dxx[i0] *= (sig*sig); Dxy[i0] *= (sig*sig); Dyy[i0] *= (sig*sig);
//        }
//    }

}

void Frangi::imgaussian(unsigned char * I, int w, int h, float sig, float* F) {

    long i0, i1;

    // gaussian filter is separated into 1D Gaussian kernel Gxy[.] that will be used for filtering along x and y
    int Lxy = ceil(3*sig);
    float Gxy[Lxy+1];
    float Gnorm = 0;
    for (int i = -Lxy; i <= Lxy; ++i) {
        Gxy[i+Lxy] = exp(-(i*i)/(2*sig*sig));
        Gnorm += Gxy[i+Lxy];
    }

    for (int i = 0; i < (2*Lxy+1); ++i) {
        Gxy[i] /= Gnorm;
    }

    float * K = new float[w*h];

    // Gx gaussian smoothing of the byte8 image array stored in I[.], result in K[.]
    for (int y = 0; y < h; ++y) {
        // x index clamping
        for (int x = 0; x < min(Lxy,w); ++x) {
            i0 = y*w+x;
            K[i0] = 0;
            for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                i1 = y*w+clamp(x1,0,w-1);
                K[i0] += I[i1] * Gxy[x1-x+Lxy];
            }
        }
        // no clamp
        for (int x = Lxy; x < (w-Lxy); ++x) {
            i0 = y*w+x;
            K[i0] = 0;
            for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                i1 = y*w+x1;
                K[i0] += I[i1] * Gxy[x1-x+Lxy];
            }
        }
        // x index clamping
        for (int x = max(w-Lxy,Lxy); x < w; ++x) {
            i0 = y*w+x;
            K[i0] = 0;
            for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                i1 = y*w+clamp(x1,0,w-1);
                K[i0] += I[i1] * Gxy[x1-x+Lxy];
            }
        }
    }

    // Gy gaussian smoothing of the float image array stored in K[.], result in F[.]
    for (int x = 0; x < w; ++x) {
        // y index clamping
        for (int y = 0; y < min(Lxy,h); ++y) {
            i0 = y*w+x;
            F[i0] = 0;
            for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                i1 = clamp(y1,0,h-1)*w+x;
                F[i0] += K[i1] * Gxy[y1-y+Lxy];
            }
        }
        // no clamp
        for (int y = Lxy; y < (h-Lxy); ++y) {
            i0 = y*w+x;
            F[i0] = 0;
            for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                i1 = y1*w+x;
                F[i0] += K[i1] * Gxy[y1-y+Lxy];
            }
        }
        // y index clamping
        for (int y = max(h-Lxy,Lxy); y < h; ++y) {
            i0 = y*w+x;
            F[i0] = 0;
            for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                i1 = clamp(y1,0,h-1)*w+x;
                F[i0] += K[i1] * Gxy[y1-y+Lxy];
            }
        }
    }
    delete [] K; K = 0;
}

void Frangi::imgaussian(unsigned char * I, int w, int h, int l, float sig, float zdist, float* F) {

    long i0, i1;

    float sigz = sig/zdist;

    // gaussian filter is separated into 1D Gaussian kernels Gxy[.] and Gz[.] if z coordinate is scaled
    int Lxy = ceil(3*sig);

    float Gxy[2*Lxy+1];

    float Gnorm = 0;
    for (int i = -Lxy; i <= Lxy; ++i) {
        Gxy[i+Lxy] = exp(-(i*i)/(2*sig*sig));
        Gnorm += Gxy[i+Lxy];
    }

    for (int i = 0; i < (2*Lxy+1); ++i) {
        Gxy[i] /= Gnorm;
    }

    Gnorm = 0;
    int Lz  = ceil(3*sigz);
    float Gz[2*Lz  +1];
    for (int i = -Lz; i <= Lz; ++i) {
        Gz[i+Lz] = exp(-(i*i)/(2*sigz*sigz));
        Gnorm += Gz[i+Lz];
    }

    for (int i = 0; i < (2*Lz+1); ++i) {
       Gz[i] /= Gnorm;
    }

    // Gx gaussian smoothing of the byte8 image array stored in I[.], result in F[.]
    for (int y = 0; y < h; ++y) {
        for (int z = 0; z < l; ++z) {
            // x index clamping
            for (int x = 0; x < min(Lxy,w); ++x) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+clamp(x1,0,w-1);
                    F[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+x1;
                    F[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }
            // x index clamping
            for (int x = max(w-Lxy,Lxy); x < w; ++x) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+clamp(x1,0,w-1);
                    F[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }

        }
    }

    // Gy gaussian smoothing of the float image array stored in F[.], result in K[.]
    float * K = new float[w*h*l];
    for (int x = 0; x < w; ++x) {
        for (int z = 0; z < l; ++z) {
            // y index clamping
            for (int y = 0; y < min(Lxy,h); ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                    K[i0] += F[i1] * Gxy[y1-y+Lxy];
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    i1 = z*w*h+y1*w+x;
                    K[i0] += F[i1] * Gxy[y1-y+Lxy];
                }
            }
            // y index clamping
            for (int y = max(h-Lxy,Lxy); y < h; ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                    K[i0] += F[i1] * Gxy[y1-y+Lxy];
                }
            }
        }
    }

    // Gz gaussain smoothing with sigz, smoothing of the float image array stored in K[.], result in F[.]
    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            // z index clamping
            for (int z = 0; z < min(Lz,l); ++z) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                    F[i0] += K[i1] * Gz[z1-z+Lz];
                }
            }
            // no clamp
            for (int z = Lz; z < (l-Lz); ++z) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                    F[i0] += K[i1] * Gz[z1-z+Lz];
                }
            }
            // z index clamping
            for (int z = max(l-Lz,Lz); z < l; ++z) {
                i0 = z*w*h+y*w+x;
                F[i0] = 0;
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                    F[i0] += K[i1] * Gz[z1-z+Lz];
                }
            }

        }
    }
    delete [] K; K = 0;
}

float Frangi::interpz(int x, int y, float z, float* img, int w, int h, int l){
    // interpolate along z only

//    int xc = (x<0)?0:x;
//    xc = (xc>=w)?(w-1):xc; // xc [0,w-1]

//    int yc = (y<0)?0:y;
//    yc = (yc>=h)?(h-1):yc; // yc [0,h-1]

    if (l==1) { // 2d
        return img[y*w+x];
    }
    else { // 3d
        int z1 = (int) z;
        z1 = (z1<0)?0:z1;
        z1 = (z1>(l-2))?(l-2):z1;       // z1 [0,l-2]
        int z2 = z1 + 1;                // z2 [1,l-1]
        float zfrac = z-z1;
        zfrac = (zfrac<0.0)?0.0:zfrac;
        zfrac = (zfrac>1.0)?1.0:zfrac;  // zfrac [0.0,1.0]

        float I1 = img[z1*w*h+y*w+x];
        float I2 = img[z2*w*h+y*w+x];

        return  (1-zfrac) * I1 + zfrac * I2;
    }

}

void Frangi::eigen_decomposition_static(double A[n][n], double V[n][n], double d[n]) {

    double e[n];
    double da[3];
    double dt, dat;
    double vet[3];
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
                V[i][j] = A[i][j];
        }
    }

    tred2(V, d, e);
    tql2(V, d, e);

    /* Sort the eigen values and vectors by abs eigen value */
    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);

    if((da[0]>=da[1])&&(da[0]>da[2])) {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
    }
    else if((da[1]>=da[0])&&(da[1]>da[2])) {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
        d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2];
    }

    if(da[0]>da[1]) {
        dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
        d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
    }

}

void Frangi::eigen_decomposition(double A[n][n], double V[n][n], double d[n]) {

    double e[n];
    double da[3];
    double dt, dat;
    double vet[3];
    int i, j;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
                V[i][j] = A[i][j];
        }
    }

    tred2(V, d, e);
    tql2(V, d, e);

    /* Sort the eigen values and vectors by abs eigen value */
    da[0]=absd(d[0]); da[1]=absd(d[1]); da[2]=absd(d[2]);

    if((da[0]>=da[1])&&(da[0]>da[2])) {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[0]; da[2]=da[0];  V[0][2] = V[0][0]; V[1][2] = V[1][0]; V[2][2] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
    }
    else if((da[1]>=da[0])&&(da[1]>da[2])) {
        dt=d[2];   dat=da[2];    vet[0]=V[0][2];    vet[1]=V[1][2];    vet[2]=V[2][2];
        d[2]=d[1]; da[2]=da[1];  V[0][2] = V[0][1]; V[1][2] = V[1][1]; V[2][2] = V[2][1];
        d[1]=dt;   da[1]=dat;    V[0][1] = vet[0];  V[1][1] = vet[1];  V[2][1] = vet[2];
    }

    if(da[0]>da[1]) {
        dt=d[1];   dat=da[1];    vet[0]=V[0][1];    vet[1]=V[1][1];    vet[2]=V[2][1];
        d[1]=d[0]; da[1]=da[0];  V[0][1] = V[0][0]; V[1][1] = V[1][0]; V[2][1] = V[2][0];
        d[0]=dt;   da[0]=dat;    V[0][0] = vet[0];  V[1][0] = vet[1];  V[2][0] = vet[2];
    }

}

/* Symmetric Householder reduction to tridiagonal form. */
void Frangi::tred2(double V[n][n], double d[n], double e[n]) {

    /*  This is derived from the Algol procedures tred2 by */
    /*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
    /*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
    /*  Fortran subroutine in EISPACK. */
        int i, j, k;
        double scale;
        double f, g, h;
        double hh;
        for (j = 0; j < n; j++) {d[j] = V[n-1][j]; }

        /* Householder reduction to tridiagonal form. */

        for (i = n-1; i > 0; i--) {
            /* Scale to avoid under/overflow. */
            scale = 0.0;
            h = 0.0;
            for (k = 0; k < i; k++) { scale = scale + fabs(d[k]); }
            if (scale == 0.0) {
                e[i] = d[i-1];
                for (j = 0; j < i; j++) { d[j] = V[i-1][j]; V[i][j] = 0.0;  V[j][i] = 0.0; }
            } else {

                /* Generate Householder vector. */

                for (k = 0; k < i; k++) { d[k] /= scale; h += d[k] * d[k]; }
                f = d[i-1];
                g = sqrt(h);
                if (f > 0) { g = -g; }
                e[i] = scale * g;
                h = h - f * g;
                d[i-1] = f - g;
                for (j = 0; j < i; j++) { e[j] = 0.0; }

                /* Apply similarity transformation to remaining columns. */

                for (j = 0; j < i; j++) {
                    f = d[j];
                    V[j][i] = f;
                    g = e[j] + V[j][j] * f;
                    for (k = j+1; k <= i-1; k++) { g += V[k][j] * d[k]; e[k] += V[k][j] * f; }
                    e[j] = g;
                }
                f = 0.0;
                for (j = 0; j < i; j++) { e[j] /= h; f += e[j] * d[j]; }
                hh = f / (h + h);
                for (j = 0; j < i; j++) { e[j] -= hh * d[j]; }
                for (j = 0; j < i; j++) {
                    f = d[j]; g = e[j];
                    for (k = j; k <= i-1; k++) { V[k][j] -= (f * e[k] + g * d[k]); }
                    d[j] = V[i-1][j];
                    V[i][j] = 0.0;
                }
            }
            d[i] = h;
        }

        /* Accumulate transformations. */

        for (i = 0; i < n-1; i++) {
            V[n-1][i] = V[i][i];
            V[i][i] = 1.0;
            h = d[i+1];
            if (h != 0.0) {
                for (k = 0; k <= i; k++) { d[k] = V[k][i+1] / h;}
                for (j = 0; j <= i; j++) {
                    g = 0.0;
                    for (k = 0; k <= i; k++) { g += V[k][i+1] * V[k][j]; }
                    for (k = 0; k <= i; k++) { V[k][j] -= g * d[k]; }
                }
            }
            for (k = 0; k <= i; k++) { V[k][i+1] = 0.0;}
        }
        for (j = 0; j < n; j++) { d[j] = V[n-1][j]; V[n-1][j] = 0.0; }
        V[n-1][n-1] = 1.0;
        e[0] = 0.0;

}

/* Symmetric tridiagonal QL algorithm. */
void Frangi::tql2(double V[n][n], double d[n], double e[n]) {

/*  This is derived from the Algol procedures tql2, by */
/*  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for */
/*  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding */
/*  Fortran subroutine in EISPACK. */

    int i, j, k, l, m;
    double f;
    double tst1;
    double eps;
    int iter;
    double g, p, r;
    double dl1, h, c, c2, c3, el1, s, s2;

    for (i = 1; i < n; i++) { e[i-1] = e[i]; }
    e[n-1] = 0.0;

    f = 0.0;
    tst1 = 0.0;
    eps = pow(2.0, -52.0);
    for (l = 0; l < n; l++) {

        /* Find small subdiagonal element */

        tst1 = MAX(tst1, fabs(d[l]) + fabs(e[l]));
        m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) { break; }
            m++;
        }

        /* If m == l, d[l] is an eigenvalue, */
        /* otherwise, iterate. */

        if (m > l) {
            iter = 0;
            do {
                iter = iter + 1;  /* (Could check iteration count here.) */
                /* Compute implicit shift */
                g = d[l];
                p = (d[l+1] - g) / (2.0 * e[l]);
                r = hypot2(p, 1.0);
                if (p < 0) { r = -r; }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                dl1 = d[l+1];
                h = g - d[l];
                for (i = l+2; i < n; i++) { d[i] -= h; }
                f = f + h;
                /* Implicit QL transformation. */
                p = d[m]; c = 1.0; c2 = c; c3 = c;
                el1 = e[l+1]; s = 0.0; s2 = 0.0;
                for (i = m-1; i >= l; i--) {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = hypot2(p, e[i]);
                    e[i+1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i+1] = h + s * (c * g + s * d[i]);
                    /* Accumulate transformation. */
                    for (k = 0; k < n; k++) {
                        h = V[k][i+1];
                        V[k][i+1] = s * V[k][i] + c * h;
                        V[k][i] = c * V[k][i] - s * h;
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                /* Check for convergence. */
            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }

    /* Sort eigenvalues and corresponding vectors. */
    for (i = 0; i < n-1; i++) {
        k = i;
        p = d[i];
        for (j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

double Frangi::hypot2(double x, double y) { return sqrt(x*x+y*y); }
