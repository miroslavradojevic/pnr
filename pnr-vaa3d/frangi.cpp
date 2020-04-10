/*
 * c++ implementation of Frangi vesselness with local vessleness direction for given unsigned char 2d/3d image. 
 * Provides with the vesselness image and the eigenvectors. 
 * [paper] http://link.springer.com/chapter/10.1007/BFb0056195
 * [inspiration] https://nl.mathworks.com/matlabcentral/fileexchange/24409-hessian-based-frangi-vesselness-filter
 */
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

unsigned char Frangi::ndirs2d = 30;
unsigned char Frangi::ndirs3d = 90;

template<typename T>
T clamp(T x, T x1, T x2) {
    T xC = (x<x1)?x1:x;
    return (xC>x2)?x2:xC;
}
#ifdef  _WIN32
double round(double r);
#endif

Frangi::Frangi(vector<float> _sigs, float _zdist, float _alpha, float _beta, float _C, float _beta_one, float _beta_two) {

    cout << "Frangi()" << endl;

    alpha   = _alpha;
    beta    = _beta;
    C       = _C;
    BetaOne = _beta_one;
    BetaTwo = _beta_two;
    zdist   = _zdist;

    sig.clear(); // initialize it
    sig = _sigs;

    // initialize directions
    cout << "Vxyz.size()=" << Vxyz.size() << endl;
    // points on sphere stored in


    blackwhite = false;

}

Frangi::~Frangi(){}

void Frangi::generate_3d_unit_directions(unsigned char Ndir, vector< vector<float> >& Vxyz) {

    cout<<"generate_3d_unit_directions()"<<endl;
    // use sphere uniform sampling
    Vxyz.clear();

    double h_k, theta_k, phi_k, phi_k_1 = 0;

    for (int k = 0; k < Ndir; k++) {

        h_k = -1 + 2 * (double)k/(Ndir-1); // range -1:1 for full sphere
        theta_k = acos(h_k);

        if(k==0 || k==(Ndir-1)){
            phi_k   = 0;
            phi_k_1 = 0;
        }
        else{
            phi_k = phi_k_1 + 3.6 / ( sqrt((double)Ndir) * sqrt(1-h_k*h_k));
            phi_k_1 = phi_k;
        }

        vector<float> p(3); // cartesian coordinates
        p[0] = sin(theta_k) * cos(phi_k);
        p[1] = sin(theta_k) * sin(phi_k);
        p[2] = cos(theta_k);

        Vxyz.push_back(p);
    }

}

void Frangi::generate_2d_unit_directions(unsigned char Ndir, vector< vector<float> >& Vxyz) {
    // circle uniform sampling
    Vxyz.clear();

    float ang;

    for (int k = 0; k < Ndir; ++k) {

        ang = k * ((2.0*3.14)/Ndir);
        vector<float> p(3);

        p[0] =  cos(ang);
        p[1] =  sin(ang);
        p[2] =  0.0;

        Vxyz.push_back(p);

    }

}

unsigned char Frangi::get_direction_idx(float vx, float vy, float vz, vector< vector<float> > Vxyz){

    unsigned char idx = 0;
    float d2_min = (vx-Vxyz[0][0])*(vx-Vxyz[0][0]) + (vy-Vxyz[0][1])*(vy-Vxyz[0][1]) + (vz-Vxyz[0][2])*(vz-Vxyz[0][2]);
    float d2; // optimiziation per arc == optimisation per chord (euclidean dist of the points on the sphere)
    for (int k = 1; k < Vxyz.size(); ++k) {
        d2 = (vx-Vxyz[k][0])*(vx-Vxyz[k][0]);
        if (d2<d2_min) {
            d2 += (vy-Vxyz[k][1])*(vy-Vxyz[k][1]);
            if (d2<d2_min) {
                d2 += (vz-Vxyz[k][2])*(vz-Vxyz[k][2]);
                if (d2<d2_min) {
                    idx = k;
                }
            }
        }
    }

    return idx;
}

unsigned char Frangi::get_direction_idx(float vx, float vy, vector< vector<float> > Vxyz){

    unsigned char idx = 0;
    float d2_min = (vx-Vxyz[0][0])*(vx-Vxyz[0][0]) + (vy-Vxyz[0][1])*(vy-Vxyz[0][1]);
    float d2; // optimiziation per arc == optimisation per chord (euclidean dist of the points on the sphere)
    for (int k = 1; k < Vxyz.size(); ++k) {
        d2 = (vx-Vxyz[k][0])*(vx-Vxyz[k][0]);
        if (d2<d2_min) {
            d2 += (vy-Vxyz[k][1])*(vy-Vxyz[k][1]);
            if (d2<d2_min) {
                idx = k;
            }
        }
    }

    return idx;
}

void Frangi::frangi3d(unsigned char* I, int w, int h, int l, float* J, float& Jmin, float& Jmax, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz) {

//    clock_t t1_frangi, t2_frangi;
//    t1_frangi = clock();
    cout << "Frangi 3d... I["<<w<<","<<h<<","<<l<<"]..."<< endl;

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
    double LambdaAbs1, LambdaAbs2, LambdaAbs3;
    double Ra, Rb, S, expRa, expRb, expS;
    double Voxel_data;
    int val;

    clock_t t1, t2;
    Jmin = FLT_MAX;
    Jmax = -FLT_MAX;
    for (int si = 0; si < sig.size(); ++si) {

//        cout << "sig=" << sig[si] << ((si<sig.size()-1)?", ":" ") << endl;
        cout << "Hessian3D[" << sig[si] << "]:" << endl;
//        t1 = clock();
        hessian3d(I, w, h, l, sig[si], zdist, Dzz, Dyy, Dyz, Dxx, Dxy, Dxz);
//        t2 = clock();
//        cout << "                                 " <<((t2-t1)/(double)CLOCKS_PER_SEC) << " sec.\n" << flush;

        t1 = clock();
        cout << "Eigen3D[" << sig[si] << "]... " << flush;
        int PLOT_INTERVAL = (w*h*l)/10;
        for (long i = 0; i < (w*h*l); ++i) {

            if (i%PLOT_INTERVAL==0) cout << ((i/PLOT_INTERVAL)*10) << "%\t" << flush;

            Ma[0][0]=Dxx[i]; Ma[0][1]=Dxy[i]; Ma[0][2]=Dxz[i];
            Ma[1][0]=Dxy[i]; Ma[1][1]=Dyy[i]; Ma[1][2]=Dyz[i];
            Ma[2][0]=Dxz[i]; Ma[2][1]=Dyz[i]; Ma[2][2]=Dzz[i];

            eigen_decomposition(Ma, Davec, Daeig);

            Lambda1 = Daeig[0];
            Lambda2 = Daeig[1];
            Lambda3 = Daeig[2];
//            Vecx = Davec[0][0]; // vx
//            Vecy = Davec[1][0]; // vy
//            Vecz = Davec[2][0]; // vz
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
            #ifdef  _WIN32
            Voxel_data = _isnan(Voxel_data)? 0 : Voxel_data ;
            #elif
            Voxel_data = isnan(Voxel_data)? 0 : Voxel_data ;
            #endif

            // add result of this scale to output
            if (si==0) {
                J[i] = Voxel_data;

                if (J[i]<Jmin) Jmin = J[i];
                if (J[i]>Jmax) Jmax = J[i];

                val = round(  ((Davec[0][0]+1)/2) * 255 ); // vx
                val = (val<0)?0:(val>255)?255:val;
                Vx[i] = (unsigned char) val;

                val = round(  ((Davec[1][0]+1)/2) * 255 ); // vy
                val = (val<0)?0:(val>255)?255:val;
                Vy[i] = (unsigned char) val;

                val = round(  ((Davec[2][0]+1)/2) * 255 ); // vz
                val = (val<0)?0:(val>255)?255:val;
                Vz[i] = (unsigned char) val;

            }
            else {
                if (Voxel_data>J[i]) {
                    J[i] = Voxel_data; // keep maximum filter response

                    if (J[i]<Jmin) Jmin = J[i];
                    if (J[i]>Jmax) Jmax = J[i];

                    val = round(  ((Davec[0][0]+1)/2) * 255 ); // vx
                    val = (val<0)?0:(val>255)?255:val;
                    Vx[i] = (unsigned char) val;

                    val = round(  ((Davec[1][0]+1)/2) * 255 ); // vy
                    val = (val<0)?0:(val>255)?255:val;
                    Vy[i] = (unsigned char) val;

                    val = round(  ((Davec[2][0]+1)/2) * 255 ); // vz
                    val = (val<0)?0:(val>255)?255:val;
                    Vz[i] = (unsigned char) val;
                }
            }
        } // pix
        t2 = clock();
        cout << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

    } // sig

    delete [] Dzz; Dzz = 0;
    delete [] Dyy; Dyy = 0;
    delete [] Dyz; Dyz = 0;
    delete [] Dxx; Dxx = 0;
    delete [] Dxy; Dxy = 0;
    delete [] Dxz; Dxz = 0;

//    t2_frangi = clock();
//    cout << ((t2_frangi-t1_frangi)/(double)CLOCKS_PER_SEC) << " sec." << endl;

}

void Frangi::hessian3d(unsigned char* I, int w, int h, int l, float sig, float zdist, float* Dzz, float* Dyy, float* Dyz, float* Dxx, float* Dxy, float* Dxz){

    clock_t t1, t2;
    t1 = clock();
    cout << "G["<<sig<<"]... " << flush;
    float* F     = new float[w*h*l];
    imgaussian(I, w, h, l, sig, zdist, F);
    t2 = clock();
    cout << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec. " << endl;

    int x,y,z;

    t1 = clock();
    // Dz
    float* DD = new float[w*h*l];
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (z==0)           DD[i] =     F[(z+1)*w*h+y*w+x]  - F[i];
        else if (z<(l-1))   DD[i] = .5*(F[(z+1)*w*h+y*w+x]  - F[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  DD[i] =     F[i]                - F[(z-1)*w*h+y*w+x];
    }
    // Dzz
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (z==0)           Dzz[i] =     DD[(z+1)*w*h+y*w+x]  - DD[i];
        else if (z<(l-1))   Dzz[i] = .5*(DD[(z+1)*w*h+y*w+x]  - DD[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dzz[i] =     DD[i]                - DD[(z-1)*w*h+y*w+x];

        Dzz[i] *= (sig*sig); // correct for scaling
    }

    cout << "Dzz["<< sig <<"]," << flush;

    // Dy
    for (long i = 0; i < (w*h*l); ++i) {
       x = i%w; z = i/(w*h); y = i/w-z*h;
       if (y==0)            DD[i] =     F[z*w*h+(y+1)*w+x] - F[i];
       else if (y<(h-1))    DD[i] = .5*(F[z*w*h+(y+1)*w+x]  - F[z*w*h+(y-1)*w+x]);
       else if (y==(h-1))   DD[i] =     F[i]                - F[z*w*h+(y-1)*w+x];
    }

    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        // Dyy
        if (y==0)           Dyy[i] =     DD[z*w*h+(y+1)*w+x] - DD[i];
        else if (y<(h-1))   Dyy[i] = .5*(DD[z*w*h+(y+1)*w+x] - DD[z*w*h+(y-1)*w+x]);
        else if (y==(h-1))  Dyy[i] =     DD[i]               - DD[z*w*h+(y-1)*w+x];

        Dyy[i] *= (sig*sig); // correct for scaling
        // Dyz
        if (z==0)           Dyz[i] =     DD[(z+1)*w*h+y*w+x] - DD[i];
        else if (z<(l-1))   Dyz[i] = .5*(DD[(z+1)*w*h+y*w+x] - DD[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dyz[i] =     DD[i]               - DD[(z-1)*w*h+y*w+x];

        Dyz[i] *= (sig*sig); // correct for scaling
    }

    cout << "Dyy["<< sig <<"]," << flush;
    cout << "Dyz["<< sig <<"]," << flush;

    // Dx
    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        if (x==0)           DD[i] =      F[z*w*h+y*w+(x+1)]  - F[i];
        else if (x<(w-1))   DD[i] =  .5*(F[z*w*h+y*w+(x+1)]  - F[z*w*h+y*w+(x-1)]);
        else if (x==(w-1))  DD[i] =      F[i]                - F[z*w*h+y*w+(x-1)];
    }

    delete [] F;   F = 0; // not necessary after the last axis derivative is done

    for (long i = 0; i < (w*h*l); ++i) {
        x = i%w; z = i/(w*h); y = i/w-z*h;
        // Dxx
        if (x==0)           Dxx[i] =      DD[z*w*h+y*w+(x+1)] - DD[i];
        else if (x<(w-1))   Dxx[i] =  .5*(DD[z*w*h+y*w+(x+1)] - DD[z*w*h+y*w+(x-1)]);
        else if (x==(w-1))  Dxx[i] =      DD[i]               - DD[z*w*h+y*w+(x-1)];

        Dxx[i] *= (sig*sig); // correct for scaling
        // Dxy
        if (y==0)           Dxy[i] =     DD[z*w*h+(y+1)*w+x] - DD[i];
        else if (y<(h-1))   Dxy[i] = .5*(DD[z*w*h+(y+1)*w+x] - DD[z*w*h+(y-1)*w+x]);
        else if (y==(h-1))  Dxy[i] =     DD[i]               - DD[z*w*h+(y-1)*w+x];

        Dxy[i] *= (sig*sig); // correct for scaling
        // Dxz
        if (z==0)           Dxz[i] =     DD[(z+1)*w*h+y*w+x] - DD[i];
        else if (z<(l-1))   Dxz[i] = .5*(DD[(z+1)*w*h+y*w+x] - DD[(z-1)*w*h+y*w+x]);
        else if (z==(l-1))  Dxz[i] =     DD[i]               - DD[(z-1)*w*h+y*w+x];

        Dxz[i] *= (sig*sig); // correct for scaling
    }

    delete [] DD;   DD = 0;

    t2 = clock();
    cout << "Dxx["<< sig <<"]," << flush;
    cout << "Dxy["<< sig <<"]," << flush;
    cout << "Dxz["<< sig <<"]... " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

}

void Frangi::frangi2d(unsigned char* I, int w, int h, int l, float* J, float& Jmin, float& Jmax, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz) {

//    generate_2d_directions(ndirs2d, Vxyz);
//    if (true) return;

    cout << "Frangi 2d... I["<<w<<","<<h<<","<<l<<"] sig = "<< flush;

    float* Dxx   = new float[w*h*l];
    float* Dxy   = new float[w*h*l];
    float* Dyy   = new float[w*h*l];

    float tmp, mag;
    float v2x, v2y, v1x, v1y;
    float mu1, mu2;
    float Lambda1, Lambda2;
    float Vecx, Vecy, Vecn;
    float Rb, S2, Ifiltered;
    int val;

    float beta  = 2*pow(BetaOne,2);
    float c     = 2*pow(BetaTwo,2);

    Jmin = FLT_MAX;
    Jmax = -FLT_MAX;
    for (int si = 0; si < sig.size(); ++si) {

        cout << sig[si] << ((si<sig.size()-1)?", ":"") << flush;

        hessian2d(I, w, h, sig[si], Dyy, Dxy, Dxx);

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
                if (J[i]<Jmin) Jmin = J[i];
                if (J[i]>Jmax) Jmax = J[i];

                Vecn = sqrt(Vecx*Vecx+Vecy*Vecy);

                val = round((((Vecx/Vecn)+1)/2)*255.0);
                val = (val<0)?0:(val>255)?255:val;
                Vx[i] = (unsigned char) val;

                val = round((((Vecy/Vecn)+1)/2)*255.0);
                val = (val<0)?0:(val>255)?255:val;
                Vy[i] = (unsigned char) val;

                Vz[i] = (unsigned char) 0;

            }
            else {
                if (Ifiltered>J[i]) {
                    J[i] = Ifiltered; // keep maximum filter response
                    if (J[i]<Jmin) Jmin = J[i];
                    if (J[i]>Jmax) Jmax = J[i];

                    Vecn = sqrt(Vecx*Vecx+Vecy*Vecy);

                    val = round((((Vecx/Vecn)+1)/2)*255.0);
                    val = (val<0)?0:(val>255)?255:val;
                    Vx[i] = (unsigned char) val;

                    val = round((((Vecy/Vecn)+1)/2)*255.0);
                    val = (val<0)?0:(val>255)?255:val;
                    Vy[i] = (unsigned char) val;

                    Vz[i] = (unsigned char) 0;
                }
            }
        }
    } // sigma loop

    delete [] Dyy; Dyy = 0;
    delete [] Dxx; Dxx = 0;
    delete [] Dxy; Dxy = 0;
}

void Frangi::hessian2d(unsigned char* I, int w, int h, float sig, float* Dyy, float* Dxy, float* Dxx) {

    float* F     = new float[w*h];
    imgaussian(I, w, h, sig, F);

    int x,y;

    // Dy
    float* DD = new float[w*h];
    for (long i = 0; i < (w*h); ++i) {
       x = i%w; y = i/w;
       if (y==0)            DD[i] =     F[(y+1)*w+x] - F[i];
       else if (y<(h-1))    DD[i] = .5*(F[(y+1)*w+x] - F[(y-1)*w+x]);
       else if (y==(h-1))   DD[i] =     F[i]         - F[(y-1)*w+x];
    }

    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        // Dyy
        if (y==0)           Dyy[i] =     DD[(y+1)*w+x] - DD[i];
        else if (y<(h-1))   Dyy[i] = .5*(DD[(y+1)*w+x] - DD[(y-1)*w+x]);
        else if (y==(h-1))  Dyy[i] =     DD[i]         - DD[(y-1)*w+x];

        Dyy[i] *= (sig*sig); // correct for scaling
    }

    // Dx
    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        if (x==0)           DD[i] =      F[y*w+(x+1)]  - F[i];
        else if (x<(w-1))   DD[i] =  .5*(F[y*w+(x+1)]  - F[y*w+(x-1)]);
        else if (x==(w-1))  DD[i] =      F[i]          - F[y*w+(x-1)];
    }

    delete [] F;   F = 0; // not necessary after the last axis derivative is done

    for (long i = 0; i < (w*h); ++i) {
        x = i%w; y = i/w;
        // Dxx
        if (x==0)           Dxx[i] =      DD[y*w+(x+1)] - DD[i];
        else if (x<(w-1))   Dxx[i] =  .5*(DD[y*w+(x+1)] - DD[y*w+(x-1)]);
        else if (x==(w-1))  Dxx[i] =      DD[i]         - DD[y*w+(x-1)];

        Dxx[i] *= (sig*sig); // correct for scaling
        // Dxy
        if (y==0)           Dxy[i] =     DD[(y+1)*w+x] - DD[i];
        else if (y<(h-1))   Dxy[i] = .5*(DD[(y+1)*w+x] - DD[(y-1)*w+x]);
        else if (y==(h-1))  Dxy[i] =     DD[i]         - DD[(y-1)*w+x];

        Dxy[i] *= (sig*sig); // correct for scaling
    }

    delete [] DD;   DD = 0;

}

void Frangi::imgaussian(unsigned char * I, int w, int h, float sig, float* F) {

    long i0, i1;

    // gaussian filter is separated into 1D Gaussian kernel Gxy[.] that will be used for filtering along x and y
    int Lxy = ceil(3*sig);
//    float Gxy[2*Lxy+1];
    vector<float> Gxy(2*Lxy+1);
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

//    float Gxy[2*Lxy+1];
    vector<float> Gxy(2*Lxy+1);

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
//    float Gz[2*Lz  +1];
    vector<float> Gz(2*Lz  +1);
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
                    i1 = z1*w*h+y*w+x;
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

void Frangi::imgaussian(unsigned char * I, int w, int h, int l, float sig) {

    long i0, i1;

    // filtering is in xy plane only
    // gaussian filter is separated into 1D Gaussian kernels Gy[Gx[.]]
    int Lxy = ceil(3*sig);

//    float Gxy[2*Lxy+1];
    vector<float> Gxy(2*Lxy+1);

    float Gnorm = 0;
    for (int i = -Lxy; i <= Lxy; ++i) {
        Gxy[i+Lxy] = exp(-(i*i)/(2*sig*sig));
        Gnorm += Gxy[i+Lxy];
    }

    for (int i = 0; i < (2*Lxy+1); ++i) {
        Gxy[i] /= Gnorm;
    }

    // Gx gaussian smoothing of the byte8 image array stored in I[.], result in K[.]
    float * K = new float[w*h*l];
    for (int y = 0; y < h; ++y) {
        for (int z = 0; z < l; ++z) {
            // x index clamping
            for (int x = 0; x < min(Lxy,w); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+clamp(x1,0,w-1);
                    K[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+x1;
                    K[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }
            // x index clamping
            for (int x = max(w-Lxy,Lxy); x < w; ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = 0;
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    i1 = z*w*h+y*w+clamp(x1,0,w-1);
                    K[i0] += I[i1] * Gxy[x1-x+Lxy];
                }
            }

        }
    }

    // Gy gaussian smoothing of the float image array stored in K[.], result back to I[.]
    for (int x = 0; x < w; ++x) {
        for (int z = 0; z < l; ++z) {
            // y index clamping
            for (int y = 0; y < min(Lxy,h); ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                    I[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = 0;
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    i1 = z*w*h+y1*w+x;
                    I[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
            // y index clamping
            for (int y = max(h-Lxy,Lxy); y < h; ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = 0;
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                    I[i0] += K[i1] * Gxy[y1-y+Lxy];
                }
            }
        }
    }

    delete [] K; K = 0;

}

void Frangi::imerode(unsigned char* I, int w, int h, int l, float rad, unsigned char* E) {

    long i0, i1;

    // erosion filter is separated into 1D erosions
    // this method does erosion in xy plane only
    int Lxy = ceil(rad);

    unsigned char* K = new unsigned char[w*h*l];
    // x erosion of the byte8 image array stored in I[.], result in K[.]
    for (int y = 0; y < h; ++y) {
        for (int z = 0; z < l; ++z) {
            // x index clamping
            for (int x = 0; x < min(Lxy,w); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]<K[i0]) K[i0] = I[i1];
                    }
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+x1;
                        if (I[i1]<K[i0]) K[i0] = I[i1];
                    }
                }
            }
            // x index clamping
            for (int x = max(w-Lxy,Lxy); x < w; ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]<K[i0]) K[i0] = I[i1];
                    }
                }
            }

        }
    }

    // y erosion of the byte8 image array stored in K[.], result in E[.]
    for (int x = 0; x < w; ++x) {
        for (int z = 0; z < l; ++z) {
            // y index clamping
            for (int y = 0; y < min(Lxy,h); ++y) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+y1*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }
                }
            }
            // y index clamping
            for (int y = max(h-Lxy,Lxy); y < h; ++y) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }
                }
            }
        }
    }

    delete [] K; K = 0;

}

void Frangi::imerode(unsigned char* I, int w, int h, int l, float rad, float zdist, unsigned char* E) {

    long i0, i1;
    // erosion filter is separated into 1D erosions
    int Lxy = ceil(rad);

    // x erosion of the byte8 image array stored in I[.], result in E[.]
    for (int y = 0; y < h; ++y) {
        for (int z = 0; z < l; ++z) {
            // x index clamping
            for (int x = 0; x < min(Lxy,w); ++x) {
                i0 = z*w*h+y*w+x;
                E[i0] = I[i0]; // I[z*w*h+y*w+clamp(x-Lxy,0,w-1)];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]<E[i0]) E[i0] = I[i1];
                    }
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = z*w*h+y*w+x;
                E[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+x1;
                        if (I[i1]<E[i0]) E[i0] = I[i1];
                    }
                }
            }
            // x index clamping
            for (int x = max(w-Lxy,Lxy); x < w; ++x) {
                i0 = z*w*h+y*w+x;
                E[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]<E[i0]) E[i0] = I[i1];
                    }
                }
            }

        }
    }

    // y erosion of the byte8 image array stored in E[.], result in K[.]
    unsigned char* K = new unsigned char[w*h*l];
    for (int x = 0; x < w; ++x) {
        for (int z = 0; z < l; ++z) {
            // y index clamping
            for (int y = 0; y < min(Lxy,h); ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = E[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (E[i1]<K[i0]) K[i0] = E[i1];
                    }
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = E[i0];
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+y1*w+x;
                        if (E[i1]<K[i0]) K[i0] = E[i1];
                    }
                }
            }
            // y index clamping
            for (int y = max(h-Lxy,Lxy); y < h; ++y) {
                i0 = z*w*h+y*w+x;
                K[i0] = E[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (E[i1]<K[i0]) K[i0] = E[i1];
                    }
                }
            }
        }
    }

    // z erosion of the byte8 image array stored in K[.], result in E[.]
    if (l==1) {
        for (long ii = 0; ii < (w*h*l); ++ii) E[ii] = K[ii]; // no erosion needed
    }
    else {

    int Lz  = ceil(rad/zdist);

    for (int x = 0; x < w; ++x) {
        for (int y = 0; y < h; ++y) {
            // z index clamping
            for (int z = 0; z < min(Lz,l); ++z) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    if (z1!=z) {
                        i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }

                }
            }
            // no clamp
            for (int z = Lz; z < (l-Lz); ++z) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    if (z1!=z) {
                        i1 = z1*w*h+y*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }
                }
            }
            // z index clamping
            for (int z = max(l-Lz,Lz); z < l; ++z) {
                i0 = z*w*h+y*w+x;
                E[i0] = K[i0];
                for (int z1 = (z-Lz); z1 <= (z+Lz); ++z1) {
                    if (z1!=z) {
                        i1 = clamp(z1,0,l-1)*w*h+y*w+x;
                        if (K[i1]<E[i0]) E[i0] = K[i1];
                    }
                }
            }

        }
    }
    } // l>1

    delete [] K; K = 0;

}

void Frangi::imdilate(unsigned char* I, int w, int h, int l, float rad) {

    long i0, i1;

    // dilatation is separated into 1D dilatations, dilates only in xy plane
    int Lxy = ceil(rad);

    unsigned char* K = new unsigned char[w*h*l];
    // x dilatation of the byte8 image array stored in I[.], result in K[.]
    for (int y = 0; y < h; ++y) {
        for (int z = 0; z < l; ++z) {
            // x index clamping
            for (int x = 0; x < min(Lxy,w); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]>K[i0]) K[i0] = I[i1];
                    }
                }
            }
            // no clamp
            for (int x = Lxy; x < (w-Lxy); ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+x1;
                        if (I[i1]>K[i0]) K[i0] = I[i1];
                    }
                }
            }
            // x index clamping
            for (int x = max(w-Lxy,Lxy); x < w; ++x) {
                i0 = z*w*h+y*w+x;
                K[i0] = I[i0];
                for (int x1 = x-Lxy; x1 <= x+Lxy; ++x1) {
                    if (x1!=x) {
                        i1 = z*w*h+y*w+clamp(x1,0,w-1);
                        if (I[i1]>K[i0]) K[i0] = I[i1];
                    }
                }
            }

        }
    }

    // y dilatation of the byte8 image array stored in K[.], result back in I[.]
    for (int x = 0; x < w; ++x) {
        for (int z = 0; z < l; ++z) {
            // y index clamping
            for (int y = 0; y < min(Lxy,h); ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = K[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (K[i1]>I[i0]) I[i0] = K[i1];
                    }
                }
            }
            // no clamp
            for (int y = Lxy; y < (h-Lxy); ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = K[i0];
                for (int y1 = (y-Lxy); y1 <= (y+Lxy); ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+y1*w+x;
                        if (K[i1]>I[i0]) I[i0] = K[i1];
                    }
                }
            }
            // y index clamping
            for (int y = max(h-Lxy,Lxy); y < h; ++y) {
                i0 = z*w*h+y*w+x;
                I[i0] = K[i0];
                for (int y1 = y-Lxy; y1 <= y+Lxy; ++y1) {
                    if (y1!=y) {
                        i1 = z*w*h+clamp(y1,0,h-1)*w+x;
                        if (K[i1]>I[i0]) I[i0] = K[i1];
                    }
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
