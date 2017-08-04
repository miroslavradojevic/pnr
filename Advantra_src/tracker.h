#ifndef TRACKER_H
#define TRACKER_H
#include "node.h"
#include "seed.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>
#include <queue>
#include <v3d_interface.h>

struct X {
    X() {x=0; y=0; z=0;     vx=0; vy=0; vz=0;   sig=0; corr=0; w=0;}
    X(float _x,float _y,float _z,float _vx,float _vy,float _vz,float _w) : x(_x),y(_y),z(_z),vx(_vx),vy(_vy),vz(_vz),w(_w) {corr=0; sig=0;}
    float x, y, z, vx, vy, vz, w, corr, sig;
};

struct X_est { // estimation, similar components as in Node class instance
    X_est() {x=0; y=0; z=0;     vx=0; vy=0; vz=0;   sig=0; corr=0;}
    X_est(float _x,float _y,float _z,float _vx,float _vy,float _vz,float _sig,float _corr) : x(_x),y(_y),z(_z),vx(_vx),vy(_vy),vz(_vz),sig(_sig),corr(_corr) {}
  float x, y, z, vx, vy, vz, sig, corr;
};

struct offxyz { // offset x,y,z
    offxyz(int _x, int _y, int _z) : x(_x),y(_y),z(_z) {}
    int x, y, z;
};

struct offvuw { // offset v,u,w, template offsets
    offvuw(int _v, int _u, int _w) : v(_v),u(_u),w(_w) {}
    int v, u, w;
};

class Tracker
{
public:

    static int    ndirs2d;
    static int    ndirs3d;

    static int NH;
    static int NH5x5[25][3];
    static int NH4x4_ul[16][3];
    static int NH4x4_ur[16][3];
    static int NH4x4_dl[16][3];
    static int NH4x4_dr[16][3];
    static int NH3x3[9][3];
    static int NH2x2_ul[4][3];
    static int NH2x2_ur[4][3];
    static int NH2x2_dl[4][3];
    static int NH2x2_dr[4][3];

    bool        is2d;
    int         ndir;
    int         sz;
    int         step;
    int         npcles;
    int         niter;
    float       kappa;

    vector<float>   sig;
    // v=(vx,vy,vz), u=(ux,uy,uz) and w=(wx,wy,wz) are aligned with vx,vy,vz of the particle orientation
    vector< vector<float> >     model_wgt;  // sig.size() array of vectors for the template weights
    vector< vector<float> >     model_img;  // sig.size() array of vectors for image values
    vector<float>               model_avg;  // sig.size() array of vectors for the template averages
    vector< vector<offvuw> >    model_vuw;  // sig.size() array of vectors for unit offsets along orthogonals
    float       ux, uy, uz; // orthogonal 1
    float       wx, wy, wz; // orhtogonal 2
    float       zDist;

    bool verbose;
    float znccth;
    int nodes_pp;
    float Kc;
    float neff_ratio;

    float** p;              // sz * 3 // predictions offsets
    float*  cws;            // sz
    float*  d;              // sz    (distances if z scaled down)
    float*  d0;             // sz    (distances if z not scaled down)
    float** u;              // sz * 3
    float** w;              // ndir * sz
    float** w_cws;          // ndir * sz (for sampling)
    float*  w0;             // sz
    float*  w0_cws;         // sz (for sampling)
    float** v;              // ndir * 3

    float*  prior;          // priors array
    float*  lhood;          // likelihoods array
    float*  res_csw;        // cummulative sum of weights array (for resampling at one iteration)
    X       x0;             // initial particle of the track
    X**     xfilt;          // filtered particles niter * npcles
    int**   idxres;         // indexes of the resamled particles niter * npcles

    // offsets used to fill up the nodemap for different radii, round(radius) gives the index
    vector< vector<offxyz> > off3; // voxel offsets for for different indexes/radiuses

    float* neff;
    X_est* xc; // centroid estimation, niter sized array (struct devoted to store the estimation)
//    int * tagg;
//    int * ovlp; // ovlp[i]=0 seed, >0 there is overlap with node idx

    int ti_limit; // last iteration index reached

    // ms clustering
    float** conv; // ms convergence
    int* labels;  // clustering labels
    bool* checked; // which labels were examined
    vector<int>* nbridxs;    // particle indexes of the neighbouring particles
    vector<int>   clst_idx; // particle indexes per cluster
    vector<float> clst_csw;  // cummulative sum of weights for sampling from the cluster
    int MAXITER;
    float EPSILON2;
    float KRAD;

    Tracker(
            float   _sigBeg,
            float   _sigStp,
            float   _sigEnd,
            int     _step,
            int     _npcles,
            int     _niter,
            float   _kappa,
            bool    _is2d,
            float   _znccth,
            float   _Kc,
            float   _neff_ratio,
            float   _zdist,
            int    _nodes_pp);

    ~Tracker();

    void generate_directions(bool is2D, float** vxyz);

    float interp(float _x, float _y, float _z,
                 unsigned char * img, int width, int height, int length);

    double bessi(int n, double x);

    static double bessi0(double x );

    static double bessi1(double x);

//    void trackNew(int _x, int _y, int _z, unsigned char * _img, vector<Node> & _nodelist, int _w, int _h, int _l, int * _nmap, int ovlp_window);

    void trackPos(seed _seed0, unsigned char* _img, vector<Node>& _nodelist, int _w, int _h, int _l, bool* _checked, unsigned char* _trc_den);

    void trackNeg(seed _seed0, unsigned char* _img, vector<Node>& _nodelist, int _w, int _h, int _l, bool* _checked, unsigned char* _trc_den);

    // sampling and resampling -- similar technique used in iter0() and iterI() for prediction
    static void sampleN(vector<float> _csw, int _N, vector<int>& _sampled_idx);

    float zncc(float _x,  float _y,   float _z,
               float _vx, float _vy,  float _vz,
               bool _return_avg,
               unsigned char * img, int img_w, int img_h, int img_l, float & _gcsstd_value);

    float zncc1(X _xp, unsigned char * img, int img_w, int img_h, int img_l, float & _gcsstd_value);

    float zncc2(X_est _xp, unsigned char * img, int img_w, int img_h, int img_l, float & _gcsstd_value);

    void extract(int _i);

    void clustering(float _dist); // labels nbridxs

    void ms(int _i, float _krad); // conv

    void runOne(float * _curr_v, float * _new_v, X * _xf, int _xfLen, float _krad);

    int getdirection(float _vx, float _vy, float _vz);

    int overlaps(float _x, float _y, float _z, int * _nodemap, int _w, int _h, int _l);

    int overlaps(X_est _x_est, int * _nodemap, int _w, int _h, int _l);

    void fillAllow(vector<Node> _nodelist, int * _map, int _w, int _h, int _l, bool _ow);

    void fill(int _nidx, vector<Node> _nodelist, int * _map, int _w, int _h, int _l, bool _ow);

    static void crossprod(float a1, float a2, float a3, float b1, float b2, float b3, float * v);

    float d2(float a1x, float a1y, float a1z,
             float a2x, float a2y, float a2z,
             float a0x, float a0y, float a0z);

    void export_off3(string savepath);
    void export_model(string savepath, bool directed=false);
//    void export_model(float vx, float vy, float vz, string savepath);
    void export_track(string savepath); // track swc reconstruction
    void export_trackcorr(string savepath); // track corr log

    void fill(X_est _xc, int _tag, int * _nmap, int _w, int _h, int _l);

    static void check1x1(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l);
    static void add1x1(float _x, float _y, float _z, unsigned char* _trc_den, int _w, int _h, int _l);

    static void check2x2(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l);

    static void check3x3(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l);

    static void check4x4(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l);

    static void check5x5(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l);

private:

    bool iter0New(unsigned char * _img, int _w, int _h, int _l); // vector<Node> & _nodelist,

    bool iterINew(int _i, unsigned char * _img, int _w, int _h, int _l); // , vector<Node> & _nodelist

};

#endif // TRACKER_H

//    void assoc0(int * _nodemap, int _w, int _h, int _l, vector<Node> & _nodelist);
//    void assocI(int _i, int * _nodemap, int _w, int _h, int _l, vector<Node> & _nodelist);
//    void nodemapFill(X_est _x_est, int _ntag, int * _nodemap, int _w, int _h, int _l, TagBuffer<int> _allow);
//    bool iterI(int _i, unsigned char * _img, vector<Node> & _nodelist, int * _nodemap, int * _trackmap, int _w, int _h, int _l);
//    float l1804(X a, X b, unsigned char * img, int w, int h, int l);
//    bool iter0(unsigned char * _img, vector<Node> & _nodelist, int * _nodemap, int * _trackmap, int _w, int _h, int _l);
//    void track(int _x, int _y, int _z, unsigned char * _img, vector<Node> & _nodelist, int * _nodemap, int * _trackmap, int _w, int _h, int _l);
