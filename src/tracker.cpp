#include "tracker.h"
#include "node.h"
#include <iostream>
#include <fstream>

#define ACC 40.0
#define BIGNO 1.0e10
#define BIGNI 1.0e-10

#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a)>(b)?(a):(b))

#ifdef MIN
#undef MIN
#endif
#define MIN(a, b) ((a)<(b)?(a):(b))

int     Tracker::ndirs2d   = 30;
int     Tracker::ndirs3d   = 50;

int     Tracker::NH = 1; // neighbourhood size 1,2,3... (1x1,2x2,3x3) used for the seed suppression (seeds are suppressed by other seeds or by trace)
int     Tracker::NH5x5[25][3] = {
                                {-2,-2,0},{-2,-1,0},{-2,0,0},{-2,1,0},{-2,2,0},
                                {-1,-2,0},{-1,-1,0},{-1,0,0},{-1,1,0},{-1,2,0},
                                {0,-2,0},{0,-1,0},{0,0,0},{0,1,0},{0,2,0},
                                {1,-2,0},{1,-1,0},{1,0,0},{1,1,0},{1,2,0},
                                {2,-2,0},{2,-1,0},{2,0,0},{2,1,0},{2,2,0}
                                };
int     Tracker::NH4x4_ul[16][3]= {
                                    {-2,-2,0},{-2,-1,0},{-2,0,0},{-2,1,0},
                                    {-1,-2,0},{-1,-1,0},{-1,0,0},{-1,1,0},
                                    {0,-2,0},{0,-1,0},{0,0,0},{0,1,0},
                                    {1,-2,0},{1,-1,0},{1,0,0},{1,1,0}
                                };

int     Tracker::NH4x4_ur[16][3]= {
                                    {-1,-2,0},{-1,-1,0},{-1,0,0},{-1,1,0},
                                    {0,-2,0},{0,-1,0},{0,0,0},{0,1,0},
                                    {1,-2,0},{1,-1,0},{1,0,0},{1,1,0},
                                    {2,-2,0},{2,-1,0},{2,0,0},{2,1,0}
                                };

int     Tracker::NH4x4_dl[16][3]= {
                                    {-2,-1,0},{-2,0,0},{-2,1,0},{-2,2,0},
                                    {-1,-1,0},{-1,0,0},{-1,1,0},{-1,2,0},
                                    {0,-1,0},{0,0,0},{0,1,0},{0,2,0},
                                    {1,-1,0},{1,0,0},{1,1,0},{1,2,0}
                                };

int     Tracker::NH4x4_dr[16][3]= {
                                    {-1,-1,0},{-1,0,0},{-1,1,0},{-1,2,0},
                                    {0,-1,0},{0,0,0},{0,1,0},{0,2,0},
                                    {1,-1,0},{1,0,0},{1,1,0},{1,2,0},
                                    {2,-1,0},{2,0,0},{2,1,0},{2,2,0}
                                };

int     Tracker::NH3x3[9][3]    = {{-1,-1,0},{0,-1,0},{1,-1,0},{-1,0,0},{0,0,0},{1,0,0},{-1,1,0},{0,1,0},{1,1,0}};
int     Tracker::NH2x2_ul[4][3] = {{0,0,0},{0,-1,0},{-1,0,0},{-1,-1,0}};
int     Tracker::NH2x2_ur[4][3] = {{0,0,0},{0,-1,0},{1,0,0},{1,-1,0}};
int     Tracker::NH2x2_dl[4][3] = {{0,0,0},{0,1,0},{-1,0,0},{-1,1,0}};
int     Tracker::NH2x2_dr[4][3] = {{0,0,0},{0,1,0},{1,0,0},{1,1,0}};

using namespace std;

template<typename T>
T clamp(T x, T x1, T x2) {T xC = (x<x1)?x1:x; return (xC>x2)?x2:xC;}

float clampf(float x, float x1, float x2) {float xC = (x<x1)?x1:x; return (xC>x2)?x2:xC;}

Tracker::Tracker(   vector<float> _sigs,
                    int    _step,
                    int    _npcles,
                    int    _niter,
                    float  _kappa,
                    bool   _is2d,
                    float  _znccth,
                    float  _Kc,
                    float  _neff_ratio,
                    float  _zdist,
                    int    _nodes_pp) {

    sig = _sigs;

    step = _step;
    npcles = _npcles;
    niter = _niter;
    kappa = _kappa;
    is2d = _is2d;
    ndir = (is2d)?ndirs2d:ndirs3d;
    zDist = _zdist;

    znccth = _znccth;
    nodespervol = _nodes_pp;
    Kc = _Kc;
    neff_ratio = _neff_ratio;

    verbose = false;

    // -----------------------------
    // clear model_* structures used
    model_avg.clear();
    for (int i = 0; i < model_wgt.size(); ++i) model_wgt[i].clear(); model_wgt.clear();
    for (int i = 0; i < model_img.size(); ++i) model_img[i].clear(); model_img.clear();
    for (int i = 0; i < model_vuw.size(); ++i) model_vuw[i].clear(); model_vuw.clear();

    for (int i = 0; i < sig.size(); ++i) {

        model_avg.push_back(0);
        vector<float> wgt;
        model_wgt.push_back(wgt);
        vector<float> img;
        model_img.push_back(img);
        vector<offvuw> vuw;
        model_vuw.push_back(vuw);

        // indexing differs in 2d and 3d
        if (is2d) {

            int V2 = round(1*sig[i]); // steps in the direciton aligned with vx,vy,vz
            int U2 = round(3*sig[i]); // orthogonal

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    float value = exp(-(uu*uu)/(2*pow(sig[i],2)));
                    model_wgt[i].push_back(value);
                    model_img[i].push_back(0); // just allocate storage for the image values
                    offvuw t(vv, uu, 0);
                    model_vuw[i].push_back(t);
                    model_avg[i] += value;
                }
            }
        }
        else {

            int V2 = round(1*sig[i]);
            int U2 = round(3*sig[i]);
            int W2 = round(3*sig[i]);

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    for (int ww = -W2; ww <= W2; ++ww) {
                        float value = exp(-((uu*uu)+(ww*ww))/(2*pow(sig[i],2)));
                        model_wgt[i].push_back(value);
                        model_img[i].push_back(0); // just allocate
                        offvuw t(vv, uu, ww);
                        model_vuw[i].push_back(t);
                        model_avg[i] += value;
                    }
                }
            }
        }

        model_avg[i] /= model_wgt[i].size(); // for each sigma

    }

//    vector< vector<float> >     model2_wgt;
//    vector< vector<float> >     model2_img;
//    vector<float>               model2_avg;
//    vector< vector<Pvuw> >      model2_vuw;

    // -----------------------------
    // clear model2_*
    model2_avg.clear();
    for (int i = 0; i < model2_wgt.size(); ++i) model2_wgt[i].clear(); model2_wgt.clear();
    for (int i = 0; i < model2_img.size(); ++i) model2_img[i].clear(); model2_img.clear();
    for (int i = 0; i < model2_vuw.size(); ++i) model2_vuw[i].clear(); model2_vuw.clear();

    model2_N = 12; // hardcoded, # samples per 3*sigma

    for (int i = 0; i < sig.size(); ++i) {

        model2_avg.push_back(0.0);
        vector<float> wgt;
        model2_wgt.push_back(wgt);
        vector<float> img;
        model2_img.push_back(img);
        vector<Pvuw> vuw;
        model2_vuw.push_back(vuw);

        // indexing differs in 2d and 3d
        if (is2d) {
            int V2 = round(1*sig[i]); // steps in the direciton aligned with vx,vy,vz
            int U2 = round(3*sig[i]); // orthogonal
            float Vs = (3.0*sig[i])/model2_N;
            Vs = (Vs<1.0)? 1.0 : Vs ;

            for (float vv = -V2; vv <= V2+FLT_MIN; vv+=Vs) {
                for (float uu = -U2; uu <= U2+FLT_MIN; uu+=Vs) {
                    float value = exp(-(uu*uu)/(2*pow(sig[i],2)));
                    model2_wgt[i].push_back(value);
                    model2_img[i].push_back(0.0); // just allocate storage for the image values
                    Pvuw t(vv, uu, 0);
                    model2_vuw[i].push_back(t);
                    model2_avg[i] += value;
                }
            }

        }
        else {
            int V2 = round(1*sig[i]);
            int U2 = round(3*sig[i]);
            int W2 = round(3*sig[i]);
            float Vs = (3.0*sig[i])/model2_N;
            Vs = (Vs<1.0)? 1.0 : Vs ;

            for (float vv = -V2; vv <= V2+FLT_MIN; vv+=Vs) {
                for (float uu = -U2; uu <= U2+FLT_MIN; uu+=Vs) {
                    for (float ww = -W2; ww <= W2+FLT_MIN; ww+=Vs) {
                        float value = exp(-((uu*uu)+(ww*ww))/(2*pow(sig[i],2)));
                        model2_wgt[i].push_back(value);
                        model2_img[i].push_back(0.0); // just allocate
                        Pvuw t(vv, uu, ww);
                        model2_vuw[i].push_back(t);
                        model2_avg[i] += value;
                    }
                }
            }

        }

        model2_avg[i] /= model2_wgt[i].size(); // for each sigma
    }

//    for (int i = 0; i < sig.size(); ++i) {
//        cout << "sig[" << i << "] = " << sig[i] << " | " << model2_vuw[i].size() << " offsets  |  " << model2_img[i].size() << endl;
//    }

//    vector<offvuw>              pattern_vuw; // offsets (pattern_vuw.size() corresponds to largest sigma)
//    vector<float>               pattern_img; // image values (pattern_img.size() corresponds to largest sigma)
//    vector< vector<float> >     pattern_wgt; // pattern weights (pattern_wgt.size() correponds to largest sigma)
//    vector< vector<int> >       pattern_sid; // sigma indexes per offset
//    vector<float>               pattern_avg; // average values (MODEL WEIGHTS)
//    vector<int>                 pattern_cnt; // number of offsets (image values) per each sig[]
//    vector<float>               pattern_ag;  // average for the IMAGE VALS per each sigma
//    vector<float>               pattern_corra; // partial result used to calculate correlation at sig[]
//    vector<float>               pattern_corrb; //
//    vector<float>               pattern_corrc; //

    pattern_avg.clear();
    pattern_cnt.clear();
    pattern_ag.clear();
    pattern_corra.clear();
    pattern_corrb.clear();
    pattern_corrc.clear();

    for (int i = 0; i < sig.size(); ++i) {
        pattern_ag.push_back(0.0);      // image values average per each sigma

        pattern_cnt.push_back(0);       // count number of offsets used per sigma
        pattern_avg.push_back(0.0);     // model values average per sigma
        pattern_corra.push_back(0.0);   // correlation auxiliaries
        pattern_corrb.push_back(0.0);   //
        pattern_corrc.push_back(0.0);   //
    }

    pattern_vuw.clear();
    pattern_img.clear();
    for (int i = 0; i < pattern_wgt.size(); ++i) pattern_wgt[i].clear(); pattern_wgt.clear();
    for (int i = 0; i < pattern_sid.size(); ++i) pattern_sid[i].clear(); pattern_sid.clear();

    if (is2d) {

        int V2 = round(1*sig.back()); // sig[] is sorted so V2, U2 define outer offset boundary
        int U2 = round(3*sig.back()); //

        for (int vv = -V2; vv <= V2; ++vv) {
            for (int uu = -U2; uu <= U2; ++uu) {

                offvuw vuw(vv, uu, 0);
                pattern_vuw.push_back(vuw);

                pattern_img.push_back(0.0); // allocate

                vector<float> wgt;
                vector<int> sid;
                for (int sig_idx = 0; sig_idx < sig.size(); ++sig_idx) {
                    if ( abs(vv)<=round(1*sig[sig_idx]) && abs(uu)<=round(3*sig[sig_idx]) ) {
                        sid.push_back(sig_idx);
                        float weight = exp(-(uu*uu)/(2*pow(sig[sig_idx],2)));
                        wgt.push_back(weight);
                        pattern_cnt[sig_idx] += 1;
                        pattern_avg[sig_idx] += weight;
                    }
                }
                pattern_wgt.push_back(wgt);
                pattern_sid.push_back(sid);
            }
        }

    }
    else { // 3d
        int V2 = round(1*sig.back());
        int U2 = round(3*sig.back());
        int W2 = round(3*sig.back());

        for (int vv = -V2; vv <= V2; ++vv) {
            for (int uu = -U2; uu <= U2; ++uu) {
                for (int ww = -W2; ww <= W2; ++ww) {

                    offvuw vuw(vv, uu, ww);
                    pattern_vuw.push_back(vuw);

                    pattern_img.push_back(0.0); // allocate

                    vector<float> wgt;
                    vector<int> sid;
                    for (int sig_idx = 0; sig_idx < sig.size(); ++sig_idx) {
                        if ( abs(vv)<=round(1*sig[sig_idx]) && abs(uu)<=round(3*sig[sig_idx]) && abs(ww)<=round(3*sig[sig_idx]) ) {
                            sid.push_back(sig_idx);
                            float weight = exp(-((uu*uu)+(ww*ww))/(2*pow(sig[sig_idx],2)));
                            wgt.push_back(weight);
                            pattern_cnt[sig_idx] += 1;
                            pattern_avg[sig_idx] += weight;
                        }
                    }
                    pattern_wgt.push_back(wgt);
                    pattern_sid.push_back(sid);
                }
            }
        }
    }

    for (int k = 0; k < pattern_cnt.size(); ++k) {
        pattern_avg[k] /= pattern_cnt[k];
//        cout << "pattern_cnt[" << k << "] = " << pattern_cnt[k] << endl;
//        cout << "pattern_avg[" << k << "] = " << pattern_avg[k] << endl;
    }

//    cout << "pattern_vuw.size()=" << pattern_vuw.size() << endl;
//    for (int k = 0; k < pattern_vuw.size(); ++k) {
//        cout << "[" << pattern_vuw[k].v << "," << pattern_vuw[k].u << "," << pattern_vuw[k].w << "] " << flush;
//        cout << "[" << flush;
//        for (int k1 = 0; k1 < pattern_sid[k].size(); ++k1) cout << pattern_sid[k][k1] << "(" << sig[pattern_sid[k][k1]] << ")  " << flush;
//        cout << "] " << flush;
//        cout << "{" << flush;
//        for (int k1 = 0; k1 < pattern_wgt[k].size(); ++k1) cout << pattern_wgt[k][k1] << "   " << flush;
//        cout << "}"<< endl;
//    }

    // ---- x filtered initialization ----
    xfilt = new X*[niter];
    for (int i = 0; i < niter; ++i) xfilt[i] = new X[npcles];

    idxres = new int*[niter];
    for (int i = 0; i < niter; ++i) idxres[i] = new int[npcles];

    for (int i = 0; i < niter; ++i) {
        for (int j = 0; j < npcles; ++j) {
            idxres[i][j] = -INT_MAX; // initialize with false index
        }
    }

    prior = new float[npcles];
    lhood = new float[npcles];
    res_csw = new float[npcles]; // allocate cummulative sum of weights for resampling

    // these could have been static allocations and called in track(),
    // then they would be allocated each time track is called,
    // while this way they are allocated once, size of these arrays does not prevent them from being statically allocated,
    // this way they have to be deleted and they're initialized with zeros
    xc = new X_est[niter];
//    tagg = new int[niter];
//    ovlp = new int[niter];
    neff = new float[niter];

    // prediction matrices
    vector<int> px, py, pz;

    int stpRg = 2 * step;

    for (int dx = -stpRg; dx <= stpRg; ++dx) {
        for (int dy = -stpRg; dy <= stpRg; ++dy) {
            if (is2d) {
                if (dx*dx+dy*dy<=stpRg*stpRg && dx*dx+dy*dy>0) {
                    px.push_back(dx); py.push_back(dy); pz.push_back(0);
                }
            }
            else {
                for (int dz = -stpRg; dz <= stpRg; ++dz) {
                    if (dx*dx+dy*dy+dz*dz<=stpRg*stpRg && dx*dx+dy*dy+dz*dz>0) {
                        px.push_back(dx); py.push_back(dy); pz.push_back(dz);
                    }
                }
            }
        }
    }

    sz = px.size();

    p = new float*[sz];         // predicted locations
    d   = new float[sz];        // norm
    d0  = new float[sz];        //
    u = new float*[sz];         // unit vector for each
    cws = new float[sz];        // cummulative weight sum
    w0  = new float[sz];        //
    w0_cws = new float[sz];

    float w0sum = 0;

    for (int i = 0; i < sz; ++i) {

        p[i] = new float[3];
        p[i][0] = px[i];
        p[i][1] = py[i];
        p[i][2] = pz[i]/zDist; // predictions are scaled down with zdist

        d[i] = sqrt(p[i][0]*p[i][0] + p[i][1]*p[i][1] + p[i][2]*p[i][2]);
        d0[i] = sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);

        u[i] = new float[3];
        u[i][0] = p[i][0]/d[i];
        u[i][1] = p[i][1]/d[i];
        u[i][2] = p[i][2]/d[i];

        cws[i] = 0;

        w0[i] = exp(-pow(d[i],2)/(2*pow(step/3.0,2)));
        w0sum += w0[i];

    }

    for (int i = 0; i < sz; ++i) {
       w0[i] /= w0sum; // normalize w0[]
       w0_cws[i] = w0[i] + ((i==0)?0:w0_cws[i-1]);
    }

    // directions v, predefined directions used to pick the corresponding oriented prior matrix
    v = new float*[ndir];
    for (int i = 0; i < ndir; ++i) v[i] = new float[3];
    generate_directions(is2d, v);

    // w per directions v
    double rad, circ, val, dotp; // radial, polar distance

    w  = new float*[ndir];
    w_cws = new float*[ndir];

    for (int i = 0; i < ndir; ++i) {
        w[i] = new float[sz];
        w_cws[i] = new float[sz];
    }

    for (int i = 0; i < ndir; i++) {

        float wsum = 0;

        for (int j = 0; j < sz; j++) {

            rad = exp(-pow(d0[j]-step, 2)/(2*pow(step/3.0,2)));

            dotp = v[i][0]*u[j][0] + v[i][1]*u[j][1] + v[i][2]*u[j][2];
            dotp = (dotp>1)? 1 : (dotp<-1)? -1 : dotp;
            circ = exp(kappa * dotp) / (2.0*3.14* bessi0(kappa)); // von mises distribution

            w[i][j] = circ * rad;

            wsum += w[i][j];

        }

        for (int j = 0; j < sz; j++) {
            w[i][j] = w[i][j]/wsum;
            w_cws[i][j] = w[i][j] + ((j==0)?0:w_cws[i][j-1]);
        }

    }

    //-------------
    // offsets
    //-------------
    for (int i = 0; i < off3.size(); ++i) off3[i].clear();
    off3.clear();


    for (int i = 0; i <= ceil(4 * sig.back()); ++i) { // 4*sig would be a very large radius, to cover all possible radiuses, sig should be sorted so sig.back() is largest

        vector<offxyz> t;     // empty vector
        off3.push_back(t);    // now off3[i] can be appended and accessed

        int rxy = i;
        int rz  = ceil(i/zDist); // scale the layer range when establishing fill-up offsets

        if (i==0) {
            offxyz curroff(0, 0, 0);
            off3[i].push_back(curroff);
        }
        else {
            for (int dx = -rxy; dx <= rxy; ++dx) {
                for (int dy = -rxy; dy <= rxy; ++dy) {
                    for (int dz = -rz; dz <= rz; ++dz) {
                        double ins=
                                ((dx*dx)/(rxy*rxy))+
                                ((dy*dy)/(rxy*rxy))+
                                ((is2d || rz==0)?0:(dz*dz)/(rz*rz));

                        if (ins<=1.0) {
                            offxyz curroff(dx, dy, dz);
                            off3[i].push_back(curroff);
                        }
                    }
                }
            }
        }
    }

    //----------------
    // ms
    //----------------
    conv     = new float*[npcles];
    for (int i = 0; i < npcles; ++i) conv[i] = new float[3];
    labels      = new int[npcles];
    checked     = new bool[npcles];
    nbridxs     = new vector<int>[npcles];
    MAXITER     = INT_MAX;
    EPSILON2    = 0.000000001;
    KRAD        = 3; // could be dependent on diam too
}

Tracker::~Tracker(){

//    cout << "~Tracker()" << endl;

    for (int i = 0; i < sz; ++i) delete [] p[i];
    delete [] p; p = 0;

    delete [] d; d = 0;
    delete [] d0; d0 = 0;

    for (int i = 0; i < sz; ++i) delete [] u[i];
    delete [] u; u = 0;

    delete [] cws; cws = 0;

    delete [] w0; w0 = 0;

    delete [] w0_cws; w0_cws = 0;

    for (int i = 0; i < ndir; ++i) {
        delete [] w[i];
        delete [] w_cws[i];
        delete [] v[i];
    }

    delete [] w; w = 0;
    delete [] w_cws; w_cws = 0;
    delete [] v; v = 0;

    for (int i = 0; i < niter; ++i) {
        delete [] xfilt[i];
        delete [] idxres[i];
    }
    delete [] xfilt; xfilt = 0;
    delete [] idxres; idxres = 0;
    delete [] prior; prior = 0;
    delete [] lhood; lhood = 0;
    delete [] res_csw; res_csw = 0;

    delete [] neff; neff = 0;
    delete [] xc; xc = 0;
//    delete [] tagg; tagg = 0;
//    delete [] ovlp; ovlp = 0;

    for (int i = 0; i < npcles; ++i) delete [] conv[i];
    delete [] conv; conv = 0;
    delete [] labels; labels = 0;
    delete [] checked; checked = 0;
    delete [] nbridxs; nbridxs = 0;

}

void Tracker::sphereXYZ(float _radius, float _zdist, vector<offxyz>& _offxyz) {
    // export set of offsets for particular radius, and z axis shrink as in off3
    // neighbourhood offsets x,y,z
    _offxyz.clear();

//    for (int i = 0; i <= ceil(4 * _sigEnd); ++i) { // 4*sig would be a very large radius, to cover all possible radiuses
//        vector<offxyz> t;     // empty vector
//        off3.push_back(t);    // now off3[i] can be appended and accessed

        _radius = (_radius<0)?0:_radius;
        int rxy = round(_radius);
        int rz  = round(_radius/_zdist); // scale the layer range when establishing fill-up offsets

//        cout << "rxy=" << rxy <<" rz=" << rz << endl;

//        if (i==0) {
//            offxyz curroff(0, 0, 0);
//            off3[i].push_back(curroff);
//        }
//        else {

        for (int dx = -rxy; dx <= rxy; ++dx) {

            for (int dy = -rxy; dy <= rxy; ++dy) {

                for (int dz = -rz; dz <= rz; ++dz) {

//                    double ins= ((dx*dx)/(double)(rxy*rxy))+((dy*dy)/(double)(rxy*rxy))+((rz==0)?0:(dz*dz)/(double)(rz*rz));
                    if (((dx*dx)/(double)(rxy*rxy))+((dy*dy)/(double)(rxy*rxy))+((rz==0)?0:(dz*dz)/(double)(rz*rz))<=1.0) {
                        offxyz o(dx, dy, dz);
                        _offxyz.push_back(o);
                    }

                }
            }
        }

//        }
//    }

}

void Tracker::export_off3(string path) { // off3[][] debug
    ofstream f;
    f.open(path.c_str(), ofstream::out | ofstream::trunc);

    float shift = 5*sig.back(); // shift to separate the vizualization
    int cnt = 1;
    for (int i = 0; i < off3.size(); ++i) {
        for (int j = 0; j < off3[i].size(); ++j) {
            f<<(cnt++)<<" "<<(i%off3.size())<<" "<<(off3[i][j].x+(i%off3.size())*shift)<<" "<<(off3[i][j].y)<<" "<<(off3[i][j].z)<<" .3 -1\n";
        }
    }

    f.close();

}

void Tracker::export_model(string savepath, bool directed) {

    ofstream f;
    f.open(savepath.c_str(), ofstream::out | ofstream::trunc);

    if (!directed) {
        float shift = 2*3*sig.back(); // separation
        int cnt = 1;
        for (int i = 0; i < model_vuw.size(); ++i) {
            for (int j = 0; j < model_vuw[i].size(); ++j) {
                f<<(cnt++)<<" "<<(i%model_vuw.size())<<" "<<(model_vuw[i][j].v+(i%model_vuw.size())*shift)<<" "<<(model_vuw[i][j].u)<<" "<<(model_vuw[i][j].w)<<" "<<(model_wgt[i][j])<<" -1\n";
            }
        }
    }
    else {

        srand (time(NULL));
        float vx = (float) (rand()%100+1);
        float vy = (float) (rand()%100+1);
        float vz = (!is2d)? (float) (rand()%100+1) : 0;


        float vn = sqrt(vx*vx+vy*vy+vz*vz);

        vx /= vn;
        vy /= vn;
        vz /= vn;

//        cout << vx << "," << vy << "," << vz << endl;

        // align models using the same rotation scheme as in zncc() with ux,uy,uz and wx,wy,wz
        float nrm = sqrt(pow(vx,2)+pow(vy,2)); // projection onto xy plane, if it is small then the vector has z component only
        if (nrm>FLT_MIN) {
            // (ux, uy, uz) denotes the unit direction of the orthogonal positioned in xy plane
            int sg = (vy<0)? -1 : 1;
            ux =  sg * (vy/nrm);
            uy = -sg * (vx/nrm);
            uz =  0;
        }
        else {
            ux = 1;
            uy = 0;
            uz = 0;
        }

        // wx, wy, wz
        if (is2d) {
            wx = 0;
            wy = 0;
            wz = 0;
        }
        else {
            wx =   uy*vz - uz*vy;
            wy = - ux*vz + uz*vx;
            wz =   ux*vy - uy*vx;
        }

        float x,y,z;

        float shift = 2*3*sig.back();
        int cnt = 1;

        for (int i = 0; i < model_vuw.size(); ++i) {
            for (int j = 0; j < model_vuw[i].size(); ++j) {

                x = model_vuw[i][j].v * (-vx) + model_vuw[i][j].u * ux + model_vuw[i][j].w * wx;
                y = model_vuw[i][j].v * (-vy) + model_vuw[i][j].u * uy + model_vuw[i][j].w * wy;
                z = model_vuw[i][j].v * (-vz) + model_vuw[i][j].u * uz + model_vuw[i][j].w * wz;

                x += (i%model_vuw.size())*shift;

                f<<(cnt++)<<" "<<(i%model_vuw.size())<<" "<<x<<" "<<y<<" "<<z<<" "<<model_wgt[i][j]<<" -1\n";

            }

        }

    }

    f.close();

}

void Tracker::export_track(string savepath) {
    ofstream f;
    f.open(savepath.c_str(), ofstream::out | ofstream::trunc);

    int cnt = 1;
    f<<(cnt++)<<" "<<Node::NOTHING<<" "<<x0.x<<" "<<x0.y<<" "<< x0.z <<" "<< x0.sig <<(-1)<<endl;

    for (int i = 0; i < ti_limit; ++i) {
        f<<cnt<<" "<<Node::BASAL_DENDRITE<<" "<<xc[i].x<<" "<<xc[i].y<<" " <<xc[i].z<<" "<<xc[i].sig<<" "<<((i==0)?-1:cnt-1)<<endl;
        cnt++;
    }
    f.close();
}

void Tracker::export_trackcorr(string savepath) {
    ofstream f;
    f.open(savepath.c_str(), ofstream::out | ofstream::trunc);

    for (int i = 0; i < ti_limit; ++i) {
        f<<xc[i].corr<<endl;
        for (int j = 0; j < npcles; ++j) {
            f<<xfilt[i][j].corr<<((j<npcles-1)?",":"")<<flush;
        }
        f<<endl;
    }

    f.close();
}

int Tracker::getdirection(float _vx, float _vy, float _vz) {

    int idx = -1;
    float maxdotp = -FLT_MAX;
    float currdotp;
    for (int i = 0; i < ndir; i++) { // go through all of the predefined directions
        currdotp = _vx * v[i][0] + _vy * v[i][1] + _vz * v[i][2];
        if (currdotp>maxdotp) {
            maxdotp = currdotp;
            idx = i;
        }
    }

    if (idx==-1) cout << "(vx,vy,vz) direction index could not be found." << endl;

    return idx;

}

void Tracker::generate_directions(bool is2D, float** vxyz){

    double h_k, theta_k, phi_k, phi_k_1 = 0;
    int nrdirs = (is2D)?ndirs2d:ndirs3d;

    for (int k = 0; k < nrdirs; k++) { // generate nrdirs directions
        if (is2D) {
            float ang1 = 0.0 + k * ((2*3.14)/(float)nrdirs);
            vxyz[k][0] = cos(ang1);
            vxyz[k][1] = sin(ang1);
            vxyz[k][2] =  0;

//            cout << "|v" << k << "|=" << sqrt(pow(vxyz[k][0],2)+pow(vxyz[k][1],2)+pow(vxyz[k][2],2)) << endl;
        }
        else {
            h_k = 1 - 2 * ((double)k/(nrdirs-1)); // 1 : -1 defines angular range
            theta_k = acos(h_k);

            if(k==0 || k==(nrdirs-1)) {
                phi_k   = 0;
                phi_k_1 = 0;
            }
            else {
                phi_k = phi_k_1 + 3.6 / ( sqrt(nrdirs) * sqrt(1 - h_k * h_k));
                phi_k_1 = phi_k;
            }

            vxyz[k][0] = (float) (sin(theta_k) * cos(phi_k));
            vxyz[k][1] = (float) (sin(theta_k) * sin(phi_k));
            vxyz[k][2] = (float)  cos(theta_k);

//            cout << "|v" << k << "|=" << sqrt(pow(vxyz[k][0],2)+pow(vxyz[k][1],2)+pow(vxyz[k][2],2)) << endl;
        }
    }

}

void Tracker::sampleN(vector<float> _csw, int _N, vector<int>& _sampled_idx){
    srand (time(NULL));
    float u1 = (_csw[_csw.size()-1]/_N) * ((float)rand()/RAND_MAX);
    _sampled_idx.clear();
    int samp = 0;
    for (int i = 0; i < _N; ++i) {
        float ui = u1 + i * (_csw[_csw.size()-1]/_N);
        while (ui > _csw[samp]) samp++;
        _sampled_idx.push_back(samp);
    }
}

void Tracker::trackNeg(seed _seed0, unsigned char* _img, vector<Node>& _nodelist, int _w, int _h, int _l, int* smap, unsigned char* _trc_den,
                       int vol, long** ioff, int* nidx_map) {
    seed seedNeg(_seed0.x, _seed0.y, _seed0.z, -_seed0.vx, -_seed0.vy, -_seed0.vz, _seed0.score, _seed0.corr);
    trackPos(seedNeg, _img, _nodelist, _w, _h, _l, smap, _trc_den, vol, ioff, nidx_map);
}

void Tracker::trackPos(seed _seed0, unsigned char* _img, vector<Node>& _nodelist, int _w, int _h, int _l, int* smap, unsigned char* _trc_den,
                       int vol, long** ioff, int* nidx_map) {

    x0.x  = _seed0.x;
    x0.y  = _seed0.y;
    x0.z  = _seed0.z;
    x0.vx = _seed0.vx;
    x0.vy = _seed0.vy;
    x0.vz = _seed0.vz;
    x0.w = 1;

    bool success;
    bool density_limit_reached;
    bool soma_reached;
    ti_limit = niter;

    if (verbose) printf("\n\t [%.2f, %.2f, %.2f] >> ", x0.vx, x0.vy, x0.vz);

    for (int i = 0; i < niter; ++i) {

        if (i==0)   success = iter0New(   _img, _w, _h, _l);
            else    success = iterINew(i, _img, _w, _h, _l);

        if (success) { // success means round(xc[i].x,xc[i].y,xc[i].z) is within the image dimensions

            long crd = (int)round(xc[i].z)*_w*_h+(int)round(xc[i].y)*_w+(int)round(xc[i].x); // get traces current voxel index

            // check if the soma was reached, if smap[crd]>0, if so, stop the trace and the node and link with the corresponding soma node
            soma_reached = smap[crd]>0;
            // check if voxel density is reched in current vox
            density_limit_reached = (int)_trc_den[crd] >= nodespervol; // nidx_map[] is filled together with _trc_den[]

            if (soma_reached) {

                if (i>0) { // linking with corresponding reached soma node
                    _nodelist[smap[crd]         ].nbr.push_back(_nodelist.size()-1);
                    _nodelist[_nodelist.size()-1].nbr.push_back(smap[crd]);
                }

                ti_limit = i; // current one is the limit

                printf("\n--%d[%d], SOMA, idx=%d", ti_limit, niter, smap[crd]);
                fflush(stdout);
                break; // stop further trace
            }
            else if (density_limit_reached) {

                if (i>0) { // linking with corresponding reached node
                    _nodelist[nidx_map[crd]     ].nbr.push_back(_nodelist.size()-1);
                    _nodelist[_nodelist.size()-1].nbr.push_back(nidx_map[crd]);
                }

                ti_limit = i;

                printf("\n--%d[%d], DENSITY, nodespervol=%d", ti_limit, niter, nodespervol);
                fflush(stdout);
                break; // stop further trace
            }
            else {
                Node nd(xc[i].x, xc[i].y, xc[i].z, xc[i].vx, xc[i].vy, xc[i].vz, xc[i].corr, xc[i].sig, ((i==0)? Node::UNDEFINED : Node::AXON));
                _nodelist.push_back(nd);

                // mark the trace node density of the neighbourhood down
                _trc_den[crd] = (unsigned char)((int)_trc_den[crd] + 1); // update node density map
                nidx_map[crd] = _nodelist.size()-1;

                if (vol>1) { // there are vol-1 neighbours to add to _trc_den[]
                    for (int j = 0; j < vol-1; ++j) {
                        _trc_den[ioff[crd][j]] = (unsigned char)((int)_trc_den[ioff[crd][j]] + 1); // update node density map
                        nidx_map[ioff[crd][j]] = _nodelist.size()-1; // update node index map
                    }
                }

                if (i>0) { // linking
                    _nodelist[_nodelist.size()-1].nbr.push_back(_nodelist.size()-2);
                    _nodelist[_nodelist.size()-2].nbr.push_back(_nodelist.size()-1);
                }

                if (verbose)
                    printf("\n\t\ti=%d\t x=[%4.1f, %4.1f, %4.1f]\t r=[%4.2f]\t zncc=%1.2f\t Neff=%1.2f[%d]", i,
                           xc[i].x, xc[i].y, xc[i].z, xc[i].sig, xc[i].corr, (neff[i]/npcles), (neff[i]/npcles<neff_ratio));

                if (i==niter-1) {// last iteration
                    printf("\n--%d[%d], TRACK LIMIT, niter=%d", ti_limit, niter, xc[ti_limit].corr, niter);
                    fflush(stdout);
                }
            }

        }
        else { // breaks the trace, trace out of dimensions or correlation below limit
            ti_limit = i;
            printf("\n--%d[%d], success=0, corr=%1.2f", ti_limit, niter, xc[ti_limit].corr);
            fflush(stdout);
            break;
        }

    }

    /* this was initially added but later on made no sense to constrain as later processing stages will take care of prunning
    if (ti_limit==1) { // there was 1 node in the trace
        for (int j = ti_limit-1; j >= 0; --j)
            _nodelist.pop_back();
    }
    */

    if (ti_limit>1)
        _nodelist.back().type = Node::END;

}

/*
void Tracker::trackNew(int _x, int _y, int _z, unsigned char * _img, vector<Node> & _nodelist, int _w, int _h, int _l, int * _nmap, int ovlp_window) {

    x0.x = _x;
    x0.y = _y;
    x0.z = _z;
    x0.vx = NAN;
    x0.vy = NAN;
    x0.vz = NAN;
    x0.w = 1;

    bool success;
    ti_limit = niter;
    int cnt_ovlp = 0;

    for (int i = 0; i < niter; ++i) {

        if (i==0)       success = iter0New(   _img, _w, _h, _l); // fills xc[0]
        else            success = iterINew(i, _img, _w, _h, _l); // fills xc[i]

//        ovlp[i] = _nmap[(int)round(xc[i].z)*(_w*_h)+(int)round(xc[i].y)*_w+(int)round(xc[i].x)];
//        if (ovlp[i]>0)  cnt_ovlp++; // count consecutive overlaps
//        else            cnt_ovlp=0; // reset counters

        if (cnt_ovlp>=ovlp_window) {
            ti_limit = i;
            printf("i=%d: cnt_ovlp>=%d\n", i, ovlp_window);
            break;
        }

        if (!success) {
            ti_limit = i;
            printf("i=%d: c=[%1.2f]\n", i, xc[i].corr);
            break;
        }
        else {
            // add node
//            tagg[i] = _nodelist.size();

            Node nd(xc[i].x, xc[i].y, xc[i].z, xc[i].sig, Node::FORK);

            _nodelist.push_back(nd);

            if (i>0) {
                _nodelist[_nodelist.size()-1].nbr.push_back(_nodelist.size()-2);
                _nodelist[_nodelist.size()-2].nbr.push_back(_nodelist.size()-1);
            }

//            if (ovlp[i]>0) {
//                _nodelist[_nodelist.size()-1].nbr.push_back(ovlp[i]);
//                _nodelist[ovlp[i]].nbr.push_back(_nodelist.size()-1);
//            }

            if (verbose)
                printf("i=%d\t x=[%4.1f, %4.1f, %4.1f]\t r=[%4.2f]\t zncc=%1.2f\t Neff=%1.2f[%d]\n", i, xc[i].x, xc[i].y, xc[i].z, xc[i].sig, xc[i].corr, (neff[i]/npcles), (neff[i]/npcles<neff_ratio));
        }
    }

    if (ti_limit>0) {
        for (int i = 0; i < ti_limit; ++i) {
//            fill(xc[i], tagg[i], _nmap, _w, _h, _l);
        }
    }

}
*/
bool Tracker::iter0New(unsigned char * _img, int _w, int _h, int _l) { // , vector<Node> & _nodelist

    srand (time(NULL));

    float wnorm_prior = 0;
    float u1 = (w0_cws[sz-1]/npcles) * ((float) rand() / RAND_MAX);
    int s = 0; // sampling index

    for (int i = 0; i < npcles; ++i) { // npcles random samples from w0 distribution, with _xp as source

        float ui = u1 + i * (w0_cws[sz-1]/npcles);

        while (ui > w0_cws[s] && s < (sz-1)) s++;

        xfilt[0][i].x = x0.x + p[s][0];
        xfilt[0][i].y = x0.y + p[s][1];
        xfilt[0][i].z = x0.z + p[s][2];

        xfilt[0][i].vx = isnan(x0.vx)?u[s][0]:x0.vx; // nan directions let the trace find the linear structure and align with it
        xfilt[0][i].vy = isnan(x0.vy)?u[s][1]:x0.vy;
        xfilt[0][i].vz = isnan(x0.vz)?u[s][2]:x0.vz;

        // prior
        prior[i] = w0[s];
        wnorm_prior += prior[i];

        // correlation
        xfilt[0][i].corr = zncc1(xfilt[0][i], _img, _w, _h, _l, xfilt[0][i].sig); // corr. with template [-1,1]
        lhood[i] = exp(Kc * xfilt[0][i].corr);

//        printf("[%f %f %f] (%f %f %f) %f\n",xfilt[0][i].x, xfilt[0][i].y, xfilt[0][i].z, xfilt[0][i].vx, xfilt[0][i].vy, xfilt[0][i].vz, xfilt[0][i].corr);

    }

    float wnorm_posterior = 0;
    for (int i = 0; i < npcles; ++i) {
       xfilt[0][i].w = (1.0/npcles) * (prior[i]/wnorm_prior) * lhood[i];
       wnorm_posterior += xfilt[0][i].w;
    }

    neff[0] = 0; // degeneracy estimate

    for (int i = 0; i < npcles; ++i) {

        xfilt[0][i].w /= wnorm_posterior; // normalize weights

        neff[0] += pow(xfilt[0][i].w,2);

        res_csw[i] = xfilt[0][i].w + ((i>0)?res_csw[i-1]:0); // cummulative sum of weights

    }

    neff[0] = 1.0/neff[0];

    // estimate (centroid) xc[0]
    xc[0].x = xc[0].y = xc[0].z = xc[0].vx = xc[0].vy = xc[0].vz = xc[0].sig = xc[0].corr = 0;
    for (int i = 0; i < npcles; ++i) {
        xc[0].x  += xfilt[0][i].w * xfilt[0][i].x;  // position
        xc[0].y  += xfilt[0][i].w * xfilt[0][i].y;
        xc[0].z  += xfilt[0][i].w * xfilt[0][i].z;
        xc[0].vx += xfilt[0][i].w * xfilt[0][i].vx; // direction
        xc[0].vy += xfilt[0][i].w * xfilt[0][i].vy;
        xc[0].vz += xfilt[0][i].w * xfilt[0][i].vz;
        xc[0].sig  += xfilt[0][i].w * xfilt[0][i].sig;  // radius
    }

    // xest has unit direction, direction is necessary for corr. (and hence likelihood)
    float vnorm = sqrt(pow(xc[0].vx,2)+pow(xc[0].vy,2)+pow(xc[0].vz,2));
    xc[0].vx /= vnorm;
    xc[0].vy /= vnorm;
    xc[0].vz /= vnorm;
    xc[0].corr = zncc2(xc[0], _img, _w, _h, _l, xc[0].sig);

    // stoppage criteria
    int x1 = round(xc[0].x);
    int y1 = round(xc[0].y);
    int z1 = round(xc[0].z);
    if (x1<0 || x1>=_w || y1<0 || y1>=_h || z1<0 || z1>=_l) return false;
    if (xc[0].corr<znccth) return false; // {printf("LOW CORR. %1.3f\n", xc[0].corr); return false;}

    // resample
    if (neff[0]/npcles < neff_ratio) {
        float u1 = (1.0/npcles) * ((float)rand()/RAND_MAX);
        int s = 0;
        for (int i = 0; i < npcles; ++i) {
            float ui = u1 + i * (1.0/npcles);
            while (ui > res_csw[s]) s++;
            idxres[0][i] = s;
        }
    }

    return true;

}

bool Tracker::iterINew(int _i, unsigned char * _img, int _w, int _h, int _l) {

    srand (time(NULL));

    xc[_i].x=xc[_i].y=xc[_i].z=xc[_i].vx=xc[_i].vy=xc[_i].vz=xc[_i].sig=xc[_i].corr=0;

    float wnorm_prior = 0;

    for (int k = 0; k < npcles; ++k) {
        // go through the particles from _i-1 iteration (previous) and predict from each
        // if resampling occured in previous iteration then predict from resampled indexes

        int k1 = (neff[_i-1]/npcles<neff_ratio)?idxres[_i-1][k]:k; // index depending on the prev. iter.

        int vi = getdirection(
                    xfilt[_i-1][k1].vx,
                    xfilt[_i-1][k1].vy,
                    xfilt[_i-1][k1].vz
                    );
        // find the corresponding direction, index of predefined v [0, ndir)

        // sample 1 random particle
        float u1 = (w_cws[vi][sz-1]) * ((float) rand() / RAND_MAX);
        int s = 0;
        while (u1>w_cws[vi][s] && s<sz-1) s++;

        xfilt[_i][k].x = xfilt[_i-1][k1].x + p[s][0];
        xfilt[_i][k].y = xfilt[_i-1][k1].y + p[s][1];
        xfilt[_i][k].z = xfilt[_i-1][k1].z + p[s][2];

        xfilt[_i][k].vx = u[s][0];
        xfilt[_i][k].vy = u[s][1];
        xfilt[_i][k].vz = u[s][2];

        // prior
        prior[k] = w[vi][s];
        wnorm_prior += prior[k];

        // likelihood
        xfilt[_i][k].corr = zncc1(xfilt[_i][k], _img, _w, _h, _l, xfilt[_i][k].sig);
        lhood[k] = exp(Kc * xfilt[_i][k].corr);

    }

    float wnorm_posterior = 0;

    for (int k = 0; k < npcles; ++k) {
       xfilt[_i][k].w = ((neff[_i-1]/npcles<neff_ratio)?(1.0/npcles):xfilt[_i-1][k].w) * (prior[k]/wnorm_prior) * lhood[k];
       wnorm_posterior += xfilt[_i][k].w;
    }

    neff[_i] = 0; // degeneracy estimate

    for (int k = 0; k < npcles; ++k) {

        xfilt[_i][k].w /= wnorm_posterior; // normalize weight

        neff[_i] += pow(xfilt[_i][k].w,2);

        res_csw[k] = xfilt[_i][k].w + ((k>0)?res_csw[k-1]:0); // cummulative sum of weights (return it)

    }

    neff[_i] = 1.0/neff[_i];

    // estimate (centroid) xc[_i]
    for (int k = 0; k < npcles; ++k) {
        xc[_i].x  += xfilt[_i][k].w * xfilt[_i][k].x;  // position
        xc[_i].y  += xfilt[_i][k].w * xfilt[_i][k].y;
        xc[_i].z  += xfilt[_i][k].w * xfilt[_i][k].z;
        xc[_i].vx += xfilt[_i][k].w * xfilt[_i][k].vx; // direction
        xc[_i].vy += xfilt[_i][k].w * xfilt[_i][k].vy;
        xc[_i].vz += xfilt[_i][k].w * xfilt[_i][k].vz;
        xc[_i].sig  += xfilt[_i][k].w * xfilt[_i][k].sig;  // radius
    }

    // xest has unit direction, direction is necessary for corr.
    float vnorm = sqrt(pow(xc[_i].vx,2)+pow(xc[_i].vy,2)+pow(xc[_i].vz,2));
    xc[_i].vx /= vnorm;
    xc[_i].vy /= vnorm;
    xc[_i].vz /= vnorm;
    xc[_i].corr = zncc2(xc[_i], _img, _w, _h, _l, xc[_i].sig);

    // stoppage criteria
    int x1 = round(xc[_i].x);
    int y1 = round(xc[_i].y);
    int z1 = round(xc[_i].z);
    if (x1<0 || x1>=_w || y1<0 || y1>=_h || z1<0 || z1>=_l) return false;
    if (xc[_i].corr<znccth) return false; // {printf("LOW CORR. %1.3f\n", xc[_i].corr); return false;}

    // resample
    if (neff[_i]/npcles < neff_ratio) {
        float u1 = (1.0/npcles) * ((float)rand()/RAND_MAX);
        int s = 0;
        for (int k = 0; k < npcles; ++k) {
            float ui = u1 + k * (1.0/npcles);
            while (ui > res_csw[s]) s++;
            idxres[_i][k] = s;
        }
    }

    return true;
}

// fill methods are abandoned, they involve off3 offset list
void Tracker::fill(int _nidx, vector<Node> _nodelist, int * _map, int _w, int _h, int _l, bool _ow){

    int x0 = round(_nodelist[_nidx].x);
    int y0 = round(_nodelist[_nidx].y);
    int z0 = round(_nodelist[_nidx].z);
    int r0 = round(_nodelist[_nidx].sig);

    int x1, y1, z1, i1;

    for (int k = 0; k < off3[r0].size(); ++k) {

        x1 = x0 + off3[r0][k].x;
        y1 = y0 + off3[r0][k].y;
        z1 = z0 + off3[r0][k].z;

        if (x1>=0 && x1<_w && y1>=0 && y1<_h && z1>=0 && z1<_l){ // boundaries
            if(_ow){
                _map[z1*(_w*_h)+y1*_w+x1] = _nidx; // enforce overwrite
            }
            else { // write if not tagged
                if(_map[z1*(_w*_h)+y1*_w+x1]==0){
                    _map[z1*(_w*_h)+y1*_w+x1] = _nidx;
                }
            }
        }

    }

}

void Tracker::fill(X_est _xc, int _tag, int * _nmap, int _w, int _h, int _l) {

    int ri = round(_xc.sig);
    int x0 = round(_xc.x);
    int y0 = round(_xc.y);
    int z0 = round(_xc.z);

    int x1,y1,z1,i1;

    for (int k = 0; k < off3[ri].size(); ++k) {

        x1 = x0 + off3[ri][k].x;
        y1 = y0 + off3[ri][k].y;
        z1 = z0 + off3[ri][k].z;

        if (x1>=0 && x1<_w && y1>=0 && y1<_h && z1>=0 && z1<_l) {
            if (_nmap[z1*(_w*_h)+y1*_w+x1]<=0) {
                _nmap[z1*(_w*_h)+y1*_w+x1] = _tag;
            }
        }
    }
}

//void Tracker::nodemapFill(X_est _x_est, int _ntag, int * _nodemap, int _w, int _h, int _l, TagBuffer<int> _allow){
//    // fill, overwrite allowed tags with the new tag + no linking
//    int ri=  round(_x_est.r);
//    int x0 = round(_x_est.x);
//    int y0 = round(_x_est.y);
//    int z0 = round(_x_est.z);

//    int x1,y1,z1,i1;

//    for (int k = 0; k < off3[ri].size(); ++k) {

//        x1 = x0 + off3[ri][k].x;
//        y1 = y0 + off3[ri][k].y;
//        z1 = z0 + off3[ri][k].z;

//        if (x1>=0 && x1<_w && y1>=0 && y1<_h && z1>=0 && z1<_l) {

//            int i1 = z1*(_w*_h)+y1*_w+x1;
//            int ntag_nbr = _nodemap[i1];

//            // linking
//            if (ntag_nbr>0) { // overwrite if >zero and allowed
//                if (_allow.isInBuffer(ntag_nbr)){
//                    _nodemap[i1] = _ntag;
//                }
//            }
//            else { // it was zero
//                _nodemap[i1] = _ntag; // overwrite it if zero
//            }
//        }
//    }
//}

//bool Tracker::iterI(int _i, unsigned char * _img, vector<Node> & _nodelist, int * _nodemap, int * _trackmap, int _w, int _h, int _l) {

//    srand (time(NULL));

//    xc[_i].x = xc[_i].y = xc[_i].z =
//            xc[_i].vx = xc[_i].vy = xc[_i].vz = xc[_i].r = xc[_i].corr = 0; // reset

//    float wnorm_prior = 0;

//    for (int k = 0; k < npcles; ++k) {
//        // go through the particles from _i-1 iteration (previous) and predict from each
//        // if resampling occured in previous iteration then predict from resampled indexes

//        int k1 = (neff[_i-1]/npcles<neff_ratio)?idxres[_i-1][k]:k; // index depending on the prev. iter.

//        int vi = getdirection(
//                    xfilt[_i-1][k1].vx,
//                    xfilt[_i-1][k1].vy,
//                    xfilt[_i-1][k1].vz
//                    );
//        // find the corresponding direction, index of predefined v [0, ndir)

//        // sample 1 random particle
//        float u1 = (w_cws[vi][sz-1]) * ((float) rand() / RAND_MAX);
//        int s = 0;
//        while (u1>w_cws[vi][s] && s<sz-1) s++;

//        xfilt[_i][k].x = xfilt[_i-1][k1].x + p[s][0];
//        xfilt[_i][k].y = xfilt[_i-1][k1].y + p[s][1];
//        xfilt[_i][k].z = xfilt[_i-1][k1].z + p[s][2];

//        xfilt[_i][k].vx = u[s][0];
//        xfilt[_i][k].vy = u[s][1];
//        xfilt[_i][k].vz = u[s][2];

//        // prior
//        prior[k] = w[vi][s];
//        wnorm_prior += prior[k];

//        // likelihood
//        xfilt[_i][k].corr = zncc1(xfilt[_i][k], _img, _w, _h, _l, xfilt[_i][k].r);
//        lhood[k] = exp(Kc * xfilt[_i][k].corr);

//    }

//    float wnorm_posterior = 0;

//    for (int k = 0; k < npcles; ++k) {
//       xfilt[_i][k].w = ((neff[_i-1]/npcles<neff_ratio)?(1.0/npcles):xfilt[_i-1][k].w) * (prior[k]/wnorm_prior) * lhood[k];
//       wnorm_posterior += xfilt[_i][k].w;
//    }

//    neff[_i] = 0; // degeneracy estimate

//    for (int k = 0; k < npcles; ++k) {

//        xfilt[_i][k].w /= wnorm_posterior; // normalize weight

//        neff[_i] += pow(xfilt[_i][k].w,2);

//        res_csw[k] = xfilt[_i][k].w + ((k>0)?res_csw[k-1]:0); // cummulative sum of weights (return it)

//    }

//    neff[_i] = 1.0/neff[_i];

//    // --------------------------------------
//    // estimate (centroid) xc[_i]
//    for (int k = 0; k < npcles; ++k) {
//        xc[_i].x  += xfilt[_i][k].w * xfilt[_i][k].x;  // position
//        xc[_i].y  += xfilt[_i][k].w * xfilt[_i][k].y;
//        xc[_i].z  += xfilt[_i][k].w * xfilt[_i][k].z;
//        xc[_i].vx += xfilt[_i][k].w * xfilt[_i][k].vx; // direction
//        xc[_i].vy += xfilt[_i][k].w * xfilt[_i][k].vy;
//        xc[_i].vz += xfilt[_i][k].w * xfilt[_i][k].vz;
//        xc[_i].r  += xfilt[_i][k].w * xfilt[_i][k].r;  // radius
//    }

//    // xest has unit direction, direction is necessary for corr.
//    float vnorm = sqrt(pow(xc[_i].vx,2)+pow(xc[_i].vy,2)+pow(xc[_i].vz,2));
//    xc[_i].vx /= vnorm;
//    xc[_i].vy /= vnorm;
//    xc[_i].vz /= vnorm;
//    // calculate corr
//    xc[_i].corr = zncc2(xc[_i], _img, _w, _h, _l, xc[_i].r);

//    // --------------------------------------
//    // stoppage criteria
//    // out of the boundaries
//    if (xc[_i].x<0 || xc[_i].x>_w-1) {printf("xout\n"); return false;}
//    if (xc[_i].y<0 || xc[_i].y>_h-1) {printf("yout\n"); return false;}
//    if (xc[_i].z<0 || xc[_i].z>_l-1) {printf("zout\n"); return false;}
//    if (xc[_i].corr<znccth)          {printf("corr<znccth %f\n",xc[_i].corr); return false;}

//    // --------------------------------------
//    // track guidance
//    ovlp[_i] = overlaps(xc[_i], _nodemap, _w, _h, _l); // read overlap at rounded locations

//    if (ovlp[_i]==0) { // no nodemap overlap, is tube

//        if(overlaps(xc[_i], _trackmap, _w, _h, _l)==0) { // no track self-overlap

//            Node nd(xc[_i].x, xc[_i].y, xc[_i].z, xc[_i].r, Node::FORK); // add it to the nodelist
//            _nodelist.push_back(nd);
//            added[_i] = _nodelist.size()-1;

//            int taken_out = allow.append(added[_i]); // first allow.sz will give NULL

//            if (taken_out!=NULL)
//                fill(taken_out, _nodelist, _trackmap, _w, _h, _l, 1); // overwrite, latest tag goes on top

//        }
//        else {
//            added[_i] = 0;
//            printf("self-overlap\n");
//            return false; // self-overlap, BREAK further trace
//        }
//    }
//    else {
//        added[_i] = 0;
//        fillAllow(_nodelist, _trackmap, _w, _h, _l, 1); // add if there is anything in the allowed buf. NOBREAK
//        allow.clear(); //tag window with allowed track tags is emptied
//    }

//    // --------------------------------------
//    // linking
//    if (_i>0) {
//        if (added[_i]>0 && ovlp[_i-1]>0) {
//            _nodelist[added[_i] ].nbr.push_back(ovlp[_i-1]);
//            _nodelist[ovlp[_i-1]].nbr.push_back(added[_i] );
//        }
//        if (ovlp[_i]>0 && added[_i-1]>0) {
//            _nodelist[ovlp[_i]   ].nbr.push_back(added[_i-1]);
//            _nodelist[added[_i-1]].nbr.push_back(ovlp[_i]);
//        }
//        if (added[_i]>0 && added[_i-1]) {
//            _nodelist[added[_i]  ].nbr.push_back(added[_i-1]);
//            _nodelist[added[_i-1]].nbr.push_back(added[_i]  );
//        }
//    }

//    // --------------------------------------
//    // resample
//    if (neff[_i]/npcles < neff_ratio) {

//        if (0) { // clustering X[] at iteration _i
//            ms(_i, KRAD); // conv  KRAD, D todo, can be warped in z direction in runOne()
//            clustering(KRAD); // labels, nbridxs
//            extract(_i); // clst_idx, clst_csw

//            // systematic resampling
//            float u1 = (clst_csw.back()/npcles) * ((float)rand()/RAND_MAX);
//            int s = 0;
//            for (int k = 0; k < npcles; ++k) {
//                float ui = u1 + k * (clst_csw.back()/npcles);
//                while (ui > clst_csw[s]) s++;
//                idxres[_i][k] = s;
//            }
//        }
//        else { // regular
//            float u1 = (1.0/npcles) * ((float)rand()/RAND_MAX);
//            int s = 0;
//            for (int k = 0; k < npcles; ++k) {
//                float ui = u1 + k * (1.0/npcles);
//                while (ui > res_csw[s]) s++;
//                idxres[_i][k] = s;
//            }
//        }

//    }

//}

//bool Tracker::iter0(unsigned char * _img, vector<Node> & _nodelist, int * _nodemap, int * _trackmap, int _w, int _h, int _l) {
//    srand (time(NULL));

//    float wnorm_prior = 0;

//    float u1 = (w0_cws[sz-1]/npcles) * ((float) rand() / RAND_MAX);

//    int s = 0; // sampling index

//    for (int i = 0; i < npcles; ++i) { // npcles random samples from w0 distribution, with _xp as source

//        float ui = u1 + i * (w0_cws[sz-1]/npcles);

//        while (ui > w0_cws[s] && s < (sz-1)) s++;

//        xfilt[0][i].x = x0.x + p[s][0];
//        xfilt[0][i].y = x0.y + p[s][1];
//        xfilt[0][i].z = x0.z + p[s][2];

//        xfilt[0][i].vx = isnan(x0.vx)?u[s][0]:x0.vx;
//        xfilt[0][i].vy = isnan(x0.vy)?u[s][1]:x0.vy;
//        xfilt[0][i].vz = isnan(x0.vz)?u[s][2]:x0.vz;

//        // prior
//        prior[i] = w0[s];
//        wnorm_prior += prior[i];

//        // correlation
//        xfilt[0][i].corr = zncc1(xfilt[0][i], _img, _w, _h, _l, xfilt[0][i].r); // corr. with template [-1,1]
//        lhood[i] = exp(Kc * xfilt[0][i].corr);

//    }

//    float wnorm_posterior = 0;
//    for (int i = 0; i < npcles; ++i) {
//       xfilt[0][i].w = (1.0/npcles) * (prior[i]/wnorm_prior) * lhood[i];
//       wnorm_posterior += xfilt[0][i].w;
//    }

//    neff[0] = 0; // degeneracy estimate

//    for (int i = 0; i < npcles; ++i) {

//        xfilt[0][i].w /= wnorm_posterior; // normalize weights

//        neff[0] += pow(xfilt[0][i].w,2);

//        res_csw[i] = xfilt[0][i].w + ((i>0)?res_csw[i-1]:0); // cummulative sum of weights

//    }

//    neff[0] = 1.0/neff[0];

//    // --------------------------------------
//    // estimate (centroid) xc[0] at first iteration
//    xc[0].x = xc[0].y = xc[0].z = xc[0].vx = xc[0].vy = xc[0].vz = xc[0].r = xc[0].corr = 0;
//    for (int i = 0; i < npcles; ++i) {
//        xc[0].x  += xfilt[0][i].w * xfilt[0][i].x;  // position
//        xc[0].y  += xfilt[0][i].w * xfilt[0][i].y;
//        xc[0].z  += xfilt[0][i].w * xfilt[0][i].z;
//        xc[0].vx += xfilt[0][i].w * xfilt[0][i].vx; // direction
//        xc[0].vy += xfilt[0][i].w * xfilt[0][i].vy;
//        xc[0].vz += xfilt[0][i].w * xfilt[0][i].vz;
//        xc[0].r  += xfilt[0][i].w * xfilt[0][i].r;  // radius
//    }

//    // xest has unit direction, direction is necessary for corr.
//    float vnorm = sqrt(pow(xc[0].vx,2)+pow(xc[0].vy,2)+pow(xc[0].vz,2));
//    xc[0].vx /= vnorm;
//    xc[0].vy /= vnorm;
//    xc[0].vz /= vnorm;
//    // calculate corr
//    xc[0].corr = zncc2(xc[0], _img, _w, _h, _l, xc[0].r);

//    // --------------------------------------
//    // stoppage criteria
//    if (xc[0].x<0 || xc[0].x>_w-1) {printf("xout\n"); return false;}
//    if (xc[0].y<0 || xc[0].y>_h-1) {printf("yout\n"); return false;}
//    if (xc[0].z<0 || xc[0].z>_l-1) {printf("zout\n"); return false;}
//    if (xc[0].corr<znccth)         {printf("corr<znccth %f\n", xc[0].corr); return false;}

//    // --------------------------------------
//    // track guidance
//    ovlp[0] = overlaps(xc[0], _nodemap, _w, _h, _l); // read overlap at rounded locations

//    if (ovlp[0]==0) { // no nodemap overlap, is tube

//        if(overlaps(xc[0], _trackmap, _w, _h, _l)==0) { // no track self-overlap

//            Node nd(xc[0].x, xc[0].y, xc[0].z, xc[0].r, Node::FORK); // add it to the nodelist
//            _nodelist.push_back(nd);
//            added[0] = _nodelist.size()-1;

//            int taken_out = allow.append(added[0]); // first allow.sz will give NULL

//            if (taken_out!=NULL)
//                fill(taken_out, _nodelist, _trackmap, _w, _h, _l, 1); // overwrite, latest tag goes on top

//        }
//        else {
//            added[0] = 0;
//            printf("self-overlap\n");
//            return false; // self-overlap, BREAK further trace
//        }
//    }
//    else {
//        added[0] = 0;
//        fillAllow(_nodelist, _trackmap, _w, _h, _l, 1); // add if there is anything in the allowed buf. NOBREAK
//        allow.clear(); //tag window with allowed track tags is emptied
//    }

//    // --------------------------------------
//    // resample
//    if (neff[0]/npcles < neff_ratio) {

//        // get indexes of the resampled
//        if (0) { // clustering X[] at iteration 0
//            ms(0, KRAD); // conv  KRAD, D todo, can be warped in z direction in runOne()
//            clustering(KRAD); // labels, nbridxs
//            extract(0); // clst_idx, clst_csw

//            // systematic resampling
//            float u1 = (clst_csw.back()/npcles) * ((float)rand()/RAND_MAX);
//            int s = 0;
//            for (int i = 0; i < npcles; ++i) {
//                float ui = u1 + i * (clst_csw.back()/npcles);
//                while (ui > clst_csw[s]) s++;
//                idxres[0][i] = s;
//            }
//        }
//        else { // regular
//            float u1 = (1.0/npcles) * ((float)rand()/RAND_MAX);
//            int s = 0;
//            for (int i = 0; i < npcles; ++i) {
//                float ui = u1 + i * (1.0/npcles);
//                while (ui > res_csw[s]) s++;
//                idxres[0][i] = s;
//            }
//        }

//    }

//}

void Tracker::extract(int _i){
    // ....
    // needs to have int[] labels calculated
    // results in vector<int> clst_idx, and vector<float> clst_csw
    // this way info on total # clusters is lost
    // ....

    for (int i = 0; i < npcles; ++i) checked[i] = 0; // reset
    clst_idx.clear();
    clst_csw.clear();
    int max_sz = -1;

    for (int i = 0; i < npcles; ++i) {
        if(!checked[i]) {

            checked[i] = 1;
            vector<int> clst;
            clst.push_back(i);

            // check the rest
            for (int j = i+1; j < npcles; ++j) {
                if(!checked[j]){
                    if(labels[j]==labels[i]){
                        clst.push_back(j);
                        checked[j] = 1;
                    }
                }
            }

            if ((int)clst.size()>max_sz) { // found max

                max_sz = clst.size();

                // update output
                clst_idx.clear();
                clst_csw.clear();

                for (int k = 0; k < max_sz; ++k) {

                    clst_idx.push_back(clst[k]); // particle indexes of largest cluster so far
                    clst_csw.push_back(xfilt[_i][clst[k]].w + ((k==0)?0:clst_csw[k-1])); // weight csw

                }

            }
        }
    }
}

void Tracker::clustering(float _dist){ // int _xfLen,

    // ....
    // needs to have float[][] conv calculated
    // results in int[] labels, vector<int>* nbridxs
    // ....

    float dist2 = pow(_dist, 2);

    for (int i = 0; i < npcles; ++i) {
        labels[i] = i;
        nbridxs[i].clear();
    }

    for (int i = 0; i < npcles; ++i) {
        for (int j = i+1; j < npcles; ++j) {
            float d2 = pow(conv[i][0]-conv[j][0], 2);
            if (d2<dist2){
                d2 += pow(conv[i][1]-conv[j][1], 2);
                if (d2<dist2){
                    d2 += pow(conv[i][2]-conv[j][2], 2);
                    if (d2<dist2){
                        nbridxs[i].push_back(j);
                        nbridxs[j].push_back(i);
                    }
                }
            }
        }
    }

    for (int i = 0; i < npcles; ++i) {
        for (int nbri = 0; nbri < nbridxs[i].size(); ++nbri) {

            int j = nbridxs[i][nbri];

            // propagate labels
            if (labels[j]!=labels[i]){

                int currLabel = labels[j];
                int newLabel = labels[i];

                labels[j] = newLabel;

                // set all that were currLabel to newLabel
                for (int k = 0; k < npcles; ++k) {
                    if(labels[k]==currLabel){
                        labels[k] = newLabel;
                    }
                }

            }

        }
    }

}

void Tracker::ms(int _i, float _krad){ // mean-shift

    for (int p = 0; p < npcles; ++p) {
        conv[p][0] = xfilt[_i][p].x;
        conv[p][1] = xfilt[_i][p].y;
        conv[p][2] = xfilt[_i][p].z;
    }

//    float* new_v = new float[3];
    float new_v[3];

    for (int i = 0; i < npcles; ++i) {

        int iter = 0;
        double d2;

        do {

            runOne(conv[i], new_v, xfilt[_i], npcles, _krad);

            d2 =
                    pow(new_v[0]-conv[i][0], 2) +
                    pow(new_v[1]-conv[i][1], 2) +
                    pow(new_v[2]-conv[i][2], 2);

            conv[i][0] = new_v[0];
            conv[i][1] = new_v[1];
            conv[i][2] = new_v[2];

            iter++;

        } while (iter < MAXITER && d2 > EPSILON2);

    }

//    delete [] new_v; new_v = 0;

}

void Tracker::runOne(float * _curr_v, float * _new_v, X * _xf, int _xfLen, float _krad) {

    _new_v[0] = 0; _new_v[1] = 0; _new_v[2] = 0;

    float x2, y2, z2;

    float sum = 0;
    for (int i = 0; i < _xfLen; ++i) {

        x2 = pow(_curr_v[0] - _xf[i].x, 2); // todo optimize the condition sequence here
        y2 = pow(_curr_v[1] - _xf[i].y, 2);
        z2 = pow(_curr_v[2] - _xf[i].z, 2);

        if(x2+y2+z2 <= pow(_krad, 2)) {
            sum++;
            _new_v[0] += _xf[i].x;
            _new_v[1] += _xf[i].y;
            _new_v[2] += _xf[i].z;
        }
    }

    if (sum>0) {
        _new_v[0] /= sum; _new_v[1] /= sum; _new_v[2] /= sum;
    }

}

int Tracker::overlaps(X_est _x_est, int * _nodemap, int _w, int _h, int _l) {
    return overlaps(_x_est.x, _x_est.y, _x_est.z, _nodemap, _w, _h, _l);
}

int Tracker::overlaps(float _x, float _y, float _z, int * _nodemap, int _w, int _h, int _l){
    return _nodemap[(int)(round(_z)*(_w*_h)+round(_y)*_w+round(_x))];
}

//// nodelist values allowed to overlap so that the returned tag is 0
//int Tracker::overlaps(X_est _x_est, TagBuffer<int> _allowed, int* _nodemap, int _w, int _h, int _l){
//    return overlaps(_x_est.x, _x_est.y, _x_est.z, _allowed, _nodemap, _w, _h, _l);
//}

//int Tracker::overlaps(float _x, float _y, float _z, TagBuffer<int> _allowed, int * _nodemap, int _w, int _h, int _l) {
//    int out = _nodemap[(int)(round(_z)*(_w*_h) + round(_y)*_w + round(_x))];

//    if (out==0) return 0;

//    if (_allowed.isInBuffer(out)) return 0;

//    return out;
//}

//int overlaps(float _x, float _y, float _z, TagBuffer<int> _allowed, int * _nodemap, int _w, int _h, int _l);

//int overlaps(X_est _x_est, TagBuffer<int> _allowed, int * _nodemap, int _w, int _h, int _l);

void Tracker::crossprod(float a1, float a2, float a3, float b1, float b2, float b3, float * v) {
    // v is cross product of (a1, a2, a3) and (b1, b2, b3)
    v[0] = a2*b3 - b2*a3;
    v[1] = -(a1*b3-b1*a3);
    v[2] = a1*b2-b1*a2;
}

float Tracker::d2(
        float a1x, float a1y, float a1z,
        float a2x, float a2y, float a2z,
        float a0x, float a0y, float a0z
        ) {

    float a21x = a2x-a1x;
    float a21y = a2y-a1y;
    float a21z = a2z-a1z;

    float a01x = a0x-a1x;
    float a01y = a0y-a1y;
    float a01z = a0z-a1z;

    float a02x = a0x-a2x;
    float a02y = a0y-a2y;
    float a02z = a0z-a2z;

    // d2=0 if it is outside
    if (a21x*a01x+a21y*a01y+a21z*a01z<=0)
        return (a01x*a01x+a01y*a01y+a01z*a01z);
    if ((-a21x)*a02x+(-a21y)*a02y+(-a21z)*a02z<=0)
        return (a02x*a02x+a02y*a02y+a02z*a02z);

    // point to line distance 3D
    // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    float a10x = a1x-a0x;
    float a10y = a1y-a0y;
    float a10z = a1z-a0z;

    float a21_x_a10[3];

    crossprod(a21x, a21y, a21z, a10x, a10y, a10z, a21_x_a10);

    return (float) (
            (pow(a21_x_a10[0],2)+pow(a21_x_a10[1],2)+pow(a21_x_a10[2],2)) / ((a21x*a21x)+(a21y*a21y)+(a21z*a21z))
    );

}

//float Tracker::l1804(X a, X b, unsigned char * img, int w, int h, int l) {
//    // sum of the intensities along the direction
//    float wxyz = 0;
//    float wsum = 0;

//    int t;
//    int x1 = round(a.x);
//    int x2 = round(b.x);
//    if (x1>x2){t=x1;x1=x2;x2=t;}
//    int y1 = round(a.y);
//    int y2 = round(b.y);
//    if (y1>y2){t=y1;y1=y2;y2=t;}
//    int z1 = round(a.z);
//    int z2 = round(b.z);
//    if (z1>z2){t=z1;z1=z2;z2=t;}

//    for (int x = x1; x <= x2; ++x) {
//        for (int y = y1; y <= y2; ++y) {
//            for (int z = z1; z <= z2; ++z) {

//                float wgt = exp(-d2(a.x,a.y,a.z,  b.x,b.y,b.z,  x,y,z)/(2*pow(.75,2)));
//                wxyz += wgt * 1;//zncc();//img[z*(w*h)+y*w+x];//interp(x, y, z, img, w, h, l);
//                wsum += wgt;

//            }
//        }
//    }

//    return wxyz/wsum;
//}

float Tracker::zncc1(X _xp, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {
//    return zncc(_xp.x, _xp.y, _xp.z, _xp.vx, _xp.vy, _xp.vz, 0, img, img_w, img_h, img_l, _sig);// 0: use max 0, 1: use average
    return znccBBB(_xp.x, _xp.y, _xp.z, _xp.vx, _xp.vy, _xp.vz, img, img_w, img_h, img_l, _sig);
}

float Tracker::zncc2(X_est _xp, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {
//    return zncc(_xp.x, _xp.y, _xp.z, _xp.vx, _xp.vy, _xp.vz, 0, img, img_w, img_h, img_l, _sig); // 0: use max, 1: use average
    return znccBBB(_xp.x, _xp.y, _xp.z, _xp.vx, _xp.vy, _xp.vz, img, img_w, img_h, img_l, _sig);
}

float Tracker::znccBBB(float _x, float _y, float _z, float _vx, float _vy, float _vz, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {
    //  ux, uy, uz
    float nrm = sqrt(pow(_vx,2)+pow(_vy,2)); // if the xy projection is close to zero, there is z component only
    if (nrm>0.0001) {
        // (ux, uy, uz) denotes the unit direction of the orthogonal positioned in xy plane
        int sg = (_vy<0)? -1 : 1;
        ux =  sg * (_vy/nrm);
        uy = -sg * (_vx/nrm);
        uz =  0;
    }
    else {
        ux = 1;
        uy = 0;
        uz = 0;
    }

    // wx, wy, wz
    if (is2d) {
        wx = 0;
        wy = 0;
        wz = 0;
    }
    else {
        wx =   uy*_vz - uz*_vy;
        wy = - ux*_vz + uz*_vx;
        wz =   ux*_vy - uy*_vx;
    }

    float x,y,z,corra, corrb, corrc, corr_val, ag;
    float out_corr = -FLT_MAX; // ensure that at least one sigma will score max
//    float out_sig = 1; // gaussian cross section sigma of the optimal solution
//    float out_corr_avg = 0; // average (test how it works), if the flag says to output it

    for (int sig_idx = 0; sig_idx < model2_vuw.size(); ++sig_idx) { // loop sigs

        ag = 0;

        // calcualte average of sampled image values
        for (int offset_idx = 0; offset_idx < model2_vuw[sig_idx].size(); ++offset_idx) {

            x = _x + model2_vuw[sig_idx][offset_idx].v * (-_vx) + model2_vuw[sig_idx][offset_idx].u * ux + model2_vuw[sig_idx][offset_idx].w * wx;
            y = _y + model2_vuw[sig_idx][offset_idx].v * (-_vy) + model2_vuw[sig_idx][offset_idx].u * uy + model2_vuw[sig_idx][offset_idx].w * wy;
            z = _z + model2_vuw[sig_idx][offset_idx].v * (-_vz) + model2_vuw[sig_idx][offset_idx].u * uz + model2_vuw[sig_idx][offset_idx].w * wz;

            model2_img[sig_idx][offset_idx] = interp(x,y,z, img, img_w, img_h, img_l);
            ag += model2_img[sig_idx][offset_idx];

        }

        ag /= model2_vuw[sig_idx].size();

        corra = 0;
        corrb = 0;
        corrc = 0;

        // calculate correlation using obtained average
        for (int offset_idx = 0; offset_idx < model2_vuw[sig_idx].size(); ++offset_idx) {

            corra += (model2_img[sig_idx][offset_idx]-ag) * (model2_wgt[sig_idx][offset_idx]-model2_avg[sig_idx]);
            corrb += pow(model2_img[sig_idx][offset_idx]-ag, 2);
            corrc += pow(model2_wgt[sig_idx][offset_idx]-model2_avg[sig_idx], 2);

        }

        corr_val = (corrb*corrc>FLT_MIN)? corra/sqrt(corrb*corrc) : 0;

        if (corr_val>out_corr) {
            out_corr = corr_val;
            _sig = sig[sig_idx];
        }
    }

    return out_corr;
}

float Tracker::znccAAA(float _x, float _y, float _z, float _vx, float _vy, float _vz, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {
    // ux, uy, uz
    float nrm = sqrt(pow(_vx,2)+pow(_vy,2)); // if the xy projection is close to zero, there is z component only
    if (nrm>0.0001) {
        // (ux, uy, uz) denotes the unit direction of the orthogonal positioned in xy plane
        int sg = (_vy<0)? -1 : 1;
        ux =  sg * (_vy/nrm);
        uy = -sg * (_vx/nrm);
        uz =  0;
    }
    else {
        ux = 1;
        uy = 0;
        uz = 0;
    }

    // wx, wy, wz
    if (is2d) {
        wx = 0;
        wy = 0;
        wz = 0;
    }
    else {
        wx =   uy*_vz - uz*_vy;
        wy = - ux*_vz + uz*_vx;
        wz =   ux*_vy - uy*_vx;
    }

    float x,y,z;

    for (int k = 0; k < sig.size(); ++k) {
        pattern_ag[k] = pattern_corra[k] = pattern_corrb[k] = pattern_corrc[k] = 0;
    }

    std::fill(pattern_ag.begin(), pattern_ag.end(), 0);

    for (int offset_idx = 0; offset_idx < pattern_vuw.size(); ++offset_idx) {

        x = _x + pattern_vuw[offset_idx].v * (-_vx) + pattern_vuw[offset_idx].u * ux + pattern_vuw[offset_idx].w * wx;
        y = _y + pattern_vuw[offset_idx].v * (-_vy) + pattern_vuw[offset_idx].u * uy + pattern_vuw[offset_idx].w * wy;
        z = _z + pattern_vuw[offset_idx].v * (-_vz) + pattern_vuw[offset_idx].u * uz + pattern_vuw[offset_idx].w * wz;

        pattern_img[offset_idx] = interp(x,y,z, img, img_w, img_h, img_l);

        for (int k = 0; k < pattern_sid[offset_idx].size(); ++k) {
            pattern_ag[ pattern_sid[offset_idx][k] ] += pattern_img[offset_idx];
        }
    }

    for (int k = 0; k < sig.size(); ++k) {
        pattern_ag[k] /= pattern_cnt[k];
    }

    // calculate correlation using obtained average
    for (int offset_idx = 0; offset_idx < pattern_vuw.size(); ++offset_idx) {

        for (int k = 0; k < pattern_sid[offset_idx].size(); ++k) {
            pattern_corra[ pattern_sid[offset_idx][k] ] +=
                    (   pattern_img[offset_idx]     - pattern_ag[ pattern_sid[offset_idx][k] ]) *
                    (   pattern_wgt[offset_idx][k]  - pattern_avg[ pattern_sid[offset_idx][k] ]);

            pattern_corrb[ pattern_sid[offset_idx][k] ] +=
                    pow(pattern_img[offset_idx]     - pattern_ag[ pattern_sid[offset_idx][k] ], 2);

            pattern_corrc[ pattern_sid[offset_idx][k] ] +=
                    pow(pattern_wgt[offset_idx][k]  - pattern_avg[ pattern_sid[offset_idx][k] ], 2);

        }

    }

    // pick the largest correlation and assign corresponding sigma
    float out_corr = -FLT_MAX; // ensure that at least one sigma will score max
    float corr_val;
    for (int i = 0; i < sig.size(); ++i) {
        corr_val = (pattern_corrb[i]*pattern_corrc[i]>FLT_MIN)? pattern_corra[i]/sqrt(pattern_corrb[i]*pattern_corrc[i]):0;
        if (corr_val>out_corr) {
            out_corr = corr_val;
            _sig = sig[i];
        }
    }
    return out_corr;
}

float Tracker::zncc(float _x, float _y, float _z, float _vx, float _vy, float _vz, bool _return_avg, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {
    //  ux, uy, uz
    float nrm = sqrt(pow(_vx,2)+pow(_vy,2)); // if the xy projection is close to zero, there is z component only
    if (nrm>0.0001) {
        // (ux, uy, uz) denotes the unit direction of the orthogonal positioned in xy plane
        int sg = (_vy<0)? -1 : 1;
        ux =  sg * (_vy/nrm);
        uy = -sg * (_vx/nrm);
        uz =  0;
    }
    else {
        ux = 1;
        uy = 0;
        uz = 0;
    }

    // wx, wy, wz
    if (is2d) {
        wx = 0;
        wy = 0;
        wz = 0;
    }
    else {
        wx =   uy*_vz - uz*_vy;
        wy = - ux*_vz + uz*_vx;
        wz =   ux*_vy - uy*_vx;
    }

    float x,y,z,corra, corrb, corrc, corr_val, ag;

    // offsets have been defined taking into accound dimensionality and scales
    //(not optimal as it will repeat sampling the same values, beneficial to implement some incremental calculation
    // say it starts with the smallest sqare and infers the image values and average using the smaller ones)
    // find correlation with corresponding template(s)
    float out_corr = -FLT_MAX; // ensure that at least one sigma will score max
    float out_sig = 1; // gaussian cross section sigma of the optimal solution
    float out_corr_avg = 0; // average (test how it works), if the flag says to output it

    for (int sig_idx = 0; sig_idx < model_vuw.size(); ++sig_idx) { // loop sigs

        ag = 0;

        // calcualte average of sampled image values
        for (int offset_idx = 0; offset_idx < model_vuw[sig_idx].size(); ++offset_idx) {

            x = _x + model_vuw[sig_idx][offset_idx].v * (-_vx) + model_vuw[sig_idx][offset_idx].u * ux + model_vuw[sig_idx][offset_idx].w * wx;
            y = _y + model_vuw[sig_idx][offset_idx].v * (-_vy) + model_vuw[sig_idx][offset_idx].u * uy + model_vuw[sig_idx][offset_idx].w * wy;
            z = _z + model_vuw[sig_idx][offset_idx].v * (-_vz) + model_vuw[sig_idx][offset_idx].u * uz + model_vuw[sig_idx][offset_idx].w * wz;

            model_img[sig_idx][offset_idx] = interp(x,y,z, img, img_w, img_h, img_l); // ((x1<0 || x2>=img_w || y1<0 || y2>=img_h || z1<0 || z2>=img_l)?0:interp(x,y,z, img, img_w, img_h, img_l));
            ag += model_img[sig_idx][offset_idx];

        }

        ag /= model_vuw[sig_idx].size();

        corra = 0;
        corrb = 0;
        corrc = 0;

        // calculate correlation using obtained average
        for (int offset_idx = 0; offset_idx < model_vuw[sig_idx].size(); ++offset_idx) {

            corra += (model_img[sig_idx][offset_idx]-ag) * (model_wgt[sig_idx][offset_idx]-model_avg[sig_idx]);
            corrb += pow(model_img[sig_idx][offset_idx]-ag, 2);
            corrc += pow(model_wgt[sig_idx][offset_idx]-model_avg[sig_idx], 2);

        }

        corr_val = (corrb*corrc>FLT_MIN)? corra/sqrt(corrb*corrc) : 0;

        if (_return_avg) out_corr_avg += corr_val;

        if (corr_val>out_corr) {
            out_corr = corr_val;
            out_sig = sig[sig_idx];
        }

    }

    if(_return_avg)
        out_corr_avg /= sig.size();

    _sig = out_sig; // side ouput
    return (_return_avg)?out_corr_avg:out_corr;

}

float Tracker::interp(float _x, float _y, float _z, unsigned char * img, int width, int height, int length) {

    float xclamp = clampf(_x, 0, width-1.001); // _x in [0,width-1) so that x1, x2 are within the boundaries
    int x1 = (int) xclamp;
    int x2 = x1 + 1;
    float xfrac = xclamp - x1;

    float yclamp = clampf(_y, 0, height-1.001);
    int y1 = (int) yclamp; // _y in [0,height-1) so that y1, y2 are within the boundaries
    int y2 = y1 + 1;
    float yfrac = yclamp - y1;

//    if (dbg) {printf("\n>>%d -- %d | %d -- %d\n", x1, x2, y1, y2);}

    if (length==1) { // _z is not necessary

//        bool isIn2D = x1>=0 && x2<width && y1>=0 && y2<height; // expell the check, avoid returning 0, just clamp the indexes
//        if(!isIn2D) {
//            printf("interp() out of boundary [%6.2f, %6.2f, %6.2f],[%d--%d],[%d--%d] M=%d, N=%d, P=%d \n", _x, _y, _z, x1, x2, y1, y2, width, height, length);
//            return 0;
//        }

//        float I11_1 = img[y1*width+x1]; // upper left
//        float I12_1 = img[y1*width+x2]; // upper right
//        float I21_1 = img[y2*width+x1]; // bottom left
//        float I22_1 = img[y2*width+x2]; // bottom right

//        if (dbg) {
//            printf("checking %d %d %d\n", width, height, length);
//            for (int i = 0; i < width*height*length; ++i)
//                if (img[i]>0) {printf("found!\n"); break;}
//            printf("took values:\t %d %f %f %f %f \n", (y1*width+x1), I11_1, I12_1, I21_1, I22_1);
//        }

        //                                 I11_1                    I12_1                                    I21_1                   I22_1
        return (1-yfrac) * ((1-xfrac)*img[y1*width+x1] + xfrac*img[y1*width+x2]) + (yfrac) * ((1-xfrac)*img[y2*width+x1] + xfrac*img[y2*width+x2]);

    }
    else {

        float zclamp = clampf(_z, 0, length-1.001);
        int z1 = (int) zclamp;
        int z2 = z1 + 1;
        float zfrac = zclamp - z1;

//        bool isIn3D = y1>=0 && y2<height && x1>=0 && x2<width && z1>=0 && z2<length;
//        if(!isIn3D) {
//            printf("interp() out of boundary [%6.2f, %6.2f, %6.2f],[%d--%d],[%d--%d],[%d--%d] M=%d, N=%d, P=%d \n", _x, _y, _z, x1, x2, y1, y2, z1, z2, width, height, length);
//            return 0;
//        }

        // take neighbourhood 3d
//        float I11_1 = img[z1*width*height+y1*width+x1]; // upper left
//        float I12_1 = img[z1*width*height+y1*width+x2]; // upper right
//        float I21_1 = img[z1*width*height+y2*width+x1]; // bottom left
//        float I22_1 = img[z1*width*height+y2*width+x2]; // bottom right

//        float I11_2 = img[z2*width*height+y1*width+x1]; // upper left
//        float I12_2 = img[z2*width*height+y1*width+x2]; // upper right
//        float I21_2 = img[z2*width*height+y2*width+x1]; // bottom left
//        float I22_2 = img[z2*width*height+y2*width+x2]; // bottom right

        //                                            I11_1                                    I12_1
        //                                            I21_1                                    I22_1
        //                                            I11_2                                    I12_2
        //                                            I21_2                                    I22_2

        return (1-zfrac) *
                (  (1-yfrac) * ((1-xfrac)*img[z1*width*height+y1*width+x1] + xfrac*img[z1*width*height+y1*width+x2]) +
                   (yfrac  ) * ((1-xfrac)*img[z1*width*height+y2*width+x1] + xfrac*img[z1*width*height+y2*width+x2]))
                +
               (zfrac) *
                (  (1-yfrac) * ((1-xfrac)*img[z2*width*height+y1*width+x1] + xfrac*img[z2*width*height+y1*width+x2]) +
                   (yfrac  ) * ((1-xfrac)*img[z2*width*height+y2*width+x1] + xfrac*img[z2*width*height+y2*width+x2]));

    }

}

double Tracker::bessi( int n, double x) {
   int j;
   double bi,bim,bip,tox,ans;


   if (n < 0)
   {
      return DBL_MAX;
   }
   if (n == 0)
      return( bessi0(x) );
   if (n == 1)
      return( bessi1(x) );


   if (x == 0.0)
      return 0.0;
   else {
      tox=2.0/fabs(x);
      bip=ans=0.0;
      bi=1.0;
      for (j=2*(n+(int) sqrt(ACC*n));j>0;j--) {
         bim=bip+j*tox*bi;
         bip=bi;
         bi=bim;
         if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi *= BIGNI;
            bip *= BIGNI;
         }
         if (j == n) ans=bip;
      }
      ans *= bessi0(x)/bi;
      return  x < 0.0 && n%2 == 1 ? -ans : ans;
   }
}

double Tracker::bessi0( double x ){
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2))))))));
   }
   return ans;
}

double Tracker::bessi1( double x){
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans *= (exp(ax)/sqrt(ax));
   }
   return x < 0.0 ? -ans : ans;
}

//void Tracker::add1x1(float _x, float _y, float _z, unsigned char* _trc_den, int _w, int _h, int _l) {

//    int x = round(_x);

//    if (x>=0 && x<_w) {
//        int y = round(_y);
//        if (y>=0 && y<_h) {
//            int z = round(_z);
//            if (z>=0 && z<_l) {
//                _trc_den[z*_w*_h+y*_w+x] = (unsigned char)((int)_trc_den[z*_w*_h+y*_w+x] + 1);
//            }
//        }
//    }
//}

void Tracker::check1x1(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l) {

    int x = round(_x);

    if (x>=0 && x<_w) {
        int y = round(_y);
        if (y>=0 && y<_h) {
            int z = round(_z);
            if (z>=0 && z<_l) {
                _smap[z*_w*_h+y*_w+x] = true;
            }
        }
    }
}

void Tracker::check2x2(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l) {

    int x = (int) _x;
    int y = (int) _y;
    int z = (int) _z;

    bool right = (_x+0.5)>=(x+1);
    bool down  = (_y+0.5)>=(y+1);

    if (right) {
        if (down) { // down right
            for (int k = 0; k < 4; ++k) {
                int x1 = x + Tracker::NH2x2_dr[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH2x2_dr[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH2x2_dr[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
        else { // up right
            for (int k = 0; k < 4; ++k) {
                int x1 = x + Tracker::NH2x2_ur[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH2x2_ur[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH2x2_ur[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
    }
    else {
        if (down) { // down left
            for (int k = 0; k < 4; ++k) {
                int x1 = x + Tracker::NH2x2_dl[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH2x2_dl[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH2x2_dl[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
        else { // up left
            for (int k = 0; k < 4; ++k) {
                int x1 = x + Tracker::NH2x2_ul[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH2x2_ul[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH2x2_ul[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
    }

}

void Tracker::check3x3(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l) {

    int x = round(_x);
    int y = round(_y);
    int z = round(_z);

    int x1, y1, z1, i1;

    for (int k = 0; k < 9; ++k) {
        x1 = x + Tracker::NH3x3[k][0];
        if (x1>=0 && x1<_w) {
            y1 = y + Tracker::NH3x3[k][1];
            if (y1>=0 && y1<_h) {
                z1 = z + Tracker::NH3x3[k][2];
                if (z1>=0 && z1<_l) {
                    _smap[z1*_w*_h+y1*_w+x1] = true;
                }
                else continue;
            }
            else continue;
        }
        else continue;
    }
}

void Tracker::check4x4(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l) {

    int x = (int) _x;
    int y = (int) _y;
    int z = (int) _z;

    bool right = (_x+0.5)>=(x+1);
    bool down  = (_y+0.5)>=(y+1);

    int Noffsets = 16;

    if (right) {
        if (down) { // down right
            for (int k = 0; k < Noffsets; ++k) {
                int x1 = x + Tracker::NH4x4_dr[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH4x4_dr[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH4x4_dr[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
        else { // up right
            for (int k = 0; k < Noffsets; ++k) {
                int x1 = x + Tracker::NH4x4_ur[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH4x4_ur[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH4x4_ur[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
    }
    else {
        if (down) { // down left
            for (int k = 0; k < Noffsets; ++k) {
                int x1 = x + Tracker::NH4x4_dl[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH4x4_dl[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH4x4_dl[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
        else { // up left
            for (int k = 0; k < Noffsets; ++k) {
                int x1 = x + Tracker::NH4x4_ul[k][0];
                if (x1>=0 && x1<_w) {
                    int y1 = y + Tracker::NH4x4_ul[k][1];
                    if (y1>=0 && y1<_h) {
                        int z1 = z + Tracker::NH4x4_ul[k][2];
                        if (z1>=0 && z1<_l) {
                            _smap[z1*_w*_h+y1*_w+x1] = true;
                        }
                        else continue;
                    }
                    else continue;
                }
                else continue;
            }
        }
    }

}

void Tracker::check5x5(float _x, float _y, float _z, bool* _smap, int _w, int _h, int _l) {

    int x = round(_x);
    int y = round(_y);
    int z = round(_z);

    int x1, y1, z1, i1;

    for (int k = 0; k < 25; ++k) {
        x1 = x + Tracker::NH5x5[k][0];
        if (x1>=0 && x1<_w) {
            y1 = y + Tracker::NH5x5[k][1];
            if (y1>=0 && y1<_h) {
                z1 = z + Tracker::NH5x5[k][2];
                if (z1>=0 && z1<_l) {
                    _smap[z1*_w*_h+y1*_w+x1] = true;
                }
                else continue;
            }
            else continue;
        }
        else continue;
    }
}

//void Tracker::track(int _x, int _y, int _z, unsigned char * _img, vector<Node> & _nodelist,
//                    int * _nodemap,
//                    int * _trackmap,
//                    int _w, int _h, int _l){
//    x0.x = _x;      x0.y = _y;      x0.z = _z;
//    x0.vx = NAN;    x0.vy = NAN;    x0.vz = NAN; // NAN: init iteration calculates directions
//    x0.w = 1;
//    bool success;
//    allow.clear();
//    ti_limit = niter;
//    for (int i = 0; i < niter; ++i) {

//        if (i==0)
//            success = iter0(_img, _nodelist, _nodemap, _trackmap, _w, _h, _l);
//        else
//            success = iterI(i, _img, _nodelist, _nodemap, _trackmap, _w, _h, _l);

//        if (!success) {
//            ti_limit = i;
//            if (verbose) printf("i=%d:\t------(x) [%1.2f | %1.2f]\n", i, xc[i].corr, znccth);
//            break;
//        }
//        else {
//            if (verbose) printf("i=%d\t x=[%4.1f, %4.1f, %4.1f]\t r=[%4.2f]\t zncc=%1.2f\t Neff=%1.2f[%d] %s\n", i, xc[i].x, xc[i].y, xc[i].z, xc[i].r, xc[i].corr, (neff[i]/npcles), (neff[i]/npcles<neff_ratio), ((ovlp[i]>0)?"X":"O"));
//        }
//    }

//    // finish with track
//    if (allow.hasItems())
//        fillAllow(_nodelist, _trackmap, _w, _h, _l, 1); // add if there is anything in the allowed buf. NOBREAK
////    allow.clear(); // redundant call

//    // merge _trackmap and _nodemap, if _nodemap>0 don't change, otherwise append to nodemap
//    if (ti_limit>0) {
//        for (long i = 0; i < (_w*_h*_l); ++i) {
//            if (_trackmap[i]>0) {

//                if (_nodemap[i]==0){
//                    _nodemap[i] = _trackmap[i]; // add track to those voxels where the nodemap did not exist
//                }

//                _trackmap[i] = 0;
//            }

//        }
//    }

//}

// procedure:
// 1-add the node
// 2-link (if added in current or at the previous iteration)
// 3-nodemap fill (if added)
// 4-update allow tag list (once it is used in fillup stage)

//void Tracker::assoc0(int * _nodemap, int _w, int _h, int _l, vector<Node> & _nodelist){
//    // xc[0] is calculated in iter0()
//    istube[0] = xc[0].corr>=znccth;
//    ovlp[0] = overlaps(xc[0], _nodemap, _w, _h, _l); // read at rounded locations

//    if (ovlp[0]>0){ // overlap with other node region
//        // same scenario whether it is a tube or not
//        added[0] = 0;     // 1
//        // 2,3
//        allow.clear(); // 4
//    }
//    else { // no overlap
//        if (istube[0]) { // is tube
//            Node node(xc[0].x, xc[0].y, xc[0].z, xc[0].r, Node::BASAL_DENDRITE);
//            _nodelist.push_back(node);
//            added[0] = _nodelist.size()-1; // 1
//            // 2 linking is in 3
//            nodemapFill(xc[0], added[0], _nodemap, _w, _h, _l, _nodelist);// 3 (only for i==0)
//            allow.append(added[0]); // 4
//        }
//        else {// is not tube
//            added[0] = 0; // 1
//            // 2,3
//            allow.clear(); // 4
//        }
//    }
//}

//void Tracker::assocI(int _i, int * _nodemap, int _w, int _h, int _l, vector<Node> & _nodelist){
//    // xc[i] is calculated in iterI()
//    istube[_i] = xc[_i].corr>=znccth;
//    ovlp[_i] = overlaps(xc[_i],allow,_nodemap,_w,_h,_l);

//    if (ovlp[_i]>0){ // overlap with other node region
//        // same scenario whether it is a tube or not
//        added[_i] = 0; // 1
//        if (added[_i-1]>0){ // 2
//            _nodelist[added[_i-1]].nbr.push_back(ovlp[_i]);
//            _nodelist[ovlp[_i]].nbr.push_back(added[_i-1]);
//        }
//        // 3
//        allow.clear();// 4
//    }
//    else { // no overlap
//        if (istube[_i]){ // is tube
//            Node node(xc[_i].x, xc[_i].y, xc[_i].z, xc[_i].r, Node::BASAL_DENDRITE);
//            _nodelist.push_back(node);
//            added[_i] = _nodelist.size()-1; // 1
//            if (added[_i-1]>0) {// 2
//                _nodelist[added[_i-1]].nbr.push_back(added[_i]  );
//                _nodelist[added[_i]  ].nbr.push_back(added[_i-1]);
//            }
//            else if (ovlp[_i-1]>0){
//                _nodelist[ovlp[_i-1]].nbr.push_back(added[_i] );
//                _nodelist[added[_i] ].nbr.push_back(ovlp[_i-1]);
//            }
//            else {
//                // non-overlapping and no tube
//                // no linking, there is no nodelist index to refer to
//                // unless linking during tag fillup (as in i==0)
//            }
//            nodemapFill(xc[_i],added[_i],_nodemap,_w,_h,_l,allow); // 3
//            allow.append(added[_i]); // 4
//        }
//        else { // is not tube
//            added[_i] = 0; // 1
//            // 2,3
//            allow.clear(); // 4
//        }
//    }
//}

//void Tracker::fillAllow(vector<Node> _nodelist, int * _map, int _w, int _h, int _l, bool _ow){
//    // use tag buffer allow accessing the nodes and fill the map
//    for (int i = 0; i < allow.buf.size(); ++i) {
//        fill(allow.buf[i], _nodelist, _map, _w, _h, _l, _ow);
//    }
//}
