#ifndef SEED_H
#define SEED_H
#include <vector>
#include <string>

using namespace std;

struct Pxy {
    Pxy(int x1, int y1) : x(x1), y(y1) {}
    int x, y;
};

struct Pxyz {
    Pxyz(int x1, int y1, int z1) : x(x1), y(y1), z(z1) {}
    int x, y, z;
};

struct Puv { // offset point in 2d cross-section
    Puv(int u1, int v1) : u(u1), v(v1) {}
    int u, v;
};

struct Puwv { // offset point in 3d cross-section
    Puwv(int u1, int w1, int v1) : u(u1), w(w1), v(v1) {}
    int u, w, v;
};

struct Ovuw { // offset v,u,w, template offsets used for correlation
    Ovuw(int _v, int _u, int _w) : v(_v),u(_u),w(_w) {}
    int v, u, w;
};

struct seed {
    seed(float _x, float _y, float _z, float _vx, float _vy, float _vz, float _score, float _corr) : x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz), score(_score), corr(_corr) {}
    float x, y, z;
    float vx, vy, vz;
    float score;
    float corr;
};

class SeedExtractor
{
public:
//    static float Luw;     // scale factor for the radius of the cross section ceil(Luw*sig)
//    static float Lv;     // scale factor for the depth  of the cross section ceil(Lv*sig)
    vector<float> sig;
    bool is2d;

    // Vx,Vy,Vz orthogonals (ux,uy,uz) and (wx,wy,wz)
    vector< vector<Puv>  > Suv;  // offset indexes that sample u=() and v=() orthogonal in 2d
    vector< vector<Puwv> > Suwv; // offset indexes that sample u=(), w=() and v=() orthogonals in 3d

    // variables used to calculate the correlation, same as in Tracker class, SeedExtractor has it's own correlation method, same as the one in Tracker
    vector< vector<float> > model_wgt; // template weights
    vector< vector<float> > model_img; // image values
    vector<float>           model_avg; // template averages
    vector< vector<Ovuw> >  model_vuw; // offsets along orthogonals used for correlation
    float ux, uy, uz;
    float wx, wy, wz;

    SeedExtractor(vector<float> _sigs, float _sig2r, bool _is2d); // float _sig_min, float _sig_stp, float _sig_max
    ~SeedExtractor();

    void extract3d( unsigned char J8th,
                    unsigned char* J8,
                    unsigned char* img,
                    int w, int h, int l,
                    unsigned char* Sc,
                    float* Vx, float* Vy, float* Vz,
                    int* smap,
                    float seed_score_min,
                    float seed_corr_min,
                    vector<seed>&  seeds,
                    vector<long>&  seedi);

    void extract2d( unsigned char J8th,
                    unsigned char* J8,
                    unsigned char* img,
                    int w, int h, int l,
                    unsigned char* Sc,
                    float* Vx, float* Vy,
                    int* smap,
                    float seed_score_min,
                    float seed_corr_min,
                    vector<seed>&  seeds,
                    vector<long>&  seedi);

    static void orthogonals(float vx, float vy, float vz, float& ux, float& uy, float& uz, float& wx, float& wy, float& wz);

    static void orthogonals(float vx, float vy, float& ux, float& uy);

    float interp(float x, float y, float z, unsigned char * img, int w, int h, int l);

    float zncc(float _x, float _y, float _z, float _vx, float _vy, float _vz, bool _return_avg, unsigned char * img, int img_w, int img_h, int img_l, float & _sig);

    void export_Suv(string savepath);
    void export_Suwv(string savepath);
    static void export_seeds(vector<seed> seeds, string swcpath, int type=0, float arrow=-1); // swc with locations
    static void export_seeds_score(vector<seed> seeds, string logpath); // log seed score values
    static void export_seeds_corr(vector<seed> seeds, string logpath); // log seed correlation values

    static void extractSeeds(double tolerance, unsigned char* J8, int w, int h, int l, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz, vector<seed>& seeds);

    static void findMaxima(unsigned char* I, int w, int h, int l, double tolerance, vector<Pxyz>& xyVector);
    static bool isWithin(int x, int y, int direction, int width, int height);
    static int DIR_X_OFFSET[8];
    static int DIR_Y_OFFSET[8];
    static unsigned char MAXIMUM;
    static unsigned char LISTED;
    static unsigned char PROCESSED;
    static unsigned char MAX_AREA;
    static unsigned char EQUAL;
    static unsigned char MAX_POINT;
    static unsigned char ELIMINATED;

};

#endif // SEED_H
