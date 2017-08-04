#ifndef SEED_H
#define SEED_H
#include <vector>

using namespace std;

struct Puv { // offset point in 2d cross-section
    Puv(int u1, int v1) : u(u1), v(v1) {}
    int u, v;
};

struct Puwv { // offset point in 3d cross-section
    Puwv(int u1, int w1, int v1) : u(u1), w(w1), v(v1) {}
    int u, w, v;
};

struct seed {
    seed(float _x, float _y, float _z, float _vx, float _vy, float _vz, float _score) : x(_x), y(_y), z(_z), vx(_vx), vy(_vy), vz(_vz), score(_score) {}
    float x, y, z;
    float vx, vy, vz;
    float score;
};

class SeedExtractor
{
public:
//    static float Luw;     // scale factor for the radius of the cross section ceil(Luw*sig)
//    static float Lv;     // scale factor for the depth  of the cross section ceil(Lv*sig)
    vector<float> sig;

    // Vx,Vy,Vz orthogonals (ux,uy,uz) and (wx,wy,wz)
//    vector< vector<int> >  offset_u;        // offset indexes (ux,uy,uz)
//    vector< vector<uw> >   offset_uw;       // offset indexes (ux,uy,uz) and (wx,wy,wz)
//    vector< vector<xyz> >  offset_xyz;      // offset indexes along x,y,z
    vector< vector<Puv>  > Suv;  // offset indexes that sample u=() and v=() orthogonal in 2d
    vector< vector<Puwv> > Suwv; // offset indexes that sample u=(), w=() and v=() orthogonals in 3d

    SeedExtractor(float _sig_min, float _sig_stp, float _sig_max, float _sig2r);
    ~SeedExtractor();

    void extract3d(
            unsigned char J8th,
            float seed_scr_th,
            unsigned char* J8,
            int w, int h, int l,
            unsigned char* Sc,
            float* Vx, float* Vy, float* Vz,
            vector<seed>& seeds);

    void extract2d(
            unsigned char J8th,
            float seed_scr_th,
            unsigned char* J8,
            int w, int h, int l,
            unsigned char* Sc,
            float* Vx, float* Vy,
            vector<seed>& seeds);

    static void orthogonals(float vx, float vy, float vz, float& ux, float& uy, float& uz, float& wx, float& wy, float& wz);

    static void orthogonals(float vx, float vy, float& ux, float& uy);

    float interp(float x, float y, float z, unsigned char * img, int w, int h, int l);

    void export_Suv(string savepath);
    void export_Suwv(string savepath);
    static void export_seedlist(vector<seed> slist, string savepath, int type=0, float vscale=5);
};

#endif // SEED_H
