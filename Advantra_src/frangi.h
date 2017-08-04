#ifndef FRANGI3D_H
#define FRANGI3D_H
#include <vector>

class Frangi
{
    public:
    std::vector<float> sig; // set of sigmas (standard dev.) of the gaussian funciton of the cross-section
    float zdist;    // gaussian smoothing along z reduced for 3d stacks
//    Frangi vesselness constants
    float alpha;    // 3d
    float beta;     // 3d
    float BetaOne;  // 2d
    float BetaTwo;  // 2d
    float C;        //

    bool blackwhite;// true: black ridges, false: white ridges

    Frangi(float _sig_min, float _sig_stp, float _sig_max, float _zdist, float _alpha, float _beta, float _C, float _beta_one, float _beta_two);

    ~Frangi();

    void  filter3d(unsigned char* I, int w, int h, int l, unsigned char* J, unsigned char* Sc, float* Vx, float* Vy, float* Vz);

    void  hessian3d(unsigned char* I, int w, int h, int l, float sig, float zdist,
                    float* Dzz, float* Dyy, float* Dyz, float* Dxx, float* Dxy, float* Dxz,
                    float* F);

    void  filter2d(unsigned char* I, int w, int h, int l, unsigned char* J, unsigned char* Sc, float* Vx, float* Vy, float* Vz);
// unsigned char* Dxx8, unsigned char* Dxy8, unsigned char* Dyy8, unsigned char* L18, unsigned char* L28
    void hessian2d(unsigned char* I, int w, int h, float sig, float* Dyy, float* Dxy, float* Dxx, float* F);

    void  imgaussian(unsigned char* I, int w, int h, int l, float sig, float zdist, float* F);
    void  imgaussian(unsigned char* I, int w, int h, float sig, float* F);

    float interpz(int x, int y, float z, float* img, int w, int h, int l);

    void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
    static void eigen_decomposition_static(double A[3][3], double V[3][3], double d[3]);
    static void tred2(double V[3][3], double d[3], double e[3]);
    static void tql2(double V[3][3], double d[3], double e[3]);
    static double hypot2(double x, double y);
    static double absd(double val){ if(val>0){ return val;} else { return -val;} }
};

#endif // FRANGI3D_H
