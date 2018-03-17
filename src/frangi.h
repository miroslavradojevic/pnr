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

    // discrete set of directions (can store up to 256, as value used to index the direction is unsigned char)
    std::vector< std::vector<float> > Vxyz; // will be assigned upon each frangi3d/frangi2d() call, array of direction indexes is returned together with vesselness score
    static unsigned char ndirs2d;
    static unsigned char ndirs3d;

    bool blackwhite;// true: black ridges, false: white ridges

    Frangi(std::vector<float> _sigs, float _zdist, float _alpha, float _beta, float _C, float _beta_one, float _beta_two);

    ~Frangi();

    void generate_3d_unit_directions(unsigned char Ndir, std::vector< std::vector<float> >& Vxyz);
    void generate_2d_unit_directions(unsigned char Ndir, std::vector< std::vector<float> >& Vxyz);
    unsigned char get_direction_idx(float vx, float vy, float vz, std::vector< std::vector<float> > Vxyz);
    unsigned char get_direction_idx(float vx, float vy,           std::vector< std::vector<float> > Vxyz);

    void  frangi3d(unsigned char* I, int w, int h, int l, float* J, float& Jmin, float& Jmax, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz);

    void  hessian3d(unsigned char* I, int w, int h, int l, float sig, float zdist,
                    float* Dzz, float* Dyy, float* Dyz, float* Dxx, float* Dxy, float* Dxz);

    void  frangi2d(unsigned char* I, int w, int h, int l, float* J, float& Jmin, float& Jmax, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz);

    void hessian2d(unsigned char* I, int w, int h, float sig, float* Dyy, float* Dxy, float* Dxx);

    static void  imgaussian(unsigned char* I, int w, int h, int l, float sig, float zdist, float* F); // 3d gaussian Gxyz
    static void  imgaussian(unsigned char* I, int w, int h, int l, float sig                       ); // 3d gaussian Gxy
    static void  imgaussian(unsigned char* I, int w, int h,        float sig,              float* F); // 2d gaussian Gxy

    static void imerode(unsigned char* I, int w, int h, int l, float rad, float zdist, unsigned char* E); // 3d
    static void imerode(unsigned char* I, int w, int h, int l, float rad, unsigned char* E); // 3d (xy plane)

    static void imdilate(unsigned char* I, int w, int h, int l, float rad); // 3d (xy plane)

    float interpz(int x, int y, float z, float* img, int w, int h, int l);

    void eigen_decomposition(double A[3][3], double V[3][3], double d[3]);
    static void eigen_decomposition_static(double A[3][3], double V[3][3], double d[3]);
    static void tred2(double V[3][3], double d[3], double e[3]);
    static void tql2(double V[3][3], double d[3], double e[3]);
    static double hypot2(double x, double y);
    static double absd(double val){ if(val>0){ return val;} else { return -val;} }
};

#endif // FRANGI3D_H
