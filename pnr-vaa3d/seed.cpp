/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

#include "seed.h"
#include <cfloat>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#ifdef MIN
#undef MIN
#endif
#define MIN(a, b) ((a)<(b)?(a):(b))
using namespace std;

class CompareSeedScore {
    vector<seed>* _s;
public:
    CompareSeedScore(vector<seed>* v) : _s(v) {}
    bool operator() (const int& a, const int& b) const { return (*_s)[a].score > (*_s)[b].score; }
};
#ifdef  _WIN32
double round(double r);
#endif

//float SeedExtractor::Luw = 1.5; // radial neighbourghood size when looking for local maxima
//float SeedExtractor::Lv  = 0.5; // axial neighbourhood size when looking for local maxima

SeedExtractor::SeedExtractor(vector<float> _sigs, float _sig2r, bool _is2d){ // float _sig_min, float _sig_stp, float _sig_max

    sig = _sigs;
    is2d = _is2d;

    /*
     *  Suv, Suwv
     */
    for (int i = 0; i < sig.size(); ++i) { // generate offests that will be used to capture the neighbourhood at each scale

        int Ruw = ceil(_sig2r*sig[i]);
        int Rv  = 1;//ceil(0.5*sig[i]); //1; // ceil(Lv*(sig[i])); // /zdist

        vector<Puv> Suv_at_sig;
        for (int u = -Ruw; u <= Ruw; ++u) {
            for (int v = -Rv; v <= Rv; ++v) {
                if (u*u>0) {Puv t(u, v); Suv_at_sig.push_back(t);}
            }
        }
        Suv.push_back(Suv_at_sig);

        vector<Puwv> Suwv_at_sig;
        for (int u = -Ruw; u <= Ruw; ++u) {
            for (int w = -Ruw; w <= Ruw; ++w) {
                for (int v = -Rv; v <= Rv; ++v) {
                    if (u*u+w*w>0 && u*u+w*w<=Ruw*Ruw) {Puwv t(u, v, w); Suwv_at_sig.push_back(t);}
                }
            }
        }
        Suwv.push_back(Suwv_at_sig);

    } // sigs

    /*
     *  model_wgt, model_img, model_avg, model_vuw
     */
    // calculate the templates
    // clear vector< vector<> > structures used
    model_avg.clear();

    for (int i = 0; i < model_wgt.size(); ++i) model_wgt[i].clear();
    model_wgt.clear();

    for (int i = 0; i < model_img.size(); ++i) model_img[i].clear();
    model_img.clear();

    for (int i = 0; i < model_vuw.size(); ++i) model_vuw[i].clear();
    model_vuw.clear();

    for (int i = 0; i < sig.size(); ++i) {

        model_avg.push_back(0);
        vector<float> wgt;
        model_wgt.push_back(wgt);
        vector<float> img;
        model_img.push_back(img);
        vector<Ovuw> vuw;
        model_vuw.push_back(vuw);

        // indexing differs in 2d and 3d
        if (is2d) {

            int V2 = ceil(1*sig[i]); // steps in the direciton aligned with vx,vy,vz
            int U2 = ceil(3*sig[i]); // orthogonal

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    float value = exp(-(uu*uu)/(2*pow(sig[i],2)));
                    model_wgt[i].push_back(value);
                    model_img[i].push_back(0); // just allocate storage for the image values
                    Ovuw t(vv, uu, 0);
                    model_vuw[i].push_back(t);
                    model_avg[i] += value;
                }
            }
        }
        else {

            int V2 = ceil(1*sig[i]);
            int U2 = ceil(3*sig[i]);
            int W2 = ceil(3*sig[i]);

            for (int vv = -V2; vv <= V2; ++vv) {
                for (int uu = -U2; uu <= U2; ++uu) {
                    for (int ww = -W2; ww <= W2; ++ww) {
                        float value = exp(-((uu*uu)+(ww*ww))/(2*pow(sig[i],2)));
                        model_wgt[i].push_back(value);
                        model_img[i].push_back(0); // just allocate
                        Ovuw t(vv, uu, ww);
                        model_vuw[i].push_back(t);
                        model_avg[i] += value;
                    }
                }
            }
        }

        model_avg[i] /= model_wgt[i].size();

    }


}

SeedExtractor::~SeedExtractor(){}

void SeedExtractor::extract3d(unsigned char J8th,
                              unsigned char* J8,
                              unsigned char* img,
                              int w, int h, int l,
                              unsigned char* Sc,
                              float* Vx, float* Vy, float* Vz,
                              int* smap,
                              float seed_score_min,
                              float seed_corr_min,
                              vector<seed>&  seeds,
                              vector<long>&  seedi) {

    vector<seed> seedsInit;
    seedsInit.clear();

    float ux, uy, uz;
    float wx, wy, wz;
    float vNeighbor;
    float score;
    float corr;
    float dummy_sigma;
    bool isMax;
    int x0, y0, z0;
    float xN, yN, zN;
    int s;

    for (long i = 0; i < w*h*l; ++i) {  // will go through all the voxels and check for local maxima seed point candidates
        if (J8[i]>J8th && smap[i]==0) { // seed considered if the vesselness is above threshold and doesn't belong to soma

            score = 0;
            orthogonals(Vx[i],Vy[i],Vz[i], ux,uy,uz, wx,wy,wz);
            s = (int) Sc[i]; // scale index
            isMax = true;

            for (int j = 0; j < Suwv[s].size(); ++j) {

                x0 = i%w;
                z0 = i/(w*h);
                y0 = i/w-z0*h;

                xN = x0 + Suwv[s][j].u * ux + Suwv[s][j].w * wx + Suwv[s][j].v * Vx[i];
                yN = y0 + Suwv[s][j].u * uy + Suwv[s][j].w * wy + Suwv[s][j].v * Vy[i];
                zN = z0 + Suwv[s][j].u * uz + Suwv[s][j].w * wz + Suwv[s][j].v * Vz[i];

                vNeighbor = interp(xN, yN, zN, J8, w, h, l);

                if (vNeighbor > (int) J8[i]) {
                    isMax = false;
                    break; // out of for loop
                }

                score += (int) J8[i] - vNeighbor;

            }

            if (isMax && score > FLT_MIN) {
                score /= Suwv[s].size(); // average
                if (score>seed_score_min) {
                    // check the correlation at the candidate (expensive calculation goes last)
                    corr = zncc(x0, y0, z0, Vx[i], Vy[i], Vz[i], 0, img, w, h, l, dummy_sigma);
                    if (corr>seed_corr_min) {
                        seed sd(x0, y0, z0, Vx[i], Vy[i], Vz[i], score, corr);
                        seedsInit.push_back(sd);
                    }
                }
            }
        }
    }

    // arrange list indices by score (start from the one with highest score)
    vector<long> indices(seedsInit.size());
    for (long i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareSeedScore(&seedsInit));

    seeds.clear();
    seedi.clear();

    for (long i = 0; i < indices.size(); ++i) {
        seeds.push_back(seedsInit[indices[i]]); // once added can be accessed as seeds[i]
        seedi.push_back((long)((int)round(seeds[i].z) * w * h + (int)round(seeds[i].y) * w + (int)round(seeds[i].x)));
    }

}

void SeedExtractor::extract2d(unsigned char J8th,
                              unsigned char *J8,
                              unsigned char* img,
                              int w, int h, int l,
                              unsigned char *Sc,
                              float* Vx, float* Vy,
                              int* smap,
                              float seed_score_min,
                              float seed_corr_min,
                              vector<seed>& seeds,
                              vector<long>& seedi) {

    vector<seed> seedsInit;
    seedsInit.clear();

    float ux, uy;
    float vNeighbor;
    float score;
    float corr;
    float dummy_sigma;
    bool isMax;
    int x0, y0, z0;
    float xN, yN;
    int s;


//    string dbg_path = "/Users/miroslav/Desktop/t.swc";
//    ofstream f;
//    f.open(dbg_path.c_str(), ofstream::out | ofstream::trunc);
//    int dbg_cnt = 1;
//    srand (time(NULL));

    for (long i = 0; i < w*h*l; ++i) {
        if (J8[i]>J8th) {

            score = 0;

            ux = -Vy[i];
            uy = Vx[i];

            s = (int) Sc[i]; // scale index
            isMax = true;

            x0 = i%w;
            z0 = 0;//i/(w*h);
            y0 = i/w;//-z0*h;

            for (int j = 0; j < Suv[s].size(); ++j) {

                xN = x0 + Suv[s][j].u * ux + Suv[s][j].v * Vx[i];
                yN = y0 + Suv[s][j].u * uy + Suv[s][j].v * Vy[i];
//                zN = 0;

                vNeighbor = interp(xN, yN, 0, J8, w, h, l);

                if (vNeighbor > (int) J8[i]) {
                    isMax = false;
                    break;
                }

                score += (int) J8[i] - vNeighbor;

            }

            if (isMax && score > FLT_MIN) {

//                int num = rand()%100+1;
//                if (num <= 2) {
//                    f<<(dbg_cnt++)<<" "<<5<<" "<<x0<<" "<<y0<<" "<<z0<<" "<<0.1<<" "<<(-1)<<endl;
//                    f<<(dbg_cnt++)<<" "<<5<<" "<<(x0+5*Vx[i])<<" "<<(y0+5*Vy[i])<<" "<<z0<<" "<<0.1<<" "<<(dbg_cnt-2)<<endl;
//                    f<<(dbg_cnt++)<<" "<<5<<" "<<(x0+5*ux)<<" "<<(y0+5*uy)<<" "<<z0<<" "<<0.1<<" "<<(dbg_cnt-3)<<endl;
//                    for (int j = 0; j < Suv[s].size(); ++j) {
//                        xN = x0 + Suv[s][j].u * ux + Suv[s][j].v * Vx[i];
//                        yN = y0 + Suv[s][j].u * uy + Suv[s][j].v * Vy[i];
//                        f<<(dbg_cnt++)<<" "<<6<<" "<<xN<<" "<<yN<<" "<<0<<" "<<0.1<<" "<<(-1)<<endl;
//                    }
//                }

                score /= Suv[s].size(); // average uint8 difference between central and the surrounding ones
                if (score>seed_score_min) {
                    corr = zncc(x0, y0, z0, Vx[i], Vy[i], 0.0, 0, img, w, h, l, dummy_sigma);
                    if (corr>seed_corr_min) {
                        seed sd(x0, y0, z0, Vx[i], Vy[i], 0, score, corr);
                        seedsInit.push_back(sd);
                    }
                }
            }
        } // if > J8th
    } // go thorough the image coordinates

//    f.close();

    // arrange list by score (start from the one with highest score)
    vector<long> indices(seedsInit.size());
    for (long i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareSeedScore(&seedsInit));

    seeds.clear();
    seedi.clear();

    for (int i = 0; i < indices.size(); ++i) {
        seeds.push_back(seedsInit[indices[i]]);
        seedi.push_back((long)((int)round(seeds[i].z) * w * h + (int)round(seeds[i].y) * w + (int)round(seeds[i].x)));
    }

}

void SeedExtractor::orthogonals(float vx, float vy, float& ux, float& uy) {
    ux = -vy;
    uy =  vx;
}

void SeedExtractor::orthogonals(float vx, float vy, float vz, float& ux, float& uy, float& uz, float& wx, float& wy, float& wz) {

    float n = sqrt(vx*vx+vy*vy);

    if (n>0.00001) {
        ux = (+1) * (vy/n);
        uy = (-1) * (vx/n);
        uz = 0;
    }
    else {
        ux = 1;
        uy = 0;
        uz = 0;
    }

    wx = uy*vz - uz*vy;
    wy = ux*vz - uz*vx;
    wz = ux*vy - uy*vx;

}

float SeedExtractor::interp(float x, float y, float z, unsigned char* img, int w, int h, int l) {

    int x1 = (int) x;
    x1 = (x1<0)?0:x1;
    x1 = (x1>(w-2))?(w-2):x1;           // x1 [0,w-2]
    int x2 = x1 + 1;                    // x2 [1,w-1]
    float xfrac = x - x1;
    xfrac = (xfrac<0.0)?0.0:xfrac;
    xfrac = (xfrac>1.0)?1.0:xfrac;      // xfrac [0.0, 1.0]

    int y1 = (int) y;
    y1 = (y1<0)?0:y1;
    y1 = (y1>(h-2))?(h-2):y1;           // y1 [0,h-2]
    int y2 = y1 + 1;                    // y2 [1,h-1]
    float yfrac = y - y1;
    yfrac = (yfrac<0.0)?0.0:yfrac;
    yfrac = (yfrac>1.0)?1.0:yfrac;      // yfrac [0.0, 1.0]

    if (l==1) { // 2d
        float I11 = img[y1*w+x1];
        float I12 = img[y1*w+x2];
        float I21 = img[y2*w+x1];
        float I22 = img[y2*w+x2];
        // bilinear interpolation in xy
        return (1-yfrac) * ((1-xfrac)*I11 + xfrac*I12) + (yfrac) * ((1-xfrac)*I21 + xfrac*I22);
    }
    else { // 3d
        int z1 = (int) z;
        z1 = (z1<0)?0:z1;
        z1 = (z1>(l-2))?(l-2):z1;       // z1 [0,l-2]
        int z2 = z1 + 1;                // z2 [1,l-1]
        float zfrac = z - z1;
        zfrac = (zfrac<0.0)?0.0:zfrac;
        zfrac = (zfrac>1.0)?1.0:zfrac;  // zfrac [0.0, 1.0]

        float I11_1 = img[z1*w*h+y1*w+x1];
        float I12_1 = img[z1*w*h+y1*w+x2];
        float I21_1 = img[z1*w*h+y2*w+x1];
        float I22_1 = img[z1*w*h+y2*w+x2];

        float I11_2 = img[z2*w*h+y1*w+x1];
        float I12_2 = img[z2*w*h+y1*w+x2];
        float I21_2 = img[z2*w*h+y2*w+x1];
        float I22_2 = img[z2*w*h+y2*w+x2];

        return  (1-zfrac) *
                (  (1-yfrac) * ((1-xfrac)*I11_1 + xfrac*I12_1) + (yfrac) * ((1-xfrac)*I21_1 + xfrac*I22_1) )   +
                zfrac    *
                (  (1-yfrac) * ((1-xfrac)*I11_2 + xfrac*I12_2) + (yfrac) * ((1-xfrac)*I21_2 + xfrac*I22_2) );

    }

}

float SeedExtractor::zncc(float _x, float _y, float _z, float _vx, float _vy, float _vz, bool _return_avg, unsigned char * img, int img_w, int img_h, int img_l, float & _sig) {

    // ux, uy, uz
    float nrm = sqrt(pow(_vx,2)+pow(_vy,2)); // projection onto xy plane, if it is small then the vector has z component only
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

void SeedExtractor::export_Suv(string path) {
    ofstream f;
    f.open(path.c_str(), ofstream::out | ofstream::trunc);
    float shift = 5*sig.back(); // shift to separate the vizualization
    int cnt = 1;
    for (int i = 0; i < Suv.size(); ++i) {
        for (int j = 0; j < Suv[i].size(); ++j) {
            f<<(cnt++)<<" "<<(i%Suv.size())<<" "<<(Suv[i][j].u+(i%Suv.size())*shift)<<" "<<(Suv[i][j].v)<<" "<<0<<" .3 -1\n";// (Suv[i][j].z)
        }
    }
    f.close();
}

void SeedExtractor::export_Suwv(string path) {
    ofstream f;
    f.open(path.c_str(), ofstream::out | ofstream::trunc);
    float shift = 5*sig.back(); // shift to separate the vizualization
    int cnt = 1;
    for (int i = 0; i < Suwv.size(); ++i) {
        for (int j = 0; j < Suwv[i].size(); ++j) {
            f<<(cnt++)<<" "<<(i%Suwv.size())<<" "<<(Suwv[i][j].u+(i%Suwv.size())*shift)<<" "<<(Suwv[i][j].w)<<" "<<(Suwv[i][j].v)<<" .3 -1\n";
        }
    }
    f.close();
}

void SeedExtractor::export_seeds(vector<seed> seeds, string swcpath, int type, float arrow) {

    ofstream f;
    f.open(swcpath.c_str(), ofstream::out | ofstream::trunc);

    long countswc = 1;
    long countseed = 0;
    for (long i = 0; i < seeds.size(); ++i) {
        countseed++;

        f<<countswc<<" "<<type<<" "<<seeds[i].x<<" "<<seeds[i].y<<" "<<seeds[i].z<<" "<<   seeds[i].corr   <<" "<<(-1)<<endl; // seeds[i].score
        countswc++;

        if (arrow>0) {
            f<<countswc<<" "<<type<<" "<<(seeds[i].x+arrow*seeds[i].vx)<<" "<<(seeds[i].y+arrow*seeds[i].vy)<<" "<<(seeds[i].z+arrow*seeds[i].vz)<<" "<<0.1<<" "<<(countswc-1)<<endl;
            countswc++;
        }

    }

    f.close();
}

void SeedExtractor::export_seeds_score(vector<seed> seeds, string logpath) {
    ofstream f;
    f.open(logpath.c_str(), ofstream::out | ofstream::trunc);
    for (long i = 0; i < seeds.size(); ++i) f<<seeds[i].score<<" ";
    f.close();
}

void SeedExtractor::export_seeds_corr(vector<seed> seeds, string logpath) {
    ofstream f;
    f.open(logpath.c_str(), ofstream::out | ofstream::trunc);
    for (long i = 0; i < seeds.size(); ++i) f<<seeds[i].corr<<" ";
    f.close();
}

void SeedExtractor::extractSeeds(double tolerance, unsigned char* J8, int w, int h, int l, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz, vector<seed>& seeds) {

    seeds.clear();

    // makeDirectionOffsets
    int shift = 0, mult=1;
    do {shift++; mult*=2;}
    while (mult < w);

    int intEncodeXMask = mult-1;
    int intEncodeYMask = ~intEncodeXMask;
    int intEncodeShift = shift;
    int dirOffset[]={-w, -w+1, +1, +w+1, +w, +w-1,   -1, -w-1};

    // extracting local maxima in 2d, per layer, as in: https://imagej.nih.gov/ij/developer/source/ij/plugin/filter/MaximumFinder.java.html
    unsigned char* types    = new unsigned char[w*h];       // will be a notepad for pixel types
    int* pList              = new int[w*h];                 // here we enter points starting from a maximum

    for (int z = 0; z < l; ++z) { // go through z layers

        for (long i = 0; i < w*h; ++i) types[i] = 0; // reset before calculating nodes for this layer

        float globalMin = FLT_MAX;
        float globalMax = -FLT_MAX;
        for (int y=0; y<h; y++) { //find local minimum/maximum now
            for (int x=0; x<w; x++) {
                float v = (int)J8[z*w*h+y*w+x];
                if (globalMin>v) globalMin = v;
                if (globalMax<v) globalMax = v;
            }
        }

        // getSortedMaxPoints()
        int nMax = 0;                       // counts local maxima
        for (int y=0; y<h; y++) {           // find local maxima now
            for (int x=0, i=x+y*w; x<w; x++, i++) {
                float v = J8[z*w*h+y*w+x];   // ip.getPixelValue(x,y);
                float vTrue = v;
                if (v==globalMin) continue;
                if (true && (x==0 || x==w-1 || y==0 || y==h-1)) continue;
                bool isMax = true;
                // check wheter we have a local maximum.
                bool isInner = (y!=0 && y!=h-1) && (x!=0 && x!=w-1); //not necessary, but faster than isWithin
                for (int d=0; d<8; d++) {                         // compare with the 8 neighbor pixels
                    if (isInner || isWithin(x, y, d, w, h)) {
                        float vNeighbor = J8[z*w*h+(y+DIR_Y_OFFSET[d])*w+(x+DIR_X_OFFSET[d])];
                        float vNeighborTrue = vNeighbor; //isEDM ? trueEdmHeight(x+DIR_X_OFFSET[d], y+DIR_Y_OFFSET[d], ip) : vNeighbor;
                        if (vNeighbor > v && vNeighborTrue > vTrue) {
                            isMax = false;
                            break;
                        }
                    }
                }
                if (isMax) {
                    types[i] = MAXIMUM;
                    nMax++;
                }
            } // for x
        } // for y

        float vFactor = (float)(2e9/(globalMax-globalMin)); //for converting float values into a 32-bit int
        long* maxPoints = new long[nMax];                  //value (int) is in the upper 32 bit, pixel offset in the lower
        for (long maxPoints_idx = 0; maxPoints_idx < nMax; ++maxPoints_idx) {
            maxPoints[maxPoints_idx] = 0;
        }
        int iMax = 0;
        for (int y=0; y<h; y++) { // enter all maxima into an array
            for (int x=0, p=x+y*w; x<w; x++, p++) {
                if (types[p]==MAXIMUM) {
                    float fValue = J8[z*w*h+y*w+x]; // ip.getPixelValue(x,y);
                    int iValue = (int)((fValue-globalMin)*vFactor); //32-bit int, linear function of float value
                    maxPoints[iMax++] = (long)iValue<<32|p;
                }
            }
        }

        sort(maxPoints, maxPoints+nMax);

        float maxSortingError = 0;

        // analyzeAndMarkMaxima()
        bool excludeEdgesNow = true;
        bool excludeOnEdges = true;
        bool displayOrCount = true;

        for (long i = 0; i < w*h; ++i) pList[i] = 0;

        for (int iMax=nMax-1; iMax>=0; iMax--) {    //process all maxima now, starting from the highest

            int offset0 = (int)maxPoints[iMax];     //type cast gets 32 lower bits, where pixel index is encoded
            //int offset0 = maxPoints[iMax].offset;
            if ((types[offset0]&PROCESSED)!=0)      //this maximum has been reached from another one, skip it
                continue;

            //we create a list of connected points and start the list at the current maximum
            int x0 = offset0 % w;
            int y0 = offset0 / w;
            float v0 = J8[z*w*h+y0*w+x0];            // ip.getPixelValue(x0,y0);
            bool sortingError;

            do {                                    // repeat if we have encountered a sortingError
                pList[0] = offset0;
                types[offset0] |= (EQUAL|LISTED);   // mark first point as equal height (to itself) and listed
                int listLen = 1;                    // number of elements in the list
                int listI = 0;                      // index of current element in the list
                bool isEdgeMaximum = (x0==0 || x0==w-1 || y0==0 || y0==h-1);
                sortingError = false;       //if sorting was inaccurate: a higher maximum was not handled so far
                bool maxPossible = true;         //it may be a true maximum
                double xEqual = x0;                 //for creating a single point: determine average over the
                double yEqual = y0;                 //  coordinates of contiguous equal-height points
                int nEqual = 1;                     //counts xEqual/yEqual points that we use for averaging

                do {                                //while neigbor list is not fully processed (to listLen)
                    int offset = pList[listI];
                    int x = offset % w;
                    int y = offset / w;
                    //if(x==18&&y==20)IJ.write("x0,y0="+x0+","+y0+"@18,20;v0="+v0+" sortingError="+sortingError);
                    bool isInner = (y!=0 && y!=h-1) && (x!=0 && x!=w-1); //not necessary, but faster than isWithin

                    for (int d=0; d<8; d++) {       //analyze all neighbors (in 8 directions) at the same level
                        int offset2 = offset+dirOffset[d];

                        if ((isInner || isWithin(x, y, d, w, h)) && (types[offset2]&LISTED)==0) {
//                            if (isEDM && edmPixels[offset2]<=0) continue;   //ignore the background (non-particles)


                            if ((types[offset2]&PROCESSED)!=0) {
                                maxPossible = false; //we have reached a point processed previously, thus it is no maximum now
                                //if(x0<25&&y0<20)IJ.write("x0,y0="+x0+","+y0+":stop at processed neighbor from x,y="+x+","+y+", dir="+d);
                                break;
                            }

                            int x2 = x+DIR_X_OFFSET[d];
                            int y2 = y+DIR_Y_OFFSET[d];
                            float v2 = J8[z*w*h+y2*w+x2]; // ip.getPixelValue(x2, y2);

                            if (v2 > v0 + maxSortingError) {
                                maxPossible = false;    //we have reached a higher point, thus it is no maximum
                                //if(x0<25&&y0<20)IJ.write("x0,y0="+x0+","+y0+":stop at higher neighbor from x,y="+x+","+y+", dir="+d+",value,value2,v2-v="+v0+","+v2+","+(v2-v0));
                                break;
                            } else if (v2 >= v0-(float)tolerance) {
                                if (v2 > v0) {          //maybe this point should have been treated earlier
                                    sortingError = true;
                                    offset0 = offset2;
                                            v0 = v2;
                                            x0 = x2;
                                            y0 = y2;

                                        }
                                        pList[listLen] = offset2;
                                        listLen++;              //we have found a new point within the tolerance
                                        types[offset2] |= LISTED;
                                        if (x2==0 || x2==w-1 || y2==0 || y2==h-1) {
                                            isEdgeMaximum = true;
                                            if (excludeEdgesNow) {
                                                maxPossible = false;
                                                break;          //we have an edge maximum;
                                            }
                                        }
                                        if (v2==v0) {           //prepare finding center of equal points (in case single point needed)
                                            types[offset2] |= EQUAL;
                                            xEqual += x2;
                                            yEqual += y2;
                                            nEqual ++;
                                        }
                                    }
                                } // if isWithin & not LISTED
                            } // for directions d
                            listI++;

                } while (listI < listLen);


                if (sortingError)  {                  //if x0,y0 was not the true maximum but we have reached a higher one
                            for (listI=0; listI<listLen; listI++)
                                types[pList[listI]] = 0;    //reset all points encountered, then retry

                } else {

                    int resetMask = ~(maxPossible?LISTED:(LISTED|EQUAL));

                    xEqual /= nEqual;
                            yEqual /= nEqual;
                            double minDist2 = 1e20;
                            int nearestI = 0;
                            for (listI=0; listI<listLen; listI++) {
                                int offset = pList[listI];
                                int x = offset % w;
                                int y = offset / w;
                                types[offset] &= resetMask;     //reset attributes no longer needed
                                types[offset] |= PROCESSED;     //mark as processed
                                if (maxPossible) {
                                    types[offset] |= MAX_AREA;
                                    if ((types[offset]&EQUAL)!=0) {
                                        double dist2 = (xEqual-x)*(double)(xEqual-x) + (yEqual-y)*(double)(yEqual-y);
                                        if (dist2 < minDist2) {
                                            minDist2 = dist2;   //this could be the best "single maximum" point
                                            nearestI = listI;
                                        }
                                    }
                                }
                            } // for listI
                            if (maxPossible) {
                                int offset = pList[nearestI];
                                types[offset] |= MAX_POINT;
                                if (displayOrCount && !(excludeOnEdges && isEdgeMaximum)) {
                                    int x = offset % w;
                                    int y = offset / w;
                                    if (true) { // roi==null || roi.contains(x, y)
                                        long si = z*w*h+y*w+x;

                                        float Ux = (((float)Vx[si]/255)*2)-1;
                                        float Uy = (((float)Vy[si]/255)*2)-1;
                                        float Uz = (((float)Vz[si]/255)*2)-1;
                                        float Un = sqrt(pow(Ux, 2) + pow(Uy, 2) + pow(Uz, 2));
                                        seed t(x, y, z, Ux/Un, Uy/Un, Uz/Un, 0, 0);
                                        seeds.push_back(t);
                                    }
                                }
                            }


                } //if !sortingError

            } while (sortingError);             //redo if we have encountered a higher maximum: handle it now.

        } // for all maxima iMax

        delete [] maxPoints; maxPoints = 0;

    } // go through z layers

    delete [] types; types = 0;
    delete [] pList; pList = 0;

}

void SeedExtractor::findMaxima(unsigned char* I, int w, int h, int l, double tolerance, vector<Pxyz>& seed_locs_xyz) {

    seed_locs_xyz.clear();

    // makeDirectionOffsets
    int shift = 0, mult=1;
    do {
        shift++; mult*=2;
    }
    while (mult < w);

    int intEncodeXMask = mult-1;
    int intEncodeYMask = ~intEncodeXMask;
    int intEncodeShift = shift;
    int dirOffset[]={-w, -w+1, +1, +w+1, +w, +w-1,   -1, -w-1};

    // allocate maps
    unsigned char* types = new unsigned char[w*h];  // will be a notepad for pixel types
    int* pList = new int[w*h];                      // here we enter points starting from a maximum

    for (int z = 0; z < l; ++z) { // go through the layers, local maxima are calculated per layer

        for (long i = 0; i < w*h; ++i) types[i] = 0;

        float globalMin = FLT_MAX;
        float globalMax = -FLT_MAX;
        for (int y=0; y<h; y++) { //find local minimum/maximum now
            for (int x=0; x<w; x++) {
                float v = (int)I[z*w*h+y*w+x];
                if (globalMin>v) globalMin = v;
                if (globalMax<v) globalMax = v;
            }
        }

        // getSortedMaxPoints()
        int nMax = 0;  //counts local maxima
        for (int y=0; y<h; y++) {         // find local maxima now
            for (int x=0, i=x+y*w; x<w; x++, i++) {
                float v = I[z*w*h+y*w+x]; // ip.getPixelValue(x,y);
                float vTrue = v;
                if (v==globalMin) continue;
                if (true && (x==0 || x==w-1 || y==0 || y==h-1)) continue;
                bool isMax = true;
                // check wheter we have a local maximum.
                bool isInner = (y!=0 && y!=h-1) && (x!=0 && x!=w-1); //not necessary, but faster than isWithin
                for (int d=0; d<8; d++) {                         // compare with the 8 neighbor pixels
                    if (isInner || isWithin(x, y, d, w, h)) {
                        float vNeighbor = I[z*w*h+(y+DIR_Y_OFFSET[d])*w+(x+DIR_X_OFFSET[d])];
                        float vNeighborTrue = vNeighbor; //isEDM ? trueEdmHeight(x+DIR_X_OFFSET[d], y+DIR_Y_OFFSET[d], ip) : vNeighbor;
                        if (vNeighbor > v && vNeighborTrue > vTrue) {
                            isMax = false;
                            break;
                        }
                    }
                }
                if (isMax) {
                    types[i] = MAXIMUM;
                    nMax++;
                }
            } // for x
        } // for y

        float vFactor = (float)(2e9/(globalMax-globalMin)); //for converting float values into a 32-bit int
        long* maxPoints = new long[nMax];                  //value (int) is in the upper 32 bit, pixel offset in the lower
        for (long maxPoints_idx = 0; maxPoints_idx < nMax; ++maxPoints_idx) {
            maxPoints[maxPoints_idx] = 0;
        }
        int iMax = 0;
        for (int y=0; y<h; y++) { // enter all maxima into an array
            for (int x=0, p=x+y*w; x<w; x++, p++) {
                if (types[p]==MAXIMUM) {
                    float fValue = I[z*w*h+y*w+x]; // ip.getPixelValue(x,y);
                    int iValue = (int)((fValue-globalMin)*vFactor); //32-bit int, linear function of float value
                    maxPoints[iMax++] = (long)iValue<<32|p;
                }
            }
        }

        sort(maxPoints, maxPoints+nMax);

        float maxSortingError = 0;

        // analyzeAndMarkMaxima()
        bool excludeEdgesNow = true;
        bool excludeOnEdges = true;
        bool displayOrCount = true;

//        int* pList = new int[w*h];       //here we enter points starting from a maximum
        for (long i = 0; i < w*h; ++i) pList[i] = 0;

        for (int iMax=nMax-1; iMax>=0; iMax--) {    //process all maxima now, starting from the highest

            int offset0 = (int)maxPoints[iMax];     //type cast gets 32 lower bits, where pixel index is encoded
            //int offset0 = maxPoints[iMax].offset;
            if ((types[offset0]&PROCESSED)!=0)      //this maximum has been reached from another one, skip it
                continue;

            //we create a list of connected points and start the list at the current maximum
            int x0 = offset0 % w;
            int y0 = offset0 / w;
            float v0 = I[z*w*h+y0*w+x0];//isEDM?trueEdmHeight(x0,y0,ip):ip.getPixelValue(x0,y0);
            bool sortingError;

            do {                                    //repeat if we have encountered a sortingError
                pList[0] = offset0;
                types[offset0] |= (EQUAL|LISTED);   //mark first point as equal height (to itself) and listed
                int listLen = 1;                    //number of elements in the list
                int listI = 0;                      //index of current element in the list
                bool isEdgeMaximum = (x0==0 || x0==w-1 || y0==0 || y0==h-1);
                sortingError = false;       //if sorting was inaccurate: a higher maximum was not handled so far
                bool maxPossible = true;         //it may be a true maximum
                double xEqual = x0;                 //for creating a single point: determine average over the
                double yEqual = y0;                 //  coordinates of contiguous equal-height points
                int nEqual = 1;                     //counts xEqual/yEqual points that we use for averaging

                do {                                //while neigbor list is not fully processed (to listLen)
                    int offset = pList[listI];
                    int x = offset % w;
                    int y = offset / w;
                    //if(x==18&&y==20)IJ.write("x0,y0="+x0+","+y0+"@18,20;v0="+v0+" sortingError="+sortingError);
                    bool isInner = (y!=0 && y!=h-1) && (x!=0 && x!=w-1); //not necessary, but faster than isWithin

                    for (int d=0; d<8; d++) {       //analyze all neighbors (in 8 directions) at the same level
                        int offset2 = offset+dirOffset[d];

                        if ((isInner || isWithin(x, y, d, w, h)) && (types[offset2]&LISTED)==0) {
//                            if (isEDM && edmPixels[offset2]<=0) continue;   //ignore the background (non-particles)


                            if ((types[offset2]&PROCESSED)!=0) {
                                maxPossible = false; //we have reached a point processed previously, thus it is no maximum now
                                //if(x0<25&&y0<20)IJ.write("x0,y0="+x0+","+y0+":stop at processed neighbor from x,y="+x+","+y+", dir="+d);
                                break;
                            }

                            int x2 = x+DIR_X_OFFSET[d];
                            int y2 = y+DIR_Y_OFFSET[d];
                            float v2 = I[z*w*h+y2*w+x2]; // isEDM ? trueEdmHeight(x2, y2, ip) : ip.getPixelValue(x2, y2);

                            if (v2 > v0 + maxSortingError) {
                                maxPossible = false;    //we have reached a higher point, thus it is no maximum
                                //if(x0<25&&y0<20)IJ.write("x0,y0="+x0+","+y0+":stop at higher neighbor from x,y="+x+","+y+", dir="+d+",value,value2,v2-v="+v0+","+v2+","+(v2-v0));
                                break;
                            } else if (v2 >= v0-(float)tolerance) {
                                if (v2 > v0) {          //maybe this point should have been treated earlier
                                    sortingError = true;
                                    offset0 = offset2;
                                            v0 = v2;
                                            x0 = x2;
                                            y0 = y2;

                                        }
                                        pList[listLen] = offset2;
                                        listLen++;              //we have found a new point within the tolerance
                                        types[offset2] |= LISTED;
                                        if (x2==0 || x2==w-1 || y2==0 || y2==h-1) {
                                            isEdgeMaximum = true;
                                            if (excludeEdgesNow) {
                                                maxPossible = false;
                                                break;          //we have an edge maximum;
                                            }
                                        }
                                        if (v2==v0) {           //prepare finding center of equal points (in case single point needed)
                                            types[offset2] |= EQUAL;
                                            xEqual += x2;
                                            yEqual += y2;
                                            nEqual ++;
                                        }
                                    }
                                } // if isWithin & not LISTED
                            } // for directions d
                            listI++;

                } while (listI < listLen);


                if (sortingError)  {                  //if x0,y0 was not the true maximum but we have reached a higher one
                            for (listI=0; listI<listLen; listI++)
                                types[pList[listI]] = 0;    //reset all points encountered, then retry

                } else {

                    int resetMask = ~(maxPossible?LISTED:(LISTED|EQUAL));

                    xEqual /= nEqual;
                            yEqual /= nEqual;
                            double minDist2 = 1e20;
                            int nearestI = 0;
                            for (listI=0; listI<listLen; listI++) {
                                int offset = pList[listI];
                                int x = offset % w;
                                int y = offset / w;
                                types[offset] &= resetMask;     //reset attributes no longer needed
                                types[offset] |= PROCESSED;     //mark as processed
                                if (maxPossible) {
                                    types[offset] |= MAX_AREA;
                                    if ((types[offset]&EQUAL)!=0) {
                                        double dist2 = (xEqual-x)*(double)(xEqual-x) + (yEqual-y)*(double)(yEqual-y);
                                        if (dist2 < minDist2) {
                                            minDist2 = dist2;   //this could be the best "single maximum" point
                                            nearestI = listI;
                                        }
                                    }
                                }
                            } // for listI
                            if (maxPossible) {
                                int offset = pList[nearestI];
                                types[offset] |= MAX_POINT;
                                if (displayOrCount && !(excludeOnEdges && isEdgeMaximum)) {
                                    int x = offset % w;
                                    int y = offset / w;
                                    if (true) { // roi==null || roi.contains(x, y)
                                        Pxyz t(x, y, z);
                                        seed_locs_xyz.push_back(t);    //addElement(new int[] {x, y});
                                    }
                                }
                            }


                } //if !sortingError

            } while (sortingError);             //redo if we have encountered a higher maximum: handle it now.

        } // for all maxima iMax

        delete [] maxPoints; maxPoints = 0;

    } // go through z layers

    delete [] types; types = 0;
    delete [] pList; pList = 0;

}

bool SeedExtractor::isWithin(int x, int y, int direction, int width, int height) {
    int xmax = width - 1;
    int ymax = height -1;
    switch(direction) {
        case 0:
            return (y>0);
        case 1:
            return (x<xmax && y>0);
        case 2:
            return (x<xmax);
        case 3:
            return (x<xmax && y<ymax);
        case 4:
            return (y<ymax);
        case 5:
            return (x>0 && y<ymax);
        case 6:
            return (x>0);
        case 7:
            return (x>0 && y>0);
    }
    return false;   //to make the compiler happy :-)
}

int SeedExtractor::DIR_X_OFFSET[8] = {  0,  1,  1,  1,  0, -1, -1, -1 };
int SeedExtractor::DIR_Y_OFFSET[8] = { -1, -1,  0,  1,  1,  1,  0, -1 };
unsigned char SeedExtractor::MAXIMUM    = (unsigned char)1;
unsigned char SeedExtractor::LISTED     = (unsigned char)2;
unsigned char SeedExtractor::PROCESSED  = (unsigned char)4;
unsigned char SeedExtractor::MAX_AREA   = (unsigned char)8;
unsigned char SeedExtractor::EQUAL      = (unsigned char)16;
unsigned char SeedExtractor::MAX_POINT  = (unsigned char)32;
unsigned char SeedExtractor::ELIMINATED = (unsigned char)64;
