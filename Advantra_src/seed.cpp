#include "seed.h"
#include <cfloat>
#include <iostream>
#include <fstream>
#include <cmath>
#ifdef MIN
#undef MIN
#endif
#define MIN(a, b) ((a)<(b)?(a):(b))
using namespace std;

class CompareBySeedScoreVal {
    vector<seed>* _s;
public:
    CompareBySeedScoreVal(vector<seed>* v) : _s(v) {}
public:
    bool operator() (const int& a, const int& b) const { return (*_s)[a].score > (*_s)[b].score; }
};

//float SeedExtractor::Luw = 1.5; // radial neighbourghood size when looking for local maxima
//float SeedExtractor::Lv  = 0.5; // axial neighbourhood size when looking for local maxima

SeedExtractor::SeedExtractor(float _sig_min, float _sig_stp, float _sig_max, float _sig2r){

    for (float s = _sig_min; s <= _sig_max+FLT_MIN; s+=_sig_stp) {
        sig.push_back(s);
    }

    // generate offests that will be used to capture the neighbourhood at each scale
    for (int i = 0; i < sig.size(); ++i) {

        int Ruw = ceil(_sig2r*sig[i]);
        int Rv  = 1; // ceil(Lv*(sig[i])); // /zdist

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

}

SeedExtractor::~SeedExtractor(){}

void SeedExtractor::extract3d(unsigned char J8th, float seed_scr_th, unsigned char* J8, int w, int h, int l, unsigned char* Sc, float* Vx, float* Vy, float* Vz, vector<seed>& seeds) {

    vector<seed> seedsInit;
    seedsInit.clear();

    float ux, uy, uz;
    float wx, wy, wz;
    float vNeighbor;
    float score;//, scoreMin = FLT_MAX, scoreMax = -FLT_MAX;
    bool isMax;
    int x0, y0, z0;
    float xN, yN, zN;
    int s;

    for (long i = 0; i < (w*h*l); ++i) {

        if (J8[i]>J8th) {

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
                    break;
                }

                score += (int) J8[i] - vNeighbor;

            }

            if (isMax && score > FLT_MIN) {
                score /= Suwv[s].size(); // average
//                if (score < scoreMin) scoreMin = score;
//                if (score > scoreMax) scoreMax = score;
                seed sd(x0, y0, z0, Vx[i], Vy[i], Vz[i], score);
                seedsInit.push_back(sd);
            }
        }
    }

//    // min-max normalize the scores
//    for (int i = 0; i < seedsInit.size(); ++i) {
//        seedsInit[i].score = (seedsInit[i].score-scoreMin) / (scoreMax-scoreMin);
//    }

    // arrange list by score (start from the one with highest score)
    vector<int> indices(seedsInit.size());
    for (int i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareBySeedScoreVal(&seedsInit));

    seeds.clear();
    for (int i = 0; i < indices.size(); ++i) {
        if (seedsInit[indices[i]].score >= seed_scr_th)
            seeds.push_back(seedsInit[indices[i]]);
    }

}

void SeedExtractor::extract2d(unsigned char J8th, float seed_scr_th, unsigned char *J8, int w, int h, int l, unsigned char *Sc, float* Vx, float* Vy, vector<seed>& seeds) {

    vector<seed> seedsInit;
    seedsInit.clear();

    float ux, uy;
    float vNeighbor;
    float score;
    bool isMax;
    int x0, y0, z0;
    float xN, yN;
    int s;

//    cout << "DEBUG" << endl;
//    string dbg_path = "/Users/miroslav/Desktop/t.swc";
//    ofstream f;
//    f.open(dbg_path.c_str(), ofstream::out | ofstream::trunc);
//    int dbg_cnt = 1;
//    srand (time(NULL));

    for (long i = 0; i < (w*h*l); ++i) {

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

//                // DEBUG
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
//                if (score < scoreMin) scoreMin = score;
//                if (score > scoreMax) scoreMax = score;
                seed sd(x0, y0, z0, Vx[i], Vy[i], 0, score);
                seedsInit.push_back(sd);
            }
        } // if > J8th
    } // go thorough the image coordinates

//    f.close();

    // arrange list by score (start from the one with highest score)
    vector<int> indices(seedsInit.size());
    for (int i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareBySeedScoreVal(&seedsInit));

    seeds.clear(); // final set of seeds
    for (int i = 0; i < indices.size(); ++i) {
        if (seedsInit[indices[i]].score >= seed_scr_th)
            seeds.push_back(seedsInit[indices[i]]);
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

void SeedExtractor::export_seedlist(vector<seed> slist, string savepath, int type, float vscale) {

    ofstream f;
    f.open(savepath.c_str(), ofstream::out | ofstream::trunc);

    int count = 1;
    for (long i = 0; i < slist.size(); ++i) {

        f<<count<<" "<<type<<" "<<slist[i].x<<" "<<slist[i].y<<" "<<slist[i].z<<" "<<slist[i].score<<" "<<(-1)<<endl;
        count++;
        f<<count<<" "<<type<<" "<<(slist[i].x+vscale*slist[i].vx)<<" "<<(slist[i].y+vscale*slist[i].vy)<<" "<<(slist[i].z+vscale*slist[i].vz)<<" "<<0.1<<" "<<(count-1)<<endl;
        count++;
//        f<<(cnt++)<<" "<<(i%Suwv.size())<<" "<<(Suwv[i][j].u+(i%Suwv.size())*shift)<<" "<<(Suwv[i][j].w)<<" "<<(Suwv[i][j].v)<<" .3 -1\n";
//        NeuronSWC n0;
//        n0.n = n0.nodeinseg_id = count;
//        n0.type = Node::AXON;
//        n0.x = seeds[i].x;
//        n0.y = seeds[i].y;
//        n0.z = seeds[i].z;
//        n0.r = seeds[i].score;
//        n0.parent = -1;
//        Sd.listNeuron.append(n0);

//        NeuronSWC n1;
//        n1.n = n1.nodeinseg_id = count;
//        n1.type = Node::AXON;
//        n1.x = seeds[i].x + 2*seeds[i].score * seeds[i].vx;
//        n1.y = seeds[i].y + 2*seeds[i].score * seeds[i].vy;
//        n1.z = seeds[i].z + 2*seeds[i].score * seeds[i].vz;
//        n1.r = 0.1;
//        n1.parent = count-1;
//        Sd.listNeuron.append(n1);

    }

//    QString vswc = PARA.inimg_file + "_Seeds.swc";
//    writeSWC_file(vswc.toStdString().c_str(), Sd);


    f.close();

//    cout << savepath << endl;
//    cout << type << endl;
//    cout << rad << endl;

}
