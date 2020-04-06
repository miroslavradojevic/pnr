/* Advantra_plugin.cpp
* Automatic neuron reconstruction from microscopy image stacks.
* 2015-8-19 : initial version
* 2016-4-4  : improved tracker implementation added trace merging mechanism
* 2016-5-24 : added frangi based prefiltering
* 2016-6-16 : added refinement module and export tools

Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

example terminal call:
~/vaa3d/vaa3d_tools/bin/vaa3d64.app/Contents/MacOS/vaa3d64 -x Advantra -f advantra_func -i $inimg_file -p 2,3 12 10 0.5 3 2 200 20 2 4 1
*/
 
#include "v3d_message.h"
#include "basic_surf_objs.h"
#include "nf_dialog.h"
#include "tracker.h"
#include "frangi.h"
#include "seed.h"
#include "node.h"
#include "toolbox.h"
#include <vector>
#include <ctime>
#include <cfloat>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <queue>
#include <climits>
#include <iomanip>
#include <ctime>
#include <math.h>

#include "Advantra_plugin.h"
Q_EXPORT_PLUGIN2(Advantra, Advantra);

using namespace std;

static V3DLONG channel = 1;      // default channel (hardcoded)


//double round(double r)
//{
//    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
//}

// input parameter default values
string  neuritesigmas = "2,3";    // list of gaussian cross-section sigmas used for filtering and tracking the neurites
int     somaradius  = 0;       // minimum soma scale (-1 if soma not detected at all)
float   tolerance   = 10;        // FLT_MIN, minimum seed score in [0-255], average of the differences towards local neighbours
float   znccth      = 0.5;      // correlation threshold (stops tracing)
float   kappa       = 2;        // von mises circular normal distribution
int     step        = 3;        // prediction step
int     ni          = 100;       // nr. iterations
int     np          = 20;       // nr. particles
float   zdist       = 2.0;      // scaling along z
int     nodepervol  = 4;        // nodes per volume limit, used for the trace suppression (suppress over-tracing)
int     vol         = 9;        // volume patterns: 1, 5, 9, 11, 19, 27 (larger it gets, suppression is )

int     nrInputParams = 11;     // as in input_PARA struct
bool    saveMidres = true;      // save midresults

float    Kc          = 20.0;    // likelihood factor
float    neff_ratio  = 0.8;     // resampling boundary (ratio of total # pcles)

float    frangi_alfa    = .5;   //
float    frangi_beta    = .5;   //
float    frangi_C       = 500;  //
float    frangi_betaone = .5;   //
float    frangi_betatwo = 15;   //

int      MAX_TRACE_COUNT    = 5000;    // INT_MAX
float    EPSILON2           = 0.0001;  // FLT_MIN

int     REFINE_ITER          = 4;

float   SIG2RADIUS          = 1.5;      // (used for the neighbourhood estimates in refinement and grouping)
float   TRACE_RSMPL         = 1.0;      // component trace resampling step
float   GROUP_RADIUS        = 2;      // node grouping radius, defines the sampling density

bool    ENFORCE_SINGLE_TREE = true;    // single largest tree as output
int     TREE_SIZE_MIN       = 10;       // trees with less than TREE_SIZE_MIN nodes are discarded (used if ENFORCE_SINGLE_TREE=false)
int     TAIL_SIZE_MIN       = 2;        // tails (junction--endpoint) with less than TAIL_SIZE_MIN are discarded

QString NAME = "Advantra";
QString COMMENT = "";                   // will be assigned after the parameters are loaded

struct input_PARA
{
    QString inimg_file;
    V3DLONG channel;
    string  neuritesigmas;      // 1 neurite scales (comma delimited list of gaussian stds)
    int     somaradius;         // 2 min soma sigma
    float   tolerance;          // 3 find maxima tolerance (as in https://imagej.nih.gov/ij/developer/source/ij/plugin/filter/MaximumFinder.java.html)
    float   znccth;             // 4 correlation threshold
    float   kappa;              // 5 von mises kappa
    int     step;               // 6 prediction step
    int     ni;                 // 7 number of iterations
    int     np;                 // 8 number of particles
    float   zdist;              // 9 the distance between layers in pixels
    int     nodepervol;         // 10 node per volume density limit
    int     vol;                // 11 volume used for nodepervol
};

//struct Pxyz {
//    Pxyz(int x1, int y1, int z1) : x(x1),y(y1),z(z1) {}
//    float x, y, z;
//};

template<typename T>
class BfsQueue  {
public:
    std::queue<T> kk;
    BfsQueue() {}
    void enqueue(T item) {this->kk.push(item);}
    T dequeue() {T output = kk.front(); kk.pop(); return output;}
    int size(){return kk.size();}
    bool hasItems(){return !kk.empty();}
};
double round(double r)
{
    return (r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5);
}

int clampi(int x, int x1, int x2) {
    int xC = (x<x1)?x1:x;
    return (xC>x2)?x2:xC;
}

int get_undiscovered(bool * disc, int nrnodes) {
    for (int i = 0; i < nrnodes; i++) {  // go indices sorted by the correlation value
        if (!disc[i])
            return i;
    }
    return -1;
}

int get_undiscovered_weighted(vector<int> indices, bool * disc) {
    for (int i = 0; i < indices.size(); i++) {  // go indices sorted by the correlation value
        if (indices[i]!=0) {
            if (!disc[indices[i]])
                return indices[i];
        }
    }
    return -1;
}

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA, bool bmenu);
 
QStringList Advantra::menulist() const
{
	return QStringList() 
		<<tr("advantra_menu")
		<<tr("about");
}

QStringList Advantra::funclist() const
{
	return QStringList()
		<<tr("advantra_func")
		<<tr("help");
}

const QString title = QObject::tr("Advantra Plugin");

void print_help(){
    printf("*** ADVANTRA usage ***\n");
    printf("vaa3d -x Advantra -f advantra_func -i <inimg_file> -p <neuritesigmas somasigma tolerance znccth kappa step ni np zdist nodepervol vol>\n");
    printf("inimg_file                  The input image.\n");
    printf("neuritesigmas               List of comma delimited Gaussian cross-section sigmas for the neurite.\n");
    printf("somaradius                  Minimum radius of the soma (-1 to skip soma detection).\n");
    printf("tolerance                   Local maxima toleance (typically 10).\n");
    printf("znccth                      Correlation threshold (0.5).\n");
    printf("kappa                       Von Mises variance (3).\n");
    printf("step                        Prediction step (2).\n");
    printf("ni                          number of trace iterations (100).\n");
    printf("np                          Number of trace particles (20).\n");
    printf("zdist                       Z-axis dist (2).\n");
    printf("nodepervol                  Nodes per volume limit (3+).\n");
    printf("vol                         Volume chunk used for nodepp [2d]: 1, 5, 9, [3d]:1, 5, 9, 11, 19, 27.\n");
    printf("outswc_file                 Will be named automatically based on the input image file name, so you don't have to specify it.\n\n");
}

void Advantra::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("advantra_menu"))
	{
        bool bmenu = true;
        input_PARA PARA;

        // take the default params
        PARA.channel = channel;
        PARA.neuritesigmas = neuritesigmas;
        PARA.somaradius = somaradius;
        PARA.tolerance = tolerance;
        PARA.znccth = znccth;
        PARA.kappa   = kappa;
        PARA.step = step;
        PARA.ni = ni;
        PARA.np = np;
        PARA.zdist = zdist;
        PARA.nodepervol = nodepervol;
        PARA.vol = vol;

        vector<string> items;
        items.push_back("neuritesigmas");
        items.push_back("somaradius");
        items.push_back("tolerance");
        items.push_back("znccth");
        items.push_back("kappa");
        items.push_back("step");
        items.push_back("ni");
        items.push_back("np");
        items.push_back("zdist");
        items.push_back("nodepervol");
        items.push_back("vol");

        // initialization
        vector<string> inits;
        inits.push_back(PARA.neuritesigmas.c_str());
        inits.push_back(QString::number(PARA.somaradius).toStdString().c_str());
        inits.push_back(QString::number(PARA.tolerance).toStdString().c_str());
        inits.push_back(QString::number(PARA.znccth).toStdString().c_str());
        inits.push_back(QString::number(PARA.kappa).toStdString().c_str());
        inits.push_back(QString::number(PARA.step).toStdString().c_str());
        inits.push_back(QString::number(PARA.ni).toStdString().c_str());
        inits.push_back(QString::number(PARA.np).toStdString().c_str());
        inits.push_back(QString::number(PARA.zdist).toStdString().c_str());
        inits.push_back(QString::number(PARA.nodepervol).toStdString().c_str());
        inits.push_back(QString::number(PARA.vol).toStdString().c_str());

        CommonDialog dialog(items, inits);
        dialog.setWindowTitle(title);
        if(dialog.exec() != QDialog::Accepted) return;

        PARA.neuritesigmas = dialog.get_para("neuritesigmas");
        dialog.get_num("somaradius", PARA.somaradius);
        dialog.get_num("tolerance", PARA.tolerance);
        dialog.get_num("znccth", PARA.znccth);
        dialog.get_num("kappa", PARA.kappa);
        dialog.get_num("step", PARA.step);
        dialog.get_num("ni", PARA.ni);
        dialog.get_num("np", PARA.np);
        dialog.get_num("zdist", PARA.zdist);
        dialog.get_num("nodepervol", PARA.nodepervol);
        dialog.get_num("vol", PARA.vol);

        // check input
//        if(PARA.neuritesigmas){v3d_msg(QObject::tr("neuritesigmas out of range")); return;}
        if(PARA.somaradius<0){v3d_msg(QObject::tr("somaradius out of range")); return;}
        if(PARA.tolerance<0){v3d_msg(QObject::tr("tolerance out of range")); return;}
        if(PARA.znccth<0 || PARA.znccth>1){v3d_msg(QObject::tr("znccth out of range")); return;}
        if(PARA.kappa<0 || PARA.kappa>5){v3d_msg(QObject::tr("kappa out of range")); return;}
        if(PARA.step<1){v3d_msg(QObject::tr("step out of range")); return;}
        if(PARA.ni<=0){v3d_msg(QObject::tr("ni out of range")); return;}
        if(PARA.np<=0){v3d_msg(QObject::tr("np out of range")); return;}
        if(PARA.zdist<1){v3d_msg(QObject::tr("zdist out of range")); return;}
        if(PARA.nodepervol<=2 || PARA.nodepervol>20){v3d_msg(QObject::tr("nodepervol out of range")); return;}
        if(!(PARA.vol==1 || PARA.vol==5 || PARA.vol==9 || PARA.vol==11 || PARA.vol==19 || PARA.vol==27)){v3d_msg(QObject::tr("vol can be 1,5,9,11,19,27")); return;}

        reconstruction_func(callback,parent,PARA,bmenu);

	}
	else
	{
        v3d_msg(tr("PNR: automatic neuron reconstruction from microscopy image stacks\n"
            "miro@braincadet.com\n"
                   "2015-8-19 first version (BigNeuron submission)\n"
                   "2016-3-29 complete re-design\n"
                   "2017-9-26 isbi 17 submission"));
	}
}

bool Advantra::dofunc(const QString & func_name, const V3DPluginArgList & input, V3DPluginArgList & output, V3DPluginCallback2 & callback,  QWidget * parent)
{
	if (func_name == tr("advantra_func"))
	{
        bool bmenu = false;
        input_PARA PARA;

        vector<char*> * pinfiles = (input.size() >= 1) ? (vector<char*> *) input[0].p : 0;
        vector<char*> * pparas = (input.size() >= 2) ? (vector<char*> *) input[1].p : 0;
        vector<char*> infiles = (pinfiles != 0) ? * pinfiles : vector<char*>();
        vector<char*> paras = (pparas != 0) ? * pparas : vector<char*>();

        if(infiles.empty())
        {
            fprintf (stderr, "Need input image. \n");
            return false;
        }
        else
            PARA.inimg_file = infiles[0];

        // constrain number of input parameters
        if (paras.size()!=nrInputParams) {
            fprintf (stderr, "\nNeeds %d input parameters.\n\n", nrInputParams);
            print_help();
            return false;
        }

        int k=0;
        PARA.channel        = channel; // hardcoded not submitted as parameter
        PARA.neuritesigmas  = (paras.size() >= k+1)   ? paras[k]                    : neuritesigmas;       k++;
        PARA.somaradius     = (paras.size() >= k+1)   ? atoi(paras[k])              : somaradius;       k++;
        PARA.tolerance      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : tolerance;    k++;
        PARA.znccth     = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : znccth;       k++;
        PARA.kappa      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : kappa;        k++;
        PARA.step       = (paras.size() >= k+1)   ? atoi(paras[k])              : step;         k++;
        PARA.ni         = (paras.size() >= k+1)   ? atoi(paras[k])              : ni;           k++;
        PARA.np         = (paras.size() >= k+1)   ? atoi(paras[k])              : np;           k++;
        PARA.zdist      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : zdist;        k++;
        PARA.nodepervol = (paras.size() >= k+1)   ? atoi(paras[k])              : nodepervol;   k++;
        PARA.vol        = (paras.size() >= k+1)   ? atoi(paras[k])              : vol;          k++;

        // check user input
//        if(PARA.neuritesigmas){v3d_msg(QObject::tr("neuritesigmas out of range")); return 0;}
        if(PARA.somaradius<0){v3d_msg(QObject::tr("somaradius out of range")); return 0;}
        if(PARA.tolerance<0){v3d_msg(QObject::tr("tolerance out of range")); return 0;}
        if(PARA.znccth<0 || PARA.znccth>1){v3d_msg(QObject::tr("znccth out of range")); return 0;}
        if(PARA.kappa<0 || PARA.kappa>5){v3d_msg(QObject::tr("kappa out of range")); return 0;}
        if(PARA.step<1){v3d_msg(QObject::tr("step out of range")); return 0;}
        if(PARA.ni<=0){v3d_msg(QObject::tr("ni out of range")); return 0;}
        if(PARA.np<=0){v3d_msg(QObject::tr("np out of range")); return 0;}
        if(PARA.zdist<1){v3d_msg(QObject::tr("zdist out of range")); return 0;}
        if(PARA.nodepervol<=2 || PARA.nodepervol>20){v3d_msg(QObject::tr("nodepervol out of range")); return 0;}
        if(!(PARA.vol==1 || PARA.vol==5 || PARA.vol==9 || PARA.vol==11 || PARA.vol==19 || PARA.vol==27)){v3d_msg(QObject::tr("vol can be 1,5,9,11,19,27")); return 0;}
        reconstruction_func(callback,parent,PARA,bmenu);
	}
    else if (func_name == tr("help"))
    {
        ////HERE IS WHERE THE DEVELOPERS SHOULD UPDATE THE USAGE OF THE PLUGIN
        print_help();
	}
	else return false;

	return true;
}

class CompareIndicesByNodeCorrVal {
    vector<Node>* nd;
public:
    CompareIndicesByNodeCorrVal(std::vector<Node>* values) : nd(values) {}
public:
    bool operator() (const int& a, const int& b) const { return (*nd)[a].corr > (*nd)[b].corr; }
};

class CompareSeedCorr {
    vector<seed>* _s;
public:
    CompareSeedCorr(vector<seed>* v) : _s(v) {}
    bool operator() (const int& a, const int& b) const { return (*_s)[a].corr > (*_s)[b].corr; }
};

/*
bool ischecked(seed _s, bool* _smap, int _w, int _h, int _l) {

    int x = round(_s.x);
    if (x<0 || x>=_w) return true;

    int y = round(_s.y);
    if (y<0 || y>=_h) return true;

    int z = round(_s.z);
    if (z<0 || z>=_l) return true;

    return _smap[z*_w*_h+y*_w+x];
}
*/

int get_undiscovered2(vector<int> dist){
    for (int i = 1; i < dist.size(); i++) {
        if (dist[i]==INT_MAX) {
           return i;
        }
    }
    return -1;
}

void bfs2(vector<Node> nlist, vector<Node>& tree, bool remove_isolated_tree_with_one_node) {

    /*
    https://en.wikipedia.org/wiki/Breadth-first_search
    1 Breadth-First-Search(Graph, root):
    2
    3     for each node n in Graph:
    4         n.distance = INFINITY
    5         n.parent = NIL
    6
    7     create empty queue Q
    8
    9     root.distance = 0
    10     Q.enqueue(root)
    11
    12     while Q is not empty:
    13
    14         current = Q.dequeue()
    15
    16         for each node n that is adjacent to current:
    17             if n.distance == INFINITY:
    18                 n.distance = current.distance + 1
    19                 n.parent = current
    20                 Q.enqueue(n)
    */

    BfsQueue<int> q;

    vector<int> dist(nlist.size());
    vector<int> nmap(nlist.size());
    vector<int> parent(nlist.size());

    for (int i = 0; i < nlist.size(); ++i) {
        dist[i] = INT_MAX;
        nmap[i] = -1; // indexing in output tree
        parent[i] = -1; // parent index in current tree
    }

    dist[0] = -1;

    Node tree0(nlist[0]); // first element of the nodelist is dummy both in input and output
    tree.clear();
    tree.push_back(tree0);
    int treecnt = 0; // independent tree counter, will be obsolete

    int seed;

    while ((seed = get_undiscovered2(dist))>0) {

        treecnt++;

        dist[seed] = 0;
        nmap[seed] = -1;
        parent[seed] = -1;
        q.enqueue(seed);

        int nodesInTree = 0;

        while (q.hasItems()) {

            // dequeue(), take from FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
            int curr = q.dequeue();

            Node n(nlist[curr]);
            n.nbr.clear();
            if (n.type!=Node::SOMA) n.type = treecnt+2; // vaa3d viz

            if (parent[curr] > 0) n.nbr.push_back(nmap[parent[curr]]);

            nmap[curr] = tree.size();
            tree.push_back(n);
            nodesInTree++;

            // for each node adjacent to current
            for (int j = 0; j < nlist[curr].nbr.size(); j++) {

                int adj = nlist[curr].nbr[j];

                if (dist[adj] == INT_MAX) {
                    dist[adj] = dist[curr] + 1;
                    parent[adj] = curr;
                    // enqueue(), add to FIFO structure, http://en.wikipedia.org/wiki/Queue_%28abstract_data_type%29
                    q.enqueue(adj);
                }

            }

            // check if there were any neighbours
            if (nodesInTree==1 && !q.hasItems() && remove_isolated_tree_with_one_node) {
                tree.pop_back();                // remove the one that was just added
                nmap[curr] = -1;                // cancel the last entry
            }

        }

    }

    cout << treecnt << " trees" << endl;

}

void save_nodelist(vector<Node> nlist, QString swcname, int type=-1, float sig2r=1, QString name="", QString comment=""){

    // NOTE: nodes have more than one parent, export will duplicate node ids to show all neighbouring connections (nodes have bidirectional links)

    NeuronTree recnodes;

    if (name.compare("")!=0) recnodes.name = name;

    if (comment.compare("")!=0) recnodes.comment = comment;

    int cnt_recnodes = 0;

    for (int i = 1; i < nlist.size(); ++i) {

        cnt_recnodes++;

        if (nlist[i].nbr.size()==0) {
            NeuronSWC n;
            n.n = n.nodeinseg_id = cnt_recnodes;
            n.type = (type==-1)?nlist[i].type:type;
            n.x = nlist[i].x;
            n.y = nlist[i].y;
            n.z = nlist[i].z;
            n.r = sig2r * nlist[i].sig;
            n.parent = -1;
            recnodes.listNeuron.append(n);
        }
        else {
            for (int j = 0; j < nlist[i].nbr.size(); ++j) { // same id will repeat as many times as there are neighbour links
                NeuronSWC n;
                n.n = n.nodeinseg_id = cnt_recnodes;
                n.type = (type==-1)?nlist[i].type:type;
                n.x = nlist[i].x;
                n.y = nlist[i].y;
                n.z = nlist[i].z;
                n.r = sig2r * nlist[i].sig;
                n.parent = nlist[i].nbr[j];
                recnodes.listNeuron.append(n);
            }
        }
    }

    writeSWC_file(swcname.toStdString().c_str(), recnodes);
}

vector<Node> compute_trees(vector<Node> nlist) {
    // convert nodelist (with bidirectional links, grouped) into the treelist (each node has one parent), both are vector<Node>
    vector<Node> ntree;
    bfs2(nlist, ntree, true); // bfs tree extraction, tree will have 0 or 1 in nbr list (it's formatted as tree)
    return ntree;
}

/*
void summarize_tree(vector<Node> ntree) {
    cout << "|ntree|=" << ntree.size() << endl;
    int root_prev=1, root_curr = 1;
    for (int i = 1; i <= ntree.size(); ++i) {
        if (i==ntree.size() || ntree[i].nbr.size()==0) { // root nodes got no neighbors
            root_prev = root_curr;
            root_curr = i;
            cout << "|t|=" << (root_curr-root_prev) << " [" <<root_prev << "," << root_curr << "] " << endl;
        }
    }
}
*/

vector<Node> extract_largest_tree(vector<Node> ntreeX) {
    // take the input TREELIST and remove all independent trees but the largest one, tag those that are to be kept
    int root_curr=1, root_prev=1;
    int tree_max_size = -INT_MAX, tree_max_beg=-INT_MAX, tree_max_end=-INT_MAX;

    for (int i = 1; i <= ntreeX.size(); ++i) {
        if (i==ntreeX.size() || ntreeX[i].nbr.size()==0) { // tree root will have zero neighbours (check bfs2 implementation for tree traveral standard)
            root_prev = root_curr;
            root_curr = i;
            if (root_curr-root_prev>tree_max_size) {
                tree_max_size = root_curr-root_prev;
                tree_max_beg = root_prev;
                tree_max_end = root_curr;
            }
        }
    }

    vector<bool> keep_map(ntreeX.size(), false);
    keep_map[0]=true;
    for (int j = tree_max_beg; j < tree_max_end; ++j) {
        keep_map[j] = true;
    }

    // go through the tree once more and concatenate together those that are kept based on the keep_map[]
    vector<int> X2Y(ntreeX.size(),-1);
    vector<Node> ntreeY;
    for (int i = 0; i < ntreeX.size(); ++i) {
        if (keep_map[i]) { // keep_map[0]=true
            X2Y[i] = ntreeY.size();
            Node nYi(ntreeX[i]);
            ntreeY.push_back(nYi); // nodes are copied, neighbours and the type stay as in ntreeX vector indexing
        }
    }

    // remap the linking
    for (int i = 1; i < ntreeY.size(); ++i) {
        for (int j = 0; j < ntreeY[i].nbr.size(); ++j) {
            ntreeY[i].nbr[j] = X2Y[     ntreeY[i].nbr[j]     ];
        }
    }

    return ntreeY;

}

vector<Node> extract_trees(vector<Node> ntreeX, int min_size) {
    // take input treelist ntreeX and remove all trees with less than min_size nodes, removed nodes belong to one connected tree, tag those that are to be removed

    vector<bool> remove_map(ntreeX.size(), false); // none is assumed to be removed at the beginning

    int root_curr=1, root_prev=1;
    for (int i = 1; i <= ntreeX.size(); ++i) { // first one is not removed, i starts from 1
        if (i==ntreeX.size() || ntreeX[i].nbr.size()==0) { // tree root will have zero neighbours (check bfs2 implementation for tree traveral standard)
            root_prev = root_curr;
            root_curr = i;
            if (root_curr-root_prev<min_size) { // current tree had less nodes than min_size => mark it for removal
                for (int j = root_prev; j < root_curr; ++j) {
                    remove_map[j] = true;
                }
            }
        }
    }

    // go through the tree once more and concatenate together those that are not removed based on the remove_map[]
    vector<int> X2Y(ntreeX.size(),-1);
    vector<Node> ntreeY;
    for (int i = 0; i < ntreeX.size(); ++i) {
        if (!remove_map[i]) { // remove_map[0]=false
            X2Y[i] = ntreeY.size();
            Node nYi(ntreeX[i]);
            ntreeY.push_back(nYi); // nodes are copied, neighbours stay as in ntreeX vector indexing
        }
    }

    // remap the linking
    for (int i = 1; i < ntreeY.size(); ++i) {
        for (int j = 0; j < ntreeY[i].nbr.size(); ++j) {
            ntreeY[i].nbr[j] = X2Y[     ntreeY[i].nbr[j]     ];
        }
    }

    return ntreeY;

}

vector<Node> remove_tails(vector<Node> ntreeX, int min_size) {
    // take input TREELIST, remove all tails with less than min_size nodes, node removal + link towards tail removal
    vector<Node> ntreeX_2d; // copy of ntreeX with bidirectional links
    ntreeX_2d = ntreeX;
    for (long i = 1; i < ntreeX_2d.size(); ++i) {
        for (int j = 0; j < ntreeX_2d[i].nbr.size(); ++j) { // existing link i -- ntreeX_2d[i].nbr[j]
            ntreeX_2d[   ntreeX_2d[i].nbr[j]   ].nbr.push_back(      i       ); // add reverse link
        }
    }

    vector<int> Nnbr(ntreeX_2d.size(), 0);
    for (long i = 1; i < ntreeX_2d.size(); ++i) {
        Nnbr[i] = ntreeX_2d[i].nbr.size();
    }

    vector<bool> remove_map(ntreeX_2d.size(), false);

    for (long i = 1; i < ntreeX_2d.size(); ++i) {
        if (ntreeX_2d[i].nbr.size()==1) { // endpoint
            vector<int> tail_i;
            tail_i.push_back(i);
            int tail_next = ntreeX_2d[i].nbr[0];

            while (ntreeX_2d[  tail_next  ].nbr.size()==2) {
                tail_i.push_back(tail_next);
                tail_next = (ntreeX_2d[tail_next].nbr[0] == tail_i[tail_i.size()-2])? ntreeX_2d[  tail_next  ].nbr[1] : ntreeX_2d[  tail_next  ].nbr[0];
            }

//            cout<<"end("<<i<<"): tail.size="<<tail_i.size()<< "(min="<< min_size<<"), reached "<< ntreeX_2d[tail_next].nbr.size()<<",idx="<<tail_next<<" "<<flush;

            // get the size of the index list
            if (ntreeX_2d[tail_next].nbr.size()>2 && tail_i.size()<min_size) { // tail if reached junction and there are less than min_size nodes
//                cout<<" Remove"<<endl;
                for (int j = 0; j < tail_i.size(); ++j) { // mark for removal
                    remove_map[tail_i[j]] = true;
                }
            }
//            else{
//                cout<<endl;
//            }
        }
    }

    // go through the tree once more and concatenate together those that are not removed based on the remove_map[]
    vector<int> X2Y(ntreeX.size(),-1);
    vector<Node> ntreeY;
    for (int i = 0; i < ntreeX.size(); ++i) {
        if (!remove_map[i]) { // remove_map[0]=false
            X2Y[i] = ntreeY.size();
            Node nYi(ntreeX[i]);
            ntreeY.push_back(nYi); // nodes are copied, neighbours stay as in ntreeX vector indexing (so the node type is unchanged)
        }
    }

    // remap the neighbor linking
    for (int i = 1; i < ntreeY.size(); ++i) {
        for (int j = ntreeY[i].nbr.size()-1; j >= 0; --j) {
            if (remove_map[ntreeY[i].nbr[j]]) {
                ntreeY[i].nbr.erase(ntreeY[i].nbr.begin()+j); // remove the link towards removed tail
            }
            else {
                ntreeY[i].nbr[j] = X2Y[ntreeY[i].nbr[j]]; // and X2Y>=0 since it was not removed
            }
        }
    }

    return ntreeY;

}

template<typename T>
void save_vector(vector<T> v, QString savename) {
    ofstream f;
    f.open(savename.toStdString().c_str(), ofstream::out | ofstream::trunc);

    for (int i = 0; i < v.size(); ++i) {
        f << v[i] << ((i<v.size()-1)?",":"\n") << flush;
    }

    f.close();
    cout << " exported:\t" << savename.toStdString() << endl;
}

void interpolate_treelist(vector<Node>& ntree, float step, int type=-1) {

    // assume onedirectional connections (1 link between 2 nodes in 1 direction)

    // interpolate all inter-node connections with the step size
    // not necessary to have bookkeeping variable, as there are no bidirectional links
    long init_size = ntree.size();

    for (long i = 1; i < init_size; ++i) {

        // change the type
        if (type>=0) { // type is not default (-1)
            // assign non-soma nodes with the given type
            if (ntree[i].type!=Node::SOMA) {
                ntree[i].type = type;
            }
        }

        // interpolation
        for (int j = 0; j < ntree[i].nbr.size(); ++j) { // there should be 0 or 1 neighbor
            long i1 = ntree[i].nbr[j];

            // interpolate between ntree[i] and ntree[i1]

//            if (ntree[i].type==Node::SOMA || ntree[i1].type==Node::SOMA) continue;

            float vnorm = sqrt(pow(ntree[i1].x-ntree[i].x,2) + pow(ntree[i1].y-ntree[i].y,2) + pow(ntree[i1].z-ntree[i].z,2));
            float vx = (ntree[i1].x-ntree[i].x)/vnorm;
            float vy = (ntree[i1].y-ntree[i].y)/vnorm;
            float vz = (ntree[i1].z-ntree[i].z)/vnorm;
            int N = ceil(vnorm/step);

            // add subsampling TODO here interpolate the tree
            for (int k = 1; k < N; ++k) {
                // add the node,only location is used in the refinement stage currently
                Node nXadd(ntree[i].x+k*(vnorm/N)*vx,
                           ntree[i].y+k*(vnorm/N)*vy,
                           ntree[i].z+k*(vnorm/N)*vz,
                           vx,vy,vz,
                           ntree[i].corr + (ntree[i1].corr-ntree[i].corr) * (k/(float)N), // 0.5*(ntree[i].corr+ntree[i1].corr),
                           ntree[i].sig  + (ntree[i1].sig -ntree[i].sig) * (k/(float)N), // 0.5*(ntree[i].sig+ntree[i1].sig),
                           (k<=N/2)?ntree[i].type:ntree[i1].type);
                ntree.push_back(nXadd);

                // link backward
                if (k==1) { // first
                    ntree[i].nbr[j] = ntree.size()-1;
//                    ntreeX[ntreeX.size()-1].nbr.push_back(i);
//                    ntreeX[i].nbr[j] = ntreeX.size()-1; // replace i1 with the link to the first addition
                }
                else { // middle
//                    ntreeX[ntreeX.size()-1].nbr.push_back(ntreeX.size()-2);
                    ntree[ntree.size()-2].nbr.push_back(ntree.size()-1);
                }
                // link forward
                if (k==N-1) { // last
                    ntree[ntree.size()-1].nbr.push_back(i1);
//                    ntreeX[i1].nbr[j1] = ntreeX.size()-1; // replace i with the link to the last addition
                }

            }
        }
    }

}

void interpolate_nodelist(vector<Node>& nX, float step) {

    // assumes bidirectional connections (2 links between 2 nodes in 2 possible directions)

    // interpolate all inter-node links with the step size
    vector< vector<bool> > chk(nX.size()); // disable accessing the same pair of nodes twice
    for (long i = 0; i < nX.size(); ++i) {

        // assume all links would be interpolated -- set to false
        vector<bool> chk1(nX[i].nbr.size(), false); // if links are 'checked' then they won't be resampled

//        for (int j = 0; j < nX[i].nbr.size(); ++j) {
//            // if one node of the link is SOMA, then set that link as checked so that it's not interpolated
//            if (nX[i].type==Node::SOMA || nX[nX[i].nbr[j]].type==Node::SOMA) {
//                chk1[j] = true; // disables from the interpolation
//            }
//        }

        chk[i] = chk1;

    }

    long init_size = nX.size();

    for (long i = 1; i < init_size; ++i) {
        for (int j = 0; j < nX[i].nbr.size(); ++j) {
            if (!chk[i][j]) {

                long i1 = nX[i].nbr[j];
                int j1 = find(nX[i1].nbr.begin(), nX[i1].nbr.end(), i) - nX[i1].nbr.begin();

                if (j1<nX[i1].nbr.size()) { // interpolate if there was existing link back

                    chk[i][j]       = true; // mark both links as checked so that they are not interpolated again
                    chk[i1][j1]     = true;

                    float vnorm = sqrt(pow(nX[i1].x-nX[i].x,2) + pow(nX[i1].y-nX[i].y,2) + pow(nX[i1].z-nX[i].z,2));

                    float vx = (nX[i1].x-nX[i].x)/vnorm;
                    float vy = (nX[i1].y-nX[i].y)/vnorm;
                    float vz = (nX[i1].z-nX[i].z)/vnorm;
                    int N = ceil(vnorm/step);

                    // add subsampling if N>1
                    for (int k = 1; k < N; ++k) {
                        // add the node,only location is used in the refinement stage currently
                        Node nXadd(nX[i].x+k*(vnorm/N)*vx,nX[i].y+k*(vnorm/N)*vy,nX[i].z+k*(vnorm/N)*vz,
                                   vx,vy,vz,
                                   nX[i].corr + (nX[i1].corr-nX[i].corr) * (k/(float)N), //  0.5*(nX[i].corr+nX[i1].corr),    // intepolate correlation
                                   nX[i].sig  + (nX[i1].sig -nX[i].sig) * (k/(float)N), // 0.5*(nX[i].sig+nX[i1].sig),      // interpolate sigma
                                   (k<=N/2)?nX[i].type:nX[i1].type);
                        nX.push_back(nXadd);

                        // link backward
                        if (k==1) {
                            // first: link nX[nX.size()-1] with nX[i]
                            nX[nX.size()-1].nbr.push_back(i);
                            nX[i].nbr[j] = nX.size()-1; // replace i1 with the link to the first addition
                        }
                        else {
                            // middle: link nX[nX.size()-1] with nX[nX.size()-2]
                            nX[nX.size()-1].nbr.push_back(nX.size()-2);
                            nX[nX.size()-2].nbr.push_back(nX.size()-1);
                        }

                        // link forward
                        if (k==N-1) {
                            // last: link nX[nX.size()-1] with nX[i1]
                            nX[nX.size()-1].nbr.push_back(i1);
                            nX[i1].nbr[j1] = nX.size()-1; // replace i with the link to the last addition
                        }

                    }

                }
            }
        }
    }

    cout << ((float)nX.size()/init_size)*100.0 << "% node # after interpolation" << endl;

}

void refine_blurring(vector<Node> nX, vector<Node>& nY, float SIG2RAD, int MAXITER, float EPSILON2) {
    // TODO blurring should work for all node types!! including soma
    int checkpoint = round(nX.size()/10.0);

    vector <vector<float> >  conv(nX.size(),vector<float>(4)); // fix allocation as it can happen that size exceeds int range
    vector <vector<float> >  next(nX.size(),vector<float>(4));
//vector<vector<float> > ray_x(ray_numbers_2d,vector<float>(100));


//    float conv[nX.size()][4]; // fix allocation as it can happen that size exceeds int range
//    float next[nX.size()][4];

    for (long i = 1; i < nX.size(); ++i) {
        conv[i][0] = nX[i].x;
        conv[i][1] = nX[i].y;
        conv[i][2] = nX[i].z;
        conv[i][3] = nX[i].sig;
    }

    float x2, y2, z2, d2, d2_max, r2;
    int iter = 0;
    int cnt = 0;

    do {

        cout << "iter = " << iter << endl;

        // each iteration will shift the whole sample = blurring
        d2_max = -FLT_MAX;

        for (long i = 1; i < nX.size(); ++i) {

            if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

            if (nX[i].type==Node::SOMA) { // skip the ms iterations for soma nodes
                next[i][0] = conv[i][0];
                next[i][1] = conv[i][1];
                next[i][2] = conv[i][2];
                next[i][3] = conv[i][3];
            }
            else { // do not use soma nodes for averaging the trace nodes

                cnt = 0;  // count neighbours

                next[i][0] = 0;
                next[i][1] = 0;
                next[i][2] = 0;
                next[i][3] = 0;

                r2 = pow(SIG2RAD * conv[i][3],2);

                for (long j = 1; j < nX.size(); ++j) {
                    if (nX[j].type!=Node::SOMA) {
                        x2 = pow(conv[j][0]-conv[i][0],2);
                        if (x2<=r2) {
                            y2 = pow(conv[j][1]-conv[i][1],2);
                            if (x2+y2<=r2) {
                                z2 = pow(conv[j][2]-conv[i][2],2);
                                if (x2+y2+z2<=r2) {
                                    next[i][0] += conv[j][0];
                                    next[i][1] += conv[j][1];
                                    next[i][2] += conv[j][2];
                                    next[i][3] += conv[j][3];
                                    cnt++;
                                }
                            }
                        }
                    }
                }

                next[i][0] /= cnt;
                next[i][1] /= cnt;
                next[i][2] /= cnt;
                next[i][3] /= cnt;

            }

            d2 = pow(next[i][0]-conv[i][0],2) + pow(next[i][1]-conv[i][1],2) + pow(next[i][2]-conv[i][2],2);
            if (d2>d2_max)
                d2_max = d2;

            conv[i][0] = next[i][0];
            conv[i][1] = next[i][1];
            conv[i][2] = next[i][2];
            conv[i][3] = next[i][3];

        }

        cout << endl;

        iter++;

    }
    while(iter<MAXITER && d2_max>EPSILON2);

    nY.clear();
    nY = nX;

    for (long i = 1; i < nY.size(); ++i) {
        nY[i].x = conv[i][0];
        nY[i].y = conv[i][1];
        nY[i].z = conv[i][2];
        nY[i].sig = conv[i][3];
    }

}

void non_blurring(vector<Node> nX, vector<Node>& nY, float SIG2RAD, int MAXITER, float EPSILON2) {

    // mean-shift (non-blurring) uses flexible neighbourhood scaled with respect to the node's sigma

    int checkpoint = round(nX.size()/10.0);

    float conv[4], next[4]; // x y z sig

    nY.clear();
    nY = nX;

    float x2, y2, z2, d2, r2;
    int iter, cnt;

    for (long i = 1; i < nY.size(); ++i) {

        if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

        // refine nX[i] node location and scale and store the result in nY[i]
//        if (nY[i].type==Node::SOMA) continue; // do not refine soma nodes

        conv[0] = nX[i].x;
        conv[1] = nX[i].y;
        conv[2] = nX[i].z;
        conv[3] = nX[i].sig;

        iter = 0;

        do {
            cnt = 0;

            next[0] = 0;        // local mean is the follow-up location
            next[1] = 0;
            next[2] = 0;
            next[3] = 0;

            r2 = pow(SIG2RAD * conv[3],2);

            for (long j = 1; j < nX.size(); ++j) {
//                if (nX[j].type==Node::SOMA) continue; // don't use soma nodes in refinement
                x2 = pow(nX[j].x-conv[0],2);
                if (x2<=r2) {
                    y2 = pow(nX[j].y-conv[1],2);
                    if (x2+y2<=r2) {
                        z2 = pow(nX[j].z-conv[2],2);
                        if (x2+y2+z2<=r2) {
                            next[0] += nX[j].x;
                            next[1] += nX[j].y;
                            next[2] += nX[j].z;
                            next[3] += nX[j].sig;
                            cnt++;
                        }
                    }
                }
            }

//            if (cnt==0) cout << "WRONG!!!" << endl;

            next[0] /= cnt; // cnt > 0, at least node location itself will be in the kernel neighbourhood
            next[1] /= cnt;
            next[2] /= cnt;
            next[3] /= cnt;

            d2 = pow(next[0]-conv[0],2) + pow(next[1]-conv[1],2) + pow(next[2]-conv[2],2);

            conv[0] = next[0]; // for the next iteration
            conv[1] = next[1];
            conv[2] = next[2];
            conv[3] = next[3];

            iter++;

        }
        while (iter<MAXITER && d2>EPSILON2);

        nY[i].x     = conv[0];
        nY[i].y     = conv[1];
        nY[i].z     = conv[2];
        nY[i].sig   = conv[3];

    } // go through nY[i], initiate with nX[i] values and refine by mean-shift averaging

    cout << endl;

}

bool is_cross_section(float n1x, float n1y, float n1z, float v1x, float v1y, float v1z, float n2x, float n2y, float n2z, float d_axial, float d_radial=-1) {

    // axial distance == dot produuct with (v1x, v1y, v1z), projection onto direction
    float da = (n2x-n1x)*v1x + (n2y-n1y)*v1y + (n2z-n1z)*v1z;
    if (fabs(da)<=d_axial) {

        if (d_radial<0) {
            // don't consider radial distance
            return true;
        }
        else {
            // consider radial distance == point to line distance
            float d_radial2 = pow(d_radial, 2);
            float dx2 = pow(n1x-n2x+d_axial*v1x,2); // x component of the distance vector
            if (dx2<=d_radial2) {
                float dy2 = pow(n1y-n2y+d_axial*v1y,2);
                if (dx2+dy2<=d_radial2) {
                    float dz2 = pow(n1z-n2z+d_axial*v1z,2);
                    if (dx2+dy2+dz2<=d_radial2) {
                        return true;
                    }
                    else
                        return false;
                }
                else
                    return false;
            }
            else
                return false;
        }
    }
    else
        return false;

}

/*
//float refine1(vector<Node> n0, vector<Node>& nI, float krad2){

//    nI.clear();
//    nI = n0;
//    int checkpoint = round(nI.size()/10.0);

//    float nbrx[n0.size()]; // allocate storage for nbr. coords
//    float nbry[n0.size()];
//    float nbrz[n0.size()];

//    float mu[3];    // c_t
//    float nconv[3]; // a_t-1
//    float nnext[3]; // a_t

//    double cov[3][3]; // eigen analysis values
//    double vec[3][3];
//    double eig[3];

//    float k, dx, dy, dz;
//    float discr2, discr2_max = -FLT_MAX;

//    for (int i = 1; i < nI.size(); ++i) {

//        if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

//        int iter = 0; // stop condition
//        float d2 = 0;

//        nconv[0] = nI[i].x; // initialize nconv with refined node coordinates
//        nconv[1] = nI[i].y;
//        nconv[2] = nI[i].z;

//        do {
//            int cnt = 0;        // count neighbours
//            float x2, y2, z2;
//            mu[0] = 0;
//            mu[1] = 0;
//            mu[2] = 0;

//            for (int j = 1; j < n0.size(); ++j) {
//                if (j!=i) {
//                    x2 = pow(n0[j].x-nconv[0],2);
//                    if (x2<=krad2) {
//                        y2 = pow(n0[j].y-nconv[1],2);
//                        if (x2+y2<=krad2) {
//                            z2 = pow(n0[j].z-nconv[2],2);
//                            if (x2+y2+z2<=krad2) {

//                                nbrx[cnt] = n0[j].x;
//                                mu[0] += n0[j].x;

//                                nbry[cnt] = n0[j].y;
//                                mu[1] += n0[j].y;

//                                nbrz[cnt] = n0[j].z;
//                                mu[2] += n0[j].z;

//                                cnt++;
//                            }
//                        }
//                    }
//                }
//            }

//            if (cnt>0) {

//                cov[0][0]=0;cov[0][1]=0;cov[0][2]=0;
//                cov[1][0]=0;cov[1][1]=0;cov[1][2]=0;
//                cov[2][0]=0;cov[2][1]=0;cov[2][2]=0;

//                for (int ii = 0; ii < cnt; ++ii) {
//                    dx = nbrx[ii]-mu[0];
//                    dy = nbry[ii]-mu[1];
//                    dz = nbrz[ii]-mu[2];
//                    cov[0][0] += dx*dx;
//                    cov[0][1] += dx*dy;
//                    cov[0][2] += dx*dz;
//                    cov[1][0] += dy*dx;
//                    cov[1][1] += dy*dy;
//                    cov[1][2] += dy*dz;
//                    cov[2][0] += dz*dx;
//                    cov[2][1] += dz*dy;
//                    cov[2][2] += dz*dz;
//                }

//                if(cov[0][0]+cov[1][1]+cov[2][2]>FLT_MIN) {
//                    Frangi::eigen_decomposition_static(cov, vec, eig);
//                    // mu[], vec[][2], and nconv[] are used to calculate nnext[]
//                    k = (nconv[0]-mu[0])*vec[0][2] + (nconv[1]-mu[1])*vec[1][2] + (nconv[2]-mu[2])*vec[2][2];
//                    nnext[0] = mu[0] + k * vec[0][2];
//                    nnext[1] = mu[1] + k * vec[1][2];
//                    nnext[2] = mu[2] + k * vec[2][2];
//                }
//                else {
//                    // all neighbours the same, no directionality
//                    nnext[0] = nconv[0];
//                    nnext[1] = nconv[1];
//                    nnext[2] = nconv[2];
//                }
//            }
//            else { // 0 nbrs
//                nnext[0] = nconv[0];
//                nnext[1] = nconv[1];
//                nnext[2] = nconv[2];
//            }

//            d2 = pow(nnext[0]-nconv[0],2) + pow(nnext[1]-nconv[1],2) + pow(nnext[2]-nconv[2],2);

//            nconv[0] = nnext[0];
//            nconv[1] = nnext[1];
//            nconv[2] = nnext[2];

//            iter++;
//        }
//        while(iter<MAXITER && d2>EPSILON2);

//        discr2 = pow(nconv[0]-nI[i].x,2) + pow(nconv[1]-nI[i].y,2) + pow(nconv[2]-nI[i].z,2);
//        if (discr2>discr2_max) discr2_max = discr2;

//        nI[i].x = nconv[0];
//        nI[i].y = nconv[1];
//        nI[i].z = nconv[2];

//    }

//    return discr2_max;

//}

//float refine(vector<Node> n0, vector<Node>& nI, bool img3d, float daxial) {

//    nI.clear();
//    nI = n0; // initialize

//    vector< vector<float> > proj(0, vector<float>((img3d)?2:1) ); // ((is2d)?vector<float>(1):vector<float>(2))
//    vector<float> proj_new((img3d)?2:1);
//    vector<float> proj_conv((img3d)?2:1);

//    float ux, uy, uz, wx, wy, wz;

//    float discr2_max = -FLT_MAX; // output score will be max squared discrepancy established in the refinement stage
//    float discr2;

//    int checkpoint = round(nI.size()/10.0);

//    for (int i = 1; i < nI.size(); ++i) { // INT_MAX

//        if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

////        proj.erase(proj.begin()+1, proj.begin()+proj.size()); // erase all but first element
//        proj.clear();

//        SeedExtractor::orthogonals(nI[i].vx, nI[i].vy, nI[i].vz, ux, uy, uz, wx, wy, wz);

//        for (int j = 1; j < n0.size(); ++j) { // go through the rest and generate projections
//            if (is_neighbour(nI[i], n0[j], daxial)) {
//                // append projections onto orthogonals u,(w)
//                vector<float> projadd(2);
//                projadd[0] = (n0[j].x-nI[i].x)*ux + (n0[j].y-nI[i].y)*uy + (n0[j].z-nI[i].z)*uz;
//                if (img3d)
//                    projadd[1] = (n0[j].x-nI[i].x)*wx + (n0[j].y-nI[i].y)*wy + (n0[j].z-nI[i].z)*wz;

//                proj.push_back(projadd);
//            }
//        }

//        if (proj.size()==1) continue; // no neighbours to refine, skip

//        // ms nI node across projections in the cross-secion plane
//        int iter = 0;
//        float d2 = 0;

//        proj_conv[0] = 0;
//        if (img3d)
//            proj_conv[1] = 0;

//        do {
//            ms(proj_conv, proj_new, proj, KRAD2, img3d);
//            d2 = pow(proj_new[0]-proj_conv[0],2) + (img3d)?pow(proj_new[1]-proj_conv[1],2):0;

//            proj_conv[0] = proj_new[0];
//            if(img3d)
//                proj_conv[1] = proj_new[1];

//            iter++;

//        }
//        while(iter<MAXITER && d2>EPSILON2);

//        float dx = proj_conv[0] * ux + (img3d)?(proj_conv[1] * wx):0;
//        float dy = proj_conv[0] * uy + (img3d)?(proj_conv[1] * wy):0;
//        float dz = proj_conv[0] * uz + (img3d)?(proj_conv[1] * wz):0;

//        cout << dx << " " << dy << " " << dz << endl;

//        // assign new x,y,z value to node in nI
//        nI[i].x += dx;
//        nI[i].y += dy;
//        nI[i].z += dz;

//        discr2 = dx*dx + dy*dy + dz*dz;

//        if (discr2>discr2_max) {
//            discr2_max = discr2;
//        }

//    }

//    return discr2_max;

//}
*/

bool is_bidirectional(vector<Node> n) {

    bool tt = true;

    for (int i = 1; i < n.size(); ++i) { // first index is dummy
        for (int j = 0; j < n[i].nbr.size(); ++j) {

            bool fnd = false;
            for (int k = 0; k < n[  n[i].nbr[j]  ].nbr.size(); ++k) {
                if (n[  n[i].nbr[j]  ].nbr[k] == i) {
                    fnd = true;
                    break;
                }
            }

//            if (!fnd) {
//                cout << "NOT FOUND" << endl;
//                cout << i << " -- ";
//                for (int i1 = 0; i1 < n[i].nbr.size(); ++i1) cout << " " << n[i].nbr[i1] << flush;
//                cout << endl;

//                cout << n[i].nbr[j] << " -- ";
//                for (int i1 = 0; i1 < n[  n[i].nbr[j]  ].nbr.size(); ++i1) cout << " " << n[  n[i].nbr[j]  ].nbr[i1] << flush;
//                cout << endl;
//            }

            tt = tt && fnd;

            if (!tt) { // loop stops first time it becomes false
                return tt;
            }
        }
    }

    return tt;
}

int get_undiscovered1(int * dist, int dist_len){

    for (int i = 0; i < dist_len; ++i) {
        if (dist[i]>0) {
            if (dist[i]==INT_MAX) {
                return i;
            }
        }
    }

    return -1;

}

//void bfs1(vector<Node> n, vector<Node>& tree) {
//    /*
//    1 Breadth-First-Search(Graph, root):
//    2
//    3     for each node n in Graph:
//    4         n.distance = INFINITY
//    5         n.parent = NIL
//    6
//    7     create empty queue Q
//    8
//    9     root.distance = 0
//    10     Q.enqueue(root)
//    11
//    12     while Q is not empty:
//    13
//    14         current = Q.dequeue()
//    15
//    16         for each node n that is adjacent to current:
//    17             if n.distance == INFINITY:
//    18                 n.distance = current.distance + 1
//    19                 n.parent = current
//    20                 Q.enqueue(n)
//    */
//    BfsQueue<int> q;
////    int * dist   = new int[n.size()]; // distances
////    int * nmap   = new int[n.size()]; // save indexing in output tree
////    int * parent = new int[n.size()]; // save parent index in current tree
//    int dist[n.size()];
//    int nmap[n.size()];
//    int parent[n.size()];
//    for (int i = 0; i < n.size(); ++i) {
//        dist[i] = INT_MAX;
//        nmap[i] = -1;
//        parent[i] = -1;
//    }
//    dist[0] = -1;
//    Node dummy_node;
//    tree.clear();
//    tree.push_back(dummy_node); // list index 0 is dummy node
//    int treecnt = 0;
//    int seed;

//    while ((seed = get_undiscovered1(dist, n.size()))>0) {

//        treecnt++;

//        dist[seed] = 0;
//        nmap[seed] = -1;
//        parent[seed] = -1;
//        q.enqueue(seed);

//        int nodesInTree = 0;

//        while (q.hasItems()) {

//            // dequeue(), take from FIFO structure
//            int curr = q.dequeue();

//            float x = n[curr].x;
//            float y = n[curr].y;
//            float z = n[curr].z;
//            float sig = n[curr].sig;

//            Node nd(x, y, z, sig, (treecnt+0)); // start from RED col
//            if (parent[curr] > 0) nd.nbr.push_back(nmap[parent[curr]]);
//            else nd.nbr.push_back(-1);
//            nmap[curr] = tree.size();
//            tree.push_back(nd);
//            nodesInTree++;
//            // for each node adjacent to current
//            for (int j = 0; j < n[curr].nbr.size(); j++) {

//                int adj = n[curr].nbr[j];

//                if (dist[adj] == INT_MAX) {
//                    dist[adj] = dist[curr] + 1;
//                    parent[adj] = curr;
//                    // enqueue(), add to FIFO structure
//                    q.enqueue(adj);
//                }
//            }
//        }
//    }
////    delete [] dist; dist = 0;
////    delete [] nmap; nmap = 0;
////    delete [] parent; parent = 0;
//}

void filter_node_density(vector<Node> nX, float wMin, vector<Node>& nY) {
    // calculate the density that at each node, remove outliers (nodes with density above wMin*mean(w))
    vector<float> nD(1); // corresonds to the nodelist, first one corresponds to the dummy node

    float nDmean; // iterative mean as density values per node are added
    int   cntnodes = 0;

    for (int i = 1; i < nX.size(); ++i) {
        float dcurr = 1;
        for (int j = 1; j < nX.size(); ++j) {
            if (j!=i) {
                float r2 = pow(3*nX[j].sig,2);
                float x2 = pow(nX[j].x-nX[i].x,2);
                if (x2<=r2) {
                    float y2 = pow(nX[j].y-nX[i].y,2);
                    if (x2+y2<=r2) {
                        float z2 = pow(nX[j].z-nX[i].z,2);
                        if (x2+y2+z2<=r2) {
                            dcurr += exp(-(x2+y2+z2)/(2*pow(nX[j].sig,2)));
                        }
                    }
                }
            }
        }

        nD.push_back(dcurr);
        cntnodes++;
        nDmean = (cntnodes==1)?dcurr:(((cntnodes-1.0)/(float) cntnodes)*nDmean+(1.0/cntnodes)*dcurr);

    }

//    int XtoY[nX.size()]; // map
    vector<int> XtoY(nX.size());
    XtoY[0] = 0;
    nY.clear();
    Node nYroot;
    nY.push_back(nYroot);
    int cntremoved = 0;
    for (int i = 1; i < nX.size(); ++i) {
        if (true || nD[i]>wMin*nDmean) {
            Node nYadd(nX[i]);
            nYadd.sig = nX[i].corr; // nD[i]/nDmean;
            XtoY[i] = nY.size();
            nY.push_back(nYadd);
        }
        else {
            XtoY[i] = -1;
            cntremoved++;
        }
    }

    // express neighbour indexes in terms of nY list indexes
    for (int i = 1; i < nY.size(); ++i) {
        for (int j = nY[i].nbr.size()-1; j >= 0; --j) {
            int nidx = nY[i].nbr[j];
            if (XtoY[nidx]==-1) {
                nY[i].nbr.erase(nY[i].nbr.begin() + j); // erase jth element from nY[i].nbr list
            }
            else {
                nY[i].nbr[j] = XtoY[nidx]; // replace jth with nY index
            }
        }
    }
}

void filter_internode_dist(vector<Node>& nX, float dth) {
    // will remove all the linkages (assumes bidirectional) that are above an euclidean distance threshold dth, removes the links only
    // links tht involve at least one soma node are not removed
    int count = 0;
    for (long i = 1; i < nX.size(); ++i) {
        if (nX[i].type!=Node::SOMA) { // can't be soma for the link to be removed
            for (int j = nX[i].nbr.size()-1; j >= 0; --j) {
                long i1 = nX[i].nbr[j];
                if (nX[i1].type!=Node::SOMA) {
                    float dist = sqrt(pow(nX[i1].x-nX[i].x,2)+pow(nX[i1].y-nX[i].y,2)+pow(nX[i1].z-nX[i].z,2));
                    if (dist>dth) {
                        count++;
                        nX[i].nbr.erase(nX[i].nbr.begin() + j);// remove element at index j, i->i1 link
                        int i2 = find(nX[i1].nbr.begin(), nX[i1].nbr.end(), i) - nX[i1].nbr.begin(); // i value at nY[i1].nbr[i2]
                        if (i2<nX[i1].nbr.size()) nX[i1].nbr.erase(nX[i1].nbr.begin()+i2); // remove i value at idx i2, i1->i link
                    }
                }
            }
        }
    }

    cout << count << " links > " << dth << endl;
}

void check_nbr(vector<Node>& nX) {
    // - ensure linkings are bidirectional
    for (long i = 1; i < nX.size(); ++i) {
        // remove double neighbourhoods from the neighbour list
        sort(nX[i].nbr.begin(), nX[i].nbr.end());
        nX[i].nbr.erase(unique(nX[i].nbr.begin(), nX[i].nbr.end()), nX[i].nbr.end());

        // remove self linkages
        int pos = find(nX[i].nbr.begin(), nX[i].nbr.end(), i) - nX[i].nbr.begin();
        if (pos>=0 && pos<nX[i].nbr.size())
            nX[i].nbr.erase(nX[i].nbr.begin()+pos); // remove at pos
    }

    // ensure linkings are bidirectional, add if not
    for (long i = 1; i < nX.size(); ++i) { // first index is dummy
        for (int j = 0; j < nX[i].nbr.size(); ++j) {
            bool fnd = false;
            for (int k = 0; k < nX[  nX[i].nbr[j]  ].nbr.size(); ++k) {
                if (nX[  nX[i].nbr[j]  ].nbr[k] == i) {
                    fnd = true;
                    break;
                }
            }

            if (!fnd) {
                // enforce link
                nX[ nX[i].nbr[j] ].nbr.push_back(i);
                cout << "enforce bidirectional link: " << nX[i].nbr[j] << " -- " << i << endl;
            }
        }
    }

}

void group1(vector<Node> nX, vector<Node>& nY, float rad) { // sphere grouping

    nX[0].corr = FLT_MAX; // so that the dummy node gets index 0 again, larges correlation
    vector<long> indices(nX.size());
    for (long i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareIndicesByNodeCorrVal(&nX));

    vector<long> X2Y(nX.size(), -1);
    X2Y[0] = 0; // first one is with max. correlation

    nY.clear();
    Node nY0(nX[0]);
    nY.push_back(nY0);

    for (long i = 1; i < nX.size(); ++i) { // add soma nodes as independent groups at the beginning
        if (nX[i].type==Node::SOMA) {
            X2Y[i] = nY.size();
            Node nYi(nX[i]);
            nYi.type = Node::SOMA;
            nY.push_back(nYi);
        }
    }

    for (long i = 1; i < indices.size(); ++i) { // add the rest of the nodes

        long ci = indices[i];

        if (X2Y[ci]!=-1) continue; // skip if it was added to a group already

        X2Y[ci] = nY.size();
        Node nYi(nX[ci]); // nX[ci] is the grabbed
        float grp_size = 1;

        float r2 = rad*rad; // sig2rad * nX[ci].sig;
        float d2;
        for (long j = 1; j < nX.size(); ++j) { // check the rest that was not groupped
            if (j!=ci && X2Y[j]==-1) {
                d2 = pow(nX[j].x-nX[ci].x,2);
                if (d2<=r2) {
                    d2 += pow(nX[j].y-nX[ci].y,2);
                    if (d2<=r2) {
                        d2 += pow(nX[j].z-nX[ci].z,2);
                        if (d2<=r2) {

                            X2Y[j]=nY.size();

                            for (int k = 0; k < nX[j].nbr.size(); ++k)
                                nYi.nbr.push_back(nX[j].nbr[k]);  // append the neighbours of the group members

                            // update local average with x,y,z,sig elements from nX[j]
                            grp_size++;
                            float a = (grp_size-1)/grp_size;
                            float b = (1.0/grp_size);
                            nYi.x   = a * nYi.x     + b * nX[j].x;
                            nYi.y   = a * nYi.y     + b * nX[j].y;
                            nYi.z   = a * nYi.z     + b * nX[j].z;
                            nYi.sig = a * nYi.sig   + b * nX[j].sig;
                            nYi.corr= a * nYi.corr  + b * nX[j].corr;
                        }
                    }
                }
            }
        }

        nYi.type = Node::AXON; // enforce type
        nY.push_back(nYi);
    }

    for (int i = 1; i < nY.size(); ++i) {
        for (int j = 0; j < nY[i].nbr.size(); ++j) {
            nY[i].nbr[j] = X2Y[ nY[i].nbr[j] ];
        }
    }

    check_nbr(nY); // remove doubles and self-linkages after grouping

}

// cylinder grouping (experimental)
void group0(vector<Node> nX, vector<Node>& nY, float d_axial=2, float sig2rad=1) {

    // grouping starts from the nodes from nX[] with highest correlation
    nX[0].corr = FLT_MAX;               // to ensure it's at index 0 after sorting by correlation
    vector<int> indices(nX.size());     // indices will guide the grouping sequence
    for (int i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareIndicesByNodeCorrVal(&nX)); // sort by correlation

    // index mapping
    vector<int> X2Y(nX.size(), -1);     // filled up during grouping - -1 means that it was not grouped yet, once it gets grouped it becomes assigned with the index in nY[] vector where the group is added X2Y[5]=1 means nY[5] is among the nodes that contributed to grouped element nY[1]
    X2Y[0] = 0;                         // dummy node mapping

    nY.clear();
    Node nY_dummy(nX[0]);
    nY.push_back(nY_dummy);             // initialize nY[0] with nX[0]

    vector<int> Nci; // indexes (nX vector) of those nodes around nX[ci], in sphere with sig2rad*nX[ci].sig radius, circular neighbourhood
    Pxyz    mu(0,0,0); // local direction using elements from N1, covariance of the local 3d coordinates
    double  cov[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
    double  vec[3][3];
    double  eig[3];

    for (int i = 1; i < indices.size(); ++i) {
        // indices[i] denotes the index in nX[] with highest correlation
        int ci = indices[i]; // group nX[ci] Node nX[ci].x, nX[ci].y, nX[ci].z...

        if (X2Y[ci] != -1) continue; // it has been assigned to the group nY[X2Y[ci]]

        // gather the nodes (their indexes) from the spherical neighbourhood used to calculate local direction
        Nci.clear(); // reset index list, nX indexing
        Nci.push_back(ci);
        float r2 = sig2rad * nX[ci].sig;
        for (int j = 1; j < nX.size(); ++j) {
            if (j!=ci) { // it's added
                float x2 = pow(nX[j].x-nX[ci].x,2);
                if (x2<=r2) {
                    float y2 = pow(nX[j].y-nX[ci].y,2);
                    if (x2+y2<=r2) {
                        float z2 = pow(nX[j].z-nX[ci].z   ,2);
                        if (x2+y2+z2<=r2) {
                            Nci.push_back(j);
                        }
                    }
                }
            }
        }

        if (Nci.size()==1) {
            X2Y[Nci[0]] = nY.size();
            Node nYi(nX[Nci[0]]);
            nY.push_back(nYi);
        }
        else { // estimate local orientation using 3d coordinates from nodes in N1 set of indexes

            // mu
            mu.x = 0; mu.y = 0; mu.z = 0;
            for (int j = 0; j < Nci.size(); ++j) {
                mu.x += nX[Nci[j]].x;
                mu.y += nX[Nci[j]].y;
                mu.z += nX[Nci[j]].z;
            }
            mu.x /= Nci.size(); mu.y /= Nci.size(); mu.z /= Nci.size();

            // cov
            cov[0][0] = 0; cov[0][1] = 0; cov[0][2] = 0;
            cov[1][0] = 0; cov[1][1] = 0; cov[1][2] = 0;
            cov[2][0] = 0; cov[2][1] = 0; cov[2][2] = 0;

            for (int j = 0; j < Nci.size(); ++j) {
                float dx = nX[Nci[j]].x-mu.x;
                float dy = nX[Nci[j]].y-mu.y;
                float dz = nX[Nci[j]].z-mu.z;
                cov[0][0] += dx*dx; cov[0][1] += dx*dy; cov[0][2] += dx*dz;
                cov[1][0] += dy*dx; cov[1][1] += dy*dy; cov[1][2] += dy*dz;
                cov[2][0] += dz*dx; cov[2][1] += dz*dy; cov[2][2] += dz*dz;
            }

            cov[0][0] /= Nci.size(); cov[0][1] /= Nci.size(); cov[0][2] /= Nci.size();
            cov[1][0] /= Nci.size(); cov[1][1] /= Nci.size(); cov[1][2] /= Nci.size();
            cov[2][0] /= Nci.size(); cov[2][1] /= Nci.size(); cov[2][2] /= Nci.size();

            // ...
            Frangi::eigen_decomposition_static(cov, vec, eig);

            // use vec[0][2], vec[1][2], vec[2][2] to extract the orthogonal cross-section elements to group them
            X2Y[Nci[0]] = nY.size();
            Node nYi(nX[Nci[0]]); // (0,0,0,     vec[0][2],vec[1][2],vec[2][2],      0, 0, Node::RED);
            float grp_size = 1;

            for (int j = 1; j < Nci.size(); ++j) { // new group: find ungroupped that were in the cross section plane centered around current node
                if (is_cross_section(nX[Nci[0]].x, nX[Nci[0]].y, nX[Nci[0]].z, vec[0][2], vec[1][2], vec[2][2],
                                 nX[Nci[j]].x, nX[Nci[j]].y, nX[Nci[j]].z, d_axial) && X2Y[ Nci[j] ]==-1) {

                    X2Y[Nci[j]] = nY.size();

                    grp_size++;
                    float a = (grp_size-1)/grp_size;
                    float b =  1.0/grp_size;
                    nYi.x       += a * nYi.x    + b * nX[Nci[j]].x;
                    nYi.y       += a * nYi.y    + b * nX[Nci[j]].y;
                    nYi.z       += a * nYi.z    + b * nX[Nci[j]].z;
                    nYi.sig     += a * nYi.sig  + b * nX[Nci[j]].sig;
                    nYi.corr    += a * nYi.corr + b * nX[Nci[j]].corr;

                    for (int k = 0; k < nX[Nci[j]].nbr.size(); ++k)
                        nYi.nbr.push_back(nX[Nci[j]].nbr[k]); // new node will accumulate neighbours from all of the members

                }
            }

            nY.push_back(nYi);

        }
    }

    // remap
    for (int i = 1; i < nY.size(); ++i) {
        for (int j = 0; j < nY[i].nbr.size(); ++j) {
            nY[i].nbr[j] = X2Y[ nY[i].nbr[j] ];
        }
    }

}

void get_node_density(  vector<Node> nX, float sig2rad, vector<float>& d) {

    d.clear();

//    float dmin =  FLT_MAX;
//    float dmax = -FLT_MAX;

    for (int i = 1; i < nX.size(); ++i) {
        float w = 1;
        for (int j = 1; j < nX.size(); ++j) {
            if (j!=i) {
                float r2 = pow(sig2rad*nX[i].sig,2);
                float x2 = pow(nX[j].x-nX[i].x,2);
                if (x2<=r2) {
                    float y2 = pow(nX[j].y-nX[i].y,2);
                    if (x2+y2<=r2) {
                        float z2 = pow(nX[j].z-nX[i].z,2);
                        if (x2+y2+z2<=r2) {
                            w += exp(-(x2+y2+z2)/(2*pow(nX[i].sig,2)));
                        }
                    }
                }
            }
        }
        d.push_back(w);
//        if (w>dmax) dmax = w;
//        if (w<dmin) dmin = w;
    }
}

void get_link_lengths(  vector<Node> nX, vector<float>& l) {

    // all internode euclidean lengths exept the ones that involve soma nodes (at least one soma node)

    vector< vector<bool> > chk(nX.size());
    for (long i = 0; i < nX.size(); ++i) {
        vector<bool> chk1(nX[i].nbr.size(), false);
//        for (int j = 0; j < nX[i].nbr.size(); ++j) {
//            if (nX[i].type==Node::SOMA || nX[nX[i].nbr[j]].type==Node::SOMA) {
//                chk1[j] = true;
//            }
//        }
        chk[i] = chk1;
    }

    l.clear();

    for (long i = 1; i < nX.size(); ++i) {

        for (long j = 0; j < nX[i].nbr.size(); ++j) {

            if (!chk[i][j]) {
                    int nidx = nX[i].nbr[j];
                    int pos = find(nX[nidx].nbr.begin(), nX[nidx].nbr.end(), i) - nX[nidx].nbr.begin();
                    if (pos<nX[nidx].nbr.size()) {
                        chk[i][j]       = true;
                        chk[nidx][pos]  = true;
                        l.push_back(sqrt(pow(nX[i].x-nX[nidx].x,2) + pow(nX[i].y-nX[nidx].y,2) + pow(nX[i].z-nX[nidx].z,2)));
                    }

            }

        }
    }

}

void get_node_corr(     vector<Node> nX, vector<float>& c) {
    c.clear();
    for (long i = 1; i < nX.size(); ++i) {
        if (nX[i].type!=Node::SOMA) {
        c.push_back(nX[i].corr);
        }
    }
}

void get_internode_dist(vector<Node> nX, vector< vector<float> > nL) {

    nL.clear();
    vector<float> t;
    nL.push_back(t);

    for (int i = 1; i < nX.size(); ++i) {
        vector<float> tt(nX[i].nbr.size());
        for (int j = 0; j < nX[i].nbr.size(); ++j) {
            int i1 = nX[i].nbr[j];
            tt[j] = sqrt(pow(nX[i].x-nX[i1].x,2)+pow(nX[i].y-nX[i1].y,2)+pow(nX[i].z-nX[i1].z,2));
        }
        nL.push_back(tt);
    }

}

void export_directionality(string swcpath, unsigned char* J, int w, int h, int l, unsigned char Jth, unsigned char* Vx, unsigned char* Vy, unsigned char* Vz) {
    ofstream f;
    f.open(swcpath.c_str(), ofstream::out | ofstream::trunc);

    int count = 1;
    for (long i = 0; i < (w*h*l); ++i) { // go through all the voxels
        if (J[i]>Jth) { // show vectors for those with tubularity above Jth [0-255]
            f<<count<<" "<<Node::OCRE_LIGHT<<" "<<i%w<<" "<<(i/w-(i/(w*h))*h)<<" "<<(i/(w*h))<<" "<<0.1<<" "<<(-1)<<endl;
            count++;
            float v = 10;//(J[i]*10.0/255);
            float ux = (((float)Vx[i]/255)*2)-1;
            float uy = (((float)Vy[i]/255)*2)-1;
            float uz = (((float)Vz[i]/255)*2)-1;
            f<<count<<" "<<Node::OCRE_LIGHT<<" "<<(i%w+v*ux)<<" "<<(i/w-(i/(w*h))*h+v*uy)<<" "<<(i/(w*h)+v*uz)<<" "<<0.1<<" "<<(count-1)<<endl;
            count++;
        }
    }

    f.close();

    cout<<"exported: "<<swcpath<<endl;
}

void parse_csv_string(string s, vector<float>& sigs) {
    std::stringstream ss(s);
    float i;
    sigs.clear();
    while (ss >> i) {
        sigs.push_back(i);
        if (ss.peek() == ',') // || ss.peek() == ' '
                ss.ignore();
    }

    std::sort (sigs.begin(), sigs.begin()+sigs.size());

}

void soma_extraction1(unsigned char* E8, unsigned char E8th, int N, int M, int P, int* smap, vector<Node>& n0) {
    // new version where soma is modelled with a sphere (x,y,z,r neuron node) and a int* smap soma map assigned to it
    unsigned char* E8bin = new unsigned char[N*M*P]; // binarized soma
    for (long i = 0; i < N*M*P; ++i) {
        E8bin[i] = (E8[i]>E8th)?UCHAR_MAX:0;
        smap[i] = 0;
    }

    vector<float> xc, yc, zc, rc;
    conn3d(E8bin, N, M, P, smap, INT_MAX, true, 0, 1, xc, yc, zc, rc);
    delete[] E8bin; E8bin = 0;

    for (int i = 0; i < xc.size(); ++i) {
        Node nd(xc[i], yc[i], zc[i], rc[i], Node::SOMA);
        n0.push_back(nd);
    }
}

void soma_extraction(unsigned char* E8, unsigned char E8th, int N, int M, int P, float Rgrp, float zdist, int* smap, vector<Node>& n0) {

    long ss = 0; // count soma voxels
    for (long i = 0; i < N*M*P; ++i) {
        smap[i] = 0; // reset soma node map smap[]
        if (E8[i]>E8th)
            ss++;
    }

    cout << "soma occupies " << ss/(float)(N*M*P) << "% vox. vol." << endl;

    long* si = new long[ss]; // global indexes of each soma voxel
    long* sk = new long[ss]; // local sv[] array indexes - used for sorting sv[] array
    unsigned char* sv = new unsigned char[ss]; // greylevel of each soma voxel

    ss = 0;
    for (long i = 0; i < N*M*P; ++i) {
        if (E8[i]>E8th) {
            si[ss] = i;
            sv[ss] = E8[i];
            sk[ss] = ss;
            ss++;
        }
        else
            E8[i] = 0;
    }

    clock_t t1, t2;
    t1 = clock();
    quicksort(sv, sk, 0, ss-1); // sort soma voxels
    t2 = clock();
    cout << " quicksort " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

    // instantiate trace node with soma nodes covering the volume that was assigned to soma after erosion+gauss+threshold+dilatation
    // E8 contains byte8 image where >0 pixels/voxels denote soma
    // add soma region to the nodelist before going further with the trace
    // establish soma nodelist index map (smap) where each 255 voxel is assigned with nodelist index
    // map is used to stop the traces and link them with nodes that belong to the soma
    vector<offxyz> Sxyz;
    Tracker::sphereXYZ(Rgrp, zdist, Sxyz);

    int scnt = 0; // soma node counter

    int x0, y0, z0, x1, y1, z1;
    long i0, i1;

    for (long i = 0; i < ss; ++i) { // go through descending soma voxels

        i0 = si[sk[i]]; // pick global index of the highest sv[i]

        if (smap[i0]==0) { // if current soma voxel location is not mapped yet

            scnt++; // node counter increment
            float xavg = 0;
            float yavg = 0;
            float zavg = 0;
            int   cnt  = 0;

            x0 = i0%N;
            y0 = i0/N-(i0/(N*M))*M;
            z0 =       i0/(N*M);

            // loop spherical neighbourhood, if E8[] value was >E8th then add it to the cluster, each added is iteratively averaged,
            for (int k = 0; k < Sxyz.size(); ++k) {

                x1 = x0 + Sxyz[k].x;
                y1 = y0 + Sxyz[k].y;
                z1 = z0 + Sxyz[k].z;

                if (x1>=0 && x1<N && y1>=0 && y1<M && z1>=0 && z1<P) {

                    i1 = z1*N*M+y1*N+x1;

                    if (E8[i1]>E8th) {
                        cnt++;
                        float A = (cnt-1)/(float)cnt;
                        float B = 1.0/cnt;
                        xavg=A*xavg+B*x1;
                        yavg=A*yavg+B*y1;
                        zavg=A*zavg+B*z1;

                        smap[i1] = scnt;

                    }
                }
            }

            float ravg = 1 + ((Rgrp-1)/(Sxyz.size()-1))*(cnt-1);

            Node nd(xavg, yavg, zavg, Rgrp, Node::SOMA);
            n0.push_back(nd);

        }
    }

    delete [] sv; sv = 0;
    delete [] sk; sk = 0;


//    cout << "ADDITIONAL check..." << endl;
//    for (long i = 0; i < N*M*P; ++i) {
//        if (E8[i]>E8th && smap[i]==0) {
//            cout << "PROBLEM" << endl;
//        }
//        if (E8[i]<=E8th && smap[i]>0){
//            cout << "PROBLEM" <<  endl;
//        }
//    }
//    cout << "done" << endl;

    // soma node linking, link those soma node tags from smap[] that were 4-connected
    for (long i = 0; i < ss; ++i) { // go through the soma voxels again (those that got E8[i]>E8th)

        i0 = si[i];
        x0 = i0%N;
        z0 = i0/(N*M);
        y0 = i0/N-z0*M;

//        if (smap[i0]==0) cout << "!!! smap[i0]==0 in soma" << endl;
//        if (E8[i0]==0)   cout << "!!! E8[i0]==0 in soma" << endl;

        x1 = x0-1;
        if (x1>=0 && x1<N) {
            i1 = z0*N*M+y0*N+x1;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

        x1 = x0+1;
        if (x1>=0 && x1<N) {
            i1 = z0*N*M+y0*N+x1;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

        y1 = y0-1;
        if (y1>=0 && y1<M) {
            i1 = z0*N*M+y1*N+x0;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

        y1 = y0+1;
        if (y1>=0 && y1<M) {
            i1 = z0*N*M+y1*N+x0;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

        z1 = z0-1;
        if (z1>=0 && z1<P) {
            i1 = z1*N*M+y0*N+x0;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

        z1 = z0+1;
        if (z1>=0 && z1<P) {
            i1 = z1*N*M+y0*N+x0;
            if (smap[i1]>0 && smap[i1]!=smap[i0]) n0[smap[i0]].nbr.push_back(smap[i1]);
        }

    }

    delete [] si; si = 0; // it was used in linking too

    cout << "done linking" << endl;

    // remove duplicate neighbourhood node links from the neighbour list
    for (int ni = 1; ni < n0.size(); ++ni) {
        sort(n0[ni].nbr.begin(), n0[ni].nbr.end());
        n0[ni].nbr.erase(unique(n0[ni].nbr.begin(), n0[ni].nbr.end()), n0[ni].nbr.end());
    }

    // remove double linkages
    for (int i = 1; i < n0.size(); ++i) {

    }

    // check if links are bidirectional... (remove later!)
    cout << "\n0 bidirectional? " << is_bidirectional(n0) << endl;

}

void reconstruct(vector<Node> n0, QString prefix, QString suffix) {

    if (saveMidres) {
        save_nodelist(n0,                   prefix+"_n0_"+suffix+".swc");
        save_nodelist(compute_trees(n0),    prefix+"_n0tree_"+suffix+".swc");
        vector<float> n0len;
        get_link_lengths(n0, n0len);
        save_vector(n0len,                  prefix+"_n0len_"+suffix+".log");
        vector<float> n0corr;
        get_node_corr(n0, n0corr);
        save_vector(n0corr,                 prefix+"_n0corr_"+suffix+".log");
    }

    interpolate_nodelist(n0, TRACE_RSMPL);

    if (saveMidres) {
        save_nodelist(n0,                   prefix+"_n0res_"+suffix+".swc");
    }

    vector<Node> n1;
//    refine_blurring(n0, n1, SIG2RADIUS, REFINE_ITER, EPSILON2);
    non_blurring(   n0, n1, SIG2RADIUS, REFINE_ITER, EPSILON2);

    if (saveMidres) {
        save_nodelist(n1,                   prefix+"_n1_"+suffix+".swc");

        vector<float> n1len;
        get_link_lengths(n1, n1len);
        save_vector(n1len,                  prefix+"_n1len_"+suffix+".log");

    }

    vector<Node> n2;
    group1(n1, n2, GROUP_RADIUS);
    n1.clear();

    if (saveMidres) {
        save_nodelist(n2,                   prefix+"_n2_"+suffix+".swc");
    }

    vector<Node> n2tree = compute_trees(n2); // assign each connected tree a new colour and keep the somas with Node::SOMA
    n2.clear();

    if (saveMidres) {
        save_nodelist(n2tree,               prefix+"_n2tree_"+suffix+".swc");
    }
    if (ENFORCE_SINGLE_TREE) { // true ||

        vector<Node> n3tree = extract_largest_tree(n2tree);

//        if (saveMidres) {
//            save_nodelist(n3tree,                prefix+"_n3tree1_"+suffix+".swc");
//        }

        interpolate_treelist(n3tree, 1.0, Node::AXON);

        save_nodelist(n3tree, prefix + "_Advantra1"+suffix+".swc", -1, 1, NAME, COMMENT);

    }

    n2tree.clear();

//    vector<Node> n3tree = (ENFORCE_SINGLE_TREE)?extract_largest_tree(n2tree):extract_trees(n2tree, TREE_SIZE_MIN);
//    n2tree.clear();

//    if (saveMidres) {
//        save_nodelist(n3tree,                prefix+"_n3tree_"+suffix+".swc");
//    }

//    vector<Node> n4tree = n3tree;
//    vector<Node> n4tree = remove_tails(n3tree, TAIL_SIZE_MIN); // expell tails (end-junction) with less than TAIL_SIZE_MIN nodes
//    n3tree.clear();
}

void reconstruction_func(V3DPluginCallback2 &callback, QWidget *parent, input_PARA &PARA, bool bmenu)
{

    unsigned char* data1d = 0;
    V3DLONG N,M,P,sc,c;
    V3DLONG in_sz[4];
    if(bmenu)
    {
        v3dhandle curwin = callback.currentImageWindow();
        if (!curwin)
        {
            QMessageBox::information(0, "", "You don't have any image open in the main window.");
            return;
        }

        Image4DSimple* p4DImage = callback.getImage(curwin);

        if (!p4DImage)
        {
            QMessageBox::information(0, "", "The image pointer is invalid. Ensure your data is valid and try again!");
            return;
        }


        data1d = p4DImage->getRawData();
        N = p4DImage->getXDim();
        M = p4DImage->getYDim();
        P = p4DImage->getZDim();
        sc = p4DImage->getCDim();

        bool ok1;

        if(sc==1)
        {
            c=1;
            ok1=true;
        }
        else
        {
            c = QInputDialog::getInteger(parent, "Channel",
                                             "Enter channel NO:",
                                             1, 1, sc, 1, &ok1);
        }

        if(!ok1)
            return;

        in_sz[0] = N;
        in_sz[1] = M;
        in_sz[2] = P;
        in_sz[3] = sc;


        PARA.inimg_file = p4DImage->getFileName();
    }
    else
    {
        int datatype = 0;
        if (!simple_loadimage_wrapper(callback,PARA.inimg_file.toStdString().c_str(), data1d, in_sz, datatype))
        {
            fprintf (stderr, "Error happens in reading the subject file [%s]. Exit. \n",PARA.inimg_file.toStdString().c_str());
            return;
        }
        if(PARA.channel < 1 || PARA.channel > in_sz[3])
        {
            fprintf (stderr, "Invalid channel number. \n");
            return;
        }
        N = in_sz[0];
        M = in_sz[1];
        P = in_sz[2];
        sc = in_sz[3];
        c = PARA.channel;
    }

    //// THIS IS WHERE THE DEVELOPERS SHOULD ADD THEIR OWN NEURON TRACING CODE

    cout<<"-------------  ADVANTRA  -------------"  <<endl;
    cout<<"neuritesigmas = " <<PARA.neuritesigmas              <<endl;
    cout<<"somaradius = "    <<PARA.somaradius              <<endl;
    cout<<"tolerance = "     <<PARA.tolerance           <<endl;
    cout<<"znccth = "    <<PARA.znccth              <<endl;
    cout<<"kappa = "     <<PARA.kappa               <<endl;
    cout<<"step = "      <<PARA.step                <<endl;
    cout<<"ni = "        <<PARA.ni                  <<endl;
    cout<<"np = "        <<PARA.np                  <<endl;
    cout<<"zdist = "     <<PARA.zdist               <<endl;
    cout<<"nodepp="      <<PARA.nodepervol          <<endl;
    cout<<"vol="         <<PARA.vol                 <<endl;
    cout<<"-------------------------------------"   <<endl;

//    if (1) {cout<<"go out"<< endl;return;}

    NAME = "Advantra";
    COMMENT =
            "email: miro@braincadet.com\n#params:\n#channel="+QString("%1").arg(PARA.channel)+
            "\n#neuritesigmas="+QString::fromStdString(PARA.neuritesigmas)+
            "\n#somaradius="+QString("%1").arg(PARA.somaradius)+
            "\n#tolerance="+QString("%1").arg(PARA.tolerance)+
            "\n#znccth="+QString("%1").arg(PARA.znccth)+
            "\n#kappa="+QString("%1").arg(PARA.kappa)+
            "\n#step="+QString("%1").arg(PARA.step)+
            "\n#ni="+QString("%1").arg(PARA.ni)+
            "\n#np="+QString("%1").arg(PARA.np)+
            "\n#zdist="+QString("%1").arg(PARA.zdist)+
            "\n#nodepervol="+QString("%1").arg(PARA.nodepervol)+
            "\n#vol="+QString("%1").arg(PARA.vol)+
            "\n#------------------------"+
            "\n#Kc="+QString("%1").arg(Kc)+
            "\n#neff_ratio="+QString("%1").arg(neff_ratio)+
            "\n#frangi_alfa="+QString("%1").arg(frangi_alfa)+
            "\n#frangi_beta="+QString("%1").arg(frangi_beta)+
            "\n#frangi_C="+QString("%1").arg(frangi_C)+
            "\n#frangi_betaone="+QString("%1").arg(frangi_betaone)+
            "\n#frangi_betatwo="+QString("%1").arg(frangi_betatwo)+
            "\n#MAX_TRACE_COUNT="+QString("%1").arg(MAX_TRACE_COUNT)+
            "\n#EPSILON2="+QString("%1").arg(EPSILON2)+
            "\n#REFINE_ITER="+QString("%1").arg(REFINE_ITER)+
            "\n#SIG2RADIUS="+QString("%1").arg(SIG2RADIUS)+
            "\n#TRACE_RSMPL="+QString("%1").arg(TRACE_RSMPL)+
            "\n#GROUP_RADIUS="+QString("%1").arg(GROUP_RADIUS)+
            "\n#ENFORCE_SINGLE_TREE="+QString("%1").arg(ENFORCE_SINGLE_TREE)+
            "\n#TREE_SIZE_MIN="+QString("%1").arg(TREE_SIZE_MIN)+
            "\n#TAIL_SIZE_MIN="+QString("%1").arg(TAIL_SIZE_MIN);

    if (0) {

        cout << "findMaxima()" << endl;
        vector<Pxyz> localMaxima;
        SeedExtractor::findMaxima(data1d, N, M, P, PARA.tolerance, localMaxima); // int w, int h, int l, double tolerance, vector<Pxy>& xyVector

        // export to swc file
        string savepath = PARA.inimg_file.toStdString() + "_findMaximaAdvantra_" + QString("%1").arg(PARA.tolerance).toStdString() + ".swc";
        ofstream f;
        f.open(savepath.c_str(), ofstream::out | ofstream::trunc);
        long cnt = 1;
        for (long i = 0; i < localMaxima.size(); ++i)
            f << fixed << setprecision(0) << (cnt++) << " " << 13 <<" "<< localMaxima[i].x <<" "<< localMaxima[i].y <<" "<< localMaxima[i].z <<" .1 -1\n";

        f.close();
        cout << localMaxima.size() << " local maxima" << endl;
        cout << savepath << endl;
        cout << "\nDONE" << endl;
        return;
    }

    long size = N * M * P; // N : width, M : height, P : nr. layers

    vector<float> sigs; // read comma delimited sigmas from PARA.sigmas
    parse_csv_string(PARA.neuritesigmas, sigs); // extract them into a sorted vector list sigs[0]<=sigs[1]<=sigs[2]...

    clock_t t1, t2;

    if (0) {
        float sig = 6;
        cout << "G["<< sig << "]... "<< flush;

        t1 = clock();
        float* G = new float[size];
        t2 = clock();
        cout << "alloc. " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

        t1 = clock();
        Frangi::imgaussian(data1d, N, M, P, sig, PARA.zdist, G);
        t2 = clock();
        cout << " proc. " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;
        delete [] G; G = 0;
        cout << "DONE" << endl;
        return;
    }

    if (0) {
        // benchmark zncc() vs. zncc_v1() vs. zncc_v2()
        Tracker t_test(sigs, PARA.step, PARA.np, PARA.ni, PARA.kappa, P==1, PARA.znccth, Kc, neff_ratio, PARA.zdist, PARA.nodepervol);
        srand (time(NULL));
        streamsize ss = cout.precision();
        float tA_avg = 0, tB_avg = 0, tC_avg = 0, t_cnt = 0, diff_AC_max = -FLT_MAX, diff_AB_max = -FLT_MAX, tA, tB, tC;
        for (long i = 0; i < N*M*P; ++i) { // random locations with random directions
            if ((int)data1d[i]>200) {

                int x1  = i%N;
                int z1  = i/(N*M);
                int y1  = i/N-z1*M;
                float vx1 = ((double) rand() / (RAND_MAX));
                float vy1 = ((double) rand() / (RAND_MAX));
                float vz1 = ((double) rand() / (RAND_MAX));
                float vn = sqrt(pow(vx1,2)+pow(vy1,2)+pow(vz1,2));
                vx1 /= vn;
                vy1 /= vn;
                vz1 /= vn;

    //            cout << fixed << setprecision(2) << "[" << x1 << "," << y1 << "," << z1 << "] ("<< vx1 << "," << vy1 << "," << vz1 <<"):" << endl;

                float dummy;

                t1 = clock();
                float znccA = t_test.zncc(x1, y1, z1, vx1, vy1, vz1, 0, data1d, N, M, P, dummy);
                t2 = clock();
                tA = (t2-t1)/(double)CLOCKS_PER_SEC;
                tA_avg += tA;
    //            cout << fixed << setprecision(20) <<"znccA=" << znccA << " [" << tA << " sec.]" << endl;

                t1 = clock();
                float znccB = t_test.znccAAA(x1, y1, z1, vx1, vy1, vz1, data1d, N, M, P, dummy);
                t2 = clock();
                tB = (t2-t1)/(double)CLOCKS_PER_SEC;
                tB_avg += tB;
    //            cout << fixed << setprecision(20) << "znccB=" << znccB << " [" << tB << " sec.]" << endl;

                t1 = clock();
                float znccC = t_test.znccBBB(x1, y1, z1, vx1, vy1, vz1, data1d, N, M, P, dummy);
                t2 = clock();
                tC = (t2-t1)/(double)CLOCKS_PER_SEC;
                tC_avg += tC;
    //            cout << fixed << setprecision(20) <<"znccC=" << znccC << " [" << tC << " sec.]" << endl;

                if (abs(znccA-znccC)>diff_AC_max) diff_AC_max = abs(znccA-znccC);
                if (abs(znccA-znccB)>diff_AB_max) diff_AB_max = abs(znccA-znccB);

                t_cnt += 1;

            }
        }

        cout << fixed << setprecision(0)  << t_cnt << " samples:" << endl;
        cout << fixed << setprecision(15) << "average exec time: A = " << (tA_avg/t_cnt) << " B = " << (tB_avg/t_cnt) << " C = " << (tC_avg/t_cnt) << endl;
        cout << fixed << setprecision(3) << "x" << ((tA_avg/t_cnt)/(tC_avg/t_cnt)) << endl;
        cout << "precision difference:" << endl;
        cout << fixed << setprecision(15) << "A-B corr diff_max=" << diff_AB_max << endl;
        cout << fixed << setprecision(15) << "A-C corr diff_max=" << diff_AC_max << endl;
        cout.precision (ss); // retrieve earlier setting
    }

    vector<Node> n0; // list of node used for the initial trace
    Node n00;
//    .push_back(0);
    n0.push_back(n00);

    int* smap = new int[size]; // soma map

    /*
     * SOMA EXTR.
     */
    if (PARA.somaradius>0) {
        // extract soma
        unsigned char*   E8   = new unsigned char[size];

        cout << "imerode("<< PARA.somaradius <<") " << flush;
        t1 = clock();
        Frangi::imerode(data1d, N, M, P, PARA.somaradius, E8);
        t2 = clock();
        cout << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

        cout << "imgaussian("<< PARA.somaradius <<") " << flush;
        t1 = clock();
        Frangi::imgaussian(E8, N, M, P, PARA.somaradius);
        t2 = clock();
        cout << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

        cout << "maxentropy_th() " << flush;
        t1 = clock();
        unsigned char E8th = maxentropy_th(E8, size); // threshold eroded & blurred image
        t2 = clock();
        cout << (int)E8th << "  " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

        soma_extraction1(E8, E8th, N, M, P, smap, n0);

        if (saveMidres) {

//            simple_saveimage_wrapper(callback, QString(PARA.inimg_file + "_Soma8.v3dpbd").toStdString().c_str(), E8, in_sz, V3D_UINT8);

            // save smap as v3dpbd
            unsigned char* smap8 = new unsigned char[size];
            int smap_min = INT_MAX, smap_max=-INT_MAX;
            for (long i = 0; i < size; ++i) {
                if (smap[i]<smap_min) smap_min = smap[i];
                if (smap[i]>smap_max) smap_max = smap[i];
            }

            for (long i = 0; i < size; ++i) {
                int val = (smap_max==smap_min)? smap_min : round(((smap[i]-smap_min)/(float)(smap_max-smap_min))*255) ;
                val = (val>255)?255:(val<0)?0:val;
                smap8[i] = (unsigned char)val;
            }

            simple_saveimage_wrapper(callback, QString(PARA.inimg_file + "_Smap.v3dpbd").toStdString().c_str(), smap8, in_sz, V3D_UINT8);

            for (long i = 0; i < size; ++i) smap8[i] = (unsigned char)((smap[i]>0)?255:0);

//            simple_saveimage_wrapper(callback, QString(PARA.inimg_file+"_Soma%1.tif").arg(PARA.somaradius).toStdString().c_str(), smap8, in_sz, V3D_UINT8);

            delete [] smap8; smap8 = 0;

//            save_nodelist(n0, PARA.inimg_file + "_Soma.swc");
        }

        delete [] E8; E8 = 0;

    }
    else { // somaradius==0 means no soma detection
        cout << "no soma detection" << endl;
        // do not extract soma, won't add any node, just reset smap
        for (long i = 0; i < size; ++i) smap[i] = 0;
    }

    Frangi frangiflt(sigs, PARA.zdist, frangi_alfa, frangi_beta, frangi_C, frangi_betaone, frangi_betatwo);

    float* J            = new float[size];
    unsigned char* Vx   = new unsigned char[size];
    unsigned char* Vy   = new unsigned char[size];
    unsigned char* Vz   = new unsigned char[size];
    float Jmin, Jmax;

    if (P>1)    frangiflt.frangi3d(data1d, N, M, P, J, Jmin, Jmax, Vx, Vy, Vz);
    else        frangiflt.frangi2d(data1d, N, M, P, J, Jmin, Jmax, Vx, Vy, Vz);

    unsigned char*   J8   = new unsigned char[size]; // float* J ---> unsigned char* J8
    // convert min-max normalized J to J8 (byte8 image)
    if (abs(Jmax-Jmin)<=FLT_MIN) {
        for (long i = 0; i < size; ++i) {
            J8[i] = (unsigned char) 0;
        }
    }
    else {
        for (long i = 0; i < size; ++i) {
            int val = round(((J[i]-Jmin)/(Jmax-Jmin)) * 255);
            val = (val<0)?0: (val>255)?255:val;
            J8[i] = (unsigned char) val;
        }
    }

    delete[] J; J = 0; // delete float* J and leave the rest to work with J8

    if (saveMidres) {
        QString of = PARA.inimg_file + "_J8.tif";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), J8, in_sz, V3D_UINT8);

        of = PARA.inimg_file + "_VxVyVz.swc";
        export_directionality(of.toStdString(), J8, N, M, P, 10, Vx, Vy, Vz); // i, frangiflt.Vxyz); // Vx, Vy, Vz); // Vx[i]==Vxyz[Vi[i]][0]
    }

    // extract seeds, pick cross-section local maxima, and calculate seed correlation
    SeedExtractor sd(sigs, SIG2RADIUS, P==1);
    Tracker t(sigs, PARA.step, PARA.np, PARA.ni, PARA.kappa, P==1, PARA.znccth, Kc, neff_ratio, PARA.zdist, PARA.nodepervol);
//    t.verbose = true;

    if (saveMidres && false) {
        string savepath = PARA.inimg_file.toStdString() + "_Suwv.swc";
        sd.export_Suwv(savepath);
        savepath = PARA.inimg_file.toStdString() + "_Suv.swc";
        sd.export_Suv(savepath);

        savepath = PARA.inimg_file.toStdString() + "_Off3.swc";
        t.export_off3(savepath);

        savepath = PARA.inimg_file.toStdString() + "_Model.swc";
        t.export_model(savepath);

        savepath = PARA.inimg_file.toStdString() + "_ModelVxyz.swc";
        t.export_model(savepath, true);
    }

    cout << "seed extraction... " << flush;
    vector<seed> seeds_init; // initial list of seeds (location+direction+score+correlation)

    t1 = clock();
    SeedExtractor::extractSeeds(PARA.tolerance, J8, N, M, P, Vx, Vy, Vz, seeds_init);
    t2 = clock();
    cout << seeds_init.size()/1000.0 << "k seeds,  " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

    delete [] J8; J8 = 0;
    delete [] Vx; Vx = 0;
    delete [] Vy; Vy = 0;
    delete [] Vz; Vz = 0;

    cout << "seed selection & sorting... " << flush; // filter those within the soma out, filter those below correlation threshold out, index backward to be able to remove
    float dummy_sig;
    t1 = clock();
    for (long i = seeds_init.size()-1; i >= 0; --i) {
        long j = (int)round(seeds_init[i].z)*N*M+(int)round(seeds_init[i].y)*N+(int)round(seeds_init[i].x);
        if (smap[j]>0)
            seeds_init.erase(seeds_init.begin()+i); // erse element at ith index
        else {
            // calculate correlation
            seeds_init[i].corr = t.znccBBB(seeds_init[i].x, seeds_init[i].y, seeds_init[i].z, seeds_init[i].vx, seeds_init[i].vy, seeds_init[i].vz, data1d, N, M, P, dummy_sig);

            if (seeds_init[i].corr<PARA.znccth)
                seeds_init.erase(seeds_init.begin()+i);

        }
    }
    t2 = clock();
    cout << seeds_init.size()/1000.0 << "k seeds, " << ((t2-t1)/(double)CLOCKS_PER_SEC) << " sec." << endl;

    // sort: seeds sorted by correlation so that the tracing starts from those with the highest corr
    // arrange list indices by score (start from the one with highest score)
    vector<long> si(seeds_init.size());
    for (long i = 0; i < si.size(); ++i) si[i] = i;
    sort(si.begin(), si.end(), CompareSeedCorr(&seeds_init));

    vector<seed> seeds;
    for (long i = 0; i < si.size(); ++i) {
        seeds.push_back(seeds_init[si[i]]);
    }

    seeds_init.clear();
    si.clear();

    if (saveMidres) {
        string savepath = PARA.inimg_file.toStdString() + "_Seeds.swc";
        SeedExtractor::export_seeds(seeds, savepath, Node::GREEN, 10);
//        savepath = PARA.inimg_file.toStdString() + "_Seeds_score.log";
//        SeedExtractor::export_seeds_score(seeds,savepath);
//        savepath = PARA.inimg_file.toStdString() + "_Seeds_corr.log";
//        SeedExtractor::export_seeds_corr(seeds, savepath);
    }

//    if (1) {cout << "DONE"<<endl; return;}

    long** ioff = new long*[size];
    for (long i = 0; i < size; ++i) {

        int x  = i%N;
        int z  = i/(N*M);
        int y  = i/N-z*M;

        if (vol==1) ioff[i] = 0; // no neighbouring voxels
        else { // 5,9,11,19,27
            ioff[i] = new long[vol-1];

            if (vol==5 || vol==9 || vol==11 || vol==19 || vol==27) { // 4+ neighbors
                ioff[i][0] = (long)(z*N*M+y                *N+clampi(x-1,0,N-1));
                ioff[i][1] = (long)(z*N*M+y                *N+clampi(x+1,0,N-1));
                ioff[i][2] = (long)(z*N*M+clampi(y-1,0,M-1)*N+                x);
                ioff[i][3] = (long)(z*N*M+clampi(y+1,0,M-1)*N+                x);
            }
            if (          vol==9 || vol==11 || vol==19 || vol==27) { // 8+ neighbors
                ioff[i][4] = (long)(z*N*M+clampi(y-1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][5] = (long)(z*N*M+clampi(y-1,0,M-1)*N+clampi(x+1,0,N-1));
                ioff[i][6] = (long)(z*N*M+clampi(y+1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][7] = (long)(z*N*M+clampi(y+1,0,M-1)*N+clampi(x+1,0,N-1));
            }
            if (                    vol==11 || vol==19 || vol==27) { // 10+ neighbors
                ioff[i][8] = (long)(clampi(z-1,0,P-1)*N*M+y*N+x);
                ioff[i][9] = (long)(clampi(z+1,0,P-1)*N*M+y*N+x);
            }
            if (                               vol==19 || vol==27) { // 18+ neighbors
                ioff[i][10] = (long)(clampi(z-1,0,P-1)*N*M+y                *N+clampi(x-1,0,N-1));
                ioff[i][11] = (long)(clampi(z-1,0,P-1)*N*M+y                *N+clampi(x+1,0,N-1));
                ioff[i][12] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y-1,0,N-1)*N+                x);
                ioff[i][13] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y+1,0,N-1)*N+                x);
                ioff[i][14] = (long)(clampi(z+1,0,P-1)*N*M+y                *N+clampi(x-1,0,N-1));
                ioff[i][15] = (long)(clampi(z+1,0,P-1)*N*M+y                *N+clampi(x+1,0,N-1));
                ioff[i][16] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y-1,0,N-1)*N+                x);
                ioff[i][17] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y+1,0,N-1)*N+                x);
            }
            if (                                          vol==27) { // 26 neighbors
                ioff[i][18] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y-1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][19] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y-1,0,M-1)*N+clampi(x+1,0,N-1));
                ioff[i][20] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y+1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][21] = (long)(clampi(z-1,0,P-1)*N*M+clampi(y+1,0,M-1)*N+clampi(x+1,0,N-1));
                ioff[i][22] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y-1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][23] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y-1,0,M-1)*N+clampi(x+1,0,N-1));
                ioff[i][24] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y+1,0,M-1)*N+clampi(x-1,0,N-1));
                ioff[i][25] = (long)(clampi(z+1,0,P-1)*N*M+clampi(y+1,0,M-1)*N+clampi(x+1,0,N-1));
            }
        }
    }

    cout << "tracing...\n" << flush;
    int trace_count = 0;
    unsigned char* npervol_map = new unsigned char[size];
    int* nidx_map = new int[size];
    for (long i = 0; i < size; ++i) {npervol_map[i] = 0; nidx_map[i] = 0;}
    bool TRACING_VERBOSE = true;
    for (long i = 0; i < seeds.size(); ++i) { // go through all the seeds and trace out of those that are above min score

        // export reconstruction and midresults each 20% of the total (presentation)
        if (false && saveMidres && i>0 && i%(seeds.size()/(int)20)==0) {

            reconstruct(n0, PARA.inimg_file, QString("%1percent").arg( 10*   (i/(seeds.size()/(int)20)) ));

            QString of = PARA.inimg_file + "_TraceDensity_"+QString("%1").arg(i/(seeds.size()/(int)10))+".tif";
            simple_saveimage_wrapper(callback, of.toStdString().c_str(), npervol_map, in_sz, V3D_UINT8);
        }

        long si = (int)round(seeds[i].z)*N*M+(int)round(seeds[i].y)*N+(int)round(seeds[i].x);
        if ((int)npervol_map[si]<PARA.nodepervol) { // trace from the one that had sufficient score (higher than PARA.scmin)
//            if (smap[si]==0) { // trace from the seed if it was out of soma
//                if (t.znccBBB(seeds[i].x, seeds[i].y, seeds[i].z, seeds[i].vx, seeds[i].vy, seeds[i].vz, data1d, N, M, P, dummy_sig)>=PARA.znccth) {

                    trace_count++; // increment trace counter

                    if (TRACING_VERBOSE) {
                        printf("\nTrace: %6d\t [%4.1f, %4.1f, %4.1f]\t sc=%6.2f\t corr=%3.2f\t |n0|=%6d[%d]\t progress %3.2f%%\t",
                               trace_count,
                               seeds[i].x, seeds[i].y, seeds[i].z, seeds[i].score, seeds[i].corr, n0.size(), INT_MAX, (100.0*i)/seeds.size());
                        fflush(stdout);
                    }

                    t.trackPos(seeds[i], data1d, n0, N, M, P, smap, npervol_map, vol, ioff, nidx_map);

                    if (false && trace_count<=1 && saveMidres) { // presentation only, this block can go out later
                        string savepath = PARA.inimg_file.toStdString()+"_Track_"+ QString("%1").arg(trace_count).toStdString()+"a.swc";
                        t.export_track(savepath);

                        savepath = PARA.inimg_file.toStdString()+"_Corr_"+ QString("%1").arg(trace_count).toStdString()+"a.log";
                        t.export_trackcorr(savepath);
                    }

                    t.trackNeg(seeds[i], data1d, n0, N, M, P, smap, npervol_map, vol, ioff, nidx_map);

                    if (false && trace_count<=1 && saveMidres) { // presentation only, this block can go out later
                        string savepath = PARA.inimg_file.toStdString()+"_Track_"+ QString("%1").arg(trace_count).toStdString()+"b.swc";
                        t.export_track(savepath);
                        savepath = PARA.inimg_file.toStdString()+"_Corr_"+ QString("%1").arg(trace_count).toStdString()+"b.log";
                        t.export_trackcorr(savepath);
                    }

                    if (trace_count>MAX_TRACE_COUNT) break;
                    if (n0.size()>=INT_MAX) break; // node references are int, use long if more are necessary

//                }
//                else cout << "NOT t.znccBBB(seeds[i].x, seeds[i].y, seeds[i].z, seeds[i].vx, seeds[i].vy, seeds[i].vz, data1d, N, M, P, dummy_sig)>=PARA.znccth" << endl;
//            }
//            else cout << "NOT smap[si]==0" << endl;
        }
    }

    cout << "\n-----\n" << ((100.0*trace_count)/seeds.size()) << "% seeds used " << endl;

    delete [] smap; smap = 0; // smap[] used to bound the traces
    for (long i = 0; i < size; ++i) {
        delete [] ioff[i];
        ioff[i] = 0;
    }
    delete [] ioff; ioff=0;

    if (saveMidres) {
        QString of = PARA.inimg_file + "_TraceDensity.tif";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), npervol_map, in_sz, V3D_UINT8);
    }

    delete [] npervol_map; npervol_map = 0;
    delete [] nidx_map; nidx_map = 0;

    reconstruct(n0, PARA.inimg_file, "");

}

//    vector<long>     seed_loci; // 3d loc. index of the seeds
//    seed_loci.clear();
//    vector<float>   seed_wght; // 3d loc. weight of the seeds
//    seed_wght.clear();
//    int * nmap     = new int[size]; // accumulated track map
//    for (long i = 0; i < size; ++i) { nmap[i] = -1;}
//    for (long i = 0; i < size; ++i) {
//        nmap[i] = -1; // false;
//        float znorm = ((float)data1d[i]-globalmin)/((float)globalmax-globalmin);
//        if (znorm>=PARA.bth) {
//            seed_loci.push_back(i);
//            seed_wght.push_back(pow(znorm,3));
//        }
//    }
//    long seed_loci_init = seed_loci.size();
//    if (saveMidres) { // save seed_loci
//        unsigned char * seedmap = new unsigned char[size];
//        for (long i = 0; i < size; ++i) seedmap[i] = 0;
//        for (long i = 0; i < seed_loci.size(); ++i) seedmap[seed_loci[i]] = 255;
//        QString of = PARA.inimg_file + "_SEED_BEG.v3dpbd";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), seedmap, in_sz, V3D_UINT8);
//        delete [] seedmap; seedmap = 0;
//    }

//    while (seed_count <= MAX_TRACE_COUNT && seed_loci.size()>0) {

//        // get the seed index by random sampling using cummulative sum of weights (csw)
//        vector<float> tc; tc.clear();
//        for (int i = 0; i < seed_wght.size(); ++i)
//            tc.push_back(seed_wght[i] + ((tc.size()>0)?tc[tc.size()-1]:0));

//        vector<int> locxyz;
//        Tracker::sampleN(tc, 1, locxyz);
//        int ti  = seed_loci[locxyz[0]]; // 'ti' is stack index
//        int tx  = ti%N;
//        int tz  = ti/(N*M);
//        int ty  = ti/N-tz*M;

//        printf(" T: %5d \t [%d,%d,%d]\t%3.2f%\t|n1|=%d\n", seed_count, tx, ty, tz, ((100.0*seed_loci.size())/seed_loci_init), n0.size());

//        t.trackNew(tx, ty, tz, data1d, n0, N, M, P, nmap, 5);

//        nmap[tz*(N*M)+ty*N+tx] = 0; // remove the seed loc from list later

//        // update seedloc/seedcsw list excluding those that were covered by nmap
//        for (int i = seed_loci.size()-1; i >= 0; --i) { // backward when deleting by indexes
//            if (nmap[seed_loci[i]]>=0) {
//                seed_loci.erase(seed_loci.begin()+i);
//                seed_wght.erase(seed_wght.begin()+i);
//            }
//        }

//        seed_count++;

//    } // while() tracks initiated with random locations

//    cout << "\nn1 bidirectional? " << is_bidirectional(n0) << endl;

//    if (saveMidres) { // save seed_loci after
//        unsigned char * seedmap = new unsigned char[size];
//        for (long i = 0; i < size; ++i) seedmap[i] = 0;
//        for (long i = 0; i < seed_loci.size(); ++i) seedmap[seed_loci[i]] = 255;
//        QString of = PARA.inimg_file + "_SEED_END.v3dpbd";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), seedmap, in_sz, V3D_UINT8);
//        delete [] seedmap; seedmap = 0;
//    }

//    if (saveMidres) { // export tmap
//        unsigned char * t = new unsigned char[size];
//        for (long i = 0; i < size; ++i) t[i] = ((nmap[i]>=0)?255:0);
//        QString of = PARA.inimg_file + "_NMAP.v3dpbd";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), t, in_sz, V3D_UINT8);
//        delete [] t; t = 0;
//    }

//    vector<Node> ntree; // ntree is the tree with the nodes obtained after clustering
//    Node dd; ntree.clear(); ntree.push_back(dd); // dummy node at i=0

//    vector<int>  nlist2ntree;
//    for (int i = 0; i < n0.size(); ++i) nlist2ntree.push_back(-1); // mark all as unlabelled

//    float dist2, d2;
//    for (int i = 1; i < n0.size(); ++i) {

//        if (nlist2ntree[i]==-1) { // pick the first unlabelled

//            nlist2ntree[i] = ntree.size(); // needs node to be added in this branch

//            dist2 = pow(PARA.sig2radius*n0[i].sig,2);

//            Node nadd(n0[i].x, n0[i].y, n0[i].z, n0[i].sig, Node::APICAL_DENDRITE);

//            int cnt = 1;
//            float alfa, beta;
//            for (int j = i+1; j < n0.size(); ++j) { // check the rest for it's neighbours
//                d2 = pow(n0[i].x-n0[j].x,2);
//                if (d2<dist2) {
//                    d2 += pow(n0[i].y-n0[j].y,2);
//                    if (d2<dist2) {
//                        d2 += pow(n0[i].z-n0[j].z,2); // zDist?
//                        if (d2<dist2) {

//                            nlist2ntree[j] = nlist2ntree[i]; // it's not going to initiate any new node
//                            cnt++;
//                            alfa = (float)(cnt-1)/cnt;
//                            beta = 1.0/cnt;
//                            nadd.x = alfa * nadd.x + beta * n0[j].x; // iterative mean of node values
//                            nadd.y = alfa * nadd.y + beta * n0[j].y;
//                            nadd.z = alfa * nadd.z + beta * n0[j].z;
//                            nadd.sig = alfa * nadd.sig + beta * n0[j].sig;

//                        }
//                    }
//                }
//            }

//            ntree.push_back(nadd); // add new node to ntree
////            cout << "|ntree|=" << ntree.size() << endl;

//        }
//    }

    // ntree linking based on the nlist2ntree mapping
//    for (int i = 1; i < n0.size(); ++i) {
//        for (int j = 0; j < n0[i].nbr.size(); ++j) {
//            int A = nlist2ntree[ i               ];
//            int B = nlist2ntree[ n0[i].nbr[j] ];
//            ntree[A].nbr.push_back(B);
//        }
//    }

//    cout << "\nntree bidirecitional? " << is_bidirectional(ntree) << endl;
//    for (int i = 1; i < ntree.size(); ++i) {
//        // remove double neighbourhoods from the neighbour list for each node
//        sort(ntree[i].nbr.begin(), ntree[i].nbr.end()); ntree[i].nbr.erase(unique(ntree[i].nbr.begin(), ntree[i].nbr.end()), ntree[i].nbr.end());
//        // remove self linkages
//        int pos = find(ntree[i].nbr.begin(), ntree[i].nbr.end(), i) - ntree[i].nbr.begin();
//        if (pos>=0 && pos<ntree[i].nbr.size()) ntree[i].nbr.erase(ntree[i].nbr.begin()+pos); // remove at pos
//    }

//    vector<Node> tree; // output tree
//    bfs1(ntree, tree); // tree will be used to generate the reconstruction
//    NeuronTree recon; // output tree reconstruction
//    recon.name = signature;

//    for (int i = 1; i < tree.size(); ++i) {
//        NeuronSWC n;
//        n.n = n.nodeinseg_id = i;
//        n.type = tree[i].type;
//        n.x = tree[i].x;
//        n.y = tree[i].y;
//        n.z = tree[i].z;
//        n.r = tree[i].sig;
//        n.parent = tree[i].nbr[0];
//        recon.listNeuron.append(n);
//    }
//    QString rec_name = PARA.inimg_file + "_Advantra.swc";
//    writeSWC_file(rec_name.toStdString().c_str(), recon);
//    delete [] nmap; nmap = 0;


/*
void non_blurring_fix(vector<Node> nX, vector<Node>& nY, float KRAD2, int MAXITER, float EPSILON2) {

    int checkpoint = round(nX.size()/10.0);

    float conv[4], next[4]; // x y z sig

    nY.clear();
    nY = nX;

    float x2, y2, z2, d2; // distance components
    int iter, cnt;

    for (int i = 1; i < nX.size(); ++i) {

        if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

        conv[0] = nX[i].x;
        conv[1] = nX[i].y;
        conv[2] = nX[i].z;
        conv[3] = nX[i].sig;

        iter = 0;

        do {

            cnt = 0;

            next[0] = 0;
            next[1] = 0;
            next[2] = 0;
            next[3] = 0;

            for (int j = 1; j < nX.size(); ++j) {
                x2 = pow(nX[j].x-conv[0],2);
                if (x2<=KRAD2) {
                    y2 = pow(nX[j].y-conv[1],2);
                    if (x2+y2<=KRAD2) {
                        z2 = pow(nX[j].z-conv[2],2);
                        if (x2+y2+z2<=KRAD2) {
                            next[0] += nX[j].x;
                            next[1] += nX[j].y;
                            next[2] += nX[j].z;
                            next[3] += nX[j].sig;
                            cnt++;
                        }
                    }
                }
            }

            if (cnt==0) cout << "WRONG!!!" << endl;

            next[0] /= cnt;
            next[1] /= cnt;
            next[2] /= cnt;
            next[3] /= cnt;

            d2 = pow(next[0]-conv[0],2) + pow(next[1]-conv[1],2) + pow(next[2]-conv[2],2);

            conv[0] = next[0];
            conv[1] = next[1];
            conv[2] = next[2];
            conv[3] = next[3];

            iter++;

        }
        while(iter<MAXITER && d2>EPSILON2);

        nY[i].x = conv[0];
        nY[i].y = conv[1];
        nY[i].z = conv[2];
        nY[i].sig = conv[3];

    }

    cout << endl;

}
*/


//    if (0) {
//        vector<Node> dummy_nodelist;
//        QString name = "Advantra";
//        QString comment =
//                "email: miro@braincadet.com\n#params:\n#channel="+QString("%1").arg(PARA.channel)+
//                "\n#sigmas="+QString::fromStdString(PARA.neuritesigmas)+       // QString::fromStdString(PARA.sigmas)+//QString("%1").arg(442)+
//                "\n#msiter="+QString("%1").arg(PARA.msiter)+
//                "\n#seedscmin="+QString("%1").arg(PARA.seedscmin)+
//                "\n#znccth="+QString("%1").arg(PARA.znccth)+
//                "\n#kappa="+QString("%1").arg(PARA.kappa)+
//                "\n#step="+QString("%1").arg(PARA.step)+
//                "\n#ni="+QString("%1").arg(PARA.ni)+
//                "\n#np="+QString("%1").arg(PARA.np)+
//                "\n#zdist="+QString("%1").arg(PARA.zdist)+
//                "\n#blurring="+QString("%1").arg(PARA.blurring)+
//                "\n#grouping="+QString("%1").arg(PARA.grouping)+
//                "\n#nodepp="+QString("%1").arg(PARA.nodepervol)+
//                "\n#vol="+QString("%1").arg(PARA.vol)+
//                "\n#------------------------"+
//                "\n#Kc="+QString("%1").arg(Kc)+
//                "\n#neff_ratio="+QString("%1").arg(neff_ratio)+
//                "\n#frangi_alfa="+QString("%1").arg(frangi_alfa)+
//                "\n#frangi_beta="+QString("%1").arg(frangi_beta)+
//                "\n#frangi_C="+QString("%1").arg(frangi_C)+
//                "\n#frangi_betaone="+QString("%1").arg(frangi_betaone)+
//                "\n#frangi_betatwo="+QString("%1").arg(frangi_betatwo)+
//                "\n#MAX_TRACE_COUNT="+QString("%1").arg(MAX_TRACE_COUNT)+
//                "\n#EPSILON2="+QString("%1").arg(EPSILON2)+
//                "\n#SIG2RADIUS="+QString("%1").arg(SIG2RADIUS)+
//                "\n#TRACE_RSMPL="+QString("%1").arg(TRACE_RSMPL)+
//                "\n#GROUP_RADIUS="+QString("%1").arg(GROUP_RADIUS)+
//                "\n#ENFORCE_SINGLE_TREE="+QString("%1").arg(ENFORCE_SINGLE_TREE)+
//                "\n#TREE_SIZE_MIN="+QString("%1").arg(TREE_SIZE_MIN)+
//                "\n#TAIL_SIZE_MIN="+QString("%1").arg(TAIL_SIZE_MIN);
//        save_nodelist(dummy_nodelist, PARA.inimg_file + "_Advantra.swc", -1, 1, name, comment);
//        return;
//    }
