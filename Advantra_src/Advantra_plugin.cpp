/* Advantra_plugin.cpp
* Tool for automatic neuron reconstruction from microscopy image stacks.
* 2015-8-19 : initial version by Miroslav Radojevic
* 2016-4-4  : new tracker implementation and track merging mechanism by Miroslav Radojevic
* 2016-5-24 : prefiltering added (frangi based)
* 2016-6-16 : added refinement module and export tools

Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
*/
 
#include "v3d_message.h"
#include "basic_surf_objs.h"
#include "nf_dialog.h"
#include "tracker.h"
#include "frangi.h"
#include "seed.h"
#include "node.h"
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

#include "Advantra_plugin.h"
Q_EXPORT_PLUGIN2(Advantra, Advantra);

using namespace std;

static V3DLONG channel = 1;      // default channel (hardcoded)

// input parameter default values
int     sigmax      = 3;     // largest gaussian cross-section standard deviation used for prefiltering and tracking
int     msiter      = 2;     // number of mean-shift iterations when blending together the trace lines
float   seedscmin   = 1;     // FLT_MIN, minimum seed score in [0-255], average of the differences towards local neighbours
float   znccth      = 0.4;   // correlation threshold (stops tracing)
float   kappa       = 2;     // von mises circular normal distribution
int     step        = 2;     // prediction step
int     ni          = 20;    // nr. iterations
int     np          = 30;    // nr. particles
float   zdist       = 3.0;   // scaling along z
int     blurring    = 0;     // 0: non-blurring, 1: blurring mean-shift implementation
int     grouping    = 1;     // 0: cylindrical, 1: spherical, spherical by default
int     nodepp      = 5;     // limit amount of nodes per pixel, used for the trace suppression (rationalize the tracing)

int     nrInputParams = 12;           // as in input_PARA struct
bool    saveMidres = false;        // save midresults

float    Kc          = 10.0;    // likelihood factor
float    neff_ratio  = 0.8;     // resampling boundary (ratio of total # pcles)
float    sig_min     = 1;       //
float    sig_stp     = 1;       //

float    frangi_alfa    = .5;   //
float    frangi_beta    = .5;   //
float    frangi_C       = 500;          //
float    frangi_betaone = .5;           //
float    frangi_betatwo = 15;           //

int      MAX_TRACE_COUNT    = INT_MAX;  // INT_MAX
float    EPSILON2           = FLT_MIN;  // FLT_MIN
//bool     SEED_SUPP          = false;    // enable usage of the trace for the suppression

float   SIG2RADIUS          = 1.5;      // (used for the neighbourhood estimates in refinement and grouping)
float   TRACE_RSMPL         = 2.0;      // component trace resampling step

bool    ENFORCE_SINGLE_TREE = false;    // 1 largest tree as output
int     TREE_SIZE_MIN       = 10;       // trees with less than TREE_SIZE_MIN nodes are discarded (used if ENFORCE_SINGLE_TREE=false)
int     TAIL_SIZE_MIN       = 2;        // tails (junction--endpoint) with less than TAIL_SIZE_MIN are discarded

struct input_PARA
{
    QString inimg_file;
    V3DLONG channel;
    int     sigmax;     // 1 scale, tube diameter
    int     msiter;     // 2 number of refinement iterations
    float   seedscmin;  // 3 seed score minimum to be taken as seed (score is average difference towards the neighbours)
    float   znccth;     // 4 correlation threshold
    float   kappa;      // 5 von mises kappa
    int     step;       // 6 prediction step
    int     ni;         // 7 number of iterations
    int     np;         // 8 number of particles
    float   zdist;      // 9 the distance between layers in pixels
    int     blurring;   // 10 averaging method
    int     grouping;   // 11 node grouping method: 0: cylinder, 1: sphere
    int     nodepp;     // 12 node per pixel density
};

struct Pxyz {
    Pxyz(int x1, int y1, int z1) : x(x1),y(y1),z(z1) {}
    float x, y, z;
};

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
    printf("---- ADVANTRA usage ----\n");
    printf("vaa3d -x Advantra -f advantra_func -i <inimg_file> -p <sigmax msiter sig2radius znccth kappa step ni np zdist>\n");
    printf("inimg_file          The input image.\n");
    printf("sigmax              Max Gaussian cross-section sigma.\n");
    printf("msiter              # refinement iterations.\n");
    printf("seedscmin           minimum seed score (uint8 greylevels 0-255).\n");
    printf("znccth              Correlation threshold.\n");
    printf("kappa               Von Mises variance.\n");
    printf("step                Prediction step.\n");
    printf("ni                  # trace iterations.\n");
    printf("np                  # trace particles.\n");
    printf("zdist               z layer dist.\n");
    printf("blurring            blurring (0-no,1-yes).\n");
    printf("grouping            grouping method (0-cylinder,1-sphere).\n");
    printf("nodepp              nodes per pixel trace density limit.\n");
    printf("outswc_file         Will be named automatically based on the input image file name, so you don't have to specify it.\n\n");
}

void Advantra::domenu(const QString &menu_name, V3DPluginCallback2 &callback, QWidget *parent)
{
	if (menu_name == tr("advantra_menu"))
	{
        bool bmenu = true;
        input_PARA PARA;

        // take the default params
        PARA.channel = channel;
        PARA.sigmax = sigmax;
        PARA.msiter = msiter;
        PARA.seedscmin = seedscmin;
        PARA.znccth = znccth;
        PARA.kappa   = kappa;
        PARA.step = step;
        PARA.ni = ni;
        PARA.np = np;
        PARA.zdist = zdist;
        PARA.blurring = blurring;
        PARA.grouping = grouping;
        PARA.nodepp = nodepp;

        vector<string> items;
        items.push_back("max_gauss_cs_std");
        items.push_back("nr_refinements");
        items.push_back("min_seed_score");
        items.push_back("min_corr");
        items.push_back("von_mises_variance");
        items.push_back("prediction_step");
        items.push_back("nr_trace_iterations");
        items.push_back("nr_trace_particles");
        items.push_back("z_layer_dist");
        items.push_back("blurring");
        items.push_back("grouping_method");
        items.push_back("node_per_pix");

        // initialization
        vector<string> inits;
        inits.push_back(QString::number(PARA.sigmax).toStdString().c_str());
        inits.push_back(QString::number(PARA.msiter).toStdString().c_str());
        inits.push_back(QString::number(PARA.seedscmin).toStdString().c_str());
        inits.push_back(QString::number(PARA.znccth).toStdString().c_str());
        inits.push_back(QString::number(PARA.kappa).toStdString().c_str());
        inits.push_back(QString::number(PARA.step).toStdString().c_str());
        inits.push_back(QString::number(PARA.ni).toStdString().c_str());
        inits.push_back(QString::number(PARA.np).toStdString().c_str());
        inits.push_back(QString::number(PARA.zdist).toStdString().c_str());
        inits.push_back(QString::number(PARA.blurring).toStdString().c_str());
        inits.push_back(QString::number(PARA.grouping).toStdString().c_str());
        inits.push_back(QString::number(PARA.nodepp).toStdString().c_str());

        CommonDialog dialog(items, inits);
        dialog.setWindowTitle(title);
        if(dialog.exec() != QDialog::Accepted) return;

        dialog.get_num("max_gauss_cs_std", PARA.sigmax);
        dialog.get_num("nr_refinements", PARA.msiter);
        dialog.get_num("min_seed_score", PARA.seedscmin);
        dialog.get_num("min_corr", PARA.znccth);
        dialog.get_num("von_mises_variance", PARA.kappa);
        dialog.get_num("prediction_step", PARA.step);
        dialog.get_num("nr_trace_iterations", PARA.ni);
        dialog.get_num("nr_trace_particles", PARA.np);
        dialog.get_num("z_layer_dist", PARA.zdist);
        dialog.get_num("blurring", PARA.blurring);
        dialog.get_num("grouping_method", PARA.grouping);
        dialog.get_num("node_per_pix", PARA.nodepp);

        // check input
        if(PARA.sigmax<=0 || PARA.sigmax>10){v3d_msg(QObject::tr("sigmax out of range")); return;}
        if(PARA.msiter<=0 || PARA.msiter>10){v3d_msg(QObject::tr("msiter out of range")); return;}
        if(PARA.seedscmin< 0 || PARA.seedscmin>255){v3d_msg(QObject::tr("seedscmin out of range")); return;}
        if(PARA.znccth<0 || PARA.znccth>1){v3d_msg(QObject::tr("znccth out of range")); return;}
        if(PARA.kappa<0 || PARA.kappa>5){v3d_msg(QObject::tr("kappa out of range")); return;}
        if(PARA.step<1){v3d_msg(QObject::tr("step out of range")); return;}
        if(PARA.ni<=0){v3d_msg(QObject::tr("ni out of range")); return;}
        if(PARA.np<=0){v3d_msg(QObject::tr("np out of range")); return;}
        if(PARA.zdist<1){v3d_msg(QObject::tr("zdist out of range")); return;}
        if(PARA.blurring<0 || PARA.blurring>1){v3d_msg(QObject::tr("blurring out of range")); return;}
        if(PARA.grouping<0 || PARA.grouping>1){v3d_msg(QObject::tr("grouping out of range")); return;}
        if(PARA.nodepp<=0 || PARA.nodepp>50){v3d_msg(QObject::tr("nodepp out of range")); return;}

        reconstruction_func(callback,parent,PARA,bmenu);

	}
	else
	{
		v3d_msg(tr("Tool for automatic neuron reconstruction from microscopy image stacks.. "
            "Developed by Miroslav Radojevic\n2015-8-19 first version (BigNeuron submission)\n2016-3-29 full re-design"));
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
        PARA.channel    = channel; // hardcoded not submitted as parameter
        PARA.sigmax     = (paras.size() >= k+1)   ? atoi(paras[k])              : sigmax;       k++;
        PARA.msiter     = (paras.size() >= k+1)   ? atoi(paras[k])              : msiter;       k++;
        PARA.seedscmin  = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : seedscmin;    k++;
        PARA.znccth     = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : znccth;       k++;
        PARA.kappa      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : kappa;        k++;
        PARA.step       = (paras.size() >= k+1)   ? atoi(paras[k])              : step;         k++;
        PARA.ni         = (paras.size() >= k+1)   ? atoi(paras[k])              : ni;           k++;
        PARA.np         = (paras.size() >= k+1)   ? atoi(paras[k])              : np;           k++;
        PARA.zdist      = (paras.size() >= k+1)   ? QString(paras[k]).toFloat() : zdist;        k++;
        PARA.blurring   = (paras.size() >= k+1)   ? atoi(paras[k])              : blurring;     k++;
        PARA.grouping   = (paras.size() >= k+1)   ? atoi(paras[k])              : grouping;     k++;
        PARA.nodepp     = (paras.size() >= k+1)   ? atoi(paras[k])              : nodepp;       k++;

        // check user input
        if(PARA.sigmax<=0 || PARA.sigmax>10){v3d_msg(QObject::tr("sigmax out of range")); return 0;}
        if(PARA.msiter<=0 || PARA.msiter>10){v3d_msg(QObject::tr("msiter out of range")); return 0;}
        if(PARA.seedscmin<0 || PARA.seedscmin>255){v3d_msg(QObject::tr("seedscmin out of range"));return 0;}
        if(PARA.znccth<0 || PARA.znccth>1){v3d_msg(QObject::tr("znccth out of range"));return 0;}
        if(PARA.kappa<0 || PARA.kappa>5){v3d_msg(QObject::tr("kappa out of range")); return 0;}
        if(PARA.step<1){v3d_msg(QObject::tr("step out of range")); return 0;}
        if(PARA.ni<=0){v3d_msg(QObject::tr("ni out of range")); return 0;}
        if(PARA.np<=0){v3d_msg(QObject::tr("np out of range")); return 0;}
        if(PARA.zdist<1){v3d_msg(QObject::tr("zdist out of range")); return 0;}
        if(PARA.blurring<0 || PARA.blurring>1){v3d_msg(QObject::tr("blurring out of range")); return 0;}
        if(PARA.grouping<0 || PARA.grouping>1){v3d_msg(QObject::tr("grouping out of range")); return 0;}
        if(PARA.nodepp<=0 || PARA.nodepp>50){v3d_msg(QObject::tr("nodepp out of range")); return 0;}

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

bool ischecked(seed _s, bool* _smap, int _w, int _h, int _l) {

    int x = round(_s.x);
    if (x<0 || x>=_w) return true;

    int y = round(_s.y);
    if (y<0 || y>=_h) return true;

    int z = round(_s.z);
    if (z<0 || z>=_l) return true;

    return _smap[z*_w*_h+y*_w+x];
}

int get_undiscovered2(int dist[], int dist_length){
    for (int i = 1; i < dist_length; i++) {
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

    int dist[nlist.size()];
    int nmap[nlist.size()];
    int parent[nlist.size()];

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

    while ((seed = get_undiscovered2(dist, nlist.size()))>0) {

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
            n.type = treecnt+10;//nicer set of colours in vaa3d viz

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

void save_nodelist(vector<Node> nlist, QString swcname, int type=-1, float sig2r=1, QString signature=""){

    // NOTE: nodes have more than one parent, export will duplicate node ids to show all neighbouring connections (nodes have bidirectional links)

    NeuronTree recnodes;

    if (signature.compare("")!=0) {
        recnodes.name = signature;
    }

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

vector<Node> extract_trees(vector<Node> ntreeX, int min_size) {
    // take input TREELIST and remove all trees with less than min_size nodes, removed nodes belong to one connected tree, tag those that are to be removed

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
    for (int i = 1; i < ntreeX_2d.size(); ++i) {
        for (int j = 0; j < ntreeX_2d[i].nbr.size(); ++j) { // existing link i -- ntreeX_2d[i].nbr[j]
            ntreeX_2d[   ntreeX_2d[i].nbr[j]   ].nbr.push_back(      i       ); // add reverse link
        }
    }

    vector<int> Nnbr(ntreeX_2d.size(), 0);
    for (int i = 1; i < ntreeX_2d.size(); ++i) {
        Nnbr[i] = ntreeX_2d[i].nbr.size();
    }

//    for (int i = 0; i < Nnbr.size(); ++i) cout << Nnbr[i] << " " << flush; cout << endl;

    vector<bool> remove_map(ntreeX_2d.size(), false);

    for (int i = 1; i < ntreeX_2d.size(); ++i) {
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
            ntreeY.push_back(nYi); // nodes are copied, neighbours stay as in ntreeX vector indexing
        }
    }

    // remap the linking
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

//void save_nodelist_tree(vector<Node> nlist, QString swcname, QString signature="", int type=-1, float sig2r=1) {
//    save_nodelist(ntree, swcname, type, sig2r, signature);
//}

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

/*
void blurring_fix(vector<Node> nX, vector<Node>& nY, float KRAD2, int MAXITER, float EPSILON2) {

    int checkpoint = round(nX.size()/10.0);

    float conv[nX.size()][4];
    float next[nX.size()][4];

    nY.clear();
    nY = nX;

    for (int i = 1; i < nX.size(); ++i) {
        conv[i][0] = nX[i].x;
        conv[i][1] = nX[i].y;
        conv[i][2] = nX[i].z;
        conv[i][3] = nX[i].sig;
    }

    float x2, y2, z2, d2, d2_max;
    int iter = 0;
    int cnt = 0;

    do {

        cout << "iter = " << iter << endl;

        // each iteration will shift the whole sample = blurring
        d2_max = -FLT_MAX;

        for (int i = 1; i < nX.size(); ++i) {

            if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

            cnt = 0;                // count neighbours

            next[i][0] = 0;
            next[i][1] = 0;
            next[i][2] = 0;
            next[i][3] = 0;

            for (int j = 1; j < nX.size(); ++j) {
                x2 = pow(conv[j][0]-conv[i][0],2);
                if (x2<=KRAD2) {
                    y2 = pow(conv[j][1]-conv[i][1],2);
                    if (x2+y2<=KRAD2) {
                        z2 = pow(conv[j][2]-conv[i][2],2);
                        if (x2+y2+z2<=KRAD2) {
                            next[i][0] += conv[j][0];
                            next[i][1] += conv[j][1];
                            next[i][2] += conv[j][2];
                            next[i][3] += conv[j][3];
                            cnt++;
                        }
                    }
                }
            }

            if (cnt==0) cout << "WRONG!!!" << endl;

            next[i][0] /= cnt;
            next[i][1] /= cnt;
            next[i][2] /= cnt;
            next[i][3] /= cnt;

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

    for (int i = 1; i < nY.size(); ++i) {
        nY[i].x = conv[i][0];
        nY[i].y = conv[i][1];
        nY[i].z = conv[i][2];
        nY[i].sig = conv[i][3];
    }

}
*/
void interpolate(vector<Node>& nX, float step) {

    // interpolate all inter-node links with the step size
    vector< vector<bool> > chk(nX.size()); // disable accessing the same pair of nodes twice
    for (int i = 0; i < nX.size(); ++i) {
        vector<bool> chk1(nX[i].nbr.size(), false);
        chk.push_back(chk1);
    }

    int check_limit = nX.size();
    float s1 = 1.0;

    for (int i = 1; i < check_limit; ++i) {
        for (int j = 0; j < nX[i].nbr.size(); ++j) {
            if (!chk[i][j]) {

                int i1 = nX[i].nbr[j];
                int j1 = find(nX[i1].nbr.begin(), nX[i1].nbr.end(), i) - nX[i1].nbr.begin();

                if (j1<nX[i1].nbr.size()) {

                    chk[i][j]       = true;
                    chk[i1][j1]     = true;

                    float vnorm = sqrt(pow(nX[i1].x-nX[i].x,2) + pow(nX[i1].y-nX[i].y,2) + pow(nX[i1].z-nX[i].z,2));

                    float vx = (nX[i1].x-nX[i].x)/vnorm;
                    float vy = (nX[i1].y-nX[i].y)/vnorm;
                    float vz = (nX[i1].z-nX[i].z)/vnorm;
                    int N = ceil(vnorm/s1);

                    // add the subsampling
                    for (int k = 1; k < N; ++k) {

                        // add the node,only location is used in the refinement stage currently
                        Node nXadd(nX[i].x+k*(vnorm/N)*vx,
                                   nX[i].y+k*(vnorm/N)*vy,
                                   nX[i].z+k*(vnorm/N)*vz,
                                   vx,
                                   vy,
                                   vz,
                                   0.5*(nX[i].corr+nX[i1].corr), // average correlation of the limit nodes
                                   0.5*(nX[i].sig+nX[i1].sig), // interpolate sig similarly
                                   nX[i].type
                                   );
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

    cout << ((float)nX.size()/check_limit)*100.0 << "% node # after interpolation" << endl;

}

void refine_blurring(vector<Node> nX, vector<Node>& nY, float SIG2RAD, int MAXITER, float EPSILON2) {

    int checkpoint = round(nX.size()/10.0);

    float conv[nX.size()][4];
    float next[nX.size()][4];

    nY.clear();
    nY = nX;

    for (int i = 1; i < nX.size(); ++i) {
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

        for (int i = 1; i < nX.size(); ++i) {

            if (i%checkpoint==0) cout << (i/checkpoint)*10 << "%  " << flush;

            cnt = 0;                // count neighbours

            next[i][0] = 0;
            next[i][1] = 0;
            next[i][2] = 0;
            next[i][3] = 0;

            r2 = pow(SIG2RAD * conv[i][3],2);

            for (int j = 1; j < nX.size(); ++j) {
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

            if (cnt==0) cout << "WRONG!!!" << endl;

            next[i][0] /= cnt;
            next[i][1] /= cnt;
            next[i][2] /= cnt;
            next[i][3] /= cnt;

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

    for (int i = 1; i < nY.size(); ++i) {
        nY[i].x = conv[i][0];
        nY[i].y = conv[i][1];
        nY[i].z = conv[i][2];
        nY[i].sig = conv[i][3];
    }

}

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

void non_blurring(vector<Node> nX, vector<Node>& nY, float SIG2RAD, int MAXITER, float EPSILON2) {

    // mean-shift (non-blurring) uses flexible neighbourhood scaled with respect to the node's sigma

    int checkpoint = round(nX.size()/10.0);

    float conv[4], next[4]; // x y z sig

    nY.clear();
    nY = nX;

    float x2, y2, z2, d2, r2;
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

            next[0] = 0;        // local mean is the follow-up location
            next[1] = 0;
            next[2] = 0;
            next[3] = 0;

            r2 = pow(SIG2RAD * conv[3],2);

            for (int j = 1; j < nX.size(); ++j) {
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

            if (cnt==0) cout << "WRONG!!!" << endl;

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
        while(iter<MAXITER && d2>EPSILON2);

        nY[i].x     = conv[0];
        nY[i].y     = conv[1];
        nY[i].z     = conv[2];
        nY[i].sig   = conv[3];

    }

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

    int XtoY[nX.size()]; // map
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

void filter_internode_dist(vector<Node>& nX, float dth) { // , vector<Node>& nY
    // will remove all the linkages (assumes bidirectional) that are above an euclidean distance threshold dth, removes the links only
    int count = 0;
    for (int i = 1; i < nX.size(); ++i) {
        for (int j = nX[i].nbr.size()-1; j >= 0; --j) {
            int i1 = nX[i].nbr[j];
            float dist = sqrt(pow(nX[i1].x-nX[i].x,2)+pow(nX[i1].y-nX[i].y,2)+pow(nX[i1].z-nX[i].z,2));
            if (dist>dth) {
                count++;
                nX[i].nbr.erase(nX[i].nbr.begin() + j);// remove element at index j, i->i1 link
                int i2 = find(nX[i1].nbr.begin(), nX[i1].nbr.end(), i) - nX[i1].nbr.begin(); // i value at nY[i1].nbr[i2]
                if (i2<nX[i1].nbr.size()) nX[i1].nbr.erase(nX[i1].nbr.begin()+i2); // remove i value at idx i2, i1->i link
            }
        }
    }

    cout << count << " links > " << dth << endl;
}

void group1(vector<Node> nX, vector<Node>& nY, float sig2rad=1) {

    nX[0].corr = FLT_MAX; // so that the dummy node gets index 0
    vector<int> indices(nX.size());
    for (int i = 0; i < indices.size(); ++i) indices[i] = i;
    sort(indices.begin(), indices.end(), CompareIndicesByNodeCorrVal(&nX));

    vector<int> X2Y(nX.size(), -1);
    X2Y[0] = 0; // first one is with max. correlation

    nY.clear();
    Node nY0(nX[0]);
    nY.push_back(nY0);

    for (int i = 1; i < indices.size(); ++i) {

        int ci = indices[i];

        if (X2Y[ci]!=-1) continue; // skip if it was added to a group already

        X2Y[ci] = nY.size();
        Node nYi(nX[ci]); // nX[ci] is the grabbed
        float grp_size = 1;

        float r2 = sig2rad * nX[ci].sig;
        for (int j = 1; j < nX.size(); ++j) { // check the rest that was not groupped
            if (j!=ci && X2Y[j]==-1) {
                float x2 = pow(nX[j].x-nX[ci].x,2);
                if (x2<=r2) {
                    float y2 = pow(nX[j].y-nX[ci].y,2);
                    if (x2+y2<=r2) {
                        float z2 = pow(nX[j].z-nX[ci].z   ,2);
                        if (x2+y2+z2<=r2) {

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

        nY.push_back(nYi);
    }

    for (int i = 1; i < nY.size(); ++i) {
        for (int j = 0; j < nY[i].nbr.size(); ++j) {
            nY[i].nbr[j] = X2Y[ nY[i].nbr[j] ];
        }
    }

}

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

void check_nbr(vector<Node>& nX) {
    // correct neighbourhood pointers
    // - ensure linkings are bidirectional
    for (int i = 1; i < nX.size(); ++i) {

        // remove double neighbourhoods from the neighbour list
        sort(nX[i].nbr.begin(), nX[i].nbr.end());
        nX[i].nbr.erase(unique(nX[i].nbr.begin(), nX[i].nbr.end()), nX[i].nbr.end());

        // remove self linkages
        int pos = find(nX[i].nbr.begin(), nX[i].nbr.end(), i) - nX[i].nbr.begin();
        if (pos>=0 && pos<nX[i].nbr.size())
            nX[i].nbr.erase(nX[i].nbr.begin()+pos); // remove at pos

    }

    // ensure linkings are bidirectional, add if not
    for (int i = 1; i < nX.size(); ++i) { // first index is dummy
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

void get_node_density(  vector<Node> nX, vector<float>& d) {

    d.clear();

//    float dmin =  FLT_MAX;
//    float dmax = -FLT_MAX;

    for (int i = 1; i < nX.size(); ++i) {
        float w = 1;
        for (int j = 1; j < nX.size(); ++j) {
            if (j!=i) {
                float r2 = pow(3*nX[j].sig,2);
                float x2 = pow(nX[j].x-nX[i].x,2);
                if (x2<=r2) {
                    float y2 = pow(nX[j].y-nX[i].y,2);
                    if (x2+y2<=r2) {
                        float z2 = pow(nX[j].z-nX[i].z,2);
                        if (x2+y2+z2<=r2) {
                            w += exp(-(x2+y2+z2)/(2*pow(nX[j].sig,2)));
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

    vector< vector<bool> > chk(nX.size());
    for (int i = 0; i < nX.size(); ++i) {
        vector<bool> chk1(nX[i].nbr.size(), false);
        chk.push_back(chk1);
    }

    l.clear();

    for (int i = 1; i < nX.size(); ++i) {
        for (int j = 0; j < nX[i].nbr.size(); ++j) {
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
    for (int i = 1; i < nX.size(); ++i) {
        c.push_back(nX[i].corr);
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

void export_directionality(string swcpath, unsigned char* J, int w, int h, int l, unsigned char Jth, float* Vx, float* Vy, float* Vz) {
    ofstream f;
    f.open(swcpath.c_str(), ofstream::out | ofstream::trunc);

    int count = 1;
    for (long i = 0; i < (w*h*l); ++i) { // go through all the voxels
        if (J[i]>Jth) { // show vectors for those with tubularity above Jth [0-255]
            f<<count<<" "<<Node::OCRE_LIGHT<<" "<<i%w<<" "<<(i/w-(i/(w*h))*h)<<" "<<(i/(w*h))<<" "<<0.1<<" "<<(-1)<<endl;
            count++;
            float v = (J[i]*10.0/255);
            f<<count<<" "<<Node::OCRE_LIGHT<<" "<<(i%w+v*Vx[i])<<" "<<(i/w-(i/(w*h))*h+v*Vy[i])<<" "<<(i/(w*h)+v*Vz[i])<<" "<<0.1<<" "<<(count-1)<<endl;
            count++;
        }
    }

    f.close();

    cout<<"exported: "<<swcpath<<endl;
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
    cout<<"sigmax = "    <<PARA.sigmax              <<endl;
    cout<<"msiter = "    <<PARA.msiter              <<endl;
    cout<<"seedscmin = " <<PARA.seedscmin           <<endl;
    cout<<"znccth = "    <<PARA.znccth              <<endl;
    cout<<"kappa = "     <<PARA.kappa               <<endl;
    cout<<"step = "      <<PARA.step                <<endl;
    cout<<"ni = "        <<PARA.ni                  <<endl;
    cout<<"np = "        <<PARA.np                  <<endl;
    cout<<"zdist = "     <<PARA.zdist               <<endl;
    cout<<"blurring="    <<PARA.blurring            <<endl;
    cout<<"grouping="    <<PARA.grouping            <<endl;
    cout<<"nodepp="      <<PARA.nodepp              <<endl;
    cout<<"-------------------------------------"   <<endl;
    cout<<"channel = "   << PARA.channel            <<endl;
    cout<<"saveMidres = "<< saveMidres              <<endl;
    cout<<"Kc = "        << Kc                      <<endl;
    cout<<"neff_ratio = "<< neff_ratio              <<endl;
    cout<<"sig_min = "   << sig_min                 <<endl;
    cout<<"sig_stp = "   << sig_stp                 <<endl;
    cout<<"-------------------------------------"   <<endl;

    long size = N * M * P; // N : width, M : height, P : nr. layers
    unsigned char*   J8   = new unsigned char[size];
    unsigned char *  Sc   = new unsigned char[size];
    float *          Vx   = new float[size];
    float *          Vy   = new float[size];
    float *          Vz   = new float[size];
//    unsigned char* Dxx8   = new unsigned char[size]; // for test purpose
//    unsigned char* Dxy8   = new unsigned char[size];
//    unsigned char* Dyy8   = new unsigned char[size];
//    unsigned char* L18   = new unsigned char[size];
//    unsigned char* L28   = new unsigned char[size];

    Frangi frangiflt(sig_min, sig_stp, PARA.sigmax, PARA.zdist, frangi_alfa, frangi_beta, frangi_C, frangi_betaone, frangi_betatwo);

    if (P>1)    frangiflt.filter3d(data1d, N, M, P, J8, Sc, Vx, Vy, Vz);
    else        frangiflt.filter2d(data1d, N, M, P, J8, Sc, Vx, Vy, Vz);

        //  Dxx8, Dxy8, Dyy8, L18, L28
//        for (long var = 0; var < size; ++var) Vz[var] = 0;
//        QString of = PARA.inimg_file + "_Dxx8.tif";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), Dxx8, in_sz, V3D_UINT8);
//        delete [] Dxx8; Dxx8 = 0;
//        of = PARA.inimg_file + "_Dyy8.tif";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), Dyy8, in_sz, V3D_UINT8);
//        delete [] Dyy8; Dyy8 = 0;
//        of = PARA.inimg_file + "_Dxy8.tif";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), Dxy8, in_sz, V3D_UINT8);
//        delete [] Dxy8; Dxy8 = 0;
//        of = PARA.inimg_file + "_L18.tif";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), L18, in_sz, V3D_UINT8);
//        delete [] L18; L18 = 0;
//        of = PARA.inimg_file + "_L28.tif";
//        simple_saveimage_wrapper(callback, of.toStdString().c_str(), L28, in_sz, V3D_UINT8);
//        delete [] L28; L28 = 0;

    if (J8==NULL) {cout<<"\nJ8==0"<<endl; return;}

    if (saveMidres) {
        QString of = PARA.inimg_file + "_J8.v3dpbd";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), J8, in_sz, V3D_UINT8);

        of = PARA.inimg_file + "_Sc.v3dpbd";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), Sc, in_sz, V3D_UINT8);

        of = PARA.inimg_file + "_VxVyVz.swc";
        export_directionality(of.toStdString(), J8, N, M, P, 20, Vx, Vy, Vz);
    }

    // extract seeds, pick cross-section local maxima
    SeedExtractor sd(sig_min, sig_stp, PARA.sigmax, SIG2RADIUS);

    if (saveMidres) {
        string savepath = PARA.inimg_file.toStdString() + "_Suwv.swc";
        sd.export_Suwv(savepath);
        savepath = PARA.inimg_file.toStdString() + "_Suv.swc";
        sd.export_Suv(savepath);
    }

    vector<seed> seeds;
//    unsigned char J8th = 0;
    if (P>1) sd.extract3d(0,PARA.seedscmin, J8, N,M,P,Sc, Vx,Vy,Vz, seeds);
    else     sd.extract2d(0,PARA.seedscmin, J8, N,M,P,Sc, Vx,Vy,    seeds);

    delete [] J8; J8 = 0;
    delete [] Sc; Sc = 0;
    delete [] Vx; Vx = 0;
    delete [] Vy; Vy = 0;
    delete [] Vz; Vz = 0;

    if (saveMidres) {
        string savepath = PARA.inimg_file.toStdString() + "_Seed.swc";
        SeedExtractor::export_seedlist(seeds, savepath, Node::GREEN);
    }

    // instantiate traces from each seed, avoid neighbours and avoid those seeds below some score threshold
    vector<Node> n0;
    Node init_node;
    n0.push_back(init_node); // list index 0 is dummy node

    Tracker t(sig_min, sig_stp, PARA.sigmax, PARA.step, PARA.np, PARA.ni, PARA.kappa, P==1, PARA.znccth, Kc, neff_ratio, PARA.zdist, PARA.nodepp);
//    t.verbose = true;

    if (saveMidres) {

        string savepath;
//        savepath = PARA.inimg_file.toStdString() + "_Off3.swc";
//        t.export_off3(savepath);

        savepath = PARA.inimg_file.toStdString() + "_Model.swc";
        t.export_model(savepath);

        savepath = PARA.inimg_file.toStdString() + "_ModelVxyz.swc";
        t.export_model(savepath, true);
    }

    // while the seed pool is not empty, sample random seeds one by one and trace from there
    // start tracing
    int seed_count = 1;
    int startidx = 0;
    int count_used_seeds = 0;
    bool* chk = new bool[size]; // used to supress the trace seeds
    unsigned char* trc_den = new unsigned char[size];
    for (long i = 0; i < size; ++i) trc_den[i] = 0;

//    // dbg find min/max
//    unsigned char tmax = 0;
//    for (long i = 0; i < size; ++i) {
//        if (trc_den[i]>tmax) tmax=trc_den[i];
//    }
//    cout << "trc_den max = " << (int)tmax << endl;


    for (long i = 0; i < size; ++i) chk[i] = false;
    bool TRACING_VERBOSE = true;

    while (seed_count<=MAX_TRACE_COUNT) {

        // take the first unchecked seed
        while (ischecked(seeds[startidx], chk, N, M, P)) startidx++;

        if (startidx>=seeds.size()) {cout<<"reached seeds.size()"<<endl; break;}

        count_used_seeds++;
        if (TRACING_VERBOSE)
            printf("\nSeed: %5d\t [%4.1f, %4.1f, %4.1f] \t sc=%1.2f% \t |n0|=%d \t %3.2f%%\t", seed_count, seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, seeds[startidx].score, n0.size(), (100.0*startidx)/seeds.size());

        //////////////////////////////////////////////////////////////////
        t.trackPos(seeds[startidx], data1d, n0, N, M, P, NULL, trc_den); // (SEED_SUPP)?chk:NULL

        if (false && saveMidres) {
            string savepath = PARA.inimg_file.toStdString()+"_Track_"+ QString("%1").arg(seed_count).toStdString()+"a.swc";
            t.export_track(savepath);

            savepath = PARA.inimg_file.toStdString()+"_Corr_"+ QString("%1").arg(seed_count).toStdString()+"a.log";
            t.export_trackcorr(savepath);
        }

        //////////////////////////////////////////////////////////////////
        t.trackNeg(seeds[startidx], data1d, n0, N, M, P, NULL, trc_den);// (SEED_SUPP)?chk:NULL

        if (false && saveMidres) {
            string savepath = PARA.inimg_file.toStdString()+"_Track_"+ QString("%1").arg(seed_count).toStdString()+"b.swc";
            t.export_track(savepath);
            savepath = PARA.inimg_file.toStdString()+"_Corr_"+ QString("%1").arg(seed_count).toStdString()+"b.log";
            t.export_trackcorr(savepath);
        }

        seed_count++;

        if (Tracker::NH==1)      Tracker::check1x1(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
        else if (Tracker::NH==2) Tracker::check2x2(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
        else if (Tracker::NH==3) Tracker::check3x3(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
        else if (Tracker::NH==4) Tracker::check4x4(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
        else if (Tracker::NH==5) Tracker::check5x5(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
        else                     Tracker::check3x3(seeds[startidx].x, seeds[startidx].y, seeds[startidx].z, chk, N, M, P);
    }

    cout << "\n-----\n" << ((100.0*count_used_seeds)/seeds.size()) << "% seeds used " << endl;

    if (saveMidres) { // export unit8 image with the map of checked
        unsigned char* checked8 = new unsigned char[size];
        for (long i = 0; i < size; ++i) checked8[i] = (chk[i])?255:0;
        QString of = PARA.inimg_file + "_Checked.v3dpbd";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), checked8, in_sz, V3D_UINT8);
        delete [] checked8; checked8 = 0;
    }

    delete [] chk; chk = 0;

    if (saveMidres) {
        QString of = PARA.inimg_file + "_TraceDensity.tif";
        simple_saveimage_wrapper(callback, of.toStdString().c_str(), trc_den, in_sz, V3D_UINT8);
    }

    delete [] trc_den; trc_den = 0;

    if (saveMidres) save_nodelist(n0, PARA.inimg_file + "_n0.swc");
    if (saveMidres) save_nodelist(compute_trees(n0), PARA.inimg_file + "_n0tree.swc");

    if (saveMidres) {
        // examine n0
        vector<float> n0len;
        get_link_lengths(n0, n0len);
        save_vector(n0len, PARA.inimg_file + "_n0len.log");
        // density log (experimental)
        vector<float> n0den;
        get_node_density(n0, n0den);
        save_vector(n0den, PARA.inimg_file + "_n0den.log");
        // correlation log
        vector<float> n0corr;
        get_node_corr(n0, n0corr);
        save_vector(n0corr, PARA.inimg_file + "_n0corr.log");
    }

    // resample
    interpolate(n0, TRACE_RSMPL);
    if (saveMidres) save_nodelist(n0, PARA.inimg_file + "_n0res.swc");

    // refinement
    vector<Node> n1;
    if (PARA.blurring==1)       refine_blurring(n0, n1, SIG2RADIUS, PARA.msiter, EPSILON2);
    else if (PARA.blurring==0)  non_blurring(n0, n1, SIG2RADIUS, PARA.msiter, EPSILON2);
    if (saveMidres) save_nodelist(n1, PARA.inimg_file + "_n1.swc");
    n0.clear(); // free memory

    //if (saveMidres) save_nodelist(n1B, PARA.inimg_file + "_n1B.swc");
//  if (saveMidres) save_nodelist(nI, PARA.inimg_file + "_nI.swc", Node::VIOLET);
//  non_blurring_fix(n0, nI, KRAD2, MAXITER, EPSILON2);     if (saveMidres) save_nodelist(nI, PARA.inimg_file + "_nI_non_blurring_fix.swc", Node::BLUE);
//  blurring_fix(n0, nI, KRAD2, MAXITER, EPSILON2);         if (saveMidres) save_nodelist(nI, PARA.inimg_file + "_nI_blurring_fix.swc", Node::GREEN);

    if (saveMidres) {
        vector<float> n1len;
        get_link_lengths(n1, n1len);
        save_vector(n1len, PARA.inimg_file + "_n1len.log");
    }

    // remove strechted out linkages after the refinement
    filter_internode_dist(n1, 3*TRACE_RSMPL);

    vector<Node> n2;
    if (PARA.grouping==0)       group0(n1, n2, 2.0, SIG2RADIUS); // hardcoded the grouping scale
    else if (PARA.grouping==1)  group1(n1, n2, 2.0);

    check_nbr(n2); // remove doubles and self-linkages
    if (saveMidres) save_nodelist(n2, PARA.inimg_file + "_n2.swc");
    n1.clear(); // free memory

    vector<Node> n2tree = compute_trees(n2);
    if (saveMidres) save_nodelist(n2tree, PARA.inimg_file + "_n2tree.swc");

//    cout<<"n2tree"<<endl; summarize_tree(n2tree);
    vector<Node> n3tree = (ENFORCE_SINGLE_TREE)?extract_largest_tree(n2tree):extract_trees(n2tree, TREE_SIZE_MIN);
//    cout<<"n3tree"<<endl; summarize_tree(n3tree);
    if (saveMidres) save_nodelist(n3tree, PARA.inimg_file + "_n3tree.swc");
    n2tree.clear();

    vector<Node> n4tree = remove_tails(n3tree, TAIL_SIZE_MIN); // expell tails (end-junction) with less than TAIL_SIZE_MIN nodes
    n3tree.clear();

    // reconstruction export
    QString signature =
            "Advantra\n#author: miro@braincadet.com\n#params:\n#channel="+QString("%1").arg(PARA.channel)+
            "\n#sigmax="+QString("%1").arg(PARA.sigmax)+
            "\n#msiter="+QString("%1").arg(PARA.msiter)+
            "\n#znccth="+QString("%1").arg(PARA.znccth)+
            "\n#kappa="+QString("%1").arg(PARA.kappa)+
            "\n#step="+QString("%1").arg(PARA.step)+
                "\n#ni="+QString("%1").arg(PARA.ni)+
                "\n#np="+QString("%1").arg(PARA.np)+
                "\n#zdist="+QString("%1").arg(PARA.zdist)+
                "\n#------"+
                "\n#Kc="+QString("%1").arg(Kc)+
                "\n#neff_ratio="+QString("%1").arg(neff_ratio)+
                "\n#sig_min="+QString("%1").arg(sig_min)+
                "\n#sig_stp="+QString("%1").arg(sig_stp)+
                "\n#frangi_alfa="+QString("%1").arg(frangi_alfa)+
                "\n#frangi_beta="+QString("%1").arg(frangi_beta)+
                "\n#frangi_C="+QString("%1").arg(frangi_C)+
                "\n#frangi_betaone="+QString("%1").arg(frangi_betaone)+
                "\n#frangi_betatwo="+QString("%1").arg(frangi_betatwo)+
                "\n#MAX_TRACE_COUNT="+QString("%1").arg(MAX_TRACE_COUNT)+
                "\n#EPSILON2="+QString("%1").arg(EPSILON2)+
                "\n#SIG2RADIUS="+QString("%1").arg(SIG2RADIUS)+
                "\n#TRACE_RSMPL="+QString("%1").arg(TRACE_RSMPL)+
                "\n#ENFORCE_SINGLE_TREE="+QString("%1").arg(ENFORCE_SINGLE_TREE)+
                "\n#TREE_SIZE_MIN="+QString("%1").arg(TREE_SIZE_MIN)+
                "\n#TAIL_SIZE_MIN="+QString("%1").arg(TAIL_SIZE_MIN);
    save_nodelist(n4tree, PARA.inimg_file + "_Advantra.swc", -1, 1,signature);

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
