#ifndef NODE_H
#define NODE_H
#include <vector>

class Node
{
public:

    float x;  // location
    float y;
    float z;

    float vx; // direction
    float vy;
    float vz;

    float corr; // correlation with the template at this location and scale
    float sig;  // standard deviation of the gaussian cross section
    int   type;

    std::vector<int> nbr; // list of indexes of the neighbouring Nodes from node list

    // types as in neuromorpho.org description
    static int NOTHING, SOMA, AXON;
    static int BASAL_DENDRITE;
    static int APICAL_DENDRITE;
    static int FORK;
    static int END;
    static int UNDEFINED;

    // some colours in vaa3d swc plot
    static int WHITE, BLACK, RED, BLUE, PINK, MAGENTA, YELLOW, GREEN, OCRE, GREEN_LIGHT,
    PINK_LIGHT, MAGENTA_LIGHT, VIOLET, PINK1, GREEN_SHARP, BLUE_LIGHT, GREEN1, OCRE_LIGHT;

    Node();
    Node(float _x, float _y, float _z, float _sig);
    Node(float _x, float _y, float _z, float _sig, int _type);
    Node(float _x, float _y, float _z, float _vx, float _vy, float _vz, float _corr, float _sig, int _type);
//    Node(Node& _n);
    ~Node();

    void print();

};

#endif // NODE_H
