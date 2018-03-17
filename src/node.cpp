/*
Copyright (C) Erasmus MC. Permission to use this software and corresponding documentation for educational, research, and not-for-profit purposes, without a fee and without a signed licensing agreement, is granted, subject to the following terms and conditions.
IT IS NOT ALLOWED TO REDISTRIBUTE, SELL, OR LEASE THIS SOFTWARE, OR DERIVATIVE WORKS THEREOF, WITHOUT PERMISSION IN WRITING FROM THE COPYRIGHT HOLDER. THE COPYRIGHT HOLDER IS FREE TO MAKE VERSIONS OF THE SOFTWARE AVAILABLE FOR A FEE OR COMMERCIALLY ONLY.
IN NO EVENT SHALL THE COPYRIGHT HOLDER BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OF ANY KIND WHATSOEVER, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
THE COPYRIGHT HOLDER SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE EXPRESS OR IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND CORRESPONDING DOCUMENTATION IS PROVIDED "AS IS". THE COPYRIGHT HOLDER HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
 */

#include "node.h"
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <cstdio>

int Node::NOTHING = 0;
int Node::SOMA = 1;
int Node::AXON = 2;
int Node::BASAL_DENDRITE = 3;
int Node::APICAL_DENDRITE = 4;
int Node::FORK = 5;
int Node::END = 6;
int Node::UNDEFINED = 7;

int Node::WHITE     = 0;
int Node::BLACK     = 1;
int Node::RED       = 2;
int Node::BLUE      = 3;
int Node::PINK      = 4;
int Node::MAGENTA   = 5;
int Node::YELLOW    = 6;
int Node::GREEN     = 7;
int Node::OCRE      = 8;
int Node::GREEN_LIGHT = 9;
int Node::PINK_LIGHT = 10;
int Node::MAGENTA_LIGHT = 11;
int Node::VIOLET    = 12;
int Node::PINK1     = 13;
int Node::GREEN_SHARP = 14;
int Node::BLUE_LIGHT = 15;
int Node::GREEN1    = 16;
int Node::OCRE_LIGHT = 17;

Node::Node() {
    x = 0;
    y = 0;
    z = 0;
    vx = 0;
    vy = 0;
    vz = 0;
    sig = 0;
    corr = -FLT_MAX;
    type = UNDEFINED;
    nbr.clear();
}

Node::Node(float _x, float _y, float _z, float _sig) {
    x = _x;
    y = _y;
    z = _z;
    vx = 0;
    vy = 0;
    vz = 0;
    sig = _sig;
    corr = -FLT_MAX;
    type = Node::UNDEFINED;
    nbr.clear();
}

Node::Node(float _x, float _y, float _z, float _sig, int _type) {
    x = _x;
    y = _y;
    z = _z;
    vx = 0;
    vy = 0;
    vz = 0;
    sig = _sig;
    corr = -FLT_MAX;
    type = _type;
    nbr.clear();
}

Node::Node(float _x, float _y, float _z, float _vx, float _vy, float _vz, float _corr, float _sig, int _type) {
    x = _x;
    y = _y;
    z = _z;
    vx = _vx;
    vy = _vy;
    vz = _vz;
    sig = _sig;
    corr = _corr;
    type = _type;
    nbr.clear();
}

//Node::Node(Node& _n) {
//    x = _n.x;
//    y = _n.y;
//    z = _n.z;
//    vx = _n.vx;
//    vy = _n.vy;
//    vz = _n.vz;
//    sig = _n.sig;
//    corr = _n.corr;
//    type = _n.type;
//    nbr.clear();
//    nbr = _n.nbr;
//}

Node::~Node(){}

void Node::print(){
    printf("[x,y,z]=(%6.2f,%6.2f,%6.2f)\t[vx,vy,vz]=(%6.2f,%6.2f,%6.2f)\tsig=%4.2f\tcorr=%3.2f\ttype=%2d\t", x, y, z, vx, vy, vz,sig,corr,type);
//    std::cout<<"[x,y,z]=("<<x<<","<<y<<","<<z<<")\t[vx,vy,vz]=("<<vx<<","<<vy<<","<<vz<<")\tsig="<<sig<<"\tcorr="<<corr<<"\ttype="<<type<<std::flush;
    std::cout<<"nbr=["<<std::flush;
    for (int i = 0; i < nbr.size(); ++i) {
        std::cout<<nbr[i]<<((i<nbr.size()-1)?",":"")<<std::flush;
    }
    std::cout<<"]"<<std::endl;
}
