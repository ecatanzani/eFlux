#ifndef SBI_H
#define SBI_H

#include <vector>

#define _trig 5
#define _msdeadtime = 3.0725

class sbi
{
public:
    sbi()  :    ntrig(_trig, 0),
                trigenable(_trig, false) {};
    ~sbi(){};

private:
    bool goodsbi;
    unsigned int run;
    unsigned int second;
    double lvtime;
    unsigned int nev;
    unsigned int firstev;
    unsigned int lastev;
    std::vector<unsigned int> ntrig;
    std::vector<bool> trigenable;
    unsigned int ntrigperiod;
    bool trigperiodenable;
    bool stkstatus;
    bool bgostatus;
    bool psdstatus;
    bool nudstatus;
    unsigned int ntracks;
    unsigned int nbgohits;
    unsigned int npsdhits;
    unsigned int nnudhits;
};

#endif