#ifndef SBI_H
#define SBI_H

#include <memory>
#include <vector>

#define _trig 5
#define _msdeadtime 3.0725
#define _att_v_val 3
#define _star_sensor_v_val 4

#include "TTree.h"
#include "TFile.h"

#include "DmpEvtHeader.h"
#include "DmpEvtAttitude.h"

#include "container.h"

class sbi
{
public:
    sbi();
    ~sbi(){};
    
    void Reset();
    void Init(
        const unsigned int _second,
        const unsigned int _firstev);
    bool SetSBIStatus(
        const DmpEvtHeader* header,
        const DmpEvtAttitude* attitude,
        const unsigned int idx);
    void Fill(
        std::shared_ptr<container> sec_info,
        std::shared_ptr<DmpEvtAttitude> attitude);
    void Write(TFile& outfile);

private:
    void sizer();
    void branch_tree();

    // Tree
	std::shared_ptr<TTree> tree;

    // SBI vars
    bool goodsbi = false;
    unsigned int run = 0;
    unsigned int second = 0;
    double lvtime = 0;
    unsigned int nev = 0;
    unsigned int firstev = 0;
    unsigned int lastev = 0;
    std::vector<unsigned int> ntrig;
    std::vector<bool> trigenable;
    unsigned int ntrigperiod = 0;
    bool trigperiodenable = false;
    bool stkstatus = false;
    bool bgostatus = false;
    bool psdstatus = false;
    bool nudstatus = false;
    unsigned int ntracks = 0;
    unsigned int nbgohits = 0;
    unsigned int npsdhits = 0;
    unsigned int nnudhits = 0;

    // Attitude vars
    double lat_geo = 0;
    double lon_geo = 0;
    double rad_geo = 0;
    double ra_zenith = 0;
    double dec_zenith = 0;
    double ra_scz = 0;
    double dec_scz = 0;
    double ra_scx = 0;
    double dec_scx = 0;
    double ra_scy = 0;
    double dec_scy = 0;
    double mjd = 0;
    double met = 0;
    double glat = 0;
    double glon = 0;

    std::vector<double> position;
    std::vector<double> velocity; 
    std::vector<double> angle;
    std::vector<double> angular_velocity;

    std::vector<double> star_sensor_A;
    std::vector<double> star_sensor_B;

    double dipoleMoment = 0;
    double L = 0;
    double bEast = 0;
    double bNorth = 0;
    double bDown = 0;
    double bAbs = 0;
    double bEquator = 0;
    double R = 0;
    double verticalRigidityCutoff = 0;
    double invariantLatitude = 0;

    std::vector<double> BVect;
    double SunRa = 0;
    double SunDec = 0;
    bool insaa = false;
};

#endif