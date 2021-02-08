#include "sbi.h"

#include <iostream>

sbi::sbi()
{
    // Sizing
    sizer();
    // Create TTree
    tree = std::make_unique<TTree>("SBItree", "SBI Tree");
    // Branch TTree
    branch_tree();
}

void sbi::sizer()
{
    ntrig = std::vector<unsigned int>(_trig, 0);
    trigenable = std::vector<bool>(_trig, false);
    position = std::vector<double>(_att_v_val, 0);
    velocity = std::vector<double>(_att_v_val, 0);
    angle = std::vector<double>(_att_v_val, 0);
    angular_velocity = std::vector<double>(_att_v_val, 0);
    star_sensor_A = std::vector<double>(_star_sensor_v_val, 0);
    star_sensor_B = std::vector<double>(_star_sensor_v_val, 0);
    BVect = std::vector<double>(_att_v_val, 0);
}

void sbi::branch_tree()
{
    tree->Branch("goodsbi", &goodsbi, "goodsbi/O");
    tree->Branch("run", &run, "run/i");
    tree->Branch("second", &second, "second/i");
    tree->Branch("lvtime", &lvtime, "lvtime/D");
    tree->Branch("nev", &nev, "nev/i");
    tree->Branch("firstev", &firstev, "firstev/i");
    tree->Branch("lastev", &lastev, "lastev/i");
    tree->Branch("ntrig", &ntrig);
    tree->Branch("trigenable", &trigenable);
    tree->Branch("ntrigperiod", &ntrigperiod, "ntrigperiod/i");
    tree->Branch("trigperiodenable", &trigperiodenable, "trigperiodenable/O");
    tree->Branch("stkstatus", &stkstatus, "stkstatus/O");
    tree->Branch("bgostatus", &bgostatus, "bgostatus/O");
    tree->Branch("psdstatus", &psdstatus, "psdstatus/O");
    tree->Branch("nudstatus", &nudstatus, "nudstatus/O");
    tree->Branch("ntracks", &ntracks, "ntracks/i");
    tree->Branch("nbgohits", &nbgohits, "nbgohits/i");
    tree->Branch("npsdhits", &npsdhits, "npsdhits/i");
    tree->Branch("nnudhits", &nnudhits, "nnudhits/i");

    tree->Branch("lat_geo", &lat_geo, "lat_geo/D");
    tree->Branch("lon_geo", &lon_geo, "lon_geo/D");
    tree->Branch("rad_geo", &rad_geo, "rad_geo/D");
    tree->Branch("ra_zenith", &ra_zenith, "ra_zenith/D");
    tree->Branch("dec_zenith", &dec_zenith, "dec_zenith/D");
    tree->Branch("ra_scz", &ra_scz, "ra_scz/D");
    tree->Branch("dec_scz", &dec_scz, "dec_scz/D");
    tree->Branch("ra_scx", &ra_scx, "ra_scx/D");
    tree->Branch("dec_scx", &dec_scx, "dec_scx/D");
    tree->Branch("ra_scy", &ra_scy, "ra_scy/D");
    tree->Branch("dec_scy", &dec_scy, "dec_scy/D");
    tree->Branch("mjd", &mjd, "mjd/D");
    tree->Branch("met", &met, "met/D");
    tree->Branch("glat", &glat, "glat/D");
    tree->Branch("glon", &glon, "glon/D");
    tree->Branch("position", &position);
    tree->Branch("velocity", &velocity);
    tree->Branch("angle", &angle);
    tree->Branch("angular_velocity", &angular_velocity);
    tree->Branch("star_sensor_A", &star_sensor_A);
    tree->Branch("star_sensor_B", &star_sensor_B);
    tree->Branch("dipoleMoment", &dipoleMoment, "dipoleMoment/D");
    tree->Branch("L", &L, "L/D");
    tree->Branch("bEast", &bEast, "bEast/D");
    tree->Branch("bNorth", &bNorth, "bNorth/D");
    tree->Branch("bDown", &bDown, "bDown/D");
    tree->Branch("bAbs", &bAbs, "bAbs/D");
    tree->Branch("bEquator", &bEquator, "bEquator/D");
    tree->Branch("R", &R, "R/D");
    tree->Branch("verticalRigidityCutoff", &verticalRigidityCutoff, "verticalRigidityCutoff/D");
    tree->Branch("invariantLatitude", &invariantLatitude, "invariantLatitude/D");
    tree->Branch("BVect", &BVect);
    tree->Branch("SunRa", &SunRa, "SunRa/D");
    tree->Branch("SunDec", &SunDec, "SunDec/D");
    tree->Branch("insaa", &insaa, "insaa/O");
}

void sbi::Reset()
{
    goodsbi = false;
    run = 0;
    second = 0;
    lvtime = 0;
    nev = 0;
    firstev = 0;
    lastev = 0;
    ntrig = std::vector<unsigned int>(_trig, 0);
    trigenable = std::vector<bool>(_trig, false);
    ntrigperiod = 0;
    trigperiodenable = false;
    stkstatus = false;
    bgostatus = false;
    psdstatus = false;
    nudstatus = false;
    ntracks = 0;
    nbgohits = 0;
    npsdhits = 0;
    nnudhits = 0;
}

void sbi::Init(
    const unsigned int _second,
    const unsigned int _firstev)
{
    Reset();
    second = _second;
    firstev = _firstev;
}

bool sbi::SetSBIStatus(
    const DmpEvtHeader *header,
    const DmpEvtAttitude *attitude,
    const unsigned int idx)
{
    bool status = true;

    // Header check
    if (!header)
    {
        std::cerr << "\nWARNING: Header obj not found at entry [" << idx << "] -- skipping\n";
        status = false;
    }
    else
        goodsbi = true;

    // Attitude check
    if (!attitude)
    {
        std::cerr << "\nWARNING: Attitude obj not found at entry [" << idx << "] -- skipping\n";
        goodsbi = false;
    }

    return status;
}
