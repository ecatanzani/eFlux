#include "sbi.h"

#include <iostream>

sbi::sbi()
{
    // Sizing
    sizer();
    // Create TTree
    tree = std::make_shared<TTree>("SBItree", "SBI Tree");
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

    lat_geo = 0;
    lon_geo = 0;
    rad_geo = 0;
    ra_zenith = 0;
    dec_zenith = 0;
    ra_scz = 0;
    dec_scz = 0;
    ra_scx = 0;
    dec_scx = 0;
    ra_scy = 0;
    dec_scy = 0;
    mjd = 0;
    met = 0;
    glat = 0;
    glon = 0;

    position = std::vector<double>(_att_v_val, 0);
    velocity = std::vector<double>(_att_v_val, 0);
    angle = std::vector<double>(_att_v_val, 0);
    angular_velocity = std::vector<double>(_att_v_val, 0);
    star_sensor_A = std::vector<double>(_star_sensor_v_val, 0);
    star_sensor_B = std::vector<double>(_star_sensor_v_val, 0);

    dipoleMoment = 0;
    L = 0;
    bEast = 0;
    bNorth = 0;
    bDown = 0;
    bAbs = 0;
    bEquator = 0;
    R = 0;
    verticalRigidityCutoff = 0;
    invariantLatitude = 0;

    BVect = std::vector<double>(_att_v_val, 0);

    SunRa = 0;
    SunDec = 0;
    insaa = false;
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

void sbi::Fill(
    std::shared_ptr<container> sec_info,
    std::shared_ptr<DmpEvtAttitude> attitude)
{
    goodsbi = sec_info->goodsbi;
    run = sec_info->run;
    second = sec_info->second;
    lvtime = sec_info->lvtime;
    nev = sec_info->nev;
    firstev = sec_info->firstev;
    lastev = sec_info->lastev;
    ntrig = sec_info->ntrig;
    trigenable = sec_info->trigenable;
    ntrigperiod = sec_info->ntrigperiod;
    trigperiodenable = sec_info->trigperiodenable;
    stkstatus = sec_info->stkstatus;
    bgostatus = sec_info->bgostatus;
    psdstatus = sec_info->psdstatus;
    nudstatus = sec_info->nudstatus;
    ntracks = sec_info->ntracks;
    nbgohits = sec_info->nbgohits;
    npsdhits = sec_info->npsdhits;
    nnudhits = sec_info->nnudhits;

    if (attitude != nullptr)
    {
        lat_geo = attitude->lat_geo;
        lon_geo = attitude->lon_geo;
        rad_geo = attitude->rad_geo;
        ra_zenith = attitude->ra_zenith;
        dec_zenith = attitude->dec_zenith;
        ra_scz = attitude->ra_scz;
        dec_scz = attitude->dec_scz;
        ra_scx = attitude->ra_scx;
        dec_scx = attitude->dec_scx;
        ra_scy = attitude->ra_scy;
        dec_scy = attitude->dec_scy;
        glat = attitude->glat;
        glon = attitude->glon;
        
        for (int idx_pos = 0; idx_pos < _att_v_val; ++idx_pos)
        {
            position[idx_pos] = attitude->Position[idx_pos];
            velocity[idx_pos] = attitude->Velocity[idx_pos];
            angle[idx_pos] = attitude->Angle[idx_pos];
            angular_velocity[idx_pos] = attitude->AngularVelocity[idx_pos];
            star_sensor_A[idx_pos] = attitude->StarSensorA[idx_pos];
            star_sensor_B[idx_pos] = attitude->StarSensorB[idx_pos];
            BVect[idx_pos] = attitude->BVect[idx_pos];
        }
        
        SunRa = attitude->SunRa;
        SunDec = attitude->SunDec;
        dipoleMoment = attitude->dipoleMoment;
        L = attitude->L;
        bEast = attitude->bEast;
        bNorth = attitude->bNorth;
        bDown = attitude->bDown;
        bAbs = attitude->bAbs;
        bEquator = attitude->bEquator;
        R = attitude->R;
        verticalRigidityCutoff = attitude->verticalRigidityCutoff;
        invariantLatitude = attitude->invariantLatitude;
        mjd = attitude->mjd;
    }

    tree->Fill();
}

void sbi::Write(TFile& outfile)
{
    outfile.cd();
    tree->Write();
}