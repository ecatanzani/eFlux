#include "acceptance_cuts.h"
#include "acceptance.h"

bool checkBGOreco(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0) || (bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0))
    {
        TVector3 orgPosition;
        orgPosition.SetX(simu_primaries->pv_x);
        orgPosition.SetY(simu_primaries->pv_y);
        orgPosition.SetZ(simu_primaries->pv_z);

        TVector3 orgMomentum;
        orgMomentum.SetX(simu_primaries->pvpart_px);
        orgMomentum.SetY(simu_primaries->pvpart_py);
        orgMomentum.SetZ(simu_primaries->pvpart_pz);

        std::vector<double> slope(2, 0);
        std::vector<double> intercept(2, 0);

        slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
        slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
        intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
        intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

        double actual_X = slope[0] * BGO_TopZ + intercept[0];
        double actual_Y = slope[1] * BGO_TopZ + intercept[1];
        double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
        double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

        int position_sensitivity = 30;

        if (fabs(actual_X - topX) > position_sensitivity || fabs(actual_Y - topY) > position_sensitivity)
            return false;
        else
            return true;
    }
    else
        return true;
}

void fillExternalMap(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    TH2D &h_noBGOenergy_real_topMap)
{
    TVector3 orgPosition;
    orgPosition.SetX(simu_primaries->pv_x);
    orgPosition.SetY(simu_primaries->pv_y);
    orgPosition.SetZ(simu_primaries->pv_z);

    TVector3 orgMomentum;
    orgMomentum.SetX(simu_primaries->pvpart_px);
    orgMomentum.SetY(simu_primaries->pvpart_py);
    orgMomentum.SetZ(simu_primaries->pvpart_pz);

    std::vector<double> slope(2, 0);
    std::vector<double> intercept(2, 0);

    slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
    slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
    intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
    intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

    double actual_X = slope[0] * BGO_TopZ + intercept[0];
    double actual_Y = slope[1] * BGO_TopZ + intercept[1];

    h_noBGOenergy_real_topMap.Fill(actual_X, actual_Y);
}

bool geometric_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
    bool passed_geometric_cut = false;

    TVector3 orgPosition;
    orgPosition.SetX(simu_primaries->pv_x);
    orgPosition.SetY(simu_primaries->pv_y);
    orgPosition.SetZ(simu_primaries->pv_z);

#if 0
    // **** Directions Cosines Method

    TVector3 dCos;
    dCos.SetX(simu_primaries->pvpart_cosx);
    dCos.SetY(simu_primaries->pvpart_cosy);
    dCos.SetZ(simu_primaries->pvpart_cosz);

    if (dCos.Z())
    {
        double ratioZ = (BGO_TopZ - orgPosition.Z()) / dCos.Z();
        double actual_X = ratioZ * dCos.X() + orgPosition.X();
        double actual_Y = ratioZ * dCos.Y() + orgPosition.Y();
        if (fabs(actual_X) < BGO_SideXY && fabs(actual_Y) < BGO_SideXY)
            passed_geometric_cut = true;
    }

#else
    // **** Moments Method

    TVector3 orgMomentum;
    orgMomentum.SetX(simu_primaries->pvpart_px);
    orgMomentum.SetY(simu_primaries->pvpart_py);
    orgMomentum.SetZ(simu_primaries->pvpart_pz);

    //auto orgMomentum_theta = orgMomentum.Theta() * TMath::RadToDeg();
    //auto orgMomentum_costheta = cos(orgMomentum.Theta());

    std::vector<double> slope(2, 0);
    std::vector<double> intercept(2, 0);

    slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
    slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
    intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
    intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

    double actual_topX = slope[0] * BGO_TopZ + intercept[0];
    double actual_topY = slope[1] * BGO_TopZ + intercept[1];

    double actual_bottomX = slope[0] * BGO_BottomZ + intercept[0];
    double actual_bottomY = slope[1] * BGO_BottomZ + intercept[1];

    if (fabs(actual_topX) < BGO_SideXY && fabs(actual_topY) < BGO_SideXY &&
        fabs(actual_bottomX) < BGO_SideXY && fabs(actual_bottomY) < BGO_SideXY)
        passed_geometric_cut = true;

#endif

    return passed_geometric_cut;
}

bool geometric_top_cut(const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries)
{
    bool passed_geometric_cut = false;

    TVector3 orgPosition;
    orgPosition.SetX(simu_primaries->pv_x);
    orgPosition.SetY(simu_primaries->pv_y);
    orgPosition.SetZ(simu_primaries->pv_z);

#if 0
    // **** Directions Cosines Method

    TVector3 dCos;
    dCos.SetX(simu_primaries->pvpart_cosx);
    dCos.SetY(simu_primaries->pvpart_cosy);
    dCos.SetZ(simu_primaries->pvpart_cosz);

    if (dCos.Z())
    {
        double ratioZ = (BGO_TopZ - orgPosition.Z()) / dCos.Z();
        double actual_X = ratioZ * dCos.X() + orgPosition.X();
        double actual_Y = ratioZ * dCos.Y() + orgPosition.Y();
        if (fabs(actual_X) < BGO_SideXY && fabs(actual_Y) < BGO_SideXY)
            passed_geometric_cut = true;
    }

#else
    // **** Moments Method

    TVector3 orgMomentum;
    orgMomentum.SetX(simu_primaries->pvpart_px);
    orgMomentum.SetY(simu_primaries->pvpart_py);
    orgMomentum.SetZ(simu_primaries->pvpart_pz);

    //auto orgMomentum_theta = orgMomentum.Theta() * TMath::RadToDeg();
    //auto orgMomentum_costheta = cos(orgMomentum.Theta());

    std::vector<double> slope(2, 0);
    std::vector<double> intercept(2, 0);

    slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
    slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
    intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
    intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

    double actual_topX = slope[0] * BGO_TopZ + intercept[0];
    double actual_topY = slope[1] * BGO_TopZ + intercept[1];

    if (fabs(actual_topX) < BGO_SideXY && fabs(actual_topY))
        passed_geometric_cut = true;

#endif

    return passed_geometric_cut;
}

void evaluateTopBottomPosition(
    const std::shared_ptr<DmpEvtSimuPrimaries> simu_primaries,
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TH1D &h_BGOrec_topX_vs_realX,
    TH1D &h_BGOrec_topY_vs_realY,
    TH1D &h_real_slopeX,
    TH1D &h_real_slopeY,
    TH1D &h_BGOrec_slopeX,
    TH1D &h_BGOrec_slopeY,
    TH1D &h_real_interceptX,
    TH1D &h_real_interceptY,
    TH1D &h_BGOrec_interceptX,
    TH1D &h_BGOrec_interceptY,
    TH2D &h_real_topMap,
    TH2D &h_BGOreco_topMap,
    TH2D &h_real_bottomMap,
    TH2D &h_BGOreco_bottomMap,
    const double energy_w)
{
    // Get the real simu position
    TVector3 orgPosition;
    orgPosition.SetX(simu_primaries->pv_x);
    orgPosition.SetY(simu_primaries->pv_y);
    orgPosition.SetZ(simu_primaries->pv_z);

    TVector3 orgMomentum;
    orgMomentum.SetX(simu_primaries->pvpart_px);
    orgMomentum.SetY(simu_primaries->pvpart_py);
    orgMomentum.SetZ(simu_primaries->pvpart_pz);

    std::vector<double> slope(2, 0);
    std::vector<double> intercept(2, 0);

    slope[0] = orgMomentum.Z() ? orgMomentum.X() / orgMomentum.Z() : -999;
    slope[1] = orgMomentum.Z() ? orgMomentum.Y() / orgMomentum.Z() : -999;
    intercept[0] = orgPosition.X() - slope[0] * orgPosition.Z();
    intercept[1] = orgPosition.Y() - slope[1] * orgPosition.Z();

    double real_topX = slope[0] * BGO_TopZ + intercept[0];
    double real_topY = slope[1] * BGO_TopZ + intercept[1];

    double real_bottomX = slope[0] * BGO_BottomZ + intercept[0];
    double real_bottomY = slope[1] * BGO_BottomZ + intercept[1];

    // Get the reco position
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    double reco_topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double reco_topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    double reco_bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
    double reco_bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

    // Fill slopes
    h_real_slopeX.Fill(slope[0], energy_w);
    h_real_slopeY.Fill(slope[1], energy_w);
    h_BGOrec_slopeX.Fill(bgoRec_slope[0], energy_w);
    h_BGOrec_slopeY.Fill(bgoRec_slope[1], energy_w);

    // Fill intercepts
    h_real_interceptX.Fill(intercept[0], energy_w);
    h_real_interceptY.Fill(intercept[1], energy_w);
    h_BGOrec_interceptX.Fill(bgoRec_intercept[0], energy_w);
    h_BGOrec_interceptY.Fill(bgoRec_intercept[1], energy_w);

    auto spread_topX = real_topX - reco_topX;
    auto spread_topY = real_topY - reco_topY;

    // Fill spreads
    h_BGOrec_topX_vs_realX.Fill(spread_topX, energy_w);
    h_BGOrec_topY_vs_realY.Fill(spread_topY, energy_w);

    // Fill maps
    h_real_topMap.Fill(real_topX, real_topY, energy_w);
    h_real_bottomMap.Fill(real_bottomX, real_bottomY, energy_w);
    h_BGOreco_topMap.Fill(reco_topX, reco_topY, energy_w);
    h_BGOreco_bottomMap.Fill(reco_bottomX, reco_bottomY, energy_w);
}












