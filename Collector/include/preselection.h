#ifndef PRESELECTION_H
#define PRESELECTION_H

struct p_cuts {
    bool maxelayer_cut {false};
    bool maxbarlayer_cut {false};
    bool bgotrack_cut {false};
    bool bgofiducial_cut {false};
    bool nbarlayer13_cut {false};
    bool maxrms_cut {false};
    bool trackselection_cut {false};
    bool psdstkmatch_cut {false};
    bool psdcharge_cut {false};
    bool stkcharge_cut {false};

    void Reset() {
        maxelayer_cut = false;
        maxbarlayer_cut = false;
        bgotrack_cut = false;
        bgofiducial_cut = false;
        nbarlayer13_cut = false;
        maxrms_cut = false;
        trackselection_cut = false;
        psdstkmatch_cut = false;
        psdcharge_cut = false;
        stkcharge_cut = false;
    }
};

class preselection {
    public:
        preselection() {};
        ~preselection() {};

        void SetMaxELayerCut(const bool cut_status) {preselection_cuts.maxelayer_cut = cut_status;};
        void SetMaxBarLayerCut(const bool cut_status) {preselection_cuts.maxbarlayer_cut = cut_status;};
        void SetBGOTrackCut(const bool cut_status) {preselection_cuts.bgotrack_cut = cut_status;};
        void setBGOFiducialCut(const bool cut_status) {preselection_cuts.bgofiducial_cut = cut_status;};
        void SetNBarLayer13Cut(const bool cut_status) {preselection_cuts.nbarlayer13_cut = cut_status;};
        void SetMaxRMSCut(const bool cut_status) {preselection_cuts.maxrms_cut = cut_status;};
        void SetTrackSelectionCut(const bool cut_status) {preselection_cuts.trackselection_cut = cut_status;};
        void SetPSDSTKMatchCut(const bool cut_status) {preselection_cuts.psdstkmatch_cut = cut_status;};
        void SetPSDChargeCut(const bool cut_status) {preselection_cuts.psdcharge_cut = cut_status;};
        void SetSTKChargeCut(const bool cut_status) {preselection_cuts.stkcharge_cut = cut_status;};
        void Reset();
        const p_cuts GetPreselectionStatus();
    private:
        p_cuts preselection_cuts;
};

#endif