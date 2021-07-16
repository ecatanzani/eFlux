# Efficiency MC/DATA computation

Facility to compute the efficiency for DATA/MC events.

## MC events

*h_geometric* -> histogram of the events that pass the geometric cut on the calorimeter

*h_geometric_trigger* -> historgam of the events that pass the geometric cut of the calorimeter and have been triggered

The ratio *h_geometric_trigger* over *h_geometric* is the **trigger efficiency**.

*h_trigger* -> histogram of the events that have been triggered

*h_maxElayer* -> histogram of the events passing the maxElayer cut

*h_maxBarlayer* -> histogram of the events passing the maxBarlayer cut

*h_BGOTrackContainment* -> histogram of the events passing the BGOTrackContainment cut

*h_bgo_fiducial* -> histogram of the events passing the BGO fiducial volume cut

The ratio between each one of those histos over the *h_trigger* represents the efficiency of the cut. 

**Each one of those cuts is independent among the others**.

*h_nbarlayer13* -> histogram of the events passing the nbarlayer13 cut and also contained into the BGO fiducial volume

*h_maxrms* -> histogram of the events passing both the nbarlayer13 and maxRms cuts

*h_track_selection* -> histogram of the events passing both the maxRms and the track selection cuts

*h_psd_stk_match* -> histogram of the events passing both the track selection and PSD/STK match cuts

*h_psd_charge* -> histogram of the events passing both the PSD/STK match and PSD charge cuts

*h_stk_charge* -> histogram of the events passing both the PSD charge and STK charge cuts

*h_all_cuts* -> histogram of the events passing all the preselection cuts

**For those cuts the efficiency needs to be obtained from the ratio with the previous one** 

## DATA events

*h_trigger* -> histogram of the events that have been triggered

*h_maxElayer* -> histogram of the events passing the maxElayer cut

*h_maxBarlayer* -> histogram of the events passing the maxBarlayer cut

*h_BGOTrackContainment* -> histogram of the events passing the BGOTrackContainment cut

*h_bgo_fiducial* -> histogram of the events passing the BGO fiducial volume cut

The ratio between each one of those histos over the *h_trigger* represents the efficiency of the cut. 

**Each one of those cuts is independent among the others**.

*h_nbarlayer13* -> histogram of the events passing the nbarlayer13 cut and also contained into the BGO fiducial volume

*h_maxrms* -> histogram of the events passing both the nbarlayer13 and maxRms cuts

*h_track_selection* -> histogram of the events passing both the maxRms and the track selection cuts

*h_psd_stk_match* -> histogram of the events passing both the track selection and PSD/STK match cuts

*h_psd_charge* -> histogram of the events passing both the PSD/STK match and PSD charge cuts

*h_stk_charge* -> histogram of the events passing both the PSD charge and STK charge cuts

*h_all_cuts* -> histogram of the events passing all the preselection cuts

**For those cuts the efficiency needs to be obtained from the ratio with the previous one** 