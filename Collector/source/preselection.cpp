#include "preselection.h"

void preselection::Reset() {
    preselection_cuts.Reset();
}

const p_cuts preselection::GetPreselectionStatus() {
    return preselection_cuts;
}