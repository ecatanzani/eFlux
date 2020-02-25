#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

extern void buildAcceptance(
                                const std::string accInputPath,
                                const bool verbose,
                                const std::vector<float> &logEBins
                            );

extern std::shared_ptr < DmpChain > aggregateEventsDmpChain(
                                                                const std::string accInputPath,
                                                                const bool verbose
                                                            );

extern std::shared_ptr < TChain > aggregateEventsTChain(
                                                            const std::string accInputPath,
                                                            const bool verbose
                                                        );

extern bool maxElater_cut(
                            std::shared_ptr < DmpEvtBgoRec > bgorec, 
                            const double egyLayerRatio, 
                            const double bgoTotalE
                        );

extern bool maxBarLayer_cut(
                                std::shared_ptr < DmpEvtBgoHits > bgohits, 
                                const int nBgoHits
                            );

extern bool BGOTrackContainment_cut(
                                        std::shared_ptr < DmpEvtBgoRec > bgorec, 
                                        bool passEvent
                                    );

#endif