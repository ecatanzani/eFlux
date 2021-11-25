#!/usr/bin/python
'''
Created on Mar 5, 2015
@author: andrii
'''

import sys

#@ do not remove - can be used for visualizing MIP events
#import MIPSelection  

printhelp = False
for arg in sys.argv:
    if arg =="-h" or arg =="-H" or arg[1:].lower() == "help" or arg[2:].lower()=="help" or arg.lower()=="help" or arg.lower()=="h": 
        printhelp = True 
        sys.argv.remove(arg)
    

from ROOT import gStyle, gSystem, TH2D, TFile, TClonesArray, TMarker, TLine, TArrow, TEllipse, kRed, kGreen,kBlue, kDashed, kDotted, TCanvas, TH1F, kTRUE, gROOT, TChain, TMath, TPaveText
gSystem.Load("libDmpEvent.so")
from ROOT import DmpEvtBgoHits, DmpEvtSimuPrimaries, DmpEvtBgoRec, DmpEvtHeader, DmpEvtPsdHits, DmpStkHkeepHeader, DmpStkTrack
from array import array
gROOT.SetBatch(kTRUE)

#gStyle.SetPalette(51) 
#gROOT.ForceStyle()

VALIDATION_PLOTS_FILE="validation_plots.root" 

MEV = 1.
GEV = 1000.
TEV = 1000000.

COLORS_FOR_TRACK_MARKERS = [kRed,kGreen,kBlue]
SIZES_FOR_TRACK_MARKERS  = [2,1,0.5]


def deg2rad(deg):
    pi = TMath.Pi()
    return (deg * pi / 180.)

def rad2deg(rad):
    pi = TMath.Pi()
    return (180.*rad/pi)

class validator:
    
    # cuts
    PSD_HIT_ENERGY_CUT     = 1.0
    
    STK_TRACKS_TOP_Z       =   480.
    BGO_DIRECTION_TOP_Z    =   450.
    BGO_DIRECTION_BOTTOM_Z =  -250.
    STK_Z_COORDINATES_X    = [-46.281957580500006, -79.0319957275, -111.38199420165, -144.1319970627, -176.8819970627, -210.0319976349]
    STK_Z_COORDINATES_Y    = [-41.96795840450001, -74.7179927368, -107.4679965515, -140.61799712369998, -173.36799664685, -206.1179980774]
    X_Y_WIDTH              = 605.
    STK_PSD_WIDTH          = 416.0 #252.0 + 164.0
    
    CORR_PLOTS_AMS_MIN   =  0
    CORR_PLOTS_AMS_MAX   =  700
    CORR_PLOTS_AMS_NBINS =  700
    
    CORR_PLOTS_STK_MIN   =  0
    CORR_PLOTS_STK_MAX   =  100
    CORR_PLOTS_STK_NBINS =  1000
    
    
    OCCUP_PLOTS_MIN_STK   = -200
    OCCUP_PLOTS_MAX_STK   =  200
    OCCUP_PLOTS_NBINS_STK = 1000
    
    OCCUP_PLOTS_MIN_AMS   = -200
    OCCUP_PLOTS_MAX_AMS   =  200
    OCCUP_PLOTS_NBINS_AMS = 1000
    
    STK_SIGNAL_HISTO_MIN   = 0
    STK_SIGNAL_HISTO_MAX   = 500
    STK_SIGNAL_HISTO_NBINS = 500
    
    AMS_SIGNAL_HISTO_MIN   = 0
    AMS_SIGNAL_HISTO_MAX   = 500
    AMS_SIGNAL_HISTO_NBINS = 500
    
    STK_POSITION_HISTO_NBINS =  800
    STK_POSITION_HISTO_MIN   = -400
    STK_POSITION_HISTO_MAX   =  400
    
    AMS_POSITION_HISTO_NBINS = 3200
    AMS_POSITION_HISTO_MIN   = 0
    AMS_POSITION_HISTO_MAX   = 640
    
    NCLSUTERS_TO_TRUNCATE_STK = 2 #2
    NCLSUTERS_TO_TRUNCATE_AMS_SSIDE = 1
    NCLSUTERS_TO_TRUNCATE_AMS_KSIDE = 1
    
    STK_SQRTSIGNAL_HISTO_NBINS = 5000
    STK_SQRTSIGNAL_HISTO_MIN   = 0
    STK_SQRTSIGNAL_HISTO_MAX   = 50
    
    AMS_NCLUSTERS_HISTO_NBINS  = 100
    AMS_NCLUSTERS_HISTO_MIN    = 0
    AMS_NCLUSTERS_HISTO_MAX    = 100
    
    STK_NCLUSTERS_HISTO_NBINS  = 100
    STK_NCLUSTERS_HISTO_MIN    = 0
    STK_NCLUSTERS_HISTO_MAX    = 100
    
    
    MIP_MAXPOSITION_SQRT = 54**0.5
    SHOWPRIMARYDIRECTION = False
    
    """
    POSITION_CUT_STK = {    #20150322-054235
                        0 : {"minx": 143.0, "maxx": 145.0, "miny": -9999999, "maxy": 9999999},
                        1 : {"minx": 143.0, "maxx": 145.0, "miny": -9999999, "maxy": 9999999},
                        2 : {"minx": 143.5, "maxx": 145.0, "miny": -9999999, "maxy": 9999999},
                        3 : {"minx": 143.5, "maxx": 145.0, "miny": -9999999, "maxy": 9999999},
                        4 : {"minx": 144.0, "maxx": 145.0, "miny": -9999999, "maxy": 9999999},
                        5 : {"minx": 142.5, "maxx": 144.0, "miny": -9999999, "maxy": 9999999}
                        }
    
    """
    
    POSITION_CUT_STK = {
                        0 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999},
                        1 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999},
                        2 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999},
                        3 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999},
                        4 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999},
                        5 : {"minx": -9999999, "maxx": 9999999, "miny": -9999999, "maxy": 9999999}
                        }
    

    
    
    """
    POSITION_CUT_STK = {  # 20150319-054016
                        0 : {"minx": 137.5, "maxx": 139.5, "miny": -9999999, "maxy": 9999999},
                        1 : {"minx": 138.0, "maxx": 140.0, "miny": -9999999, "maxy": 9999999},
                        2 : {"minx": 138.0, "maxx": 140.0, "miny": -9999999, "maxy": 9999999},
                        3 : {"minx": 138.5, "maxx": 140.0, "miny": -9999999, "maxy": 9999999},
                        4 : {"minx": 138.0, "maxx": 140.5, "miny": -9999999, "maxy": 9999999},
                        5 : {"minx": 138.5, "maxx": 140.5, "miny": -9999999, "maxy": 9999999}
                        }
    """
    
    
    
    #MAX_EVENTS = 10000
    
    
    def __init__(self, amsconffile, datafilename,maxevents, showprimarydirection):
        #self.canvases = []
        self.drawclustersfortrack = True
        self.amsconffile = amsconffile
        self.SHOWPRIMARYDIRECTION = showprimarydirection
        #self.datafile = TFile(datafilename, "READ")
        #self.tree = self.datafile.Get("CollectionTree")
        self.tree = TChain("CollectionTree")
        for f in datafilename:
            print "Using file: "+f
            self.tree.Add(f) 
        self.tree.SetBranchStatus("*",1)

        if self.amsconffile is not None:
            self.__get_ams_parameters__()
        else:
            print "    WARNING: no AMS congig file set"
        
        self.MAX_EVENTS = maxevents
        self.branchnames = map(lambda branch: branch.GetName(), self.tree.GetListOfBranches())
        
        self.amsclustersfound = False
        self.stktracks    = None
        self.bgorec       = None        
        self.amsclusters  = None
        self.stkclusters  = TClonesArray("DmpStkSiCluster")
        self.bgohits      = None
        self.psdhits      = None
        self.eventheader  = None
        self.mcprimaries  = None
        self.truthtrajectory = None
        self.stkbesttrack = None
        
        self.tree.SetBranchAddress("StkClusterCollection",self.stkclusters)
        self.b_stkclusters = self.tree.GetBranch("StkClusterCollection")
        if "DmpEvtBgoHits" in self.branchnames:
            self.bgohits      = DmpEvtBgoHits()
            self.tree.SetBranchAddress("DmpEvtBgoHits", self.bgohits)
            self.b_bgohits  = self.tree.GetBranch("DmpEvtBgoHits")
        else:
            print "    WARNING: no BGO hits collection found"
        
        self.b_mcprimaries = None
        if "DmpEvtSimuPrimaries" in self.branchnames:
            self.mcprimaries  = DmpEvtSimuPrimaries()
            #TClonesArray("DmpEvtMCPrimaryParticle")
            self.tree.SetBranchAddress("DmpEvtSimuPrimaries",self.mcprimaries)
            self.b_mcprimaries = self.tree.GetBranch("DmpEvtSimuPrimaries")
        
        if "DmpTruthTrajectoriesCollection" in self.branchnames:
            self.truthtrajectory = TClonesArray("DmpSimuTrajectory")
            self.tree.SetBranchAddress("DmpTruthTrajectoriesCollection", self.truthtrajectory)
            self.b_truthtrajectory = self.tree.GetBranch("DmpTruthTrajectoriesCollection")
            
        psdname = [name for name in ['DmpPsdHits'] if name in self.branchnames ]
        if psdname:
            psdname = psdname[0]
            self.psdhits  = DmpEvtPsdHits()
            self.tree.SetBranchAddress(psdname, self.psdhits)        
            self.b_psdhits = self.tree.GetBranch(psdname)
        else:
            print "    WARNING: no PSD hits collection found"
        
        if "StkKalmanTracks" in self.branchnames:
            self.stktracks = TClonesArray("DmpStkTrack")
            self.tree.SetBranchAddress("StkKalmanTracks", self.stktracks)
            self.b_stktracks = self.tree.GetBranch("StkKalmanTracks")
            
        else:
            print "    WARNING: no STK tracks collection found"
        
        if "DmpEvtBgoRec" in self.branchnames:
            self.bgorec = DmpEvtBgoRec()
            self.tree.SetBranchAddress("DmpEvtBgoRec", self.bgorec)
            self.b_bgorec = self.tree.GetBranch("DmpEvtBgoRec")
        else:
            print "    WARNING: no BGO shower direction collection found"
            
        if "AmsClusterCollection" in self.branchnames:
            self.amsclusters = TClonesArray("Cluster")
            self.tree.SetBranchAddress("AmsClusterCollection", self.amsclusters)
            self.amsclustersfound = True
        else:
            print "    WARNING: no AMS cluster collection"
        
        if "EventHeader" in  self.branchnames:
            self.eventheader  = DmpEvtHeader()
            self.tree.SetBranchAddress("EventHeader", self.eventheader)
            self.b_eventheader = self.tree.GetBranch("EventHeader")
         
        if self.amsconffile is not None:   
            self.all_ams_z = map(lambda x: self.__get_ams_zposition__(x["ladderid"]), self.amsparameters)
            max_ams_z = max(self.all_ams_z)
            min_ams_z = min(self.all_ams_z)
            self.AMS_TOP_Z    = max_ams_z + (max_ams_z - min_ams_z)* 0.2
            self.AMS_BOTTOM_Z = min_ams_z - (max_ams_z - min_ams_z)* 0.2
        
        if "STKBestTrack" in self.branchnames:
            self.stkbesttrack = DmpStkTrack()
            self.tree.SetBranchAddress("STKBestTrack", self.stkbesttrack)
    
    def __get_ams_parameters__(self):
        self.amsparameters = []
        f = open(self.amsconffile, "r")
        l = f.readlines()
        l = map(lambda line: line.strip(), l)
        for line in l:
            if not line: continue
            if line[0] == "#": continue
            values = line.split()
            self.amsparameters.append({
                                       "ladderid":   int(values[0]), 
                                       "distance": float(values[1]),
                                       "offsetx" : float(values[2]),
                                       "offsety" : float(values[3]),
                                       })
    
            
    def __get_ams_zposition__(self, ladderid):
            ladder = filter(lambda x: x["ladderid"] == ladderid, self.amsparameters)[0]
            return self.__get_ams_z__(ladder)
        
    def __get_ams_z__(self, ladder):
            return - ladder["distance"] - self.STK_PSD_WIDTH
            
    def __is_ams_cluster_x__(self,cluster):
        if cluster.address < 640: 
            return True
        return False
    
    
    def __is_ams_cluster_sside__(self,cluster):
        if cluster.side == 0:
            return True
        return False
    
        
    def __get_ams_coordinate_x__(self,cluster):
        ladder = self.__get_ams_ladder__(cluster)
        return ladder["offsetx"] - cluster.GetCoG() * 0.110
    
    
    def __get_ams_coordinate_y__(self,cluster):
        ladder = self.__get_ams_ladder__(cluster)
        return ladder["offsety"] + (cluster.GetCoG() - 640.) * 0.208
    
    
    def __get_ams_coordinate_z__(self,cluster):
        ladder = self.__get_ams_ladder__(cluster)
        return self.__get_ams_z__(ladder)
    
    
    def __get_ams_cluster_cog_x__(self,cluster):
        return cluster.GetCoG()
    
    def __get_ams_cluster_cog_y__(self,cluster):
        return cluster.GetCoG() - 640.
    
    def __get_ams_ladder__(self, cluster):
        return filter(lambda x: x["ladderid"] == cluster.ladder, self.amsparameters)[0] 
    
    
    def __get_ams_cluster_energy__(self,cluster):
        return sum ([cluster.Signal[i] for i in xrange(cluster.length)])
            
    
    def __get_stk_ladder_centroid__(self,cluster):
        return cluster.getLadderCentroid()
    
    
    def __get_stk_layer__(self,cluster):
        if cluster.isX():
            return self.STK_Z_COORDINATES_X.index(cluster.GetZ())
        else:
            return self.STK_Z_COORDINATES_Y.index(cluster.GetZ())
        
        
    def __get_stk_cluster_energy__(self,cluster):
        return cluster.getEnergy()
        
            
    def __plot_event_diplay__(self, canvasname, event_index):

        
        
        
        
        if self.amsclustersfound:
            c = TCanvas(canvasname, canvasname, 1210, 870)
            c.Divide(2,2)
            DAMPE_X_CANVAS = 1
            DAMPE_Y_CANVAS = 3
        else:
            c = TCanvas(canvasname, canvasname, 1210, 435);
            c.Divide(2,1)
            DAMPE_X_CANVAS = 1
            DAMPE_Y_CANVAS = 2
        

        if self.bgorec:
            en = self.bgorec.GetTotalEnergy()
            if en > TEV:
                en = en / TEV
                en = "   E=%.3f TeV"%(en)
            elif en > GEV:
                en = en / GEV
                en = "   E=%.3f GeV"%(en)
            else:
                en = en / MEV
                en = "   E=%.3f MeV"%(en)
        else:
            en = ""

        
        histo_bgo_x = TH2D("xz", "xz",   44,-self.X_Y_WIDTH,self.X_Y_WIDTH,28,-333,479);
        histo_bgo_y = TH2D("yz","yz",    44,-self.X_Y_WIDTH,self.X_Y_WIDTH,28,-333,479);
        #histo_bgo_box_x = TH2D("DAMPE XZ","DAMPE XZ"+en,44,-self.X_Y_WIDTH,self.X_Y_WIDTH,30,-391,479);
        #histo_bgo_box_y = TH2D("DAMPE YZ","DAMPE YZ"+en,44,-self.X_Y_WIDTH,self.X_Y_WIDTH,30,-391,479);


        histo_bgo_box_x = TH2D("DAMPE XZ","DAMPE XZ"+en,44,-self.X_Y_WIDTH,self.X_Y_WIDTH,30,-391,479);
        histo_bgo_box_y = TH2D("DAMPE YZ","DAMPE YZ"+en,44,-self.X_Y_WIDTH,self.X_Y_WIDTH,30,-391,479);
        
        """
        histo_psd_x_1 = TH2D("XZ", "XZ", 21, -420, 420, 1, -291.5, -277.5)
        histo_psd_x_2 = TH2D("XZ", "XZ", 20, -400, 400, 1, -305.5, -291.5)
        histo_psd_y_1 = TH2D("YZ", "YZ", 21, -420, 420, 1, -317.7, -303.7)
        histo_psd_y_2 = TH2D("YZ", "YZ", 20, -400, 400, 1, -331.7, -317.7)
        histo_psd_box_x_1 = TH2D("XZ", "XZ", 21, -420, 420, 1, -291.5, -277.5)
        histo_psd_box_x_2 = TH2D("XZ", "XZ", 20, -400, 400, 1, -305.5, -291.5)
        histo_psd_box_y_1 = TH2D("YZ", "YZ", 21, -420, 420, 1, -317.7, -303.7)
        histo_psd_box_y_2 = TH2D("YZ", "YZ", 20, -400, 400, 1, -331.7, -317.7)
        """
        bins1 = array("d", [-414 + 28 * ((i+1)/2) + 12 * (i/2) for i in xrange(42)])
        bins2 = array("d", [-394 + 28 * ((i+1)/2) + 12 * (i/2) for i in xrange(40)])
        histo_psd_x_1 = TH2D(canvasname+"XZ_psd_1", "XZ", 41, bins1, 1, -291.5, -277.5)
        histo_psd_x_2 = TH2D(canvasname+"XZ_psd_2", "XZ", 39, bins2, 1, -305.5, -291.5)
        histo_psd_y_1 = TH2D(canvasname+"YZ_psd_1", "YZ", 41, bins1, 1, -317.7, -303.7)
        histo_psd_y_2 = TH2D(canvasname+"YZ_psd_2", "YZ", 39, bins2, 1, -331.7, -317.7)
        histo_psd_box_x_1 = TH2D(canvasname+"XZ_bgo_1", "XZ", 41, bins1, 1, -291.5, -277.5)
        histo_psd_box_x_2 = TH2D(canvasname+"XZ_bgo_2", "XZ", 39, bins2, 1, -305.5, -291.5)
        histo_psd_box_y_1 = TH2D(canvasname+"YZ_bgo_1", "YZ", 41, bins1, 1, -317.7, -303.7)
        histo_psd_box_y_2 = TH2D(canvasname+"YZ_bgo_2", "YZ", 39, bins2, 1, -331.7, -317.7)
        
        
        
                
        histo_bgo_box_x.SetStats(False)
        histo_bgo_box_y.SetStats(False)
        histo_bgo_x.SetStats(False)
        histo_bgo_y.SetStats(False)
        
        histo_psd_box_x_1.SetStats(False)
        histo_psd_box_x_2.SetStats(False)
        histo_psd_box_y_1.SetStats(False)
        histo_psd_box_y_2.SetStats(False)
        
        for x in xrange(12, histo_bgo_x.GetNbinsX()-10):
            for y in xrange(17, histo_bgo_x.GetNbinsY()+2, 2):
                histo_bgo_box_x.SetBinContent(x,y,0.001)
                histo_bgo_box_y.SetBinContent(x,y-1,0.001)
        """
        for i in xrange(1,histo_psd_box_x_1.GetNbinsX()+1):
            histo_psd_box_x_1.SetBinContent(i,1,0.001)
            histo_psd_box_y_1.SetBinContent(i,1,0.001)
            if i < histo_psd_box_x_2.GetNbinsX()+1:
                histo_psd_box_x_2.SetBinContent(i,1,0.001)
                histo_psd_box_y_2.SetBinContent(i,1,0.001)
        """
        for i in xrange(1,histo_psd_box_x_1.GetNbinsX()+1,2):
            histo_psd_box_x_1.SetBinContent(i,1,0.001)
            histo_psd_box_y_1.SetBinContent(i,1,0.001)
        
        for i in xrange(1,histo_psd_box_x_2.GetNbinsX()+1,2):
            histo_psd_box_x_2.SetBinContent(i,1,0.001)
            histo_psd_box_y_2.SetBinContent(i,1,0.001)
        
        
        
        if self.bgohits:
            for i in xrange(len(self.bgohits.fEnergy)):
                if self.bgohits.GetLayerID(i)%2:
                    histo_bgo_x.Fill(self.bgohits.GetHitX(i),self.bgohits.GetHitZ(i),self.bgohits.fEnergy[i])
                else:
                    histo_bgo_y.Fill(self.bgohits.GetHitY(i),self.bgohits.GetHitZ(i),self.bgohits.fEnergy[i])
        
        stkpointsx = []
        stkpointsy = []
        for i in xrange(self.stkclusters.GetLast()+1):
            cluster = self.stkclusters.ConstructedAt(i)
            if cluster.isX():
                stkpointsx.append(TMarker(cluster.GetX(),cluster.GetZ(),3))
            else:
                stkpointsy.append(TMarker( cluster.GetY(), cluster.GetZ(),3 ))
                
        amspointsx = []
        amspointsy = []

        if self.amsclustersfound and self.amsconffile:
            for i in xrange(self.amsclusters.GetLast()+1):
                cluster = self.amsclusters.ConstructedAt(i)
                if self.__is_ams_cluster_x__(cluster):
                    amspointsx.append(TMarker(self.__get_ams_coordinate_x__(cluster),self.__get_ams_coordinate_z__(cluster),3))
                else:
                    amspointsy.append(TMarker(self.__get_ams_coordinate_y__(cluster),self.__get_ams_coordinate_z__(cluster),3))
                
        
        stktracksx = []
        stktracksy = []
        hitsfortrackx = []
        hitsfortracky = []
        if self.stktracks:
            for i in xrange(self.stktracks.GetLast()+1):
                track = self.stktracks.ConstructedAt(i)
                x = track.getImpactPoint().x()
                y = track.getImpactPoint().y()
                z = track.getImpactPoint().z()
                incl_x = track.getTrackParams().getSlopeX()
                incl_y = track.getTrackParams().getSlopeY()
                stktracksx.append(TLine(x , z, incl_x * (self.STK_TRACKS_TOP_Z - z) +x, self.STK_TRACKS_TOP_Z))
                stktracksy.append(TLine(y , z, incl_y * (self.STK_TRACKS_TOP_Z - z) +y, self.STK_TRACKS_TOP_Z))

                # clusters for track
                if i >= len(COLORS_FOR_TRACK_MARKERS): continue
                hitsfortrackx.append([])
                hitsfortracky.append([])
                for p in xrange(track.GetNPoints()):
                    track_clx = track.GetClusterX(p,self.stkclusters)
                    if not track_clx: continue
                    clustertmp = TMarker(track_clx.GetX(),track_clx.GetZ(),3)
                    clustertmp.SetMarkerColor(COLORS_FOR_TRACK_MARKERS[i])
                    clustertmp.SetMarkerSize(SIZES_FOR_TRACK_MARKERS[i])
                    hitsfortrackx[-1].append(clustertmp)
                for p in xrange(track.GetNPoints()):
                    track_cly = track.GetClusterY(p,self.stkclusters)
                    if not track_cly: continue
                    clustertmp = TMarker(track_cly.GetY(),track_cly.GetZ(),3)
                    clustertmp.SetMarkerColor(COLORS_FOR_TRACK_MARKERS[i])
                    clustertmp.SetMarkerSize(SIZES_FOR_TRACK_MARKERS[i])
                    hitsfortracky[-1].append(clustertmp)
                
        #### Read the external file for best tracks
        x = self.stkbesttrack.getImpactPoint().x()
        y = self.stkbesttrack.getImpactPoint().y()
        z = self.stkbesttrack.getImpactPoint().z()
        incl_x = self.stkbesttrack.getTrackParams().getSlopeX()
        incl_y = self.stkbesttrack.getTrackParams().getSlopeY()
        stktracksx.append(TLine(x , z, incl_x * (self.STK_TRACKS_TOP_Z - z) +x, self.STK_TRACKS_TOP_Z))
        stktracksy.append(TLine(y , z, incl_y * (self.STK_TRACKS_TOP_Z - z) +y, self.STK_TRACKS_TOP_Z))
        stktracksx[-1].SetLineColor(kRed)
        stktracksy[-1].SetLineColor(kRed)
        stktracksx[-1].SetLineWidth(2)
        stktracksy[-1].SetLineWidth(2)
                
        #### add true direction here #####        
        trueDirection_x = None
        trueDirection_y = None
        trueEnergy = None
        if self.mcprimaries:
            x = self.mcprimaries.pv_x
            y = self.mcprimaries.pv_y
            z = self.mcprimaries.pv_z

            px = self.mcprimaries.pvpart_px
            py = self.mcprimaries.pvpart_py
            pz = self.mcprimaries.pvpart_pz

            trueTheta =  self.mcprimaries.pv_theta
            truePhi   =  self.mcprimaries.pv_phi
            vertexRadius = self.mcprimaries.pv_r # vertexRadius
            trueEnergy = self.mcprimaries.pvpart_ekin # in MeV
            trueDirection_x = TMarker(x, z, 15)
            trueDirection_y = TMarker(y, z, 15)
            trueDirection_x.SetMarkerColor(kGreen)
            trueDirection_y.SetMarkerColor(kGreen)
            trueDirection_x.SetMarkerStyle(34)
            trueDirection_y.SetMarkerStyle(34)
            
            #trueDirection_x = TLine(x , z, (px+x), (pz+z))
            #trueDirection_y = TLine(y , z, (py+y), (pz+z))
            if not pz: pz = 0.0000001
            trueDirection_x = TLine(x, z, x + (self.STK_TRACKS_TOP_Z - z)*px/pz, self.STK_TRACKS_TOP_Z)
            trueDirection_y = TLine(y, z, y + (self.STK_TRACKS_TOP_Z - z)*py/pz, self.STK_TRACKS_TOP_Z)
        
        true_trajectoryx = []
        true_trajectoryy = []
        if self.truthtrajectory:
            for i in xrange(self.truthtrajectory.GetLast()+1):
                traj = self.truthtrajectory.ConstructedAt(i)
                start_x = traj.start_x
                start_y = traj.start_y
                start_z = traj.start_z
                stop_x = traj.stop_x
                stop_y = traj.stop_y
                stop_z = traj.stop_z
                if i == 0:
                    true_trajectoryx.append(TMarker(start_x, start_z,4))
                    true_trajectoryy.append(TMarker(start_y, start_z,4))
                    true_trajectoryx[-1].SetMarkerColor(kRed)
                    true_trajectoryy[-1].SetMarkerColor(kRed)
                    true_trajectoryx[-1].SetMarkerStyle(34)
                    true_trajectoryy[-1].SetMarkerStyle(34)
                true_trajectoryx.append(TLine(start_x, start_z, stop_x, stop_z))
                true_trajectoryy.append(TLine(start_y, start_z, stop_y, stop_z))
                
        bgodirection_x = None
        bgodirection_y = None
        if self.bgorec:
            x = self.bgorec.GetTrajectoryLocation2D().x()
            y = self.bgorec.GetTrajectoryLocation2D().y()
            z = self.bgorec.GetTrajectoryLocation2D().z()
            incl_x = self.bgorec.GetSlopeXZ()
            incl_y = self.bgorec.GetSlopeYZ()
            bgodirection_x = TLine(x , z, incl_x * (self.BGO_DIRECTION_TOP_Z - z) +x, self.BGO_DIRECTION_TOP_Z)
            #print "bgoshower_x = TLine(",x," , ",z,", ",incl_x*(self.BGO_DIRECTION_TOP_Z - z) +x,",", self.BGO_DIRECTION_TOP_Z,")"
            #print "xz(200) = ", incl_x*(200-z+1./incl_x*x)
            bgodirection_y = TLine(y , z, incl_y * (self.BGO_DIRECTION_TOP_Z - z) +y, self.BGO_DIRECTION_TOP_Z)
            #print "bgoshower_y = TLine(",y," , ",z,", ",incl_y*(self.BGO_DIRECTION_TOP_Z - z) +y,",", self.BGO_DIRECTION_TOP_Z,")"
            #print "yz(200) = ", incl_y*(200-z+1./incl_y*y)
            chi2_x = self.bgorec.GetChi2XZ()
            chi2_y = self.bgorec.GetChi2YZ()
            bgodirection_x.SetLineWidth(2)
            bgodirection_y.SetLineWidth(2)
            bgodirection_x.SetLineStyle(2)
            bgodirection_y.SetLineStyle(2)
            bgodirection_x.SetLineColor(8);
            bgodirection_y.SetLineColor(8);
        # NOT AVAILABLE IN RELEASE 5-1-2
        """
        if self.bgodirection:
            bgodirection_y = TLine(
                                   self.bgodirection.fTrajectoryLocation.y() + (self.BGO_DIRECTION_TOP_Z  -self.bgodirection.fTrajectoryLocation.z()) * self.bgodirection.fTrajectoryDirection.y() / self.bgodirection.fTrajectoryDirection.z(),
                                   self.BGO_DIRECTION_TOP_Z,
                                   self.bgodirection.fTrajectoryLocation.y() + (self.BGO_DIRECTION_BOTTOM_Z - self.bgodirection.fTrajectoryLocation.z()) * self.bgodirection.fTrajectoryDirection.y() / self.bgodirection.fTrajectoryDirection.z(),
                                   self.BGO_DIRECTION_TOP_Z
                                   )
            bgodirection_x = TLine(
                                   self.bgodirection.fTrajectoryLocation.x() + (self.BGO_DIRECTION_TOP_Z -self.bgodirection.fTrajectoryLocation.z()) * self.bgodirection.fTrajectoryDirection.x()/self.bgodirection.fTrajectoryDirection.z(),
                                   self.BGO_DIRECTION_TOP_Z,
                                   self.bgodirection.fTrajectoryLocation.x() + (self.BGO_DIRECTION_BOTTOM_Z -self.bgodirection.fTrajectoryLocation.z()) * self.bgodirection.fTrajectoryDirection.x()/self.bgodirection.fTrajectoryDirection.z(),
                                   self.BGO_DIRECTION_BOTTOM_Z)
            bgodirection_y.SetLineColor(kGreen)
            bgodirection_x.SetLineColor(kGreen)
        
        """ 

        """
        psdhitsx = []
        psdhitsy = []
        for i in xrange(len(self.psdhits.fEnergy)):
            if self.psdhits.fGlobalBarID[i] > 3000:
                psdhitsx.append(TMarker(self.psdhits.GetHitX(i),self.psdhits.GetHitZ(i),3))
                psdhitsx[-1].SetMarkerColor(kRed)
            else:
                psdhitsy.append(TMarker(self.psdhits.GetHitY(i),self.psdhits.GetHitZ(i),3))
                psdhitsy[-1].SetMarkerColor(kRed)
        """
        

        if self.psdhits:
            for i in xrange(len(self.psdhits.fEnergy)):
                if self.psdhits.fEnergy[i] < self.PSD_HIT_ENERGY_CUT: continue
                if self.psdhits.fGlobalBarID[i] > 3000:
                    histo_psd_x_1.Fill(self.psdhits.GetHitX(i),self.psdhits.GetHitZ(i),self.psdhits.fEnergy[i])
                    histo_psd_x_2.Fill(self.psdhits.GetHitX(i),self.psdhits.GetHitZ(i),self.psdhits.fEnergy[i])
                else:
                    histo_psd_y_1.Fill(self.psdhits.GetHitY(i),self.psdhits.GetHitZ(i),self.psdhits.fEnergy[i])
                    histo_psd_y_2.Fill(self.psdhits.GetHitY(i),self.psdhits.GetHitZ(i),self.psdhits.fEnergy[i])
                
                
                
        linesx = map(lambda z: TLine(-400, z, 400, z), self.STK_Z_COORDINATES_X)
        linesy = map(lambda z: TLine(-400, z, 400, z), self.STK_Z_COORDINATES_Y)
        
        if self.amsconffile:
            linesxamsx = map(lambda z: TLine(-400, z, 400, z), self.all_ams_z)
            linesxamsy = map(lambda z: TLine(-400, z, 400, z), self.all_ams_z)
        
        if self.mcprimaries:
            paveText = TPaveText(.706954,.522399,.891487,.880786,"NDC")
            paveText.AddLine(.0,.75,1.,.75)
            paveText.AddText("MC truth information")
            paveText.AddText("E_{true} = %1.2f GeV"%float(trueEnergy/1.e3))
            paveText.AddText("R_{vertex} = %1.2f m"%float(vertexRadius/1.e3))
            paveText.AddText("#theta = %1.2f #phi = %1.2f"%(trueTheta,truePhi))
            if self.truthtrajectory:
                paveText.AddText("# trajectories %i"%(len(true_trajectoryx)-1))
        #gStyle.SetPalette(51) 
        c.cd(DAMPE_X_CANVAS)
        #gStyle.SetPalette(51) # kDeepSea
        histo_bgo_box_x.Draw("box")
        histo_psd_box_x_1.Draw("box same")
        histo_psd_box_x_2.Draw("box same")
        #histo_bgo_x.Draw("box same")
        histo_bgo_x.Draw("colz same")

        c.Update()
        palette_bgo_x = histo_bgo_x.GetListOfFunctions().FindObject("palette")
        palx1 = 0.90
        palx2 = 0.93
        paly1 = 0.50
        paly2 = 0.90
        try:
            palette_bgo_x.SetX1NDC(palx1)
            palette_bgo_x.SetX2NDC(palx2)
            palette_bgo_x.SetY1NDC(paly1)
            palette_bgo_x.SetY2NDC(paly2)
        except:
            print "failde to change z palette for bgo"
        
        maxpsdz = max(self.psdhits.fEnergy) if len(self.psdhits.fEnergy) else self.PSD_HIT_ENERGY_CUT

        histo_psd_x_1.GetZaxis().SetRangeUser(0,maxpsdz)
        histo_psd_x_2.GetZaxis().SetRangeUser(0,maxpsdz)
        histo_psd_x_1.Draw("col same")
        histo_psd_x_2.Draw("colz same")

        c.Update()
        palette_psd_x = histo_psd_x_2.GetListOfFunctions().FindObject("palette")
        paly3 = 0.10
        paly4 = 0.30
        try:
            palette_psd_x.SetX1NDC(palx1)
            palette_psd_x.SetX2NDC(palx2)
            palette_psd_x.SetY1NDC(paly3)
            palette_psd_x.SetY2NDC(paly4)
        except:
            print "failde to change z palette for psd"
        #palette.SetLabelSize(0.01)

        #map(lambda a: a.Draw("same"), psdhitsx)
        map(lambda a: a.Draw("same"), stktracksx)
        map(lambda a: a.Draw("same"), stkpointsx)
        map(lambda a: a.Draw("same"), linesx)
        
        map(lambda a: a.SetLineColor(kRed) if isinstance(a,TLine) else a.SetMarkerColor(kRed), true_trajectoryx)
        map(lambda a: a.SetLineStyle(kDashed) if isinstance(a,TLine) else a.SetMarkerStyle(34), true_trajectoryx)
        map(lambda a: a.Draw("same"), true_trajectoryx)

        if self.drawclustersfortrack:
            [p.Draw("same") for hitsfortrack in hitsfortrackx for p in hitsfortrack]
        
        if bgodirection_x:
            bgodirection_x.Draw("same")
        if trueDirection_x:
            trueDirection_x.SetLineColor(kBlue)
            trueDirection_x.SetLineStyle(kDashed)
            if self.SHOWPRIMARYDIRECTION: trueDirection_x.Draw("same")
        if self.mcprimaries:
            paveText.Draw("same")
            #ell_x = TEllipse(0.,0.,vertexRadius,vertexRadius)
            #ell_x.SetLineColor(kRed)
            #ell_x.SetLineStyle(kDotted)
            #ell_x.SetFillStyle(0)
            #ell_x.Draw("lsame")

        c.cd(DAMPE_Y_CANVAS)
        #gStyle.SetPalette(51) # kDeepSea
        histo_bgo_box_y.Draw("box")
        histo_psd_box_y_1.Draw("box same")
        histo_psd_box_y_2.Draw("box same")
        histo_bgo_y.SetMarkerColor(kRed)
        #histo_bgo_y.Draw("box same")
        histo_bgo_y.Draw("colz same")

        c.Update()
        palette_bgo_y = histo_bgo_y.GetListOfFunctions().FindObject("palette")
        try:
            palette_bgo_y.SetX1NDC(palx1)
            palette_bgo_y.SetX2NDC(palx2)
            palette_bgo_y.SetY1NDC(paly1)
            palette_bgo_y.SetY2NDC(paly2)
        except:
            print "failde to change z palette for bgo"


        histo_psd_y_1.GetZaxis().SetRangeUser(0,maxpsdz)
        histo_psd_y_2.GetZaxis().SetRangeUser(0,maxpsdz)
        histo_psd_y_1.Draw("col same")
        histo_psd_y_2.Draw("colz same")

        c.Update()
        palette_psd_y = histo_psd_y_2.GetListOfFunctions().FindObject("palette")
        try:
            palette_psd_y.SetX1NDC(palx1)
            palette_psd_y.SetX2NDC(palx2)
            palette_psd_y.SetY1NDC(paly3)
            palette_psd_y.SetY2NDC(paly4)
        except:
            print "failde to change z palette for psd"

        #map(lambda a: a.Draw("same"), psdhitsy)
        map(lambda a: a.Draw("same"), stktracksy)
        map(lambda a: a.Draw("same"), stkpointsy)
        map(lambda a: a.Draw("same"), linesy) ### STK tracks
        
        map(lambda a: a.SetLineColor(kRed) if isinstance(a,TLine) else a.SetMarkerColor(kRed), true_trajectoryy)
        map(lambda a: a.SetLineStyle(kDashed) if isinstance(a,TLine) else a.SetMarkerStyle(34), true_trajectoryy)

        map(lambda a: a.Draw("same"), true_trajectoryy)

        #for hitsfortrack in hitsfortracky:
        #    for p in hitsfortrack:
        #        print "test"
        #        print "  ",p.GetX()
        #        print "  ",p.GetY()
        #        p.Draw("same")

        if self.drawclustersfortrack:
            [p.Draw("same") for hitsfortrack in hitsfortracky for p in hitsfortrack ]

        if bgodirection_y:
            bgodirection_y.Draw("same")
        if trueDirection_y:
            trueDirection_y.SetLineColor(kBlue)
            trueDirection_y.SetLineStyle(kDashed)
            if self.SHOWPRIMARYDIRECTION: trueDirection_y.Draw("same")
        if self.amsconffile:    
            c.cd(2)
            histoamsframex = TH2D("AMS XZ","AMS XZ", 10, - self.X_Y_WIDTH, self.X_Y_WIDTH, 10, self.AMS_BOTTOM_Z, self.AMS_TOP_Z )
            histoamsframex.SetStats(False)
            histoamsframex.Draw()
            map(lambda a: a.Draw("same"), linesxamsx)
            map(lambda a: a.Draw("same"), amspointsx)
                
            c.cd(4)
            histoamsframey = TH2D("AMS YZ","AMS YZ", 10, - self.X_Y_WIDTH, self.X_Y_WIDTH, 10, self.AMS_BOTTOM_Z, self.AMS_TOP_Z )
            histoamsframey.SetStats(False)
            histoamsframey.Draw()
            map(lambda a: a.Draw("same"), linesxamsy)
            map(lambda a: a.Draw("same"), amspointsy)
        if self.mcprimaries:
            paveText.Draw("same")
            #ell_y = TEllipse(0.,0.,vertexRadius,vertexRadius)
            #ell_y.SetLineColor(kRed)
            #ell_y.SetLineStyle(kDotted)
            #ell_y.SetFillStyle(0)
            #ell_y.Draw("lsame")

        #histo_bgo_y.Draw("colz")
        #a=raw_input()
       
        #self.canvases.append(c)
        #c.Draw()
        #c.SaveAs("test3.png")
        #a=raw_input()
        c.Write()
        #c.SaveAs("event_%d.pdf"%self.current_entry_id)
        c.SaveAs("event.pdf")
        
        
    def __get_random_stk_coordinate__(self, layer, xdirection = True):
        for i in xrange(self.stkclusters.GetLast()+1):
            cluster = self.stkclusters.ConstructedAt(i)
            if xdirection:
                if not cluster.isX(): continue
            else:
                if not cluster.isY(): continue
            curlayer = self.__get_stk_layer__(cluster)
            if curlayer != layer: continue
            
            #return coordinate
            return self.__get_stk_ladder_centroid__(cluster)
        return None

    
    def __get_random_ams_coordinate__(self, ladder, xdirection = True):
        for i in xrange(self.amsclusters.GetLast()+1):
            cluster = self.amsclusters.ConstructedAt(i)
            if cluster.ladder != ladder: continue
            if xdirection:
                if not self.__is_ams_cluster_x__(cluster): continue
                coordinate =  self.__get_ams_cluster_cog_x__(cluster)
            else:
                if self.__is_ams_cluster_x__(cluster): continue
                coordinate =  self.__get_ams_cluster_cog_y__(cluster)
            return coordinate
        return None

    
        
    def __plot_stk_ams_correlations__(self, topdir, amsladderstoplot = [0,1,4], stklayersplot = [3,4,5]):
        cx = TCanvas("Plots_x", "Plots_x", 1600, 1200)
        cy = TCanvas("Plots_y", "Plots_y", 1600, 1200)
        cx.Divide(len(amsladderstoplot),len(stklayersplot))
        cy.Divide(len(amsladderstoplot),len(stklayersplot))
        histosx = dict( ("stk%d_ams%d"%(i,j), TH2D("X_stk%d_ams%d"%(i,j),"X_stk%d_ams%d"%(i,j), 
                                                   self.CORR_PLOTS_STK_NBINS, 
                                                   self.CORR_PLOTS_STK_MIN, 
                                                   self.CORR_PLOTS_STK_MAX, 
                                                   self.CORR_PLOTS_AMS_NBINS, 
                                                   self.CORR_PLOTS_AMS_MIN, 
                                                   self.CORR_PLOTS_AMS_MAX)
                         )  for i  in  stklayersplot for j in amsladderstoplot)
        histosy = dict( ("stk%d_ams%d"%(i,j), TH2D("Y_stk%d_ams%d"%(i,j),"Y_stk%d_ams%d"%(i,j), 
                                                   self.CORR_PLOTS_STK_NBINS, 
                                                   self.CORR_PLOTS_STK_MIN, 
                                                   self.CORR_PLOTS_STK_MAX, 
                                                   self.CORR_PLOTS_AMS_NBINS, 
                                                   self.CORR_PLOTS_AMS_MIN, 
                                                   self.CORR_PLOTS_AMS_MAX)
                         )  for i  in  stklayersplot for j in amsladderstoplot)
        
        #test_entry = 0
        for entry in xrange(min(self.tree.GetEntries(),self.MAX_EVENTS)):
            
            #test_entry+=1
            #print "test_entry = ",test_entry
            #if self.stkclusters.GetLast()+1> 100: continue
            
            self.tree.GetEntry(entry)
            for i in stklayersplot:
                for j in amsladderstoplot:
                    stk_x_coordinate = self.__get_random_stk_coordinate__(i, xdirection = True)
                    ams_x_coordinate = self.__get_random_ams_coordinate__(j, xdirection = True)
                    stk_y_coordinate = self.__get_random_stk_coordinate__(i, xdirection = False)
                    ams_y_coordinate = self.__get_random_ams_coordinate__(j, xdirection = False)
                    
                    """
                    print "\n "
                    print "stk_x_coordinate: ", stk_x_coordinate
                    print "stk_y_coordinate: ", stk_y_coordinate
                    print "ams_x_coordinate: ", ams_x_coordinate
                    print "ams_y_coordinate: ", ams_y_coordinate
                    """
                    
                    if stk_x_coordinate and ams_x_coordinate:
                        histosx["stk%d_ams%d"%(i,j)].Fill(stk_x_coordinate, ams_x_coordinate)
                    if stk_y_coordinate and ams_y_coordinate:
                        histosy["stk%d_ams%d"%(i,j)].Fill(stk_y_coordinate, ams_y_coordinate)
        for i in xrange(len(stklayersplot)):
            for j in xrange(len(amsladderstoplot)):
                padnumber = 1+i*len(amsladderstoplot)+j 
                cx.cd(padnumber)
                histoname = "stk%d_ams%d"%(stklayersplot[i],amsladderstoplot[j]) 
                histosx[histoname].Draw()
                cy.cd(padnumber)
                histosy[histoname].Draw()
        
        topdir.cd()
        cx.Write()
        cy.Write()
        
        
    def __plot_signal_all_ladders__(self, plotdir, nstkclusterslayer):
        cx_stk = TCanvas("Signal_STK_xladders", "Signal_STK_xladders", 1600, 1200)
        cy_stk = TCanvas("Signal_STK_yladders", "Signal_STK_yladders", 1600, 1200)
        cx_stk.Divide(len(self.STK_Z_COORDINATES_X)/2,2)
        cy_stk.Divide(len(self.STK_Z_COORDINATES_Y)/2,2)
        
        cx_stk_position = TCanvas("Position_STK_xladders", "Position_STK_xladders", 1600, 1200)
        cy_stk_position = TCanvas("Position_STK_yladders", "Position_STK_yladders", 1600, 1200)
        cx_stk_position.Divide(len(self.STK_Z_COORDINATES_X)/2,2)
        cy_stk_position.Divide(len(self.STK_Z_COORDINATES_Y)/2,2)
        
        cx_stk_nclusters = TCanvas("Nclusters_STK_xladders", "Nclusters_STK_xladders", 1600, 1200)
        cy_stk_nclusters = TCanvas("Nclusters_STK_yladders", "Nclusters_STK_yladders", 1600, 1200)
        cx_stk_nclusters.Divide(len(self.STK_Z_COORDINATES_X)/2,2)
        cy_stk_nclusters.Divide(len(self.STK_Z_COORDINATES_Y)/2,2)
        
        cx_sqrt_stk = TCanvas("SqrtSignal_STK_xladders", "SqrtSignal_STK_xladders", 1600, 1200)
        cy_sqrt_stk = TCanvas("SqrtSignal_STK_yladders", "SqrtSignal_STK_yladders", 1600, 1200)
        cx_sqrt_stk.Divide(len(self.STK_Z_COORDINATES_X)/2,2)
        cy_sqrt_stk.Divide(len(self.STK_Z_COORDINATES_Y)/2,2)
        
        
        histosstksignal_all = TH1F("histosstksignal_all", "histosstksignal_all", self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX)
        #histosstksignal_all = TH1F("histosstksignal_all_nstrips", "histosstksignal_all", self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX)
                                   
        
        histosstksignal_all_2 = TH1F("stk_sqrtsignal_all", "stk_sqrtsignal_all", self.STK_SQRTSIGNAL_HISTO_NBINS, self.STK_SQRTSIGNAL_HISTO_MIN, self.STK_SQRTSIGNAL_HISTO_MAX)

        
        
        histosstksignal_x = dict( (i, TH1F("stk_signal_xside_layer%d"%i,
                                           "stk_signal_xside_layer%d"%i,
                                           self.STK_SIGNAL_HISTO_NBINS,
                                           self.STK_SIGNAL_HISTO_MIN,
                                           self.STK_SIGNAL_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_X)) )
        
        histosstksignal_x_2 = dict( (i, TH1F("stk_sqrtsignal_xside_layer%d"%i,
                                             "stk_sqrtsignal_xside_layer%d"%i,
                                             self.STK_SQRTSIGNAL_HISTO_NBINS,
                                             self.STK_SQRTSIGNAL_HISTO_MIN,
                                             self.STK_SQRTSIGNAL_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_X)) )

        
        
        histosstkposition_x = dict( (i, TH1F("stk_position_xside_layer%d"%i,
                                           "stk_position_xside_layer%d"%i,
                                           self.STK_POSITION_HISTO_NBINS,
                                           self.STK_POSITION_HISTO_MIN,
                                           self.STK_POSITION_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_X)) )
        
        
        histosstknclusters_x = dict( (i, TH1F("stk_nclusters_xside_layer%d"%i,
                                              "stk_nclusters_xside_layer%d"%i,
                                              self.STK_NCLUSTERS_HISTO_NBINS,
                                              self.STK_NCLUSTERS_HISTO_MIN,
                                              self.STK_NCLUSTERS_HISTO_MAX)
                                      ) for i in xrange(len(self.STK_Z_COORDINATES_X)) )
        
        histosstknclusters_y = dict( (i, TH1F("stk_nclusters_yside_layer%d"%i,
                                              "stk_nclusters_yside_layer%d"%i,
                                              self.STK_NCLUSTERS_HISTO_NBINS,
                                              self.STK_NCLUSTERS_HISTO_MIN,
                                              self.STK_NCLUSTERS_HISTO_MAX)
                                      ) for i in xrange(len(self.STK_Z_COORDINATES_Y)) )
        
        
        histosstksignal_y = dict( (i, TH1F("stk_signal_xside_layer%d"%i,
                                           "stk_signal_xside_layer%d"%i,
                                           self.STK_SIGNAL_HISTO_NBINS,
                                           self.STK_SIGNAL_HISTO_MIN,
                                           self.STK_SIGNAL_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_Y)) )
        
        histosstkposition_y = dict( (i, TH1F("stk_position_yside_layer%d"%i,
                                           "stk_position_yside_layer%d"%i,
                                           self.STK_POSITION_HISTO_NBINS,
                                           self.STK_POSITION_HISTO_MIN,
                                           self.STK_POSITION_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_X)) )
        
        histosstksignal_y_2 = dict( (i, TH1F("stk_sqrtsignal_yside_layer%d"%i,
                                             "stk_sqrtsignal_yside_layer%d"%i,
                                             self.STK_SQRTSIGNAL_HISTO_NBINS,
                                             self.STK_SQRTSIGNAL_HISTO_MIN,
                                             self.STK_SQRTSIGNAL_HISTO_MAX)
                                   ) for i in xrange(len(self.STK_Z_COORDINATES_Y)) )

        
        
        hisosallstksignals = {}
        histotruncatedstk = TH1F("truncated_stk_energy_allladders_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                 "truncated_stk_energy_allladders_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                 self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX*2)
        
        histotruncatedstk_sqrtsignal_x = TH1F("truncated_stk_sqrtsignal_allladdersx_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                              "truncated_stk_sqrtsignal_allladdersx_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                              self.STK_SQRTSIGNAL_HISTO_NBINS, self.STK_SQRTSIGNAL_HISTO_MIN, self.STK_SQRTSIGNAL_HISTO_MAX)
        
        histotruncatedstk_sqrtsignal_y = TH1F("truncated_stk_sqrtsignal_allladdersy_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                              "truncated_stk_sqrtsignal_allladdersy_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                              self.STK_SQRTSIGNAL_HISTO_NBINS, self.STK_SQRTSIGNAL_HISTO_MIN, self.STK_SQRTSIGNAL_HISTO_MAX)
        
        histotruncatedstk_signal_x = TH1F("truncated_stk_signal_allladdersx_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                          "truncated_stk_signal_allladdersx_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                          self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX)
        
        
        histotruncatedstk_signal_y = TH1F("truncated_stk_signal_allladdersy_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                          "truncated_stk_signal_allladdersy_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_STK,
                                          self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX)
         
        
        if self.amsclustersfound:
            cx_ams = TCanvas("Signal_AMS_xladders", "Signal_AMS_xladders", 1600, 1200)
            cy_ams = TCanvas("Signal_AMS_yladders", "Signal_AMS_yladders", 1600, 1200)
            cx_ams.Divide((len(self.amsparameters)+1)/2,2)
            cy_ams.Divide((len(self.amsparameters)+1)/2,2)
            
            cx_pos_ams = TCanvas("Position_AMS_xladders", "Position_AMS_xladders", 1600, 1200)
            cy_pos_ams = TCanvas("Position_AMS_yladders", "Position_AMS_yladders", 1600, 1200)
            cx_pos_ams.Divide((len(self.amsparameters)+1)/2,2)
            cy_pos_ams.Divide((len(self.amsparameters)+1)/2,2)
            
            cx_ncluster_ams = TCanvas("Nclusters_AMS_xladders", "Nclusters_AMS_xladders", 1600, 1200)
            cy_ncluster_ams = TCanvas("Nclusters_AMS_yladders", "Nclusters_AMS_yladders", 1600, 1200)
            cx_ncluster_ams.Divide((len(self.amsparameters)+1)/2,2)
            cy_ncluster_ams.Divide((len(self.amsparameters)+1)/2,2)
        
            histosamssignal_x = dict( (self.amsparameters[i]['ladderid'], TH1F("ams_signal_xside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                               "ams_signal_xside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                               self.AMS_SIGNAL_HISTO_NBINS,
                                                                               self.AMS_SIGNAL_HISTO_MIN,
                                                                               self.AMS_SIGNAL_HISTO_MAX)
                                       ) for i in xrange(len(self.amsparameters)) )
            
            histosamssignal_y = dict( (self.amsparameters[i]['ladderid'], TH1F("ams_signal_yside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                               "ams_signal_yside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                               self.AMS_SIGNAL_HISTO_NBINS,
                                                                               self.AMS_SIGNAL_HISTO_MIN,
                                                                               self.AMS_SIGNAL_HISTO_MAX)
                                       ) for i in xrange(len(self.amsparameters)) )
            hisosallamssignals_sside = dict( (self.amsparameters[i]['ladderid'], TH1F("ams_signal_sside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                                      "ams_signal_sside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                                      self.AMS_SIGNAL_HISTO_NBINS,
                                                                                      self.AMS_SIGNAL_HISTO_MIN,
                                                                                      self.AMS_SIGNAL_HISTO_MAX)
                                              ) for i in xrange(len(self.amsparameters)) )
            
            hisosallamssignals_kside = dict( (self.amsparameters[i]['ladderid'], TH1F("ams_signal_kside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                                      "ams_signal_kside_ladder%d"%self.amsparameters[i]['ladderid'],
                                                                                      self.AMS_SIGNAL_HISTO_NBINS,
                                                                                      self.AMS_SIGNAL_HISTO_MIN,
                                                                                      self.AMS_SIGNAL_HISTO_MAX)
                                              ) for i in xrange(len(self.amsparameters)) )
            
            
            histosamsposition_x = dict( (i, TH1F("ams_position_xside_layer%d"%i,
                                                 "ams_position_xside_layer%d"%i,
                                                 self.AMS_POSITION_HISTO_NBINS,
                                                 self.AMS_POSITION_HISTO_MIN,
                                                 self.AMS_POSITION_HISTO_MAX)
                                         ) for i in xrange(len(self.amsparameters)) )
            
            histosamsposition_y = dict( (i, TH1F("ams_position_yside_layer%d"%i,
                                                 "ams_position_yside_layer%d"%i,
                                                 self.AMS_POSITION_HISTO_NBINS,
                                                 self.AMS_POSITION_HISTO_MIN,
                                                 self.AMS_POSITION_HISTO_MAX)
                                         ) for i in xrange(len(self.amsparameters)) )
            
            histosamsnclusters_x = dict( (i, TH1F("ams_ncluster_xside_layer%d"%i,
                                                 "ams_nclusters_xside_layer%d"%i,
                                                 self.AMS_NCLUSTERS_HISTO_NBINS,
                                                 self.AMS_NCLUSTERS_HISTO_MIN,
                                                 self.AMS_NCLUSTERS_HISTO_MAX)
                                         ) for i in xrange(len(self.amsparameters)) )
            
            histosamsnclusters_y = dict( (i, TH1F("ams_ncluster_yside_layer%d"%i,
                                                 "ams_nclusters_yside_layer%d"%i,
                                                 self.AMS_NCLUSTERS_HISTO_NBINS,
                                                 self.AMS_NCLUSTERS_HISTO_MIN,
                                                 self.AMS_NCLUSTERS_HISTO_MAX)
                                         ) for i in xrange(len(self.amsparameters)) )
            
            
            
            
            histotruncatedams_sside = TH1F("truncated_ams_energy_allladders_sside_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_AMS_SSIDE,
                                           "truncated_ams_energy_allladders_sside_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_AMS_SSIDE,
                                           self.STK_SIGNAL_HISTO_NBINS, 
                                           self.STK_SIGNAL_HISTO_MIN, 
                                           self.STK_SIGNAL_HISTO_MAX*2)
            
            histotruncatedams_kside = TH1F("truncated_ams_energy_allladders_kside_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_AMS_KSIDE,
                                           "truncated_ams_energy_allladders_kside_minus_%d"%self.NCLSUTERS_TO_TRUNCATE_AMS_KSIDE,
                                           self.STK_SIGNAL_HISTO_NBINS, 
                                           self.STK_SIGNAL_HISTO_MIN, 
                                           self.STK_SIGNAL_HISTO_MAX*2)
            
        if self.stktracks and self.bgohits:
            stk_x_vs_bgo_energy = TH2D("stk_x_vs_bgo_energy", "stk_x_vs_bgo_energy", 800, -400, 400, 1000, 0, 100000)
            stk_y_vs_bgo_energy = TH2D("stk_y_vs_bgo_energy", "stk_y_vs_bgo_energy", 800, -400, 400, 1000, 0, 100000)
            
            
        nentries = min(self.tree.GetEntries(),self.MAX_EVENTS)
        for entry in xrange(nentries):
            sys.stdout.write("\r Processing file:  %10d / %10d ( %2d%% )"%(entry,nentries, 100 * entry/nentries))
            sys.stdout.flush()
            truncatedclusters = []
            truncatedclusters_x = []
            truncatedclusters_y = []
            self.tree.GetEntry(entry)
            
            
            if self.stktracks and self.bgohits:
                for i in xrange(self.stktracks.GetLast()+1):
                    track = self.stktracks.ConstructedAt(i)
                    if track.GetNPoints() < 6: continue
                    x = track.getImpactPoint().x()
                    y = track.getImpactPoint().y()
                    bgoenergy = sum(self.bgohits.fEnergy)
                    stk_x_vs_bgo_energy.Fill(x,bgoenergy)
                    stk_y_vs_bgo_energy.Fill(y,bgoenergy)
            
            
            
            #print "test:", self.stkclusters.GetLast()+1
            
            
            #if self.stkclusters.GetLast()+1> 100: continue
            
            if nstkclusterslayer:
                clusters_ladder_x = {}
                clusters_ladder_y = {}
                badevent = False
                #MAXCLUSTERS_PERLAYER = 4
                for i in xrange(self.stkclusters.GetLast()+1):
                    cluster = self.stkclusters.ConstructedAt(i)
                    #if cluster.isY(): continue
                    ladder = self.__get_stk_layer__(cluster)
                    
                    if cluster.isX():
                        if self.POSITION_CUT_STK[ladder]["minx"] > cluster.GetX(): continue
                        if self.POSITION_CUT_STK[ladder]["maxx"] < cluster.GetX(): continue
                        clusters_ladder_x[ladder] =  clusters_ladder_x[ladder]+1  if ladder in clusters_ladder_x.keys() else 1
                        if clusters_ladder_x[ladder] > nstkclusterslayer:
                            badevent = True
                            break
                    else:
                        if self.POSITION_CUT_STK[ladder]["miny"] > cluster.GetY(): continue
                        if self.POSITION_CUT_STK[ladder]["maxy"] < cluster.GetY(): continue
                        clusters_ladder_y[ladder] =  clusters_ladder_y[ladder]+1  if ladder in clusters_ladder_y.keys() else 1
                        if clusters_ladder_y[ladder] > nstkclusterslayer:
                            badevent = True
                            break
                if badevent:
                    continue
                
            
            
            
            
            neventsstk_x_layer = dict( [i,0] for i in xrange(len(self.STK_Z_COORDINATES_X)) )
            neventsstk_y_layer = dict( [i,0] for i in xrange(len(self.STK_Z_COORDINATES_X)) )
            for i in xrange(self.stkclusters.GetLast()+1):
                cluster = self.stkclusters.ConstructedAt(i)
                layer = self.__get_stk_layer__(cluster)
                energy = self.__get_stk_cluster_energy__(cluster)
                energy_05 = energy**0.5
                cluster.energy = energy
                if cluster.isX():
                    if self.POSITION_CUT_STK[layer]["minx"] > cluster.GetX(): continue
                    if self.POSITION_CUT_STK[layer]["maxx"] < cluster.GetX(): continue
                    
                                        
                    histosstksignal_x[layer].Fill(energy)
                    histosstksignal_x_2[layer].Fill(energy_05/self.MIP_MAXPOSITION_SQRT)
                    histosstkposition_x[layer].Fill(cluster.GetX())
                    neventsstk_x_layer[layer]+=1
                    truncatedclusters_x.append(cluster)
                else:
                    if self.POSITION_CUT_STK[layer]["miny"] > cluster.GetY(): continue
                    if self.POSITION_CUT_STK[layer]["maxy"] < cluster.GetY(): continue
                    histosstksignal_y[layer].Fill(energy)
                    histosstksignal_y_2[layer].Fill(energy_05/self.MIP_MAXPOSITION_SQRT)
                    histosstkposition_y[layer].Fill(cluster.GetY())
                    neventsstk_y_layer[layer]+=1
                    truncatedclusters_y.append(cluster)
                
                ladedrhardwareid = cluster.getLadderHardware()
                histoname = "signal_ladder_%d"%ladedrhardwareid
                if histoname not in hisosallstksignals.keys():
                    hisosallstksignals[histoname] = TH1F(histoname,histoname, self.STK_SIGNAL_HISTO_NBINS, self.STK_SIGNAL_HISTO_MIN, self.STK_SIGNAL_HISTO_MAX)
                hisosallstksignals[histoname].Fill(energy)
                truncatedclusters.append(cluster)
                histosstksignal_all.Fill(energy)
                histosstksignal_all_2.Fill(energy_05/self.MIP_MAXPOSITION_SQRT)
                
            truncatedclusters.sort(key = lambda x: x.energy)
            truncatedclusters_x.sort(key = lambda x: x.energy)
            truncatedclusters_y.sort(key = lambda x: x.energy)
            truncatedclusters = truncatedclusters[:-self.NCLSUTERS_TO_TRUNCATE_STK]
            truncatedclusters_x = truncatedclusters_x[:-self.NCLSUTERS_TO_TRUNCATE_STK]
            truncatedclusters_y = truncatedclusters_y[:-self.NCLSUTERS_TO_TRUNCATE_STK]
            if truncatedclusters:
                histotruncatedstk.Fill(sum(cluster.energy for cluster in truncatedclusters)/len(truncatedclusters))
                
            if truncatedclusters_x:
                energy_2 = sum(cluster.energy**0.5 for cluster in truncatedclusters_x)
                energy   = sum(cluster.energy for cluster in truncatedclusters_x)
                length = len(truncatedclusters_x)
                histotruncatedstk_sqrtsignal_x.Fill(energy_2/length/self.MIP_MAXPOSITION_SQRT)
                histotruncatedstk_signal_x.Fill(energy/length)
                
            if truncatedclusters_y:
                energy_2 = sum(cluster.energy**0.5 for cluster in truncatedclusters_y)
                energy   = sum(cluster.energy for cluster in truncatedclusters_y)
                length = len(truncatedclusters_y)
                histotruncatedstk_sqrtsignal_y.Fill(energy_2/length/self.MIP_MAXPOSITION_SQRT)
                histotruncatedstk_signal_y.Fill(energy/length)
                
                
            for layer, value in neventsstk_x_layer.items():
                histosstknclusters_x[layer].Fill(value)
                
            for layer, value in neventsstk_y_layer.items():
                histosstknclusters_y[layer].Fill(value)
            #for cluster in truncatedclusters:
            #    histotruncatedstk.Fill(cluster.energy)
            
                
                
            if not self.amsclustersfound: continue
            truncatedamsclusters_sside = []
            truncatedamsclusters_kside = []
            namsclustersx_ladder = dict( [ parameter["ladderid"],0] for parameter in self.amsparameters)
            namsclustersy_ladder = dict( [ parameter["ladderid"],0] for parameter in self.amsparameters)
            for i in xrange(self.amsclusters.GetLast()+1):
                cluster = self.amsclusters.ConstructedAt(i)
                energy = self.__get_ams_cluster_energy__(cluster)
                cluster.energy = energy
                if self.__is_ams_cluster_x__(cluster):
                    histosamssignal_x[cluster.ladder].Fill(energy)
                    histosamsposition_x[cluster.ladder].Fill(self.__get_ams_cluster_cog_x__(cluster))
                    namsclustersx_ladder[cluster.ladder]+=1
                else:
                    histosamssignal_y[cluster.ladder].Fill(energy)
                    histosamsposition_y[cluster.ladder].Fill(self.__get_ams_cluster_cog_y__(cluster))
                    namsclustersy_ladder[cluster.ladder]+=1
                    
                if self.__is_ams_cluster_sside__(cluster):
                    hisosallamssignals_sside[cluster.ladder].Fill(energy)
                    truncatedamsclusters_sside.append(cluster)
                else:
                    hisosallamssignals_kside[cluster.ladder].Fill(energy)
                    truncatedamsclusters_kside.append(cluster)
            truncatedamsclusters_sside.sort(key = lambda x: x.energy)
            truncatedamsclusters_kside.sort(key = lambda x: x.energy)
            truncatedamsclusters_sside = truncatedamsclusters_sside[:-self.NCLSUTERS_TO_TRUNCATE_AMS_SSIDE]
            truncatedamsclusters_kside = truncatedamsclusters_sside[:-self.NCLSUTERS_TO_TRUNCATE_AMS_KSIDE]
            if truncatedamsclusters_sside:
                histotruncatedams_sside.Fill(sum(cluster.energy for cluster in truncatedamsclusters_sside)/len(truncatedamsclusters_sside))
            
            if truncatedamsclusters_kside:
                histotruncatedams_kside.Fill(sum(cluster.energy for cluster in truncatedamsclusters_kside)/len(truncatedamsclusters_kside))
            
            for ladder, value in namsclustersx_ladder.items():
                histosamsnclusters_x[ladder].Fill(value)
                
            for ladder, value in namsclustersy_ladder.items():
                histosamsnclusters_y[ladder].Fill(value)
                
            
            
        
        plotdir.cd()
        stkdir = plotdir.mkdir("STK")
        stkdir.cd()            
        tmp = stkdir.mkdir("all_ladders")
        tmp.cd()
        [histo.Write() for histo in hisosallstksignals.values()]
        for i in xrange(len(histosstksignal_x)):
            cx_stk.cd(i+1)
            histosstksignal_x.values()[i].Draw()
            
        for i in xrange(len(histosstksignal_y)):
            cy_stk.cd(i+1)
            histosstksignal_y.values()[i].Draw()
            
        for i in xrange(len(histosstkposition_x)):
            cx_stk_position.cd(i+1)
            histosstkposition_x.values()[i].Draw()
            
        for i in xrange(len(histosstkposition_y)):
            cy_stk_position.cd(i+1)
            histosstkposition_y.values()[i].Draw()
            
        for i in xrange(len(histosstksignal_x_2)):
            cx_sqrt_stk.cd(i+1)
            histosstksignal_x_2.values()[i].Draw()
            
        for i in xrange(len(histosstksignal_y_2)):
            cy_sqrt_stk.cd(i+1)
            histosstksignal_y_2.values()[i].Draw()
            
        for i in xrange(len(histosstknclusters_x)):
            cx_stk_nclusters.cd(i+1)
            histosstknclusters_x.values()[i].Draw()
            
        for i in xrange(len(histosstknclusters_y)):
            cy_stk_nclusters.cd(i+1)
            histosstknclusters_y.values()[i].Draw()
            
            
            
        if self.amsclustersfound:
            plotdir.cd()
            amsdir = plotdir.mkdir("AMS")
            amsdir.cd()
            for i in xrange(len(histosamssignal_x)):
                cx_ams.cd(i+1)
                histosamssignal_x.values()[i].Draw()
                
            for i in xrange(len(histosamssignal_y)):
                cy_ams.cd(i+1)
                histosamssignal_y.values()[i].Draw()
                
            for i in xrange(len(histosamsposition_x)):
                cx_pos_ams.cd(i+1)
                histosamsposition_x.values()[i].Draw()
                
            for i in xrange(len(histosamsposition_y)):
                cy_pos_ams.cd(i+1)
                histosamsposition_y.values()[i].Draw()
                
            for i in xrange(len(histosamsnclusters_x)):
                cx_ncluster_ams.cd(i+1)
                histosamsnclusters_x.values()[i].Draw()
                
            for i in xrange(len(histosamsnclusters_y)):
                cy_ncluster_ams.cd(i+1)
                histosamsnclusters_y.values()[i].Draw()
        
                
            tmp = amsdir.mkdir("all_ladders_sside_kside")
            tmp.cd()
            [histo.Write() for histo in hisosallamssignals_sside.values()]
            [histo.Write() for histo in hisosallamssignals_kside.values()]
            
        
        stkdir.cd()
        cx_stk.Write()
        cy_stk.Write()
        cx_stk_nclusters.Write()
        cy_stk_nclusters.Write()
        cx_stk_position.Write()
        cy_stk_position.Write()
        cx_sqrt_stk.Write()
        cy_sqrt_stk.Write()
        histotruncatedstk.Write()
        histotruncatedstk_sqrtsignal_y.Write()
        histotruncatedstk_sqrtsignal_x.Write()
        histotruncatedstk_signal_y.Write()
        histotruncatedstk_signal_x.Write()
        histosstksignal_all.Write()
        histosstksignal_all_2.Write()
        
        if self.stktracks and self.bgohits:
            stk_x_vs_bgo_energy.Write()
            stk_y_vs_bgo_energy.Write()
        
        if self.amsclustersfound:
            amsdir.cd()
            cx_ams.Write()
            cy_ams.Write()
            cx_pos_ams.Write()
            cy_pos_ams.Write()
            cx_ncluster_ams.Write()
            cy_ncluster_ams.Write()      
            histotruncatedams_sside.Write()
            histotruncatedams_kside.Write()                 
             
        
    
    
    """
        auxdir = topdir.mkdir("aux")
        auxdir.cd()
        [histo.Write() for histo in histosx.values()]
        [histo.Write() for histo in histosy.values()]
    """
        
    """
    def __plot_cluster_distributions__(self):
        hisots_stk_x = [TH1F("clusterx_layerfrombgo_%d"%i, "clusterx_layerfrombgo_%d"%i, 
                             self.OCCUP_PLOTS_NBINS_STK,
                             self.OCCUP_PLOTS_MIN_STK, 
                             self.OCCUP_PLOTS_MAX_STK 
                             ) for i in xrange(len(self.STK_Z_COORDINATES_X))]
        hisots_stk_y = [TH1F("clustery_layerfrombgo_%d"%i, "clustery_layerfrombgo_%d"%i, 
                             self.OCCUP_PLOTS_NBINS_STK,
                             self.OCCUP_PLOTS_MIN_STK, 
                             self.OCCUP_PLOTS_MAX_STK 
                             ) for i in xrange(len(self.STK_Z_COORDINATES_Y))]
        hisots_ams_x = [TH1F("clustery_layerfrombgo_%d"%i, "clustery_layerfrombgo_%d"%i, 
                             self.OCCUP_PLOTS_NBINS_STK,
                             self.OCCUP_PLOTS_MIN_STK, 
                             self.OCCUP_PLOTS_MAX_STK 
                             ) for i in xrange(len(self.STK_Z_COORDINATES_Y))]
    """
        
        
                
            
    def create_plots(self, nevents, skipevents, nstkclusterslayer, displaysonly, sec, msec, mintracks, maxtracks, minbgoenergy, maxbgoenergy):
        
        #Some sainity checks
        if (sec or msec) and not self.eventheader:
            print "ERROR: EventHeader object is not available in the ROOT file, while [sec]/[msec] is specified"
            raise SystemExit

        #Event displays
        self.fout = TFile(VALIDATION_PLOTS_FILE, "RECREATE")
        self.eventdiplaydir = self.fout.mkdir("EventDisplays")
        self.eventdiplaydir.cd()
        processedevents = 0
        nentries =  self.tree.GetEntries()
        for i in xrange(skipevents, nentries):
            #self.current_entry_id = i
            if processedevents>=self.MAX_EVENTS or processedevents>=nevents: break
            sys.stdout.write("\r Processing file: %2d%% "%(100 * (i-skipevents)/nentries))
            sys.stdout.flush()
            
            # header selection
            if sec is not None or msec is not None:
                self.b_eventheader.GetEntry(i)
#mintracks=args["mintracks"],maxtracks=args["maxtracks"],minbgoenergy=args["minbgoenergy"],maxbgoenergy=args["maxbgoenergy"]
                if sec is not None  and self.eventheader.GetSecond()!= sec: continue
                if msec is not None and self.eventheader.GetMillisecond()!= msec: continue

            # track selection 
            if mintracks is not None or maxtracks is not None:
                self.b_stktracks.GetEntry(i)
                if mintracks is not None and self.stktracks.GetLast()+1 < mintracks: continue
                if maxtracks is not None and self.stktracks.GetLast()+1 > maxtracks: continue

            # bgo selection
            if minbgoenergy is not None or maxbgoenergy is not None:
                self.b_bgorec.GetEntry(i)
                if minbgoenergy is not None and self.bgorec.GetTotalEnergy() < minbgoenergy: continue
                if maxbgoenergy is not None and self.bgorec.GetTotalEnergy() > maxbgoenergy: continue

            # mc branches
            if self.b_mcprimaries:
                self.b_mcprimaries.GetEntry(i) # read out the MC branch            

            
            self.tree.GetEntry(i)
            s = self.eventheader.GetSecond()
            ms = self.eventheader.GetMillisecond()
            date = DmpStkHkeepHeader.TimeCode2TimeStamp_PMO(s)
            self.__plot_event_diplay__("%d_event_%ds_%dms-%s"%(processedevents,s,ms,date), i)
            processedevents+=1
            
        
        if not displaysonly:
            #Correlation plots
            if self.amsclustersfound:
                self.corrdir = self.fout.mkdir("StkAmsCorrelations")
                self.corrdir.cd()
                self.__plot_stk_ams_correlations__(self.corrdir)
        
        
            #Signal
            self.amssignaldir = self.fout.mkdir("Signal")
            self.amssignaldir.cd()
            self.__plot_signal_all_ladders__(self.amssignaldir, nstkclusterslayer)
        
        
        
        
        #Cluster occupancies, energy, number of x,y clusters per AMS ladder
        #self.datafile.Close()
        self.fout.Close()        
        print 
        print "SUCCESS!  File with the plots produced: ",VALIDATION_PLOTS_FILE 


def __get_module_name__():
    return str(sys.modules[__name__]).split('from')[-1].split("'")[1].split("/")[-1]


def __print_help__():
    print """
    Usage:
        %s MyDampeAmsDataFile.root [--amsconfig=MyAmsConfigFile.txt] [--neventdisplays=N] [--nstkcluslayer=N] [--displaysonly] [--nskipeventdisplays=N]
 
            MyDampeAmsDataFile.root  - name of ROOT file with reconstructed data 
          
            --amsconfig            :  (optional) name of ams configuration file, holding the AMS ladder indices, z poistions, alignment constants, etc.
            --neventdisplays       :  (optional) number of event diplays to store in the output validation file
            --nstkcluslayer        :  (optional) maximum number of stk cluster per layer
            --displaysonly         :  (optional) if specified, only eventdisplays are created
            --nskipeventdisplays   :  (optional) number of event displays to skip in the beginning of file
            --sec                  :  (optional) timestamp (in seconds) of event to display, e.g. --sec=75053758
            --msec                 :  (optional) timestamp (in miliseconds) of event to display, e.g. --sec=603    
            --mintracks            :  (optional) minimal number of STK tracks, e.g. --mintracks=1
            --maxtracks            :  (optional) maximum number of STK tracks
            --minbgoenergy         :  (optional) minimal BGO total energy (MeV), e.g. --minbgoenergy=1000. 
            --maxbgoenergy         :  (optional) maximum BGO total energy (MeV)
            --showprimarydirection :  (optional) show direction of primary particle
            
    Example:
        %s  DAMPE_AMS_ANC_20141108-174841.root --amsconfig=ams_config_2014PS.txt --neventdisplays=50"
            
    """%(__get_module_name__(),__get_module_name__())

def run_script():
    args = {"neventdisplays":9999999,
            #"maxevents":99}
            "maxevents":9999999,
            "amsconfig":None,
            "nskipeventdisplays":0,
            "timestampsec":None,
            "timestampmsec":None,
            "rootfile":[],
            "mintracks":None,     
            "maxtracks":None,
            "minbgoenergy":None,
            "maxbgoenergy":None,
            "showprimarydirection":False
            }
    
    nstkclusterslayer = None
    displaysonly      = False
    
    #print sys.argv
    if printhelp:
        __print_help__()
        return
    
    for item in sys.argv[1:]:
        if "-H" in item:
            __print_help__()
            return
        elif "--amsconfig" in item:
            args["amsconfig"] = item.split("=")[1].strip()
        #elif "--rootfile" in item:
        #    args["rootfile"] =  item.split("=")[1].strip()
        elif ".root" in item.lower():
            args["rootfile"].append(item)
        elif "--neventdisplays" in item:
            args["neventdisplays"] = int(item.split("=")[1])
        elif "--maxevents" in item:
            args["maxevents"] = int(item.split("=")[1])
        elif "--nstkcluslayer" in item:
            nstkclusterslayer = int(item.split("=")[1])
        elif "--displaysonly" in item:
            displaysonly = True
        elif "--nskipeventdisplays" in item:
            args["nskipeventdisplays"] = int(item.split("=")[1])
        elif "--sec" in item:
            args["timestampsec"] = int(item.split("=")[1])
        elif "--msec" in item:
            args["timestampmsec"] = int(item.split("=")[1])
        elif "--mintracks" in item:
            args["mintracks"] = int(item.split("=")[1])
        elif "--maxtracks" in item:
            args["maxtracks"] = int(item.split("=")[1])
        elif "--minbgoenergy" in item:
            args["minbgoenergy"] = float(item.split("=")[1])
        elif "--maxbgoenergy" in item:
            args["maxbgoenergy"] = float(item.split("=")[1])
        elif "--showprimarydirection" in item:
            args["showprimarydirection"] = True
        elif "--" in item or "-" in item:
            print "Unrecognized parameter: ", item, "==> see <%s -H> for more information"%__get_module_name__()
            return
    
    if "amsconfig" not in args.keys():
        print "\n<amsconfig> parameter is not set ==> see <%s -H> for more information\n"%__get_module_name__()
        return 
        
    if not args["rootfile"]:
        print "\n<rootfile> parameter is not set ==> see <%s -H> for more information\n"%__get_module_name__()
        return
    
    thevalidator = validator(args["amsconfig"], args["rootfile"], args["maxevents"], args["showprimarydirection"])
    thevalidator.create_plots(skipevents = args["nskipeventdisplays"], nevents = args["neventdisplays"], nstkclusterslayer = nstkclusterslayer, displaysonly = displaysonly, sec = args["timestampsec"] , msec = args["timestampmsec"],mintracks=args["mintracks"],maxtracks=args["maxtracks"],minbgoenergy=args["minbgoenergy"],maxbgoenergy=args["maxbgoenergy"])
    
run_script()            