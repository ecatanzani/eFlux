#include "Dmp/DmpStkContainer.h"

#include "TMath.h"

void DmpStkContainer::scanSTKHits(
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const std::vector<double> bgoRec_slope,
    const std::vector<double> bgoRec_intercept)
{
    // Fund the number of clusters per layer
    for (int idx=0; idx<stkclusters->GetEntries(); ++idx)
        ++clusters_on_plane[static_cast<DmpStkSiCluster*>(stkclusters->ConstructedAt(idx))->getPlane()];
    
    // Find the number of clusters and their energy within 1 Moliere radius
    clusters_1_rm(stkclusters, stktracks, bgoRec_slope, bgoRec_intercept);
}

void DmpStkContainer::clusters_1_rm(
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const std::vector<double> bgoRec_slope,
    const std::vector<double> bgoRec_intercept)
    {
        double mindAngleTrackBgoRec 					{185};
        double mindDrTopTrackBgoRec						{999};
        DmpStkTrack *mindAngleTrack 					{nullptr};
        std::vector<double> mindAngleTrack_slope 		(2, -999);
        std::vector<double> mindAngleTrack_intercept 	(2, -999);
        std::vector<int> LadderToLayer					(nSTKladders, -1);
        TVector3 mindAngleTrackDirection;
        TVector3 bgoRecEntrance;
        TVector3 bgoRecDirection;
        
        const bool removeTrack3Clus 					{true};
        const bool removeTrack2Holes 					{true};
        const bool useTrack3clu 						{true};

        bool passTrackMatch			{false};
        bool passTrackMatch_3clu	{false};

        link_ladders(LadderToLayer);
        fill_BGO_vectors(
            bgoRecEntrance,
            bgoRecDirection,
            bgoRec_slope,
            bgoRec_intercept);

        // Loop over tracks
        for (unsigned int itrk=0; itrk<stktracks->GetLast()+1; ++itrk)
        {	
            std::vector<int> track_nHoles(2, 0);

            auto track = static_cast<DmpStkTrack*>(stktracks->ConstructedAt(itrk));

            get_track_points(
                stkclusters,
                track,
                LadderToLayer,
                track_nHoles);

            // filter the number of X and Y hits
            if ((track->getNhitX() <= 3 || track->getNhitY() <= 3) && removeTrack3Clus)
                continue;

            // filter the number of X and Y holes
            if ((track_nHoles[0] > 1 || track_nHoles[1] > 1) && removeTrack2Holes)
                continue;

            // Get the track direction
            auto trackDirection = (track->getDirection()).Unit();
            // Get the angular distance between the track and the shower axis
            auto dAngleTrackBgoRec = trackDirection.Angle(bgoRecDirection) * TMath::RadToDeg();
            // Extrapolate to the top of the BGO
            double topX 	{track->getTrackParams().getSlopeX() * BGO_TopZ + track->getTrackParams().getInterceptX()};
            double topY		{track->getTrackParams().getSlopeY() * BGO_TopZ + track->getTrackParams().getInterceptY()};
            double dxTop	{topX - bgoRecEntrance[0]};
            double dyTop	{topY - bgoRecEntrance[1]};
            double drTop	{sqrt(pow(dxTop, 2) + pow(dyTop, 2))};

            if (drTop < mindDrTopTrackBgoRec && dAngleTrackBgoRec < 10)
            {
                mindAngleTrack = static_cast<DmpStkTrack*>(track);
                mindDrTopTrackBgoRec = drTop;
                mindAngleTrackDirection = trackDirection;
                mindAngleTrackBgoRec = dAngleTrackBgoRec;
                mindAngleTrack_slope[0] = track->getTrackParams().getSlopeX();
                mindAngleTrack_slope[1] = track->getTrackParams().getSlopeY();
                mindAngleTrack_intercept[0] = track->getTrackParams().getInterceptX();
                mindAngleTrack_intercept[1] = track->getTrackParams().getInterceptY();
            }
        }

        passTrackMatch = mindAngleTrackBgoRec > 10 ? false : true;

        if (!passTrackMatch && useTrack3clu)
        {
            // Reainitialize variables
            mindAngleTrackBgoRec = 185;
            mindDrTopTrackBgoRec = 999;
            mindAngleTrack = nullptr;
            mindAngleTrack_slope = std::vector<double> (2, -999);
            mindAngleTrack_intercept = std::vector<double> (2, -999);
            mindAngleTrackDirection = TVector3();

            // Loop over tracks
            for (int itrk=0; itrk<stktracks->GetLast()+1; ++itrk)
            {	
                std::vector<int> track_nHoles(2, 0);

                auto track = static_cast<DmpStkTrack*>(stktracks->ConstructedAt(itrk));

                get_track_points(
                    stkclusters,
                    track,
                    LadderToLayer,
                    track_nHoles);

                // filter the number of X and Y hits
                if (track->getNhitX() != 3 && track->getNhitY() != 3)
                    continue;

                // filter the number of X and Y holes
                if (track->getNhitX() == 3 && track_nHoles[0] > 0)
                    continue;
                if (track->getNhitY() == 3 && track_nHoles[1] > 0)
                    continue;
                if (track_nHoles[0] > 1 || track_nHoles[1] > 1)
                    continue;
                if (track->getNhitXY() < 3)
                    continue;
                if (track->getNhitX() + track->getNhitY() - 2 * track->getNhitXY() > 1)
                    continue;

                // Get the track direction
                auto trackDirection = (track->getDirection()).Unit();
                // Get the angular distance between the track and the shower axis
                auto dAngleTrackBgoRec = trackDirection.Angle(bgoRecDirection) * TMath::RadToDeg();
                // Extrapolate to the top of the BGO
                double topX 	{track->getTrackParams().getSlopeX() * BGO_TopZ + track->getTrackParams().getInterceptX()};
                double topY		{track->getTrackParams().getSlopeY() * BGO_TopZ + track->getTrackParams().getInterceptY()};
                double dxTop	{topX - bgoRecEntrance[0]};
                double dyTop	{topY - bgoRecEntrance[1]};
                double drTop	{sqrt(pow(dxTop, 2) + pow(dyTop, 2))};
            
                if (drTop < mindDrTopTrackBgoRec && dAngleTrackBgoRec < 10)
                {
                    mindAngleTrack = static_cast<DmpStkTrack*>(track);
                    mindDrTopTrackBgoRec = drTop;
                    mindAngleTrackDirection = trackDirection;
                    mindAngleTrackBgoRec = dAngleTrackBgoRec;
                    mindAngleTrack_slope[0] = track->getTrackParams().getSlopeX();
                    mindAngleTrack_slope[1] = track->getTrackParams().getSlopeY();
                    mindAngleTrack_intercept[0] = track->getTrackParams().getInterceptX();
                    mindAngleTrack_intercept[1] = track->getTrackParams().getInterceptY();
                }
            
            }

            passTrackMatch_3clu = mindAngleTrackBgoRec > 40 ? false : true;
        }

        // Find the number of clusters (and their energy) within 1 Molier radius within the minimum distance track
		TVector3 trackStkTopPoint;
		
        trackStkTopPoint[0] = mindAngleTrack_slope[0] * STK_TopZ + mindAngleTrack_intercept[0];
        trackStkTopPoint[1] = mindAngleTrack_slope[1] * STK_TopZ + mindAngleTrack_intercept[1];
		trackStkTopPoint[2] = STK_TopZ;
		
		TVector3 XbarDirection(0, 1, 0);
  		TVector3 YbarDirection(1, 0, 0);

		const double Rm_W {9.327};

		std::vector<int> LadderToIsN (nSTKladders, -1);

		for (int ilad = 0; ilad < nSTKladders; ++ilad)
		{
			int iTRB = ilad / 24;
			int isN = iTRB / 4;
			LadderToIsN[ilad] = isN;
		}

		// Loop over STK clusters
		for (int iclu=0; iclu<stkclusters->GetLast()+1; ++iclu)
		{
			auto cluster = static_cast<DmpStkSiCluster*>(stkclusters->ConstructedAt(iclu));

			int isX 		{cluster->isX()};
			double hitX		{cluster->GetX()};
			double hitY		{cluster->GetY()};
			double hitZ		{cluster->GetZ()};
			int hardID		{cluster->getLadderHardware()};
			double hitE		{cluster->getEnergy()};
			int isN			{LadderToIsN[hardID]};
			//int iLay		{LadderToLayer[hardID]};

			int k 				{isX ? 0 : 1};
        	int m 				{isX ? 1 : 0};
       		double thisCoord 	{isX ? hitX : hitY};

			TVector3 mindAnglePoint(mindAngleTrack_slope[0] * hitZ + mindAngleTrack_intercept[0],
                                mindAngleTrack_slope[1] * hitZ + mindAngleTrack_intercept[1], hitZ);
        	TVector3 bgoRecPoint(bgoRec_slope[0] * hitZ + bgoRec_intercept[0],
                             	bgoRec_slope[1] * hitZ + bgoRec_intercept[1], hitZ);

			int number_inputs {!passTrackMatch && !passTrackMatch_3clu ? 1 : 2};

			for (int input=0; input < number_inputs; ++input)
			{
				TVector3 thisPoint;
          		TVector3 thisDirection;

				if (!input)
				{
					thisPoint = bgoRecPoint;
					thisDirection = bgoRecDirection;
				}
				if (input)
				{
					thisPoint = mindAnglePoint;
					thisDirection = mindAngleTrackDirection;
				}

				double predToCenter {0};
				if ((isN && thisPoint[m] > 0) || (!isN && thisPoint[m] < 0))
					predToCenter = fabs(thisPoint[m]);
				
				thisPoint[k] = thisCoord;
				TVector3 trackStkTopPointToThis 	{thisPoint - trackStkTopPoint};
				//double longiLength 				{thisDirection * trackStkTopPointToThis};
          		//double transLength 				{(trackStkTopPointToThis - longiLength * thisDirection).Mag()};
				TVector3 stripDirection 			{isX ? XbarDirection : YbarDirection};
         		double dxy 							{stripDirection * thisDirection};
				TVector3 crossStripTrack 			{stripDirection.Cross(thisDirection)};
				TVector3 perpStripToTrack 			{thisDirection - dxy * stripDirection};
				double transLength 					{fabs(trackStkTopPointToThis * crossStripTrack) / crossStripTrack.Mag()};
				//double diffTrans 					{(transLength - transLength) / Rm_W};

				double dRmOther 					{predToCenter / Rm_W};
          		double dRm 							{transLength / Rm_W};
				
				if (dRm < 1 && dRmOther < 1)
				{
					stkEcore1Rm[input] += hitE;
					++nStkClu1Rm[input];
				}
			}
		}
    }

void DmpStkContainer::link_ladders(std::vector<int> &LadderToLayer)
{
	for (int ilad = 0; ilad < nSTKladders; ++ilad)
	{
		int iTRB = ilad / 24;
		int iladTRB = ilad % 24;
		int iPlane = 5 - iladTRB / 4;
		int isY = (iTRB / 2 + 1) % 2;
		int iLay = iPlane * 2 + isY;
		LadderToLayer[ilad] = iLay;
	}
}

void DmpStkContainer::fill_BGO_vectors(
	TVector3 &bgoRecEntrance,
	TVector3 &bgoRecDirection,
	const std::vector<double> bgoRec_slope,
	const std::vector<double> bgoRec_intercept)
    {
        // Build bgoRecDirection TVector3
        TVector3 vec_s0_a(bgoRec_intercept[0], bgoRec_intercept[1], 0.);
        TVector3 vec_s1_a(bgoRec_intercept[0] + bgoRec_slope[0], bgoRec_intercept[1] + bgoRec_slope[1], 1.);
        bgoRecDirection = (vec_s1_a - vec_s0_a).Unit(); //uni vector pointing from front to back

        // Build bgoRecEntrance TVector3
        double topZ = BGO_TopZ;
        double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
        double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

        if (fabs(topX) > BGO_SideXY || fabs(topY) > BGO_SideXY)
        {
            // possibly enter from the x-sides
            if (fabs(topX) > BGO_SideXY)
            {
                if (topX > 0)
                    topX = BGO_SideXY;
                else
                    topX = -BGO_SideXY;
                topZ = (topX - bgoRec_intercept[0]) / bgoRec_slope[0];
                topY = bgoRec_slope[1] * topZ + bgoRec_intercept[1];
                // possibly enter from the y-sides
                if (fabs(topY) > BGO_SideXY)
                {
                    if (topY > 0)
                        topY = BGO_SideXY;
                    else
                        topY = -BGO_SideXY;
                    topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
                    topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
                }
            }
            //enter from the y-sides
            else if (fabs(topY) > BGO_SideXY)
            {
                if (topY > 0)
                    topY = BGO_SideXY;
                else
                    topY = -BGO_SideXY;
                topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
                topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
            }
        }

        bgoRecEntrance[0] = topX;
        bgoRecEntrance[1] = topY;
        bgoRecEntrance[2] = topZ;
    }

void DmpStkContainer::get_track_points(
    const std::shared_ptr<TClonesArray> stkclusters,
	DmpStkTrack *track,
	const std::vector<int> LadderToLayer,
	std::vector<int> &track_nHoles)
    {
        std::vector<int> prevHole(2, -2);
        std::vector<int> firstLayer(2, -1);
        std::vector<int> lastLayer(2, -1);
        std::vector<int> lastPoint(2, -1);
        std::vector<unsigned int> track_nHoles_cont(2, 0);

        // Loop on track points to find last layer values
        for (int ip = track->GetNPoints() - 1; ip >= 0; --ip)
        {
            if (lastLayer[0] == -1)
            {
                if (track->getHitMeasX(ip) > -99999)
                {
                    lastPoint[0] = ip;
                    DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                    auto hardID = cluster->getLadderHardware();
                    lastLayer[0] = LadderToLayer[hardID];
                }
            }
            if (lastLayer[1] == -1)
            {
                if (track->getHitMeasY(ip) > -99999)
                {
                    lastPoint[1] = ip;
                    DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
                    auto hardID = cluster->getLadderHardware();
                    lastLayer[1] = LadderToLayer[hardID];
                }
            }
        }

        // Found the number of holes on both X and Y
        for (int ip = 0; ip <= lastPoint[0]; ++ip)
        {
            if (track->getHitMeasX(ip) > -99999)
            {
                DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                if (firstLayer[0] == -1)
                    firstLayer[0] = LadderToLayer[hardID];
            }
            else
            {
                if (firstLayer[0] != -1)
                    ++track_nHoles[0];
                if (ip == prevHole[0] + 1)
                    ++track_nHoles_cont[0];
                prevHole[0] = ip;
            }
        }

        for (int ip = 0; ip <= lastPoint[1]; ++ip)
        {
            if (track->getHitMeasY(ip) > -99999)
            {
                DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                if (firstLayer[1] == -1)
                    firstLayer[1] = LadderToLayer[hardID];
            }
            else
            {
                if (firstLayer[1] != -1)
                    ++track_nHoles[1];
                if (ip == prevHole[1] + 1)
                    ++track_nHoles_cont[1];
                prevHole[1] = ip;
            }
        }
    }

const std::vector<int> DmpStkContainer::GetNPlaneClusters()
{
    return clusters_on_plane;
}

const std::vector<double> DmpStkContainer::GetStkEcore1Rm()
{
    return stkEcore1Rm;
}

const std::vector<unsigned int> DmpStkContainer::GetNStkClu1Rm()
{
    return nStkClu1Rm;
}

void DmpStkContainer::Reset()
{
    clusters_on_plane = std::vector<int>(DAMPE_stk_planes, 0);
    stkEcore1Rm = std::vector<double> (2, 0);
    nStkClu1Rm = std::vector<unsigned int> (2, 0);
}