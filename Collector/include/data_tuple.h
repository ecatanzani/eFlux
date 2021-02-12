#ifndef DATA_TUPLE_H
#define DATA_TUPLE_H

#include "tuple.h"
#include "tmpstruct.h"
#include "data_tmpstruct.h"

#include "DmpEvtAttitude.h"

class data_tuple : public ntuple
{
public:
	data_tuple(const active_cuts &acuts)
	{
		init(acuts);
	};
	~data_tuple(){};

	void Fill(
		const std::shared_ptr<_tmp_filter> _filter_res,
		const std::vector<int> _stk_clusters_on_plane,
		const std::shared_ptr<_tmp_bgo> _bgo_res,
		const std::shared_ptr<_tmp_energy_data> _energy_res,
		const std::shared_ptr<DmpEvtAttitude> attitude,
		const std::shared_ptr<DmpEvtHeader> header,
		const std::shared_ptr<_tmp_nud> _nud_res);
		
	void Reset();
	
private:
	void init(const active_cuts &acuts);
	void branch_tree();
	void fill_filter_info(const filter_output &output);
	void fill_attitude_info(const std::shared_ptr<DmpEvtAttitude> &attitude);

	// Event time
	unsigned int second = 0;
	// Attitude
	double glat = -999;
	double glon = -999;
	double geo_lat = -999;
	double geo_lon = -999;
	double ra_zenith = -999;
	double dec_zenith = -999;
	double ra_scz = -999;
	double dec_scz = -999;
	double ra_scx = -999;
	double dec_scx = -999;
	double ra_scy = -999;
	double dec_scy = -999;
	double cutoff = -999;
	// Filters
	bool evtfilter_evt_in_saa = false;
};

#endif