#include "trigger.h"

bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
    auto trigger_mip1 = evt_header->GeneratedTrigger(1);
	auto trigger_mip2 = evt_header->GeneratedTrigger(2);
	auto trigger_HET = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
	auto trigger_LET = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
	auto trigger_MIP = trigger_mip1 || trigger_mip2;
	auto trigger_general = trigger_MIP || trigger_HET || trigger_LET;
    return trigger_general;
}

bool check_HE_trigger(const std::shared_ptr<DmpEvtHeader> evt_header) {
	return evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
}