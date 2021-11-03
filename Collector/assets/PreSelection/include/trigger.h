#ifndef TRIGGER_H
#define TRIGGER_H

#include <memory>

#include "DmpEvtHeader.h"

extern bool check_trigger(const std::shared_ptr<DmpEvtHeader> evt_header);
extern bool check_HE_trigger(const std::shared_ptr<DmpEvtHeader> evt_header);

#endif