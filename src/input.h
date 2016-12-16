#ifndef INPUT_H_
#define INPUT_H_

#include <string>
#include "moves.h"
#include "system.h"
#include "insert.h"
#include "delete.h"
#include "translate.h"
#include "swap.h"

// JSON interface from local distro of rapidjson
#include "rapidjson/include/rapidjson/document.h"
#include "rapidjson/include/rapidjson/writer.h"
#include "rapidjson/include/rapidjson/stringbuffer.h"
#include "rapidjson/include/rapidjson/filereadstream.h"
#include "rapidjson/include/rapidjson/prettywriter.h"

void checkBounds (simSystem &sys);
void initializeSystemBarriers (simSystem &sys, const rapidjson::Document &doc);
void setPairPotentials (simSystem &sys, const rapidjson::Document &doc);
void setup (simSystem &sys, const std::string filename);
simSystem initialize (const std::string filename);

#endif
