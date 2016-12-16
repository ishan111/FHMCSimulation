#ifndef INPUT_H_
#define INPUT_H_

#include "system.h"
#include <string>

// JSON interface from local distro of rapidjson
#include "../src/rapidjson/include/rapidjson/document.h"
#include "../src/rapidjson/include/rapidjson/writer.h"
#include "../src/rapidjson/include/rapidjson/stringbuffer.h"
#include "../src/rapidjson/include/rapidjson/filereadstream.h"
#include "../src/rapidjson/include/rapidjson/prettywriter.h"

void initializeSystemBarriers (simSystem &sys, const rapidjson::Document &doc);
void parse_json (const std::string filename, simSystem &sys);

#endif
