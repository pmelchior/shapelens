#include "../include/ShapeLens.h"
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>

using namespace shapelens;
using namespace std;

bool Config::VERBOSITY = false;
bool Config::CHECK_OBJECT = true;
data_t Config::ADD_BORDER = 0.5;
bool Config::USE_WCS = false;

Config::Config() {
}

Config::Config(string filename) {
  // open config file
  ifstream configfile (filename.c_str());
  if (configfile.fail()) {
    cerr << "Config: configuration file does not exists!" << endl;
    terminate();
  }
  // read in config file
  string line;
  while(getline(configfile, line)) {
    typedef boost::tokenizer<boost::char_separator<char> > Tok;
    // split entries at tabs
    boost::char_separator<char> sep("\t");
    Tok tok(line, sep);
    // first of all we copy the token into string vector
    // though this is not too sophisticated it does not hurt and is more 
    // convenient
    std::vector<std::string> column;
    for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter)
      column.push_back(*tok_iter);
    // exclude empty and comment lines
    if (column.size() >= 2 && column[0] != "#") {
      if (column[0] == "VERBOSITY")
	VERBOSITY = boost::lexical_cast<bool>(column[1].c_str());
      if (column[0] == "CHECK_OBJECT")
        CHECK_OBJECT = boost::lexical_cast<bool>(column[1].c_str());
      if (column[0] == "ADD_BORDER")
        ADD_BORDER = boost::lexical_cast<data_t>(column[1].c_str());
      if (column[0] == "USE_WCS")
	USE_WCS = boost::lexical_cast<bool>(column[1].c_str());
    }
  }
}

