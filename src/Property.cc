#include "../include/Property.h"
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <stdexcept>

namespace shapelens {

Property::Property() : std::map<std::string, variant_t>() {}

variant_t& Property::operator[](const std::string& key) {
  return std::map<std::string, variant_t>::operator[](key);
}

const variant_t& Property::operator[](const std::string& key) const {
  Property::const_iterator iter = Property::find(key);
  if (iter != Property::end())
    return iter->second;
  else
    throw std::invalid_argument("Property: Key " + key + " does not exist!");
}


// Helper class for performing write()
class OutFunctor : public boost::static_visitor<void> {
public:
  OutFunctor(std::ostream& o, std::string linesep) : out(o),ls(linesep) {
  }
  void operator()(int e) const {
    out << "\tI\t" << e;
    if (ls.size() == 0)
      out << std::endl;
    else
      out << ls;
  }
  void operator()(data_t e) const {
    out << "\tD\t" << e;
    if (ls.size() == 0)
      out << std::endl;
    else
      out << ls;
  }
  void operator()(const std::string& e) const {
    out << "\tS\t" << e;
    if (ls.size() == 0)
      out << std::endl;
    else
      out << ls;
  }
  void operator()(const std::vector<int>& e) const {
    out << "\tVI\t";
    writeVector(e);
  }
  void operator()(const std::vector<data_t>& e) const {
    out << "\tVD\t";
    writeVector(e);
  }
  void operator()(const std::vector<std::string>& e) const {
    out << "\tVS\t";
    writeVector(e);
  }
private:
  template<class T>
  void writeVector(const std::vector<T>& v) const {
    for (typename std::vector<T>::const_iterator iter = v.begin(); iter != v.end(); iter++) {
      out << *iter;
      if (iter != --v.end())
	out  << ",";
    }
    if (ls.size() == 0)
      out << std::endl;
    else
      out << ls;
  }
  std::ostream & out;
  std::string  ls;
};

void Property::write(std::ostream& o,std::string linesep) const {
  OutFunctor out(o,linesep);
  for (Property::const_iterator iter = std::map<std::string, variant_t>::begin(); iter != std::map<std::string, variant_t>::end(); iter++) {
    o << iter->first << "\t";
    boost::apply_visitor(out, iter->second);
  }
}

void Property::read(std::istream& in) {
  std::map<std::string, variant_t>::clear();
  std::string line, name, type, value;
  typedef boost::tokenizer<boost::char_separator<char> > Tok;
  boost::char_separator<char> fieldsep("\t");
  boost::char_separator<char> vectorsep(",");
  while(getline(in, line)) {
    if (line[0] != '#' && line != "") {
      Tok tok(line, fieldsep);
      int i=0;
      for(Tok::iterator tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter) {
        switch(i) {
        case 0: name = *tok_iter; break;
        case 1: type = *tok_iter; break;
        case 2: value = *tok_iter; break; 
        }
        i++;
        if (i>2)
          break;
      } 
      // name, type, value filled
      if (type[0] == 'V') { // vector type
	Tok vtok(value, vectorsep);
	if (type[1] == 'I') {
	  std::vector<int> vi;
	  for(Tok::iterator tok_iter = vtok.begin(); tok_iter != vtok.end(); ++tok_iter) 
	    vi.push_back(boost::lexical_cast<int>(*tok_iter));
	  std::map<std::string, variant_t>::insert(std::make_pair(name,vi));
	} else if (type[1] == 'D') {
	  std::vector<data_t> vd;
	  for(Tok::iterator tok_iter = vtok.begin(); tok_iter != vtok.end(); ++tok_iter) 
	    vd.push_back(boost::lexical_cast<data_t>(*tok_iter));
	  std::map<std::string, variant_t>::insert(std::make_pair(name,vd));
	} else if (type[1] == 'S') {
	  std::vector<std::string> vs;
	  for(Tok::iterator tok_iter = vtok.begin(); tok_iter != vtok.end(); ++tok_iter) 
	    vs.push_back(*tok_iter);
	  std::map<std::string, variant_t>::insert(std::make_pair(name,vs));
	}
      } else {
	if (type[0] == 'S')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,value));
	else if (type[0] == 'I')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,boost::lexical_cast<int>(value)));
	else if (type[0] == 'D')
	  std::map<std::string, variant_t>::insert(std::make_pair(name,boost::lexical_cast<data_t>(value)));
      }
    }
  }
}

} // end namespace
