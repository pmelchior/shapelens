#ifndef SHAPELENS_HISTORY_H
#define SHAPELENS_HISTORY_H

#include <sstream>
#include <iostream>
#include <string>
#include "ShapeLens.h"

namespace shapelens {

  /// History class.
  /// The class stores the history of processing steps.
  /// New text is added by using the \p operator<< (borrowed from 
  /// std::ostringstream).\n
  /// If Config::VERBOSITY is set to \p true, everything appended to History
  /// is also printed to \p stdout.

  class History {
  public:
    /// Default constructor.
    History() {
      s.clear();
      silent = 0;
    }
    /// Copy constructor
    History(const History& h) {
      s.clear();
      s << h.s.str();
      silent = h.silent;
    }
    /// Copy operator.
    History& operator=(const History& h) {
      s.clear();
      s << h.s.str();
      silent = h.silent;
      return *this;
    }
    /// Addition operator.
    History& operator+=(const History& h) {
      s << h.s.str();
      return *this;
    }
    /// Overloaded \p operator<<.
    /// With this operator, History behaves like a \p std::ostringstream.
    template <typename T> History& operator<<(T t) {
      if (Config::VERBOSITY && !silent)
	std::cout << t;
      s << t; 
      return *this; 
    }

    History& operator<<(std::ostream& (*func)(std::ostream&)) {
      if (Config::VERBOSITY && !silent)
	std::cout << func;
      s << func;
      return *this;
    } 

    History& operator<<(std::ios& (*func)(std::ios&)) {
      if (Config::VERBOSITY && !silent)
	std::cout << func;
      s << func;
      return *this;
    }

    History& operator<<(std::ios_base& (*func)(std::ios_base&)) {
      if (Config::VERBOSITY && !silent)
	std::cout << func;
      s << func;
      return *this;
    }
    /// Return the history string.
    std::string str() const {
      return s.str();
    }
    /// Clear the history.
    void clear() {
      s.str("");
    }
    /// Set the verbosity of \p operator<<.
    /// If set to 1, all texts added to History are printed to stdout.
    /// The default value is 0.
    static void setVerbosity(bool v) {
      Config::VERBOSITY = v;
    }
    /// Get the actual verbosity.
    bool getVerbosity() const {
      return Config::VERBOSITY;
    }
    /// Turn verbose mode off temporarily.
    /// This does not change Config::VERBOSITY, but turns
    /// the output off in cases where it is unwanted.
    void setSilent() {
      silent = 1;
    }
    /// Undo setSilent().
    void unsetSilent() {
      silent = 0;
    }
    /// Whether history is empty.
    bool isEmpty() const {
      return !(bool) s.str().size();
    }
  private:
    std::ostringstream s;
    bool silent;
  }; 
} // end namespace
#endif
