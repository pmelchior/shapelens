#include "../include/FITS.h"
#include <boost/tokenizer.hpp>
#include <fstream>

namespace shapelens {

  fitsfile* FITS::openFile(const std::string& filename, bool write) {
    int status = 0;
    fitsfile* outfptr;
    fits_open_file(&outfptr, filename.c_str(), (int) write, &status);
    if (status != 0)
      throw std::runtime_error("FITS: Cannot open " + filename + ": " + getErrorMessage(status));
    return outfptr;
  }

  fitsfile* FITS::openTable(const std::string& filename, bool write) {
    int status = 0;
    fitsfile* outfptr;
    fits_open_table(&outfptr, filename.c_str(), (int) write, &status);
    if (status != 0)
      throw std::runtime_error("FITS: Cannot open " + filename + ": " + getErrorMessage(status));
    return outfptr;
  }
			      
  fitsfile* FITS::createFile(const std::string& filename) {
    int status = 0;
    fitsfile *outfptr;
    std::string newfilename = "!"+filename; // overwrite existing file if necessary
    // create fits file
    fits_create_file(&outfptr,newfilename.c_str(), &status);
    if (status != 0)
      throw std::runtime_error("FITS: Cannot create " + filename + ": " + getErrorMessage(status));
    return outfptr;
  }

  void FITS::closeFile(fitsfile* fptr) {
    int status = 0;
    fits_close_file(fptr, &status);
    if (status != 0)
	  throw std::runtime_error("FITS: Cannot close FITS file " + getFileName(fptr) + ": " + getErrorMessage(status));
  }

  std::string FITS::getFileName(fitsfile *fptr) {
    int status = 0;
    char header[FLEN_FILENAME];
    fits_file_name(fptr, header, &status);
    if (status != 0)
      throw std::invalid_argument("FITS: Cannot get filename pointer: " + getErrorMessage(status));
    return std::string(header);
  }

  std::string FITS::getErrorMessage(int status) {
    char err[FLEN_STATUS];
    fits_get_errstatus(status, err);
    std::ostringstream note;
    note << err << " (" << status << ")";
    return note.str();
  }

  void FITS::moveToExtension(fitsfile* fptr, unsigned int i) {
    int status = 0;
    fits_movabs_hdu(fptr, i, NULL, &status);
    if (status != 0) {
      std::ostringstream note;
      note << "FITS: Cannot move to extension " << i << " in " + getFileName(fptr) << ": " << getErrorMessage(status);
      throw std::runtime_error(note.str());
    }
  }

  void FITS:: moveToExtension(fitsfile* fptr, const std::string& name) {
    int status = 0, hdutype = ANY_HDU, extver = 0;
    fits_movnam_hdu(fptr, hdutype, const_cast<char*>(name.c_str()), extver, &status);
    if (status != 0) {
      std::ostringstream  note;
      note << "FITS: Cannot move to extension " << name << " in " << getFileName(fptr) << ": " << getErrorMessage(status);
      throw std::runtime_error(note.str());
    }
  }

  void FITS::appendHistory(fitsfile *fptr, const std::string& history) {
    // since it is too long the be saved in one shot
    // split it line by line
    // and for each line insert a HISTORY line in the header
    boost::char_separator<char> sep("\n");
    boost::tokenizer<boost::char_separator<char> > tok(history, sep);
    boost::tokenizer<boost::char_separator<char> >::iterator tok_iter;
    int status = 0, spaces = 8;
    for(tok_iter = tok.begin(); tok_iter != tok.end(); ++tok_iter) {
      // replace tab by spaces because tab cannot be saved in FITS header card
      // in order to mimick the tab style, an adaptive amount of spaces is inserted
      std::string card = (*tok_iter);
      char tabReplacement(' ');
      while ((card).find("\t") != std::string::npos) {
	int index = (card).find("\t");
	card.replace(index,1,spaces - (index%spaces), tabReplacement);
      }
      fits_write_history (fptr,const_cast<char*>(card.c_str()), &status);
    }
    if (status != 0)
      throw std::runtime_error("FITS: Cannot append FITS history to " + getFileName(fptr) + ": " + getErrorMessage(status));
  }

  void FITS::readKeyCards(fitsfile *fptr, const std::string& key, std::string& value) {
    int status = 0, nkeys, keylen=0;
    char* comment = NULL;
    char card[FLEN_CARD], keyword[FLEN_KEYWORD];
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    for (int i = 1; i <= nkeys; i++) {
      fits_read_record(fptr, i, card, &status);
      fits_get_keyname(card,keyword,&keylen,&status);
      if (std::string(keyword) == key) {
	// append all but the first 8 characters
	value += card + 8;
	value += "\n";
      }
    }
    if (status != 0)
      throw std::invalid_argument("FITS: Cannot read FITS keycards from " + getFileName(fptr) + ": " + getErrorMessage(status));
  }

  long FITS::getTableRows(fitsfile* fptr) {
    int status = 0;
    long nrows;
    fits_get_num_rows(fptr, &nrows, &status);
    if (status != 0)
      throw std::runtime_error("FITS: Cannot find number of rows in table in " + getFileName(fptr) + ": " + getErrorMessage(status));
    return nrows;
  }

  int FITS::getTableColumnNumber(fitsfile* fptr, const std::string& name) {
    int status = 0, colnum;
    fits_get_colnum(fptr, 0, const_cast<char*>(name.c_str()), &colnum, &status);
    if (status == COL_NOT_UNIQUE) {
      std::ostringstream note;
      note << "FITS: Column " << name + " in FITS table is not unique in " << getFileName(fptr) << ": " << getErrorMessage(status);
      throw std::runtime_error(note.str());
    }
    else if (status != 0) {
      std::ostringstream note;
      note << "FITS: Cannot find column " << name << " in FITS table in " << getFileName(fptr)  << ": " << getErrorMessage(status);
      throw std::invalid_argument(note.str());
    }
    return colnum;
  }
  
  int FITS::getTableColumnType(fitsfile* fptr, int colnr) {
    int status = 0, typecode;
    long repeat, width;
    fits_get_coltype(fptr, colnr, &typecode, &repeat, &width, &status);
    if (status != 0) {
      std::ostringstream note;
      note << "FITS: Cannot find type of column " << colnr << " in FITS table in " << getFileName(fptr) << ": " << getErrorMessage(status);
      throw std::invalid_argument(note.str());
    }
    return typecode;
  }
  
} // end namespace
