/*!
  \file main.cc
  \brief pull in reference allele frequencies for results files
  \copyright Released under the MIT License.
  Copyright 2020 Cameron Palmer.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "annotate_frequency/cargs.h"
#include "finter/finter.h"

bool cicompare(const std::string &s1, const std::string &s2) {
  if (s1.size() != s2.size()) return false;
  for (unsigned i = 0; i < s1.size(); ++i) {
    if (toupper(s1.at(i)) != toupper(s2.at(i))) return false;
  }
  return true;
}

class freq_handler {
 public:
  freq_handler()
      : _filename(""),
        _metadata(""),
        _input(0),
        _input_metadata(0),
        _chr(0),
        _pos(0),
        _ref(""),
        _target_index(0) {}
  explicit freq_handler(const std::string &filename,
                        const std::string &metadata)
      : _filename(filename),
        _metadata(metadata),
        _input(0),
        _input_metadata(0),
        _chr(0),
        _pos(0),
        _ref(""),
        _target_index(0) {}
  freq_handler(const freq_handler &obj)
      : _filename(""),
        _metadata(""),
        _input(0),
        _input_metadata(0),
        _chr(0),
        _pos(0),
        _ref(""),
        _target_index(0) {
    throw std::domain_error("copy constructor not permitted");
  }
  ~freq_handler() throw() {
    if (_input) delete _input;
    if (_input_metadata) delete _input_metadata;
  }
  void initialize(const std::string &supercontinent) {
    _input = finter::reconcile_reader(_filename);
    if (!_metadata.empty()) {
      _input_metadata = finter::reconcile_reader(_metadata);
    }
    std::string line = "", id = "";
    if (!_input->getline(&line))
      throw std::domain_error("no header available from file \"" + _filename +
                              "\"");
    if (_metadata.empty()) {
      std::istringstream strm1(line);
      for (unsigned i = 0; i < 6; ++i) {
        if (!(strm1 >> id))
          throw std::domain_error("cannot parse header intro info");
      }
      unsigned counter = 1;
      while (true) {
        if (!(strm1 >> id)) {
          break;
        }
        if (cicompare(id, supercontinent)) {
          _target_index = counter + 6;
          break;
        }
        ++counter;
      }
      if (!_target_index)
        throw std::domain_error("unable to locate header \"" + supercontinent +
                                "\" from frequency file \"" + _filename + "\"");
    } else {
      _input_metadata->getline(&line);
    }
  }

  bool align(unsigned chr, unsigned pos, const std::string &id);

  bool find(unsigned chr, unsigned pos, const std::string &id,
            const std::string &a1, const std::string &a2, double *freq);

 private:
  std::string _filename;
  std::string _metadata;
  finter::finter_reader *_input;
  finter::finter_reader *_input_metadata;
  unsigned _chr;
  unsigned _pos;
  std::string _ref;
  std::vector<std::string> _alt;
  std::vector<double> _freq;
  unsigned _target_index;
};

bool freq_handler::align(unsigned chr, unsigned pos, const std::string &id) {
  std::string line = "", freqline = "", catcher = "", ref = "", alt = "";
  unsigned next_chr = 0, next_pos = 0;
  double freq = 0.0;
  bool dup_position_failure = false;
  if (!_input) throw std::domain_error("align called on null input pointer");
  if (!_metadata.empty() && !_input_metadata)
    throw std::domain_error("align called on null metadata pointer");
  while (next_chr < chr || (next_chr == chr && next_pos <= pos)) {
    if (_input_metadata) {
      if (!_input_metadata->getline(&line)) break;
    } else {
      if (!_input->getline(&line)) break;
    }
    std::istringstream strm1(line);
    if (!(strm1 >> catcher >> next_chr >> next_pos >> ref >> alt)) {
      throw std::domain_error("cannot parse frequency file \"" + _filename +
                              "\" line \"" + line + "\"");
    }
    if (!catcher.compare(id.substr(3))) {
      _chr = next_chr;
      _pos = next_pos;
      _ref = ref;
      _alt.resize(1);
      _alt.at(0) = alt;
      _freq.resize(1);
      if (_input_metadata) {
        if (!_input->getline(&freqline)) break;
        std::istringstream strm2(freqline);
        if (!(strm2 >> freq))
          throw std::domain_error("cannot parse frequency split file \"" +
                                  _filename + "\" line \"" + freqline + "\"");
      } else {
        strm1 >> catcher;
        for (unsigned i = 6; i < _target_index; ++i) {
          if (!(strm1 >> freq))
            throw std::domain_error("cannot parse frequency from line \"" +
                                    line + "\"");
        }
      }
      _freq.at(0) = freq;
      return true;
    } else if (_input_metadata) {
      _input->getline(&freqline);
    }
    if (next_chr == chr && next_pos == pos) {
      dup_position_failure = true;
    }
  }
  /*
    this condition was added to detect situations in which multiallelics were
    split and ordered in such a way that you could reasonably miss one of the
    alternates because your input file's ordering differed from that of the
    reference file. this then promptly never came up as an issue during the
    entire first imputation round with PLCO. it has however come up as an issue
    for the second round. i don't think the severity of the issue will be
    very high? for the moment, going to change this to emit a warning message
    to better evaluate the number of possible instances of this being a problem

    I think overriding the error is fine because the logic seems to indicate
    that the issue will just introduce sporadic NA frequencies where in fact
    some data were available. but compared to the total number of variants,
    this really isn't much of a problem; and the code should still run...
  */
  if (dup_position_failure)
    std::cerr << "warning: possible multiallelic sorting issue for chr " << chr
              << " pos " << pos << std::endl;
  return false;
}

bool freq_handler::find(unsigned chr, unsigned pos, const std::string &id,
                        const std::string &a1, const std::string &a2,
                        double *freq) {
  if (!freq)
    throw std::domain_error(
        "freq_handler::find: called with null freq pointer");
  *freq = -1.0;
  std::string alt_target = "";
  if (align(chr, pos, id)) {
    if (cicompare(a1, _ref)) {
      alt_target = a2;
    } else if (cicompare(a2, _ref)) {
      alt_target = a1;
    } else {
      return false;
    }
    for (unsigned i = 0; i < _alt.size(); ++i) {
      if (cicompare(alt_target, _alt.at(i))) {
        *freq = !alt_target.compare(a2) ? _freq.at(i) : 1.0 - _freq.at(i);
        return true;
      }
    }
  }
  return false;
}

void process_file(const std::string &input_filename,
                  const std::string &freq_filename,
                  const std::string &freq_metadata,
                  const std::string &supercontinent,
                  const std::string &output_filename) {
  finter::finter_reader *input = 0;
  finter::finter_writer *output = 0;
  freq_handler freq(freq_filename, freq_metadata);
  freq.initialize(supercontinent);
  std::string line = "", id = "", a1 = "", a2 = "", catcher = "", n = "",
              beta = "", se = "", p = "";
  unsigned chr = 0, pos = 0, total_input = 0, mapped_input = 0,
           unmapped_input = 0;
  double updated_freq = 0.0;
  try {
    input = finter::reconcile_reader(input_filename);
    output = finter::reconcile_writer(output_filename);
    input->getline(&line);
    output->writeline(line);
    while (input->getline(&line)) {
      ++total_input;
      std::istringstream strm1(line);
      if (!(strm1 >> chr >> pos >> id >> a1 >> a2 >> catcher >> beta >> se >>
            p >> n))
        throw std::domain_error("cannot parse file \"" + input_filename +
                                "\" line \"" + line + "\"");
      std::ostringstream o;
      if (freq.find(chr, pos, id, a2, a1, &updated_freq)) {
        ++mapped_input;
        if (!(o << chr << '\t' << pos << '\t' << id << '\t' << a1 << '\t' << a2
                << '\t' << updated_freq << '\t' << beta << '\t' << se << '\t'
                << p << '\t' << n)) {
          throw std::domain_error("cannot format output");
        }
        while (strm1 >> catcher) {
          if (!(o << '\t' << catcher))
            throw std::domain_error("cannot format output");
        }
        output->writeline(o.str());
      } else {
        ++unmapped_input;
        output->writeline(line);
      }
      if (total_input % 100000 == 0) {
        std::cout << "\tprocessed " << total_input << "; mapped "
                  << mapped_input << ", unmapped " << unmapped_input
                  << std::endl;
      }
    }
    input->close();
    delete input;
    input = 0;
    output->close();
    delete output;
    output = 0;
  } catch (...) {
    if (input) delete input;
    if (output) delete output;
    throw;
  }
}

int main(int argc, char **argv) {
  try {
    // parse command line input
    annotate_frequency::cargs ap(argc, argv);
    // if help is requested or no flags specified
    if (ap.help() || argc == 1) {
      // print a help message and exist
      ap.print_help(std::cout);
      return 0;
    }
    std::string input_filename = ap.get_input_filename();
    std::string freq_filename = ap.get_frequency_filename();
    std::string freq_metadata = ap.get_frequency_metadata_filename();
    std::string supercontinent = ap.get_supercontinent();
    std::string output_filename = ap.get_output_filename();
    std::cout << "starting processing of \"" << input_filename << "\" with \""
              << freq_filename << "\"" << std::endl;
    process_file(input_filename, freq_filename, freq_metadata, supercontinent,
                 output_filename);
    std::cout << "all done: \"" << output_filename << "\"" << std::endl;
    return 0;
  } catch (const std::domain_error &e) {
    std::cerr << "error: " << e.what() << std::endl;
    return 1;
  }
}
