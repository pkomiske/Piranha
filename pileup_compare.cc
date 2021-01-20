// C++ standard library
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// OpenMP
#include <omp.h>

// Boost
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

// FastJet
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/tools/GridMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/Transformer.hh"
#include "fastjet/contrib/ConstituentSubtractor.hh"
#include "fastjet/contrib/IteratedVoronoiSubtractor.hh"

// HepPID - for particle charges
#include "HepPID/ParticleIDMethods.hh"

using fastjet::PseudoJet;

class EventSource {
public:

  EventSource() {}
  virtual ~EventSource() {}

  virtual bool append_next(std::vector<PseudoJet> & particles) = 0;

  // this version uses the internal particles vector
  bool next() {
    particles_.clear();
    return append_next(particles_);
  }

  const std::vector<PseudoJet> & particles() const { return particles_; }

private:

  std::vector<PseudoJet> particles_;

}; // EventSource


class FileEventSource : public EventSource {
public:

  FileEventSource(const std::string & filepath) :
    gzipped_(filepath.size() >= 3 && 0 == filepath.compare(filepath.size() - 3, 3, ".gz"))
  {
    if (!filepath.empty()) {
      if (gzipped_) {
        file_.open(filepath, std::ios_base::in | std::ios_base::binary);
        filtered_stream_.push(boost::iostreams::gzip_decompressor());
      }
      else file_.open(filepath);

      if (!file_.good())
        throw std::runtime_error("Failed to open " + filepath);

      filtered_stream_.push(file_);  
    }
  }
  virtual ~FileEventSource() {}

  boost::iostreams::filtering_istream & file() { return filtered_stream_; }

private:

  bool gzipped_;
  std::ifstream file_;
  boost::iostreams::filtering_istream filtered_stream_;

}; // FileEventSource


class PileupInfo : public PseudoJet::UserInfoBase {
public:

  PileupInfo(int pdgid, int vertex = 0, int barcode = 0) :
    pdgid_(pdgid),
    three_charge_(HepPID::threeCharge(pdgid_)),
    vertex_(vertex),
    barcode_(barcode)
  {}

  // access functions
  int pdgid() const { return pdgid_; }
  int three_charge() const { return three_charge_; }
  int vertex() const { return vertex_; }
  int barcode() const { return barcode_; }

  bool is_charged() const { return three_charge_ != 0; }
  double charge() const { return three_charge_ / 3.0; }

private:

  int pdgid_, three_charge_, vertex_, barcode_;

}; // PileupInfo


class CERN2014PileupWorkshopFileEventSource : public FileEventSource {
public:

  CERN2014PileupWorkshopFileEventSource(const std::string & filepath) :
    FileEventSource(filepath),
    vertex_(0)
  {}

  bool append_next(std::vector<PseudoJet> & particles) {

    std::string line;
    double px, py, pz, m, E;
    int pdgid;

    size_t orig_size(particles.size());

    while (std::getline(file(), line)) {

      // ignore blank lines and comments
      if (line.length() == 0 || line[0] == '#') continue;

      // end of event
      if (line[0] == 'e' && line.substr(0, 3) == "end") break;

      // read particle and load into back of the vector
      std::istringstream s(line);
      s >> px >> py >> pz >> m >> pdgid;
      E = std::sqrt(px*px + py*py + pz*pz + m*m);

      particles.emplace_back(px, py, pz, E);
      particles.back().set_user_info(new PileupInfo(pdgid, vertex_));
    }

    return particles.size() != orig_size;
  }

  void set_vertex(int v) { vertex_ = v; }

private:

  int vertex_;

}; // CERN2014PileupWorkshopFileEventSource


class MasslessTransformer : public fastjet::Transformer {
public:

  MasslessTransformer() {}
  std::string description() const { return ""; }
  PseudoJet result(const PseudoJet & pj) const {
    PseudoJet newpj(pj);
    newpj.reset_momentum_PtYPhiM(pj.pt(), pj.rap(), pj.phi());
    return newpj;
  }

}; // MasslessTransformer


class CERN2014PileupWorkshopEventMixer {
public:

  CERN2014PileupWorkshopEventMixer(const po::variables_map & vm) :
    hard_filepath_(vm["hard"].as<std::string>()),
    pileup_filepath_(vm["pileup"].as<std::string>()),
    poisson_(vm["poisson"].as<bool>()),
    mu_(vm["mu"].as<int>()),
    chs_factor_(vm["chs-factor"].as<double>()),
    massless_(vm["massless"].as<bool>()),

    // initialize randomness
    rng_(),
    poisson_dist_(mu_ > 0 ? mu_ : 1.0),

    // setup event sources
    hard_events_(hard_filepath_),
    pileup_events_(pileup_filepath_)
  {
    if (mu_ > 0 && pileup_filepath_.empty())
      throw std::runtime_error("requested pileup but did not provide a file");
  }

  bool next() {

    // reset particles
    particles_.clear();

    // determine npu for this event
    npu_ = (mu_ == 0 || !poisson_ ? mu_ : poisson_dist_(rng_));

    // first get the hard event
    if (!hard_events_.append_next(particles_)) return false;
    size_t hard_size(particles_.size());

    // now append pileup events
    for (int i = 1; i <= npu_; i++) {
      pileup_events_.set_vertex(i);
      if (!pileup_events_.append_next(particles_)) return false;
    }

    // make massless if requested
    if (massless_)
      particles_ = MasslessTransformer()(particles_);

    // apply CHS
    if (chs_factor_ != 1)
      for (size_t i = hard_size; i < particles_.size(); i++)
        if (particles_[i].user_info<PileupInfo>().is_charged())
          particles_[i] *= chs_factor_;

    return true;
  }

  void clear() {
    particles_.clear();
  }

  // particles for last mixed event
  const std::vector<PseudoJet> & particles() const { return particles_; }

  // number of pileup events in last mixed event
  int npu() const { return npu_; }

private:

  // parameters in vm
  std::string hard_filepath_, pileup_filepath_;
  bool poisson_;
  int mu_;
  double chs_factor_;
  bool massless_;

  // randomness
  std::default_random_engine rng_;
  std::poisson_distribution<int> poisson_dist_;

  // event sources
  CERN2014PileupWorkshopFileEventSource hard_events_, pileup_events_;

  // particles of the current mixed event
  std::vector<PseudoJet> particles_;
  int npu_;

}; // CERN2014PileupWorkshopEventMixer


class SelectorWorkerVertex : public fastjet::SelectorWorker {
public:
  
  SelectorWorkerVertex(int vertex) : vertex_(vertex) {}

  bool pass(const PseudoJet & particle) const {
    // we check that the user_info_ptr is non-zero so as to make
    // sure that explicit ghosts don't cause the selector to fail
    return (particle.user_info_ptr() != 0 && 
            particle.user_info<PileupInfo>().vertex() == vertex_);
  }
  
  std::string description() const {
    std::ostringstream ostr;
    ostr << "vertex == " << vertex_;
    return ostr.str();
  }

private:

  int vertex_;

}; // SelectorWorkerVertex


fastjet::Selector SelectorVertex(int vertex) {
  return new SelectorWorkerVertex(vertex);
}

fastjet::Selector SelectorIsHard() { return SelectorVertex(0); }


int main(int argc, char** argv) {

  // variables
  int nev, max_print, num_threads, print_every;
  unsigned njet;
  double jet_R, rapmax, match;
  bool rap_rescale, massless;
  std::string outputfp;

  // catch bad input
  po::variables_map vm;
  try {

    po::options_description od("Options for pileup removal comparison");

    // add options in a chained manner
    od.add_options()
      ("hard",        po::value<std::string>()->default_value(""), "Path to file with hard events")
      ("pileup",      po::value<std::string>()->default_value(""), "Path to file with pileup events")
      ("output,o",    po::value<std::string>(&outputfp)->default_value("out.txt"), "Path to output file")
      ("nev,n",       po::value<int>(&nev)->default_value(1), "Number of events")
      ("R",           po::value<double>(&jet_R)->default_value(0.4), "Jet radius")
      ("njet",        po::value<unsigned>(&njet)->default_value(2), "Number of jets to keep")
      ("rapmax",      po::value<double>(&rapmax)->default_value(4), "Maximum rapidity to consider")
      ("match",       po::value<double>(&match)->default_value(1), "Matching factor between jets (multiplied by R)")
      ("rap-rescale", po::bool_switch(&rap_rescale), "Includes rescaling by rapidity")
      ("poisson,p",   po::bool_switch(), "Draws the number of pileup events from a Poisson distribution with the given mean")
      ("mu",          po::value<int>()->default_value(0), "The pileup amount")
      ("chs-factor",  po::value<double>()->default_value(1e-60), "Amount to scale charged particles by")
      ("massless",    po::bool_switch(&massless), "Whether to make particles massless")
      ("max-print",   po::value<int>(&max_print)->default_value(100), "Number of events to print out")
      ("num-threads,t", po::value<int>(&num_threads)->default_value(omp_get_max_threads()), "Number of threads to use")
      ("print-every", po::value<int>(&print_every)->default_value(100), "How often to print updates")
      ("help,h",      "Prints the help message")
    ; // ends add_options()

    // parse command line
    po::store(po::command_line_parser(argc, argv).options(od).run(), vm);

    // print help if requested
    if (vm.count("help")) {
      std::cout << od << '\n';
      return 1;
    }

    // check things
    po::notify(vm);
  }

  // catch exceptions from po::notify
  catch (std::exception & e) {
    std::cerr << e.what() << std::endl;
    throw e;
  }

  // setup jet finding
  fastjet::ClusterSequence::print_banner();

  // setup events
  CERN2014PileupWorkshopEventMixer mixer(vm);

  // setup output
  std::ofstream of(outputfp);
  of << std::setprecision(10);

  // parallelized event section
  int iEvent(0);
  bool end(false);
  #pragma omp parallel num_threads(num_threads) \
                       default(none) \
                       firstprivate(jet_R, njet, rapmax, massless, rap_rescale, max_print, print_every) \
                       shared(std::cout, iEvent, nev, mixer, end, of)
  {
    int thread_i(omp_get_thread_num());

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_R);
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts);
    fastjet::Selector jet_sel(fastjet::SelectorAbsRapMax(rapmax - jet_R) * fastjet::SelectorNHardest(njet));

    // setup area subtraction
    auto rescaling = new fastjet::BackgroundRescalingYPolynomial(1.1685397, 0, -0.0246807, 0, 5.94119e-05);
    fastjet::GridMedianBackgroundEstimator gmbge(rapmax, 0.55);
    if (rap_rescale) {
      gmbge.set_rescaling_class(rescaling);
      gmbge.set_compute_rho_m(!massless);  
    }
    fastjet::Subtractor area_sub(&gmbge);
    area_sub.set_use_rho_m(!massless);

    // setup other pileup subtraction
    fastjet::contrib::piranha::IteratedVoronoiSubtractorCylinder ivsub(rapmax);
    ivsub.set_background_estimator(&gmbge);
    fastjet::contrib::ConstituentSubtractor cssub;

    cssub.set_distance_type(fastjet::contrib::ConstituentSubtractor::deltaR);
    cssub.set_max_distance(0.3);
    cssub.set_alpha(0);
    cssub.set_ghost_area(0.01);
    cssub.set_max_eta(rapmax);
    cssub.set_background_estimator(&gmbge);
    cssub.initialize();

    #pragma omp single
    {
      std::cout << ivsub.description() << std::endl
                << cssub.description() << std::endl;  
    }

    std::vector<PseudoJet> event;
    int npu, iEvent_i;

    // loop forever
    while (true) {

      // get another event
      #pragma omp critical (read)
      {
        if (!end && mixer.next() && ((iEvent_i = iEvent++) < nev || nev == -1)) {
          event = mixer.particles();
          npu = mixer.npu();
          mixer.clear();
          if (iEvent < max_print)
            std::cout << "Event " << iEvent_i
                      << ", npu = " << npu
                      << ", mult = " << event.size() << '\n';
        }
        else end = true;

        if (iEvent % print_every == 0)
          std::cout << "Done with " << iEvent << " events" << std::endl;
      }
      if (end) break;

      // get hard event and pileup event
      std::vector<PseudoJet> hard_event, pileup_event;
      SelectorIsHard().sift(event, hard_event, pileup_event);

      // cluster jets
      fastjet::ClusterSequence cs_hard(hard_event, jet_def);
      fastjet::ClusterSequenceArea cs_full(event, jet_def, area_def);
      std::vector<PseudoJet> hard_jets(fastjet::sorted_by_pt(jet_sel(cs_hard.inclusive_jets()))),
                             full_jets(fastjet::sorted_by_pt(jet_sel(cs_full.inclusive_jets())));

      // do area subtraction
      gmbge.set_particles(event);
      std::vector<PseudoJet> area_sub_jets(fastjet::sorted_by_pt(area_sub(full_jets)));

      // do IVS subtraction
      std::vector<PseudoJet> iv_sub_event(ivsub(event));

      // find IVS jets
      fastjet::ClusterSequence cs_iv(iv_sub_event, jet_def);
      std::vector<PseudoJet> iv_sub_jets(fastjet::sorted_by_pt(jet_sel(cs_iv.inclusive_jets())));

      // constituent subtraction
      std::vector<PseudoJet> cs_sub_event(cssub.subtract_event(event));
      fastjet::ClusterSequence cs_cs(cs_sub_event, jet_def);
      std::vector<PseudoJet> cs_sub_jets(fastjet::sorted_by_pt(jet_sel(cs_cs.inclusive_jets())));

      // print pt of leading jet
      #pragma omp critical (write)
      {
        if (iEvent_i < max_print)
          std::cout << "Event " << iEvent_i << ", pt_hard_jet0 = " << hard_jets[0].pt()
                                            << ", pt_full_jet0 = " << full_jets[0].pt() 
                                            << ", pt_corr_jet0 = " << area_sub_jets[0].pt()
                                            << ", pt_ivcr_jet0 = " << iv_sub_jets[0].pt() 
                                            << ", pt_cscr_jet0 = " << cs_sub_jets[0].pt() 
                                            << std::endl;
        of << npu                   << ' '
           << hard_jets[0].pt()     << ' '
           << full_jets[0].pt()     << ' '
           << area_sub_jets[0].pt() << ' '
           << cs_sub_jets[0].pt()   << ' '
           << iv_sub_jets[0].pt()   << ' '
           << hard_jets[0].m()      << ' '
           << full_jets[0].m()      << ' '
           << area_sub_jets[0].m()  << ' '
           << cs_sub_jets[0].m()    << ' '
           << iv_sub_jets[0].m()    << ' '
           << hard_jets[1].pt()     << ' '
           << full_jets[1].pt()     << ' '
           << area_sub_jets[1].pt() << ' '
           << cs_sub_jets[1].pt()   << ' '
           << iv_sub_jets[1].pt()   << ' '
           << hard_jets[1].m()      << ' '
           << full_jets[1].m()      << ' '
           << area_sub_jets[1].m()  << ' '
           << cs_sub_jets[1].m()    << ' '
           << iv_sub_jets[1].m()    
           << '\n';
      }
    }

    #pragma omp barrier
  }

  std::cout << "I've finished my pileup removal!" << std::endl;
}