#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Math/Vector4.hh"
#include "fastjet/contrib/IFNPlugin.hh"
#include "fastjet/contrib/CMPPlugin.hh"
#include "fastjet/contrib/GHSAlgo.hh"
#include "fastjet/contrib/SDFlavPlugin.hh"

//#define DebugLog

using namespace fastjet;
using namespace fastjet::contrib;

/// \class Angularity
/// definition of angularity
///
class Angularity :
//public FunctionOfPseudoJet<double>,
public PseudoJet{
public:
  /// ctor
  Angularity(double alpha, double jet_radius, double kappa=1.0, Selector constitCut=SelectorPtMin(0.)) : _alpha(alpha), _radius(jet_radius), _kappa(kappa), _constitCut(constitCut) {}

  Angularity() = default;
  /// description
  std::string description() const{
    std::ostringstream oss;
    oss << "Angularity with alpha=" << _alpha;
    return oss.str();
  }

  /// computation of the angularity itself
  //double result(const PseudoJet &jet) const{
  double operator()(const PseudoJet &jet) const{
    // get the jet constituents
    std::vector<PseudoJet> constits = jet.constituents();

    // get the reference axis
    PseudoJet reference_axis = _get_reference_axis(jet);

    // do the actual coputation
    double numerator = 0.0, denominator = 0.0;
    unsigned int num = 0;
    for (const auto &c : constits){
      if (!_constitCut.pass(c)) continue;
      double pt = c.pt();
      // Note: better compute (dist^2)^(alpha/2) to avoid an extra square root
      numerator   += pow(pt, _kappa) * pow(c.squared_distance(reference_axis), 0.5*_alpha);
      denominator += pt;
      num += 1;
    }
    if (denominator == 0) return -1;
    // the formula is only correct for the the typical angularities which satisfy either kappa==1 or alpha==0.
    else return numerator/(pow(denominator, _kappa)*pow(_radius, _alpha));
  }

protected:
  PseudoJet _get_reference_axis(const PseudoJet &jet) const{
    if (_alpha>1) return jet;

    fastjet::Recluster recluster(JetDefinition(antikt_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme));
    return recluster(jet);
  }

  double _alpha, _radius, _kappa;
  Selector _constitCut;
};


namespace Rivet {


  class FinalSPPartons : public FinalState {
    /// Final state accessing the intermediate partonic finals state
    /// that can be traced back to the hard process
    /// based on code form arXiv:2012.09574, arXiv:2112.09545
  public:

    /// Constructor
    FinalSPPartons(const Cut& c=Cuts::open())
       : FinalState(c) { }

    /// Clone method
    DEFAULT_RIVET_PROJ_CLONE(FinalSPPartons);

    /// Do the calculation
    void project(const Event& e) override;

  protected:
    /// Cut-applying method overload
    bool accept(const Particle& p) const override;
  };


  class FlavAlgAnalysis : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(FlavAlgAnalysis);

    /// Book histograms and initialise projections before the run
    void init() {

#ifdef DebugLog
      // set optionally the verbosity for the internal Rivet message system
      getLog().setLevel(0);
#endif

      if ( getOption("ALG") == "IFN" )       flavAlg = IFN;
      else if ( getOption("ALG") == "CMP" )  flavAlg = CMP;
      else if ( getOption("ALG") == "CMP2" ) flavAlg = CMP2;
      else if ( getOption("ALG") == "GHS" )  flavAlg = GHS;
      else if ( getOption("ALG") == "SDF" )  flavAlg = SDF;
      else if ( getOption("ALG") == "AKT" )  flavAlg = AKT;
      else if ( getOption("ALG") == "TAG" )  flavAlg = TAG;
      else if ( getOption("ALG") == "OTAG" ) flavAlg = OTAG;
      else if ( getOption("ALG") == "CONE" ) flavAlg = CONE;
      else {
        cout<<"unkown flavour algorithm '"<<getOption("ALG")<<"'.";
        exit(1);
      }

      if ( getOption("ALG2") == "IFN" )       flavAlg2 = IFN;
      else if ( getOption("ALG2") == "CMP" )  flavAlg2 = CMP;
      else if ( getOption("ALG2") == "CMP2" ) flavAlg2 = CMP2;
      else if ( getOption("ALG2") == "GHS" )  flavAlg2 = GHS;
      else if ( getOption("ALG2") == "SDF" )  flavAlg2 = SDF;
      else if ( getOption("ALG2") == "AKT" )  flavAlg2 = AKT;
      else if ( getOption("ALG2") == "TAG" )  flavAlg2 = TAG;
      else if ( getOption("ALG2") == "OTAG" ) flavAlg2 = OTAG;
      else if ( getOption("ALG2") == "CONE" ) flavAlg2 = CONE;
      else                                    flavAlg2 = NONE;

      if ( getOption("TAGPID") == "4" ) tagPID = 4;
      else if ( getOption("TAGPID") == "5" ) tagPID = 5;
      else tagPID = 5; // b-jet tagging is the default

      if ( getOption("DEBUG") == "1" ) debug = true;
      else debug = false;

      FinalState fs; ///< @todo No cuts?

      HeavyHadrons HHs(Cuts::pT > 5*GeV);
      declare(HHs, "HeavyHadrons");

      ZFinder zeeFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zmumuFinder, "ZmumuFinder");

      double R = 0.5;

      #ifdef hepmc3
        cout << " THAT WORKED " << endl;
        if( getOption("LEVEL","HADRON") == "HADRON") {
          VetoedFinalState jetConstits = VetoedFinalState(VisibleFinalState(fs));
          jetConstits.addVetoOnThisFinalState(zeeFinder);
          jetConstits.addVetoOnThisFinalState(zmumuFinder);
          declare(jetConstits, "jetConstits");
          FastJets akt05Jets(jetConstits, FastJets::ANTIKT, R);
          declare(akt05Jets, "AntiKt05Jets");
        }
        else {
          VetoedFinalState jetConstits = VetoedFinalState(FinalSPPartons());
          jetConstits.addVetoOnThisFinalState(zeeFinder);
          jetConstits.addVetoOnThisFinalState(zmumuFinder);
          declare(jetConstits, "jetConstits");
          FastJets akt05Jets(jetConstits, FastJets::ANTIKT, R);
          declare(akt05Jets, "AntiKt05Jets");
        }
      #else
        VetoedFinalState jetConstits = VetoedFinalState(VisibleFinalState(fs));
        jetConstits.addVetoOnThisFinalState(zeeFinder);
        jetConstits.addVetoOnThisFinalState(zmumuFinder);
        declare(jetConstits, "jetConstits");
        FastJets akt05Jets(jetConstits, FastJets::ANTIKT, R);
        declare(akt05Jets, "AntiKt05Jets");
      #endif

      // we start with a base jet definition (should be either
      // antikt_algorithm or cambridge_algorithm, or their e+e- variants)
      base_jet_def = JetDefinition(antikt_algorithm, R);
      // enable it to track flavours (default is net flavour)
      base_jet_def.set_recombiner(&flav_recombiner);

      if ( flavAlg == IFN ) {
        // And then we set up the IFN_Plugin that builds on the base_jet_def
        // The main free parameter, alpha, in the uij distance,
        //   uij = max(pt_i, pt_j)^alpha min(pt_i, pt_j)^(2-alpha) Omega_ij
        double alpha = 2.0;
        // The parameter that sets the nature of the Omega rapidity term;
        // only change the default of 3-alpha if you are sure you know what you are doing
        double omega = 3.0 - alpha;
        // The flavour summation scheme; should be one of
        //   - FlavRecombiner::net
        //   - FlavRecombiner::modulo_2
        FlavRecombiner::FlavSummation flav_summation = FlavRecombiner::net;
        // then construct the IFNPlugin jet definition
        flav_jet_def = JetDefinition(new IFNPlugin(base_jet_def, alpha, omega, flav_summation));
        flav_jet_def.delete_plugin_when_unused();
      }
      else if ( flavAlg == CMP || flavAlg == CMP2 ) {
        // CMP parameters:
        // CMP 'a' parameter in
        //   kappa_ij = 1/a * (kT_i^2 + kT_j^2) / (2*kT_max^2)
        double CMP_a = 0.1;
        // correction to original CMP algo: do not change this if you want IRC safety!
        CMPPlugin::CorrectionType CMP_corr = CMPPlugin::CorrectionType::OverAllCoshyCosPhi_a2;
        // Dynamic definition of ktmax
        CMPPlugin::ClusteringType CMP_clust = CMPPlugin::ClusteringType::DynamicKtMax;
        // CMP plugin
        flav_jet_def = JetDefinition(new CMPPlugin(R, CMP_a, CMP_corr, CMP_clust));
        // enable it to track flavours (default is net flavour)
        flav_jet_def.set_recombiner(&flav_recombiner);
        flav_jet_def.delete_plugin_when_unused();
      }
      else if ( flavAlg == GHS ) {
        GHS_alpha = 1.0; // < flav-kt distance parameter alpha
        GHS_omega = 2.0; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)
        GHS_ptcut = 15.0; // < overall ptcut
      }
      else if ( flavAlg == SDF ) {
        double zcut = 0.1;
        double beta = 1;
        sdFlavCalc = SDFlavourCalc(beta,zcut,R);
      }
      else if ( flavAlg == AKT ) {
        flav_jet_def = base_jet_def;
      }

      if ( flavAlg2 == IFN ) {
        // And then we set up the IFN_Plugin that builds on the base_jet_def
        // The main free parameter, alpha, in the uij distance,
        //   uij = max(pt_i, pt_j)^alpha min(pt_i, pt_j)^(2-alpha) Omega_ij
        double alpha = 2.0;
        // The parameter that sets the nature of the Omega rapidity term;
        // only change the default of 3-alpha if you are sure you know what you are doing
        double omega = 3.0 - alpha;
        // The flavour summation scheme; should be one of
        //   - FlavRecombiner::net
        //   - FlavRecombiner::modulo_2
        FlavRecombiner::FlavSummation flav_summation = FlavRecombiner::net;
        // then construct the IFNPlugin jet definition
        flav_jet_def2= JetDefinition(new IFNPlugin(base_jet_def, alpha, omega, flav_summation));
        flav_jet_def2.delete_plugin_when_unused();
      }
      else if ( flavAlg2 == CMP || flavAlg2 == CMP2) {
        // CMP parameters:
        // CMP 'a' parameter in
        //   kappa_ij = 1/a * (kT_i^2 + kT_j^2) / (2*kT_max^2)
        double CMP_a = 0.1;
        // correction to original CMP alg2o: do not change this if you want IRC safety!
        CMPPlugin::CorrectionType CMP_corr = CMPPlugin::CorrectionType::OverAllCoshyCosPhi_a2;
        // Dynamic definition of ktmax
        CMPPlugin::ClusteringType CMP_clust = CMPPlugin::ClusteringType::DynamicKtMax;
        // CMP plugin
        flav_jet_def2 = JetDefinition(new CMPPlugin(R, CMP_a, CMP_corr, CMP_clust));
        // enable it to track flavours (default is net flavour)
        flav_jet_def2.set_recombiner(&flav_recombiner);
        flav_jet_def2.delete_plugin_when_unused();
      }
      else if ( flavAlg2 == GHS ) {
        GHS_alpha = 1.0; // < flav-kt distance parameter alpha
        GHS_omega = 2.0; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)
        GHS_ptcut = 15.0; // < overall ptcut
      }
      else if ( flavAlg2 == SDF ) {
        double zcut = 0.1;
        double beta = 1;
        sdFlavCalc = SDFlavourCalc(beta,zcut,R);
      }
      else if ( flavAlg2 == AKT ) {
        flav_jet_def2 = base_jet_def;
      }

      //Histograms booking

      book(_h_fidxsec,"fidxsec",1,0.,1.);
      book(_h_fidxsec_b,"fidxsec_b",1,0.,1.);

      book(_h_first_bjet_pt_b ,1,1,1);
      book(_h_first_bjet_abseta_b ,3,1,1);
      book(_h_Z_pt_b ,5,1,1);
      book(_h_HT_b ,7,1,1);
      book(_h_Dphi_Zb_b ,9,1,1);

      book(_h_first_jet_pt_ratio ,2,1,1);
      book(_h_first_jet_abseta_ratio ,4,1,1);
      book(_h_Z_pt_ratio ,6,1,1);
      book(_h_HT_ratio ,8,1,1);
      book(_h_Dphi_Zj_ratio ,10,1,1);

      book(_h_first_jet_pt, "first_jet_pt", refData(1,1,1) ); // (*_h_first_bjet_pt_b);
      book(_h_first_jet_abseta, "first_jet_abseta", refData(3,1,1) ); // (*_h_first_bjet_abseta_b);
      book(_h_Z_pt, "Z_pt", refData(5,1,1) ); // (*_h_Z_pt_b);
      book(_h_HT, "HT", refData(7,1,1) ); // (*_h_HT_b);
      book(_h_Dphi_Zj, "Dphi_Zj", refData(9,1,1) ); // (*_h_Dphi_Zb_b);

      book(_h_first_bjet_pt_bb ,11,1,1);
      book(_h_second_bjet_pt_bb ,12,1,1);
      book(_h_Z_pt_bb ,13,1,1);
      book(_h_bb_mass_bb ,14,1,1);
      book(_h_Zbb_mass_bb ,15,1,1);
      book(_h_Dphi_bb ,16,1,1);
      book(_h_DR_bb ,17,1,1);
      book(_h_DR_Zbmin_bb ,18,1,1);
      book(_h_A_DR_Zb_bb ,19,1,1);

      book(_h_bjet_multiplicity ,20,1,1);

      book(_h_compare,"compare_algs",6,-0.5,5.5);

      book(_h_first_jet_pt_by_cat[LEAD_TAGGED_BOTH], "jet_pt_LEAD_TAGGED_BOTH", refData(1,1,1) );
      book(_h_first_jet_pt_by_cat[LEAD_TAGGED_ALG1], "jet_pt_LEAD_TAGGED_ALG1", refData(1,1,1) );
      book(_h_first_jet_pt_by_cat[LEAD_TAGGED_ALG2], "jet_pt_LEAD_TAGGED_ALG2", refData(1,1,1) );
      book(_h_first_jet_pt_by_cat[SUBLEAD_TAG_AGREE], "jet_pt_LEAD_TAG_AGREE", refData(1,1,1) );
      book(_h_first_jet_pt_by_cat[SUBLEAD_TAG_DISAGREE], "jet_pt_LEAD_TAG_DISAGREE", refData(1,1,1) );
      book(_h_first_jet_pt_by_cat[NO_TAG_BOTH], "jet_pt_NO_TAG_BOTH", refData(1,1,1) );

      book(_h_first_bjet_pt_by_cat[LEAD_TAGGED_BOTH], "bjet_pt_LEAD_TAGGED_BOTH", refData(1,1,1) );
      book(_h_first_bjet_pt_by_cat[LEAD_TAGGED_ALG1], "bjet_pt_LEAD_TAGGED_ALG1", refData(1,1,1) );
      book(_h_first_bjet_pt_by_cat[LEAD_TAGGED_ALG2], "bjet_pt_LEAD_TAGGED_ALG2", refData(1,1,1) );
      book(_h_first_bjet_pt_by_cat[SUBLEAD_TAG_AGREE], "bjet_pt_LEAD_TAG_AGREE", refData(1,1,1) );
      book(_h_first_bjet_pt_by_cat[SUBLEAD_TAG_DISAGREE], "bjet_pt_LEAD_TAG_DISAGREE", refData(1,1,1) );
      book(_h_first_bjet_pt_by_cat[NO_TAG_BOTH], "bjet_pt_NO_TAG_BOTH", refData(1,1,1) );

      book( _h_ang05, "ang05",  20, 0, 1);
      book( _h_ang10, "ang10",  20, 0, 1);
      book( _h_ang20, "ang20",  20, 0, 1);
      book( _h_mass, "mass",  20, 0, 1);
      book( _h_ang05_b, "ang05_b",  20, 0, 1);
      book( _h_ang10_b, "ang10_b",  20, 0, 1);
      book( _h_ang20_b, "ang20_b",  20, 0, 1);
      book( _h_mass_b, "mass_b",  20, 0, 1);

      book(_h_bbcorr, "bb_correlations",     5,0,5, 5,0,0.5);
      book(_h_bbcorr_2, "bb_correlations_2", 5,0,5, 5,0,0.5);

      book(_h_pTb1_highPT, "bjet_pT_highPT",  40, 200., 1200.);
      book(_h_pTj1_highPT, "jet_pT_highPT" , 40, 200., 1200.);

      book(_h_pTb1_nb , "pTb1_nb" , 10,200.,1200.,20,0,20.);
      book(_h_pTb1_nb2, "pTb1_nb2", 10,200.,1200.,20,0,20.);

      ang05 = Angularity(0.5, R);
      ang10 = Angularity(1.0, R);
      ang20 = Angularity(2.0, R);

      if (debug) {
        std::cout << "Jet descr \n";
        std::cout << "jet def descr = \n" << flav_jet_def.description() << "\n";
        std::cout << "jet def descr = \n" << base_jet_def.description() << "\n";
      }

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zeeFS = applyProjection<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmumuFS = applyProjection<ZFinder>(event, "ZmumuFinder");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      // We did not find exactly one Z. No good.
      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      //event identification depending on mass window
      bool ee_event=false;
      bool mm_event=false;

      if (zees.size() == 1) { ee_event = true; }
      if (zmumus.size() == 1) { mm_event = true; }
      const Particles& theLeptons = zees.size() ? zeeFS.constituents() : zmumuFS.constituents();

      if (debug) {
        std::cout<<"~~~~~~~~~~~~Event in rivet \n";
        for (unsigned int i=0; i<event.allParticles().size(); i++) {
          std::cout<<event.allParticles()[i]<<std::endl;
        }
      }

      Jets goodjets;
      Jets jb_final;
      Jets goodjets2;
      Jets jb_final2;
      bool alg1_lead_is_btagged = false;
      bool alg2_lead_is_btagged = false;
      int alg1_first_btagged = -1;
      int alg2_first_btagged = -1;
      FourMomentum bTagMom;
      bool bTagInJet = false;

      if (flavAlg != TAG && flavAlg != OTAG && flavAlg != CONE) { // Flavoured jet algorithms

        // NB. Veto has already been applied on leptons and photons used for dressing

        const FinalState& jetConstits_flav= applyProjection<FinalState>(event, "jetConstits");

        if (debug) {
          std::cout<<"~~~~~~~~~~~~~Projection \n";
          for (unsigned int i=0; i< jetConstits_flav.particles().size(); i++) {
            std::cout << jetConstits_flav.particles()[i] << std::endl;
          }
        }

        PseudoJets fj_flav = FastJets::mkClusterInputs(jetConstits_flav.particles());
        for (unsigned int i=0; i<  fj_flav.size(); i++) {
          if (debug) std::cout<<fj_flav[i].description()<<"\n";
          const int pdgid = jetConstits_flav.particles()[i].pid();
          fastjet::contrib::FlavInfo flav_info_init(pdgid);
          if (flavAlg == CMP2) flav_info_init.reset_all_but_flav(tagPID);
          fj_flav[i].set_user_info(new fastjet::contrib::FlavHistory(const_cast<const fastjet::contrib::FlavInfo&>(flav_info_init)));
        }

        if (debug) {
          std::cout<<"Convert FS into pseudojets \n";
          std::cout<<"~~~~~~~~~ FS \n";

          for (unsigned int i=0; i<  fj_flav.size(); i++) {
            std::cout<<fj_flav[i].description()<<"\n";
            std::cout<<"pseudo jet rap="<<fj_flav[i].rap()<<" pT="<<fj_flav[i].perp() <<" flav= "<< FlavHistory::current_flavour_of(fj_flav[i]).description()<<"\n";
          }

          std::cout<<"user info flav done \n";
          std::cout<<"jet def descr = \n"<< flav_jet_def.description()<<"\n";
          std::cout<<"jet def descr = \n"<< base_jet_def.description()<<"\n";
        }

        vector<PseudoJet> base_jets = sorted_by_pt(base_jet_def(fj_flav));
        vector<PseudoJet> flav_pseudojets;
        if (flavAlg == IFN || flavAlg == CMP || flavAlg == CMP2) {
          flav_pseudojets = sorted_by_pt(flav_jet_def(fj_flav));
        }
        else if (flavAlg == GHS) {
          flav_pseudojets = run_GHS(base_jets, GHS_ptcut,
                                    GHS_alpha, GHS_omega, flav_recombiner);
        }
        else if (flavAlg == SDF) {
          flav_pseudojets = base_jets;
          sdFlavCalc(flav_pseudojets);
        }
        else if (flavAlg == AKT) {
          flav_pseudojets = base_jets;
        }

        const Jets& jets = FastJets::mkJets(flav_pseudojets, jetConstits_flav.particles());

        if (debug) {
          std::cout<<"The set of clustered pseudo jets \n";

            for (const PseudoJet& j: flav_pseudojets) {
            std::cout<<"jet rap="<<j.rap()<<" pT="<<j.perp() <<" flav= "<< FlavHistory::current_flavour_of(j).description()<<"\n";
          }
          std::cout<<"The set of clustered jets \n";

            for (const Jet& j: jets) {
            std::cout<<"jet rap="<<j.rap()<<" pT="<<j.perp() <<" flav= "<< FlavHistory::current_flavour_of(j).description()<<"\n";
          }
        }
        for (const Jet& j: jets) {
          if (j.perp()>30 && std::abs(j.eta())<2.4) goodjets.push_back(j);
        }
        if (debug) {
          std::cout<<"The set of good jets \n";

            for (const Jet& j: goodjets) {
            std::cout<<"jet rap="<<j.rap()<<" pT="<<j.perp() <<" flav= "<< FlavHistory::current_flavour_of(j).description()<<"\n";
          }
        }

        //identification of bjets/cjets
        for (unsigned int i=0; i<goodjets.size(); i++ ) {
          const bool btagged =  std::abs(FlavHistory::current_flavour_of(goodjets[i])[tagPID]%2) == 1;
          if (btagged) {
            if(jb_final.size() == 0) {
              for(const Particle& j: goodjets[i].constituents()) {
                if(std::abs(j.pid()) == tagPID) {
                  if(!bTagInJet || j.pT() > bTagMom.pT()) bTagMom = j.momentum();
                  bTagInJet = true;
                }
              }
            }
            jb_final.push_back(goodjets[i]);
            if (alg1_first_btagged < 0) alg1_first_btagged = i;
          }
          if (i==0) alg1_lead_is_btagged = btagged;
      }
        if (debug) {
          std::cout<<"The set of good b jets \n";

            for (const Jet& j: jb_final) {
            std::cout<<"jet rap="<<j.rap()<<" pT="<<j.perp() <<" flav= "<< FlavHistory::current_flavour_of(j).description()<<"\n";
          }
        }
    } else {

        const FastJets fj = applyProjection<FastJets>(event, "AntiKt05Jets");
        goodjets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);

        //ATLAS STYLE TRUTH TAGGING
        if (flavAlg == CONE) {

          const HeavyHadrons& HHs = applyProjection<HeavyHadrons>(event, "HeavyHadrons");
          Particles bHadrons;
          if (tagPID == 5) bHadrons = HHs.bHadrons(Cuts::pT > 5*GeV);
          else if (tagPID == 4) bHadrons = HHs.cHadrons(Cuts::pT > 5*GeV);

          Particles matchedBs;

          for (const Jet& j : goodjets) {
            Jet closest_j;
            Particle closest_b;
            double minDR_j_b = 10;

            for (const Particle& b : bHadrons) {
              bool alreadyMatched = false;

              for (const Particle& matchedB : matchedBs) {
                alreadyMatched = matchedB.isSame(b);
              }
              if(alreadyMatched) continue;
              double DR_j_b = deltaR(j, b);

              if (DR_j_b < 0.3 && DR_j_b < minDR_j_b) {
                minDR_j_b = DR_j_b;
                closest_j = j;
                closest_b = b;
              }
            }
            if (minDR_j_b < 0.3) {
              jb_final.push_back(closest_j);
              matchedBs.push_back(closest_b);
            }
          }

        } else if (flavAlg==TAG) {
          //CMS STYLE TAGGING
          for (const Jet& j : goodjets) {
            if (tagPID == 5 && j.bTagged() ) { jb_final.push_back(j); }
            if (tagPID == 4 && j.cTagged() ) { jb_final.push_back(j); }
          }
        } else if (flavAlg==OTAG) {
          //CMS STYLE TAGGING, but requiring an odd number of btags
          for (const Jet& j : goodjets) {
            if( tagPID == 5 && j.bTagged()){
              const int btags = j.bTags().size();
              if(btags%2 == 1)  jb_final.push_back(j);
            }
            if( tagPID == 4 && j.cTagged()){
              const int ctags = j.cTags().size();
              if(ctags%2 == 1)  jb_final.push_back(j);
            }
          }
        }
      }

      // Algorithm-algorithm correlations
      if(flavAlg2 != NONE && flavAlg2 != TAG && flavAlg2 != OTAG && flavAlg2 != CONE){

        // NB. Veto has already been applied on leptons and photons used for dressing

        const FinalState& jetConstits_flav= applyProjection<FinalState>(event, "jetConstits");

        PseudoJets fj_flav = FastJets::mkClusterInputs(jetConstits_flav.particles());
        for (unsigned int i=0; i<  fj_flav.size(); i++) {
          if (debug) std::cout << fj_flav[i].description() << "\n";
          const int pdgid = jetConstits_flav.particles()[i].pid();
          fj_flav[i].set_user_info(new  fastjet::contrib::FlavHistory(pdgid));
        }

        vector<PseudoJet> base_jets = sorted_by_pt(base_jet_def(fj_flav));
        vector<PseudoJet> flav_pseudojets;

        if (flavAlg2 == IFN || flavAlg2 == CMP || flavAlg2 == CMP2) {
          flav_pseudojets = sorted_by_pt(flav_jet_def2(fj_flav));
        }
        else if (flavAlg2 == GHS) {
          flav_pseudojets = run_GHS(base_jets, GHS_ptcut,
                                    GHS_alpha, GHS_omega, flav_recombiner);
        }
        else if (flavAlg2 == SDF) {
          flav_pseudojets = base_jets;
          sdFlavCalc(flav_pseudojets);
        }
        else if (flavAlg2 == AKT) {
          flav_pseudojets = base_jets;
        }

        const Jets& jets = FastJets::mkJets(flav_pseudojets, jetConstits_flav.particles());

        for( auto j: jets){
          if (j.perp()>30 && std::abs(j.eta())<2.4) goodjets2.push_back(j);
        }


        //identification of bjets
        for (unsigned int i=0; i<goodjets2.size(); i++ ) {
          const bool btagged =  std::abs(FlavHistory::current_flavour_of(goodjets2[i])[tagPID]%2) ==1 ;
          if (btagged) {
            jb_final2.push_back(goodjets2[i]);
            if(alg2_first_btagged < 0) alg2_first_btagged = i;
          }
          if(i==0) alg2_lead_is_btagged = btagged;
        }
      } /// at this point we have the flav alg to compare too

      int cat = -1;
      if (goodjets.size() != 0 && goodjets2.size() != 0) {
        if (alg1_lead_is_btagged && alg2_lead_is_btagged)   cat = LEAD_TAGGED_BOTH;   //_h_compare->fill(LEAD_TAGGED_BOTH);
        else if (alg1_lead_is_btagged)                    cat = LEAD_TAGGED_ALG1; //_h_compare->fill(LEAD_TAGGED_ALG1);
        else if (alg2_lead_is_btagged)                    cat = LEAD_TAGGED_ALG2; //_h_compare->fill(LEAD_TAGGED_ALG2);
        else {
          if(alg1_first_btagged < 0 && alg2_first_btagged < 0) cat = NO_TAG_BOTH; // _h_compare->fill(NO_TAG_BOTH);
          else if(alg1_first_btagged == alg2_first_btagged)    cat = SUBLEAD_TAG_AGREE; //_h_compare->fill(SUBLEAD_TAG_AGREE);
          else                                                 cat = SUBLEAD_TAG_DISAGREE; //_h_compare->fill(SUBLEAD_TAG_DISAGREE);
        }
      }

      //Event weight
      const double w = 0.5;
      _h_compare->fill(cat,w);

      const HeavyHadrons& HHs = applyProjection<HeavyHadrons>(event, "HeavyHadrons");
      unsigned nb_event = 0;
      if (tagPID == 5) for (auto &el : HHs.bHadrons()) nb_event++;
      else if (tagPID == 4) for (auto &el : HHs.cHadrons()) nb_event++;

      Particles tagIdHadrons;
      if (tagPID == 5) tagIdHadrons = HHs.bHadrons(Cuts::pT > 5*GeV);
      else if (tagPID == 4) tagIdHadrons = HHs.cHadrons(Cuts::pT > 5*GeV);

      //histogram filling
      // double z_pt;
      // if ( ee_event ) z_pt = zees[0].pt();
      // if ( mm_event ) z_pt = zmumus[0].pt();
      // Do we need a z_pt > 20*GeV cut here?

      if ((ee_event || mm_event) && goodjets.size() > 0) {
        _h_fidxsec->fill(0.5,w);
        double Ht = 0;

        for (const Jet& j : goodjets) {
          Ht += j.pT();
        }

        FourMomentum j1(goodjets[0].momentum());

        _h_first_jet_pt->fill(j1.pt(),w);
        _h_pTj1_highPT->fill(j1.pt(),w);
        _h_first_jet_pt_by_cat[cat]->fill(j1.pt(),w);
        _h_first_jet_abseta->fill(fabs(j1.eta()),w);
        if ( ee_event ) _h_Z_pt->fill(zees[0].pt(),w);
        if ( mm_event ) _h_Z_pt->fill(zmumus[0].pt(),w);
        _h_HT->fill(Ht,w);
        if ( ee_event ) _h_Dphi_Zj->fill(deltaPhi(zees[0], j1),w);
        if ( mm_event ) _h_Dphi_Zj->fill(deltaPhi(zmumus[0], j1),w);

        _h_ang05->fill(ang05(goodjets[0].pseudojet()),w);
        _h_ang10->fill(ang10(goodjets[0].pseudojet()),w);
        _h_ang20->fill(ang20(goodjets[0].pseudojet()),w);
        _h_mass->fill(j1.mass()/j1.pt(),w);

        if ( jb_final.size() > 0 ) {
          _h_fidxsec_b->fill(0.5,w);

          FourMomentum b1(jb_final[0].momentum());

          _h_bjet_multiplicity->fill(1.,w);

          _h_first_bjet_pt_b->fill(b1.pt(),w);
          _h_pTb1_highPT->fill(b1.pt(),w);
          _h_first_bjet_pt_by_cat[cat]->fill(j1.pt(),w);
          _h_first_bjet_abseta_b->fill(fabs(b1.eta()),w);
          if ( ee_event ) _h_Z_pt_b->fill(zees[0].pt(),w);
          if ( mm_event ) _h_Z_pt_b->fill(zmumus[0].pt(),w);
          _h_HT_b->fill(Ht,w);
          if ( ee_event ) _h_Dphi_Zb_b->fill(deltaPhi(zees[0], b1.phi()),w);
          if ( mm_event ) _h_Dphi_Zb_b->fill(deltaPhi(zmumus[0], b1.phi()),w);

          _h_ang05_b->fill(ang05(jb_final[0].pseudojet()),w);
          _h_ang10_b->fill(ang10(jb_final[0].pseudojet()),w);
          _h_ang20_b->fill(ang20(jb_final[0].pseudojet()),w);
          _h_mass_b->fill(b1.mass()/b1.pt(),w);

          if(bTagInJet && tagIdHadrons.size() == 2) {
            const double dR = deltaR(tagIdHadrons[0],tagIdHadrons[1]);
            const double relPt = bTagMom.pT()/b1.pT();
            _h_bbcorr->fill(dR,relPt,w);
            if (alg1_lead_is_btagged && !alg2_lead_is_btagged)
              _h_bbcorr_2->fill(dR,relPt,w);
          }

          _h_pTb1_nb->fill(b1.pT(),nb_event,w);
          unsigned nb_bj1 = 0;
          for (const auto &el : jb_final[0].constituents())
            if (tagPID ==5 && el.hasBottom()) nb_bj1++;
            else if (tagPID == 4 && el.hasCharm()) nb_bj1++;

          _h_pTb1_nb2->fill(b1.pT(),nb_bj1,w);

          if ( jb_final.size() > 1 ) {

            FourMomentum b2(jb_final[1].momentum());

            _h_bjet_multiplicity->fill(2.,w);

            _h_first_bjet_pt_bb->fill(b1.pt(),w);
            _h_second_bjet_pt_bb->fill(b2.pt(),w);
            if ( ee_event ) _h_Z_pt_bb->fill(zees[0].pt(),w);
            if ( mm_event ) _h_Z_pt_bb->fill(zmumus[0].pt(),w);

            FourMomentum bb = add(b1,b2);
            FourMomentum Zbb;
            if (ee_event) Zbb = add(zees[0],bb);
            if (mm_event) Zbb = add(zmumus[0],bb);

            _h_bb_mass_bb->fill(bb.mass(),w);
            _h_Zbb_mass_bb->fill(Zbb.mass(),w);

            _h_Dphi_bb->fill(deltaPhi(b1,b2),w);
            if (deltaR(b1,b2)>0.5) {
              _h_DR_bb->fill(deltaR(b1,b2),w);
            }

            double DR_Z_b1(0.), DR_Z_b2(0.);
            if ( ee_event ) {
              DR_Z_b1 = deltaR(zees[0],b1);
              DR_Z_b2 = deltaR(zees[0],b2);
            }
            if ( mm_event ) {
              DR_Z_b1 = deltaR(zmumus[0],b1);
              DR_Z_b2 = deltaR(zmumus[0],b2);
            }

            double DR_Zb_min = DR_Z_b1;
            double DR_Zb_max = DR_Z_b2;
            if ( DR_Zb_min > DR_Zb_max ) {
              DR_Zb_min = DR_Z_b2;
              DR_Zb_max = DR_Z_b1;
            }
            double A_Zbb = (DR_Zb_max - DR_Zb_min)/(DR_Zb_max + DR_Zb_min);

            _h_DR_Zbmin_bb->fill(DR_Zb_min,w);
            _h_A_DR_Zb_bb->fill(A_Zbb,w);

          }

        }

      }

    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/picobarn/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      scale( _h_first_bjet_pt_b, 100. );
      scale( _h_first_bjet_abseta_b, 100. );
      scale( _h_Z_pt_b, 100. );
      scale( _h_HT_b, 100. );
      scale( _h_Dphi_Zb_b, 100. );

      divide( _h_first_bjet_pt_b , _h_first_jet_pt , _h_first_jet_pt_ratio );
      divide( _h_first_bjet_abseta_b , _h_first_jet_abseta , _h_first_jet_abseta_ratio );
      divide( _h_Z_pt_b , _h_Z_pt , _h_Z_pt_ratio );
      divide( _h_HT_b , _h_HT , _h_HT_ratio );
      divide( _h_Dphi_Zb_b , _h_Dphi_Zj , _h_Dphi_Zj_ratio );

      scale( _h_first_bjet_pt_b, norm/100. );
      scale( _h_first_bjet_abseta_b, norm/100. );
      scale( _h_Z_pt_b, norm/100. );
      scale( _h_HT_b, norm/100. );
      scale( _h_Dphi_Zb_b, norm/100. );

      scale( _h_fidxsec, norm);
      scale( _h_fidxsec_b, norm);

      scale( _h_pTb1_highPT, norm);
      scale( _h_pTj1_highPT, norm);

      scale( _h_first_jet_pt, norm);
      scale( _h_first_jet_abseta, norm);
      scale( _h_Z_pt, norm);
      scale( _h_HT, norm);
      scale( _h_Dphi_Zj, norm);

      scale( _h_first_bjet_pt_bb, norm);
      scale( _h_second_bjet_pt_bb, norm);
      scale( _h_Z_pt_bb, norm);
      scale( _h_bb_mass_bb, norm);
      scale( _h_Zbb_mass_bb, norm);
      scale( _h_Dphi_bb, norm);
      scale( _h_DR_bb, norm);
      scale( _h_DR_Zbmin_bb, norm);
      scale( _h_A_DR_Zb_bb, norm);

      scale( _h_bjet_multiplicity, norm );

      scale( _h_compare, norm);
      scale( _h_ang05, norm);
      scale( _h_ang10, norm);
      scale( _h_ang20, norm);
      scale( _h_mass, norm);
      scale( _h_ang05_b, norm);
      scale( _h_ang10_b, norm);
      scale( _h_ang20_b, norm);
      scale( _h_mass_b, norm);

      scale( _h_bbcorr, norm);
      scale( _h_bbcorr_2, norm);

      scale( _h_pTb1_nb, norm);
      scale( _h_pTb1_nb2, norm);

      for ( auto h : _h_first_jet_pt_by_cat ) {
        scale(h, norm);
      }
      for ( auto h : _h_first_bjet_pt_by_cat ) {
        scale(h, norm);
      }

    }

  private:

    /// @name Histograms

    Histo1DPtr     _h_fidxsec, _h_fidxsec_b;

    Histo1DPtr     _h_first_jet_pt, _h_first_bjet_pt_b, _h_pTj1_highPT, _h_pTb1_highPT;
    Histo1DPtr     _h_first_jet_abseta, _h_first_bjet_abseta_b;
    Histo1DPtr     _h_Z_pt, _h_Z_pt_b;
    Histo1DPtr     _h_HT, _h_HT_b;
    Histo1DPtr     _h_Dphi_Zj, _h_Dphi_Zb_b;

    Scatter2DPtr   _h_first_jet_pt_ratio;
    Scatter2DPtr   _h_first_jet_abseta_ratio;
    Scatter2DPtr   _h_Z_pt_ratio;
    Scatter2DPtr   _h_HT_ratio;
    Scatter2DPtr   _h_Dphi_Zj_ratio;

    Histo1DPtr     _h_first_bjet_pt_bb, _h_second_bjet_pt_bb;
    Histo1DPtr     _h_Z_pt_bb;
    Histo1DPtr     _h_bb_mass_bb, _h_Zbb_mass_bb;
    Histo1DPtr     _h_Dphi_bb, _h_DR_bb, _h_DR_Zbmin_bb, _h_A_DR_Zb_bb;

    Histo1DPtr     _h_bjet_multiplicity;

    Histo1DPtr     _h_ang05, _h_ang10, _h_ang20, _h_mass;
    Histo1DPtr     _h_ang05_b, _h_ang10_b, _h_ang20_b, _h_mass_b;
    Angularity     ang05, ang10, ang20;

    Histo2DPtr     _h_bbcorr;
    Histo2DPtr     _h_bbcorr_2;

    Histo2DPtr     _h_pTb1_nb;
    Histo2DPtr     _h_pTb1_nb2;

    enum {
      LEAD_TAGGED_BOTH = 0,
      LEAD_TAGGED_ALG1 = 1,
      LEAD_TAGGED_ALG2 = 2,
      SUBLEAD_TAG_AGREE = 3,
      SUBLEAD_TAG_DISAGREE = 4,
      NO_TAG_BOTH = 5,
    };
    Histo1DPtr     _h_compare;

    std::array<Histo1DPtr,6>  _h_first_jet_pt_by_cat;
    std::array<Histo1DPtr,6>  _h_first_bjet_pt_by_cat;

    JetDefinition base_jet_def;
    JetDefinition flav_jet_def;
    JetDefinition flav_jet_def2;
    FlavRecombiner flav_recombiner;

    bool debug = false;

    // GHS parameters:
    double GHS_alpha; // < flav-kt distance parameter alpha
    double GHS_omega; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)
    double GHS_ptcut; // < overall ptcut

    SDFlavourCalc sdFlavCalc;

    int flavAlg;
    int flavAlg2;
    enum {
      IFN = 0,
      CMP = 1,
      GHS = 2,
      SDF = 3,
      AKT = 4,
      TAG = 5,
      OTAG = 6,
      CONE = 7,
      NONE = 8,
      CMP2 = 9
    };

    int tagPID = 5;

    std::string flavAlgName() {
      if(flavAlg == IFN) return "IFN";
      else if(flavAlg == CMP) return "CMP";
      else if(flavAlg == CMP2) return "CMP2";
      else if(flavAlg == GHS) return "GHS";
      else if(flavAlg == SDF) return "SDF";
      else if(flavAlg == AKT) return "AKT";
      else if(flavAlg == TAG) return "TAG";
      else if(flavAlg == OTAG) return "OTAG";
      else if(flavAlg == CONE) return "CONE";
      else return "unknown";
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(FlavAlgAnalysis);

  #ifdef hepmc3
    bool FinalSPPartons::accept(const Particle& p) const {
      // Reject if *not* a parton
      if (!isParton(p))
        return false;
      if (p.genParticle()->end_vertex() and p.genParticle()->production_vertex()) {
        // Accept partons if they end on a standard hadronization vertex
        if (p.genParticle()->end_vertex()->status() == 5) {
          auto pv = p.genParticle()->production_vertex();
          // Accept if some ancenstor starts on SP vertex
          for(auto pp: p.ancestors(Cuts::OPEN,false)) {
            if(pp.genParticle() and pp.genParticle()->production_vertex()) {
              if(pp.genParticle()->production_vertex()->status() == 1) {
                return _cuts->accept(p);
              }
            }
          }
        }
      }
      return false;
    }


    void FinalSPPartons::project(const Event& e) {
      _theParticles.clear();
      for (auto gp : HepMCUtils::particles(e.genEvent())) {
        if (!gp) continue;
        const Particle p(gp);
        if (accept(p)) {
          _theParticles.push_back(p);
        }
      }
    }
  #endif
}
