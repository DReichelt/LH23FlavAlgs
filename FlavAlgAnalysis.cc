#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Math/Vector4.hh"
#include "IFNPlugin/IFNPlugin.hh"
#include "CMPPlugin/CMPPlugin.hh"
#include "GHSAlgo/GHSAlgo.hh"
#include "SDFlavPlugin/SDFlavPlugin.hh"

//#define DebugLog

using namespace fastjet;
using namespace fastjet::contrib;

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

      if ( getOption("ALG") == "IFN" ) flavAlg = IFN;
      else if ( getOption("ALG") == "CMP" ) flavAlg = CMP;
      else if ( getOption("ALG") == "GHS" ) flavAlg = GHS;
      else if ( getOption("ALG") == "SDF" ) flavAlg = SDF;
      else if ( getOption("ALG") == "AKT" )  flavAlg = AKT;
      else if ( getOption("ALG") == "TAG" )  flavAlg = TAG;
      else if ( getOption("ALG") == "OTAG" )  flavAlg = OTAG;
      else if ( getOption("ALG") == "CONE" ) flavAlg = CONE;
      else {
        cout<<"unkown flavour algorithm '"<<getOption("ALG")<<"'.";
        exit(1);
      }

      FinalState fs; ///< @todo No cuts?
      
      HeavyHadrons HHs(Cuts::pT > 5*GeV);
      declare(HHs, "HeavyHadrons");

      ZFinder zeeFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71.0*GeV, 111.0*GeV, 0.1 );
      declare(zmumuFinder, "ZmumuFinder");

      double R = 0.5;

      #ifdef hepmc3
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

      if(flavAlg == IFN) {
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
        flav_jet_def= JetDefinition(new IFNPlugin(base_jet_def, alpha, omega, flav_summation));
        flav_jet_def.delete_plugin_when_unused();
      }
      else if(flavAlg == CMP) {
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
      else if(flavAlg == GHS) {
        GHS_alpha = 1.0; // < flav-kt distance parameter alpha
        GHS_beta  = 1.0; // < SoftDrop parameter beta for flavour clusters
        GHS_zcut  = 0.1; // < SoftDrop zcut for flavour clusters
        GHS_Rcut  = 0.1; // < Rcut parameter
        GHS_omega = 0.0; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)
        GHS_ptcut = 15.0; // < overall ptcut
      }
      else if(flavAlg == SDF) {
        double zcut = 0.1;
        double beta = 1;
        sdFlavCalc = SDFlavourCalc(beta,zcut,R);
      }
      else if(flavAlg == AKT) {
        flav_jet_def = base_jet_def;
      }


      //Histograms booking

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



      if(debug){
	std::cout<<"Jet descr \n";
	std::cout<<"jet def descr = \n"<< flav_jet_def.description()<<"\n";
	std::cout<<"jet def descr = \n"<< base_jet_def.description()<<"\n";
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

      if(debug){
	std::cout<<"~~~~~~~~~~~~Event in rivet \n";
	for (unsigned int i=0; i<event.allParticles().size(); i++){
	  std::cout<<event.allParticles()[i]<<std::endl;
	}
      }


      Jets goodjets;
      Jets jb_final;
      double Ht = 0;

      if(flavAlg != TAG && flavAlg != OTAG && flavAlg != CONE){

        // NB. Veto has already been applied on leptons and photons used for dressing

        const FinalState& jetConstits_flav= applyProjection<FinalState>(event, "jetConstits");

        if(debug){
          std::cout<<"~~~~~~~~~~~~~Projection \n";
          for (unsigned int i=0; i< jetConstits_flav.particles().size(); i++){
            std::cout<< jetConstits_flav.particles()[i]<<std::endl;
          }
        }



        PseudoJets fj_flav = FastJets::mkClusterInputs(jetConstits_flav.particles());
        for (unsigned int i=0; i<  fj_flav.size(); i++){
          //	std::cout<<fj_flav[i].description()<<"\n";
          const int pdgid = jetConstits_flav.particles()[i].pid();
          fj_flav[i].set_user_info(new  fastjet::contrib::FlavHistory(pdgid));
        }


        if(debug){
          std::cout<<"Convert FS into pseudojets \n";
          std::cout<<"~~~~~~~~~ FS \n";

          for (unsigned int i=0; i<  fj_flav.size(); i++){
            //	std::cout<<fj_flav[i].description()<<"\n";
            std::cout<<"pseudo jet rap="<<fj_flav[i].rap()<<" pT="<<fj_flav[i].perp() <<" flav= "<< FlavHistory::current_flavour_of(fj_flav[i]).description()<<"\n";
          }

          std::cout<<"user info flav done \n";
          std::cout<<"jet def descr = \n"<< flav_jet_def.description()<<"\n";
          std::cout<<"jet def descr = \n"<< base_jet_def.description()<<"\n";
        }


        vector<PseudoJet> base_jets = sorted_by_pt(base_jet_def(fj_flav));
        vector<PseudoJet> flav_pseudojets;
        if(flavAlg == IFN || flavAlg == CMP) {
          flav_pseudojets = sorted_by_pt(flav_jet_def(fj_flav));
        }
        else if(flavAlg == GHS) {
          flav_pseudojets = run_GHS(base_jets, GHS_ptcut,
                                    GHS_beta, GHS_zcut, GHS_Rcut, GHS_alpha, GHS_omega);
        }
        else if(flavAlg == SDF) {
          flav_pseudojets = base_jets;
          sdFlavCalc(flav_pseudojets);
        }
        else if(flavAlg == AKT) {
          flav_pseudojets = base_jets;
        }

        const Jets& jets = FastJets::mkJets(flav_pseudojets, jetConstits_flav.particles());
        //const Jets& jets= jets_unordered.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30);

        // Perform lepton-jet overlap and HT calculation




        for( auto j: jets){
          if(j.perp()>30 && std::abs(j.eta())<2.4) goodjets.push_back(j);
        }


        //identification of bjets
        for (unsigned int i=0; i<goodjets.size(); i++ ) {
          Ht += goodjets[i].pT();
          const bool btagged =  std::abs(FlavHistory::current_flavour_of(goodjets[i])[5]%2) ==1 ;
          if (btagged) jb_final.push_back(goodjets[i]);
          //if ( j.bTagged() ) { jb_final.push_back(j); }
          //if  FlavHistory::current_flavour_of(jets[i]);
      }
    }else{

	const FastJets fj = applyProjection<FastJets>(event, "AntiKt05Jets");
	goodjets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);
	
	
	//ATLAS STYLE TRUTH TAGGING
	if(flavAlg == CONE) {
        
	  const HeavyHadrons& HHs = applyProjection<HeavyHadrons>(event, "HeavyHadrons");
	  const Particles& bHadrons = HHs.bHadrons(Cuts::pT > 5*GeV);
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
	
	}else if(flavAlg==TAG){ 
	  //CMS STYLE TAGGING
	  for (const Jet& j : goodjets) {
	    if ( j.bTagged() ) { jb_final.push_back(j); }
	  }
	}else if(flavAlg==OTAG){ 
	  //CMS STYLE TAGGING, but requiring an odd number of btags
	  for (const Jet& j : goodjets) {
	    if( j.bTagged()){
	      const int btags = j.bTags().size();
	      if(btags%2 ==1)  jb_final.push_back(j);
	    }
	  }
	}
	
      }
       //cout<<flavAlgName()<<" "<<"goodjets: "<<goodjets.size()<<", btagged: "<<jb_final.size()<<"\n";


      // if(jb_final.size() >0){

      //   std::cout<<"Found bjet!\n";
      //   for (unsigned int i=0; i<event.allParticles().size(); i++){
      //     std::cout<<event.allParticles()[i]<<std::endl;
      //   }
      // }

      for (const Jet& j : goodjets) {
	      Ht += j.pT();
      }




      //Event weight
      const double w = 0.5;

      //histogram filling

      if ((ee_event || mm_event) && goodjets.size() > 0) {

        FourMomentum j1(goodjets[0].momentum());

        _h_first_jet_pt->fill(j1.pt(),w);
        _h_first_jet_abseta->fill(fabs(j1.eta()),w);
        if ( ee_event ) _h_Z_pt->fill(zees[0].pt(),w);
        if ( mm_event ) _h_Z_pt->fill(zmumus[0].pt(),w);
        _h_HT->fill(Ht,w);
        if ( ee_event ) _h_Dphi_Zj->fill(deltaPhi(zees[0], j1),w);
        if ( mm_event ) _h_Dphi_Zj->fill(deltaPhi(zmumus[0], j1),w);

        if ( jb_final.size() > 0 ) {

          FourMomentum b1(jb_final[0].momentum());

          _h_bjet_multiplicity->fill(1.,w);

          _h_first_bjet_pt_b->fill(b1.pt(),w);
          _h_first_bjet_abseta_b->fill(fabs(b1.eta()),w);
          if ( ee_event ) _h_Z_pt_b->fill(zees[0].pt(),w);
          if ( mm_event ) _h_Z_pt_b->fill(zmumus[0].pt(),w);
          _h_HT_b->fill(Ht,w);
          if ( ee_event ) _h_Dphi_Zb_b->fill(deltaPhi(zees[0], b1.phi()),w);
          if ( mm_event ) _h_Dphi_Zb_b->fill(deltaPhi(zmumus[0], b1.phi()),w);

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

    }


  private:

    /// @name Histograms

    Histo1DPtr     _h_first_jet_pt, _h_first_bjet_pt_b;
    Histo1DPtr     _h_first_jet_abseta, _h_first_bjet_abseta_b;
    Histo1DPtr     _h_Z_pt, _h_Z_pt_b;
    Histo1DPtr     _h_HT, _h_HT_b;
    Histo1DPtr     _h_Dphi_Zj, _h_Dphi_Zb_b;

    Scatter2DPtr     _h_first_jet_pt_ratio;
    Scatter2DPtr     _h_first_jet_abseta_ratio;
    Scatter2DPtr     _h_Z_pt_ratio;
    Scatter2DPtr     _h_HT_ratio;
    Scatter2DPtr     _h_Dphi_Zj_ratio;

    Histo1DPtr     _h_first_bjet_pt_bb, _h_second_bjet_pt_bb;
    Histo1DPtr     _h_Z_pt_bb;
    Histo1DPtr     _h_bb_mass_bb, _h_Zbb_mass_bb;
    Histo1DPtr     _h_Dphi_bb, _h_DR_bb, _h_DR_Zbmin_bb, _h_A_DR_Zb_bb;

    Histo1DPtr     _h_bjet_multiplicity;


    JetDefinition base_jet_def;
    JetDefinition flav_jet_def;
    FlavRecombiner flav_recombiner;

    bool debug = false;

    // GHS parameters:
    double GHS_alpha; // < flav-kt distance parameter alpha
    double GHS_beta; // < SoftDrop parameter beta for flavour clusters
    double GHS_zcut; // < SoftDrop zcut for flavour clusters
    double GHS_Rcut; // < Rcut parameter
    double GHS_omega; // < omega parameter for GHS_Omega (omega = 0 uses DeltaR_ij^2)
    double GHS_ptcut; // < overall ptcut

    SDFlavourCalc sdFlavCalc;

    int flavAlg;
    enum {
      IFN = 0,
      CMP = 1,
      GHS = 2,
      SDF = 3,
      AKT = 4,
      TAG = 5,
      OTAG = 6,
      CONE = 7,
    };

    std::string flavAlgName() {
      if(flavAlg == IFN) return "IFN";
      else if(flavAlg == CMP) return "CMP";
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
