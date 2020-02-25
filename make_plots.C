#include <fstream>
#include <sstream>

// STVs (for .root-file labeling)
std::string p_mu("mc_truth_mu_mom");
std::string phi_mu("mc_truth_mu_phi");
std::string theta_mu("mc_truth_mu_theta");
std::string costheta_mu("mc_truth_mu_costheta");

std::string p_p("mc_truth_leading_p_mom");
std::string phi_p("mc_truth_leading_p_phi");
std::string theta_p("mc_truth_leading_p_theta");
std::string costheta_p("mc_truth_leading_p_costheta");

std::string p_p2("mc_truth_2nd_leading_p_mom");
std::string phi_p2("mc_truth_2nd_leading_p_phi");
std::string theta_p2("mc_truth_2nd_leading_p_theta");
std::string costheta_p2("mc_truth_2nd_leading_p_costheta");

std::string p_pi("mc_truth_leading_pi_mom");
std::string phi_pi("mc_truth_leading_pi_phi");
std::string theta_pi("mc_truth_leading_pi_theta");
std::string costheta_pi("mc_truth_leading_pi_costheta");

std::string p_n("mc_truth_leading_n_mom");
std::string phi_n("mc_truth_leading_n_phi");
std::string theta_n("mc_truth_leading_n_theta");
std::string costheta_n("mc_truth_leading_n_costheta");

std::string p_T("delta_p_T");
std::string alpha_T("delta_alpha_T");
std::string phi_T("delta_phi_T");
std::string p_T_2p("delta_p_T_2p");
std::string p_T_nN("delta_p_T_nN");

std::string pL("delta_p_L");
std::string delta_p("delta_p");

std::string pL_noE("delta_p_L_noE");
std::string delta_p_noE("delta_p_noE");

std::string E_kin("E_kin");
std::string E_nu("E_nu");

std::string variable;
TLatex gvname;

constexpr int HIST_LINE_WIDTH = 5;
int HIST_LINE_STYLE;

constexpr double AXIS_LABEL_SIZE = 0.035;
constexpr double AXIS_TITLE_SIZE = 0.051;
constexpr double X_AXIS_TITLE_OFFSET = 1.15;
constexpr double Y_AXIS_TITLE_OFFSET_XSEC = 0.9;
constexpr double Y_AXIS_TITLE_OFFSET_EVENTS = 0.9;
constexpr double FLUX_INTEGRAL_RELATIVE_TOLERANCE = 1e-5;
const double PI = std::acos(-1.);

constexpr double MICROBOONE_FIDUCIAL_MASS = 60e6; // grams
constexpr double ATOMIC_MASS_NATURAL_ARGON = 6.6335209e-23; // grams
constexpr double NUM_TARGETS = MICROBOONE_FIDUCIAL_MASS / ATOMIC_MASS_NATURAL_ARGON;

constexpr int A_TARGET = 40; // 40Ar mass number

// Number of GiBUU runs used to make the samples
// This is set by the job card parameter "num_runs_SameEnergy"
constexpr int NUM_GIBUU_RUNS = 17;

// Global cut that will be applied to all plots. For the 2D
// sliced plots, the appropriate additions will be made to
// this basic cut.

int cat_index = 1;
std::string legend_attr;

// For CC events, include 'cc', for QE events, include 'cat == 1' and so on
//const std::string CUT_TO_USE( "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 1" ); // CCQE and threshold
//  " && nf_pi_above_threshold == 0" );
std::string CUT_TO_USE;
//const std::string CUT_TO_USE( "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 2" ); // CC MEC and threshold
//const std::string CUT_TO_USE( "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 3" ); // CC RES and threshold

// String that will be prepended to all generated plot titles
//const std::string PLOT_TITLE_PREFIX( "1#mu + Np" ); // turn on plot title later
const std::string PLOT_TITLE_PREFIX( "" );

constexpr int DEFAULT_NUM_BINS = 50;
constexpr int DEFAULT_NUM_CUT_BINS = 4;

// Alternative to std::to_string that allows us to control the precision.
// See https://stackoverflow.com/a/16606128/4081973
std::string double_to_string(double x, int precision = 2) {
  std::ostringstream out;
  out.precision( precision );
  out << std::fixed << x;
  return out.str();
}

// Struct that stores the branch name, plot title, etc. for a particular
// variable in the STV TTree
struct VarInfo {
  VarInfo(double axis_min, double axis_max, const std::string& hist_name,
    const std::string& title1, const std::string& title2,
    const std::string& display_name, int num_bins = DEFAULT_NUM_BINS,
    int num_bins_cut = DEFAULT_NUM_CUT_BINS) : fAxisMin( axis_min ), fAxisMax( axis_max ),
    fNumBins( num_bins ), fNumBinsCut( num_bins_cut ), fHistNameBase( hist_name ),
    fHistTitlePart1( title1 ), fHistTitlePart2( title2 ),
    fVarDisplayName( display_name ) {}

  double fAxisMin;
  double fAxisMax;
  int fNumBins;
  int fNumBinsCut;
  std::string fHistNameBase;
  std::string fHistTitlePart1;
  std::string fHistTitlePart2;
  std::string fVarDisplayName;
};

// Type used to instruct STVPlotMaker to create plots in terms of differential
// cross sections or in terms of expected event rates in MicroBooNE
enum class STVPlotMode {
  CrossSections,
  FiducialEvents,
  MCProbDensity
};

class STVPlotMaker {

  public:

    // The member variables we need to make plots are re-initialized
    // with every call to make_plots(), so don't bother to do anything
    // with them in the constructor
    STVPlotMaker( STVPlotMode mode = STVPlotMode::CrossSections )
      : fPlotMode( mode )
    {
      initialize_variable_info();
    }

    std::pair< std::string, std::vector<TH1*> > make_plots(
      const std::string& generator,
      const std::string& flux_filename, const std::string& spline_filename,
      const std::string& stv_tree_filename, const std::string& gst_tree_filename,
      int plot_color, std::string& plot_legend)
    {
      // By default, assume that we're using unweighted events
      fUsingWeightedEvents = false;

      // Read my_plots_input_*.txt
      TFile flux_file( flux_filename.c_str(), "read" ); // Flux file
      TFile spline_file( spline_filename.c_str(), "read" ); // Spline file (relevant for GENIE, else input 'dummy.root')
      TFile gst_tree_file( gst_tree_filename.c_str(), "read" ); // generator output .root-file
      TFile stv_tree_file( stv_tree_filename.c_str(), "read" ); // stv_calculation output .root-file
/*
      // For legend
      if ( cat_index == 1 ) {
      legend_attr = "Total CC <FSIs on>";
      }

      if ( cat_index == 2 ) {
      legend_attr = "CC QE <FSIs on>";
      }

      if ( cat_index == 3 ) {
          if ( generator == "GiBUU" ) {
          legend_attr = "CC 2p2h <FSIs on>";
          }
          else legend_attr = "CC MEC <FSIs on>";
      }

      if ( cat_index == 4 ) {
      legend_attr = "Other CC <FSIs on>";
      }

      if ( cat_index == 5 ) {
      legend_attr = "Total CC <FSIs off>";
      }
*/
	// GENIE
      if ( generator == "GENIE" ) {
//      TFile flux_file( flux_filename.c_str(), "read" ); // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu", flux_hist );

//      TFile spline_file( spline_filename.c_str(), "read" ); // Spline file (relevant for GENIE, else input 'dummy.root')
      TGraph* spline_graph = nullptr;
      spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

//      TFile gst_tree_file( gst_tree_filename.c_str(), "read" ); // generator output .root-file (need it only for cc cut)
      fGSTTree = nullptr;
      gst_tree_file.GetObject( "gst", fGSTTree );
      assert( fGSTTree );

//      TFile stv_tree_file( stv_tree_filename.c_str(), "read" ); // stv_calculation output .root-file
      fSTVTree = nullptr;
      stv_tree_file.GetObject( "stv_tree", fSTVTree );
      assert( fSTVTree );
      fSTVTree->AddFriend( fGSTTree );

      fTotalXSecAvg = flux_averaged_total_xsec(*flux_hist, *spline_graph);

      TNamed* GENIE_tune = nullptr;
      stv_tree_file.GetObject( "GENIE_tune", GENIE_tune );

//      fGenieTuneName = Form( "%s_%d", GENIE_tune->GetTitle(), cat_index );
      }

      // NuWro
      if ( generator == "NuWro" ) {

      // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu", flux_hist );

      // Dummy for spline file
//      TGraph* spline_graph = nullptr;
//      spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

      // This is the NuWro output file (I need it only to get the xsections (bool cc is in stv_tree))
      TH1D* xsec_hist = nullptr;
      gst_tree_file.GetObject( "xsections", xsec_hist );
//      assert( fGSTTree );

//      TFile stv_tree_file( stv_tree_filename.c_str(), "read" ); // stv_calculation output .root-file
      fSTVTree = nullptr;
      stv_tree_file.GetObject( "stv_tree", fSTVTree );
//      assert( fSTVTree );
//      fSTVTree->AddFriend( fGSTTree );

//      fTotalXSecAvg = 1.01987e-38;
      fTotalXSecAvg = flux_averaged_total_xsec_NuWro(*xsec_hist);

      TNamed* NuWro_tune = nullptr;
      stv_tree_file.GetObject( "NuWro_tune", NuWro_tune );

//      fGenieTuneName = Form( "%s_%s", NuWro_tune->GetTitle(), legend_attr.c_str() );
      }

      // GiBUU
      if ( generator == "GiBUU" ) {

      // GiBUU uses weighted MC events
      fUsingWeightedEvents = true;

      // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu", flux_hist );

      // Dummy for spline file
//      TGraph* spline_graph = nullptr;
//      spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

      // This is the GiBUU output file (I need the GiBUU output file only to get the xsections (bool cc is in stv_tree and in GiBUU all events are CC events as it is specified in jobCard already))
      fGSTTree = nullptr;
      gst_tree_file.GetObject( "RootTuple", fGSTTree );
      fGSTTree->Draw("weight>>weight_hist", "", "goff");
      TH1D *weight_hist = (TH1D*)gDirectory->Get("weight_hist");
//	weight_hist->Draw();

//      assert( fGSTTree );

//      TFile stv_tree_file( stv_tree_filename.c_str(), "read" ); // stv_calculation output .root-file
      fSTVTree = nullptr;
      stv_tree_file.GetObject( "stv_tree", fSTVTree );
//      assert( fSTVTree );
//      fSTVTree->AddFriend( fGSTTree );

      fTotalXSecAvg = 1.; // GiBUU uses weighted events, so we don't need this factor

      TNamed* GiBUU_tune = nullptr;
      stv_tree_file.GetObject( "GiBUU_tune", GiBUU_tune );

//      fGenieTuneName = GiBUU_tune->GetTitle();
//      fGenieTuneName = Form( "%s", legend_attr.c_str() );
//	if (cat_index == 1 ) fGenieTuneName = Form( "%s", "Total CC <FSIs on>" );
//        if (cat_index == 2) fGenieTuneName = "Total CC <FSI off>";
/*
      if ( cat_index == 3 ) {
          if ( generator == "GiBUU" ) {
          legend_attr = "CC 2p2h (FSIs on)";
          }
          else legend_attr = "CC MEC (FSIs on)";
      }

      if ( cat_index == 4 ) {
      legend_attr = "Other CC (FSIs on)";
      }

      if ( cat_index == 5 ) {
      legend_attr = "Total CC (FSIs off)";
      }
*/
//      fGenieTuneName = Form( "%s", "dddddddd d" );
      }

      // NEUT
      if ( generator == "NEUT" ) {

      // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu", flux_hist );

      // Dummy for spline file
//      TGraph* spline_graph = nullptr;
//      spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

      // This is the NEUT output file (I need the NEUT output file only to get the xsections
      TH1D* evtrt_hist = nullptr;
      gst_tree_file.GetObject( "evtrt_numu", evtrt_hist );
//      fGSTTree->Draw("weight>>weight_hist", "", "goff");
//      TH1D *weight_hist = (TH1D*)gDirectory->Get("weight_hist");
//      weight_hist->Draw();

//      assert( fGSTTree );

//      TFile stv_tree_file( stv_tree_filename.c_str(), "read" ); // stv_calculation output .root-file
      fSTVTree = nullptr;
      stv_tree_file.GetObject( "stv_tree", fSTVTree );
//      assert( fSTVTree );
//      fSTVTree->AddFriend( fGSTTree );

//      fTotalXSecAvg = 1.01987e-38;
      fTotalXSecAvg = flux_averaged_total_xsec_NEUT( *flux_hist, *evtrt_hist);

      TNamed* NEUT_tune = nullptr;
      stv_tree_file.GetObject( "NEUT_tune", NEUT_tune );

//      fGenieTuneName = NEUT_tune->GetTitle();
//      fGenieTuneName = Form( "%s", legend_attr.c_str() );
      }


      // The gspl2root files use units of 1e-38 cm^2 for the cross section
      // splines, so make this adjustment before continuing if we're
      // making plots in event mode
      if ( fPlotMode == STVPlotMode::FiducialEvents ) {
        // fTotalXSecAvg is now in cm^2
        fTotalXSecAvg *= 1e-38;
      }


      if ( cat_index == 1 ) {
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1"; // all channels, no cuts, Npi
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 1"; // all channels, no cuts, 1pi
      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 0"; // all channels, no cuts, 0pi
      HIST_LINE_STYLE = 1; // solid line
      fGenieTuneName = Form( "%s", "Total CC <FSIs on>" );
      }

      if ( cat_index == 2 ) {
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1"; // all channels, no cuts, Npi
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 1"; // all channels, no cuts, 1pi
      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 0"; // all channels, no cuts, 0pi
      HIST_LINE_STYLE = 2; // dotted line
      fGenieTuneName = "Total CC <FSIs off>";
      }

      if ( cat_index == 3 ) {
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 1"; // QE, Npi
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 1"; // all channels, no cuts, 1pi
      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 1 && nf_pi_above_threshold == 0"; // QE, 0pi
      HIST_LINE_STYLE = 9; // dashed line
      fGenieTuneName = "CCQE <FSI on>";
      }

      if ( cat_index == 4 ) {
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 2"; // MEC, Npi
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 1"; // all channels, no cuts, 1pi
      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat == 2 && nf_pi_above_threshold == 0"; // MEC, 0pi
      HIST_LINE_STYLE = 10; // dashed dotted line
          if ( generator == "GiBUU" || generator == "NEUT" ) {
          fGenieTuneName = "CC2p2h <FSIs on>";
          }
          else fGenieTuneName = "CCMEC <FSIs on>";
      }

      if ( cat_index == 5 ) {
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat > 2"; // Other channels, Npi
//      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && nf_pi_above_threshold == 1"; // all channels, no cuts, 1pi
      CUT_TO_USE = "cc && STVs_ok && nf_p_above_threshold >= 1 && cat > 2 && nf_pi_above_threshold == 0"; // Other channels, 0pi
      HIST_LINE_STYLE = 8; // different dotted line
      fGenieTuneName = "Other CC <FSIs on>";
      }

      if ( cat_index == 2 ) {
      legend_attr = "CC QE <FSIs on>";
      }

      if ( generator == "GiBUU" ) {
        // If we're working with GiBUU events (which are weighted), apply the appropriate
        // weights to the output histograms by including a factor in the cut
        CUT_TO_USE = "weight * (" + CUT_TO_USE + ')';
      }

      cout << "Generator is " << generator << endl;
//      cout << "Tune name is " << fGenieTuneName << endl;
      cout << "CUT_TO_USE is " << CUT_TO_USE << endl;

      cout << "Total xsec_averaged is " << fTotalXSecAvg << endl;

      fNumEvents = fSTVTree->GetEntries();
/*
      TNamed* GENIE_tune = nullptr;
      stv_tree_file.GetObject( "GENIE_tune", genie_tune );

      fGenieTuneName = GENIE_tune->GetTitle();
*/
      fPlotColor = plot_color;
      fPlotLegend = plot_legend;

      std::vector<TH1*> histograms = this->build_histograms();

      // If we're working with GiBUU events, apply an extra scaling
      // factor needed to correct the usual approach
      if ( generator == "GiBUU" ) {
        for ( auto* hist : histograms ) {
          hist->Scale( A_TARGET * fNumEvents
            / static_cast<double>( NUM_GIBUU_RUNS ) );
        }
      }

      cout << "cat_index is " << cat_index << endl;
      cat_index += 1;

      return std::pair< std::string, std::vector<TH1*> >(
        fGenieTuneName, histograms );
    }

  protected:

    STVPlotMode fPlotMode;
    TTree* fSTVTree;
    TTree* fGSTTree;
    double fTotalXSecAvg;
    double fFluxIntegral;
    int fNumEvents;
    std::string fGenieTuneName;
    int fPlotColor;
    std::string fPlotLegend;
    bool fUsingWeightedEvents;

    std::map<std::string, VarInfo> fVariableInfo;

    // Sets up the plot variable info map for the requested plot mode
    void initialize_variable_info();

    // Computes the flux-averaged total cross section given a histogram
    // of the NEUTrino flux and a total cross section TGraph (both with
    // arbitrary units)
    // For GENIE
    double flux_averaged_total_xsec(TH1& flux_hist,
      const TGraph& xsec_spline)
    {
      fFluxIntegral = flux_hist.Integral( "width" );

      // Integrate up to the smaller of either the final point on the
      // total cross section spline or the left edge of the overflow bin
      // of the flux histogram
      // TODO: Add more careful checking here. This assumes that you're doing
      // something reasonable.
      double spline_max_energy = xsec_spline.GetX()[ xsec_spline.GetN() - 1 ];
      double flux_max_energy = flux_hist.GetBinLowEdge( flux_hist.GetNbinsX() + 1 );
      double max_energy = std::min( spline_max_energy, flux_max_energy );

      // Function to use for numerical integration. Element x[0] is the
      // NEUTrino energy
      TF1 xsec_weighted_flux_func("temp_func", [&](double* x, double*)
        {
          int flux_bin = flux_hist.FindBin( x[0] );
          double flux = flux_hist.GetBinContent( flux_bin );

          double total_xsec = xsec_spline.Eval( x[0] );

          return flux * total_xsec / fFluxIntegral;
        }, 0., max_energy, 0);

      return xsec_weighted_flux_func.Integral(0., max_energy,
        FLUX_INTEGRAL_RELATIVE_TOLERANCE);
    }

    // For NuWro
    double flux_averaged_total_xsec_NuWro( TH1D& xsec_hist )
    {
      return xsec_hist.Integral() * 1e+38 * 40; // cross section in NuWro is given per nucleon (we have 40Ar)
    }

    // For GiBUU
    double flux_averaged_total_xsec_GiBUU( TH1D& weight_hist )
    {
      // The cross-section is the mean of the weight times the integral
      return weight_hist.GetMean() * weight_hist.Integral();
    }

    // For NEUT
    double flux_averaged_total_xsec_NEUT( TH1D& flux_hist, TH1D& evtrt_hist )
    {

      double flux_numu = flux_hist.GetSumOfWeights();
      double evtrt_numu = evtrt_hist.GetSumOfWeights();

//      if ( cat_index == 2 ) return 40 * evtrt_numu /flux_numu/2; // FSI off: consider factor 2 due to hadd
     /* else*/ return 40 * evtrt_numu /flux_numu;
    }

    void xsec_normalize_event_hist(TH1& hist)
    {
      bool has_y_axis = hist.GetDimension() > 1;
      bool has_z_axis = hist.GetDimension() > 2;

      int overflow_x = hist.GetNbinsX() + 1;
      int overflow_y = has_y_axis ? hist.GetNbinsY() + 1 : 0;
      int overflow_z = has_z_axis ? hist.GetNbinsZ() + 1 : 0;

      // Get global index of the last overflow bin
      int last_overflow_bin_index = hist.GetBin(overflow_x, overflow_y, overflow_z);

      for (int bin = 0; bin <= last_overflow_bin_index; ++bin) {

        // Get axis-specific bin indices for the current bin (with global
        // bin index "bin")
        int bin_x_index, bin_y_index, bin_z_index;
        hist.GetBinXYZ( bin, bin_x_index, bin_y_index, bin_z_index );

        // Get the bin's width in the x, y, and z dimensions. If the histogram
        // is less than 3D, set the widths to unity (so that dividing by them
        // below will have no effect)
        double x_width = hist.GetXaxis()->GetBinWidth( bin_x_index );
        double y_width = hist.GetDimension() > 1 ? hist.GetYaxis()->GetBinWidth( bin_y_index ) : 1.;
        double z_width = hist.GetDimension() > 2 ? hist.GetZaxis()->GetBinWidth( bin_z_index ) : 1.;

        // Renormalize the current bin to represent the MC estimate of the
        // appropriate differential cross section. Also set the bin error to the
        // appropriate value from the binomial distribution.
        double old_bin_count = hist.GetBinContent( bin );

        // MC statistical uncertainty on the old bin count
        double error_old_bin_count = 0.;
        if ( !fUsingWeightedEvents ) {
          // Note: unweighted event bin counts follow a binomial distribution
          error_old_bin_count = std::sqrt( (fNumEvents - old_bin_count)
            * old_bin_count / fNumEvents );
        }
        else {
          // We approximate the error on weighted event bin counts using a
          // Poisson distribution. When histograms are filled with weighted
          // events (and TH1::Sumw2() has been called), ROOT sets the bin error
          // like this automatically. See
          // https://www.pp.rhul.ac.uk/~cowan/stat/notes/weights.pdf
          // TODO: revisit this, come up with a better estimator
          error_old_bin_count = hist.GetBinError( bin );
        }

        double new_bin_count, new_bin_error;

        if ( fPlotMode == STVPlotMode::CrossSections ) {
          new_bin_count = old_bin_count * fTotalXSecAvg
            / ( fNumEvents * x_width * y_width * z_width );
          new_bin_error = fTotalXSecAvg * error_old_bin_count / ( fNumEvents * x_width * y_width * z_width );
        }
        else if ( fPlotMode == STVPlotMode::FiducialEvents ) {
          new_bin_count = old_bin_count * fTotalXSecAvg * fFluxIntegral * NUM_TARGETS / fNumEvents;

          new_bin_error = error_old_bin_count * fTotalXSecAvg * fFluxIntegral * NUM_TARGETS / fNumEvents;
        }
        else {
          // fPlotMode == STVPlotMode::MCProbDensity
          new_bin_count = old_bin_count / ( fNumEvents * x_width * y_width * z_width );

          new_bin_error = error_old_bin_count / ( fNumEvents * x_width * y_width * z_width );
        }

        hist.SetBinContent(bin, new_bin_count);
        hist.SetBinError(bin, new_bin_error);
      }
    }

    void fill_histogram(TH1& hist, const std::string& var_exp,
      const std::string& cuts = "STVs_ok", const std::string& options = "")
    {
      fSTVTree->Project( hist.GetName(), var_exp.c_str(), cuts.c_str(), options.c_str() );

      xsec_normalize_event_hist( hist );

      hist.SetDirectory( nullptr );
      hist.SetStats( false );
      hist.SetLineColor( fPlotColor );
      hist.SetLineWidth( HIST_LINE_WIDTH );
      hist.SetLineStyle( HIST_LINE_STYLE );

      hist.GetXaxis()->SetLabelSize( AXIS_LABEL_SIZE );
      hist.GetXaxis()->SetTitleSize( AXIS_TITLE_SIZE );
      hist.GetXaxis()->SetTitleOffset( X_AXIS_TITLE_OFFSET );

      hist.GetYaxis()->SetLabelSize( AXIS_LABEL_SIZE );
      hist.GetYaxis()->SetTitleSize( AXIS_TITLE_SIZE );

      double y_offset = (fPlotMode == STVPlotMode::CrossSections) ?
        Y_AXIS_TITLE_OFFSET_XSEC : Y_AXIS_TITLE_OFFSET_EVENTS;
      hist.GetYaxis()->SetTitleOffset( y_offset );

    }

    std::vector<TH1*> build_histograms()
    {
      std::vector<TH1*> hists;

      // Make single-variable histograms
      for (const auto& pair : fVariableInfo) {
        const std::string& var_name = pair.first;

        // Create the strings needed to generate a unique histogram
        // name and create the plot and axis titles
        const auto& info = pair.second;
        std::string hist_name = info.fHistNameBase + fGenieTuneName;
        std::string hist_title = PLOT_TITLE_PREFIX + ' ' + info.fHistTitlePart1
          + ';' + info.fHistTitlePart2;

        // Create the histogram
        TH1D* temp_hist = new TH1D(hist_name.c_str(), hist_title.c_str(), // with histo_title
          info.fNumBins, info.fAxisMin, info.fAxisMax);

        // If we're using weighted events, then configure the histogram to
        // compute bin errors by summing the squares of the weights.
        if ( fUsingWeightedEvents ) temp_hist->Sumw2( true );

        this->fill_histogram( *temp_hist, var_name.c_str(), CUT_TO_USE.c_str() );
        hists.push_back( temp_hist );
/*
        // Create 2D plots cutting on the current variable
        for (const auto& pair2 : fVariableInfo) {
          const std::string& plot_var_name = pair2.first;
          const VarInfo& plot_var_info = pair2.second;
          if ( plot_var_name == var_name ) continue;

          double cut_step = (info.fAxisMax - info.fAxisMin) / info.fNumBinsCut;
          for ( int cb = 0; cb < info.fNumBinsCut; ++cb ) {
            double cut_min = cb*cut_step;
            double cut_max = cut_min + cut_step;

            std::string hist_cut_title = PLOT_TITLE_PREFIX + ' '
              + plot_var_info.fHistTitlePart1 + " for " + info.fVarDisplayName
              + " #in  [" + double_to_string(cut_min) + ", " + double_to_string(cut_max) + "];"
              + plot_var_info.fHistTitlePart2;

            TH1D* temp_cut_hist = new TH1D( (hist_name + "_cut" + std::to_string(cb)).c_str(),
              hist_cut_title.c_str(), plot_var_info.fNumBins, plot_var_info.fAxisMin,
              plot_var_info.fAxisMax );

            this->fill_histogram( *temp_cut_hist, plot_var_name,
              (CUT_TO_USE + " && " + std::to_string(cut_min) + " <= " + var_name + " && "
              + std::to_string(cut_max) + " >= " + var_name).c_str() );
//	      hists.push_back( temp_cut_hist ); // if commented, no 2D plots will be produced

          }
        }*/
      }

      return hists;
    }

};


// Specify settings for plotting each of the variables of interest
void STVPlotMaker::initialize_variable_info() {
  if ( fPlotMode == STVPlotMode::CrossSections ) {
    fVariableInfo = {



      { "mc_truth_mu_mom", VarInfo(0., 2.0, "mumom", "",
        "p^{truth}_{#mu} (GeV); d#sigma/p^{truth}_{#mu} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "p^{truth}_{#mu}") },

      { "mc_truth_mu_phi", VarInfo(-PI, PI, "muphi", "",
        "#phi^{truth}_{#mu} (GeV); d#sigma/#phi^{truth}_{#mu} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#phi^{truth}_{#mu}") },

      { "mc_truth_mu_theta", VarInfo(0, PI, "mutheta", "",
        "#theta^{truth}_{#mu} (GeV); d#sigma/#theta^{truth}_{#mu} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#theta^{truth}_{#mu}") },

      { "mc_truth_mu_costheta", VarInfo(-1., 1., "cosmutheta", "",
        "cos(#theta)^{truth}_{#mu} (GeV); d#sigma/cos(#theta)^{truth}_{#mu} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "cos(#theta)^{truth}_{#mu}") },

      { "mc_truth_leading_p_mom", VarInfo(0., 1.9, "pmom", "",
        "p^{truth}_{p_{lead}} (GeV); d#sigma/p^{truth}_{p_{lead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "p^{truth}_{p_{lead}}") },

      { "mc_truth_leading_p_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{p_{lead}} (GeV); d#sigma/#phi^{truth}_{p_{lead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#phi^{truth}_{p_{lead}}") },

      { "mc_truth_leading_p_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{p_{lead}} (GeV); d#sigma/#theta^{truth}_{p_{lead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#theta^{truth}_{p_{lead}}") },

      { "mc_truth_leading_p_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{p_{lead}} (GeV); d#sigma/cos(#theta)^{truth}_{p_{lead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "cos(#theta)^{truth}_{p_{lead}}") },

      { "mc_truth_2nd_leading_p_mom", VarInfo(0., 1.2, "pmom", "",
        "p^{truth}_{p_{2ndlead}} (GeV); d#sigma/p^{truth}_{p_{2ndlead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "p^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{p_{2ndlead}} (GeV); d#sigma/#phi^{truth}_{p_{2ndlead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#phi^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{p_{2ndlead}} (GeV); d#sigma/#theta^{truth}_{p_{2ndlead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#theta^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{p_{2ndlead}} (GeV); d#sigma/cos(#theta)^{truth}_{p_{2ndlead}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "cos(#theta)^{truth}_{p_{2ndlead}}") },

      { "mc_truth_leading_pi_mom", VarInfo(0., 1.8, "pmom", "",
        "p^{truth}_{#pi} (GeV); d#sigma/p^{truth}_{#pi} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "p^{truth}_{#pi}") },

      { "mc_truth_leading_pi_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{#pi} (GeV); d#sigma/#phi^{truth}_{#pi} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#phi^{truth}_{#pi}") },

      { "mc_truth_leading_pi_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{#pi} (GeV); d#sigma/#theta^{truth}_{#pi} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#theta^{truth}_{#pi}") },

      { "mc_truth_leading_pi_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{#pi} (GeV); d#sigma/cos(#theta)^{truth}_{#pi} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "cos(#theta)^{truth}_{#pi}") },

      { "mc_truth_leading_n_mom", VarInfo(0., 1.8, "pmom", "",
        "p^{truth}_{n} (GeV); d#sigma/p^{truth}_{n} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "p^{truth}_{n}") },

      { "mc_truth_leading_n_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{n} (GeV); d#sigma/#phi^{truth}_{n} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#phi^{truth}_{n}") },

      { "mc_truth_leading_n_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{n} (GeV); d#sigma/#theta^{truth}_{n} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#theta^{truth}_{n}") },

      { "mc_truth_leading_n_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{n} (GeV); d#sigma/cos(#theta)^{truth}_{p} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "cos(#theta)^{truth}_{n}") },

//      { "delta_p_T", VarInfo(0., 1.2, "pT", "d#sigma/#deltap_{T}",
      { "delta_p_T", VarInfo(0., 1.2, "pT", "",
        "#deltap_{T} (GeV); d#sigma/#deltap_{T} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#deltap_{T}") },

//      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "d#sigma/#delta#alpha_{T}",
      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "",
        "#delta#alpha_{T}; d#sigma/#delta#alpha_{T} (10^{-38} cm^{2} / ^{40}Ar)",
        "#delta#alpha_{T}") },

//      { "delta_phi_T", VarInfo(0., PI, "phiT", "d#sigma/#delta#phi_{T}",
      { "delta_phi_T", VarInfo(0., PI, "phiT", "",
        "#delta#phi_{T}; d#sigma/#delta#phi_{T} (10^{-38} cm^{2} / ^{40}Ar)",
        "#delta#phi_{T}") },

      { "delta_p_T_2p", VarInfo(0., 1.4, "pT_2p", "",
        "#deltap_{T_{2p}} (GeV); d#sigma/#deltap_{T_{2p}} (10^{-38} cm^{2} / ^{40}Ar)",
        "#deltap_{T_{2p}}") },

      { "delta_p_T_nN", VarInfo(0., 1.0, "pT_nN", "",
        "#deltap_{T_{nN}} (GeV); d#sigma/#deltap_{T_{nN}} (10^{-38} cm^{2} / ^{40}Ar)",
        "#deltap_{T_{nN}}") },

      { "delta_p_L", VarInfo(0., 0.6, "pL", "",
        "#deltap_{L} (GeV); d#sigma/#deltap_{L} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#deltap_{L}") },

      { "delta_p", VarInfo(0., 1.3, "deltap", "",
        "#deltap (GeV); d#sigma/#deltap (10^{-38} cm^{2} / ^{40}Ar / GeV)", "#deltap") },

      { "delta_p_L_noE", VarInfo(0., 0.6, "pL_noE", "",
        "#deltap_{L_{noE}} (GeV); d#sigma/#deltap_{L_{noE}} (10^{-38} cm^{2} / ^{40}Ar / GeV)",
        "#deltap_{L_{noE}}") },

      { "delta_p_noE", VarInfo(0., 1.3, "deltap_noE", "",
        "#deltap_{noE} (GeV); d#sigma/#deltap_{L_{noE}} (10^{-38} cm^{2} / ^{40}Ar / GeV)", "#deltap_{noE}") },

      { "E_kin_leading_p", VarInfo(0., 2.1, "E_kin_leading_p", "",
        "E_{kin_{p}} (GeV); d#sigma/E_{kin_{p}} (10^{-38} cm^{2} / ^{40}Ar / GeV)", "E_{kin_p}") },

      { "E_nu", VarInfo(0., 3.5, "E_nu", "",
        "E_{#nu} (GeV); d#sigma/E_{#nu} (10^{-38} cm^{2} / ^{40}Ar / GeV)", "E_{#nu}") }

    };
  }
  else if ( fPlotMode == STVPlotMode::FiducialEvents ) {
    fVariableInfo = {

      { "delta_p_T", VarInfo(0., 1.5, "pT", "#deltap_{T} events",
        "#deltap_{T} (GeV); fiducial events / 10^{20} POT", "#deltap_{T}") },

      { "delta_phi_T", VarInfo(0., PI, "phiT", "#delta#phi_{T} events",
        "#delta#phi_{T}; fiducial events / 10^{20} POT", "#delta#phi_{T}") },

      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "#delta#alpha_{T} events",
        "#delta#alpha_{T}; fiducial events / 10^{20} POT", "#delta#alpha_{T}") }

    };

  }
  else {
    // fPlotMode == STVPlotMode::MCProbDensity
    fVariableInfo = {
      { "delta_p_T", VarInfo(0., 1.5, "pT", "#deltap_{T} distribution",
        "#deltap_{T} (GeV); probability density", "#deltap_{T}") },

      { "delta_phi_T", VarInfo(0., PI, "phiT", "#delta#phi_{T} distribution",
        "#delta#phi_{T}; probability density", "#delta#phi_{T}") },

      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "#delta#alpha_{T} distribution",
        "#delta#alpha_{T}; probability density", "#delta#alpha_{T}") }

    };


  }

}


// Reads in a list of input files and makes histograms for each set (which will
// correspond to a particular GENIE tune). This function then draws the
// histograms together in comparison plots
void make_plots( )
{

// Choose generator and corresponding files
  const std::string input_filename( "input_GENIE_fsi.txt" );
//  const std::string input_filename( "input_NuWro_fsi.txt" );
//  const std::string input_filename( "input_GiBUU_fsi.txt" );
//  const std::string input_filename( "input_NEUT_fsi.txt" );

  // Keys are GENIE tune names, values are vectors of histograms
  std::map< std::string, std::vector<TH1*> > plots_map;

//  STVPlotMaker stv_plots( STVPlotMode::MCProbDensity );
  STVPlotMaker stv_plots( STVPlotMode::CrossSections );

  // Generate the histograms for each set of files. Store them indexed
  // by GENIE tune name in the plots map
  std::ifstream input_file( input_filename );
  std::string generator_name, flux_filename, spline_filename, gst_tree_filename, stv_tree_filename;
  int plot_color = 1;
  std::string plot_legend = "Total";

	// Create file
        TFile *f1 = new TFile("output.root","UPDATE");

  while ( input_file >> generator_name >> flux_filename >> spline_filename >> gst_tree_filename
    >> stv_tree_filename )
  {

    plots_map.emplace(
      stv_plots.make_plots( generator_name, flux_filename, spline_filename, stv_tree_filename,
        gst_tree_filename, plot_color, plot_legend)
    );

//    ++plot_color;
    // Skip yellow (too hard to see)
//    if ( plot_color == 5 ) ++plot_color;
    // This shade of green is too hard to distinguish from an earlier one
//    if ( plot_color == 8 ) plot_color += 3;

    if ( cat_index == 1 ) {
	plot_color = 1;
	plot_legend = "total";
    }

    if ( cat_index == 2 ) plot_color = 15; // grey
    if ( cat_index == 3 ) plot_color = 2; // red
    if ( cat_index == 4 ) plot_color = 4; // blue
    if ( cat_index == 5 ) plot_color = 8; // green

  }

  // The same number of histograms will be generated for all tunes, so we can
  // use this to figure out how many plots to make
  size_t num_plots = plots_map.cbegin()->second.size();

  for ( size_t p = 0; p < num_plots; ++p ) { // Loop over STV

    double y_max = 0.;
    int max_bin = 1;
    double x_max = 0;
    int left = 0;

    TCanvas* canvas = new TCanvas;
    canvas->SetBottomMargin(0.15);
    canvas->SetLeftMargin(0.12);

    if ( p == 0 ) variable = p_mu;
    if ( p == 1 ) variable = phi_mu;
    if ( p == 2 ) variable = theta_mu;
    if ( p == 3 ) variable = costheta_mu;
    if ( p == 4 ) variable = p_p;
    if ( p == 5 ) variable = phi_p;
    if ( p == 6 ) variable = theta_p;
    if ( p == 7 ) variable = costheta_p;
    if ( p == 8 ) variable = p_p2;
    if ( p == 9 ) variable = phi_p2;
    if ( p == 10 ) variable = theta_p2;
    if ( p == 11 ) variable = costheta_p2;
    if ( p == 12 ) variable = p_pi;
    if ( p == 13 ) variable = phi_pi;
    if ( p == 14 ) variable = theta_pi;
    if ( p == 15 ) variable = costheta_pi;
    if ( p == 16 ) variable = p_n;
    if ( p == 17 ) variable = phi_n;
    if ( p == 18 ) variable = theta_n;
    if ( p == 19 ) variable = costheta_n;
    if ( p == 20 ) variable = p_T;
    if ( p == 21 ) variable = alpha_T;
    if ( p == 22 ) variable = phi_T;
    if ( p == 23 ) variable = p_T_2p;
    if ( p == 24 ) variable = p_T_nN;
    if ( p == 25 ) variable = pL;
    if ( p == 26 ) variable = delta_p;
    if ( p == 27 ) variable = pL_noE;
    if ( p == 28 ) variable = delta_p_noE;
    if ( p == 29 ) variable = E_kin;
    if ( p == 30 ) variable = E_nu;

    canvas->SetName( ( generator_name + "_" + variable + "_" + CUT_TO_USE ).c_str() );
//    canvas->SetGrid();

//    TLegend* legend = new TLegend(0.60, 0.56, 0.86, 0.79); // TLegend position 1 (GiBUU and NEUT)
//    TLegend* legend = new TLegend(0.60, 0.58, 0.88, 0.81); // TLegend position 1 (GENIE)
    TLegend* legend = new TLegend(0.40, 0.65, 0.66, 0.89); // TLegend position 2 (for inlet)
//    TLegend* legend = new TLegend(0.37, 0.65, 0.63, 0.89); // TLegend position 2 (for inlet NuWro)

    legend->SetBorderSize( 0 );

    int num_plots_plots = 0;

    for ( auto iter = plots_map.cbegin(); iter != plots_map.cend(); ++iter ) { // Loop over tune

    num_plots_plots += 1;

      TH1* hist = iter->second.at(p);
      std::string s = iter->first;
      std::replace( s.begin(), s.end(), '<', '(' );
      std::replace( s.begin(), s.end(), '>', ')' );
      legend->AddEntry( hist, s.c_str(), "l" );
//      legend->AddEntry( hist, iter->first.c_str(), "l" );
//      legend->AddEntry( hist, plot_legend.c_str(), "l" );
//      cout << "Number of plots in plot: " << num_plots_plots << endl;
/*
      if ( num_plots_plots == 1 ) {
      legend->AddEntry( hist, "Total CC (FSIs on)", "l" );
      }

      if ( num_plots_plots == 2 ) {
      legend->AddEntry( hist, "Total CC (FSIs off)", "l" );
      }
*/

///*
      if ( iter == plots_map.cbegin() ) hist->Draw("hist e");
      else hist->Draw("hist e same");
//*/

/*
      if ( iter == plots_map.cbegin() ) hist->Draw("hist c");
      else hist->Draw("hist c same");
*/


      int bin = hist->GetMaximumBin();
      double y = hist->GetBinContent( bin );
      x_max = hist->GetXaxis()->GetXmax();
      if ( y > y_max ) {
        y_max = y;
        max_bin = bin;
      }
    }


    // Set the maximum y value to include the full range
    TH1* first_hist = plots_map.cbegin()->second.at(p);
    first_hist->GetYaxis()->SetRangeUser(0., /*y_max*1.02 */ y_max*1.12 );
    cout << "y_max " << y_max << endl;
    first_hist->Draw("hist e same");
//    first_hist->Draw("hist c same");

    // Move the legend to the left-hand side if needed
    if ( max_bin > first_hist->GetNbinsX() / 2 ) {
      legend->SetX1( 0.15 );
      legend->SetX2( 0.40 );

      left = 1;

      // If the distribution is nearly flat and close to the
      // maximum, move the legend to near the bottom of the plot
      if ( first_hist->GetBinContent(1) > 0.7*y_max ) {
        legend->SetY1( 0.15 );
        legend->SetY2( 0.40 );
      }
    }

    legend->Draw("same");

    gvname.SetTextSize(0.055);
    if ( generator_name == "GiBUU" && left == 0 ) {
//    gvname.DrawLatex( x_max - 0.365*x_max , y_max - 0.09*y_max, ( generator_name + " 2019" ).c_str() ); // Position 1
    gvname.DrawLatex( x_max - 0.25*x_max , y_max - 0.06*y_max, ( generator_name + " 2019" ).c_str() ); // Position 2 (for inlet)
    }
    else if ( generator_name == "GiBUU" ) {
    gvname.DrawLatex( 0. + 0.07*x_max , y_max - 0.09*y_max, ( generator_name + " 2019" ).c_str() );
    }

    if ( generator_name == "NuWro" && left == 0 ) {
    gvname.DrawLatex( x_max - 0.365*x_max , y_max - 0.09*y_max, ( generator_name + " 19.02.1" ).c_str() ); // Position 1
//    gvname.DrawLatex( x_max - 0.32*x_max , y_max - 0.06*y_max, ( generator_name + "  19.02.1" ).c_str() ); // Position 2 (for inlet)
    }
    else if ( generator_name == "NuWro" ) {
    gvname.DrawLatex( 0. + 0.07*x_max , y_max - 0.09*y_max, ( generator_name + " 19.02.1" ).c_str() );
    }

    if ( generator_name == "GENIE" && left == 0 ) {
//    gvname.DrawLatex( x_max - 0.365*x_max , y_max - 0.09*y_max, ( generator_name + " 3.0.6" ).c_str() ); // Position 1
    gvname.DrawLatex( x_max - 0.27*x_max , y_max - 0.03*y_max, ( generator_name + " 3.0.6" ).c_str() ); // Position 2 (for inlet)
    }
    else if ( generator_name == "GENIE" ) {
//    gvname.DrawLatex( 0. + 0.07*x_max , y_max - 0.03*y_max, ( generator_name + " 3.0.6" ).c_str() );
    gvname.DrawLatex( 0. + 0.37*x_max , y_max - 0.03*y_max, ( generator_name + " 3.0.6" ).c_str() );
    }

    if ( generator_name == "NEUT" && left == 0 ) {
    gvname.DrawLatex( x_max - 0.365*x_max , y_max - 0.09*y_max, ( generator_name + " 5.4.0" ).c_str() ); // Position 1
//    gvname.DrawLatex( x_max - 0.27*x_max , y_max - 0.06*y_max, ( generator_name + " 5.4.0" ).c_str() ); // Position 2 (for inlet)
    }
    else if ( generator_name == "NEUT" ) {
    gvname.DrawLatex( 0. + 0.07*x_max , y_max - 0.06*y_max, ( generator_name + " 5.4.0" ).c_str() );
    }

//    canvas->SaveAs( ("Plots/plot_" + generator_name + "_" + variable + "_" + CUT_TO_USE + ".pdf").c_str() );
    canvas->SaveAs( ("Plots/" + generator_name + "/plot" + std::to_string(p) + ".jpg").c_str() );

    f1->cd();
    canvas->Write();

    delete canvas;
    delete legend;
  }
}
