#include <fstream>
#include <sstream>

// Set cuts
std::string CUTS = "cc";

// Vars (for .root-file labeling)
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

std::string p_ni("mc_truth_initial_n_mom");

std::string decay_ang_mec("decay_ang_mec");

std::string variable;
TLatex gvname;

constexpr int HIST_LINE_WIDTH = 3;
int HIST_LINE_STYLE;

constexpr double AXIS_LABEL_SIZE = 0.035;
constexpr double AXIS_TITLE_SIZE = 0.051;
constexpr double X_AXIS_TITLE_OFFSET = 1.15;
constexpr double Y_AXIS_TITLE_OFFSET_XSEC = 1.2; 
constexpr double Y_AXIS_TITLE_OFFSET_EVENTS = 1.2;
constexpr double FLUX_INTEGRAL_RELATIVE_TOLERANCE = 1e-5;
const double PI = std::acos(-1.);

constexpr double FLUX_NORMALISATION_POT = 1.1e21;
constexpr double FLUX_NORMALISATION_CONVERSION_TO_CM2 = 1E-4;

//constexpr double MICROBOONE_FIDUCIAL_MASS = 60e6; // grams
//constexpr double DUNE_ND_FIDUCIAL_MASS = 67e6; // grams https://arxiv.org/pdf/2103.13910, p. 2-37
constexpr double DUNE_ND_FIDUCIAL_MASS = 6.02831e6;
constexpr double ATOMIC_MASS_NATURAL_ARGON = 6.6335209e-23; // grams
//constexpr double NUM_TARGETS = MICROBOONE_FIDUCIAL_MASS / ATOMIC_MASS_NATURAL_ARGON; // MicroBooNE
constexpr double NUM_TARGETS = DUNE_ND_FIDUCIAL_MASS / ATOMIC_MASS_NATURAL_ARGON; // DUNE

constexpr int A_TARGET = 40; // 40Ar mass number

// Number of GiBUU runs used to make the samples
// This is set by the job card parameter "num_runs_SameEnergy"
constexpr int NUM_GIBUU_RUNS = 50;

// Global cut that will be applied to all plots. For the 2D
// sliced plots, the appropriate additions will be made to
// this basic cut.

int cat_index = 1;
std::string legend_attr;

// For CC events, include 'cc', for QE events, include 'cat == 1' and so on
std::string CUT_TO_USE;

// String that will be prepended to all generated plot titles
//const std::string PLOT_TITLE_PREFIX( "1#mu + Np" ); // optional
const std::string PLOT_TITLE_PREFIX( "" );

constexpr int DEFAULT_NUM_BINS = 30;
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

// Type used to instruct VarPlotMaker to create plots in terms of differential
// cross sections or in terms of expected event rates in MicroBooNE
enum class VarPlotMode {
  CrossSections,
  FiducialEvents,
  MCProbDensity
};

class VarPlotMaker {

  public:

    // The member variables we need to make plots are re-initialized
    // with every call to make_plots(), so don't bother to do anything
    // with them in the constructor
    // VarPlotMaker( VarPlotMode mode = VarPlotMode::FiducialEvents )
    // VarPlotMaker( VarPlotMode mode = VarPlotMode::MCProbDensity )
    VarPlotMaker( VarPlotMode mode = VarPlotMode::CrossSections )
      : fPlotMode( mode )
    {
      initialize_variable_info();
    }

    std::pair< std::string, std::vector<TH1*> > make_plots(
      const std::string& generator,
      const std::string& flux_filename, const std::string& spline_filename,
      const std::string& var_tree_filename, const std::string& gst_tree_filename,
      int plot_color, std::string& plot_legend)
    {
      // By default, assume that we're using unweighted events
      fUsingWeightedEvents = false;

      cout << "A" << endl;

      // Read my_plots_input_*.txt
      TFile flux_file( flux_filename.c_str(), "read" ); // Flux file
      TFile spline_file( spline_filename.c_str(), "read" ); // Spline file (relevant for GENIE, else input 'dummy.root')
      TFile gst_tree_file( gst_tree_filename.c_str(), "read" ); // generator output .root-file
      TFile var_tree_file( var_tree_filename.c_str(), "read" ); // var_calculation output .root-file

      cout << "B" << endl;

    // GENIE
    if ( generator == "GENIE" ) {
      TH1D* flux_hist = nullptr; // For uBooNE flux, use TH1F, for DUNE flux, it is TH1D!! (otherwise gives Segmentation error)
      flux_file.GetObject( "numu_NDFHC_flux", flux_hist );

    if (!flux_hist) {
      std::cerr << "Error: flux_hist is null!" << std::endl;
    }

      TGraph* spline_graph = nullptr;
      spline_file.GetObject("nu_mu_Ar40/tot_cc", spline_graph); // CC inclusive

      fGSTTree = nullptr;
      gst_tree_file.GetObject( "gst", fGSTTree );
      assert( fGSTTree );

      fGRWGHTTree = nullptr;
      var_tree_file.GetObject( "NormCCMEC", fGRWGHTTree );
      //var_tree_file.GetObject( "DecayAngMEC", fGRWGHTTree );
      //var_tree_file.GetObject( "DecayAngMECLegendre", fGRWGHTTree );
      //var_tree_file.GetObject( "DecayAngMECLegendre2", fGRWGHTTree );
      //var_tree_file.GetObject( "twk_DecayAngMECLegendre", fGRWGHTTree );

      assert( fGRWGHTTree ); 

      fVarTree = nullptr;
      var_tree_file.GetObject( "var_tree", fVarTree );
      assert( fVarTree );
      fVarTree->AddFriend( fGSTTree );
      fVarTree->AddFriend( fGRWGHTTree );

      fTotalXSecAvg = flux_averaged_total_xsec(*flux_hist, *spline_graph);

      TNamed* GENIE_tune = nullptr;
      var_tree_file.GetObject( "GENIE_tune", GENIE_tune );
      }

      // NuWro
      if ( generator == "NuWro" ) {

      // Flux file
      TH1F* flux_hist = nullptr;
      flux_file.GetObject( "numu_NDFHC_flux", flux_hist );

      // Dummy for spline file

      // This is the NuWro output file (I need it only to get the xsections (bool cc is in var_tree))
      TH1D* xsec_hist = nullptr;
      gst_tree_file.GetObject( "xsections", xsec_hist );

      fVarTree = nullptr;
      var_tree_file.GetObject( "var_tree", fVarTree );

      // fTotalXSecAvg = 1.01987e-38;
      fTotalXSecAvg = flux_averaged_total_xsec_NuWro(*xsec_hist);

      TNamed* NuWro_tune = nullptr;
      var_tree_file.GetObject( "NuWro_tune", NuWro_tune );
      }

      // GiBUU
      if ( generator == "GiBUU" ) {

      // GiBUU uses weighted MC events
      fUsingWeightedEvents = true;

      // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu_NDFHC_flux", flux_hist );

      // Dummy for spline file

      // This is the GiBUU output file (I need the GiBUU output file only to get the xsections (bool cc is in var_tree and in GiBUU all events are CC events as it is specified in jobCard already))
      fGSTTree = nullptr;
      gst_tree_file.GetObject( "RootTuple", fGSTTree );
      fGSTTree->Draw("weight>>weight_hist", "", "goff");
      TH1D *weight_hist = (TH1D*)gDirectory->Get("weight_hist");

      fVarTree = nullptr;
      var_tree_file.GetObject( "var_tree", fVarTree );

      fTotalXSecAvg = 1.; // GiBUU uses weighted events, so we don't need this factor

      TNamed* GiBUU_tune = nullptr;
      var_tree_file.GetObject( "GiBUU_tune", GiBUU_tune );
      }

      // NEUT
      if ( generator == "NEUT" ) {

      // Flux file
      TH1D* flux_hist = nullptr;
      flux_file.GetObject( "numu_NDFHC_flux", flux_hist );

      // Dummy for spline file

      // This is the NEUT output file (I need the NEUT output file only to get the xsections
      TH1D* evtrt_hist = nullptr;
      gst_tree_file.GetObject( "evtrt_numu", evtrt_hist );

      fVarTree = nullptr;
      var_tree_file.GetObject( "var_tree", fVarTree );

      fTotalXSecAvg = flux_averaged_total_xsec_NEUT( *flux_hist, *evtrt_hist);

      TNamed* NEUT_tune = nullptr;
      var_tree_file.GetObject( "NEUT_tune", NEUT_tune );
      }

      // The gspl2root files use units of 1e-38 cm^2 for the cross section
      // splines, so make this adjustment before continuing if we're
      // making plots in event mode
      if ( fPlotMode == VarPlotMode::FiducialEvents ) {
        // fTotalXSecAvg is now in cm^2
        fTotalXSecAvg *= 1e-38;
      }

      std::string fGRWGHTTreeName = fGRWGHTTree->GetName();
      fGRWGHTTree->Show(1); // Show the first entry in the tree

      std::cout << "cat_index is " << cat_index << std::endl;

      if ( cat_index == 1 ) {
      CUT_TO_USE = "cc * weights->At(0)";
      std::cout << "fGRWGHTTreeName is " << fGRWGHTTreeName << std::endl;
      fGenieTuneName = Form( "%s", ( fGRWGHTTreeName + "  = 0 - default" ).c_str() );
      }
 
      if ( cat_index == 2 ) {
      CUT_TO_USE = "cc * weights->At(1)";
      HIST_LINE_STYLE = 2; // dotted line
      fGenieTuneName = fGRWGHTTreeName + " = 1";
      }

      if ( cat_index == 3 ) {
      CUT_TO_USE = "cc * weights->At(2)";
      HIST_LINE_STYLE = 9; // dashed line    canvas->SaveAs( ("Plots/" + generator_name + "/plot" + std::to_string(p) + ".png").c_str() );
      fGenieTuneName = fGRWGHTTreeName + " = 2";
      }

      if ( cat_index == 4 ) {
      CUT_TO_USE = "cc * weights->At(3)";
      HIST_LINE_STYLE = 10; // dashed dotted line
          if ( generator == "GiBUU" || generator == "NEUT" ) {
          fGenieTuneName = "";
          }
          else fGenieTuneName = fGRWGHTTreeName + " = 3";
      }

      if ( cat_index == 5 ) {
      CUT_TO_USE = "cc * weights->At(4)";
      HIST_LINE_STYLE = 8; // different dotted line
      fGenieTuneName = fGRWGHTTreeName + " = 4";
      }

      if ( cat_index == 6 ) {
      CUT_TO_USE = "cc * weights->At(5)";
      HIST_LINE_STYLE = 6; // different dotted line
      fGenieTuneName = fGRWGHTTreeName + " = 5";
      }

      if ( cat_index == 7 ) {
      CUT_TO_USE = "cc * weights->At(6)";
      HIST_LINE_STYLE = 4; // different dotted line
      fGenieTuneName = fGRWGHTTreeName + " = 6";
      }

      if ( cat_index == 8 ) {
      CUT_TO_USE = "cc * weights->At(7)";
      HIST_LINE_STYLE = 8; // different dotted line
      fGenieTuneName = fGRWGHTTreeName + " = 7";
      }

      if ( generator == "GiBUU" ) {
        // If we're working with GiBUU events (which are weighted), apply the appropriate
        // weights to the output histograms by including a factor in the cut
        CUT_TO_USE = "weight * (" + CUT_TO_USE + ')';
      }

      cout << "Generator is " << generator << endl;
      //cout << "Tune name is " << fGenieTuneName << endl;
      cout << "CUT_TO_USE is " << CUT_TO_USE << endl;

      cout << "Total xsec_averaged is " << fTotalXSecAvg << endl;

      fNumEvents = fVarTree->GetEntries();
      cout << "Number of events: " << fNumEvents << endl;

      fPlotColor = plot_color;
      fPlotLegend = plot_legend;

      std::vector<TH1*> histograms = this->build_histograms();

      // If we're working with GiBUU events, apply an extra scaling
      // factor needed to correct the usual approach
      if ( generator == "GiBUU" && fPlotMode == VarPlotMode::CrossSections ) {
        for ( auto* hist : histograms ) {
          hist->Scale( A_TARGET * fNumEvents
            / static_cast<double>( NUM_GIBUU_RUNS ) ); // per Ar nucleus
        }
      }

      cout << "new_bin_count = old_bin_count * fTotalXSecAvg * fFluxIntegral * NUM_TARGETS / fNumEvents;" << std::endl;
      cout << "fTotalXSecAvg = " << fTotalXSecAvg << std::endl;
      cout << "fScaleFactor = fTotalXSecAvg/fNumEvents = " << fTotalXSecAvg/fNumEvents << std::endl;
      cout << "fFluxIntegral = " << fFluxIntegral << std::endl;
      cout << "NUM_TARGETS = " << NUM_TARGETS << std::endl;
      cout << "fNumEvents = " << fNumEvents << std::endl;

      // xsec per nucleon (leave out if you want per nucleus) 
      for ( auto* hist : histograms ) {
        //hist->Scale( 1.
        //    / A_TARGET );
        // Following three lines in combination with FiducialEvents mode and new_bin_count = old_bin_count; below gives geniev3_validation_plots normalisation; It differs to FiducialEvents normalisation by additional factors that are num_POT and Fiducial_mass in tonnes; test four lines below here (use FiducialEvents mode normally and scale hist->Scale(5.98*1.1E21); where this is num_POT and Fiducial_mass in tonnes).
      }

      cout << "cat_index is " << cat_index << endl;
      cat_index += 1;

      return std::pair< std::string, std::vector<TH1*> >(
        fGenieTuneName, histograms );
    }

  protected:

    VarPlotMode fPlotMode;
    TTree* fVarTree;
    TTree* fGRWGHTTree;
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

      return 40 * evtrt_numu /flux_numu;
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

        if ( fPlotMode == VarPlotMode::CrossSections ) {
          new_bin_count = old_bin_count;
          new_bin_error = fTotalXSecAvg * error_old_bin_count / ( fNumEvents * x_width * y_width * z_width );
        }
        else if ( fPlotMode == VarPlotMode::FiducialEvents ) {
          new_bin_count = old_bin_count * fTotalXSecAvg * fFluxIntegral * NUM_TARGETS / fNumEvents * FLUX_NORMALISATION_POT * FLUX_NORMALISATION_CONVERSION_TO_CM2;
          new_bin_error = error_old_bin_count * fTotalXSecAvg * fFluxIntegral * NUM_TARGETS / fNumEvents;
        }
        else {
          // fPlotMode == VarPlotMode::MCProbDensity
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
      fVarTree->Project( hist.GetName(), var_exp.c_str(), cuts.c_str(), options.c_str() );

      xsec_normalize_event_hist( hist );

      hist.SetDirectory( nullptr );
      hist.SetStats( false );
      hist.SetLineColor( fPlotColor );
      hist.SetLineWidth( HIST_LINE_WIDTH );
      hist.SetLineStyle( HIST_LINE_STYLE );

      hist.GetXaxis()->SetLabelSize( AXIS_LABEL_SIZE );
      hist.GetXaxis()->SetTitleSize( AXIS_TITLE_SIZE );
      hist.GetXaxis()->SetTitleOffset( X_AXIS_TITLE_OFFSET );
      hist.GetXaxis()->CenterTitle(true);
      hist.GetYaxis()->CenterTitle(true);

      hist.GetYaxis()->SetLabelSize( AXIS_LABEL_SIZE );
      hist.GetYaxis()->SetTitleSize( AXIS_TITLE_SIZE );

      double y_offset = (fPlotMode == VarPlotMode::CrossSections) ?
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
        }
      }

      return hists;
    }

};


// Specify settings for plotting each of the variables of interest
void VarPlotMaker::initialize_variable_info() {
  if ( fPlotMode == VarPlotMode::CrossSections ) {
    fVariableInfo = {

      { "mc_truth_mu_mom", VarInfo(0., 2.0, "mumom", "",
        "p^{truth}_{#mu} [GeV/c]; #LT #frac{d#sigma}{dp^{truth}_{#mu}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{truth}_{#mu}") },
/*
      { "mc_truth_mu_phi", VarInfo(-PI, PI, "muphi", "",
        "#phi^{truth}_{#mu} [rad]; #LT #frac{d#sigma}{d#phi^{truth}_{#mu}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#phi^{truth}_{#mu}") },

      { "mc_truth_mu_theta", VarInfo(0, PI, "mutheta", "",
        "#theta^{truth}_{#mu} [rad]; #LT #frac{d#sigma}{d#theta^{truth}_{#mu}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#theta^{truth}_{#mu}") },    

      { "mc_truth_mu_costheta", VarInfo(-1., 1., "cosmutheta", "",
        "cos(#theta)^{truth}_{#mu}; #LT #frac{d#sigma}{dcos(#theta)^{truth}_{#mu}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta)^{truth}_{#mu}") },
*/
      { "mc_truth_leading_p_mom", VarInfo(0., 1.9, "pmom", "",
        "p^{truth}_{p_{lead}} [GeV/c]; #LT #frac{d#sigma}{dp^{truth}_{p_{lead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{truth}_{p_{lead}}") },
/*
      { "mc_truth_leading_p_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{p_{lead}} [rad]; #LT #frac{d#sigma}{d#phi^{truth}_{p_{lead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#phi^{truth}_{p_{lead}}") },

      { "mc_truth_leading_p_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{p_{lead}} [rad]; #LT #frac{d#sigma}{d#theta^{truth}_{p_{lead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#theta^{truth}_{p_{lead}}") },

      { "mc_truth_leading_p_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{p_{lead}}; #LT #frac{d#sigma}{dcos(#theta)^{truth}_{p_{lead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta)^{truth}_{p_{lead}}") },

      { "mc_truth_2nd_leading_p_mom", VarInfo(0., 1.2, "pmom", "",
        "p^{truth}_{p_{2ndlead}} [GeV/c]; #LT #frac{d#sigma}{dp^{truth}_{p_{2ndlead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{p_{2ndlead}} [rad]; #LT #frac{d#sigma}{d#phi^{truth}_{p_{2ndlead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#phi^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{p_{2ndlead}} [rad]; #LT #frac{d#sigma}{d#theta^{truth}_{p_{2ndlead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#theta^{truth}_{p_{2ndlead}}") },

      { "mc_truth_2nd_leading_p_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{p_{2ndlead}}; #LT #frac{d#sigma}{dcos(#theta)^{truth}_{p_{2ndlead}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta)^{truth}_{p_{2ndlead}}") },

      { "mc_truth_leading_pi_mom", VarInfo(0., 1.8, "pmom", "",
        "p^{truth}_{#pi} [GeV/c]; #LT #frac{d#sigma}{dp^{truth}_{#pi}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{truth}_{#pi}") },

      { "mc_truth_leading_pi_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{#pi} [rad]; #LT #frac{d#sigma}{d#phi^{truth}_{#pi}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#phi^{truth}_{#pi}") },

      { "mc_truth_leading_pi_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{#pi} [rad]; #LT #frac{d#sigma}{d#theta^{truth}_{#pi}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#theta^{truth}_{#pi}") },

      { "mc_truth_leading_pi_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{#pi}; #LT #frac{d#sigma}{dcos(#theta)^{truth}_{#pi}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta)^{truth}_{#pi}") },

      { "mc_truth_leading_n_mom", VarInfo(0., 1.8, "pmom", "",
        "p^{truth}_{n} [GeV/c]; #LT #frac{d#sigma}{dp^{truth}_{n}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{truth}_{n}") },

      { "mc_truth_leading_n_phi", VarInfo(-PI, PI, "pphi", "",
        "#phi^{truth}_{n} [rad]; #LT #frac{d#sigma}{d#phi^{truth}_{n}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#phi^{truth}_{n}") },

      { "mc_truth_leading_n_theta", VarInfo(0, PI, "ptheta", "",
        "#theta^{truth}_{n} [rad]; #LT #frac{d#sigma}{d#theta^{truth}_{n}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#theta^{truth}_{n}") },

      { "mc_truth_leading_n_costheta", VarInfo(-1., 1., "cosptheta", "",
        "cos(#theta)^{truth}_{n}; #LT #frac{d#sigma}{dcos(#theta)^{truth}_{p}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta)^{truth}_{n}") },

      { "delta_p_T", VarInfo(0., 2.0, "pT", "",
        "#deltap_{T} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{T}} #GT_{flux} #left[#frac{10^{-38} cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}GeV/c} #right]",
        "#deltap_{T}") },

      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "",
        "#delta#alpha_{T} [rad]; #LT #frac{d#sigma}{d#delta#alpha_{T}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#delta#alpha_{T}") },

      { "delta_phi_T", VarInfo(0., PI, "phiT", "",
        "#delta#phi_{T} [rad]; #LT #frac{d#sigma}{d#delta#phi_{T}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} rad} #right]",
        "#delta#phi_{T}") },

      { "delta_p_T_2p", VarInfo(0., 1.4, "pT_2p", "",
        "#deltap_{T_{2p}} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{T_{2p}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "#deltap_{T_{2p}}") },

      { "decay_ang_mec", VarInfo(-1., 1., "decayangmec", "",
        "cos(#theta_{DecayAngMEC}); #LT #frac{d#sigma}{cos(#theta_{DecayAngMEC})} #GT_{flux} #left[ #frac{10^{-39} cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}}} #right]",
        "cos(#theta_{DecayAngMEC})") },

      { "delta_p_T_nN", VarInfo(0., 1.0, "pT_nN", "",
        "#deltap_{T_{nN}} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{T_{nN}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "#deltap_{T_{nN}}") },

      { "delta_p_L", VarInfo(0., 0.6, "pL", "",
        "#deltap_{L} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{L}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "#deltap_{L}") },

      { "delta_p", VarInfo(0., 1.3, "deltap", "",
        "#deltap [GeV/c]; #LT #frac{d#sigma}{d#deltap} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "#deltap") },

      { "delta_p_L_noE", VarInfo(0., 0.6, "pL_noE", "",
        "#deltap_{L_{noE}} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{L_{noE}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "#deltap_{L_{noE}}") },

      { "delta_p_noE", VarInfo(0., 1.3, "deltap_noE", "",
        "#deltap_{noE} [GeV/c]; #LT #frac{d#sigma}{d#deltap_{L_{noE}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]", "#deltap_{noE}") },

      { "E_kin_leading_p", VarInfo(0., 2.1, "E_kin_leading_p", "",
        "E_{kin_{p}} [GeV]; #LT #frac{d#sigma}{dE_{kin_{p}}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV} #right]", "E_{kin_p}") },

      { "E_nu", VarInfo(0., 3.5, "E_nu", "",
        "E_{#nu} [GeV]; #LT #frac{d#sigma}{dE_{#nu}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV} #right]", "E_{#nu}") },

      { "mc_truth_initial_n_mom", VarInfo(0., 1.0, "pni", "Fermi Motion",
        "p^{initial}_{N} [GeV/c]; #LT #frac{d#sigma}{dp^{initial}_{N}} #GT_{flux} #left[ #frac{cm^{2}}{^{^{40}Ar}nucleon^{{}^{{}^{{}^{}}}} GeV/c} #right]",
        "p^{initial}_{N}") }
*/
    };
  }
  else if ( fPlotMode == VarPlotMode::FiducialEvents ) {
    fVariableInfo = {

      { "delta_p_T", VarInfo(0., 2.0, "pT", "#deltap_{T} events",
        "#deltap_{T} [GeV/c]; fiducial events / 1.1#times 10^{21} POT", "#deltap_{T}") },

      { "delta_phi_T", VarInfo(0., PI, "phiT", "#delta#phi_{T} events",
        "#delta#phi_{T}; fiducial events / 1.1#times 10^{21} POT", "#delta#phi_{T}") },

      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "#delta#alpha_{T} events",
        "#delta#alpha_{T}; fiducial events / 1.1#times 10^{21} POT", "#delta#alpha_{T}") },

      { "mc_truth_initial_n_mom", VarInfo(0., 1.0, "pni", "Fermi Motion",
        "p_{n}; fiducial events / 1.1#times 10^{21} POT", "p_{n}") },

      { "mc_truth_leading_p_mom", VarInfo(0., 2.0, "pmom", "",
        "p^{truth}_{p_{lead}} [GeV/c]; fiducial events / 1.1#times 10^{21} POT", "p^{truth}_{p_{lead}}") },

      { "decay_ang_mec", VarInfo(-1., 1., "decayangmec", "",
        "cos(#theta_{DecayAngMEC}); fiducial events / 1.1#times 10^{21} POT", "cos(#theta_{DecayAngMEC})") }

    };
  }
  else {
    // fPlotMode == VarPlotMode::MCProbDensity
    fVariableInfo = {
      { "delta_p_T", VarInfo(0., 1.5, "pT", "#deltap_{T} distribution",
        "#deltap_{T} [GeV/c]; probability density", "#deltap_{T}") },

      { "delta_phi_T", VarInfo(0., PI, "phiT", "#delta#phi_{T} distribution",
        "#delta#phi_{T}; probability density", "#delta#phi_{T}") },

      { "delta_alpha_T", VarInfo(0., PI, "alphaT", "#delta#alpha_{T} distribution",
        "#delta#alpha_{T}; probability density", "#delta#alpha_{T}") },

      { "mc_truth_initial_n_mom", VarInfo(0., 1.0, "pni", "Fermi Motion",
        "p_{n}; probability density", "p_{n}") }

    };
  }
}


// Reads in a list of input files and makes histograms for each set (which will
// correspond to a particular GENIE tune). This function then draws the
// histograms together in comparison plots
void make_plots( )
{

  // Choose generator and corresponding files
  const std::string input_filename( "input_files/input_GENIEReweight_NormCCMEC.txt" );

  // Keys are GENIE tune names, values are vectors of histograms
  std::map< std::string, std::vector<TH1*> > plots_map;

  //VarPlotMaker var_plots( VarPlotMode::FiducialEvents );
  //VarPlotMaker var_plots( VarPlotMode::MCProbDensity );
  VarPlotMaker var_plots( VarPlotMode::CrossSections );

  // Generate the histograms for each set of files. Store them indexed
  // by GENIE tune name in the plots map
  std::ifstream input_file( input_filename );
  std::string generator_name, flux_filename, spline_filename, gst_tree_filename, var_tree_filename;
  int plot_color = kOrange-3;
  std::string plot_legend = "Total";

  // Create file
  TFile *f1 = new TFile("output.root","UPDATE");

  while ( input_file >> generator_name >> flux_filename >> spline_filename >> gst_tree_filename
    >> var_tree_filename )
  {

    plots_map.emplace(
      var_plots.make_plots( generator_name, flux_filename, spline_filename, var_tree_filename,
        gst_tree_filename, plot_color, plot_legend)
    );

    if ( cat_index == 1 ) {
	plot_color = kOrange-3;
	plot_legend = "total";
    }

    if ( cat_index == 2 ) plot_color = kRed; 
    if ( cat_index == 3 ) plot_color = kViolet -4;
    if ( cat_index == 4 ) plot_color = 15;
    if ( cat_index == 5 ) plot_color = kGreen + 2;
    if ( cat_index == 6 ) plot_color = kOrange;
    if ( cat_index == 7 ) plot_color = kViolet + 2;
    if ( cat_index == 8 ) plot_color = kAzure -2;
  }


  // Make legend labels
  std::string l1 = "";
  std::string l2 = "";
  std::string l3 = ""; 
  std::string l4 = "";
  std::string l5 = "";
  std::string l6 = "";
  std::string l7 = "";
  std::string l8 = "";
  std::string l9 = "";
  std::string l10 = "";
  std::string l11 = "";
  std::string l12 = "";
  std::string l13 = "";
  std::string l15 = "(";
  std::string l16 = ")";

  // The same number of histograms will be generated for all tunes, so we can
  // use this to figure out how many plots to make
  size_t num_plots = plots_map.cbegin()->second.size();

  for ( size_t p = 0; p < num_plots; ++p ) { // Loop over Var

    double y_max = 0.;
    int max_bin = 1;
    double x_max = 0;
    int left = 0;

    TCanvas* canvas = new TCanvas;
    canvas->SetBottomMargin(0.15);
    canvas->SetRightMargin(0.05);
    //canvas->SetLeftMargin(0.11); // for one-line label
    canvas->SetLeftMargin(0.15); // for other label

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
    if ( p == 30 ) variable = p_ni;
    if ( p == 31 ) variable = decay_ang_mec;

    canvas->SetName( ( generator_name + "_" + variable + "_" + CUT_TO_USE ).c_str() );
    //canvas->SetGrid();

    TLegend* legend = new TLegend(0.61, 0.58, 0.94, 0.89); // TLegend position 2 (for inlet)
    legend->SetBorderSize( 0 );

    int num_plots_plots = 0;

    for ( auto iter = plots_map.cbegin(); iter != plots_map.cend(); ++iter ) { // Loop over tune

    num_plots_plots += 1;

      TH1* hist = iter->second.at(p);
      std::string s = iter->first;
      std::replace( s.begin(), s.end(), '<', '(' );
      std::replace( s.begin(), s.end(), '>', ')' );

      if ( s == "a" ) s.replace(s.begin(), s.end(), l1);
      if ( s == "b" ) s.replace(s.begin(), s.end(), l2);
      if ( s == "c" ) s.replace(s.begin(), s.end(), l3);
      if ( s == "d" ) s.replace(s.begin(), s.end(), l4);
      if ( s == "d2" ) s.replace(s.begin(), s.end(), l9);
      if ( s == "e" ) s.replace(s.begin(), s.end(), l5);
      if ( s == "f" ) s.replace(s.begin(), s.end(), l6);
      if ( s == "g" ) s.replace(s.begin(), s.end(), l7);
      if ( s == "h" ) s.replace(s.begin(), s.end(), l8);

      if ( s == "k" ) s.replace(s.begin(), s.end(), l15);
      if ( s == "l" ) s.replace(s.begin(), s.end(), l16);

      if ( s == "GENIE" ) s.replace(s.begin(), s.end(), l10);
      if ( s == "NuWro" ) s.replace(s.begin(), s.end(), l11);
      if ( s == "GiBUU" ) s.replace(s.begin(), s.end(), l12);
      if ( s == "NEUT" ) s.replace(s.begin(), s.end(), l13);

      std::replace( s.begin(), s.end(), '>', ')' );

      legend->AddEntry( hist, s.c_str(), "l" );

      if ( iter == plots_map.cbegin() ) hist->Draw("hist e");
      else hist->Draw("hist e same");

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
    //cout << "y_max " << y_max << endl;
    //cout << "0.7*y_max " << 0.7*y_max << endl;
    //cout << "Bin content max bin " << first_hist->GetBinContent(1) << endl;

    // Move the legend to the left-hand side if needed
    if ( max_bin > first_hist->GetNbinsX() / 2 ) {
      legend->SetX1( 0.15 ); // one-line label
      legend->SetX1( 0.17 ); // other label
      legend->SetX2( 0.54 ); 

      left = 1;
    }
      // If the distribution is nearly flat and close to the
      // maximum, move the legend to near the bottom of the plot
     // cout << first_hist->GetBinContent(1) << " && " << first_hist->GetBinContent(max_bin) << " && " << first_hist->GetBinContent( first_hist->GetNbinsX() - 1 ) << " 0.7 y_max " << 0.7*y_max << endl;

      if ( ( first_hist->GetBinContent( 1 ) > 0.7*y_max ) && ( first_hist->GetBinContent( first_hist->GetNbinsX() - 1 ) ) > 0.7*y_max  ) {
        legend->SetY1( 0.40 );
        legend->SetY2( 0.70 );
      }

    TLatex* left_title; // = new TLatex(0.12, 0.92, "GENIE 3.0.6");
    if ( generator_name == "GENIE" ) left_title = new TLatex(0.806, 0.92, "GENIE 3.4.0"); 
    if ( generator_name == "NuWro" ) left_title = new TLatex(0.12, 0.92, "NuWro 19.02.1");
    if ( generator_name == "GiBUU" ) left_title = new TLatex(0.12, 0.92, "GiBUU 2019");
    if ( generator_name == "NEUT" ) left_title = new TLatex(0.12, 0.92, "NEUT 5.4.0");
    left_title->SetTextFont(62);
    left_title->SetTextColor(kGray+2);
    left_title->SetNDC();
    left_title->SetTextSize(1/25.);
    left_title->SetTextAlign(10);//left adjusted
    left_title->Draw();


    TLatex* prelim = new TLatex(0.94,0.93, "MicroBooNE Preliminary");
    prelim->SetTextFont(62);
    prelim->SetTextColor(kGray+2);
    prelim->SetNDC();
    prelim->SetTextSize(1/30.);
    prelim->SetTextAlign(32);
    prelim->SetTextSize(0.04631579);
    //prelim->Draw();

    legend->Draw("same");

    //canvas->SaveAs( ("Plots/plot_" + generator_name + "_" + variable + "_" + CUT_TO_USE + ".pdf").c_str() );
    canvas->SaveAs( ("Plots/" + generator_name + "/plot" + std::to_string(p) + ".png").c_str() );
    canvas->SaveAs( ("Plots/" + generator_name + "/plot" + std::to_string(p) + ".root").c_str() );

    f1->cd();
    canvas->Write();

    delete canvas;
    delete legend;
  }
}
