// Plots single transverse variables using entries from a GENIE-generated
// GST TTree
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

constexpr int MUON = 13;
constexpr int PROTON = 2212;
constexpr int PI_PLUS = 211;
constexpr int PI_MINUS = -211;
constexpr int PI_ZERO = 111;
constexpr int NEUTRON = 2112;
constexpr double BOGUS = -1e30;

const double MUON_MASS = TDatabasePDG::Instance()->GetParticle( MUON )->Mass();
const double PROTON_MASS = TDatabasePDG::Instance()->GetParticle( PROTON )->Mass();
const double PI_PLUS_MASS = TDatabasePDG::Instance()->GetParticle( PI_PLUS )->Mass();
const double PI_ZERO_MASS = TDatabasePDG::Instance()->GetParticle( PI_PLUS )->Mass();
const double PI_MINUS_MASS = TDatabasePDG::Instance()->GetParticle( PI_PLUS )->Mass();
const double NEUTRON_MASS = TDatabasePDG::Instance()->GetParticle( NEUTRON )->Mass();
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV (from GENIE v3.0.6)
constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV (from GENIE v3.0.6)

// Nominal threshold values taken from
// https://microboone.fnal.gov/wp-content/uploads/MICROBOONE-NOTE-1024-PUB.pdf
// Thresholds on particle kinetic energy
constexpr double THRESHOLD_KE_MUON = 0.; // GENIE uses GeV
constexpr double THRESHOLD_KE_PROTON = 0.; // 0.082; // GeV
//constexpr double THRESHOLD_KE_PION = 0.037; // GeV

// Thresholds on particle momenta
constexpr double THRESHOLD_P_MUON = .0;
constexpr double THRESHOLD_P_PROTON = .300; // GeV
constexpr double THRESHOLD_P_PION = 0.000; // GeV
constexpr double THRESHOLD_P_NEUTRON = 0.000; // GeV

// Thresholds on particle angle deviated from neutrino direction
//constexpr double THRESHOLD_COSTHETA_MUON = -0.6;
//constexpr double THRESHOLD_COSTHETA_PROTON = 0.4;


void compute_stv_GENIE(const std::string& input_file_name,
  const std::string& genie_tune_name)
{

  TFile in_file(input_file_name.c_str(), "read");
  TTree* gst = NULL;

  in_file.GetObject("gst", gst);
  if ( !gst ) return;

  int neu; // neutrino PDG code

  // Flag to mark charged-current events (true if CC event)
  bool cc; // true if CC event (will be used in make_plots.C)

  double Ev, El, pxl, pyl, pzl; // outgoing lepton 4-momentum

  int nf; // Number of final hadron(s)
  int ni; // Number of particles in hadronic system befor FSIs
  int nfp; // Number of final state p and anti-p after intranuclear scattering
  int pdgf[1000]; // PDG code for final hadron(s)
  int pdgi[1000]; // PDG code for final hadron(s) before FSIs
  double Ef[1000]; // Total energy of final state hadron
  double pxf[1000]; // x-component of final state hadron momentum
  double pyf[1000]; // y-component of final state hadron momentum
  double pzf[1000]; // z-component of final state hadron momentum

  double pxi[1000]; // x-component of final state hadron momentum befor FSIs
  double pyi[1000]; // y-component of final state hadron momentum befor FSIs
  double pzi[1000]; // z-component of final state hadron momentum befor FSIs

  bool qel, mec, res, dis, coh;

  gst->SetBranchAddress("neu", &neu); // Neutrion PDG code
  gst->SetBranchAddress("cc", &cc); // Flags charged current events

  gst->SetBranchAddress("Ev", &Ev); // Incoming neutrino energy
  gst->SetBranchAddress("El", &El); // Final state primary lepton energy
  gst->SetBranchAddress("pxl", &pxl); // Final state primary lepton px
  gst->SetBranchAddress("pyl", &pyl); // Final state primary lepton py
  gst->SetBranchAddress("pzl", &pzl); // Final state primary lepton pz

  gst->SetBranchAddress("nf", &nf); // Number of final state particles
  gst->SetBranchAddress("ni", &ni); // Number of before FSIs / intranuclear rescattering)
  gst->SetBranchAddress("nfp", &nfp); // Number of final state p and p̄ (after intranuclear rescattering)
  gst->SetBranchAddress("pdgf", &pdgf); // PDG code of k-th final state particle in hadronic system
  gst->SetBranchAddress("pdgi", &pdgi); // PDG code of k-th particle in ‘primary’ hadronic system
  gst->SetBranchAddress("Ef", &Ef); // Energy of k-th final state particle in hadronic system
  gst->SetBranchAddress("pxf", &pxf); // x-component of k-th final state particle in hadronic system
  gst->SetBranchAddress("pyf", &pyf); // y-component of k-th final state particle in hadronic system
  gst->SetBranchAddress("pzf", &pzf); // z-component of k-th final state particle in hadronic system

  gst->SetBranchAddress("qel", &qel); // quasi-elastic
  gst->SetBranchAddress("mec", &mec); // meson exchange current
  gst->SetBranchAddress("res", &res); // resonance neutrino-production
  gst->SetBranchAddress("dis", &dis); // deep-inelastic scattering
  gst->SetBranchAddress("coh", &coh); // coherent meson production

  // Make output .root file
  TFile out_file( ("stv_files/stv_out_" + genie_tune_name + ".root").c_str(), "recreate");
  TTree* out_tree = new TTree("stv_tree", "GENIE STV tree");

  // Kinematic variables
  double        mc_truth_mu_mom, mc_truth_mu_phi, mc_truth_mu_theta, mc_truth_mu_costheta,
                mc_truth_leading_p_mom, mc_truth_leading_p_phi, mc_truth_leading_p_theta, mc_truth_leading_p_costheta,
                mc_truth_2nd_leading_p_mom, mc_truth_2nd_leading_p_phi, mc_truth_2nd_leading_p_theta, mc_truth_2nd_leading_p_costheta,
                mc_truth_leading_pi_mom, mc_truth_leading_pi_phi, mc_truth_leading_pi_theta, mc_truth_leading_pi_costheta,
                mc_truth_leading_n_mom, mc_truth_leading_n_phi, mc_truth_leading_n_theta, mc_truth_leading_n_costheta;

  // Single Transverse Variables
  double pTmag, delta_alpha_T, delta_phi_T, pTmag_2p, pTmag_nN;

  // Initial neutron momentum and variables needed to calculate it
  double pLmag, delta_p, pLmag_noE, delta_p_noE;

  // Incoming neutrino and leading proton kinetic energy
  double E_nu, E_kin_leading_p;

  // Integers to count particles
  int nf_mu_above_threshold, nf_p_above_threshold, nf_pi_above_threshold, nf_neutrons;

  // Flag process channel (QE, 2p2h/MEC, RES, DIS, COH, Other)
  int cat;

  // Flag to mark whether STVs could be calculated for the current event
  bool STVs_ok = false;

  // Output TTree variables
  // Muon kinematic variables
  out_tree->Branch("mc_truth_mu_mom", &mc_truth_mu_mom, "mc_truth_mu_mom/D");
  out_tree->Branch("mc_truth_mu_phi", &mc_truth_mu_phi, "mc_truth_mu_phi/D");
  out_tree->Branch("mc_truth_mu_theta", &mc_truth_mu_theta, "mc_truth_mu_theta/D");
  out_tree->Branch("mc_truth_mu_costheta", &mc_truth_mu_costheta, "mc_truth_mu_costheta/D");

  // Momentum leading proton kinematic variables
  out_tree->Branch("mc_truth_leading_p_mom", &mc_truth_leading_p_mom, "mc_truth_leading_p_mom/D");
  out_tree->Branch("mc_truth_leading_p_phi", &mc_truth_leading_p_phi, "mc_truth_leading_p_phi/D");
  out_tree->Branch("mc_truth_leading_p_theta", &mc_truth_leading_p_theta, "mc_truth_leading_p_theta/D");
  out_tree->Branch("mc_truth_leading_p_costheta", &mc_truth_leading_p_costheta, "mc_truth_leading_p_costheta/D");

  // 2nd momentum leading proton kinematic variables
  out_tree->Branch("mc_truth_2nd_leading_p_mom", &mc_truth_2nd_leading_p_mom, "mc_truth_2nd_leading_p_mom/D");
  out_tree->Branch("mc_truth_2nd_leading_p_phi", &mc_truth_2nd_leading_p_phi, "mc_truth_2nd_leading_p_phi/D");
  out_tree->Branch("mc_truth_2nd_leading_p_theta", &mc_truth_2nd_leading_p_theta, "mc_truth_2nd_leading_p_theta/D");
  out_tree->Branch("mc_truth_2nd_leading_p_costheta", &mc_truth_2nd_leading_p_costheta, "mc_truth_2nd_leading_p_costheta/D");

  // Momentum leading pion kinematic variables
  out_tree->Branch("mc_truth_leading_pi_mom", &mc_truth_leading_pi_mom, "mc_truth_leading_pi_mom/D");
  out_tree->Branch("mc_truth_leading_pi_phi", &mc_truth_leading_pi_phi, "mc_truth_leading_pi_phi/D");
  out_tree->Branch("mc_truth_leading_pi_theta", &mc_truth_leading_pi_theta, "mc_truth_leading_pi_theta/D");
  out_tree->Branch("mc_truth_leading_pi_costheta", &mc_truth_leading_pi_costheta, "mc_truth_leading_pi_costheta/D");

  // Momentum leading neutron kinematic variables
  out_tree->Branch("mc_truth_leading_n_mom", &mc_truth_leading_n_mom, "mc_truth_leading_n_mom/D");
  out_tree->Branch("mc_truth_leading_n_phi", &mc_truth_leading_n_phi, "mc_truth_leading_n_phi/D");
  out_tree->Branch("mc_truth_leading_n_theta", &mc_truth_leading_n_theta, "mc_truth_leading_n_theta/D");
  out_tree->Branch("mc_truth_leading_n_costheta", &mc_truth_leading_n_costheta, "mc_truth_leading_n_costheta/D");

  // Single Transverse Kinematic Variables
  out_tree->Branch("delta_p_T", &pTmag, "delta_p_T/D");
  out_tree->Branch("delta_alpha_T", &delta_alpha_T, "delta_alpha_T/D");
  out_tree->Branch("delta_phi_T", &delta_phi_T, "delta_phi_T/D");
  out_tree->Branch("delta_p_T_2p", &pTmag_2p, "delta_p_T_2p/D");
  out_tree->Branch("delta_p_T_nN", &pTmag_nN, "delta_p_T_nN/D");

  // Initial neutron momentum and variables needed to calculate it
  out_tree->Branch("delta_p_L", &pLmag, "delta_p_L/D");
  out_tree->Branch("delta_p", &delta_p, "delta_p/D");
  out_tree->Branch("delta_p_L_noE", &pLmag_noE, "delta_p_L_noE/D");
  out_tree->Branch("delta_p_noE", &delta_p_noE, "delta_p_noE/D");

  // Incoming neutrino energy
  out_tree->Branch("E_nu", &E_nu, "E_nu/D");

  // Leading proton kinetic energy
  out_tree->Branch("E_kin_leading_p", &E_kin_leading_p, "E_kin_leading_p/D");

  // Integers to count particles
  out_tree->Branch("nf_mu_above_threshold", &nf_mu_above_threshold, "nf_mu_above_threshold/I");
  out_tree->Branch("nf_p_above_threshold", &nf_p_above_threshold, "nf_p_above_threshold/I");
  out_tree->Branch("nf_pi_above_threshold", &nf_pi_above_threshold, "nf_pi_above_threshold/I");
  out_tree->Branch("nf_neutrons", &nf_neutrons, "nf_neutrons/I");

  // Tag process channel (CC or not)
  out_tree->Branch("cc", &cc, "cc/O");

  // Tag process channel (QE, 2p2h/MEC, RES, DIS, COH, Other)
  out_tree->Branch("cat", &cat, "cat/I");

  // STVs ok?
  out_tree->Branch("STVs_ok", &STVs_ok, "STVs_ok/O");

  // Attach histogram with leading proton kinetic energy to output file 
  Double_t bins[] = { 0., THRESHOLD_KE_PROTON, 5.5 };
  Int_t  binnum = sizeof(bins)/sizeof(Double_t) - 1;
  TH1D * h_Ekin_leading_p = new TH1D ("h_Ekin_leading_p", ";Kinetic energy of leading proton [GeV];Entries", binnum , bins );
  h_Ekin_leading_p->SetCanExtend(TH1::kXaxis);
  h_Ekin_leading_p->SetTitle( ( genie_tune_name ).c_str() ) ;


  // Loop over events
  for (int e = 0; e < gst->GetEntries(); ++e) {

//    cout << endl;
//    cout << "Event " << e << "                                                           <- New event" << endl;

    if ( e % 1000 == 0 ) std::cout << "Entry " << e << '\n';

    // Get current event
    gst->GetEntry( e );

    // Find indices of leading final state protons, pions and neutrons
    int index_leading_N = 0; // Momentum leading proton
    int index_leading_N2 = BOGUS; // Momentum 2nd leading proton
    int index_leading_pi = 0; // Index pion
    int index_leading_n = 0; // Index neutron
    int index_N_multi_proton[1000]; // Indices of all final state protons in one event
    int index_N_multi_pion[1000]; // Indices of all final state pions in one event
    int index_N_multi_neutron[1000]; // Indices of all final state neutrons in one event

    // Reset
    STVs_ok = false;

//    cout << "Event " << e << " has " << ni << " particles before FSIs and " << nf << " particles that leave nucleus (after FSIs) and " << nf_p << " proton after intranuclear rescattering" << '\n';

    double pN_mag_p = 0.; // magnitude of final state nucleon (proton) mom
    double pN_mag_p2; // magnitude of 2nd leading final state nucleon (proton) mom
    double pN_mag_pi = 0.; // magnitude of final state meson (pion) mom
    double pN_mag_n = 0.; // magnitude of final state nucleon (neutron) mom

    nf_mu_above_threshold = 0; // number of muons above threshold
    nf_p_above_threshold = 0; // number of final state protons above threshold
    nf_pi_above_threshold = 0; // number of final state pions above threshold 
    nf_neutrons = 0; // number of final state neutrons

    if ( cc ) { // Beginning of muon check (if nu_mu CC event there must be a muon in the final state )

      // Get muon information
      // Check whether this muon is above the kinetic energy or momentum threshold
      double KE_mu = El - MUON_MASS;

      // Calculate kinematic variables of muon
     mc_truth_mu_mom = std::sqrt( pxl*pxl + pyl*pyl + pzl*pzl );
      mc_truth_mu_costheta = pzl / std::sqrt( pxl*pxl + pyl*pyl + pzl*pzl );
      mc_truth_mu_theta = std::acos( mc_truth_mu_costheta );
      mc_truth_mu_phi = std::atan2( pyl, pxl );

      // Check threshold on energy, momentum and/or angle
//    if ( KE_mu > THRESHOLD_KE_MUON ) {
      if ( cc && mc_truth_mu_mom > THRESHOLD_P_MUON ) {
//    if ( mc_truth_mu_mom > THRESHOLD_P_MUON && costheta_mu > THRESHOLD_COSTHETA_MUON ) {
	++nf_mu_above_threshold;
      }

      cout << "E_kin of muon track is " << KE_mu;
      cout << " and costheta " << mc_truth_mu_costheta << '\n';
      cout << " pxl is " << pxl << '\n';
      cout << " pyl is " << pyl << '\n';
      cout << " pzl is " << pzl << '\n';
      cout << "Number of muons above threshold of event " << e << " is " << nf_mu_above_threshold << endl;

    } // End of muon check

    // Loop over particles in one event
    for (int f = 0; f < nf; ++f) { // Beginning of loop over particles in one event

//      cout << "Before FSIs: Particle " << f+1 << " has PDG " << pdgi[f] << " and momentum = " << std::sqrt( pxi[f] * pxi[f] + pyi[f] * pyi[f] + pzi[f] * pzi[f] ) << ";  " << '\n';
//      cout << "After FSIs: Particle " << f+1 << " has PDG " << pdgf[f] << " and momentum = " << std::sqrt( pxf[f] * pxf[f] + pyf[f] * pyf[f] + pzf[f] * pzf[f] ) << '\n';

      // Get proton information
      if ( pdgf[f] == PROTON ) {

        // Save all proton indices
        index_N_multi_proton[ nf_p_above_threshold ] = f;

	double pNx = pxf[f];
	double pNy = pyf[f];
	double pNz = pzf[f];
       
	// Check whether this proton is above the kinetic energy or momentum threshold
        double KE = Ef[f] - PROTON_MASS;
        double new_pN_mag_p = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
        double costheta = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

//	cout << " and costheta " << costheta << '\n';

        // Check threshold on energy, momentum and/or angle
//	if ( KE > THRESHOLD_KE_PROTON ) {
	if ( new_pN_mag_p > THRESHOLD_P_PROTON ) {
//	if ( new_pN_mag_p > THRESHOLD_P_PROTON && costheta > THRESHOLD_COSTHETA_PROTON ) {
	  ++nf_p_above_threshold;
          STVs_ok = true; // There is a leading proton with momentum exceeding threshold -> STVs can be calculated
	}

	// Get leading and second leading proton momenta
	if ( pN_mag_p < new_pN_mag_p ) {
	  pN_mag_p2 = pN_mag_p; // first save second to leading proton momentum...
	  index_leading_N2 = index_leading_N; // ...and second to leading proton index 
	  pN_mag_p = new_pN_mag_p; // then overwrite leading proton momentum...
	  index_leading_N = f; // ...and its index
	}

	if ( new_pN_mag_p > THRESHOLD_P_PROTON && new_pN_mag_p < pN_mag_p && new_pN_mag_p > pN_mag_p2 ) {
	  pN_mag_p2 = new_pN_mag_p; // first save second to leading proton momentum...
	  index_leading_N2 = f; // ...and second to leading proton index          
	}

//	cout << "and total momentum " << pN_mag_p << '\n';

     } // End of proton check 


      // Get pion information
      else if ( std::abs(pdgf[f]) == PI_PLUS || std::abs(pdgf[f]) == PI_MINUS || std::abs(pdgf[f]) == PI_ZERO ) {

	// Save all pion indices
	index_N_multi_pion[ nf_pi_above_threshold ] = f;

	double pNx = pxf[f];
	double pNy = pyf[f];
	double pNz = pzf[f];

	double KE;

	// Check whether this pion is above the kinetic energy or momentum threshold
	if ( std::abs(pdgf[f]) == PI_PLUS ) KE = Ef[f] - PI_PLUS_MASS;
        if ( std::abs(pdgf[f]) == PI_ZERO ) KE = Ef[f] - PI_ZERO_MASS;
        if ( std::abs(pdgf[f]) == PI_MINUS ) KE = Ef[f] - PI_MINUS_MASS;

	double new_pN_mag_pi = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//	double costheta = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

//	if ( KE > THRESHOLD_KE_PION ) {
	if ( new_pN_mag_pi > THRESHOLD_P_PION ) {
//	if ( pN_mag_pi > THRESHOLD_P_PION && costheta_pi > THRESHOLD_COSTHETA_PION ) {
	  ++nf_pi_above_threshold;
	}

	// Get leading pion momentum
	if ( pN_mag_pi < new_pN_mag_pi ) {
	  pN_mag_pi = new_pN_mag_pi; // leading pion momentum
	  index_leading_pi = f;
	}

	  cout << "Track belongs to a pion (" << std::abs(pdgf[f]) << ")" << endl;
	  cout << "E_kin of track  " << f+1 << " is " << KE << endl;
	  cout << "nf_pi_above_threshold of event " << e << " is " << nf_pi_above_threshold << endl;

	} // End of pion check


      // Get neutron information
      else if ( std::abs(pdgf[f]) == NEUTRON ) { // Beginning of neutron check

	// Save all neutron indices
	index_N_multi_neutron[ nf_neutrons ] = f;

        double pNx = pxf[f];
        double pNy = pyf[f];
        double pNz = pzf[f];

	// Check whether this pion is above the kinetic energy or momentum threshold
	double KE = Ef[f] - NEUTRON_MASS;
	double new_pN_mag_n = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//	double costheta_n = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

        // Check threshold on energy, momentum and/or angle
//      if ( KE > THRESHOLD_KE_NEUTRON ) {
        if ( new_pN_mag_n > THRESHOLD_P_NEUTRON ) {
//      if ( pN_mag_n > THRESHOLD_P_NEUTRON && costheta_n > THRESHOLD_COSTHETA_NEUTRON ) {
          ++nf_neutrons; // THRESHOLD_P_NEUTRON currently set to 0. GeV
        }

        // Get leading neutron momentum
        if ( pN_mag_n < new_pN_mag_n ) {
        pN_mag_n = new_pN_mag_n; // leading neutron momentum
        index_leading_n = f;
        }

        cout << "Track belongs to a neutron (" << std::abs(pdgf[f]) << ")" << endl;
        cout << "E_kin of track  " << f+1 << " is " << KE << endl;
        cout << "nf_neutrons of event " << e << " is " << nf_neutrons << endl;

      } // End of neutron check

    } // End of loop over particles in one event

//    cout << endl;
//    cout << "STVs ok? " << STVs_ok << endl;
//    cout << "Number of protons above threshold of event " << e << " is " << nf_p_above_threshold << endl;
//    cout << "Number of pions above threshold of event " << e << " is " << nf_pi_above_threshold << endl;
//    cout << "Number of neutrons above threshold of event " << e << " is " << nf_neutrons << endl;
//    cout << endl;

    // Momentum leading proton kinetic energy
    if ( nf_p_above_threshold > 0 ) {
      h_Ekin_leading_p->Fill( Ef[index_leading_N] - PROTON_MASS );
      E_kin_leading_p = Ef[index_leading_N] - PROTON_MASS;
    }

    // Leading proton momentum components
    double pN_lead_x = pxf[ index_leading_N ];
    double pN_lead_y = pyf[ index_leading_N ];
    double pN_lead_z = pzf[ index_leading_N ];

//    cout << "pN_lead_x is " << pN_lead_x << endl;
//    cout << "pN_lead_y is " << pN_lead_y << endl;
//    cout << "pN_lead_z is " << pN_lead_z << endl;

    // Momentum leading proton kinematic variables
    mc_truth_leading_p_mom = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y + pN_lead_z*pN_lead_z );
    mc_truth_leading_p_phi = std::atan2( pN_lead_y, pN_lead_x );
    mc_truth_leading_p_costheta = pN_lead_z / std::sqrt( mc_truth_leading_p_mom );
    mc_truth_leading_p_theta = std::acos( mc_truth_leading_p_costheta );

    // 2nd momentum leading proton momentum components
    double pN_lead_x2, pN_lead_y2, pN_lead_z2;

    if ( nf_p_above_threshold >= 2 ) {
    pN_lead_x2 = pxf[ index_leading_N2 ];
    pN_lead_y2 = pyf[ index_leading_N2 ];
    pN_lead_z2 = pzf[ index_leading_N2 ];

//    cout << "pN_lead_x2 is " << pN_lead_x2 << endl;
//    cout << "pN_lead_y2 is " << pN_lead_y2 << endl;
//    cout << "pN_lead_z2 is " << pN_lead_z2 << endl;

    // 2nd momentum leading proton kinematic variables
    mc_truth_2nd_leading_p_mom = std::sqrt( pN_lead_x2*pN_lead_x2 + pN_lead_y2*pN_lead_y2 + pN_lead_z2*pN_lead_z2 );
    mc_truth_2nd_leading_p_phi = std::atan2( pN_lead_y2, pN_lead_x2 );
    mc_truth_2nd_leading_p_costheta = pN_lead_z2 / std::sqrt( mc_truth_2nd_leading_p_mom );
    mc_truth_2nd_leading_p_theta = std::acos( mc_truth_2nd_leading_p_costheta );
    }
    else {
    pN_lead_x2 = BOGUS;
    pN_lead_y2 = BOGUS;
    pN_lead_z2 = BOGUS;
    mc_truth_2nd_leading_p_mom = BOGUS;
    mc_truth_2nd_leading_p_phi = BOGUS;
    mc_truth_2nd_leading_p_costheta = BOGUS;
    mc_truth_2nd_leading_p_theta = BOGUS;
    }

    // Pion momentum components
    double pN_lead_xpi, pN_lead_ypi, pN_lead_zpi;

    if ( nf_pi_above_threshold > 0 ) {
    pN_lead_xpi = pxf[ index_leading_pi ];
    pN_lead_ypi = pyf[ index_leading_pi ];
    pN_lead_zpi = pzf[ index_leading_pi ];

    // Pion kinematic variables
    mc_truth_leading_pi_mom = std::sqrt( pN_lead_xpi*pN_lead_xpi + pN_lead_ypi*pN_lead_ypi + pN_lead_zpi*pN_lead_zpi );
    mc_truth_leading_pi_phi = std::atan2( pN_lead_ypi, pN_lead_xpi );
    mc_truth_leading_pi_costheta = pN_lead_zpi / std::sqrt( mc_truth_leading_pi_mom );
    mc_truth_leading_pi_theta = std::acos( mc_truth_leading_pi_costheta );
    }
    else {
    pN_lead_xpi = BOGUS;
    pN_lead_ypi = BOGUS;
    pN_lead_zpi = BOGUS;
    mc_truth_leading_pi_mom = BOGUS;
    mc_truth_leading_pi_phi = BOGUS;
    mc_truth_leading_pi_costheta = BOGUS;
    mc_truth_leading_pi_theta = BOGUS;
    }

    // Neutron momentum components
    double pN_lead_xn, pN_lead_yn, pN_lead_zn;

    if ( nf_neutrons > 0 ) {
    pN_lead_xn = pxf[ index_leading_n ];
    pN_lead_yn = pyf[ index_leading_n ];
    pN_lead_zn = pzf[ index_leading_n ];

    // Neutron kinematic variables
    mc_truth_leading_n_mom = std::sqrt( pN_lead_xn*pN_lead_xn + pN_lead_yn*pN_lead_yn + pN_lead_zn*pN_lead_zn );
    mc_truth_leading_n_phi = std::atan2( pN_lead_yn, pN_lead_xn );
    mc_truth_leading_n_costheta = pN_lead_zn / std::sqrt( mc_truth_leading_n_mom );
    mc_truth_leading_n_theta = std::acos( mc_truth_leading_n_costheta );
    }
    else {
    pN_lead_xn = BOGUS;
    pN_lead_yn = BOGUS;
    pN_lead_zn = BOGUS;
    mc_truth_leading_n_mom = BOGUS;
    mc_truth_leading_n_phi = BOGUS;
    mc_truth_leading_n_costheta = BOGUS;
    mc_truth_leading_n_theta = BOGUS;
    }


    // Transverse momentum of all selected hadrons
    double pTx_np = 0.;
    double pTy_np = 0.;

    double pTx_npi = 0.;
    double pTy_npi = 0.;

    double pTx_nn = 0.;
    double pTy_nn = 0.;

    // Proton
    if ( nf_p_above_threshold >= 1 ) {
      for ( int i = 0; i < nf_p_above_threshold ; ++i) {

        // Add all final (above threshold) proton momenta
        pTx_np += pxf[ index_N_multi_proton[i] ];
        pTy_np += pxf[ index_N_multi_proton[i] ];
      }
    }

    // Pions
    if ( nf_pi_above_threshold >= 1 ) {
      for ( int i = 0; i < nf_pi_above_threshold ; ++i) {

        // Add all final (above threshold) pion momenta
        pTx_npi += pxf[ index_N_multi_pion[i] ];
        pTy_npi += pxf[ index_N_multi_pion[i] ];
      }
    }

    // Neutrons
    if ( nf_neutrons >= 1 ) {
      for ( int i = 0; i < nf_neutrons ; ++i) {

        // Add all final (above threshold) proton momenta
        pTx_nn += pxf[ index_N_multi_neutron[i] ];
        pTy_nn += pyf[ index_N_multi_neutron[i] ];
      }
    }


    // Calculate Single Transverse Variables

    // Calculate delta_pT
    double pTx = pxl + pN_lead_x;
    double pTy = pyl + pN_lead_y;
    pTmag = std::sqrt( pTx*pTx + pTy*pTy );

    // Calculate delta_alpha_T
    double pTl = std::sqrt( pxl*pxl + pyl*pyl );
    delta_alpha_T = std::acos( ( -pxl*pTx - pyl*pTy ) / ( pTl * pTmag ) );

    // Calculate delta_phi_T
    double pTn = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y );
    delta_phi_T = std::acos( ( -pxl*pN_lead_x - pyl*pN_lead_y ) / ( pTl * pTn ) );

    // Calculate delta_pT_2p (Single Transverse Momentum Imbalance with two final state protons)
    double pTx_2p = pxl + pN_lead_x + pN_lead_x2;
    double pTy_2p = pyl + pN_lead_y + pN_lead_y2;
    if ( nf_p_above_threshold > 1 ) {
    pTmag_2p = std::sqrt( pTx_2p*pTx_2p + pTy_2p*pTy_2p );
    }
    else pTmag_2p = BOGUS;

    // delta_pT_nN (Single Transverse Momentum Imbalance with N final state hadrons (protons, pions, neutrons))
    double pTx_nN = pxl/* + pTx_np + pTx_npi + pTx_nn*/;
    double pTy_nN = pyl /*+ pTy_np + pTy_npi + pTy_nn*/;
    pTmag_nN = std::sqrt( pTx_nN*pTx_nN + pTy_nN*pTy_nN );

//    std::cout << "Event " << e << " has 1mu" << nf_p_above_threshold << "p summed transverse momentum " << pTmag_nN <<  '\n'; // without pions
//    std::cout << "Event " << e << " has 1mu" << nf_p_above_threshold << "p" << nf_pi_above_threshold << "pi" << " summed transverse momentum " << pTmag_nN <<  '\n'; // with pions

    // STK 4-momentum-imbalance delta p (see arXiv:1805.05486) 
    // There are two ways of calculating delta p: 1) when neutrino energy E_nu is known and 2) when it is unknown
    // 1)

    // Incoming neutrino energy
    E_nu = Ev;
    cout << "Incoming Neutrino Energy " << Ev << endl;

    // delta pL
    pLmag = pzl + pN_lead_z - Ev;

//    cout << "Ev " << Ev << " pzl " << pzl << " pN_lead_z " << pN_lead_z << " pLmag " << pLmag <<endl;

    // delta_p (STK 4-momentum-imbalance)
    delta_p = std::sqrt( pTmag*pTmag + pLmag*pLmag );

    // 2)

    double R = TARGET_MASS + pzl + pN_lead_z - El - ( E_kin_leading_p + PROTON_MASS );

    // Estimated mass of the final remnant nucleus (CCQE assumption)
    double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;

    // delta pL without knowing E_nu
    pLmag_noE = 0.5 * R - ( std::pow( mf, 2 ) + std::pow( pTmag, 2 ) ) / ( 2. * R );

    // delta p without knowing E_nu
    delta_p_noE = std::sqrt( std::pow( pLmag_noE, 2 ) + std::pow( pTmag, 2 ) );

    // Flag processes
    if ( qel ) cat = 1; // Quasielastic
    else if ( mec ) cat = 2; // MEC
    else if ( res ) cat = 3; // Resonance production
    else if ( dis ) cat = 4; // Deep Inelastic Scattering
    else if ( coh ) cat = 5; // Coherent meson production
    else cat = 6; // Other

//    if ( !STVs_ok || !muon_above_thresholds || !proton_above_thresholds ) {
    if ( !STVs_ok || !cc ) {
      mc_truth_mu_mom = BOGUS;
      mc_truth_mu_phi = BOGUS;
      mc_truth_mu_theta = BOGUS;
      mc_truth_mu_costheta = BOGUS;
      mc_truth_leading_p_mom = BOGUS;
      mc_truth_leading_p_phi = BOGUS;
      mc_truth_leading_p_theta = BOGUS;
      mc_truth_leading_p_costheta = BOGUS;
      mc_truth_2nd_leading_p_mom = BOGUS;
      mc_truth_2nd_leading_p_phi = BOGUS;
      mc_truth_2nd_leading_p_theta = BOGUS;
      mc_truth_2nd_leading_p_costheta = BOGUS;
      mc_truth_leading_pi_mom = BOGUS;
      mc_truth_leading_pi_phi = BOGUS;
      mc_truth_leading_pi_theta = BOGUS;
      mc_truth_leading_pi_costheta = BOGUS;
      mc_truth_leading_n_mom = BOGUS;
      mc_truth_leading_n_phi = BOGUS;
      mc_truth_leading_n_theta = BOGUS;
      mc_truth_leading_n_costheta = BOGUS;
      pTmag = BOGUS;
      delta_alpha_T = BOGUS;
      delta_phi_T = BOGUS;
      pTmag_2p = BOGUS;
      pTmag_nN = BOGUS;
      pLmag = BOGUS;
      delta_p = BOGUS;
      pLmag_noE = BOGUS;
      delta_p_noE = BOGUS;
    }

    if ( nf_p_above_threshold == 0 ) {
      E_kin_leading_p = BOGUS;
    }

     out_tree->Fill();

  } // End of event loop

  out_tree->Write();

  // Write the tune name to the ROOT file as metadata for later plotting
  TNamed temp_named( "genie_tune", genie_tune_name.c_str() );
  temp_named.Write();


  out_file.Write(); 
  out_file.Close();
}
