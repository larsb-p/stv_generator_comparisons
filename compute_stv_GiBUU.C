// Plots single transverse variables using entries from a GiBUU-generated RootTuple
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"

constexpr int MUON = 13;
constexpr int PROTON = 2212;
constexpr int PI_PLUS = 211;
constexpr int PI_ZERO = 111;
constexpr int PI_MINUS = -211;
constexpr int NEUTRON = 2112;
constexpr double BOGUS = -1e30;

const double MUON_MASS = TDatabasePDG::Instance()->GetParticle( MUON )->Mass();
const double PROTON_MASS = TDatabasePDG::Instance()->GetParticle( PROTON )->Mass();
const double PI_PLUS_MASS = TDatabasePDG::Instance()->GetParticle( PI_PLUS )->Mass();
const double PI_ZERO_MASS = TDatabasePDG::Instance()->GetParticle( PI_ZERO )->Mass();
const double PI_MINUS_MASS = TDatabasePDG::Instance()->GetParticle( PI_MINUS )->Mass();
const double NEUTRON_MASS = TDatabasePDG::Instance()->GetParticle( NEUTRON )->Mass();
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV (from GENIE v3.0.6)
constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV (from GENIE v3.0.6)

// Nominal threshold values taken from
// https://microboone.fnal.gov/wp-content/uploads/MICROBOONE-NOTE-1024-PUB.pdf
// Thresholds on particle kinetic energy
//constexpr double THRESHOLD_KE_PROTON = 0.082; // GeV
//constexpr double THRESHOLD_KE_PION = 0.037; // GeV
//constexpr double THRESHOLD_KE_MUON = 0.; // GeV

// Thresholds on particle momenta
constexpr double THRESHOLD_P_MUON = .0; // GeV
constexpr double THRESHOLD_P_PROTON = 0.300; // GeV
constexpr double THRESHOLD_P_PION = 0.000; // GeV
constexpr double THRESHOLD_P_NEUTRON = 0.000; // GeV

// Thresholds on particle angle deviated from neutrino direction
//constexpr double THRESHOLD_COSTHETA_MUON = -0.6;
//constexpr double THRESHOLD_COSTHETA_PROTON = 0.4;

void compute_stv_GiBUU(	const std::string& input_file_name,
 	            	const std::string& GiBUU_tune_name)
{

  // Input file
  TFile in_file(input_file_name.c_str(), "read");
//  TTree *ntuple = NULL;
  TTree *ntuple = (TTree* )in_file.Get("RootTuple");

//  in_file.GetObject("ntuple", RootTuple);
//  if ( !ntuple ) return;

  std::vector<int> *barcode = 0;
  std::vector<double> *Px = 0;
  std::vector<double> *Py = 0;
  std::vector<double> *Pz = 0;
  std::vector<double> *E = 0;
  int evType;
  double weight, lepIn_E, lepIn_Px, lepIn_Py, lepIn_Pz, lepOut_E, lepOut_Px, lepOut_Py, lepOut_Pz, nuc_E, nuc_Px, nuc_Py, nuc_Pz; // outgoing lepton 4-momentum

  ntuple->SetBranchAddress("weight", &weight); // weight
  ntuple->SetBranchAddress("barcode", &barcode); // PDG of particle
  ntuple->SetBranchAddress("Px", &Px); // final state particle Px
  ntuple->SetBranchAddress("Py", &Py); // final state particle Py
  ntuple->SetBranchAddress("Pz", &Pz); // final state particle Pz
  ntuple->SetBranchAddress("E", &E); // final state particle energy E
  ntuple->SetBranchAddress("evType", &evType); // Channel (1=QE, 2-31=res ID, 32,33=1pi, 34=DIS, 35,36=2p2h, 37=2pi)
  ntuple->SetBranchAddress("lepIn_E", &lepIn_E); // Incoming lepton (neutrino) energy
  ntuple->SetBranchAddress("lepIn_Px", &lepIn_Px); // Incoming lepton (neutrino) x-momentum
  ntuple->SetBranchAddress("lepIn_Py", &lepIn_Py); // Incoming lepton (neutrino) y-momentum
  ntuple->SetBranchAddress("lepIn_Pz", &lepIn_Pz); // Incoming lepton (neutrino) z-momentum
  ntuple->SetBranchAddress("lepOut_E", &lepOut_E); // Outgoing lepton (muon) energy
  ntuple->SetBranchAddress("lepOut_Px", &lepOut_Px); // Outgoing lepton (muon) x-momentum
  ntuple->SetBranchAddress("lepOut_Py", &lepOut_Py); // Outgoing lepton (muon) y-momentum
  ntuple->SetBranchAddress("lepOut_Pz", &lepOut_Pz); // Outgoing lepton (muon) z-momentum
  ntuple->SetBranchAddress("nuc_E", &nuc_E); // initial state nucleon (neutron) energy
  ntuple->SetBranchAddress("nuc_Px", &nuc_Px); // initial state nucleon (neutron) x-momentum
  ntuple->SetBranchAddress("nuc_Py", &nuc_Py); // initial state nucleon (neutron) y-momentum
  ntuple->SetBranchAddress("nuc_Pz", &nuc_Pz); // initial state nucleon (neutron) z-momentum

  // Output file
  TFile out_file( ("stv_files/stv_out_" + GiBUU_tune_name + ".root").c_str(), "recreate");
  TTree* out_tree = new TTree("stv_tree", "GiBUU STV tree");


  // Kinematic variables
  double	mc_truth_mu_mom, mc_truth_mu_phi, mc_truth_mu_theta, mc_truth_mu_costheta,
 		mc_truth_leading_p_mom, mc_truth_leading_p_phi, mc_truth_leading_p_theta, mc_truth_leading_p_costheta,
                mc_truth_2nd_leading_p_mom, mc_truth_2nd_leading_p_phi, mc_truth_2nd_leading_p_theta, mc_truth_2nd_leading_p_costheta,
                mc_truth_leading_pi_mom, mc_truth_leading_pi_phi, mc_truth_leading_pi_theta, mc_truth_leading_pi_costheta,
                mc_truth_leading_n_mom, mc_truth_leading_n_phi, mc_truth_leading_n_theta, mc_truth_leading_n_costheta;

  // Single Transverse Kinematic Variables
  double pTmag, delta_alpha_T, delta_phi_T, pTmag_2p, pTmag_nN;

  // Initial neutron momentum and variables needed to calculate it
  double pLmag, delta_p, pLmag_noE, delta_p_noE;

  // Incoming neutrino and leading proton kinetic energy
  double E_nu, E_kin_leading_p;

  // Integers to count particles
  int nf_mu_above_threshold, nf_p_above_threshold, nf_pi_above_threshold, nf_neutrons;

  // Flag to mark charged-current events (true if CC event)
  bool cc = false;

  // Tag process channel (QE, 2p2h/MEC, RES, DIS, COH, Other)
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

  // Tag process channel (CC or not
  out_tree->Branch("cc", &cc, "cc/O");

  // Tag process channel (QE, 2p2h/MEC, RES, DIS, COH, Other)
  out_tree->Branch("cat", &cat, "cat/I");

  // STVs ok?
  out_tree->Branch("STVs_ok", &STVs_ok, "STVs_ok/O");

  // Event "perweight"
  out_tree->Branch("weight", &weight, "weight/D");

  // Loop over events
//  for (int e = 0; e < ntuple->GetEntries() - ntuple->GetEntries() + 100000; ++e) {
  for (int e = 0; e < ntuple->GetEntries(); ++e) {

    cout << endl;
    cout << "Event " << e << "                                                           <- New event" << endl;

    if ( e % 1000 == 0 ) std::cout << "Entry " << e << '\n';

    // Get current event
    ntuple->GetEntry(e);

    // Find indices of leading final state protons, pions and neutrons
    int index_leading_N = 0; // Momentum leading proton
    int index_leading_N2 = BOGUS; // Momentum 2nd leading proton
    int index_leading_pi = 0; // Index pion
    int index_leading_n = 0; // Index neutron
    int index_N_multi_proton[1000]; // Indices of all final state protons in one event
    int index_N_multi_pion[1000]; // Indices of all final state pions in one event
    int index_N_multi_neutron[1000]; // Indices of all final state neutrons in one event

    // Find out what kind of event and reset quantities
    cc = false;
    STVs_ok = false;

    double pN_mag_p = 0.; // magnitude of final state leading nucleon (proton) mom
    double pN_mag_p2 = 0.; // magnitude of 2nd leading final state nucleon (proton) mom
    double pN_mag_pi = 0.; // magnitude of final state meson (pion) mom
    double pN_mag_n = 0.; // magnitude of final state nucleon (neutron) mom

    nf_mu_above_threshold = 0; // number of muons above threshold
    nf_p_above_threshold = 0; // number of final state protons above threshold
    nf_pi_above_threshold = 0; // number of final state pions above threshold
    nf_neutrons = 0; // number of final state neutrons

    cout << "Tracks in event " << e << " is " << barcode->size() << endl;

    cc = true; // I want bool variable 'cc' in output file in analogy to NuWro. In NuWro cut is made here, because there not only CC events. In GiBUU one has to specify in jobCard what kind of events. Therefore here I only have CC events. How can I check cc events in GiBUU?

    if ( cc ) { // Beginning of muon check

      // Get muon information
      double plx = lepOut_Px;
      double ply = lepOut_Py;
      double plz = lepOut_Pz;

      // Check whether this muon is above the kinetic energy or momentum threshold
      double KE_mu = lepOut_E - MUON_MASS;

      // Calculate kinematic variables of muon
      mc_truth_mu_mom = std::sqrt( plx*plx + ply*ply + plz*plz );
      mc_truth_mu_costheta = plz / std::sqrt( plx*plx + ply*ply + plz*plz );
      mc_truth_mu_theta = std::acos( mc_truth_mu_costheta );
      mc_truth_mu_phi = std::atan2( ply, plx );

      // Check threshold on energy, momentum and/or angle
//      if ( KE_mu > THRESHOLD_KE_MUON ) {
      if ( mc_truth_mu_mom > THRESHOLD_P_MUON ) {
//      if ( mc_truth_mu_mom > THRESHOLD_P_MUON && costheta_mu > THRESHOLD_COSTHETA_MUON ) {
	++nf_mu_above_threshold;
      }

      cout << "E_kin of muon track is " << KE_mu << endl;
      cout << "Number of muons above threshold of event " << e << " is " << nf_mu_above_threshold << endl;

    } // End of muon check

    // Loop over particles in one event
    for (int f = 0; f < barcode->size(); ++f) { // Beginning of loop over particles in one event

      cout << "--------------------------" << barcode->at(f) << "--------------------------" << endl;

      // compare output with RootTuple->Scan("Px:Py:Pz:E:barcode")
      cout << "Barcode of track " << f+1 << " is " << barcode->at(f) << endl;
      cout << "Px of track " << f+1 << " is " << Px->at(f) << endl;
      cout << "Py of track " << f+1 << " is " << Py->at(f) << endl;
      cout << "Pz of track " << f+1 << " is " << Pz->at(f) << endl;
      cout << "E = E_kin + m of track " << f+1 << " is " << E->at(f) << endl;

      // Get proton information
      if ( barcode->at(f) == PROTON ) { // Beginning of proton check

	// Save all proton indices
	index_N_multi_proton[ nf_p_above_threshold ] = f;

        double pNx = Px->at(f);
        double pNy = Py->at(f);
        double pNz = Pz->at(f);

        // Check whether this proton is above the kinetic energy and/or momentum threshold
        double KE = E->at(f) - PROTON_MASS;
        double new_pN_mag_p = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
        double costheta = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

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

        cout << "Track belongs to a proton (" << barcode->at(f) << ")" << endl;
        cout << "E_kin of track  " << f+1 << " is " << KE << endl;
        cout << "nf_p_above_threshold of event " << e << " is " << nf_p_above_threshold << endl;

      } // End of proton check


      // Get pion information
      else if ( barcode->at(f) == PI_PLUS || barcode->at(f) == PI_ZERO || barcode->at(f) == PI_MINUS ) { // Beginning of pion check

        // Save all pion indices
        index_N_multi_pion[ nf_pi_above_threshold ] = f;

        double pNx = Px->at(f);
        double pNy = Py->at(f);
        double pNz = Pz->at(f);

	double KE;

	// Check whether this pion is above the kinetic energy or momentum threshold
	if (  barcode->at(f) == PI_PLUS  ) KE = E->at(f) - PI_PLUS_MASS;
        if (  barcode->at(f) == PI_ZERO  ) KE = E->at(f) - PI_ZERO_MASS;
        if (  barcode->at(f) == PI_MINUS  ) KE = E->at(f) - PI_MINUS_MASS;

        double new_pN_mag_pi = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//	double costheta_pi = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

	// Check threshold on energy, momentum and/or angle
//	if ( KE > THRESHOLD_KE_PION ) {
        if ( new_pN_mag_pi > THRESHOLD_P_PION ) {
//      if ( pN_mag_pi > THRESHOLD_P_PION && costheta_pi > THRESHOLD_COSTHETA_PION ) {
	  ++nf_pi_above_threshold; // THRESHOLD_P_PION currently set to 0. GeV
	}

        // Get leading pion momentum
        if ( pN_mag_pi < new_pN_mag_pi ) {
	pN_mag_pi = new_pN_mag_pi; // leading pion momentum
	index_leading_pi = f;
	}

	cout << "Track belongs to a pion (" << barcode->at(f) << ")" << endl;
	cout << "E_kin of track  " << f+1 << " is " << KE << endl;
	cout << "nf_pi_above_threshold of event " << e << " is " << nf_pi_above_threshold << endl;

       }


      // Get neutron information
      else if ( barcode->at(f) == NEUTRON ) { // Beginning of neutron check

        // Save all neutron indices
        index_N_multi_neutron[ nf_neutrons ] = f;

        double pNx = Px->at(f);
        double pNy = Py->at(f);
        double pNz = Pz->at(f);

        // Check whether this pion is above the kinetic energy or momentum threshold
	double KE = E->at(f) - NEUTRON_MASS;
        double new_pN_mag_n = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//      double costheta_n = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

        // Check threshold on energy, momentum and/or angle
//	if ( KE > THRESHOLD_KE_NEUTRON ) {
        if ( new_pN_mag_n > THRESHOLD_P_NEUTRON ) {
//      if ( pN_mag_n > THRESHOLD_P_NEUTRON && costheta_n > THRESHOLD_COSTHETA_NEUTRON ) {
	  ++nf_neutrons; // THRESHOLD_P_NEUTRON currently set to 0. GeV
	}

        // Get leading neutron momentum
        if ( pN_mag_n < new_pN_mag_n ) {
	pN_mag_n = new_pN_mag_n; // leading neutron momentum
        index_leading_n = f;
	}

	cout << "Track belongs to a neutron (" << barcode->at(f) << ")" << endl;
	cout << "E_kin of track  " << f+1 << " is " << KE << endl;
	cout << "nf_neutrons of event " << e << " is " << nf_neutrons << endl;

      } // End of neutron check

    } // End of loop over particles in one event

    cout << endl;
    cout << "STVs ok? " << STVs_ok << endl;
    cout << "Number of protons above threshold of event " << e << " is " << nf_p_above_threshold << endl;
    cout << "Number of pions above threshold of event " << e << " is " << nf_pi_above_threshold << endl;
    cout << "Number of neutrons above threshold of event " << e << " is " << nf_neutrons << endl;
    cout << endl;

    // Momentum leading proton kinetic energy
    if ( nf_p_above_threshold > 0 ) {
      E_kin_leading_p = E->at( index_leading_N ) - PROTON_MASS;
    }


    // Leading proton momentum components
    double pN_lead_x = Px->at( index_leading_N );
    double pN_lead_y = Py->at( index_leading_N );
    double pN_lead_z = Pz->at( index_leading_N );

    // Momentum leading proton kinematic variables
    mc_truth_leading_p_mom = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y + pN_lead_z*pN_lead_z );
    mc_truth_leading_p_phi = std::atan2( pN_lead_y, pN_lead_x );
    mc_truth_leading_p_costheta = pN_lead_z / std::sqrt( mc_truth_leading_p_mom );
    mc_truth_leading_p_theta = std::acos( mc_truth_leading_p_costheta );

    // 2nd momentum leading proton momentum components
    double pN_lead_x2, pN_lead_y2, pN_lead_z2;

    if ( nf_p_above_threshold >= 2 ) {
    pN_lead_x2 = Px->at( index_leading_N2 );
    pN_lead_y2 = Py->at( index_leading_N2 );
    pN_lead_z2 = Pz->at( index_leading_N2 );

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
    pN_lead_xpi = Px->at( index_leading_pi );
    pN_lead_ypi = Py->at( index_leading_pi );
    pN_lead_zpi = Pz->at( index_leading_pi );

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
    pN_lead_xn = Px->at( index_leading_n );
    pN_lead_yn = Py->at( index_leading_n );
    pN_lead_zn = Pz->at( index_leading_n );

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
        pTx_np += Px->at( index_N_multi_proton[i] );
        pTy_np += Py->at( index_N_multi_proton[i] );
      }
    }

    // Pions
    if ( nf_pi_above_threshold >= 1 ) {
      for ( int j = 0; j < nf_pi_above_threshold ; ++j) {

        // Add all final (above threshold) pion momenta
        pTx_npi += Px->at( index_N_multi_pion[j] );
        pTy_npi += Py->at( index_N_multi_pion[j] );
      }
    }

    // Neutrons
    if ( nf_neutrons >= 1 ) {
      for ( int j = 0; j < nf_neutrons ; ++j) {

        // Add all final (above threshold) neutrons momenta
        pTx_nn += Px->at( index_N_multi_neutron[j] );
        pTy_nn += Py->at( index_N_multi_neutron[j] );
      }
    }


    // Calculate Single Transverse Variables

    // Muon momentum components
    double plx = lepOut_Px;
    double ply = lepOut_Py;
    double plz = lepOut_Pz;

    // Transverse momentum components of leading proton
    double pTx = lepOut_Px + pN_lead_x;
    double pTy = lepOut_Py + pN_lead_y;

    // Calculate delta_p_T
    pTmag = std::sqrt( pTx*pTx + pTy*pTy );

    // Calculate delta_alpha_T
    double pTl = std::sqrt( lepOut_Px*lepOut_Px + lepOut_Py*lepOut_Py );
    delta_alpha_T = std::acos( ( -lepOut_Px*pTx - lepOut_Py*pTy ) / ( pTl * pTmag ) );

    // Calculate delta_phi_T
    double pTn = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y );
    delta_phi_T = std::acos( ( -lepOut_Px*pN_lead_x - lepOut_Py*pN_lead_y ) / ( pTl * pTn ) );

    // delta pT_2p (Single Transverse Momentum Imbalance with two final state protons)
    double pTx_2p = plx + pN_lead_x + pN_lead_x2;
    double pTy_2p = ply + pN_lead_y + pN_lead_y2;
    if ( nf_p_above_threshold > 1 ) {
    pTmag_2p = std::sqrt( pTx_2p*pTx_2p + pTy_2p*pTy_2p );
    }
    else pTmag_2p = BOGUS;

    // delta_pT_nN (Single Transverse Momentum Imbalance with N final state hadrons (protons, pions, neutrons))
    double pTx_nN = plx + pTx_np + pTx_npi + pTx_nn;
    double pTy_nN = ply + pTy_np + pTy_npi + pTy_nn;
    pTmag_nN = std::sqrt( pTx_nN*pTx_nN + pTy_nN*pTy_nN );


    // STK 4-momentum-imbalance delta p (see arXiv:1805.05486)
    // There are two ways of calculating delta p: 1) when neutrino energy E_nu is known and 2) when it is unknown
    // 1)

    // Incoming neutrino energy
    E_nu = lepIn_E;
    cout << "Incoming Neutrino Energy " << lepIn_E << endl;

    // delta pL
    pLmag = plz + pN_lead_z - lepIn_E;

    // delta p
    delta_p = std::sqrt( pTmag*pTmag + pLmag*pLmag );

    // 2)

    double R = TARGET_MASS + plz + pN_lead_z - lepOut_E - E->at( index_leading_N );

    // Estimated mass of the final remnant nucleus (CCQE assumption)
    double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;

    // delta pL without knowing E_nu
    pLmag_noE = 0.5 * R - ( std::pow( mf, 2 ) + std::pow( pTmag, 2 ) ) / ( 2. * R );

    // delta p without knowing E_nu
    delta_p_noE = std::sqrt( std::pow( pLmag_noE, 2 ) + std::pow( pTmag, 2 ) );


//    cc = true; // I want bool variable 'cc' in output file in analogy to NuWro. In NuWro cut is made here, because there not only CC events. In GiBUU one has to specify in jobCard what kind of events. Therefore here I only have CC events.

    // Flag processes
    if ( evType == 1 ) cat = 1; // Quasielastic
    else if ( evType == 35 || evType == 36 ) cat = 2; // 2p2h
    else if ( 1 < evType && evType < 32 ) cat = 3; // res ID (Resonance production)
    else if ( evType == 34 ) cat = 4; // Deep inelastic Scattering
//    else if ( evType == 1 ) cat = 5; // Coherent meson production
    else cat = 6; // Other

    // Do not consider events where STVs shall not be calculated and consider only CC events
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
  } // End of loop over events

  out_tree->Write();

  // Write the tune name to the ROOT file as metadata for later plotting
  TNamed temp_named( "GiBUU_tune", GiBUU_tune_name.c_str() );
  temp_named.Write();

  out_file.Close();
}
