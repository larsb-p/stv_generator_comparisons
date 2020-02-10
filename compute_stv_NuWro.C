// Plots single transverse variables using entries from a NuWro-generated
// TTree called 'treeout'
//#include <iostream>
//#include <string>

//#include "TFile.h"
//#include "TTree.h"

//#include "src/vec.h" // (for 3-D vectors)
//#include "src/vect.h" // (for 4-D vectors)
//#include "src/particle.h"
//#include "src/params.h"

constexpr int MUON = 13;
constexpr int PROTON = 2212;
constexpr int PI_PLUS = 211;
constexpr int PI_ZERO = 111;
constexpr int PI_MINUS = -211;
constexpr int NEUTRON = 2112;
constexpr double BOGUS = -1e30;

//const double MUON_MASS = TDatabasePDG::Instance()->GetParticle( MUON )->Mass();
const double PROTON_MASS = TDatabasePDG::Instance()->GetParticle( PROTON )->Mass();
//const double PI_PLUS_MASS = TDatabasePDG::Instance()->GetParticle( PI_PLUS )->Mass();
const double NEUTRON_MASS = TDatabasePDG::Instance()->GetParticle( NEUTRON )->Mass();
constexpr double TARGET_MASS = 37.215526; // 40Ar, GeV (from GENIE v3.0.6)
constexpr double BINDING_ENERGY = 0.0295; // 40Ar, GeV (from GENIE v3.0.6)

// Nominal threshold values taken from
// https://microboone.fnal.gov/wp-content/uploads/MICROBOONE-NOTE-1024-PUB.pdf
// Thresholds on particle kinetic energy
constexpr double THRESHOLD_KE_MUON = 0.; // MeV
constexpr double THRESHOLD_KE_PROTON = 82.; // NuWro uses MeV
constexpr double THRESHOLD_KE_PION = 37.; // MeV

// Thresholds on particle momenta
constexpr double THRESHOLD_P_MUON = 250.;
//constexpr double THRESHOLD_P_PROTON = 450.; // NuWro uses MeV and c = 1
constexpr double THRESHOLD_P_PROTON = 300.; // NuWro uses MeV and c = 1
constexpr double THRESHOLD_P_PION = 0.000; // GeV
constexpr double THRESHOLD_P_NEUTRON = 0.000; // GeV

// Thresholds on particle angle deviated from neutrino direction
//constexpr double THRESHOLD_COSTHETA_MUON = -0.6;
//constexpr double THRESHOLD_COSTHETA_PROTON = 0.4;


void compute_STVs(const std::string& input_file_name,
                  const std::string& nuwro_tune_name)
{

  TFile in_file(input_file_name.c_str(), "read");
  TTree* treeout = NULL;

  // Declare input used from NuWro event generated tree. event is a class with members, e. g. the vector <particle> in. particle is another class (subclass of vect) and has a function/method p(), that gives out the 3-mom vector by using the class vec .
  in_file.GetObject("treeout", treeout);
  if ( !treeout ) return;

  event *e = new event();

  treeout->SetBranchAddress("e", &e);

  // Make output .root file
  TFile out_file( ("stv_out_" + nuwro_tune_name + ".root").c_str(), "recreate");
  TTree* out_tree = new TTree("stv_tree", "NuWro STV tree");

  // Variables that will be written to out_tree (variables used for calculation of STVs will be initiated in loop.)

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

  // Flag to mark charged-current events (true if CC event)
  bool cc = false;

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


  // Loop over events
  for (int k = 0; k < treeout->GetEntries(); ++k ) {

    cout << endl;
    cout << "Event " << k << "                                                           <- New event" << endl;

    if ( k % 1000 == 0 ) std::cout << "Entry " << k << '\n';

    // Get current event  
    treeout->GetEntry( k );

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

    double pN_mag_p = 0.; // magnitude of final state nucleon (proton) mom
    double pN_mag_p2 = 0.; // magnitude of 2nd leading final state nucleon (proton) mom
    double pN_mag_pi = 0.; // magnitude of final state meson (pion) mom
    double pN_mag_n = 0.; // magnitude of final state nucleon (neutron) mom

    nf_mu_above_threshold = 0; // number of muons above threshold
    nf_p_above_threshold = 0; // number of final state protons above threshold
    nf_pi_above_threshold = 0; // number of final state pions above threshold 
    nf_neutrons = 0; // number of final state neutrons

    cout << "Event " << k << " has " << e->f() << " particles that leave nucleus. " << e->nof(PROTON) << " proton(s) come out of primary vertex (before FSIs) and " << e->fof(PROTON) << " proton(s) come out of nucleus (after FSIs)" << '\n';

    // Check if event is CC event
    if ( e->flag.cc ) cc = true;

    if ( cc && std::abs(e->post[0].pdg) == MUON ) { // Beginning of muon check (if nu_mu CC event there must be a muon in the final state )

      // Get muon information
      double plx = e->post[0].p().x;
      double ply = e->post[0].p().y;
      double plz = e->post[0].p().z;

      // Check whether this muon is above the kinetic energy or momentum threshold
      double KE_mu = e->post[0].Ek();

      // Calculate kinematic variables of muon
      mc_truth_mu_mom = std::sqrt( plx*plx + ply*ply + plz*plz ) / 1000.;
      mc_truth_mu_costheta = plz / std::sqrt( plx*plx + ply*ply + plz*plz );
      mc_truth_mu_theta = std::acos( mc_truth_mu_costheta );
      mc_truth_mu_phi = std::atan2( ply, plx );

      // Check threshold on energy, momentum and/or angle
//    if ( KE_mu > THRESHOLD_KE_MUON ) {
      if ( mc_truth_mu_mom > THRESHOLD_P_MUON / 1000. ) {
//    if ( mc_truth_mu_mom > THRESHOLD_P_MUON && costheta_mu > THRESHOLD_COSTHETA_MUON ) {
	++nf_mu_above_threshold;
      }


      cout << "E_kin of muon track is " << KE_mu;
      cout << " and costheta " << mc_truth_mu_costheta << '\n';
      cout << " plx is " << plx << '\n';
      cout << " ply is " << ply << '\n';
      cout << " plz is " << plz << '\n';
      cout << "Number of muons above threshold of event " << k << " is " << nf_mu_above_threshold << endl;
    
    } // End of muon check


    // Loop over particles in one event
    for ( int l = 0; l < e->f(); ++l ) { // f() gives out number of particles that leave nucleus

      cout << "Particle " << l+1 << " has PDG " << e->post[l].pdg << " and E_kin = " << e->post[l].Ek() << " and momentum = " << std::sqrt( e->post[l].p().x * e->post[l].p().x + e->post[l].p().y * e->post[l].p().y + e->post[l].p().z * e->post[l].p().z ) << '\n';

      // Get proton information
      if ( e->post[l].pdg == PROTON ) { // Beginning of proton check

	// Save all proton indices
	index_N_multi_proton[ nf_p_above_threshold ] = l;

	double pNx = e->post[l].p().x;
	double pNy = e->post[l].p().y;
	double pNz = e->post[l].p().z;

	// Check whether this proton is above the kinetic energy and/or momentum threshold
	double KE = e->post[l].Ek();
	double new_pN_mag_p = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
	double costheta = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

	cout << "and costheta " << costheta << '\n';
	cout << "Energy check: Total Enery  " << e->post[l].E() << " and mass " << e->post[l].m() << '\n';

	// Check threshold on energy, momentum and/or angle
//      if ( KE > THRESHOLD_KE_PROTON ) {
	if ( new_pN_mag_p > THRESHOLD_P_PROTON ) {
//      if ( new_pN_mag_p > THRESHOLD_P_PROTON && costheta > THRESHOLD_COSTHETA_PROTON ) {
	  ++nf_p_above_threshold;
          STVs_ok = true; // There is a leading proton with momentum exceeding threshold -> STVs can be calculated
	}

	// Get leading and second leading proton momenta
	if ( pN_mag_p < new_pN_mag_p ) {
	  pN_mag_p2 = pN_mag_p; // first save second to leading proton momentum...
	  index_leading_N2 = index_leading_N; // ...and second to leading proton index
	  pN_mag_p = new_pN_mag_p; // then overwrite leading proton momentum...
	  index_leading_N = l; // ...and its index
	}

	if ( new_pN_mag_p > THRESHOLD_P_PROTON && new_pN_mag_p < pN_mag_p && new_pN_mag_p > pN_mag_p2 ) {
          pN_mag_p2 = new_pN_mag_p; // first save second to leading proton momentum...
          index_leading_N2 = l; // ...and second to leading proton index	  
	}

	cout << "and total momentum " << pN_mag_p << '\n';

      } // End of proton check
 
//      else if ( std::abs(e->post[l].pdg) == MUON ) { // Beginning of muon check (already done above)
//      }


      // Get pion information
      else if ( std::abs(e->post[l].pdg) == PI_PLUS || std::abs(e->post[l].pdg) == PI_ZERO  || std::abs(e->post[l].pdg) == PI_MINUS ) { // Beginning of pion check

	// Save all pion indices
	index_N_multi_pion[ nf_pi_above_threshold ] = l;

	double pNx = e->post[l].p().x;
	double pNy = e->post[l].p().y;
	double pNz = e->post[l].p().z;

	// Check whether this pion is above the kinetic energy or momentum threshold
	double KE = e->post[l].Ek();

	double new_pN_mag_pi = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//	double costheta_pi = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

 
//	if ( KE > THRESHOLD_KE_PION ) {
	if ( new_pN_mag_pi > THRESHOLD_P_PION ) {
//	if ( pN_mag_pi > THRESHOLD_P_PION && costheta_pi > THRESHOLD_COSTHETA_PION ) {
	  ++nf_pi_above_threshold;
	}

	// Get leading pion momentum
	if ( pN_mag_pi < new_pN_mag_pi ) {
	  pN_mag_pi = new_pN_mag_pi; // leading pion momentum
	  index_leading_pi = l;
	}

	cout << "Track belongs to a pion (" << std::abs(e->post[l].pdg) << ")" << endl;
	cout << "E_kin of track  " << l+1 << " is " << KE << endl;
	cout << "nf_pi_above_threshold of event " << k << " is " << nf_pi_above_threshold << endl;

      } // End of pion check


      // Get neutron information
      else if ( std::abs(e->post[l].pdg) == NEUTRON ) { // Beginning of neutron check

	// Save all neutron indices
	index_N_multi_neutron[ nf_neutrons ] = l;

        double pNx = e->post[l].p().x;
        double pNy = e->post[l].p().y;
        double pNz = e->post[l].p().z;

        // Check whether this pion is above the kinetic energy or momentum threshold
        double KE = e->post[l].Ek();
        double new_pN_mag_n = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );
//      double costheta_n = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

        // Check threshold on energy, momentum and/or angle
//      if ( KE > THRESHOLD_KE_NEUTRON ) {
        if ( new_pN_mag_n > THRESHOLD_P_NEUTRON ) {
//      if ( pN_mag_n > THRESHOLD_P_NEUTRON && costheta_n > THRESHOLD_COSTHETA_NEUTRON ) {
	  ++nf_neutrons; // THRESHOLD_P_NEUTRON currently set to 0. GeV
	}

        // Get leading neutron momentum
        if ( pN_mag_n < new_pN_mag_n ) {
        pN_mag_n = new_pN_mag_n; // leading neutron momentum
        index_leading_n = l;
        }

        cout << "Track belongs to a neutron (" << std::abs( e->post[l].pdg ) << ")" << endl;
        cout << "E_kin of track  " << l+1 << " is " << KE << endl;
        cout << "nf_neutrons of event " << k << " is " << nf_neutrons << endl;

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
      E_kin_leading_p = e->post[ index_leading_N ].Ek() / 1000.;
    }


    // Leading proton momentum components
    double pN_lead_x = e->post[ index_leading_N ].p().x / 1000.;
    double pN_lead_y = e->post[ index_leading_N ].p().y / 1000.;
    double pN_lead_z = e->post[ index_leading_N ].p().z / 1000.;

    cout << "pN_lead_x is " << pN_lead_x << endl;
    cout << "pN_lead_y is " << pN_lead_y << endl;
    cout << "pN_lead_z is " << pN_lead_z << endl;

    // Momentum leading proton kinematic variables
    mc_truth_leading_p_mom = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y + pN_lead_z*pN_lead_z );
    mc_truth_leading_p_phi = std::atan2( pN_lead_y, pN_lead_x );
    mc_truth_leading_p_costheta = pN_lead_z / std::sqrt( mc_truth_leading_p_mom );
    mc_truth_leading_p_theta = std::acos( mc_truth_leading_p_costheta );

    // 2nd momentum leading proton momentum components
    double pN_lead_x2, pN_lead_y2, pN_lead_z2;

    if ( nf_p_above_threshold >= 2 ) {
    pN_lead_x2 = e->post[ index_leading_N2 ].p().x / 1000.;
    pN_lead_y2 = e->post[ index_leading_N2 ].p().y / 1000.;
    pN_lead_z2 = e->post[ index_leading_N2 ].p().z / 1000.;


    cout << "pN_lead_x2 is " << pN_lead_x2 << endl;
    cout << "pN_lead_y2 is " << pN_lead_y2 << endl;
    cout << "pN_lead_z2 is " << pN_lead_z2 << endl;

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
    pN_lead_xpi = e->post[ index_leading_pi ].p().x / 1000.;
    pN_lead_ypi = e->post[ index_leading_pi ].p().y / 1000.;
    pN_lead_zpi = e->post[ index_leading_pi ].p().z / 1000.;

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
    pN_lead_xn = e->post[ index_leading_n ].p().x / 1000.;
    pN_lead_yn = e->post[ index_leading_n ].p().y / 1000.;
    pN_lead_zn = e->post[ index_leading_n ].p().z / 1000.;

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
        pTx_np += e->post[ index_N_multi_proton[i] ].p().x / 1000.;
        pTy_np += e->post[ index_N_multi_proton[i] ].p().y / 1000.;
      }
    }

    // Pions
    if ( nf_pi_above_threshold >= 1 ) {
      for ( int i = 0; i < nf_pi_above_threshold ; ++i) {

        // Add all final (above threshold) pion momenta
        pTx_npi += e->post[ index_N_multi_pion[i] ].p().x / 1000.;
        pTy_npi += e->post[ index_N_multi_pion[i] ].p().y / 1000.;
      }
    }

    // Neutrons
    if ( nf_neutrons >= 1 ) {
      for ( int i = 0; i < nf_neutrons ; ++i) {

        // Add all final (above threshold) proton momenta
        pTx_nn += e->post[ index_N_multi_neutron[i] ].p().x / 1000.;
        pTy_nn += e->post[ index_N_multi_neutron[i] ].p().y / 1000.;
      }
    }


    // Calculate Single Transverse Variables

    // Muon momentum components
    double plx = e->post[ 0 ].p().x / 1000.;
    double ply = e->post[ 0 ].p().y / 1000.;
    double plz = e->post[ 0 ].p().z / 1000.;

    // Transverse momentum components
    double pTx = plx + pN_lead_x; 
    double pTy = ply + pN_lead_y;

    // Calculate delta_p_T
    pTmag = std::sqrt( pTx*pTx + pTy*pTy );
//    pTmag = pTmag; // divided by 1000., so that I have delta_p_T in analogy to GENIE and GiBUU in GeV and not in MeV

    // Calculate delta_alpha_T
    double pTl_mag = std::sqrt( plx*plx + ply*ply );
    delta_alpha_T = std::acos( ( -plx*pTx - ply*pTy ) / ( pTl_mag * pTmag ) );  

    // Calculate delta_phi_T
    double pTN_mag = std::sqrt( pN_lead_x * pN_lead_x + pN_lead_y * pN_lead_y );
    delta_phi_T = std::acos( ( -plx*pN_lead_x - ply*pN_lead_y ) / ( pTl_mag * pTN_mag ) );

    // delta pT_2p (Single Transverse Momentum Imbalance with two final state protons)
    double pTx_2p = plx + pN_lead_x + pN_lead_x2;
    double pTy_2p = ply + pN_lead_y + pN_lead_y2;
    if ( nf_p_above_threshold > 1 ) {
    pTmag_2p = std::sqrt( pTx_2p*pTx_2p + pTy_2p*pTy_2p );
//    pTmag_2p = pTmag_2p; // MeV -> GeV
    }
    else pTmag_2p = BOGUS;

    // delta_pT_nN (Single Transverse Momentum Imbalance with N final state hadrons (protons, pions, neutrons))
    double pTx_nN = plx + pTx_np + pTx_npi + pTx_nn;
    double pTy_nN = ply + pTy_np + pTy_npi + pTy_nn;
    pTmag_nN = std::sqrt( pTx_nN*pTx_nN + pTy_nN*pTy_nN );
//    pTmag_nN = pTmag_nN; // MeV -> GeV

    // STK 4-momentum-imbalance delta p (see arXiv:1805.05486) 
    // There are two ways of calculating delta p: 1) when neutrino energy E_nu is known and 2) when it is unknown
    // 1)

    // Incoming neutrino energy
    E_nu = e->in[ 0 ].p().z / 1000.;
    cout << "Incoming Neutrino Energy " << e->in[ 0 ].p().z / 1000. << endl;

    // delta pL
    pLmag = plz + pN_lead_z - E_nu;

    // delta p
    delta_p = std::sqrt( pTmag*pTmag + pLmag*pLmag );

    // 2)

    double R = TARGET_MASS + plz + pN_lead_z - e->post[ 0 ].p().z / 1000. - ( E_kin_leading_p + PROTON_MASS );

    // Estimated mass of the final remnant nucleus (CCQE assumption)
    double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;

    // delta pL without knowing E_nu
    pLmag_noE = 0.5 * R - ( std::pow( mf, 2 ) + std::pow( pTmag, 2 ) ) / ( 2. * R );

    // delta p without knowing E_nu
    delta_p_noE = std::sqrt( std::pow( pLmag_noE, 2 ) + std::pow( pTmag, 2 ) );

    // Flag processes
    if ( e->flag.qel ) cat = 1; // Quasielastic
    else if ( e->flag.mec ) cat = 2; // MEC
    else if ( e->flag.res ) cat = 3; // Resonance production
    else if ( e->flag.dis ) cat = 4; // Deep Inelastic Scattering
    else if ( e->flag.coh ) cat = 5; // Coherent meson production
    else cat = 6; // Other

//    if ( !STVs_ok || nf_mu_above_threshold == 0 || nf_p_above_thresholds == 0 ) {
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
  TNamed temp_named( "nuwro_tune", nuwro_tune_name.c_str() );
  temp_named.Write();

  out_file.Write();
  out_file.Close();
}
