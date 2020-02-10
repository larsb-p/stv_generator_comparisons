#include <iostream>
//#include "neut_files/src/neutclass/neutvtx.h"
//#include "neut_files/src/neutclass/neutvect.h"

// Load libraries needed to access NEUT output .root-file content
R__LOAD_LIBRARY(neut_files/src/neutclass/neutvtx.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutpart.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutfsipart.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutfsivert.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutnucfsistep.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutnucfsivert.so)
R__LOAD_LIBRARY(neut_files/src/neutclass/neutvect.so)

void load_libraries() {
  gSystem->Load("neut_files/src/neutclass/neutvtx.so");
  gSystem->Load("neut_files/src/neutclass/neutpart.so");
  gSystem->Load("neut_files/src/neutclass/neutfsipart.so");
  gSystem->Load("neut_files/src/neutclass/neutfsivert.so");
  gSystem->Load("neut_files/src/neutclass/neutnucfsistep.so");
  gSystem->Load("neut_files/src/neutclass/neutnucfsivert.so");
  gSystem->Load("neut_files/src/neutclass/neutvect.so");
}


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

// Thresholds on particle kinetic energy
constexpr double THRESHOLD_KE_MUON = 0.;
constexpr double THRESHOLD_KE_PROTON = 0.; // 82.; // MeV
constexpr double THRESHOLD_KE_PION = 0.; // 37.; // MeV

// Thresholds on particle momenta
constexpr double THRESHOLD_P_MUON = 0.; // 250.;
//constexpr double THRESHOLD_P_PROTON = 450.; // NEUT uses whatever is set in Card. Here MeV and  c = 1
constexpr double THRESHOLD_P_PROTON = 300.; // Threshold on proton momentum
constexpr double THRESHOLD_P_PION = 0.; // 300.; // Threshold on pion momentum
constexpr double THRESHOLD_P_NEUTRON = 0.000; // GeV

// Thresholds on particle angle deviated from neutrino direction
constexpr double THRESHOLD_COSTHETA_MUON = 0.; // -0.6;
constexpr double THRESHOLD_COSTHETA_PROTON = 0.; // 0.4;


void compute_stv_NEUT( 	const std::string& input_file_name,
			const std::string& neut_tune_name)
{

  // Load libraries
  load_libraries();

  // Initialize
  Int_t i, j;
  Int_t nevents;

  TTree  *tn;

  // Make object with class
  NeutVtx *nvtx;
  NeutVect *nvect;

  TFile f( input_file_name.c_str() );
  tn = (TTree *)(f.Get("neuttree"));
  if ( !tn ) return;

  nvtx = new NeutVtx();
  tn->SetBranchAddress("vertexbranch",&nvtx);

  nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

  // Get number of events
  nevents = tn->GetEntries();
  cout << "Number of events: " << nevents << endl;
 
  // Make output .root-file
  TFile out_file( ( "stv_files/stv_out_" + neut_tune_name + ".root" ).c_str(), "recreate" );
  TTree* out_tree = new TTree( "stv_tree", "NEUT STV tree" );

  // Initiate variables 
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

  // Output branches
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
//  for ( j = 0 ; j < nevents - nevents + 10000 ; j++ ){ // test
  for ( j = 0 ; j < nevents; j++ ){

    cout << endl;
    cout << "Event " << j << "                                                           <- New event" << endl;

    if ( j % 1000 == 0 ) std::cout << "Entry " << j << '\n';

    cout << "---------------------------------------------" << "\n";

    // Get event one by one
    tn->GetEntry(j);

        /***************************************************************/
	// Print event information
        cout << "Event #        :" << nvect->EventNo << "\n";
        cout << "Target A       :" << nvect->TargetA << "\n";
        cout << "Target Z       :" << nvect->TargetA << "\n";
        cout << "VNuclIni       :" << nvect->VNuclIni << "\n";
        cout << "VNuclFin       :" << nvect->VNuclFin << "\n";
        cout << "PF Surface     :" << nvect->PFSurf   << "\n";
        cout << "PF Maximum     :" << nvect->PFMax    << "\n";
        cout << "Flux ID        :" << nvect->FluxID   << "\n";

        cout << "Intr. mode     :" << nvect->Mode   << "\n";

    // Find indices of leading final state protons, pions and neutrons
    int index_leading_N = 0; // Momentum leading proton 
    int index_leading_N2 = BOGUS; // Momentum 2nd leading proton
    int index_muon = 0; // Muon index
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

    cout << "Tracks in event " << j << " is " << nvect->Npart() << endl;

    // Loop over particles in one event
    for ( i = 0 ; i < nvect->Npart() ; i++ ){ // Beginning of loop over particles in one event

      cout << "--------------------------" << nvect->PartInfo(i)->fPID << "--------------------------" << endl;

      // For FSIs
      if ((nvect->PartInfo(i))->fStatus != 0 || (nvect->PartInfo(i))->fIsAlive != 1 ) continue;

      // Print particle information
      cout << "i=" << i << "\n";
      cout << "Vertex         =" << nvect->VertexID(i) << "\n";
      cout << "Parent Index   =" << nvect->ParentIdx(i) << "\n";

      cout << "Particle Code  = " << (nvect->PartInfo(i))->fPID   << "\n";
      cout << "Particle Mass  = " << (nvect->PartInfo(i))->fMass / 1000.   << "\n";
      cout << "Particle Mom.  =(" << (nvect->PartInfo(i))->fP.Px() / 1000. << ","
       		   << (nvect->PartInfo(i))->fP.Py() / 1000. << ","
                   << (nvect->PartInfo(i))->fP.Pz() / 1000. << ","
                   << (nvect->PartInfo(i))->fP.E() / 1000.  << ")"
                   << "\n";
      cout << "Particle Flag  = " << (nvect->PartInfo(i))->fIsAlive << "\n";
      cout << "Particle Stat. = " << (nvect->PartInfo(i))->fStatus  << "\n";
      cout << "Particle Pos(1)=(" << (nvect->PartInfo(i))->fPosIni.X() << ","
                   << (nvect->PartInfo(i))->fPosIni.Y() << ","
                   << (nvect->PartInfo(i))->fPosIni.Z() << ","
                   << (nvect->PartInfo(i))->fPosIni.T()  << ")"
                   << "\n";

      // Get proton information
      if ( ( nvect->PartInfo(i) )->fPID == PROTON ) { // Beginning of proton check

	// Save all proton indices
	index_N_multi_proton[ nf_p_above_threshold ] = i;
 
	// Get proton momentum components
	double pNx = ( nvect->PartInfo(i) )->fP.Px();
	double pNy = ( nvect->PartInfo(i) )->fP.Py();
	double pNz = ( nvect->PartInfo(i) )->fP.Pz();


	// Calculate energy, magnitude of proton momentum and angle (for (momentum) leading and second leading proton )
	double KE = ( nvect->PartInfo(i) )->fP.E() - PROTON_MASS; // Get kinetic energy of proton
	double new_pN_mag_p = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz ); // Calculate momentum of proton
	double costheta = pNz / std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz ); // Calculate cos(theta)

	// Is proton above thresholds (set above)?
//	if ( KE > THRESHOLD_KE_PROTON ) {
	if ( new_pN_mag_p > THRESHOLD_P_PROTON ) {
//	if ( new_pN_mag_p > THRESHOLD_P_PROTON && costheta > THRESHOLD_COSTHETA_PROTON ) {
	
	  // Count proton above thresholds
	  ++nf_p_above_threshold;
	  STVs_ok = true; // There is a leading proton with momentum exceeding threshold -> STVs can be calculated
	}

	// Find leading and second leading proton in event
	if ( pN_mag_p < new_pN_mag_p ) {

	  // Save prior momentum (will be second leading proton momentum)
          pN_mag_p2 = pN_mag_p; // first save second to leading proton momentum...
          index_leading_N2 = index_leading_N; // ...and second to leading proton index
	  pN_mag_p = new_pN_mag_p; // then overwrite leading proton momentum...
	  index_leading_N = i; // ...and its index
	}

        if ( new_pN_mag_p > THRESHOLD_P_PROTON && new_pN_mag_p < pN_mag_p && new_pN_mag_p > pN_mag_p2 ) {
          pN_mag_p2 = new_pN_mag_p; // first save second to leading proton momentum...
          index_leading_N2 = i; // ...and second to leading proton index          
        }

        cout << "Track belongs to a proton (" << ( nvect->PartInfo(i) )->fPID << ")" << endl;
        cout << "E_kin of track  " << i+1 << " is " << KE / 1000. << endl;
        cout << "nf_p_above_threshold of event " << j << " is " << nf_p_above_threshold << endl;

      } // End of proton


      // Find muon
      else if ( (  nvect->PartInfo(i) )->fPID == MUON ) { // Beginning of muon check

	// Get muon momentum components
	double plx = ( nvect->PartInfo(i) )->fP.Px();
	double ply = ( nvect->PartInfo(i) )->fP.Py();
	double plz = ( nvect->PartInfo(i) )->fP.Pz();

        // Get kinetic energy of muon
        double KE_mu = ( nvect->PartInfo(i) )->fP.E() - MUON_MASS;

	// Calculate magnitude of muon momentum
	double pN_mag_p_muon = std::sqrt( plx*plx + ply*ply + plz*plz );

	// Calculate costheta
	double costheta_mu = plz / std::sqrt( pN_mag_p_muon );

	// Is muon above thresholds (set above)?
//	if ( KE_mu > THRESHOLD_KE_MUON ) {
	if ( pN_mag_p_muon > THRESHOLD_P_MUON ) {
//	if ( pN_mag_p_muon > THRESHOLD_P_MUON && costheta_mu > THRESHOLD_COSTHETA_MUON ) {

	  // Count muon above thresholds
	  ++nf_mu_above_threshold;
	}

	// Save index of muon in event
	index_muon = i;

      cout << "E_kin of muon track is " << KE_mu / 1000. << endl;
      cout << "Number of muons above threshold of event " << j << " is " << nf_mu_above_threshold << endl;

      } // End of muon check


      // Find pions
      else if ( (  nvect->PartInfo(i) )->fPID == PI_PLUS || ( nvect->PartInfo(i) )->fPID == PI_ZERO || ( nvect->PartInfo(i) )->fPID == PI_MINUS ) { // Beginning of pion check

        // Save all pion indices
        index_N_multi_pion[ nf_pi_above_threshold ] = i;

	// Get pion momentum components
	double pNx = ( nvect->PartInfo(i) )->fP.Px();
	double pNy = ( nvect->PartInfo(i) )->fP.Py();
	double pNz = ( nvect->PartInfo(i) )->fP.Pz();

	double KE;

        // Get kinetic energy of pion
        if ( (  nvect->PartInfo(i) )->fPID == PI_PLUS ) KE = ( nvect->PartInfo(i) )->fP.E() - PI_PLUS_MASS;
        if ( (  nvect->PartInfo(i) )->fPID == PI_ZERO ) KE = ( nvect->PartInfo(i) )->fP.E() - PI_ZERO_MASS;
        if ( (  nvect->PartInfo(i) )->fPID == PI_MINUS ) KE = ( nvect->PartInfo(i) )->fP.E() - PI_MINUS_MASS;

	// Calculate magnitude of pion momentum
	double new_pN_mag_pi = std::sqrt( pNx*pNx + pNy*pNy + pNz*pNz );

	// Calculate costheta
	double costheta = pNz / std::sqrt( new_pN_mag_pi );  

	// Is pion above thresholds (set above)?
//	if ( KE > THRESHOLD_KE_PION ) {
	if ( new_pN_mag_pi > THRESHOLD_P_PION ) {
//	if ( pN_mag_pi > THRESHOLD_P_PION && costheta_pi > THRESHOLD_COSTHETA_PION ) {

	  // Count proton above thresholds
	  ++nf_pi_above_threshold; // THRESHOLD_P_PION currently set to 0. GeV
        }

	// Get leading pion momentum
	if ( pN_mag_pi < new_pN_mag_pi ) {
	  pN_mag_pi = new_pN_mag_pi; // leading pion momentum
	  index_leading_pi = i;
	}

        cout << "Track belongs to a pion (" << ( nvect->PartInfo(i) )->fPID << ")" << endl;
        cout << "E_kin of track  " << i+1 << " is " << KE / 1000. << endl;
        cout << "nf_pi_above_threshold of event " << j << " is " << nf_pi_above_threshold << endl;

      } // End of pion check


      // Get neutron information
      else if ( ( nvect->PartInfo(i) )->fPID == NEUTRON ) { // Beginning of neutron check

        // Save all neutron indices
        index_N_multi_neutron[ nf_neutrons ] = i;

        // Get neutron momentum components
        double pNx = ( nvect->PartInfo(i) )->fP.Px();
        double pNy = ( nvect->PartInfo(i) )->fP.Py();
        double pNz = ( nvect->PartInfo(i) )->fP.Pz();

        // Check whether this pion is above the kinetic energy or momentum threshold
        double KE = ( nvect->PartInfo(i) )->fP.E() - NEUTRON_MASS;
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
        index_leading_n = i;
        }

        cout << "Track belongs to a neutron (" << ( nvect->PartInfo(i) )->fPID << ")" << endl;
        cout << "E_kin of track  " << i+1 << " is " << KE / 1000. << endl;
        cout << "nf_neutrons of event " << j << " is " << nf_neutrons << endl;

      } // End of neutron check

    } // End of particle loop

    std::cout << "Event #" << nvtx->EventNo << "\n";
    cout << endl;
    cout << "STVs ok? " << STVs_ok << endl;
    cout << "Number of protons above threshold of event " << j << " is " << nf_p_above_threshold << endl;
    cout << "Number of pions above threshold of event " << j << " is " << nf_pi_above_threshold << endl;
    cout << "Number of neutrons above threshold of event " << j << " is " << nf_neutrons << endl;
    cout << endl;

    // Momentum leading proton kinetic energy
    if ( nf_p_above_threshold > 0 ) {
      E_kin_leading_p = ( ( nvect->PartInfo( index_leading_N ) )->fP.E() - PROTON_MASS ) / 1000.;
    }

    // Print vertex position of event
    for (i = 0 ; i < nvtx->Nvtx() ; i++){
      cout << "i=" << i << "\n";

      cout << "Vertex Pos(1)=(" << (nvtx->Pos(i))->X() << ","
                   << (nvtx->Pos(i))->Y() << ","
                   << (nvtx->Pos(i))->Z() << ","
                   << (nvtx->Pos(i))->T()  << ")"
                   << "\n";
    } // End of vertex position loop


    // Get momentum components of this event's leading proton
    double pN_lead_x = ( nvect->PartInfo( index_leading_N ) )->fP.Px() / 1000.;
    double pN_lead_y = ( nvect->PartInfo( index_leading_N ) )->fP.Py() / 1000.;
    double pN_lead_z = ( nvect->PartInfo( index_leading_N ) )->fP.Pz() / 1000.;

    // Save leading proton kinematic variables
    mc_truth_leading_p_mom = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y + pN_lead_z*pN_lead_z ); // Calculate momentum for leading proton
    mc_truth_leading_p_phi = std::atan2( pN_lead_y, pN_lead_x ); // Calculate angle phi for leading proton
    mc_truth_leading_p_costheta = pN_lead_z / mc_truth_leading_p_mom; // Calculate cos(theta) for leading proton
    mc_truth_leading_p_theta = std::acos( mc_truth_leading_p_costheta ); // Calculate angle theta for leading proton

    cout << "Index of leading proton: " << index_leading_N << endl;
    cout << "Leading proton momentum: " << mc_truth_leading_p_mom << endl;


    // 2nd momentum leading proton momentum components
    double pN_lead_x2, pN_lead_y2, pN_lead_z2;

    if ( nf_p_above_threshold >= 2 ) {
    pN_lead_x2 = ( nvect->PartInfo( index_leading_N2 ) )->fP.Px() / 1000.;
    pN_lead_y2 = ( nvect->PartInfo( index_leading_N2 ) )->fP.Py() / 1000.;
    pN_lead_z2 = ( nvect->PartInfo( index_leading_N2 ) )->fP.Pz() / 1000.;

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


    // Get components of this event's muon
    double p_mu_x = ( nvect->PartInfo( index_muon ) )->fP.Px() / 1000.;
    double p_mu_y = ( nvect->PartInfo( index_muon ) )->fP.Py() / 1000.;
    double p_mu_z = ( nvect->PartInfo( index_muon ) )->fP.Pz() / 1000.;

    // Save kinematic variables of MC muon
    mc_truth_mu_mom = std::sqrt( p_mu_x*p_mu_x + p_mu_y*p_mu_y + p_mu_z*p_mu_z ); // Calculate muon momentum
    mc_truth_mu_costheta = p_mu_z / std::sqrt( mc_truth_mu_mom ); // Calculate muon cos(theta)
    mc_truth_mu_theta = std::acos( mc_truth_mu_costheta ); // Calculate muon deviation angle theta
    mc_truth_mu_phi = std::atan2( p_mu_y, p_mu_x ); // Calculate muon deviation angle phi

    cout << "Index of muon: " << index_muon << endl;
    cout << "Muon momentum: " << mc_truth_mu_mom << endl;


    // Pion momentum components
    double pN_lead_xpi, pN_lead_ypi, pN_lead_zpi;

    if ( nf_pi_above_threshold > 0 ) {
    pN_lead_xpi = ( nvect->PartInfo( index_leading_pi ) )->fP.Px() / 1000.;
    pN_lead_ypi = ( nvect->PartInfo( index_leading_pi ) )->fP.Py() / 1000.;
    pN_lead_zpi = ( nvect->PartInfo( index_leading_pi ) )->fP.Pz() / 1000.;

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
    pN_lead_xn = ( nvect->PartInfo( index_leading_n ) )->fP.Px() / 1000.;
    pN_lead_yn = ( nvect->PartInfo( index_leading_n ) )->fP.Py() / 1000.;
    pN_lead_zn = ( nvect->PartInfo( index_leading_n ) )->fP.Pz() / 1000.;

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
        pTx_np += ( nvect->PartInfo( index_N_multi_proton[i] ) )->fP.Px() / 1000.;
        pTy_np += ( nvect->PartInfo( index_N_multi_proton[i] ) )->fP.Py() / 1000.;
      }
    }

    // Pions
    if ( nf_pi_above_threshold >= 1 ) {
      for ( int j = 0; j < nf_pi_above_threshold ; ++j) {

        // Add all final (above threshold) pion momenta
        pTx_npi += ( nvect->PartInfo( index_N_multi_pion[j] ) )->fP.Px() / 1000.;
        pTy_npi += ( nvect->PartInfo( index_N_multi_pion[j] ) )->fP.Py() / 1000.;
      }
    }

    // Neutrons
    if ( nf_neutrons >= 1 ) {
      for ( int j = 0; j < nf_neutrons ; ++j) {

        // Add all final (above threshold) neutrons momenta
        pTx_nn += ( nvect->PartInfo( index_N_multi_neutron[j] ) )->fP.Px() / 1000.;
        pTy_nn += ( nvect->PartInfo( index_N_multi_neutron[j] ) )->fP.Py() / 1000.;
      }
    }

    // Calculate single transverse momentum components of leading proton
    double pTx = p_mu_x + pN_lead_x;
    double pTy = p_mu_y + pN_lead_y;

    // Calculate single transverse momentum magnitude
    pTmag = std::sqrt( pTx*pTx + pTy*pTy );

    // Calculate magnitude of transverse lepton momentum
    double pTl_mag = std::sqrt( p_mu_x*p_mu_x + p_mu_y*p_mu_y );

    // Calculate single transverse angle delta_alpha_T
    delta_alpha_T = std::acos( ( -p_mu_x*pTx - p_mu_y*pTy ) / ( pTl_mag * pTmag ) );

    // Calculate magnitude of leading proton momentum
    double pTN_mag = std::sqrt( pN_lead_x*pN_lead_x + pN_lead_y*pN_lead_y );	

    // Calculate single transverse angle delta_T
    delta_phi_T = std::acos( ( -p_mu_x*pN_lead_x - p_mu_y*pN_lead_y ) / ( pTl_mag * pTN_mag ) );	

    // delta pT_2p (Single Transverse Momentum Imbalance with two final state protons)
    double pTx_2p = p_mu_x + pN_lead_x + pN_lead_x2;
    double pTy_2p = p_mu_y + pN_lead_y + pN_lead_y2;
    if ( nf_p_above_threshold > 1 ) {
    pTmag_2p = std::sqrt( pTx_2p*pTx_2p + pTy_2p*pTy_2p );
    }
    else pTmag_2p = BOGUS;

    // delta_pT_nN (Single Transverse Momentum Imbalance with N final state hadrons (protons, pions, neutrons))
    double pTx_nN = p_mu_x + pTx_np + pTx_npi + pTx_nn;
    double pTy_nN = p_mu_y + pTy_np + pTy_npi + pTy_nn;
    pTmag_nN = std::sqrt( pTx_nN*pTx_nN + pTy_nN*pTy_nN );


    // STK 4-momentum-imbalance delta p (see arXiv:1805.05486) 
    // There are two ways of calculating delta p: 1) when neutrino energy E_nu is known and 2) when it is unknown
    // 1)

    // Incoming neutrino energy
    E_nu = (nvect->PartInfo(0))->fP.E() / 1000.;
    cout << "Incoming Neutrino Energy " << E_nu << endl;

    // delta pL
    pLmag = p_mu_z + pN_lead_z - E_nu;

    // delta p
    delta_p = std::sqrt( pTmag*pTmag + pLmag*pLmag );

    // 2)

    double R = TARGET_MASS + p_mu_z + pN_lead_z - ( nvect->PartInfo( index_muon ) )->fP.E() / 1000. - (nvect->PartInfo( index_leading_N ) )->fP.E() / 1000.;

    // Estimated mass of the final remnant nucleus (CCQE assumption)
    double mf = TARGET_MASS - NEUTRON_MASS + BINDING_ENERGY;

    // delta pL without knowing E_nu
    pLmag_noE = 0.5 * R - ( std::pow( mf, 2 ) + std::pow( pTmag, 2 ) ) / ( 2. * R );

    // delta p without knowing E_nu
    delta_p_noE = std::sqrt( std::pow( pLmag_noE, 2 ) + std::pow( pTmag, 2 ) );


    // Flag event according to interaction mode (channel)
    if ( index_muon != 0 ) cc = true; // CC

    // Flag processes
    if ( nvect->Mode == 1 ) cat = 1; // Quasielastic
    else if ( nvect->Mode == 2 ) cat = 2; // 2p2h (MEC included)
    else if ( nvect->Mode > 10 && nvect->Mode < 30 && nvect->Mode != 16 && nvect->Mode != 26 ) cat = 3; // 1p (Resonance production)
    else if ( nvect->Mode == 26 ) cat = 4; // Deep Inelastic Scattering
    else if ( nvect->Mode == 16 ) cat = 5; // Coherent meson production
    else cat = 6; // Other

    // Only fill tree with calculated values if it is a 1muNp event that satisfies set thresholds
    if ( !STVs_ok || !cc ) {
      mc_truth_mu_mom = BOGUS;
      mc_truth_mu_phi = BOGUS;
      mc_truth_mu_theta = BOGUS;
      mc_truth_mu_costheta = BOGUS;
      mc_truth_leading_p_mom = BOGUS;
      mc_truth_leading_p_phi = BOGUS;
      mc_truth_leading_p_costheta = BOGUS;
      mc_truth_leading_p_theta = BOGUS;
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

  TNamed temp_named( "neut_tune", neut_tune_name.c_str() );
  temp_named.Write();

  out_file.Write();
  out_file.Close();
} // End of void loop

