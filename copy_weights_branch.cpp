void copy_weights_branch( const char* destFilePath, const char* sourceFilePath, const char* systParam ) {
    // Open the source file (read mode)
    TFile *sourceFile = TFile::Open(sourceFilePath, "READ");
    if (!sourceFile || sourceFile->IsZombie()) {
        std::cerr << "Error opening source file!" << std::endl;
        return;
    }

    // Open the destination file (update mode)
    TFile *destFile = TFile::Open(destFilePath, "UPDATE");
    if (!destFile || destFile->IsZombie()) {
        std::cerr << "Error opening destination file!" << std::endl;
        return;
    }

    // Get the tree from the source file
    TTree *sourceTree = (TTree*)sourceFile->Get(systParam);
    if (!sourceTree) {
        std::cerr << "Error: Tree " << systParam << " not found in source file!" << std::endl;
        return;
    }

    // Create a new tree in the destination file, copying structure but no entries
    TTree *newTree = sourceTree->CloneTree(0);  // 0 means no entries copied

    // Get the address of the weights branch
    TArrayF *weights = nullptr;
    sourceTree->SetBranchAddress("weights", &weights);

    // Copy the entries from the source tree to the new tree in the destination file
    for (int i = 0; i < sourceTree->GetEntries(); ++i) {
        sourceTree->GetEntry(i);
        newTree->Fill();  // Fill the new tree with data
    }

    // Write the new tree (overwriting if already exists) into the destination file
    newTree->Write(systParam, TObject::kOverwrite);

    // Close both files
    sourceFile->Close();
    destFile->Close();
}

