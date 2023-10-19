#include <iostream>
#include <vector>
#include <string>
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <ROOT/RDataFrame.hxx>
#include <string>
#include <utility>

std::vector<std::string> particle_types = {"protonp", "kaonp", "pionp", "protonm", "kaonm", "pionm"};

void create_lorentz_vectors(ROOT::RDataFrame df, const std::vector<std::string>& particle_types) {
    for (const auto& particle : particle_types) {
        *df = df->Define(particle, "TLorentzVector(" + particle + "_PX, " + particle + "_PY, " + particle + "_PZ, " + particle + "_PE)");
    }
}

void all_combos() {
    auto d = ROOT::RDataFrame("B2LcLcpiOS/DecayTree", {"/data/lhcb01/mwhitehead/LcLcpi_2018_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2018_MD.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MU.root","/data/lhcb01/mwhitehead/LcLcpi_2017_MD.root"});
    auto df = d.Filter("Lambdacp_p_ProbNNp > 0.7 && Lambdacp_K_ProbNNk > 0.6 && Lambdacp_pi_ProbNNpi > 0.4 && Lambdacm_p_ProbNNp > 0.7 && Lambdacm_K_ProbNNk > 0.6 && Lambdacm_pi_ProbNNpi > 0.4");
    create_lorentz_vectors(&df, particle_types);
}

