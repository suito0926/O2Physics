// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
///
/// \brief This task is an empty skeleton that fills a simple eta histogram.
///        it is meant to be a blank page for further developments.
/// \author everyone

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/ASoAHelpers.h"
#include "EMCALBase/Geometry.h"
#include "EMCALCalib/BadChannelMap.h"
#include "PWGJE/DataModel/EMCALClusters.h"
#include "DataFormatsEMCAL/Cell.h"
#include "DataFormatsEMCAL/Constants.h"
#include "DataFormatsEMCAL/AnalysisCluster.h"
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA, 
o2::aod::pidTOFPi, o2::aod::pidTOFEl, o2::aod::pidTPCPi, o2::aod::pidTPCEl, o2::aod::pidTOFmass, o2::aod::pidTOFbeta>>;

struct myExampleTask {

  Preslice<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalclustercell::emcalclusterId;
  
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> nBinsPt{"nBinsPt", 5000, "N bins in pT histo"};
  Configurable<int> nBinsNSigma{"nBinsNSigma", 120, "N bins in NSigma histo"};
  Configurable<int> nBinsNSigma_1{"NBinsNSigma_1", 5000, "N bins NSigma histo"};

  Configurable<float> dcaxy_cut{"dcaxy_cut", 2.0f, ""};
  Filter dcaxyfilter = (nabs(aod::track::dcaXY) < dcaxy_cut);

  void init(InitContext const&)
  {
    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisEta{50, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{nBinsPt, 0, 10, "p_{T}(GeV/c)"};
    const AxisSpec axisPt_2{nBinsPt, 0, 100, "p_{T}(GeV/c)"};
    const AxisSpec axisP{nBinsPt, 0, 10, "p(GeV/c)"};
    const AxisSpec axisNSigmaTPCEl{nBinsNSigma, -6.0, 6.0, "N_{#sigma}^{TPC}(e)"};
    const AxisSpec axisNSigmaTOFEl{nBinsNSigma, -6.0, 6.0, "N_{#sigma}^{TOF}(e)"};
    const AxisSpec axisNSigmaTPCEl_1{nBinsNSigma_1, -10.0, 10.0, "N_{#sigma}^{TPC}(e)"};
    const AxisSpec axisNSigmaTOFEl_1{nBinsNSigma_1, -6.0, 6.0, "N_{#sigma}^{TOF}(e)"};
    const AxisSpec axisdca{500, -6.0, 6.0, "cm"};
    const AxisSpec axisChi2{500, 0.0, 50.0, "#chi^{2}"};
    const AxisSpec axisCluster{100, 0.0, 200.0, "counts"};
    const AxisSpec axisITSNCls{20, 0.0, 20, "counts"};
    const AxisSpec axisEP{30, 0, 3.0, "E/p"};
    const AxisSpec axisE{500, 0, 5.0, "GeV"};
    const AxisSpec axisM02{500, 0.0, 5.0, "M02"};
    const AxisSpec axisdPhi{200, -0.5, 0.5, "dPhi"};
    const AxisSpec axisdEta{200, -0.5, 0.5, "dEta"};

    histos.add("hEventCounter", "hEventCounter", kTH1F, {axisCounter});
    histos.add("etaHistogram", "etaHistogram", kTH1F, {axisEta});
    histos.add("PtHistogram", "PtHistogram", kTH1F, {axisPt});
    histos.add("PHist", "PHist", kTH1F, {axisP});
    histos.add("N_sigma_TPC_eHist", "N_{#sigma}^{TPC}(e)Hist", kTH2F, {{axisP}, {axisNSigmaTPCEl}});
    histos.add("N_sigma_TOF_eHist", "N_{#sigma}^{TOF}(e)Hist", kTH2F, {{axisP}, {axisNSigmaTOFEl}});
    histos.add("NSigmaTPCElHist/wocut", "N_{#sigma}^{TPC}(e) w/oCut", kTH1F, {axisNSigmaTPCEl_1});
    histos.add("NSigmaTPCElHist/wcut", "N_{#sigma}^{TPC}(e) wCut", kTH1F, {axisNSigmaTPCEl_1});
    histos.add("NSigmaTPCElHist/wpcut", "N_{#sigma}^{TPC}(e) wPCut", kTH1F, {axisNSigmaTPCEl_1});
    histos.add("only_N_Sigma_TOF_eHist", "N_{#sigma}^{TOF}(e)", kTH1F, {axisNSigmaTOFEl_1});
    histos.add("DCA_XY_Hist", "DCAxyHist", kTH1F, {axisdca});
    histos.add("DCA_Z_Hist", "DCAzHist", kTH1F, {axisdca});
    histos.add("ITS_Chi2_Hist", "ITS #chi^{2} Hist", kTH1F, {axisChi2});
    histos.add("TPC_Chi2_Hist", "TPC #chi^{2} Hist", kTH1F, {axisChi2});
    histos.add("TPC_NCls_Hist", "TPC_NCls_Hist", kTH1F, {axisCluster});
    histos.add("ITS_NCls_Hist", "ITS_NCls_Hist", kTH1F, {axisITSNCls});
    histos.add("TPC_NClsCrossedRows_Hist", "TPC_NClsCrossedRows_Hist", kTH1F, {axisCluster});
    histos.add("E_hist", "E_hist", kTH1F, {axisE});
    histos.add("M02hist", "M02hist", kTH1F, {axisM02});

    histos.add("E_PwoGM/onlyE_Pwocut", "E/P w/ocut", kTH1F, {axisEP});
    histos.add("E_PwoGM/E_PvPtwocut", "E/p Hist w/ocut", kTH2F, {{axisPt_2}, {axisEP}});
    histos.add("E_PwoGM/onlyE_Pwcut", "E/P w/cut", kTH1F, {axisEP});
    histos.add("E_PwoGM/E_PvPtwcut", "E/p Hist w/cut", kTH2F, {{axisPt_2}, {axisEP}});
    histos.add("E_PwoGM/onlyE_Pwpcut", "E/P w/pcut", kTH1F, {axisEP});
    histos.add("E_PwoGM/E_PvPtwpcut", "E/p Hist w/pcut", kTH2F, {{axisPt_2}, {axisEP}});

    histos.add("E_P/onlyE_Pwocut", "E/P w/ocut", kTH1F, {axisEP});
    histos.add("E_P/E_PvPtwocut", "E/p Hist w/ocut", kTH2F, {{axisPt_2}, {axisEP}});
    histos.add("E_P/onlyE_Pwcut", "E/P w/cut", kTH1F, {axisEP});
    histos.add("E_P/E_PvPtwcut", "E/p Hist w/cut", kTH2F, {{axisPt_2}, {axisEP}});
    histos.add("E_P/onlyE_Pwpcut", "E/P w/pcut", kTH1F, {axisEP});
    histos.add("E_P/E_PvPtwpcut", "E/p Hist w/pcut", kTH2F, {{axisPt_2}, {axisEP}});

    histos.add("TrMatchHistogram", "TrMatchHistogram", kTH2F, {{axisdPhi}, {axisdEta}});
    histos.add("TrMatchHistogram_wcut", "TrMatchHistogram_wcut", kTH2F, {{axisdPhi}, {axisdEta}});
  }

  void process(aod::Collision const& collision, TrackEle const& tracks, o2::aod::EMCALClusters const& clusters,
                o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      histos.fill(HIST("etaHistogram"), track.eta());
      histos.fill(HIST("PtHistogram"), track.pt());
      histos.fill(HIST("PHist"), track.p());
      histos.fill(HIST("DCA_XY_Hist"), track.dcaXY());
      histos.fill(HIST("DCA_Z_Hist"), track.dcaZ());
      histos.fill(HIST("ITS_Chi2_Hist"), track.itsChi2NCl());
      histos.fill(HIST("TPC_Chi2_Hist"), track.tpcChi2NCl());
      histos.fill(HIST("TPC_NCls_Hist"), track.tpcNClsFound());
      histos.fill(HIST("ITS_NCls_Hist"), track.itsNCls());
      histos.fill(HIST("TPC_NClsCrossedRows_Hist"), track.tpcNClsCrossedRows());
      histos.fill(HIST("NSigmaTPCElHist/wocut"), track.tpcNSigmaEl());
      if (track.tpcNClsCrossedRows() < 100) continue;
      if (fabs(track.dcaXY()) > 2.4) continue;
      if (fabs(track.dcaZ()) > 2.0) continue;
      if (fabs(track.eta()) > 0.6) continue;
      if (track.itsChi2NCl() > 15) continue;
      if (track.tpcChi2NCl() > 4) continue;
      if (track.tpcNClsFound() < 100) continue;
      if (track.itsNCls() < 2) continue;
      histos.fill(HIST("NSigmaTPCElHist/wcut"), track.tpcNSigmaEl());
      if (track.p() > 2.0) continue;
      if (track.p() < 0.5) continue;
      histos.fill(HIST("N_sigma_TPC_eHist"), track.p(), track.tpcNSigmaEl());
      histos.fill(HIST("N_sigma_TOF_eHist"), track.p(), track.tofNSigmaEl());
      histos.fill(HIST("NSigmaTPCElHist/wpcut"), track.tpcNSigmaEl());
      histos.fill(HIST("only_N_Sigma_TOF_eHist"), track.tofNSigmaEl());
    }

    for (const auto& cluster : clusters) {
      if (cluster.energy() < 0.5) continue; // 低エネルギークラスターを除外

      histos.fill(HIST("E_hist"), cluster.energy());
      histos.fill(HIST("M02hist"), cluster.m02());

      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, cluster.globalIndex());

      if (tracksofcluster.size() > 0) {


        for (const auto& match : tracksofcluster) {
          double eop = cluster.energy()/match.track_as<TrackEle>().p();
          double dPhi = match.track_as<TrackEle>().phi() - cluster.phi();
          double dEta = match.track_as<TrackEle>().eta() - cluster.eta();
          double deltaR = sqrt(dPhi * dPhi + dEta * dEta);

          histos.fill(HIST("TrMatchHistogram"), dPhi, dEta);

          histos.fill(HIST("E_PwoGM/onlyE_Pwocut"), eop);
          histos.fill(HIST("E_PwoGM/E_PvPtwocut"), match.track_as<TrackEle>().pt(), eop);

          if(match.track_as<TrackEle>().tpcNClsCrossedRows()>100 && 
          match.track_as<TrackEle>().itsNCls()>2 && 
          TMath::Abs(match.track_as<TrackEle>().dcaXY())<2.4 && 
          match.track_as<TrackEle>().tpcNSigmaEl()>-1.0 && 
          match.track_as<TrackEle>().tpcNSigmaEl()<3.0)
          {
            histos.fill(HIST("E_PwoGM/onlyE_Pwcut"), eop);
            histos.fill(HIST("E_PwoGM/E_PvPtwcut"), match.track_as<TrackEle>().pt(), eop);
          }

          if(match.track_as<TrackEle>().tpcNClsCrossedRows()>100 && 
          match.track_as<TrackEle>().itsNCls()>2 && 
          TMath::Abs(match.track_as<TrackEle>().dcaXY())<2.4 && 
          match.track_as<TrackEle>().tpcNSigmaEl()>-1.0 && 
          match.track_as<TrackEle>().tpcNSigmaEl()<3.0 &&
          match.track_as<TrackEle>().pt()>1.0)
          {
            histos.fill(HIST("E_PwoGM/onlyE_Pwpcut"), eop);
            histos.fill(HIST("E_PwoGM/E_PvPtwpcut"), match.track_as<TrackEle>().pt(), eop);
          }

          // マッチング条件を厳格化
          if (fabs(dPhi) < 0.1 && fabs(dEta) < 0.1){
            histos.fill(HIST("E_P/onlyE_Pwocut"), eop);
            histos.fill(HIST("E_P/E_PvPtwocut"), match.track_as<TrackEle>().pt(), eop);

            histos.fill(HIST("TrMatchHistogram_wcut"), dPhi, dEta);
            
            if(match.track_as<TrackEle>().tpcNClsCrossedRows()<100)continue;
            if(match.track_as<TrackEle>().itsNCls()<2)continue;
            if(TMath::Abs(match.track_as<TrackEle>().dcaXY())>2.4)continue;
            if(cluster.m02()<0.1 || cluster.m02()>0.3) continue;
            if(match.track_as<TrackEle>().tpcNSigmaEl()<-1.0 || match.track_as<TrackEle>().tpcNSigmaEl()>3.0) continue;
            histos.fill(HIST("E_P/onlyE_Pwcut"), eop);
            histos.fill(HIST("E_P/E_PvPtwcut"), match.track_as<TrackEle>().pt(), eop);

            if(match.track_as<TrackEle>().pt()<1.0)continue;
            
            histos.fill(HIST("E_P/onlyE_Pwpcut"), eop);
            histos.fill(HIST("E_P/E_PvPtwpcut"), match.track_as<TrackEle>().pt(), eop);
          };

        }

      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<myExampleTask>(cfgc)};
}