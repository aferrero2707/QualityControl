// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

///
/// \file    MergeMCH.cxx
/// \author  Sebastien Perrin
///
/// \brief DPL Wotkflow for MCH Merger

#include "Framework/RootSerializationSupport.h"
#include "Mergers/MergerBuilder.h"

#include <Framework/CompletionPolicy.h>

using namespace o2::framework;
using namespace o2::mergers;

void customize(std::vector<CompletionPolicy>& policies)
{
  MergerBuilder::customizeInfrastructure(policies);
}

#include "Framework/runDataProcessing.h"
#include "Mergers/MergerInfrastructureBuilder.h"
#include "MCH/CustomMergeableTH2Quotient.h"

#include <fairmq/FairMQLogger.h>
#include <TH1F.h>
#include <memory>
#include <random>

using namespace std::chrono;
using namespace o2::quality_control_modules::muonchambers;

// clang-format off
WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  WorkflowSpec specs;

  // one 1D histo
  {
//    WorkflowSpec specs; // enable comment to disable the workflow

    size_t producersAmount = 8;
    Inputs mergersInputs;
    for (size_t p = 0; p < producersAmount; p++) {
      mergersInputs.push_back({ "mo",               "MCH",
                                "preclusters",            static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1),
                                Lifetime::Timeframe });
      DataProcessorSpec producer{
        "producer-histo" + std::to_string(p), Inputs{},
        Outputs{ { { "mo" },
                   "MCH",
                   "preclusters",
                   static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1),
                   Lifetime::Timeframe } },
        AlgorithmSpec{(AlgorithmSpec::ProcessCallback)
                      [ p, producersAmount ](ProcessingContext & processingContext) mutable {

                        //            usleep(100000 + (rand() % 10000) - 5000);
                        usleep(100000);

            static int i = 0;
            if (i++ >= 1000) { return; }

            auto subspec = static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1);
            TH1F& histo = processingContext.outputs().make<TH1F>(Output{ "MCH", "preclusters", subspec });
            histo.Fill(p / (double) producersAmount);
          }
        }
      };
      specs.push_back(producer);
    }

    MergerInfrastructureBuilder mergersBuilder;
    mergersBuilder.setInfrastructureName("histos");
    mergersBuilder.setInputSpecs(mergersInputs);
    mergersBuilder.setOutputSpec({{ "main" }, "MCH", "preclusters", 0 });
    MergerConfig config;
    config.inputObjectTimespan = { InputObjectsTimespan::LastDifference };
    config.publicationDecision = { PublicationDecision::EachNSeconds, 5 };
    config.mergedObjectTimespan = { MergedObjectTimespan::FullHistory };
    config.topologySize = { TopologySize::NumberOfLayers, 2 };
    mergersBuilder.setConfig(config);

    mergersBuilder.generateInfrastructure(specs);

    DataProcessorSpec printer{
      "printer-bins",
      Inputs{
        { "histo", "MCH", "preclusters", 0 }
      },
      Outputs{},
      AlgorithmSpec{
        (AlgorithmSpec::InitCallback) [](InitContext&) {
          return (AlgorithmSpec::ProcessCallback) [](ProcessingContext& processingContext) mutable {
//            LOG(INFO) << "printer invoked";
            auto histo = processingContext.inputs().get<TH1F*>("histo");
            std::string bins = "BINS:";
            for (int i = 1; i <= histo->GetNbinsX(); i++) {
              bins += " " + std::to_string((int) histo->GetBinContent(i));
            }
            LOG(INFO) << bins;
          };
        }
      }
    };
    specs.push_back(printer);
  }

  // custom merge
  {
//    WorkflowSpec specs; // enable comment to disable the workflow

    size_t producersAmount = 4;
    Inputs mergersInputs;
    for (size_t p = 0; p < producersAmount; p++) {
      mergersInputs.push_back({ "mo",               "MCH",
                                "preclusters",           static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1),
                                Lifetime::Timeframe });
      DataProcessorSpec producer{ "producer-custom" + std::to_string(p), Inputs{},
                                  Outputs{ { { "mo" },
                                             "MCH",
                                             "preclusters",
                                             static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1),
                                             Lifetime::Timeframe } },
                                  AlgorithmSpec{(AlgorithmSpec::ProcessCallback)[ p ](
                                    ProcessingContext & processingContext) mutable { usleep(100000);

            static int i = 0;
            if (i++ >= 1000) { return; }

            auto histo = std::make_unique<CustomMergeableTH2Quotient>();
            auto subspec = static_cast<o2::header::DataHeader::SubSpecificationType>(p + 1);
            processingContext.outputs().snapshot(OutputRef{ "mo", subspec }, *histo);
          }
        }
      };
      specs.push_back(producer);
    }

    MergerInfrastructureBuilder mergersBuilder;
    mergersBuilder.setInfrastructureName("custom");
    mergersBuilder.setInputSpecs(mergersInputs);
    mergersBuilder.setOutputSpec({{ "main" }, "MCH", "preclusters", 0 });
    MergerConfig config;
    config.inputObjectTimespan = { InputObjectsTimespan::LastDifference };
    config.publicationDecision = { PublicationDecision::EachNSeconds, 5 };
    config.mergedObjectTimespan = { MergedObjectTimespan::FullHistory };
    config.topologySize = { TopologySize::NumberOfLayers, 1 };
    mergersBuilder.setConfig(config);

    mergersBuilder.generateInfrastructure(specs);

    DataProcessorSpec printer{
      "printer-custom",
      Inputs{
        { "custom", "MCH", "preclusters", 0 }
      },
      Outputs{},
      AlgorithmSpec{
        (AlgorithmSpec::InitCallback) [](InitContext&) {
          return (AlgorithmSpec::ProcessCallback) [](ProcessingContext& processingContext) mutable {
            auto obj = processingContext.inputs().get<CustomMergeableTH2Quotient*>("custom");
           // LOG(INFO) << "SECRET:" << obj->getSecret();
          };
        }
      }
    };
    specs.push_back(printer);
  }

  return specs;
}
// clang-format on
