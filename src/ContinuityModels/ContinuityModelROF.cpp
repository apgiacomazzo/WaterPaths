//
// Created by bernardo on 1/26/17.
//

#include <iostream>
#include "ContinuityModelROF.h"

ContinuityModelROF::ContinuityModelROF(const vector<WaterSource *> &water_source,
                                       const vector<vector<int>> &water_source_connectivity_matrix,
                                       const vector<Utility *> &utilities,
                                       const vector<vector<int>> &water_source_utility_connectivity_matrix,
                                       const int rof_type, const int realization_id) : ContinuityModel(water_source,
                                                                             water_source_connectivity_matrix,
                                                                             utilities,
                                                                             water_source_utility_connectivity_matrix),
                                                             rof_type(rof_type), realization_id(realization_id) {
    // calculate total combined storage volumes (here, it equals sum of reservoirs capacities.
    for (Utility *u : this->utilities) {
        u->updateTotalStoredVolume();
    }
}

/**
 * Runs one the full rof calculations for realization #realization_id for a given week.
 * @param week for which rof is to be calculated.
 */
vector<double> ContinuityModelROF::calculateROF(int week) {

    // vector where risks of failure will be stored.
    vector<double> risk_of_failure(utilities.size(), 0.0);
    vector<double> year_failure(utilities.size(), 0.0);

    // short-term rof calculations.
    if (rof_type == SHORT_TERM_ROF) {

        // perform a continuity simulation for 50 yearly realization.
        for (int r = 0; r < NUMBER_REALIZATIONS_ROF; ++r) {
            for (int w = week; w < week + WEEKS_ROF_SHORT_TERM; ++w) {

                // one week continuity time-step.
                this->continuityStep(w, r);

                // check total available storage for each utility and, if smaller than the fail ration,
                // increase the rof of that utility to 1 / (50 [yearly realizations] * 52 [short-term rof])
                for (int i = 0; i < utilities.size(); ++i) {
                    if (utilities[i]->getStorageToCapacityRatio() <= STORAGE_CAPACITY_RATIO_FAIL)
                        year_failure[i] = FAILURE;
                }
            }

            for (int j = 0; j < utilities.size(); ++j) {
                risk_of_failure[j] += year_failure[j];
                year_failure[j] = NON_FAILURE;
            }

            // reset reservoirs' and utilities' storage and combined storage, respectively, they currently
            // have in the corresponding realization simulation.
            this->resetUtilitiesAndReservoirs();
        }
    }

    for (int i = 0; i < utilities.size(); ++i) {
        risk_of_failure[i] /= NUMBER_REALIZATIONS_ROF;
    }

    return risk_of_failure;
}

/**
 * reset reservoirs' and utilities' storage and last release, and combined storage, respectively, they
 * currently have in the corresponding realization simulation.
 */
void ContinuityModelROF::resetUtilitiesAndReservoirs() {

    // update water sources info.
    for (int i = 0; i < water_sources.size(); ++i) {
        water_sources[i]->setAvailable_volume(water_sources_realization[i]->getAvailable_volume());
        water_sources[i]->setOutflow_previous_week(water_sources_realization[i]->getOutflow_previous_week());
    }

    // update utilities combined storage.
    for (Utility *u : this->utilities) {
        u->updateTotalStoredVolume();
    }
}

/**
 * Pass to the rof continuity model the locations of the utilities of the realization it calculated rofs
 * for. This is done so that the rof calculations can reset the reservoir and utilities states for each of
 * the rofs' year realizations.
 * @param water_sources_realization
 */
void ContinuityModelROF::setWater_sources_realization(const vector<WaterSource *> &water_sources_realization) {
    ContinuityModelROF::water_sources_realization = water_sources_realization;
}
