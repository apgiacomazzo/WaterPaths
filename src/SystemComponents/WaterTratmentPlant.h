//
// Created by Bernardo on 5/3/2020.
//

#ifndef WATERPATHS_WATERTRATMENTPLANT_H
#define WATERPATHS_WATERTRATMENTPLANT_H

#include <vector>
#include <Bonds/Base/Bond.h>

using namespace std;

class WaterTreatmentPlant {
private:
    WaterTreatmentPlant(const int id, int connectedUtility, int connectedSource,
                        double totalTreatmentCapacity);

    const int id;
    vector<int> connected_sources;
    vector<int> connected_utilities;
    vector<double> allocated_treatment_capacities;
    vector<double> allocated_treatment_fractions;

    double total_treatment_capacity;
    vector<Bond *> bond;

public:
    WaterTreatmentPlant(int id, const vector<int> &connectedSources,
                        const vector<int> &connectedUtilities,
                        const vector<double> &allocatedTreatmentFractions,
                        double totalTreatmentCapacity);

    WaterTreatmentPlant(int id, const vector<int> &connectedSources,
                        const vector<int> &connectedUtilities,
                        const vector<double> &allocatedTreatmentFractions,
                        double totalTreatmentCapacity,
                        vector<Bond *> bond);
};


#endif //WATERPATHS_WATERTRATMENTPLANT_H
