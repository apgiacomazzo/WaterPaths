//
// Created by Bernardo on 5/3/2020.
//

#include "WaterTratmentPlant.h"

WaterTreatmentPlant::WaterTreatmentPlant(const int id,
                                         const vector<int> &connectedSources,
                                         const vector<int> &connectedUtilities,
                                         const vector<double> &allocatedTreatmentFractions,
                                         double totalTreatmentCapacity) : id(
        id), connected_sources(connectedSources), connected_utilities(
        connectedUtilities), allocated_treatment_fractions(
        allocatedTreatmentFractions), total_treatment_capacity(
        totalTreatmentCapacity) {}

WaterTreatmentPlant::WaterTreatmentPlant(const int id,
                                         const vector<int> &connectedSources,
                                         const vector<int> &connectedUtilities,
                                         const vector<double> &allocatedTreatmentFractions,
                                         double totalTreatmentCapacity,
                                         const vector<Bond *> bond) :
        id(id), connected_sources(connectedSources),
        connected_utilities(connectedUtilities),
        allocated_treatment_fractions(allocatedTreatmentFractions),
        total_treatment_capacity(totalTreatmentCapacity),
        bond(bond) {}

WaterTreatmentPlant::WaterTreatmentPlant(const int id, int connectedUtility, int connectedSource, double totalTreatmentCapacity):
    id(id), connected_sources(vector<int>(1, connectedSource)),
    connected_utilities(vector<int>(1, connectedUtility)),
    total_treatment_capacity(totalTreatmentCapacity) {}