//
// Created by bernardoct on 5/3/17.
//

#ifndef TRIANGLEMODEL_WATERREUSE_H
#define TRIANGLEMODEL_WATERREUSE_H


#include "Base/WaterSource.h"

class WaterReuse : public WaterSource {
private:
    double treated_volume;
public:
    WaterReuse(const char *name, const int id, const double capacity);

    WaterReuse(
            const char *name, const int id, const double capacity,
            const double construction_rof_or_demand,
            const vector<double> &construction_time_range,
            double permitting_period,
            double construction_cost_of_capital, double bond_term,
            double bond_interest_rate);

    void applyContinuity(
            int week, double upstream_source_inflow,
            vector<double> *demand_outflow) override;

    WaterReuse &operator=(const WaterReuse &water_reuse);

    double getReused_volume() const;

};

#endif //TRIANGLEMODEL_WATERREUSE_H
