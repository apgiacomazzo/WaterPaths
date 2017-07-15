//
// Created by bernardo on 1/22/17.
//

#ifndef TRIANGLEMODEL_WATERSOURCE_H
#define TRIANGLEMODEL_WATERSOURCE_H

#include <string>
#include "../../Catchment.h"
#include "../../../Utils/Constants.h"

using namespace std;
using namespace Constants;
const int BOND_INTEREST_PAYMENTS_PER_YEAR = 1;

class WaterSource {
protected:
    double available_volume = 0;
    double total_outflow = 0;
    double upstream_source_inflow = 0;
    double upstream_catchment_inflow = 0;
    double demand = 0;
    double upstream_min_env_inflow;
    double capacity;
    double *available_allocated_volumes = new double();;
    bool online;
    vector<Catchment *> catchments;
    double min_environmental_outflow;

    virtual void applyContinuity(
            int week, double upstream_source_inflow,
            vector<double> *demand_outflow) = 0;

    void bypass(int week, double upstream_source_inflow);

public:
    const int id;
    const char *name;
    const int source_type;
    const double raw_water_main_capacity;
    const double construction_rof;
    const double construction_cost_of_capital;
    const double construction_time;
    const double bond_term;
    const double bond_interest_rate;

    WaterSource(
            const char *name, const int id,
            const vector<Catchment *> &catchments, const double capacity,
            const double raw_water_main_capacity, const int source_type);

    WaterSource(
            const char *source_name, const int id,
            const vector<Catchment *> &catchments, const double capacity,
            const double raw_water_main_capacity, const int source_type,
            const double construction_rof,
            const vector<double> construction_time_range,
            double construction_cost_of_capital, double bond_term,
            double bond_interest_rate);

    WaterSource(const WaterSource &water_source);

    virtual ~WaterSource();

    WaterSource &operator=(const WaterSource &water_source);

    bool operator<(const WaterSource *other);

    bool operator>(const WaterSource *other);

    void continuityWaterSource(
            int week, double upstream_source_inflow,
            vector<double> *demand_outflow);

    virtual void setOnline();

    virtual double getAvailableAllocatedVolume(int utility_id);

    double getMin_environmental_outflow() const;

    void setMin_environmental_outflow(double min_environmental_outflow);

    double getAvailable_volume() const;

    double getTotal_outflow() const;

    bool isOnline() const;

    virtual void setFull();

    void setOutflow_previous_week(double outflow_previous_week);

    virtual double getCapacity();

    virtual void addCapacity(double capacity);

    double getUpstream_source_inflow() const;

    double getDemand() const;

    double getUpstreamCatchmentInflow() const;

    double calculateNetPresentConstructionCost(
            int week, double
    discount_rate, double *level_debt_service_payment) const;

    static bool compare(WaterSource *lhs, WaterSource *rhs);

    virtual void setRealization(unsigned long r);

    virtual void removeWater(int allocation_id, double volume);

    virtual double getAllocatedCapacity(int utility_id);

    virtual double getAllocatedFraction(int utility_id);
};


#endif //TRIANGLEMODEL_WATERSOURCE_H