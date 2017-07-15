//
// Created by bernardo on 2/6/17.
//

#ifndef TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H
#define TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H

#include "../../SystemComponents/Utility.h"
#include "../../Utils/Constants.h"
#include "../../Utils/Graph/Graph.h"
#include "../../Utils/Matrices.h"

class DroughtMitigationPolicy {
protected:
    DroughtMitigationPolicy(const DroughtMitigationPolicy &drought_mitigation_policy);

    vector<int> utilities_ids;
    vector<Utility *> realization_utilities;
    const Matrix3D<double> *storage_to_rof_table_;

public:
    const int id;
    const int type;

    DroughtMitigationPolicy(const int id, const int type);

    virtual void applyPolicy(int week)= 0;

    virtual void addSystemComponents(vector<Utility *> utilities, vector<WaterSource *> water_sources)= 0;

    const vector<int> &getUtilities_ids() const;

    bool operator<(const DroughtMitigationPolicy *other);

    bool operator>(const DroughtMitigationPolicy *other);

    virtual ~DroughtMitigationPolicy();

    void setStorage_to_rof_table_(const Matrix3D<double> *storage_to_rof_table_);

};


#endif //TRIANGLEMODEL_DROUGHTMITIGATIONPOLICY_H