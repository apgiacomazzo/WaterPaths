//
// Criado por Bernardo em 23/11/2017 e adaptado por Andressa em 21/08/2019.
//
//DISTRITO FEDERAL MODEL
#ifndef TRIANGLEMODEL_CAESB_H
#define TRIANGLEMODEL_CAESB_H

#ifdef PARALLEL
#include "../../Borg/borgms.h"
#endif
#include "Base/Problem.h"
#include "../Simulation/Simulation.h"

using namespace std;

class Caesb : public Problem {
private:
    const int n_utilities = 1;
//obs: não inseri os streamflows do subsistema gama (provavelmente será excluído da análise)
    vector<vector<double>> streamflows_chacara89; //descoberto
    vector<vector<double>> streamflows_chapadinha; //descoberto
    vector<vector<double>> streamflows_olaria; //descoberto
    vector<vector<double>> streamflows_rodeador; //descoberto
    vector<vector<double>> streamflows_capao_comp; //descoberto
    vector<vector<double>> streamflows_pedras; //descoberto
    vector<vector<double>> streamflows_descoberto; //descoberto
    vector<vector<double>> streamflows_santa_maria; //tortoSM
    vector<vector<double>> streamflows_milho_cozido; //tortoSM
    vector<vector<double>> streamflows_vargem_grande; //tortoSM
    vector<vector<double>> streamflows_tortoSM; //tortoSM
    vector<vector<double>> streamflows_bananal; //subsistema bananal e lago paranoa
    vector<vector<double>> streamflows_torto; //subsistema lago norte e lago paranoa
    vector<vector<double>> streamflows_riacho_fundo; //lago paranoa
    vector<vector<double>> streamflows_gama; //lago paranoa
    vector<vector<double>> streamflows_cabeca_veado; //lago paranoa
    vector<vector<double>> streamflows_paranoa; //lago paranoa
    vector<vector<double>> streamflows_areias; //corumba
    vector<vector<double>> streamflows_engenho_lajes; //corumba
    vector<vector<double>> streamflows_alagado; //corumba
    vector<vector<double>> streamflows_fazenda_beira; //corumba
    vector<vector<double>> streamflows_campo_limpo; //corumba
    vector<vector<double>> streamflows_corumbaIV; //corumba
    vector<vector<double>> demand_caesb_descoberto;
    vector<vector<double>> demand_caesb_tortoSM;
    vector<vector<double>> demand_caesb_paranoa;
    vector<vector<double>> demand_caesb_corumbaIV;
    vector<vector<double>> evap_descoberto;
    vector<vector<double>> evap_tortoSM;
    vector<vector<double>> evap_paranoa;
    vector<vector<double>> evap_corumba;
    vector<vector<double>> demand_to_wastewater_fraction_caesb_descoberto;
    vector<vector<double>> demand_to_wastewater_fraction_caesb_tortoSM;
    vector<vector<double>> caesbDescobertoDemandClassesFractions;
    vector<vector<double>> caesbTortoSMDemandClassesFractions;
    vector<vector<double>> caesbUserClassesWaterPrices;
    vector<vector<double>> caesbPriceRestrictionMultipliers;

public:
    Caesb(unsigned long n_weeks, int import_export_rof_table);

    ~Caesb();

#ifdef PARALLEL
    void setProblemDefinition(BORG_Problem &problem);
#endif

    int functionEvaluation(double *vars, double *objs, double *consts) override;

    int simulationExceptionHander(const std::exception &e, Simulation *s,
                                  double *objs, const double *vars);

    void readInputData();

};


#endif //TRIANGLEMODEL_CAESB_H
