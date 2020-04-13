// Criado por Bernardo em 23/11/2017 e adaptado por Andressa em 12/08/2019.
//Os arquivos .cpp contêm as fontes das implementações em c++. Já os arquivos .h são os chamados headers, utilizados para extrair a declaração de funções, classes e outras declarações, que podem ser compartilhados entre vários arquivos com o código fonte. São, por isso, reutilizáveis. Esses arquivos possuem códigos que o compilador precisa para compilar outras partes

#include <algorithm> //biblioteca que define uma coleção de funções a serem usadas em uma gama de elementos
#include <numeric> //biblioteca que descreve um conjunto de algoritmos para executar determinadas operações em sequências de valores numéricos
#include <iostream> //biblioteca que define os objetos padrão de fluxo de entrada / saída
#include <iterator> //iterador é qualquer objeto que, apontando para algum elemento de um intervalo de elementos (como uma matriz ou um contêiner), tem a capacidade de iterar através dos elementos desse intervalo usando um conjunto de operadores (por exemplo, o operador de incremento ++)
#include <fstream> //biblioteca que fornece classes de fluxos de arquivos
#include <omp.h> //biblioteca que permite a programação multi-processo de memória compartilhada em múltiplas plataformas, de modo a acrescentar simultaneidade aos programas escritos em C++
#include <vector>

#ifdef  PARALLEL
#include <mpi.h> // biblioteca que permite a comunicação de dados em computação paralela
#endif

#include "Caesb.h" //include é uma diretiva de compilação
#include "../Controls/SeasonalMinEnvFlowControl.h"
#include "../Controls/StorageMinEnvFlowControl.h"
#include "../Controls/InflowMinEnvFlowControl.h"
#include "../Controls/FixedMinEnvFlowControl.h"
#include "../SystemComponents/WaterSources/ReservoirExpansion.h"
#include "../SystemComponents/WaterSources/SequentialJointTreatmentExpansion.h"
#include "../Simulation/Simulation.h"
#include "../SystemComponents/Bonds/LevelDebtServiceBond.h"
#include "../SystemComponents/Bonds/BalloonPaymentBond.h"
#include "../DroughtMitigationInstruments/Transfers.h"
#include "../SystemComponents/WaterSources/AllocatedReservoir.h"

/**
#include "../../../../../AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/usr/include/c++/7/cmath"
#include "../../../../../AppData/Local/Packages/CanonicalGroupLimited.Ubuntu18.04onWindows_79rhkp1fndgsc/LocalState/rootfs/usr/include/c++/7/vector"
**/


#ifdef PARALLEL //#ifdef é uma diretiva que permite que uma seção de um programa seja compilado somente se a macro especificada como parâmetro tiver sido definida, não importando qual seja o seu valor.
void Caesb::setProblemDefinition(BORG_Problem &problem) //void = vazio. O tipo void permite fazer funções que não retornam nada e funções que não têm parâmetros. Algumas funções não precisam retornar nenhum valor para funcionar, apenas realizar alguma ação.
{ //BORG é o algoritmo de otimização.(problem, número de identificação da variável de decisão, limite inferior da variável, limite superior da variável)
    // The parameter bounds are the same for all formulations
    BORG_Problem_set_bounds(problem, 0, 0.001, 1.0); //Gatilho para acionar a restrição de uso da água - descoberto
    BORG_Problem_set_bounds(problem, 1, 0.001, 1.0); //Gatilho para acionar a restrição de uso da água - tortoSM
    BORG_Problem_set_bounds(problem, 2, 0.001, 1.0); //Diferença entre um estágio de maior restrição para o estágio de menor restrição - descoberto
    BORG_Problem_set_bounds(problem, 3, 0.001, 1.0); //Diferença entre um estágio de maior restrição para o estágio de menor restrição - tortoSM
    BORG_Problem_set_bounds(problem, 4, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - descoberto
    BORG_Problem_set_bounds(problem, 5, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - tortoSM
    BORG_Problem_set_bounds(problem, 6, 0.0, 0.1); // Percentual da receita anual alocada para o fundo de contingência da companhia 0. O limite superior representa 10% da receita anual.
    BORG_Problem_set_bounds(problem, 7, 0.0, 0.1); // Percentual da receita anual alocada para o fundo de contingência da companhia 1. O limite superior representa 10% da receita anual.
    BORG_Problem_set_bounds(problem, 8, 0.001, 1.0); //Gatilho para acionar construção de infraestrutura pela Caesb - descoberto
    BORG_Problem_set_bounds(problem, 9, 0.001, 1.0); //Gatilho para acionar construção de infraestrutura pela Caesb - tortoSM
    BORG_Problem_set_bounds(problem, 10, 0.0, 1.0); //Ordem de "construção" do upgrade 1 da ETA Paranoá Sul (ampliação da capacidade de produção - 700 l/s)
    BORG_Problem_set_bounds(problem, 11, 0.0, 1.0); //Ordem de "construção" do upgrade 2 da ETA Paranoá Sul (ampliação da capacidade de produção - 700 l/s)
    BORG_Problem_set_bounds(problem, 12, 0.0, 1.0); //Ordem de "construção" do upgrade 3 da ETA Paranoá (ampliação da capacidade de produção - 700 l/s)
    BORG_Problem_set_bounds(problem, 13, 0.0, 1.0); //Ordem de "construção" do upgrade 1 da ETA Corumbá (implantação da ETA - 1400 l/s)
    BORG_Problem_set_bounds(problem, 14, 0.0, 1.0); //Ordem de "construção" do upgrade 2 da ETA Corumbá (ampliação da capacidade de produção - 1400 l/s)
    BORG_Problem_set_bounds(problem, 15, 0.0, 1.0); //Ordem de "construção" do upgrade 3 da ETA Corumbá (implantação de + 1.200 l/s só pra CAESB)
    BORG_Problem_set_bounds(problem, 16, 0.0, 1.0); //Ordem de "construção" da elevação do nível da barragem do Descoberto
    BORG_Problem_set_bounds(problem, 17, 0.0, 20.0); //Buffer de infraestrutura por parte da Caesb - descoberto
    BORG_Problem_set_bounds(problem, 18, 0.0, 20.0); //Buffer de infraestrutura por parte da Caesb - tortoSM

    // Set epsilons for objectives //(problem, n° de identificação da função objetivo, valor do epsilon). O valor do epsilon indica a precisão das funções objetivo.
    BORG_Problem_set_epsilon(problem, 0, 0.002);
    BORG_Problem_set_epsilon(problem, 1, 0.02);
    BORG_Problem_set_epsilon(problem, 2, 10.);
    BORG_Problem_set_epsilon(problem, 3, 0.02);
    BORG_Problem_set_epsilon(problem, 4, 0.01);
}
#endif


/**
 * Rodar o problema do Distrito Federal
 * @param vars
 * @param n_realizations
 * @param n_weeks
 * @param sol_number
 * @param output_directory

 */
int Caesb::functionEvaluation(double *vars, double *objs, double *consts) {
    double realization_start = omp_get_wtime();

//    cout << "Building Caesb Problem." << endl;
    // ===================== SET UP DECISION VARIABLES  =====================

    Simulation *s = nullptr; //nullptr: ponteiro nulo //*operador de dereferencia: a dereferencia busca o VALOR que está contida no endereço gravado no ponteiro
//    try {
    //throw invalid_argument("Test error");

    //VARIÁVEIS DE DECISÃO

    double caesb_descoberto_restriction_trigger = vars[0]; //gatilho para acionar restrição no Descoberto
    double caesb_tortoSM_restriction_trigger = vars[1]; //gatilho para acionar restrição no TortoSM
    double delta_descoberto_restriction_trigger = vars[2];
    double delta_tortoSM_restriction_trigger = vars[3];
    double caesb_descoberto_transfer_trigger = vars[4]; //gatilho para acionar transferência de água para o Descoberto
    double caesb_tortoSM_transfer_trigger = vars[5]; //gatilho para acionar transferência de água para o Descoberto
    double caesb_descoberto_annual_payment = vars[6]; // pagamento anual ao fundo de contingência. O valor é constante (igual para todo ano).
    double caesb_tortoSM_annual_payment = vars[7]; // pagamento anual ao fundo de contingência. O valor é constante (igual para todo ano).
    double caesb_descoberto_inftrigger = vars[8]; //gatilho para acionar a construção de nova infraestrutura por parte da Companhia Descoberto
    double caesb_tortoSM_inftrigger = vars[9]; //gatilho para acionar a construção de nova infraestrutura por parte da Companhia Torto/SM
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        caesb_descoberto_inftrigger = 1.1;
        caesb_tortoSM_inftrigger = 1.1;
    }
    double ETA_paranoaSul_upgrade1_ranking = vars[10]; // implantação de nova ETA no Paranoá Sul. É como se fossem o "low" e o "high" do estudo de caso da Carolina do Norte.
    double ETA_paranoaSul_upgrade2_ranking = vars[11]; // ampliação da capacidade da ETA Paranoá Sul
    double ETA_paranoaSul_upgrade3_ranking = vars[12]; // ampliação da capacidade da ETA Paranoá Sul
    double ETA_corumba_upgrade1_ranking = vars[13]; // ampliação da ETA Corumbá (+ 1400 l/s)
    double ETA_corumba_upgrade2_ranking = vars[14]; // ampliação da ETA Corumbá (+ 1400 l/s)
    double ETA_corumba_upgrade3_ranking = vars[15]; // ampliação da ETA Corumbá (+ 1200 l/s)
    double descoberto_expansao_ranking = vars[16]; // expansão da capacidade de armazenamento do reservatório do Descoberto
    //double tortoSM_descoberto_dupli_adut_ranking = vars[10]; // 3 // dupli adut = duplicação da adutora para aumentar a capacidade de transferência entre os sistemas
    double caesb_descoberto_inf_buffer = vars[17];
    double caesb_tortoSM_inf_buffer = vars[18];

    //ANALISAR POSSIBILIDADE DE INCLUIR O RIO DO SAL COMO OPÇÃO DE AMPLIAÇÃO DA INFRAESTRUTURA DE OFERTA

    //IDENTIFICADOR DE CADA INFRAESTRUTURA FUTURA. Obs: as infraestruturas já existentes devem ser numeradas antes, começando do 0.

    vector<infraRank> caesb_descoberto_infra_order_raw = { //A Companhia Caesb Descoberto abrange os reservatórios do Descoberto e de Corumbá IV
            infraRank(6, ETA_corumba_upgrade1_ranking),
            infraRank(7, ETA_corumba_upgrade2_ranking),
            infraRank(8, ETA_corumba_upgrade3_ranking),
            infraRank(12, descoberto_expansao_ranking)
    };

    vector<infraRank> caesb_tortoSM_infra_order_raw = { //A Companhia Caesb TortoSM abrange os reservatórios do TortoSM e do Lago Paranoá
            infraRank(9, ETA_paranoaSul_upgrade1_ranking),
            infraRank(10, ETA_paranoaSul_upgrade2_ranking),
            infraRank(11, ETA_paranoaSul_upgrade3_ranking)
    };

    // GET INFRASTRUCTURE CONSTRUCTION ORDER BASED ON DECISION VARIABLES
    sort(caesb_descoberto_infra_order_raw.begin(),
         caesb_descoberto_infra_order_raw.end(),
         by_xreal());
    sort(caesb_tortoSM_infra_order_raw.begin(),
         caesb_tortoSM_infra_order_raw.end(),
         by_xreal());

    vector<int> rof_triggered_infra_order_caesb_descoberto =
            vecInfraRankToVecInt(caesb_descoberto_infra_order_raw);
    vector<int> rof_triggered_infra_order_caesb_tortoSM =
            vecInfraRankToVecInt(caesb_tortoSM_infra_order_raw);

    // Create vectors with each utility's long-term ROF values assigned to all
    // infrastructure options.
    vector<double> rofs_infra_caesb_descoberto = vector<double>
            (rof_triggered_infra_order_caesb_descoberto.size(),
             caesb_descoberto_inftrigger);
    vector<double> rofs_infra_caesb_tortoSM = vector<double>
            (rof_triggered_infra_order_caesb_tortoSM.size(),
             caesb_tortoSM_inftrigger);

    // Remove small expansions being built after big expansions that would encompass the small expansions.

    // ==================== SET UP RDM FACTORS ============================
    // RDM factors são os fatores de grande incerteza (DU factors)
    if (utilities_rdm.empty()) {
        /// All matrices below have dimensions n_realizations x nr_rdm_factors
        utilities_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(5, 1.));
        water_sources_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(27,
                                               1.)); //n corresponde a: 1 + 2 * N° de infraestruturas.
        policies_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(2, 1.));
    }

    // ===================== SET UP PROBLEM COMPONENTS =====================

    //cout << "BEGINNING TRIANGLE TEST" << endl << endl;
    // cout << "Using " << omp_get_num_threads() << " cores." << endl;
    //    cout << getexepath() << endl;

    /// READ STREAMFLOWS
    int streamflow_n_weeks = 52 * (50 +
                                   40); //52 semanas em um ano / 50 anos de dados históricos / 40 anos de simulação futura.

    /// In case a vector containing realizations numbers to be calculated is passed, set
    /// number of realizations to number of realizations in that vector.

    //    vector<double> sewageFractions = Utils::parse1DCsvFile(
    //            io_directory + "/TestFiles/sewageFractions.csv");

    //SÉRIES DE EVAPORAÇÃO DE CADA RESERVATÓRIO

    EvaporationSeries evaporation_descoberto(&evap_descoberto,
                                             streamflow_n_weeks);
    EvaporationSeries evaporation_tortoSM(
            &evap_tortoSM,
            streamflow_n_weeks);
    EvaporationSeries evaporation_paranoa(
            &evap_paranoa,
            streamflow_n_weeks);
    EvaporationSeries evaporation_corumba(&evap_corumba,
            streamflow_n_weeks); //evaporação obtida por meio da evap do Paranoá (multiplicada por 4.478)

    // CRIAÇÃO DOS VETORES RELACIONADOS ÀS VAZÕES AFLUENTES DE CADA RESERVATÓRIO

    // Descoberto
    Catchment descoberto_inflows(&streamflows_descoberto,
                                 streamflow_n_weeks); //as vazões dos 6 afluentes foram somadas em um único arquivo

    vector<Catchment *> bacia_descoberto;
    bacia_descoberto.push_back(&descoberto_inflows);

    // Santa Maria
    Catchment tortoSM_inflows(&streamflows_tortoSM,
                              streamflow_n_weeks); //são 3 afluentes: Santa Maria, Milho Cozido e Vargem Grande

    vector<Catchment *> bacia_tortoSM;
    bacia_tortoSM.push_back(&tortoSM_inflows);

    //Sistema Bananal e Torto (Reforçam o Santa Maria) - CAP.TOR.001 (Torto)
    Catchment bananal_torto_inflows(&streamflows_bananal_torto,
                             streamflow_n_weeks); //esse subsistema reforça o Torto/Santa Maria -> a água captada
                                                    //no riberião Bananal e no Torto é bombeada para a ETA Brasília
    vector<Catchment *> sistema_bananal_torto;
    sistema_bananal_torto.push_back(&bananal_torto_inflows);

//    //Captação no Torto (Torto reforça o Santa Maria) - Captação CAP.TOR.001
//    Catchment torto_inflows(&streamflows_torto,
//                            streamflow_n_weeks);
//    vector<Catchment *> sistema_torto;
//    sistema_torto.push_back(&torto_inflows);

    // Lago Paranoá (Obs: o Subsistema Lago Norte capta água no Lago Paranoá, precisamente no braço do Torto)
    Catchment paranoa_inflows(&streamflows_paranoa,
                              streamflow_n_weeks); //são 5 afluentes: Ribeirão do Torto, Bananal, Riacho Fundo, Ribeirão Gama, Cabeça de Veado

    vector<Catchment *> bacia_paranoa;
    bacia_paranoa.push_back(&paranoa_inflows);

    // Corumbá IV
    Catchment corumbaIV_inflows(&streamflows_corumbaIV,
                                streamflow_n_weeks); //são 5 afluentes: Areias, Engenho das Lajes, Alagado, Fazenda Beira e Campo Limpo

    vector<Catchment *> bacia_corumba;
    bacia_corumba.push_back(&corumbaIV_inflows);

    // CURVAS VOLUME X ÁREA DOS RESERVATÓRIOS

    //Curva do Descoberto - baseado no arquivo de batimetria passado pela CAESB (volume útil em hm³)
    vector<double> descoberto_storage = {0,
                                         4.508 * table_gen_storage_multiplier,
                                         9.930 * table_gen_storage_multiplier,
                                         16.269 * table_gen_storage_multiplier,
                                         23.436 * table_gen_storage_multiplier,
                                         31.282 * table_gen_storage_multiplier,
                                         40.002 * table_gen_storage_multiplier,
                                         49.843 * table_gen_storage_multiplier,
                                         60.988 * table_gen_storage_multiplier,
                                         73.357 * table_gen_storage_multiplier,
                                         86.694 * table_gen_storage_multiplier,
                                         100.862 * table_gen_storage_multiplier}; //dados do volume (hm³) do reservatório do Descoberto
    vector<double> descoberto_area = {408.953, 494.839, 587.373, 678.283,
                                      751.075, 821.850, 923.983, 1049.369,
                                      1178.847, 1290.198, 1373.824,
                                      1461.072}; //dados da área (hm²) do reservatório do Descoberto (correspondente a cada volume acima)

    //Curva do Santa Maria - baseado no arquivo de batimetria passado pela CAESB (volume útil em hm³)
    vector<double> tortoSM_storage = {0, 3.584 * table_gen_storage_multiplier,
                                      7.507 * table_gen_storage_multiplier,
                                      11.758 * table_gen_storage_multiplier,
                                      16.322 * table_gen_storage_multiplier,
                                      21.228 * table_gen_storage_multiplier,
                                      24.878 * table_gen_storage_multiplier,
                                      32.122 * table_gen_storage_multiplier,
                                      38.122 * table_gen_storage_multiplier,
                                      44.504 * table_gen_storage_multiplier,
                                      51.299 * table_gen_storage_multiplier,
                                      58.507 * table_gen_storage_multiplier,
                                      66.090 * table_gen_storage_multiplier,
                                      73.954 * table_gen_storage_multiplier}; //dados do volume (hm³) do reservatório de Santa Maria

    vector<double> tortoSM_area = {340.830, 375.330, 409.180, 440.600, 472.690,
                                   508.480, 534.380, 580.740, 619.150, 658.450,
                                   700.400, 740.800, 774.640,
                                   793.160}; //dados da area (hm²) do reservatório de Santa Maria (correspondente a cada volume acima)

    //Curva do Paranoá - baseado na batimetria da CAESB, realizada em 2003 (volume total do lago)
    vector<double> paranoa_storage = {0, 0.0278 * table_gen_storage_multiplier,
                                      0.325 * table_gen_storage_multiplier,
                                      1.581 * table_gen_storage_multiplier,
                                      7.979 * table_gen_storage_multiplier,
                                      26.169 * table_gen_storage_multiplier,
                                      62.251 * table_gen_storage_multiplier,
                                      118.494 * table_gen_storage_multiplier,
                                      193.662 * table_gen_storage_multiplier,
                                      288.837 * table_gen_storage_multiplier,
                                      460.490 * table_gen_storage_multiplier}; //dados do volume (hm³) do reservatório do Paranoá

    vector<double> paranoa_area = {0, 2.757, 15.154, 81.656, 324.785, 763.788,
                                   1321.336, 1874.680, 2429.789, 3006.939,
                                   3881.452}; //dados da area (hm²) do reservatório do Paranoá (correspondente a cada volume acima)


    DataSeries descoberto_storage_area(&descoberto_storage,
                                       &descoberto_area); //aqui ele está juntando o vetor storage e o vetor area criado para cada reservatório acima
    DataSeries tortoSM_storage_area(&tortoSM_storage,
                                    &tortoSM_area);
    DataSeries paranoa_storage_area(&paranoa_storage,
                                    &paranoa_area);


    /// REGRAS RELACIONADAS ÀS VAZÕES REMANESCENTES (RESTRIÇÕES AMBIENTAIS)
    // Valores de vazão em hm³/semana

    // Vazão remanescente do Descoberto - baseado no artigo de Rocha e Cézar (2015)
    FixedMinEnvFlowControl descoberto_min_env_control(0, 0.6e-6 * 3600 * 24 * 7);

    // Vazão remanescente do Torto/Santa Maria - não tem, apenas verte a água
    FixedMinEnvFlowControl tortoSM_min_env_control(1, 0);

    // Vazão remanescente de Corumbá IV - baseado no EIA da UHE de Corumbá IV e na ação civil pública
    FixedMinEnvFlowControl corumba_min_env_control(2, 5.3e-6 * 3600 * 24 * 7);

    // Vazão remanescente do Paranoá - baseado na resolução n° 33/2018 da ADASA
    vector<int> paranoa_weeks = {0, 18, 44,
                                 53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    vector<double> paranoa_releases = {(1.2e-6 * 3600 * 24 * 7), (0.7e-6 * 3600 * 24 * 7),
                                       (1.2e-6 * 3600 * 24 * 7)}; // mínimo de 0,7 m³/s no período de estiagem e de 1,2 m³/s no período chuvoso
    SeasonalMinEnvFlowControl paranoa_min_env_control(3, paranoa_weeks,
                                                      paranoa_releases);

    //Vazão remanescente do Ribeirão Bananal e do Torto - combinação entre as duas water sources
    vector<int> bananal_torto_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44, 48,
                                        53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    vector<double> bananal_torto_releases = {(0.416e-6 * 3600 * 24 * 7), (0.503e-6 * 3600 * 24 * 7), (0.462e-6 * 3600 * 24 * 7),
                                             (0.437e-6 * 3600 * 24 * 7), (0.333e-6 * 3600 * 24 * 7), (0.249e-6 * 3600 * 24 * 7),
                                             (0.216e-6 * 3600 * 24 * 7), (0.593e-6 * 3600 * 24 * 7), (0.529e-6 * 3600 * 24 * 7),
                                             (0.179e-6 * 3600 * 24 * 7), (0.242e-6 * 3600 * 24 * 7), (0.341e-6 * 3600 * 24 * 7)};
    SeasonalMinEnvFlowControl bananal_torto_min_env_control(4, bananal_torto_weeks,
                                                      bananal_torto_releases);

//    // Vazão remanescente do Ribeirão Bananal - baseado no ato de outorga da ADASA para a captação no Bananal
//    vector<int> bananal_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44, 48,
//                                 53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
//    vector<double> bananal_releases = {(0.422e-6 * 3600 * 24 * 7), (0.446e-6 * 3600 * 24 * 7), (0.474e-6 * 3600 * 24 * 7),
//                                       (0.464e-6 * 3600 * 24 * 7), (0.406e-6 * 3600 * 24 * 7), (0.358e-6 * 3600 * 24 * 7),
//                                       (0.322e-6 * 3600 * 24 * 7), (0.286e-6 * 3600 * 24 * 7), (0.258e-6 * 3600 * 24 * 7),
//                                       (0.258e-6 * 3600 * 24 * 7), (0.314e-6 * 3600 * 24 * 7), (0.372e-6 * 3600 * 24 * 7)};
//    SeasonalMinEnvFlowControl bananal_min_env_control(4, bananal_weeks,
//                                                      bananal_releases);
//
//    // Vazão remanescente do Ribeirão do Torto - dado fornecido pela Caesb - exigência do PGIRH/DF (os valores correspondem a 20% da vazão média mínima mensal)
//    vector<int> torto_weeks = {0, 5, 9, 13, 18, 22, 26, 31, 35, 39, 44, 48,
//                               53}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
//    vector<double> torto_releases = {(0.410e-6 * 3600 * 24 * 7), (0.560e-6 * 3600 * 24 * 7), (0.450e-6 * 3600 * 24 * 7),
//                                     (0.410e-6 * 3600 * 24 * 7), (0.260e-6 * 3600 * 24 * 7), (0.140e-6 * 3600 * 24 * 7),
//                                     (0.110e-6 * 3600 * 24 * 7), (0.90e-6 * 3600 * 24 * 7), (0.80e-6 * 3600 * 24 * 7),
//                                     (0.100e-6 * 3600 * 24 * 7), (0.170e-6 * 3600 * 24 * 7), (0.310e-6 * 3600 * 24 * 7)};
//    SeasonalMinEnvFlowControl torto_min_env_control(5, torto_weeks,
//                                                      torto_releases);

    //FixedMinEnvFlowControl captacoes_gama();

    vector<MinEnvFlowControl *> min_env_flow_controls; //criação do vetor min_env_flow_controls
    min_env_flow_controls.push_back(
            &descoberto_min_env_control); //todos esses elementos (tortoSM, paranoa..) estão sendo jogados para dentro do vetor min_env_flow_controls
    min_env_flow_controls.push_back(&tortoSM_min_env_control);
    min_env_flow_controls.push_back(&paranoa_min_env_control);
    min_env_flow_controls.push_back(&corumba_min_env_control);
    min_env_flow_controls.push_back(&bananal_torto_min_env_control);
    //min_env_flow_controls.push_back(&torto_min_env_control);
    //min_env_flow_controls.push_back(&captacoes_gama_min_env_control);


    // CRIAÇÃO DOS RESERVATÓRIOS E VETORES CORRESPONDENTES
    vector<double> construction_time_interval = {3.0,
                                                 5.0}; //o período de construção das novas infraestruturas varia entre 3 e 5 anos
    vector<double> city_infrastructure_rof_triggers = {
            caesb_descoberto_inftrigger,
            caesb_tortoSM_inftrigger};
    //todas essas são variáveis de decisão definidas lá no começo do código, gatilho para acionar a construção de nova infraestrutura

    // EXISTING SOURCES //descrição das fontes de água já existentes

    Reservoir descoberto("Descoberto",//nome do reservatório
                         0,//número de identificação
                         bacia_descoberto,//vetor criado lá em cima - contém as vazões referentes aos afluentes do Descoberto
                         100.862 * table_gen_storage_multiplier, //capacidade de armazenamento do reservatório (hm³)
                         5.0e-6 * 3600 * 24 * 7, //capacidade máxima de tratamento da ETA Descoberto (hm³/semana) - e-mail da ADASA
                         evaporation_descoberto,                //obs: outorga da represa de Sta Maria é de 1.478 l/s (PDSB, 2017)
                         &descoberto_storage_area);

    Reservoir tortoSM("Torto / Santa Maria", 1,
                      bacia_tortoSM,
                      73.954 * table_gen_storage_multiplier, //hm³
                      1.1e-6 * 3600 * 24 * 7, //capacidade máxima de tratamento da ETA Brasília (hm³/semana)
                      evaporation_tortoSM,
                      &tortoSM_storage_area);

    // Reservatório de Corumbá IV
    double cIV_supply_caesb_capacity = 12.849 * table_gen_storage_multiplier; //valor obtido baseado na proporção das vazões retiradas
                                                                                // (conferir arquivo excel)
    double cIV_supply_saneago_capacity = 12.849 * table_gen_storage_multiplier;
    double cIV_energy_capacity = 745.702 * table_gen_storage_multiplier;
    double cIV_wq_capacity = 2936.6 * table_gen_storage_multiplier;
    double cIV_storage_capacity = cIV_wq_capacity + cIV_energy_capacity +
                                  cIV_supply_saneago_capacity +
                                  cIV_supply_caesb_capacity; //o armazenamento de água total é igual a soma da parte destinada à abastecimento,
                                                            // destinada à energia e da parte destinada à preservação ambiental

            //Curva de Corumbá IV - baseado nos dados do portal da ANA (volume útil)
    vector<double> corumba_storage = {0,
                                      2936.6 * table_gen_storage_multiplier, //hm³
                                      cIV_storage_capacity * table_gen_storage_multiplier};
    vector<double> corumba_area = {0, 13712.0, //hm²
                                   17330.0};
    DataSeries corumba_storage_area(&corumba_storage, &corumba_area);

    vector<int> cIV_allocations_ids = {0,
                                       WATER_QUALITY_ALLOCATION}; //0 é a id da companhia do Descoberto
    vector<double> cIV_allocation_fractions = {
            cIV_supply_caesb_capacity / cIV_storage_capacity,
            (cIV_wq_capacity + cIV_supply_saneago_capacity +
             cIV_energy_capacity) / cIV_storage_capacity};
    vector<double> cIV_treatment_allocation_fractions = {1.0};  //A companhia descoberto trata água do Corumbá IV. A companhia TortoSM não trata nada.

    AllocatedReservoir corumba("Corumba IV",
                               2,
                               bacia_corumba,
                               cIV_storage_capacity * table_gen_storage_multiplier,
                               0,
                               evaporation_corumba,
                               &corumba_storage_area,
                               &cIV_allocations_ids,
                               &cIV_allocation_fractions,
                               &cIV_treatment_allocation_fractions);

    // Lago Paranoá parameters
    double lp_supply_capacity = 36.966 *
                                table_gen_storage_multiplier; //volume destinado a abastecimento em hm³
    double lp_wq_capacity = 423.524 *
                            table_gen_storage_multiplier; //volume destinado a qualidade da água do lago em hm³
    double lp_storage_capacity = lp_wq_capacity + lp_supply_capacity;
    vector<int> lp_allocations_ids = {1,
                                      WATER_QUALITY_ALLOCATION}; //1 é a id da companhia do TortoSM
    vector<double> lp_allocation_fractions = {
            lp_supply_capacity / lp_storage_capacity,
            lp_wq_capacity / lp_storage_capacity};
    vector<double> lp_treatment_allocation_fractions = {1.0}; //A companhia descoberto não trata nenhuma água do Lago Paranoá. A companhia TortoSM trata.

    AllocatedReservoir paranoa("Lago Paranoa", 3,
                               bacia_paranoa,
                               lp_storage_capacity,
                               0.7e-6 * 3600 * 24 * 7, //capacidade de tratamento da ETA Lago Norte atual (hm³/semana)
                               evaporation_paranoa,
                               &paranoa_storage_area,
                               &lp_allocations_ids,
                               &lp_allocation_fractions,
                               &lp_treatment_allocation_fractions);

    Intake ribeirao_bananal_torto("Captacao no Ribeirao Bananal e Ribeiro do Torto",
                            4,                                                          //Obs: a série de vazão utilizada referente ao Bananal
                            sistema_bananal_torto,                                          //foi retirada de uma estação fluviométrica localizada
                            1.7e-6 * 3600 * 24 * 7); // hm³/semana      //a justante do ponto de captação. Não há problema,
                                                                                            // pois a captação começou apenas ao final de 2017,
                                                                                            // então a série é basicamente composta pela vazão natural do ribeirão.
//    Intake ribeirao_torto("Captacao no Ribeirao do Torto",
//                            5,
//                            sistema_torto,
//                            0.95e-6 * 3600 * 24 * 7); //950 m³/s corresponde à captação média no Torto.

    LevelDebtServiceBond dummy_bond(5, 1., 1, 1., vector<int>(1, 0));
    Reservoir dummy_endpoint("Dummy Node", 5, vector<Catchment *>(), 1., 0,
                             evaporation_corumba, 1,
                             construction_time_interval, 0, dummy_bond);

    // OPÇÕES DE INFRAESTRUTURAS FUTURAS - Descrição das opções futuras de infraestrutura para ampliar a oferta de água

    // The capacities listed here for expansions are what additional capacity is gained relative to existing capacity,
    // NOT the total capacity after the expansion is complete. For infrastructure with a high and low option, this means
    // the capacity for both is relative to current conditions - if Lake Michie is expanded small it will gain 2.5BG,
    // and a high expansion will provide 7.7BG more water than current. if small expansion happens, followed by a large
    // expansion, the maximum capacity through expansion is 7.7BG, NOT 2.5BG + 7.7BG.

    //OBS: Assumiu-se que todas as infraestruturas já tinham licença para serem construídas (permitting period = 0)
    /* Tipos de empréstimos definidos no WaterPaths:
        - Level Debt Service
        - Ballon Payment
        - Variable-Interest Bonds */

    //Captação em Corumbá IV - Etapa 1 (1.4 m³/s) + Etapa 2 (+ 1.4 m³/s) + Etapa 3 (+ 1.2 m³/s)

    vector<double> capacity_ETA_corumba_upgrade_1 = {1.4e-6 * 3600 * 24 * 7,
                                                     0}; //introdução da primeira ETA em Corumbá IV (1,4 m³/s). A unidade está em hm³/semana.
    vector<double> capacity_ETA_corumba_upgrade_2 = {1.4e-6 * 3600 * 24 * 7,
                                                     0}; //ampliação da capacidade de produção (+ 1,4 m³/s).
    vector<double> capacity_ETA_corumba_upgrade_3 = {1.2e-6 * 3600 * 24 * 7,
                                                     0}; //ampliação da capacidade de produção (+ 1,2 m³/s)

            //Empréstimo para Expansão da ETA Corumbá (Sistema Corumbá)

    vector<Bond *> debendure_expansao_ETA_corumba_1 = {
            new BalloonPaymentBond(11, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(6, 276521740.0, 20, 0.07,
                                     vector<int>(1, 0))};
    vector<Bond *> debendure_expansao_ETA_corumba_2 = {
            new BalloonPaymentBond(12, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(7, 222066142.8, 20, 0.07,
                                     vector<int>(1, 0))};
    vector<Bond *> debendure_expansao_ETA_corumba_3 = {
            new BalloonPaymentBond(13, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(8, 251383400, 20, 0.07,
                                     vector<int>(1, 0))};

    SequentialJointTreatmentExpansion ETA_corumba_etapa1(
            "Etapa 1 de Corumba IV", 6, 2, 0, {6, 7, 8},
            capacity_ETA_corumba_upgrade_1,
            debendure_expansao_ETA_corumba_1, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2020

    SequentialJointTreatmentExpansion ETA_corumba_etapa2(
            "Etapa 2 de Corumba IV", 7, 2, 1, {6, 7, 8},
            capacity_ETA_corumba_upgrade_2,
            debendure_expansao_ETA_corumba_2, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2030 a 2033

    SequentialJointTreatmentExpansion ETA_corumba_etapa3(
            "Etapa 3 de Corumba IV", 8, 2, 2, {6, 7, 8},
            capacity_ETA_corumba_upgrade_3,
            debendure_expansao_ETA_corumba_3, construction_time_interval,
            5 * WEEKS_IN_YEAR); //previsão: depois de 2037


    //Sistema Paranoá - Construção da ETA Paranoá Sul (0.7 m³/s), sua primeira ampliação (upgrade 2, com + 0.7 m³/s), segunda ampliação (upgrade 3, com + 0.35 m³/s) e
    // ampliação da ETA Lago Norte (upgrade 3, com + 0.35 m³/s))

    vector<double> capacities_ETA_paranoaSul_upgrade_1 = {0,
                                                          0.7e-6 * 3600 * 24 * 7}; //capacidade de produção da ETA paranoá Sul na sua etapa 1 = 0.7 m³/s
    vector<double> capacities_ETA_paranoaSul_upgrade_2 = {0,
                                                          0.7e-6 * 3600 * 24 * 7}; //aumento da capacidade de produção da ETA Paranoá Sul na sua etapa 2 = 0.7 m³/s
    vector<double> capacities_ETAs_paranoa_upgrade_3 = {0,
                                                          0.7e-6 * 3600 * 24 * 7}; // aumento da capacidade de produção da ETA paranoá Sul na sua etapa 3 = 0.350 m³/s
                                                                                    // + ampliação da ETA paranoá Norte em 0.350 m³/s

            // Empréstimos para a implantação e ampliação da ETA Paranoá Sul e Norte (Norte: apenas upgrade 3)

    vector<Bond *> debendure_expansao_ETA_paranoa_1 = {
            new BalloonPaymentBond(14, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(9, 60332016, 20, 0.07, vector<int>(1, 0))}; //bond term = 25 anos
    vector<Bond *> debendure_expansao_ETA_paranoa_2 = {
            new BalloonPaymentBond(15, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(10, 60332016, 20, 0.07, vector<int>(1, 0))};
    vector<Bond *> debendure_expansao_ETA_paranoa_3 = {
            new BalloonPaymentBond(16, 0, 20, 0.07, vector<int>(1, 0)),
            new LevelDebtServiceBond(11, 60332016, 20, 0.07, vector<int>(1, 0))};


    SequentialJointTreatmentExpansion etapa1_ETA_paranoaSul(
            "Etapa 1 da ETA Paranoa Sul ", 9, 3, 0, {9, 10, 11},
            capacities_ETA_paranoaSul_upgrade_1,
            debendure_expansao_ETA_paranoa_1, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2020

    SequentialJointTreatmentExpansion etapa2_ETA_paranoaSul(
            "Etapa 2 da ETA Paranoa Sul", 10, 3, 1, {9, 10, 11},
            capacities_ETA_paranoaSul_upgrade_2,
            debendure_expansao_ETA_paranoa_2, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2022

    SequentialJointTreatmentExpansion etapa3_ETAs_paranoa(
            "Etapa 3 da ETA Paranoa Sul e Norte", 11, 3, 2, {9, 10, 11},
            capacities_ETAs_paranoa_upgrade_3,
            debendure_expansao_ETA_paranoa_3, construction_time_interval,
            0 * WEEKS_IN_YEAR); //previsão: 2034 a 2037

    // Expansão do Reservatório do Descoberto

    LevelDebtServiceBond descoberto_exp_bond(12, 7541502, 20, 0.07, vector<int>(1,
                                                                                0)); //Elevação do nível d'água da barragem do descoberto (aumento da capacidade de armazenamento em 25%)
    ReservoirExpansion descoberto_expansion(
            "Expansao da capacidade de armazenamento do Descoberto", 12, 0,
            25.216, //25.216 = aumento em hm³ da capacidade de armazenamento do Descoberto (25%)
            construction_time_interval,
            0 * WEEKS_IN_YEAR,
            descoberto_exp_bond); //previsão: 2022 //acréscimo de 0.4 m³/s na vazão captável (PDSB, 2017).

    vector<WaterSource *> water_sources; //water_sources é um vetor comum, que comportará todas as
                                        // opções descritas acima de ampliação da infraestrutura de abastecimento
    water_sources.push_back(&descoberto);
    water_sources.push_back(&tortoSM);
    water_sources.push_back(&corumba);
    water_sources.push_back(&paranoa);
    water_sources.push_back(&ribeirao_bananal_torto);
//    water_sources.push_back(&ribeirao_bananal);
//    water_sources.push_back(&ribeirao_torto);

    water_sources.push_back(&ETA_corumba_etapa1);
    water_sources.push_back(&ETA_corumba_etapa2);
    water_sources.push_back(&ETA_corumba_etapa3);
    water_sources.push_back(&etapa1_ETA_paranoaSul);
    water_sources.push_back(&etapa2_ETA_paranoaSul);
    water_sources.push_back(&etapa3_ETAs_paranoa);
    water_sources.push_back(&descoberto_expansion);

    water_sources.push_back(&dummy_endpoint);

    //water_sources.push_back(&dummy_endpoint);

    /*
     * System connection diagram (water
     * flows from top to bottom)
     * Potential projects and expansions
     * of existing sources in parentheses
     *
     *                                  1 (Santa Maria)
     *      0(13)                       |
     *       \                          |
     *        \                         |  4(Ban)   4 (Torto)
     *         \                        |   \      /
     *          \                       |    \ __/
     *           \                      |     3 (Paranoá) (9, 10, 11)
     *            \                     |     |
     *             \                    |    |
     *              \                   |   |
     *               \                  |  |
     *                \                 |_|
     *                 \               /
     *                  \             /
     *                   \           /
     *                2(6, 7, 8)_   /
     *                           \ /
     *                           \/
     *                           5 (dummy)
     *
     */

    Graph g(6); //graph é uma forma de estruturar dados que consiste em dois componentes:
                    // um conjunto finito de vértices denominado de nós; e um cojunto finito de pares ordenados (x,y), denominados de edges.
    g.addEdge(0, 2); //essa conexão indica que existe uma aresta do vértice 0 ao vértice 2.
    g.addEdge(2, 5);
    g.addEdge(1, 5);
    g.addEdge(4, 3);
    //g.addEdge(5, 3);
    g.addEdge(3, 5);

    auto demand_n_weeks = (int) std::round(40 *
                                           WEEKS_IN_YEAR); //40 é o número de anos a serem simulados. A fç auto serve para
                                                            // declarar variáveis cujo tipo vai ser inferido pelo compilador a partir da inicialização delas

    // Criação dos vetores de vazão de retorno, relacionada aos efluentes lançados em corpos hídricos utilizados também para abastecimento urbano.

    vector<int> caesb_descoberto_ws_return_id = {
            3}; //ID do reservatório onde o efluente será lançado (montante desse reservatório) - 3 (ID do Paranoá)
    WwtpDischargeRule wwtp_discharge_caesb_descoberto( //em demand_to_waste_water_fraction, deve-se colocar a % de agua tratada (que sai da ETA - demanda) que vai voltar como vazão efluente para cada corpo receptor
            demand_to_wastewater_fraction_caesb_descoberto, // n° de séries (de demanda de vazão efluente) correspondentes ao n° de reservatórios onde os efluentes são lançados. Obs: se no meu estudo de caso, apenas o Paranoá for receptor, esse arquivo csv terá uma linha e 53 colunas (vazão efluente de cada semana)
            caesb_descoberto_ws_return_id); // id de cada reservatório/corpo receptor de efluente

    vector<int> caesb_tortoSM_ws_return_id = {
            3}; //ID do reservatório onde o efluente será lançado (montante desse reservatório) - 3 (ID do Paranoá)
    WwtpDischargeRule wwtp_discharge_caesb_tortoSM(
            demand_to_wastewater_fraction_caesb_tortoSM,
            caesb_tortoSM_ws_return_id); // id de cada reservatório/corpo receptor de efluente

    //Criação das companhias de água. A descrição de cada termo está no arquivo .doc.

    Utility caesb_descoberto((char *) "CAESB Descoberto", 0,
                             demand_caesb_descoberto, demand_n_weeks,
                             caesb_descoberto_annual_payment,
                             caesbDescobertoDemandClassesFractions,
                             caesbUserClassesWaterPrices,
                             wwtp_discharge_caesb_descoberto,
                             caesb_descoberto_inf_buffer,
                             rof_triggered_infra_order_caesb_descoberto,
                             vector<int>(), rofs_infra_caesb_descoberto, 0.04, //taxa de desconto
                             20, 0.07);

    Utility caesb_tortoSM((char *) "CAESB Torto/Santa Maria", 1,
                          demand_caesb_tortoSM, demand_n_weeks,
                          caesb_tortoSM_annual_payment,
                          caesbTortoSMDemandClassesFractions,
                          caesbUserClassesWaterPrices,
                          wwtp_discharge_caesb_tortoSM,
                          caesb_tortoSM_inf_buffer,
                          rof_triggered_infra_order_caesb_tortoSM,
                          vector<int>(), rofs_infra_caesb_tortoSM, 0.04, //0.04 = taxa de desconto
                          20, 0.07);


    vector<Utility *> utilities; //vetor que junta as companhias criadas acima
    utilities.push_back(&caesb_descoberto);
    utilities.push_back(&caesb_tortoSM);

    /// Water-source-utility connectivity matrix (each row corresponds to a utility and numbers are water
    /// sources IDs.
    vector<vector<int>> reservoir_utility_connectivity_matrix = {
            {0, 2, 6, 7, 8, 12},  //CAESB Descoberto
            {1, 3, 4, 9, 10, 11}  //CAESB Torto/Santa Maria
    };

//    @TODO: verificar se há necessidade de corrigir volumes de reservatórios construídos.
//    // O que table_storage_shift representa? O que são esses números (2000, 5000...) [3] [17]
    auto table_storage_shift = vector<vector<double>>(4,
                                                      vector<double>(25, 0.));
//    table_storage_shift[3][17] = 2000.; //tem a ver com a RdF
//    table_storage_shift[3][8] = 5000.;
//    table_storage_shift[1][14] = 100.;
//    table_storage_shift[1][20] = 500.;
//    table_storage_shift[1][21] = 500.;
//    table_storage_shift[1][15] = 700.;
//    table_storage_shift[1][9] = 700.;

    vector<DroughtMitigationPolicy *> drought_mitigation_policies; //criação do vetor das políticas de mitigação de seca (racionamento + transferência)

    // POLÍTICA DE RESTRIÇÃO

    vector<double> initial_restriction_triggers = {
            caesb_descoberto_restriction_trigger,
            caesb_tortoSM_restriction_trigger};
    //variáveis de decisão que representam os gatilhos para acionar o racionamento em cada companhia

    vector<double> restriction_stage_multipliers_caesb_descoberto = {0.98, 0.94,
                                                                     0.86}; //São 3 estágios de racionamento. Os fatores 0.98, 0.96, 0.88 são as restrições da demanda.
                                                                            // 0.98 significa que a demanda será restringida em 2% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_descoberto = {
            caesb_descoberto_restriction_trigger, //estágio 1 - 2% de redução (campanha de conscientização)
            caesb_descoberto_restriction_trigger +
            delta_descoberto_restriction_trigger, // estágio 2 - 4% de redução (tarifa de conting) + 2% do primeiro
            caesb_descoberto_restriction_trigger + 2 *
            delta_descoberto_restriction_trigger}; //estágio 3 - 12% de redução (racionamento) + 2% do primeiro
    //Obs: o 1° estágio corresponde à campanha, o 2° à tarifa de contingência e o 3° ao racionamento.
    //Acumulou-se a campanha com os outros dois estágios.
    //Como a tarifa de contingência é uma medida bastante penosa, bem como o racionamento,
    //optou-se por não acumular ambas medidas.

    vector<double> restriction_stage_multipliers_caesb_tortoSM = {0.98, 0.94,
                                                                  0.86};

    vector<double> restriction_stage_triggers_caesb_tortoSM = {
            caesb_tortoSM_restriction_trigger, //estágio 1
            caesb_tortoSM_restriction_trigger +
            delta_tortoSM_restriction_trigger, // estágio 2
            caesb_tortoSM_restriction_trigger +
            2 * delta_tortoSM_restriction_trigger}; //estágio 3

    // Criação das políticas de restrição de uso da água (puxa o código Restrictions.h, que contém os comandos dessa política)

    Restrictions restrictions_descoberto(0,
            restriction_stage_multipliers_caesb_descoberto,
            restriction_stage_triggers_caesb_descoberto,
            &caesbDescobertoDemandClassesFractions,
            &caesbUserClassesWaterPrices,
            &caesbPriceRestrictionMultipliers);

    Restrictions restrictions_tortoSM(1,
            restriction_stage_multipliers_caesb_tortoSM,
            restriction_stage_triggers_caesb_tortoSM,
            &caesbTortoSMDemandClassesFractions,
            &caesbUserClassesWaterPrices,
            &caesbPriceRestrictionMultipliers);

    drought_mitigation_policies = {&restrictions_descoberto, &restrictions_tortoSM};

    // POLÍTICA DE TRANSFERÊNCIA

    Graph transfer_graph_tortoSM_descoberto(2);
    transfer_graph_tortoSM_descoberto.addEdge(1,
                                              0); // Água do tortoSM para o Descoberto
    Transfers transfer_tortoSM_descoberto(0, 1, 1, 0.1, {0},
                                          {0.7e-6 * 3600 * 24 * 7}, //hm³/semana
                                          {caesb_descoberto_transfer_trigger}, //TortoSM transfere até 0.7 m³/s para o Descoberto
                                          transfer_graph_tortoSM_descoberto,
                                          vector<double>(), vector<int>());

    Graph transfer_graph_descoberto_tortoSM(2);
    transfer_graph_descoberto_tortoSM.addEdge(0,
                                              1); // Água do Descoberto para o TortoSM
    Transfers transfer_descoberto_tortoSM(1, 0, 0, 0.1, {1},
                                          {0.5e-6 * 3600 * 24 * 7},
                                          {caesb_tortoSM_transfer_trigger}, //Descoberto transfere até 0.5 m³/s para o TortoSM
                                          transfer_graph_descoberto_tortoSM,
                                          vector<double>(), vector<int>());

    drought_mitigation_policies.push_back(&transfer_tortoSM_descoberto);
    drought_mitigation_policies.push_back(&transfer_descoberto_tortoSM);

    /// Creates simulation object depending on use (or lack thereof) ROF tables
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run,
                           rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run,
                           rof_tables,
                           table_storage_shift,
                           rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    } else {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm,
                           water_sources_rdm,
                           policies_rdm,
                           n_weeks,
                           realizations_to_run);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, vars);
    }
//    double end_time = omp_get_wtime();
//	printf("Function evaluation time: %f s\n", end_time - start_time);

    double realization_end = omp_get_wtime();
    std::cout << "Simulation took  " << realization_end - realization_start
          << "s" << std::endl;

    /// Calculate objectives and store them in Borg decision variables array.
#ifdef  PARALLEL
    objectives = calculateAndPrintObjectives(false);

        int i = 0;
        objs[i] = min(objectives[i], objectives[5 + i]) - 1.;
        for (i = 1; i < 5; ++i) {                       //são 5 objetivos que serão otimizados (0 a 4)
            objs[i] = max(objectives[i], objectives[5 + i]);
        }

        if (s != nullptr) {	 // != significa "diferente de"
            delete s;
    }
    s = nullptr;
#endif
//    } catch (const std::exception& e) {
//        simulationExceptionHander(e, s, objs, vars);
//	return 1;
//    }

    delete s;

    return 0;
}

int Caesb::simulationExceptionHander(const std::exception &e,
                                     Simulation *s, // :: significa "resolução de escopo"
                                     double *objs, const double *vars) {
    int num_dec_var = 19; //número de variáveis desse estudo de caso
//        printf("Exception called during calculations. Decision variables are below:\n");
    ofstream sol;
    int world_rank;

#ifdef  PARALLEL
    // int mpi_initialized;
    // MPI_Initialized(&mpi_initialized);
    // if (mpi_initialized)
 //            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    // else
        world_rank = 0;
#else
    world_rank = 0;
#endif
    string error_file = "sol_error_rank_" + to_string(world_rank) + ".csv";
    sol.open(error_file.c_str());
    for (int i = 0; i < num_dec_var; ++i) {
        sol << vars[i] << ",";
    }
    sol << flush;
    sol.close();
    printf("Error. Decision variables printed in %s\n", error_file.c_str());

#ifdef PARALLEL
    objs[0] = 0.;
    objs[1] = 1.1;
    objs[2] = 1000;
    objs[3] = 5.;
    objs[4] = 5.;

    if (s != nullptr) {
        delete s;
        s = nullptr;
    }
#else
    Utils::print_exception(e);
#endif

    return 1;
}


Caesb::~Caesb() = default;

Caesb::Caesb(unsigned long n_weeks, int import_export_rof_table)
        : Problem(n_weeks) {
    if (import_export_rof_table == EXPORT_ROF_TABLES) {
        table_gen_storage_multiplier = BASE_STORAGE_CAPACITY_MULTIPLIER;
    } else {
        table_gen_storage_multiplier = 1.;
    }
}


void
Caesb::readInputData() { //A partir dessa linha serão inseridos os dados de entrada para que o modelo possa funcionar
    cout << "Reading input data." << endl;
    string data_dir = DEFAULT_DATA_DIR + BAR;

#pragma omp parallel default(none) num_threads(omp_get_thread_num())
    {
#pragma omp single
        streamflows_descoberto = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "descoberto_inflows.csv", n_realizations);

#pragma omp single
        streamflows_tortoSM = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "tortoSM_inflows.csv", n_realizations);
#pragma omp single
        streamflows_bananal_torto = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "bananal_torto_inflows.csv", n_realizations);
//#pragma omp single
//        streamflows_torto = Utils::parse2DCsvFile( //inserção dos dados de vazão
//                io_directory + DEFAULT_DATA_DIR + "inflows" +
//                evap_inflows_suffix +
//                BAR + "torto_inflows.csv", n_realizations);
        // }
#pragma omp single
        streamflows_paranoa = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "paranoa_inflows.csv", n_realizations);

#pragma omp single
        streamflows_corumbaIV = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" +
                evap_inflows_suffix +
                BAR + "corumbaIV_inflows.csv", n_realizations);

// };
        //cout << "Reading evaporations." << endl;
#pragma omp single
        evap_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "descoberto_evap.csv",
                n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "tortoSM_evap.csv",
                n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_paranoa = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "paranoa_evap.csv",
                n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_corumba = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" +
                evap_inflows_suffix +
                BAR + "corumba_evap.csv",
                n_realizations); //inserção dos dados de evaporação


        //cout << "Reading demands." << endl;
#pragma omp single
        demand_caesb_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" +
                evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_descoberto_demand.csv", n_realizations);

        demand_caesb_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" +
                evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_tortoSM_demand.csv", n_realizations);


        //cout << "Reading others." << endl;
#pragma omp single
        {
            demand_to_wastewater_fraction_caesb_descoberto = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "demand_to_wastewater_fraction_caesb_descoberto.csv"); //demanda de efluentes da caesb

            demand_to_wastewater_fraction_caesb_tortoSM = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "demand_to_wastewater_fraction_caesb_tortoSM.csv"); //demanda de efluentes da caesb

            caesbDescobertoDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbDescobertoDemandClassesFractions.csv"); //demanda de cada categoria de usuário da caesb

            caesbTortoSMDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbTortoSMDemandClassesFractions.csv"); //demanda de cada categoria de usuário da caesb

            caesbUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbUserClassesWaterPrices.csv"); //tarifa de água para cada categoria de usuário da caesb
            //OBS: inserir no valor a tarifa de esgoto também na última coluna desse arquivo (100% da de água)
            caesbPriceRestrictionMultipliers = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR +
                    "caesbPriceRestrictionMultipliers.csv"); //% de aumento da tarifa para cada categoria durante o racionamento
        }
//    cout << "Done reading input data." << endl;
    }

}
