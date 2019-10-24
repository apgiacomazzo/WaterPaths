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
    BORG_Problem_set_bounds(problem, 2, 0.001, 1.0); //Gatilho para acionar a restrição de uso da água - paranoa
    BORG_Problem_set_bounds(problem, 3, 0.001, 1.0); //Gatilho para acionar a restrição de uso da água - corumbaIV
    BORG_Problem_set_bounds(problem, 4, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - descoberto
    BORG_Problem_set_bounds(problem, 5, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - tortoSM
    BORG_Problem_set_bounds(problem, 6, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - paranoa
    BORG_Problem_set_bounds(problem, 7, 0.001, 1.0); //Gatilho para acionar a transferência de água entre sistemas - corumba IV
    BORG_Problem_set_bounds(problem, 8, 0.0, 0.1); // Percentual da receita anual alocada para o fundo de contingência da Caesb. O limite superior representa 10% da receita anual.
    BORG_Problem_set_bounds(problem, 9, 0.001, 1.0); //Gatilho para acionar construção de infraestrutura pela Caesb - descoberto
    BORG_Problem_set_bounds(problem, 10, 0.001, 1.0); //Gatilho para acionar construção de infraestrutura pela Caesb - paranoa
    BORG_Problem_set_bounds(problem, 11, 0.001, 1.0); //Gatilho para acionar construção de infraestrutura pela Caesb - corumba IV
    BORG_Problem_set_bounds(problem, 12, 0.0, 1.0); //Ordem de "construção" do upgrade 1 da ETA Paranoá Sul (ampliação da capacidade de produção - 700 l/s)
    BORG_Problem_set_bounds(problem, 13, 0.0, 1.0); //Ordem de "construção" do upgrade 2 da ETA Paranoá Sul (ampliação da capacidade de produção - 700 l/s)
    BORG_Problem_set_bounds(problem, 14, 0.0, 1.0); //Ordem de "construção" do upgrade 3 da ETA Paranoá Sul (ampliação da capacidade de produção - 350 l/s)
    BORG_Problem_set_bounds(problem, 15, 0.0, 1.0); //Ordem de "construção" do upgrade da ETA Lago Norte (ampliação da capacidade de produção - 350 l/s)
    BORG_Problem_set_bounds(problem, 16, 0.0, 1.0); //Ordem de "construção" do upgrade 1 da ETA Valparaíso de Corumbá IV (implantação de + 1.400 l/s só pra CAESB)
    BORG_Problem_set_bounds(problem, 17, 0.0, 1.0); //Ordem de "construção" do upgrade 2 da ETA Valparaíso de Corumbá IV (implantação de + 1.200 l/s só pra CAESB)
    BORG_Problem_set_bounds(problem, 18, 0.0, 1.0); //Ordem de "construção" da elevação do nível da barragem do Descoberto
    //BORG_Problem_set_bounds(problem, , 0.0, 1.0); //Ordem de "construção" da duplicação da adutora que conecta os sistemas Torto/Santa Maria ao Descoberto
    //BORG_Problem_set_bounds(problem, , 0.0, 0.5); //Percentual de água de Corumbá IV alocada para a caesb
    //BORG_Problem_set_bounds(problem, , 0.001, 1.0); //Fração da ETA de Corumbá IV que vai para a Caesb
    BORG_Problem_set_bounds(problem, 19, 0.0, 20.0); //Buffer de infraestrutura por parte da Caesb

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
//    cout << "Building Caesb Problem." << endl;
    // ===================== SET UP DECISION VARIABLES  =====================

    Simulation *s = nullptr; //nullptr: ponteiro nulo //*operador de dereferencia: a dereferencia busca o VALOR que está contida no endereço gravado no ponteiro
//    try {
    //throw invalid_argument("Test error");
    double caesb_descoberto_restriction_trigger = vars[0];
    double caesb_tortoSM_restriction_trigger = vars[1];
    double caesb_paranoa_restriction_trigger = vars[2];
    double caesb_corumbaIV_restriction_trigger = vars[3];
    double caesb_descoberto_transfer_trigger = vars[4];
    //double caesb_tortoSM_transfer_trigger = vars[5]; //excluir?
    //double caesb_paranoa_transfer_trigger = vars[6]; //excluir?
    //double caesb_corumbaIV_transfer_trigger = vars[7];//excluir?
    double caesb_annual_payment = vars[8]; // pagamento anual ao fundo de contingência
    double caesb_descoberto_inftrigger = vars[9];
    double caesb_tortoSM_inftrigger = vars[10];
    //double caesb_corumbaIV_inftrigger = vars[11];
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        caesb_descoberto_inftrigger = 1.1;
        caesb_tortoSM_inftrigger = 1.1;
        //caesb_corumbaIV_inftrigger = 1.1;
    }
    double ETA_paranoaSul_upgrade1_ranking = vars[12]; // 6 - essas etapas se referem a ampliação da capacidade de produção de água (implantação de nova ETA ou ampliação de alguma já existente). É como se fossem o "low" e o "high" do estudo de caso da Carolina do Norte.
    double ETA_paranoaSul_upgrade2_ranking = vars[13]; // 7
    double ETA_paranoaSul_upgrade3_ranking = vars[14]; // 8
    double ETA_lagoNorte_upgrade_ranking = vars[15]; // 9 - ampliação da ETA Lago Norte (Captação no Lago Paranoá)
    double ETA_corumba_upgrade1_ranking = vars[16]; // 11 - ampliação da ETA Valparaíso (+ 1400 l/s)
    double ETA_corumba_upgrade2_ranking = vars[17]; // 12 - ampliação da ETA Valparaíso (+ 1200 l/s)
    double descoberto_expansao_ranking = vars[18]; // 10
    //double tortoSM_descoberto_dupli_adut_ranking = vars[10]; // 3 // dupli adut = duplicação da adutora para aumentar a capacidade de transferência entre os sistemas
    //double caesb_corumba = vars[14];
    //double eta_corumba_frac_caesb = vars[15]; // Não sei se essa entra como uma variável de decisão, porque acredito que esse valor já seja definido (50% da água tratada pela ETA vai pra caesb e os outros 50% vai para Corumbá.
    double caesb_inf_buffer = vars[19];

    //ANALISAR POSSIBILIDADE DE INCLUIR O RIO DO SAL COMO OPÇÃO DE AMPLIAÇÃO DA INFRAESTRUTURA DE OFERTA

    //Identificador de cada infraestrutura. As infraestruturas já existentes devem ser numeradas antes, começando do 0.
    vector<infraRank> caesb_descoberto_infra_order_raw = {
            infraRank(10, descoberto_expansao_ranking),
            infraRank(11, ETA_corumba_upgrade1_ranking),
            infraRank(12, ETA_corumba_upgrade2_ranking)
            };

    vector<infraRank> caesb_tortoSM_infra_order_raw = {
            infraRank(6, ETA_paranoaSul_upgrade1_ranking),
            infraRank(7, ETA_paranoaSul_upgrade2_ranking),
            infraRank(8, ETA_paranoaSul_upgrade3_ranking),
            infraRank(9, ETA_lagoNorte_upgrade_ranking)
            };

            //infraRank(, tortoSM_descoberto_dupli_adut_ranking), //deletar - acredito que essa medida não resulte em ampliação da oferta

    //unidade oficial: l/s
    /*
    double added_production_ETA_paranoaSul_upgrade1 = 700; // ampliação da capacidade de produção da ETA Paranoá (+ 700 l/s) na primeira etapa
    double added_production_ETA_paranoaSul_upgrade2 = 700;
    double added_production_ETA_paranoaSul_upgrade3 = 350;
    double added_production_ETA_lagoNorte_upgrade = 350;
    double added_production_ETA_corumba_upgrade1 = 1400; //1400 = só pra caesb; 2800 = caesb + saneago
    double added_production_ETA_corumba_upgrade2 = 1200; //1200 = só pra caesb; 2400 = caesb + saneago
    */

    /// get infrastructure construction order based on decision variables
    sort(caesb_descoberto_infra_order_raw.begin(), //durham foi substituído por caesb (excluí os correspondentes a owasa e raleigh)
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
    vector<double> rofs_infra_caesb_descoberto = vector<double> //durham foi substituído por caesb (excluí os correspondentes a owasa e raleigh)
            (rof_triggered_infra_order_caesb_descoberto.size(), caesb_descoberto_inftrigger);
    vector<double> rofs_infra_caesb_tortoSM = vector<double>
            (rof_triggered_infra_order_caesb_tortoSM.size(), caesb_tortoSM_inftrigger);

    /// Remove small expansions being built after big expansions that would
    /// encompass the small expansions.
    /*
    added_production_paranoa_etapa4 =
            checkAndFixInfraExpansionHighLowOrder(
                    &rof_triggered_infra_order_caesb,
                    &rofs_infra_caesb,
                    7,
                    8,
                    9,
                    added_production_paranoa_etapa2,
                    added_production_paranoa_etapa3,
                    added_production_paranoa_etapa4);

    added_production_corumba_etapa3 =
            checkAndFixInfraExpansionHighLowOrder(
                    &rof_triggered_infra_order_caesb,
                    &rofs_infra_caesb,
                    10,
                    11,
                    added_production_corumba_etapa2,
                    added_production_corumba_etapa3);
    */

    // ==================== SET UP RDM FACTORS ============================
    // RDM factors são os fatores de grande incerteza
    if (utilities_rdm.empty()) {
        /// All matrices below have dimensions n_realizations x nr_rdm_factors
        utilities_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(1, 1.));
        water_sources_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(51, 1.)); //alterar número "51" para (1 + 2 * N° de novas infraestruturas), caso eu opte por não fazer as análises de deep uncertainty.
        policies_rdm = std::vector<vector<double>>(
                n_realizations, vector<double>(1, 1.));
    }

    // ===================== SET UP PROBLEM COMPONENTS =====================

    //cout << "BEGINNING TRIANGLE TEST" << endl << endl;
    // cout << "Using " << omp_get_num_threads() << " cores." << endl;
    //    cout << getexepath() << endl;

    /// Read streamflows
    int streamflow_n_weeks = 52 * (50 + 40); //52 semanas em um ano / 50 anos de dados históricos / 40 anos de simulação futura.

    /// In case a vector containing realizations numbers to be calculated is passed, set
    /// number of realizations to number of realizations in that vector.

    //    vector<double> sewageFractions = Utils::parse1DCsvFile(
    //            io_directory + "/TestFiles/sewageFractions.csv");

    EvaporationSeries evaporation_descoberto(&evap_descoberto, streamflow_n_weeks); //substituí durham por descoberto
    EvaporationSeries evaporation_tortoSM(
            &evap_tortoSM, //substituí jordan_lake por tortoSM
            streamflow_n_weeks);
    EvaporationSeries evaporation_paranoa(
            &evap_paranoa, //substituí falls_lake por paranoa
            streamflow_n_weeks);
    EvaporationSeries evaporation_corumba(&evap_corumba, streamflow_n_weeks); //substituí owasa por corumba

    // Create catchments and corresponding vectors

    // Descoberto
    /*Catchment descoberto_chacara89(&streamflows_chacara89, streamflow_n_weeks); //
    Catchment chapadinha_aviario(&streamflows_chapadinha, streamflow_n_weeks);
    Catchment olaria(&streamflows_olaria, streamflow_n_weeks);
    Catchment ribeirao_rodeador(&streamflows_rodeador, streamflow_n_weeks);
    Catchment capao_comprido(&streamflows_capao_comp, streamflow_n_weeks);
    Catchment ribeirao_pedras(&streamflows_pedras, streamflow_n_weeks); */
    Catchment descoberto_inflows(&streamflows_descoberto, streamflow_n_weeks); //as vazões dos 6 afluentes foram somadas em um único arquivo

    vector<Catchment *> bacia_descoberto;
    /*bacia_descoberto.push_back(descoberto_chacara89);
    bacia_descoberto.push_back(chapadinha_aviario);
    bacia_descoberto.push_back(olaria);
    bacia_descoberto.push_back(ribeirao_rodeador);
    bacia_descoberto.push_back(capao_comprido);
    bacia_descoberto.push_back(ribeirao_pedras);*/
    bacia_descoberto.push_back(&descoberto_inflows);

    //Subsistema Gama - reforça o Descoberto - Contribui com uma vazão de 320 l/s para o Descoberto (adotar valor fixo - mais simples)
    // Verificar a possibilidade de somar todas as vazões desses córregos e deixar um único arquivo apenas.
    // Problema: as medições do alagado começam em 2007, as demais em 1993 e 1985.
    /*Catchment ponte_terra_2(&streamflows_ponte_terra_2, streamflow_n_weeks); //ponte terra 2 e 3 são captações no mesmo córrego, apenas em ptos distintos
    Catchment ponte_terra_3(&streamflows_ponte_terra_3, streamflow_n_weeks);
    Catchment olhos_dagua(&streamflows_olhos_dagua, streamflow_n_weeks);
    Catchment alagado_gama(&streamflows_alagado_gama, streamflow_n_weeks);
    Catchment crispim(&streamflows_crispim, streamflow_n_weeks);*/

    /*vector<Catchment *> subsistema_gama;
    subsistema_gama.pushback(ponte_terra_2);
    subsistema_gama.pushback(ponte_terra_3);
    subsistema_gama.pushback(olhos_dagua);
    subsistema_gama.pushback(alagado_gama);
    subsistema_gama.pushback(crispim); */

    // Torto/Santa Maria
    //Catchment santa_maria(&streamflows_santa_maria, streamflow_n_weeks);
    //Catchment milho_cozido(&streamflows_milho_cozido, streamflow_n_weeks);
    //Catchment vargem_grande(&streamflows_vargem_grande, streamflow_n_weeks);
    Catchment tortoSM_inflows(&streamflows_tortoSM, streamflow_n_weeks); //são 3 afluentes: Santa Maria, Milho Cozido e Vargem Grande

    vector<Catchment *> bacia_tortoSM;
    /*bacia_tortoSM.push_back(santa_maria);
    bacia_tortoSM.push_back(milho_cozido);
    bacia_tortoSM.push_back(vargem_grande);*/
    bacia_tortoSM.push_back(&tortoSM_inflows);

    //Subsistema Bananal
    Catchment bananal_inflow(&streamflows_bananal, streamflow_n_weeks); //esse subsistema reforça o Torto/Santa Maria -> a água captada
                                                                    //no riberião Bananal é bombeada para a ETA Brasília
    vector<Catchment *> subsistema_bananal;
    subsistema_bananal.push_back(&bananal_inflow);

    // Lago Paranoá (Obs: o Subsistema Lago Norte capta água no Lago Paranoá, precisamente no braço do Torto)
    /*Catchment ribeirao_torto(&streamflows_torto, streamflow_n_weeks);
    Catchment bananal(&streamflows_bananal, streamflow_n_weeks);
    Catchment riacho_fundo(&streamflows_riacho_fundo, streamflow_n_weeks);
    Catchment ribeirao_gama(&streamflows_gama, streamflow_n_weeks);
    Catchment cabeca_veado(&streamflows_cabeca_veado, streamflow_n_weeks);*/
    Catchment paranoa_inflows(&streamflows_paranoa, streamflow_n_weeks); //são 5 afluentes: Ribeirão do Torto, Bananal, Riacho Fundo, Ribeirão Gama, Cabeça de Veado

    vector<Catchment *> bacia_paranoa;
    /*bacia_paranoa.push_back(ribeirao_torto);
    bacia_paranoa.push_back(bananal);
    bacia_paranoa.push_back(riacho_fundo);
    bacia_paranoa.push_back(ribeirao_gama);
    bacia_paranoa.push_back(cabeca_veado);*/
    bacia_paranoa.push_back(&paranoa_inflows);

    // Corumbá IV
    /*Catchment areias(&streamflows_areias, streamflow_n_weeks);
    Catchment engenho_lajes(&streamflows_engenho_lajes, streamflow_n_weeks);
    Catchment alagado(&streamflows_alagado, streamflow_n_weeks);
    Catchment fazenda_beira(&streamflows_fazenda_beira, streamflow_n_weeks);
    Catchment campo_limpo(&streamflows_campo_limpo, streamflow_n_weeks);*/
    Catchment corumbaIV_inflows(&streamflows_corumbaIV, streamflow_n_weeks); //são 5 afluentes: Areias, Engenho das Lajes, Alagado, Fazenda Beira e Campo Limpo

    vector<Catchment *> bacia_corumba;
    /*bacia_corumba.push_back(areias);
    bacia_corumba.push_back(engenho_lajes);
    bacia_corumba.push_back(alagado);
    bacia_corumba.push_back(fazenda_beira);
    bacia_corumba.push_back(campo_limpo);*/
    bacia_corumba.push_back(&corumbaIV_inflows);


    // Storage vs. area reservoir curves.
    vector<double> descoberto_storage = {0, 7230000 * table_gen_storage_multiplier, 14460000 * table_gen_storage_multiplier, 21690000 * table_gen_storage_multiplier,
                                         28920000 * table_gen_storage_multiplier, 36150000 * table_gen_storage_multiplier, 43370000 * table_gen_storage_multiplier,
                                         50600000 * table_gen_storage_multiplier, 57830000 * table_gen_storage_multiplier, 65060000 * table_gen_storage_multiplier,
                                         72290000 * table_gen_storage_multiplier}; //dados do volume (m³) do reservatório do Descoberto
    vector<double> descoberto_area = {0, 4918367, 5436090, 5768617, 6050209, 6319930, 6601218, 6893733, 7210723, 7547564, 7909190}; //dados da area (m²) do reservatório do Descoberto (correspondente a cada volume acima)
    vector<double> tortoSM_storage = {0, 6130000 * table_gen_storage_multiplier, 12260000 * table_gen_storage_multiplier, 18390000 * table_gen_storage_multiplier,
                                      24520000 * table_gen_storage_multiplier, 30650000 * table_gen_storage_multiplier, 36780000 * table_gen_storage_multiplier,
                                      42920000 * table_gen_storage_multiplier, 49050000 * table_gen_storage_multiplier, 55180000 * table_gen_storage_multiplier,
                                      61310000 * table_gen_storage_multiplier}; //dados do volume (m³) do reservatório de Santa Maria
    vector<double> tortoSM_area = {0, 3980519, 4227586, 4399522, 4549165, 4693721, 4845850, 5002331, 5168599, 5346899, 5533394}; //dados da area (m²) do reservatório de Santa Maria (correspondente a cada volume acima)
    vector<double> paranoa_storage = {0, 27773 * table_gen_storage_multiplier, 325307 * table_gen_storage_multiplier, 1580840 * table_gen_storage_multiplier,
                                      7979477 * table_gen_storage_multiplier, 26169300 * table_gen_storage_multiplier, 62251080 * table_gen_storage_multiplier,
                                      118494265 * table_gen_storage_multiplier, 193662316 * table_gen_storage_multiplier, 288836784 * table_gen_storage_multiplier,
                                      460489632 * table_gen_storage_multiplier}; //dados do volume (m³) do reservatório do Paranoá
    vector<double> paranoa_area = {0, 27568, 151539, 816557, 3247850, 7637880, 13213364, 18746798, 24297885, 30069390, 38814515}; //dados da area (m²) do reservatório do Paranoá (correspondente a cada volume acima)
    vector<double> corumba_storage = {0, 2936600000 * table_gen_storage_multiplier, 3708000000 * table_gen_storage_multiplier}; //só obtive esses dados
    vector<double> corumba_area = {0, 137120000, 173300000}; //só obtive esses dados

    DataSeries descoberto_storage_area(&descoberto_storage, &descoberto_area); //aqui ele está juntando o vetor storage e o vetor area criado para cada reservatório acima
    DataSeries tortoSM_storage_area(&tortoSM_storage,
                                           &tortoSM_area);
    DataSeries paranoa_storage_area(&paranoa_storage,
                                 &paranoa_area);
    DataSeries corumba_storage_area(&corumba_storage,
                                         &corumba_area);


    /// Minimum environmental flow rules (controls) //se houver dúvidas, olhar o arquivo do Triangle.cpp inalterado

    FixedMinEnvFlowControl descoberto_min_env_control(0, 0.6); //0.6 m³/s
    FixedMinEnvFlowControl tortoSM_min_env_control(1, 0); //o Santa Maria não tem vazão remanescente (apenas verte)

    vector<int> paranoa_weeks = {18, 44}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    vector<double> paranoa_releases = {0.7, 1.2}; // mínimo de 0,7 m³/s no período de estiagem e de 1,2 m³/s no período chuvoso
    SeasonalMinEnvFlowControl paranoa_min_env_control(4, paranoa_weeks, paranoa_releases);

    FixedMinEnvFlowControl corumba_min_env_control(2, 16.8); //encontrei apenas essa informação de 16,8 m³/s - EIA da UHE de Corumbá IV

    vector<int> bananal_weeks = {18, 44}; // período de estiagem (maio - week 18, a outubro) e período chuvoso (novembro - week 44 a abril)
    vector<double> bananal_releases = {0.36, 0.5}; //vazões baseadas no gráfico de vazões remanescentes do córrego Bananal no PGIRH
    SeasonalMinEnvFlowControl bananal_min_env_control(3, bananal_weeks, bananal_releases);

    //FixedMinEnvFlowControl captacoes_gama();


    vector<MinEnvFlowControl *> min_env_flow_controls; //criação do vetor min_env_flow_controls
    min_env_flow_controls.push_back(&descoberto_min_env_control); //todos esses elementos (durham, falls, wheeler e etc) estão sendo jogados para dentro do vetor min_env_flow_controls
    min_env_flow_controls.push_back(&tortoSM_min_env_control);
    min_env_flow_controls.push_back(&paranoa_min_env_control);
    min_env_flow_controls.push_back(&corumba_min_env_control);
    min_env_flow_controls.push_back(&bananal_min_env_control);
    //min_env_flow_controls.push_back(&captacoes_gama_min_env_control);


    /// Create reservoirs and corresponding vector
    vector<double> construction_time_interval = {3.0, 5.0}; //o período de construção das novas infraestruturas varia entre 3 e 5 anos
    vector<double> city_infrastructure_rof_triggers = {caesb_descoberto_inftrigger,
                                                       caesb_tortoSM_inftrigger
                                                       }; //todas essas são variáveis de decisão definidas lá no começo do código, gatilho para acionar a construção de nova infraestrutura

    vector<double> bond_term = {25}; //tempo até o empréstimo ser pago na totalidade
    vector<double> bond_rate = {0.05};//taxa de juros do empréstimo
    double discount_rate = 0.05;


    // Existing Sources //descrição das fontes de água já existentes
    Reservoir descoberto("Descoberto",//nome do reservatório
                                0,//número de identificação
                                bacia_descoberto,//vetor criado lá em cima - contém as vazões referentes aos afluentes do Descoberto
                                72290000 * table_gen_storage_multiplier,//capacidade de armazenamento do reservatório
                                ILLIMITED_TREATMENT_CAPACITY,
                                evaporation_descoberto,
                                &descoberto_storage_area);

    Reservoir tortoSM("Torto / Santa Maria", 1,
                                bacia_tortoSM,
                                61310000 * table_gen_storage_multiplier,
                                ILLIMITED_TREATMENT_CAPACITY,
                                evaporation_tortoSM,
                                &tortoSM_storage_area);

    // Corumbá IV parameters
    double cIV_supply_caesb_capacity = 14924.0 * table_gen_storage_multiplier;
    double cIV_supply_saneago_capacity = 14924.0 * table_gen_storage_multiplier;
    double cIV_energy_capacity = 14924.0 * table_gen_storage_multiplier;
    double cIV_wq_capacity = 30825.0 * table_gen_storage_multiplier;
    double cIV_storage_capacity = cIV_wq_capacity + cIV_energy_capacity + cIV_supply_saneago_capacity + cIV_supply_caesb_capacity; //o armazenamento de água total é igual a soma da parte destinada a abastecimento, destinada a energia e e da parte destinada a preservação ambiental
    vector<int> cIV_allocations_ids = {0, WATER_QUALITY_ALLOCATION}; //0 é a id da companhia do Descoberto
    vector<double> cIV_allocation_fractions = {
            cIV_supply_caesb_capacity / cIV_storage_capacity,
            cIV_supply_saneago_capacity / cIV_storage_capacity,
            cIV_energy_capacity / cIV_storage_capacity,
            cIV_wq_capacity / cIV_storage_capacity};
    vector<double> cIV_treatment_allocation_fractions = {1.0, 0.0};  //A companhia descoberto água do Corumbá IV. A companhia TortoSM não trata nada.

    AllocatedReservoir corumba("Corumba IV", 2, //colocar como alocated reservoir
                                bacia_corumba,
                                cIV_storage_capacity,
                                1.4, //capacidade de tratamento da ETA Corumbá atualmente (1.4 m³/s)
                                evaporation_corumba,
                                &corumba_storage_area,
                                &cIV_allocations_ids,
                                &cIV_allocation_fractions,
                                &cIV_treatment_allocation_fractions);

    /// Lago Paranoá parameters
    double lp_supply_capacity = 36965656 * table_gen_storage_multiplier; //volume destinado a abastecimento em m³
    double lp_wq_capacity = 423523976 * table_gen_storage_multiplier; //volume destinado a qualidade da água do lago em m³
    double lp_storage_capacity = lp_wq_capacity + lp_supply_capacity;
    vector<int> lp_allocations_ids = {1, WATER_QUALITY_ALLOCATION}; //1 é a id da companhia do TortoSM
    vector<double> lp_allocation_fractions = {
            lp_supply_capacity / lp_storage_capacity,
            lp_wq_capacity / lp_storage_capacity};
    vector<double> lp_treatment_allocation_fractions = {0.0, 1.0}; //A companhia descoberto não trata nenhuma água do Lago Paranoá. A companhia TortoSM trata.

    AllocatedReservoir paranoa("Lago Paranoa", 4,
                                bacia_paranoa,
                                lp_storage_capacity,
                                0.7, //capacidade de tratamento da ETA Lago Norte atual = 0,7 m³/s
                                evaporation_paranoa,
                                &paranoa_storage_area,
                                &lp_allocations_ids,
                                &lp_allocation_fractions,
                                &lp_treatment_allocation_fractions);

    Intake ribeirao_bananal("Captacao no Ribeirao Bananal", 3,
                            subsistema_bananal,
                            0.73); // 0,73 m³/s

    /*Intake captacoes_gama("Subsistema Gama", 5,
                                subsistema_gama,
                                0.31, //capacidade de produção de água (m³/s)
                                ILLIMITED_TREATMENT CAPACITY);*/


    // Potential Sources
    // The capacities listed here for expansions are what additional capacity is gained relative to existing capacity,
    // NOT the total capacity after the expansion is complete. For infrastructure with a high and low option, this means
    // the capacity for both is relative to current conditions - if Lake Michie is expanded small it will gain 2.5BG,
    // and a high expansion will provide 7.7BG more water than current. if small expansion happens, followed by a large
    // expansion, the maximum capacity through expansion is 7.7BG, NOT 2.5BG + 7.7BG.

    //OBS: os custos estão em unidade de milhão. - ALTERAR OS PERMITTING PERIOD

    LevelDebtServiceBond descoberto_exp_bond(10, 7.5, 25, 0.05, vector<int>(1, 0)); //Elevação do nível d'água da barragem do descoberto (aumento da capacidade de armazenamento em 25%)
    ReservoirExpansion descoberto_expansion("Expansao da capacidade de armazenamento do Descoberto", 10, 0, 0.4, construction_time_interval, //22775000 = aumento em m³ da capacidade de armazenamento do Descoberto
                                     17 * WEEKS_IN_YEAR, descoberto_exp_bond); //previsão: 2022 //acréscimo de 0.4 m³/s na vazão captável (PDSB, 2017).

    //Empréstimo para Expansão da ETA Corumbá (Sistema Corumbá)

    double cost_ETA_corumba_upgrade_1 = 220.844;
    double cost_ETA_corumba_upgrade_2 = 150.2; //dado aleatório, inserido apenas pare fins de teste

    vector<Bond *> bonds_ETA_corumba_upgrade_1; //empréstimo para ampliação 1 da ETA (menor capacidade)
    int uid = 0;
    for (double &cost : cost_ETA_corumba_upgrade_1) {
        bonds_ETA_corumba_upgrade_1.emplace_back(new LevelDebtServiceBond(11 + uid, cost, 25, 0.05, vector<int>(1, 0))); //alterar 25, 0.05
        uid++;
    }
    vector<Bond *> bonds_ETA_corumba_upgrade_2; //empréstimo para ampliação 2 da ETA (maior capacidade)
    uid = 0;
    for (double &cost : cost_ETA_corumba_upgrade_2) {
        bonds_ETA_corumba_upgrade_2.emplace_back(new LevelDebtServiceBond(12 + uid, cost, 25, 0.05, vector<int>(1, 0))); //alterar 25, 0.05
        uid++;
    }
    /// Expansão da ETA Corumbá (Sistema Corumbá)
    double capacity_ETA_corumba_upgrade_1 = 1.4; //ampliação da capacidade de produção (+ 1,4 m³/s)
    double capacity_ETA_corumba_upgrade_2 = 1.2; //ampliação da capacidade de produção (+ 1,2 m³/s)

    SequentialJointTreatmentExpansion ETA_corumba_etapa2("Etapa 2 de Corumba IV", 11, 2, 0, {11, 12}, capacity_ETA_corumba_upgrade_1,
                                                bonds_ETA_corumba_upgrade_1, construction_time_interval, 12 * WEEKS_IN_YEAR); //previsão: 2030 a 2033

    SequentialJointTreatmentExpansion ETA_corumba_etapa3("Etapa 3 de Corumba IV", 12, 2, 1, {11, 12}, capacity_ETA_corumba_upgrade_2,
                                                  bonds_ETA_corumba_upgrade_2, construction_time_interval, 12 * WEEKS_IN_YEAR); //previsão: depois de 2037
    //Obs: No arquivo SequentialJointTreatmentExpansion.h, o argumento vector<double> sequential_treatment_capacity foi alterado para
    //double sequential_treatment_capacity (visto que, para esse estudo de caso, há apenas uma companhia).

    //Sistema Paranoá - Construção da ETA Paranoá Sul (0.7 m³/s), sua primeira ampliação (upgrade 2, com + 0.7 m³/s), segunda ampliação (upgrade 3, com + 0.35 m³/s) e
    // ampliação da ETA Lago Norte (upgrade 1, com + 0.35 m³/s))

    double capacities_ETA_paranoaSul_upgrade_1 = 0.7; //capacidade de produção da ETA paranoá Sul na sua etapa 1 = 0.7 m³/s

    double cost_ETA_paranoaSul_upgrade_1 = 60; //custo do upgrade 1 da ETA Paranoá Sul = 60 milhões

    double capacities_ETA_paranoaSul_upgrade_2 = (1.4 - 0.7); //aumento da capacidade de produção da ETA Paranoá Sul na sua etapa 2 = 0.7 m³/s

    double cost_ETA_paranoaSul_upgrade_2 = (120 - 60); //custo do upgrade 2 da ETA Paranoá Sul = 60 milhões

    double capacities_ETA_paranoaSul_upgrade_3 = 0.35; // aumento da capacidade de produção da ETA paranoá Sul na sua etapa 3 = 0.350 m³/s

    double cost_ETA_paranoaSul_upgrade_3 = (150 - 120); //custo do upgrade 3 da ETA Paranoá Sul = 30 milhões


    /// Bonds ETA Paranoá Sul
    vector<Bond *> ETA_paranoaSul_bonds_capacity_1; //criação do vetor relacionado ao empréstimo para implantação da ETA Paranoá Sul com capacidade 1
    uid = 0;
    for (double &cost : cost_ETA_paranoaSul_upgrade_1) {
        ETA_paranoaSul_bonds_capacity_1.emplace_back(new LevelDebtServiceBond(6 + uid, cost, 25, 0.05, vector<int>(1, 0)));
        uid++;
    }
    vector<Bond *> ETA_paranoaSul_bonds_capacity_2; //criação do vetor relacionado ao empréstimo para implantação da ETA Paranoá Sul com capacidade 2
    uid = 0;
    for (double &cost : cost_ETA_paranoaSul_upgrade_2) {
        ETA_paranoaSul_bonds_capacity_2.emplace_back(new LevelDebtServiceBond(7 + uid, cost, 25, 0.05, vector<int>(1, 0)));
        uid++;
    }
    vector<Bond *> ETA_paranoaSul_bonds_capacity_3; //criação do vetor relacionado ao empréstimo para implantação da ETA Paranoá Sul com capacidade 2
    uid = 0;
    for (double &cost : cost_ETA_paranoaSul_upgrade_3) {
        ETA_paranoaSul_bonds_capacity_3.emplace_back(new LevelDebtServiceBond(8 + uid, cost, 25, 0.05, vector<int>(1, 0)));
        uid++;
    }
    /// ETA Paranoá Sul (Sistema Paranoá)
    SequentialJointTreatmentExpansion etapa1_ETA_paranoaSul("Etapa 1 da ETA Paranoa Sul ", 6, 4, 0, {6, 7, 8}, capacities_ETA_paranoaSul_upgrade_1,
                                                 ETA_paranoaSul_bonds_capacity_1, construction_time_interval, 12 * WEEKS_IN_YEAR);

    SequentialJointTreatmentExpansion etapa2_ETA_paranoaSul("Etapa 2 da ETA Paranoa Sul", 7, 4, 1, {6, 7, 8}, capacities_ETA_paranoaSul_upgrade_2,
                                                  ETA_paranoaSul_bonds_capacity_2, construction_time_interval, 12 * WEEKS_IN_YEAR);

    SequentialJointTreatmentExpansion etapa3_ETA_paranoaSul("Etapa 3 da ETA Paranoa Sul", 8, 4, 2, {6, 7, 8}, capacities_ETA_paranoaSul_upgrade_3,
                                                ETA_paranoaSul_bonds_capacity_3, construction_time_interval, 12 * WEEKS_IN_YEAR);

    /// Bonds da expansão da ETA Lago Norte (Sistema Paranoá)
    double cost_ETA_lagoNorte_upgrade = 30;

    vector<Bond *> bonds_ETA_lagoNorte_upgrade; //empréstimo para ampliação 1 da ETA (menor capacidade)
    uid = 0;
    for (double &cost : cost_ETA_lagoNorte_upgrade) {
        bonds_ETA_lagoNorte_upgrade.emplace_back(new LevelDebtServiceBond(9 + uid, cost, 25, 0.05, vector<int>(1, 0)));
        uid++;
    }

    /// Expansão da ETA Lago Norte (Sistema Paranoá)
    double capacity_ETA_lagoNorte_upgrade = 0.35; //aumento da capacidade de produção da ETA Lago Norte = + 0.35 m³/s

    SequentialJointTreatmentExpansion ETA_lagoNorte_upgrade("Expansao da ETA Lago Norte", 9, 4, 0, {4, 9}, //4 é a ID do Lago Paranoá
                                                      capacity_ETA_lagoNorte_upgrade, bonds_ETA_lagoNorte_upgrade,
                                                      construction_time_interval, 12 * WEEKS_IN_YEAR); //previsão: 2034 a 2037


    /*LevelDebtServiceBond dummy_bond(13, 1., 1, 1., vector<int>(1, 0)); // O QUE É DUMMY ENDPOINT?
    Reservoir dummy_endpoint("Dummy Node", 13, vector<Catchment *>(), 1., 0, evaporation_durham, 1,
                             construction_time_interval, 0, dummy_bond);*/

    vector<WaterSource *> water_sources; //water_sources é um vetor comum, que comportará todas as opções descritas acima de ampliação da infraestrutura de abastecimento
    water_sources.push_back(&descoberto);
    water_sources.push_back(&tortoSM);
    water_sources.push_back(&corumba);
    water_sources.push_back(&paranoa);
    water_sources.push_back(&ribeirao_bananal);
    //water_sources.push_back(&subsistema_gama);

    water_sources.push_back(&descoberto_expansion);
    water_sources.push_back(&ETA_corumba_etapa2);
    water_sources.push_back(&ETA_corumba_etapa3);
    water_sources.push_back(&etapa1_ETA_paranoaSul);
    water_sources.push_back(&etapa2_ETA_paranoaSul);
    water_sources.push_back(&etapa3_ETA_paranoaSul);
    water_sources.push_back(&ETA_lagoNorte_upgrade);

    //water_sources.push_back(&dummy_endpoint);

    /*
     * System connection diagram (water
     * flows from top to bottom)
     * Potential projects and expansions
     * of existing sources in parentheses
     *
     *                         1
     *      0(10)               \
     *       \                   \
     *        \                   3
     *         \                   \______4(6, 7, 8, 9)
     *          \                     /
     *           \                   /
     *            \                 /
     *             \               /
     *              \             /
     *               \           /
     *                \         /
     *                 \       /
     *                  \     /
     *                   \   /
     *                     5
     *                     /
     *                    /
     *                   /
     *                  /
     *                 /
     *                /
     *         ------
     *        /
     *   2(11, 12)
     *
     */

    Graph g(6); //graph é uma forma de estruturar dados que consiste em dois componentes: um conjunto finito de vértices denominado de nós; e um cojunto finito de pares ordenados (x,y), denominados de edges.
    g.addEdge(0, 5); //essa conexão indica que existe uma aresta do vértice 0 ao vértice 5.
    g.addEdge(4, 5);
    g.addEdge(1, 3);
    g.addEdge(3, 4);
    g.addEdge(5, 2);

    auto demand_n_weeks = (int) std::round(40 * WEEKS_IN_YEAR); //40 é o número de anos? A fç auto serve para declarar variáveis cujo tipo vai ser inferido pelo compilador a partir da inicialização delas

   // Criação dos vetores de vazão de retorno, relacionada aos efluentes lançados em corpos hídricos utilizados também para abastecimento urbano.

    vector<int> caesb_ws_return_id = {4}; //ID do reservatório onde o efluente será lançado (montante desse reservatório) - 4 (ID do Paranoá)
    WwtpDischargeRule wwtp_discharge_caesb( //em demand_to_waste_water_fraction, deve-se colocar a % de agua tratada (que sai da ETA - demanda) que vai voltar como vazão efluente para cada corpo receptor
            demand_to_wastewater_fraction_caesb, // n° de séries (de demanda de vazão efluente) correspondentes ao n° de reservatórios onde os efluentes são lançados. Obs: se no meu estudo de caso, apenas o Paranoá for receptor, esse arquivo csv terá uma linha e 53 colunas (vazão efluente de cada semana)
            caesb_ws_return_id); // id de cada reservatório/corpo receptor de efluente

    Utility caesb_descoberto((char *) "CAESB Descoberto", 0, demand_caesb_descoberto, demand_n_weeks, caesb_annual_payment,//criação das companhias de água. A descrição de cada termo está no arquivo .doc.
                  caesbDemandClassesFractions, caesbUserClassesWaterPrices, wwtp_discharge_caesb,
                  caesb_inf_buffer, rof_triggered_infra_order_caesb_descoberto,
                  vector<int>(), rofs_infra_caesb_descoberto, discount_rate, bond_term[0], bond_rate[0]);

    Utility caesb_tortoSM((char *) "CAESB Torto/Santa Maria", 1, demand_caesb_tortoSM, demand_n_weeks, caesb_annual_payment,//criação das companhias de água. A descrição de cada termo está no arquivo .doc.
                             caesbDemandClassesFractions, caesbUserClassesWaterPrices, wwtp_discharge_caesb,
                             caesb_inf_buffer, rof_triggered_infra_order_caesb_tortoSM,
                             vector<int>(), rofs_infra_caesb_tortoSM, discount_rate, bond_term[0], bond_rate[0]);


    vector<Utility *> utilities; //vetor que junta as companhias criadas acima
    utilities.push_back(&caesb_descoberto);
    utilities.push_back(&caesb_tortoSM);

    /// Water-source-utility connectivity matrix (each row corresponds to a utility and numbers are water
    /// sources IDs.
    vector<vector<int>> reservoir_utility_connectivity_matrix = {
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}  //CAESB
    };
    // O que table_storage_shift representa? O que são esses números (2000, 5000...) [3] [17]
    auto table_storage_shift = vector<vector<double>>(4, vector<double>(25, 0.));
    table_storage_shift[3][17] = 2000.; //tem a ver com a RdF
    table_storage_shift[3][8] = 5000.;
    table_storage_shift[1][14] = 100.;
    table_storage_shift[1][20] = 500.;
    table_storage_shift[1][21] = 500.;
    table_storage_shift[1][15] = 700.;
    table_storage_shift[1][9] = 700.;

    vector<DroughtMitigationPolicy *> drought_mitigation_policies; //criação do vetor das políticas de mitigação de seca (racionamento + transferência)
    /// Restriction policies
    vector<double> initial_restriction_triggers = {caesb_descoberto_restriction_trigger,
                                                   caesb_tortoSM_restriction_trigger,
                                                   caesb_paranoa_restriction_trigger,
                                                   caesb_corumbaIV_restriction_trigger}; //variáveis de decisão que representam os gatilhos para acionar o racionamento em cada companhia

    vector<double> restriction_stage_multipliers_caesb_descoberto = {0.9, 0.8, 0.7, 0.6}; //São 4 estágios de racionamento. Os fatores 0.9, 0.8, 0.7, 0.6 são as restrições da demanda. 0.9 significa que a demanda será restringida em 10% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_descoberto = {initial_restriction_triggers[0], //estágio 1
                                                                    initial_restriction_triggers[0] + 0.15f, // estágio 2 -> 0.15f está relacionado ao limite da métrica de risco utilizado para acionar a implementação de racionamento do estágio 2 (mais severo do que o estágio 1)
                                                                    initial_restriction_triggers[0] + 0.35f, //estágio 3
                                                                    initial_restriction_triggers[0] + 0.6f}; //estágio 4

    vector<double> restriction_stage_multipliers_caesb_tortoSM = {0.9, 0.8, 0.7, 0.6}; //São 4 estágios de racionamento. Os fatores 0.9, 0.8, 0.7, 0.6 são as restrições da demanda. 0.9 significa que a demanda será restringida em 10% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_tortoSM = {initial_restriction_triggers[0], //estágio 1
                                                                  initial_restriction_triggers[0] + 0.15f, // estágio 2 -> 0.15f está relacionado ao limite da métrica de risco utilizado para acionar a implementação de racionamento do estágio 2 (mais severo do que o estágio 1)
                                                                  initial_restriction_triggers[0] + 0.35f, //estágio 3
                                                                  initial_restriction_triggers[0] + 0.6f}; //estágio 4

    vector<double> restriction_stage_multipliers_caesb_paranoa = {0.9, 0.8, 0.7, 0.6}; //São 4 estágios de racionamento. Os fatores 0.9, 0.8, 0.7, 0.6 são as restrições da demanda. 0.9 significa que a demanda será restringida em 10% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_paranoa = {initial_restriction_triggers[0], //estágio 1
                                                                  initial_restriction_triggers[0] + 0.15f, // estágio 2 -> 0.15f está relacionado ao limite da métrica de risco utilizado para acionar a implementação de racionamento do estágio 2 (mais severo do que o estágio 1)
                                                                  initial_restriction_triggers[0] + 0.35f, //estágio 3
                                                                  initial_restriction_triggers[0] + 0.6f}; //estágio 4

    vector<double> restriction_stage_multipliers_caesb_corumbaIV = {0.9, 0.8, 0.7, 0.6}; //São 4 estágios de racionamento. Os fatores 0.9, 0.8, 0.7, 0.6 são as restrições da demanda. 0.9 significa que a demanda será restringida em 10% e assim por diante.
    vector<double> restriction_stage_triggers_caesb_corumbaIV = {initial_restriction_triggers[0], //estágio 1
                                                                  initial_restriction_triggers[0] + 0.15f, // estágio 2 -> 0.15f está relacionado ao limite da métrica de risco utilizado para acionar a implementação de racionamento do estágio 2 (mais severo do que o estágio 1)
                                                                  initial_restriction_triggers[0] + 0.35f, //estágio 3
                                                                  initial_restriction_triggers[0] + 0.6f}; //estágio 4

    Restrictions restrictions_d(0, // Criação das políticas de restrição de uso da água (puxa o código Restrictions.h, que contém os comandos dessa política)
                              restriction_stage_multipliers_caesb_descoberto,
                              restriction_stage_triggers_caesb_descoberto,
                              &caesbDemandClassesFractions,
                              &caesbUserClassesWaterPrices,
                              &caesbPriceRestrictionMultipliers);

    Restrictions restrictions_t(1, // Criação das políticas de restrição de uso da água (puxa o código Restrictions.h, que contém os comandos dessa política)
                              restriction_stage_multipliers_caesb_tortoSM,
                              restriction_stage_triggers_caesb_tortoSM,
                              &caesbDemandClassesFractions,
                              &caesbUserClassesWaterPrices,
                              &caesbPriceRestrictionMultipliers);

    Restrictions restrictions_p(4, // Criação das políticas de restrição de uso da água (puxa o código Restrictions.h, que contém os comandos dessa política)
                              restriction_stage_multipliers_caesb_paranoa,
                              restriction_stage_triggers_caesb_paranoa,
                              &caesbDemandClassesFractions,
                              &caesbUserClassesWaterPrices,
                              &caesbPriceRestrictionMultipliers);

    Restrictions restrictions_c(2, // Criação das políticas de restrição de uso da água (puxa o código Restrictions.h, que contém os comandos dessa política)
                              restriction_stage_multipliers_caesb_corumbaIV,
                              restriction_stage_triggers_caesb_corumbaIV,
                              &caesbDemandClassesFractions,
                              &caesbUserClassesWaterPrices,
                              &caesbPriceRestrictionMultipliers);

    drought_mitigation_policies = {&restrictions_d, &restrictions_t, &restrictions_p, &restrictions_c};

    // TRANSFER POLICY

    int buyer_id = 0; //apenas descoberto solicitará água por enquanto

    double buyer_transfer_capacity = 0.7; //a capacidade de transferência de água do TortoSM para o Descoberto é de até 700 l/s
    double buyer_transfer_trigger = caesb_descoberto_transfer_trigger;

    Graph transfer_graph(2);
    transfer_graph.addEdge(1, 0); // Água do tortoSM para o Descoberto

    Transfers transfer_intersystem(0, 1, 1, ILLIMITED_TREATMENT_CAPACITY, buyer_id,
                                   buyer_transfer_capacity, buyer_transfer_trigger,
                                   transfer_graph, vector<double>(), vector<int>());

    drought_mitigation_policies = {&transfer_intersystem};

    /// Creates simulation object depending on use (or lack thereof) ROF tables
    double start_time = omp_get_wtime();
    if (import_export_rof_tables == EXPORT_ROF_TABLES) {
        s = new Simulation(water_sources,
                           g,
                           reservoir_utility_connectivity_matrix,
                           utilities,
                           drought_mitigation_policies,
                           min_env_flow_controls,
                           utilities_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                           water_sources_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                           policies_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                           n_weeks,
                           realizations_to_run,
                           rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, nullptr);
    } else if (import_export_rof_tables == IMPORT_ROF_TABLES) {
        s = new Simulation (water_sources,
                            g,
                            reservoir_utility_connectivity_matrix,
                            utilities,
                            drought_mitigation_policies,
                            min_env_flow_controls,
                            utilities_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                            water_sources_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                            policies_rdm, //comentar essa linha (até eu gerar esse arquivo - se refere aos fatores de DU)
                            n_weeks,
                            realizations_to_run,
                            rof_tables,
                            table_storage_shift,
                            rof_tables_directory);
        //realization_start = omp_get_wtime();
        this->master_data_collector = s->runFullSimulation(n_threads, nullptr);
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
        this->master_data_collector = s->runFullSimulation(n_threads, nullptr);
    }
    double end_time = omp_get_wtime();
//	printf("Function evaluation time: %f s\n", end_time - start_time);

    //double realization_end = omp_get_wtime();
    //std::cout << "Simulation took  " << realization_end - realization_start
    //      << "s" << std::endl;

    /// Calculate objectives and store them in Borg decision variables array.
#ifdef  PARALLEL
    objectives = calculateAndPrintObjectives(false);

        int i = 0;
        objs[i] = min(min(objectives[i], objectives[5 + i]),
        		   min(objectives[10 + i], objectives[15 + i])) - 1.;
        for (i = 1; i < 5; ++i) {                       //são 5 objetivos que serão otimizados (0 a 4)
            objs[i] = max(max(objectives[i], objectives[5 + i]),
      	                  max(objectives[10 + i], objectives[15 + i]));
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

int Caesb::simulationExceptionHander(const std::exception &e, Simulation *s, // :: significa "resolução de escopo"
                                        double *objs, const double *vars) {
    int num_dec_var = 11; //número de variáveis desse estudo de caso - alterar para o valor do meu estudo
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


void Caesb::readInputData() { //A partir dessa linha serão inseridos os dados de entrada para que o modelo possa funcionar
    cout << "Reading input data." << endl;
    string data_dir = DEFAULT_DATA_DIR + BAR;

#pragma omp parallel num_threads(omp_get_thread_num())
    {
#pragma omp single
        streamflows_descoberto = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "descoberto_inflows.csv", n_realizations);
/*#pragma omp single
        streamflows_subsistema_gama = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "subsistema_gama_inflows.csv", n_realizations);*/
#pragma omp single
        streamflows_tortoSM = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "tortoSM_inflows.csv", n_realizations);
#pragma omp single
        streamflows_bananal = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "subsistema_bananal_inflows.csv", n_realizations);
        // }
#pragma omp single
        streamflows_paranoa = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "paranoa_inflows.csv", n_realizations);

#pragma omp single
        streamflows_corumbaIV = Utils::parse2DCsvFile( //inserção dos dados de vazão
                io_directory + DEFAULT_DATA_DIR + "inflows" + evap_inflows_suffix +
                BAR + "corumbaIV_inflows.csv", n_realizations);

// };
        //cout << "Reading evaporations." << endl;
#pragma omp single
        evap_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "descoberto_evap.csv", n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "tortoSM_evap.csv", n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_paranoa = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "paranoa_evap.csv", n_realizations); //inserção dos dados de evaporação
#pragma omp single
        evap_corumba = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "evaporation" + evap_inflows_suffix +
                BAR + "corumba_evap.csv", n_realizations); //inserção dos dados de evaporação


        //cout << "Reading demands." << endl;
#pragma omp single
        demand_caesb_descoberto = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_descoberto_demand.csv", n_realizations);

        demand_caesb_tortoSM = Utils::parse2DCsvFile(
                io_directory + DEFAULT_DATA_DIR + "demands" + evap_inflows_suffix + //inserção dos dados de demanda da caesb
                BAR + "caesb_tortoSM_demand.csv", n_realizations);


        //cout << "Reading others." << endl;
#pragma omp single
        {
            demand_to_wastewater_fraction_caesb = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "demand_to_wastewater_fraction_caesb.csv"); //demanda de efluentes da caesb

            caesbDemandClassesFractions = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "caesbDemandClassesFractions.csv"); //demanda de cada categoria de usuário da caesb

            caesbUserClassesWaterPrices = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "caesbUserClassesWaterPrices.csv"); //tarifa de água para cada categoria de usuário da caesb
                                                                                                   //OBS: inserir no valor a tarifa de esgoto também na última coluna desse arquivo (100% da de água)
            caesbPriceRestrictionMultipliers = Utils::parse2DCsvFile(
                    io_directory + DEFAULT_DATA_DIR + "caesbPriceRestrictionMultipliers.csv"); //% de aumento da tarifa para cada categoria durante o racionamento
        }
//    cout << "Done reading input data." << endl;
    }

}