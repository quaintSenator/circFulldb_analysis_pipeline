//
// Created by tail on 2024/3/11.
//

#ifndef HUMAN_ALLBSJS_FINAL_JSON_METHODS_H
#define HUMAN_ALLBSJS_FINAL_JSON_METHODS_H
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <sys/types.h>
#include <sys/unistd.h>
#include <unistd.h>
//#include "vcpkg/packages/jsoncpp_x64-windows/include/json/json.h"
#define MILLION 1000000
class isoform_info{
public:
    int iso_len;
    char strand;
    isoform_info(int il, char str):iso_len(il),strand(str){}
    isoform_info(const isoform_info&)=default;
};
class bsj_map_copy{
public:
    std::string strand;
    std::string gene_id;
    std::string gene_name;
    int isoform_count;
    bsj_map_copy(std::string &strnd, std::string &gid, std::string &gname, int iso_count):
    strand(strnd),gene_id(gid),gene_name(gname),isoform_count(iso_count){}
    bsj_map_copy(){
        isoform_count = 0;
    }
};
void treat_humaninfo();
void mainFun4MouseBSJ_Treatment();
void get_humanBSJallInfo();
void get_humanIso();
void get_humanExpression();
void get_human_miRNAbinding();
void get_human_IRES();
void get_human_TIS();
void get_human_m6A();
void long_isoID_encoder();
void all_seq_split_into_5k();
void split_then_send2miRanda();
void IRES_result_preprocess();
void batch_RBPmap_launcher();
void specific_expression_analysis();
int decideTissueCollection(int);
void get_most_specific_isos();
int decideChrFromID(std::string &s);
void get_RBPdata_from_all_batch(int batchNo);
void merge_Human_Mouse_BSJs();
void mouse_id_recoder(std::string&);
void judge_Human_Mouse_BSJs_has_collapse();
void merge_Human_Mouse_isos();
void merge_Human_Mouse_exps();
void count_Human_BSJs();
void count_Mouse_BSJs();
void get_Mouse_length_chromo_distribution_matlab_data();
void get_Mouse_read_chromo_distribution_matlab_data();
void get_ORFs();
void generate_mouse_bsjs_without_function();
void automatic_scp_and_process_RBP_batch();
void count_raw_RBP_batch_count();
void get_Mouse_IRES();
void get_Mouse_miRNA();
void get_Mouse_m6a_from_flat_text();
void reconstruct_RBP_report_from_all_batches();
void reconstruct_RBP_report_from_batch(std::vector<std::string> &article, std::map<int,std::string> &seqNoMap, std::ofstream &outputFile);
void correct_RBP_output();

void report_our_identification_results();
void get_Mouse_chromo_micro_chromosome_distribution_matlab_data();
std::string interval_style_ID_convert2_array_style(std::string &ID, bool containBSJID);
std::string interval_style_ID_convert2_array_style(std::string &ID, bool containBSJID, char strand);

int judgeORFvalid(int circLen, int l, int r);

std::string substr_by_be(std::string& str, int b, int e);
void split_line_into_strvec(char cutByChar, bool needQuote, std::vector<std::string> &lineStrVec, std::string &line);

void get_Mouse_chromo_RBP_chromosome_distribution_matlab_data();
#endif //HUMAN_ALLBSJS_FINAL_JSON_METHODS_H
