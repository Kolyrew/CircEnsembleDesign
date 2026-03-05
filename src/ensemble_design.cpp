#include <iomanip>
#include "optimizer.cpp"
#include "Utils/reader.h"


#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

using namespace EnsembleDesign;

int main(int argc, char** argv) {

    int beam_size = 0;
    int num_epochs = 30;
    double learning_rate = 0.03;
    double epsilon = 1.0;
    bool is_verbose = true;
    unsigned int rand_seed = 0;
    string init_solution = "";
    string prefix = "";
    string suffix = "";
    int ires_end = 0;
    double ires_orf_lambda = 1.0;
    double cross_pair_prob_threshold = 0.01;
    int max_cross_pairs = 5;

    string CODON_TABLE = "./codon_usage_freq_table_human.csv";

    if (argc > 1) beam_size = atoi(argv[1]);
    if (argc > 2) num_epochs = atoi(argv[2]);
    if (argc > 3) learning_rate = atof(argv[3]);
    if (argc > 4) epsilon = atof(argv[4]);
    if (argc > 5) rand_seed = atoi(argv[5]);
    if (argc > 6) init_solution = argv[6];

    bool has_named_args = false;
    for (int i = 7; i < argc; i++) {
        if (strncmp(argv[i], "--", 2) == 0) { has_named_args = true; break; }
    }

    if (has_named_args) {
        for (int i = 7; i < argc; i++) {
            string arg = argv[i];
            if (arg == "--prefix" && i+1 < argc) { prefix = argv[++i]; }
            else if (arg == "--suffix" && i+1 < argc) { suffix = argv[++i]; }
            else if (arg == "--iresEND" && i+1 < argc) { ires_end = atoi(argv[++i]); }
            else if (arg == "--lambda" && i+1 < argc) { ires_orf_lambda = atof(argv[++i]); }
            else if (arg == "--cross_pair_prob_threshold" && i+1 < argc) { cross_pair_prob_threshold = atof(argv[++i]); }
            else if (arg == "--max_cross_pairs" && i+1 < argc) { max_cross_pairs = atoi(argv[++i]); }
        }
    } else {
        if (argc > 7) prefix = argv[7];
        if (argc > 8) ires_orf_lambda = atof(argv[8]);
        if (argc > 9) cross_pair_prob_threshold = atof(argv[9]);
        if (argc > 10) max_cross_pairs = atoi(argv[10]);
        ires_end = prefix.length();
    }

    Codon codon(CODON_TABLE);
    std::unordered_map<std::string, Lattice> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, 0);


    // main loop
    string aa_seq, aa_tri_seq;
    vector<string> aa_seq_list, aa_name_list;
    // load input
    for (string seq; getline(cin, seq);){
        if (seq.empty()) continue;
        if (seq[0] == '>'){
            aa_name_list.push_back(seq); // sequence name
            if (!aa_seq.empty())
                aa_seq_list.push_back(aa_seq);
            aa_seq.clear();
            continue;
        }else{
            rtrim(seq);
            aa_seq += seq;
        }
    }
    if (!aa_seq.empty())
        aa_seq_list.push_back(aa_seq);


    // start design
    for(int i = 0; i < aa_seq_list.size(); i++){
        if (aa_name_list.size() > i)
            cout << aa_name_list[i] << endl;
        auto& aa_seq = aa_seq_list[i];
        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);
        aa_tri_seq.clear();

        if (!ReaderTraits<Fasta>::cvt_to_seq(aa_seq, aa_tri_seq)) 
            continue;

        Optimizer parser(beam_size, num_epochs, learning_rate, epsilon, init_solution, is_verbose, rand_seed, prefix, prefix.length(), ires_orf_lambda, cross_pair_prob_threshold, max_cross_pairs, suffix, ires_end);

        auto protein = util::split(aa_tri_seq, ' ');
        auto dfa = get_dfa(aa_graphs_with_ln_weights, util::split(aa_tri_seq, ' '));
        parser.optimize(dfa, codon, aa_seq, protein,
                       aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit,
                       aa_graphs_with_ln_weights);


    }/**/
    return 0;
}
