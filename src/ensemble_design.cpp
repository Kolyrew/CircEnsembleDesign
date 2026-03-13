#include <iomanip>
#include "optimizer.cpp"
#include "Utils/reader.h"


#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif

using namespace EnsembleDesign;

void prepend_prefix_to_dfa(WDFA& dfa, const string& prefix) {
    int prefix_len = prefix.length();
    if (prefix_len == 0) return;

    WDFA new_dfa;
    int orig_nodes_count = dfa.nodes.size();
    int new_nodes_count = orig_nodes_count + prefix_len;
    int new_edges_count = dfa.edges.size() + prefix_len + 10;
    new_dfa.reserve_capacity(new_nodes_count, new_edges_count);

    for (int k = 0; k < prefix_len; k++) {
        NodeType node = make_pair(static_cast<IndexType>(k), static_cast<NumType>(0));
        new_dfa.add_node(node);
    }

    for (int k = 0; k < prefix_len; k++) {
        NodeType n1 = make_pair(static_cast<IndexType>(k), static_cast<NumType>(0));
        NodeType n2 = make_pair(static_cast<IndexType>(k + 1), static_cast<NumType>(0));
        NucType nuc = GET_ACGU_NUC(prefix[k]);
        new_dfa.add_edge(n1, n2, nuc, 0.0, util::value_min<ScoreType>(), 0.0);
    }

    for (int pos = 0; pos < orig_nodes_count; pos++) {
        for (auto& node : dfa.nodes[pos]) {
            NodeType new_node = make_pair(
                static_cast<IndexType>(node.first + prefix_len),
                node.second);
            new_dfa.add_node(new_node);
        }
    }

    for (auto& edge : dfa.edges) {
        NodeType old_n1 = get<0>(edge);
        NodeType old_n2 = get<1>(edge);
        NucType nuc = get<2>(edge);
        Parameter* param = get<3>(edge);

        NodeType new_n1 = make_pair(
            static_cast<IndexType>(old_n1.first + prefix_len), old_n1.second);
        NodeType new_n2 = make_pair(
            static_cast<IndexType>(old_n2.first + prefix_len), old_n2.second);

        new_dfa.add_edge(new_n1, new_n2, nuc, param->weight, param->gradient, param->cai_score);
    }

    dfa = std::move(new_dfa);
}

void append_suffix_to_dfa(WDFA& dfa, const string& suffix) {
    int suffix_len = suffix.length();
    if (suffix_len == 0) return;

    WDFA new_dfa;
    int orig_nodes_count = dfa.nodes.size();
    int end_pos = orig_nodes_count - 1;
    int new_nodes_count = orig_nodes_count + suffix_len;
    int extra_edges = suffix_len + static_cast<int>(dfa.nodes[end_pos].size()) + 20;
    int new_edges_count = dfa.edges.size() + extra_edges;
    new_dfa.reserve_capacity(new_nodes_count, new_edges_count);

    for (int pos = 0; pos < orig_nodes_count; pos++) {
        for (auto& node : dfa.nodes[pos]) {
            new_dfa.add_node(node);
        }
    }

    for (auto& edge : dfa.edges) {
        NodeType n1 = get<0>(edge);
        NodeType n2 = get<1>(edge);
        NucType nuc = get<2>(edge);
        Parameter* param = get<3>(edge);
        new_dfa.add_edge(n1, n2, nuc, param->weight, param->gradient, param->cai_score);
    }

    for (int k = 1; k <= suffix_len; k++) {
        NodeType node = make_pair(static_cast<IndexType>(end_pos + k), static_cast<NumType>(0));
        new_dfa.add_node(node);
    }

    NucType first_nuc = GET_ACGU_NUC(suffix[0]);
    for (const auto& node_at_end : dfa.nodes[end_pos]) {
        NodeType n2 = make_pair(static_cast<IndexType>(end_pos + 1), static_cast<NumType>(0));
        new_dfa.add_edge(node_at_end, n2, first_nuc, 0.0, util::value_min<ScoreType>(), 0.0);
    }
    for (int k = 1; k < suffix_len; k++) {
        NodeType n1 = make_pair(static_cast<IndexType>(end_pos + k), static_cast<NumType>(0));
        NodeType n2 = make_pair(static_cast<IndexType>(end_pos + k + 1), static_cast<NumType>(0));
        NucType nuc = GET_ACGU_NUC(suffix[k]);
        new_dfa.add_edge(n1, n2, nuc, 0.0, util::value_min<ScoreType>(), 0.0);
    }

    dfa = std::move(new_dfa);
}

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
    if (has_named_args && argc > 6 && strncmp(argv[6], "--", 2) != 0) {
        init_solution = argv[6];
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
        if (!prefix.empty()) {
            prepend_prefix_to_dfa(dfa, prefix);
        }
        parser.optimize(dfa, codon, aa_seq, protein,
                       aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit,
                       aa_graphs_with_ln_weights);


    }/**/
    return 0;
}
