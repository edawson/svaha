#include <iostream>
#include <set>
#include <functional>
#include <algorithm>
#include <getopt.h>
#include <vector>
#include <unordered_map>
#include <zlib.h>
#include <fstream>
#include <omp.h>

#include "rodeo.hpp"
#include "wrangler.hpp"
#include "kseq.hpp"
#include "gfakluge.hpp"
//#include "Variant.h"
#include "vcfparse.hpp"
#include <map>


using namespace std;
using namespace gfak;
using namespace rodeo;
namespace rodeo{
    KSEQ_INIT(gzFile, gzread);

    struct Node{
        int id; // Must be unique within a contig given how we're going to
        // abuse the id fields later.
        std::string sequence;
        std::string path;
        vector<Node*> next;
        vector<Node*> prev;
    };

    struct Edge{
        bool is_reverse;
        Node* from;
        Node* to;
    };

    /**
      void parse_variants(vector<string>& v_files,
      vector<vcflib::Variant>& variants){
      for (auto v : v_files){
      vcflib::VariantCallFile variant_file;
      variant_file.open(v);
      if (!variant_file.is_open()) {
      cerr << "Error: couldn't open vcf file " << v << endl;
      exit(1);
      }




      }
      }



      void make_graph(vector<vcflib::Variant>& vars,
      map<string, char*>& seqs,
      map<string, int>& lengths){

      }
      */
}


int main(int argc, char** argv){
    vector<char*> ref_files;
    vector<string> var_files;
    int threads;
    int max_node_size = 1000;


    int c;
    int optind = 2;

    if (argc <= 2){
        cerr << "Usage: " << argv[0] << " [options] -r <REF.fa> -v <VAR.vcf>" << endl;
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {"max-node-size", required_argument, 0, 'm'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hr:v:t:m:", long_options, &option_index);
        if (c == -1){
            break;
        }

        switch (c){
            case 't':
                threads = atoi(optarg);
                break;
            case 'r':
                ref_files.push_back(optarg);
                break;
            case 'v':
                var_files.push_back(optarg);
                break;
            case 'm':
                max_node_size = atoi(optarg);
                break;
            default:
                abort();
        }
    }

    if (ref_files.size() != 1){
        cerr << "One and exactly one reference must be provided." << endl;
        exit(1);
    }

    if (var_files.size() != 1){
        cerr << "One and exactly one variant call file must be provided." << endl;
        exit(1);
    }

    map<string, char*> name_to_seq;
    map<string, int> name_to_length;
    map<string, vector<vcfparse::Variant> > contig_to_variants;
    map<string, map<int, vector<vcfparse::Variant> > > contig_to_breakpoint_to_variants;
    map<string, vector<int> > contig_to_breakpoints;

    std::function<int(int)> pos_to_node = [&](int i){
        return 1;
    };

    std::function<pair<int, int>(vcfparse::Variant)> var_to_bp = [&](vcfparse::Variant var){
        int front = 0;
        int back = 0;
        if ( !var.info["SVTYPE"].empty() && var.info["SVLEN"].empty()){
            cerr << "No end info. Skipping sv." << endl;
            return make_pair(-1, -1);
        }
        if (var.info["SVTYPE"] == "DEL"){
            front = var.pos;
            back = var.pos + (uint32_t) stoi(var.info["SVLEN"]);
        }
        else if (var.info["SVTYPE"] == "INS"){
            front = var.pos;
            back = var.pos + (uint32_t) stoi(var.info["SVLEN"]);
        }
        else if (var.info["SVTYPE"] == "DUP"){
            front = var.pos;
            back = var.pos + (uint32_t) stoi(var.info["SVLEN"]);

        }
        else if (var.info["SVTYPE"] == "INV"){
            front = var.pos;
            back = var.pos + (uint32_t) stoi(var.info["SVLEN"]);
        }
        else{
            front = var.pos;
            back = var.pos + 1;
        }
        return std::make_pair(front, back);
    };

    // parse fasta reference files using kseq
    vector<vcfparse::Variant> variants;
    vcfparse::Variant var;
    vector<vcfparse::Variant> insertions;

    std::ifstream infile(var_files[0]);
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line[0] != '#'){
            vcfparse::Variant vv = vcfparse::parse_line(line);
            variants.push_back(vv);
            pair<int, int> brks = var_to_bp(vv);
            if (brks.first != -1){
                contig_to_variants[vv.seq].push_back(vv);
                contig_to_breakpoints[vv.seq].push_back(brks.first);
                contig_to_breakpoints[vv.seq].push_back(brks.second);
                contig_to_breakpoint_to_variants[vv.seq][brks.first].push_back(vv);
                contig_to_breakpoint_to_variants[vv.seq][brks.second].push_back(vv);
            }
            if (vv.info["SVTYPE"].empty()){
                insertions.push_back(vv);
            }
            else if (vv.info["SVTYPE"] == "INS"){
                insertions.push_back(vv);
            }
        }
    }

    map<int64_t, Node*> bp_to_node;


    for (auto f : ref_files){
        gzFile fp;
        kseq_t *seq;
        int l;
        fp = gzopen(f, "r");
        seq = kseq_init(fp);
        uint64_t currID = 0;
        while ((l = kseq_read(seq)) >= 0) {
            //char * x = new char[seq->seq.l];
            //lengths[seq->name.s] = seq->seq.l;
            //

            GFAKluge og;
            og.set_version();

            string contig = seq->name.s;
            char* cseq = seq->seq.s;
            int clen = seq->seq.l;

            vector<Node*> contig_nodes;
            vector<Edge*> contig_edges;
            int curr_id = 0;
            vector<int64_t> ins_nodes;
            if (contig_to_variants.find(contig) != contig_to_variants.end()){

                // our current id
                int c_id = 0;

                // our current position, used for substring generation.
                int c_pos = 0;

                // maps from a breakpoint to the node right before that breakpoint.
                map<int, Node*> bp_to_node;

                // insert breakpoints
                set<int> s_bps = set<int>(contig_to_breakpoints[contig].begin(), contig_to_breakpoints[contig].end());
                for (int i = max_node_size; i < clen; i+=max_node_size){
                    s_bps.insert(i);
                }
                contig_to_breakpoints[contig] = vector<int>(s_bps.begin(), s_bps.end());
                std::sort(contig_to_breakpoints[contig].begin(), contig_to_breakpoints[contig].end());
                for (int i = 0; i < contig_to_breakpoints[contig].size(); ++i){

                    // offset due to inserted sequence
                    int ins_offset = 0;

                    Node* n = new Node();
                    int brk = contig_to_breakpoints[contig][i]; 
                    n->id = ++c_id;        
                    //cerr << n-> id << endl;
                    //cerr << brk << "\t" << c_pos << endl;
                    n->sequence.assign(cseq + c_pos, brk - c_pos);
                    n->path = contig;

                    if (i > 0){
                        n->prev.push_back(contig_nodes.back());


                        // SNPS / insertions are off the ref path
                        if (contig_nodes.back()->path.empty()){
                            int prev_ref = contig_nodes.size() - 1;
                            while (contig_nodes[prev_ref]->path.empty()){
                                --prev_ref; 

                                Edge* secondary_edge = new Edge();
                                secondary_edge->from = contig_nodes[prev_ref];
                                secondary_edge->to = n;
                                contig_edges.push_back(secondary_edge);
                            }
                        }

                        Edge* e = new Edge();
                        e->from = contig_nodes.back();
                        e->to = n;
                        e->is_reverse = false;
                        contig_edges.push_back(e);

                        (contig_nodes.back())->next.push_back(n);
                    }

                    contig_nodes.push_back(n);
                    bp_to_node[ contig_to_breakpoints[contig][i] ] = n;

                    // Add in extra nodes if needed.
                    vector<vcfparse::Variant> bp_vars = contig_to_breakpoint_to_variants[contig][brk];
                    for (auto bpv : bp_vars){
                        if (bpv.info["SVTYPE"].empty()){
                            for (int altp = 0; altp < bpv.alts.size(); ++altp){
                                Node* snp_n = new Node();
                                snp_n->sequence = string(bpv.alts[altp]);
                                //cerr << bpv.alts[altp] << endl;
                                snp_n->id = ++c_id;
                                n->next.push_back(snp_n);
                                contig_nodes.push_back(snp_n);

                                snp_n->prev.push_back(n);

                                Edge* ie = new Edge();
                                ie->from = n;
                                ie->to = snp_n;
                                ie->is_reverse = false;
                            }
                            ins_offset = max(ins_offset, 1);
                        }
                        else if (bpv.info["SVTYPE"] == "INS"){
                            Node* ins_n = new Node();
                            ins_n->sequence.assign(bpv.info["SVSEQ"]);
                            ins_n->id = ++c_id;
                            n->next.push_back(ins_n);

                            contig_nodes.push_back(ins_n);
                            ins_n->prev.push_back(n);
                            Edge* ie = new Edge();
                            ie->from = n;
                            ie->to = ins_n;
                            ie->is_reverse = false;
                            ins_offset = max(ins_offset, (int) ins_n->sequence.length());
                        }
                        else{
                            continue;
                        }
                    }

                    c_pos = contig_to_breakpoints[contig][i];
                }

                // finally, wire up the edges for DEL / INV
                for (auto vvar : contig_to_variants[contig]){

                    Edge* e_from = new Edge();
                    Edge* e_to = new Edge();
                    if (vvar.info["SVTYPE"] == "DEL"){
                        e_from->from = bp_to_node[ vvar.pos ];
                        e_from->to = bp_to_node [ vvar.pos + stoi(vvar.info["SVLEN"]) ];
                    }
                    else if (var.info["SVTYPE"] == "INS"){

                    }
                    else if (var.info["SVTYPE"] == "DUP"){

                    }
                    else if (var.info["SVTYPE"] == "INV"){
                        e_from->from = bp_to_node[ vvar.pos ];
                        e_from->to = bp_to_node [ vvar.pos + stoi(vvar.info["SVLEN"]) ];
                        e_from->is_reverse = true;

                        e_to->from = bp_to_node [ vvar.pos ];
                        e_to->to = bp_to_node [ vvar.pos + stoi(vvar.info["SVLEN"]) ];
                        e_to->is_reverse = true;
                    }
                    else{

                    }

                }


                // Label paths for contig
                //   and variant ID

                // format GFA and output it.
                // Delete all our business as we go
                int prank = 0;
                for (auto cnode : contig_nodes){
                    sequence_elem s;
                    s.sequence = cnode->sequence;
                    s.name = to_string(cnode->id);
                    og.add_sequence(s);
                    if (!cnode->path.empty()){
                    path_elem p;
                    p.source_name = s.name;
                    p.name = cnode->path;
                    p.is_reverse = false;
                    p.rank = ++prank;
                    p.cigar = to_string(s.sequence.length()) + "M";
                    og.add_path(s.name, p);
                    }
                }
                for (auto cedge : contig_edges){
                    link_elem l;
                    l.source_name = to_string(cedge->from->id);
                    l.sink_name = to_string(cedge->to->id);
                    l.source_orientation_forward = !cedge->is_reverse;
                    l.sink_orientation_forward = !cedge->is_reverse;
                    l.cigar = "";
                    og.add_link(l.source_name, l);
                    delete cedge;
                }
                for (auto cnode : contig_nodes){
                    delete cnode;
                }

            }
            cout << og.to_string();
        }
        gzclose(fp);



        return 0;
    }
}
