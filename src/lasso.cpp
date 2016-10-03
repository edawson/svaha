#include <iostream>
#include <functional>
#include <algorithm>
#include <getopt.h>
#include <vector>
#include <unordered_map>
#include <tuple>
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
    KSEQ_INIT(gzFile, gzread)

        void parse_fastas(vector<char*>& refs,
                map<string, char*>& seqs,
                map<string, int>& lengths){
            for (auto f : refs){
                gzFile fp;
                kseq_t *seq;
                int l;
                fp = gzopen(f, "r");
                seq = kseq_init(fp);
                while ((l = kseq_read(seq)) >= 0) {
                    char * x = new char[seq->seq.l];
                    memcpy(x, seq->seq.s, seq->seq.l);
                    seqs[seq->name.s] = x;
                    lengths[seq->name.s] = seq->seq.l;
                }
                gzclose(fp);
            }
        }

    struct Node{
        int id; // Must be unique within a contig given how we're going to
        // abuse the id fields later.
        std::string sequence;
        std::string path;
        vector<Node*> next;
        vector<Node*> prev;
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

void print_help(char** argv){
    cerr << "lasso: make structura variation graphs." << endl
        << "usage: " << argv[0] << " [options]" << endl
        << "-r/--reference  the reference genome to add SVs to." << endl
        << " -v/--vcf       a VCF file containing SVs." << endl
        << endl;
}


int main(int argc, char** argv){
    vector<char*> ref_files;
    vector<string> var_files;
    int threads;


    int c;
    int optind = 2;

    if (argc <= 2){
        cerr << "Usage: lasso -r <REF.fa> -v <VAR.vcf>" << endl;
        exit(1);
    }

    while (true){
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"reference", required_argument, 0, 'r'},
            {"vcf", required_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {0,0,0,0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hr:v:t:", long_options, &option_index);
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

    int global_max_id = 0;
    // parse fasta reference files using kseq
    parse_fastas(ref_files, name_to_seq, name_to_length); 


    //map<string, map<int, vector<pair<int, pair<bool, bool> > > > > contig_to_source_to_sink_tostart_fromend;
    map<string, vector<int> > contig_to_breakpoints;
    map<string, vector<vcfparse::Variant> > contig_to_variants;
    map<string, map<int, vector<tuple< int, bool, bool> > > > contig_to_source_to_sink_tostart_fromend;

    map<string, map<int, vector<string> > > contig_to_pos_to_seq;

    map<string, int> contig_max_id;


    std::ifstream infile(var_files[0]);
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line[0] != '#'){
            vcfparse::Variant v = vcfparse::parse_line(line);
            contig_to_variants[v.seq].push_back(v);
        }
    }
    /**
     * At this point we have a pseudo graph of Nodes,
     * one per base in the sequence.
     *
     * We also have a vector of VCF entries in a simple string format.
     *
     * Now, we need to loop over those variants.
     * We'll switch on SVTYPE, and for DEL / INV / (DUP??)
     * we will pull out the end and then length. We can check these against the
     * confidence intervals if needed.
     * *** WE WILL ALSO NEED THE POSITION FROM THE VCF RECORD ***
     * Then, if we have a deletion, we simply insert node (start + len) into
     * the next member of node start.
     *
     * if we have an inversion, we insert node (start + len - 1) into the NEXT member
     * if start AND we insrt node ( start + 1 ) into the PREV member of node ( start + len ).
     *
     */
    for (auto c_v : contig_to_variants){
        vector<vcfparse::Variant> variants = c_v.second;
        for (int i = 0; i < variants.size(); i++){
            vcfparse::Variant var = variants[i];
            std::string::size_type sz;

            if (!var.info["SVLEN"].empty()){
                if (var.info["SVTYPE"] == "DEL"){
                    // Add two breakpoints
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 );
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 + (unsigned int) stoi(var.info["SVLEN"], &sz) );

                    //Add a single edge to represent the deletion
                    contig_to_source_to_sink_tostart_fromend[var.seq][ var.pos -1 ].push_back(std::make_tuple(var.pos - 1 + stoi(var.info["SVLEN"], &sz), true, true));
                }
                else if (var.info["SVTYPE"] == "INV"){
                    // Add two breakpoints
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 );
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 + stoi(var.info["SVLEN"], &sz) );

                    // Add two edges
                    contig_to_source_to_sink_tostart_fromend[var.seq][ var.pos - 1 ].push_back(std::make_tuple( var.pos - 1, true, false));
                    contig_to_source_to_sink_tostart_fromend[var.seq][ var.pos - 1  + stoi(var.info["SVLEN"], &sz) ].push_back(std::make_tuple(var.pos - 1 + stoi(var.info["SVLEN"], &sz) ,false, true));
                }
                else if (var.info["SVTYPE"] == "INS"){
                    // Add two breakpoints
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 );
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 + stoi(var.info["SVLEN"], &sz) );

                    // Add a new node with the inserted sequence

                    // Add two edges

                }
                else if (var.info["SVTYPE"] == "DUP"){
                    // Add a single cyclic edge and two breakpoints
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 );
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 + stoi(var.info["SVLEN"], &sz) );
                    contig_to_source_to_sink_tostart_fromend[var.seq][var.pos - 1].push_back(std::make_tuple(var.pos -1, true, true));
                }
                else if (var.info["SVTYPE"] == "BND"){
                    // TODO I'm just not sure how to handle breakends.
                    // Their natural representation is, at its simplest, just "breakpoints" in the graph, but there
                    // should be corresponding edges, too.

                }
                else{
                    // it's a SNP?
                    // add two breakpoints, a single-base node and two edges
                    //contig_to_breakpoints[var.seq].push_back( var.pos - 1 - 1);
                    //contig_to_breakpoints[var.seq].push_back( var.pos - 1 + 1);
                    //for (auto i_alt : var.alts){
                    //    contig_to_pos_to_seq[var.seq][var.pos - 1].push_back(i_alt);
                    //}
                }
            }
        else{
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 );
                    contig_to_breakpoints[var.seq].push_back( var.pos - 1 + 1);
                    int count = 0;
                    for (auto i_alt : var.alts){
                        contig_to_pos_to_seq[var.seq][var.pos - 1].push_back(i_alt);
                        //contig_to_source_to_sink_tostart_fromend[var.seq][ var.pos - 1 ].push_back(std::make_tuple(var.pos - 1, true, true));S
                        ++count;
                    }

        }
        }
        contig_to_breakpoints[c_v.first].push_back(name_to_length[c_v.first]);
        if (!contig_to_breakpoints[c_v.first].empty()){
            contig_max_id[c_v.first] = contig_to_breakpoints[c_v.first].back();
        }
    }

    map<int, Node*> id_to_node;

    map<string, vector<Node*> >cont_to_nodes;


    GFAKluge gg;
    gg.set_version();


    //map<string, vector<Node*> > cont_nodes;
    map<string, map<int, vector<tuple<int, bool, bool> > > >::iterator it;
    for (it = contig_to_source_to_sink_tostart_fromend.begin(); it != contig_to_source_to_sink_tostart_fromend.end(); it++){
        char* seq = name_to_seq[it->first];
        vector<int> breakpoints = contig_to_breakpoints[it->first];
        std::sort(breakpoints.begin(), breakpoints.end());
        // gives us the ability to get the previous node, which we need for all types
        // of SV pretty much.
        vector<int> con;
        int n_start = 0;
        for (int i = 0; i < breakpoints.size(); i++){

            sequence_elem sq;
            sq.sequence.assign(seq + n_start, breakpoints[i] - n_start);
            sq.name = std::to_string(n_start + 1);
            gg.add_sequence(sq);

            vector<string> seqs = contig_to_pos_to_seq[it->first][n_start]; 
            if (seqs.size() > 0){

                for (int i_seq= 0; i_seq < seqs.size(); i_seq++){
                    sequence_elem insert_seq;
                    insert_seq.sequence.assign(seqs[i_seq]);
                    insert_seq.name = std::to_string(++contig_max_id[it->first]); //+ "." + std::to_string(char ( i_seq ));

                    gg.add_sequence(insert_seq);

                    link_elem pre_ll;
                    // TODO could lead to Out of Bounds error
                    cerr << breakpoints[ i ] << endl;
                    pre_ll.source_name = std::to_string( breakpoints[i - 2] + 1 ); //std::to_string(n_start + 1);
                    pre_ll.sink_name = insert_seq.name;
                    pre_ll.source_orientation_forward = true;
                    pre_ll.sink_orientation_forward = true;
                    pre_ll.cigar = "0M";

                    gg.add_link(pre_ll.source_name, pre_ll);
                    cerr << pre_ll.source_name << " " << pre_ll.sink_name << endl;

                    link_elem post_ll;
                    post_ll.source_name = insert_seq.name;
                    post_ll.sink_name = std::to_string( breakpoints[i] + 1);
                    post_ll.source_orientation_forward = true;
                    post_ll.sink_orientation_forward = true;
                    post_ll.cigar = "0M";

                    gg.add_link(post_ll.source_name, post_ll);

                }
            }

            if (i < breakpoints.size() - 1){
                link_elem ref_l;
                ref_l.source_name = std::to_string(n_start + 1);
                ref_l.sink_name = std::to_string(breakpoints[i] + 1);
                ref_l.source_orientation_forward = true;
                ref_l.sink_orientation_forward = true;
                ref_l.cigar = "0M";
                gg.add_link(ref_l.source_name, ref_l);
            }

            for (auto x : it->second[n_start]){
                link_elem ll;
                ll.source_name = std::to_string(con[con.size() - 1] + 1);
                ll.sink_name = std::to_string( std::get<0>(x) + 1);
                ll.source_orientation_forward = std::get<1>(x);
                ll.sink_orientation_forward = std::get<2>(x);
                ll.cigar = "0M";

                gg.add_link(ll.source_name, ll);

            }

            // Add the reference path, as these are reference nodes.
            path_elem pe;
            pe.name = it->first;
            pe.rank = i+1;
            pe.source_name = sq.name;
            pe.is_reverse = false;
            stringstream p_cig;
            p_cig << sq.sequence.length() << "M";
            pe.cigar = p_cig.str();
            gg.add_path(pe.source_name, pe);

            con.push_back(n_start);
            n_start = breakpoints[i];
            global_max_id = global_max_id > contig_max_id[it->first] ? global_max_id : contig_max_id[it->first];
        }


        vector<vcfparse::Variant> variants = contig_to_variants[it->first];
        for (int i = 0; i < variants.size(); i++){
            vcfparse::Variant var = variants[i];
            std::string::size_type sz;
            if (var.info["SVTYPE"] == "DEL"){

            }
            else if (var.info["SVTYPE"] == "INV"){
            }
            else if (var.info["SVTYPE"] == "INS"){
            }
            else if (var.info["SVTYPE"] == "DUP"){
            }
            else if (var.info["SVTYPE"] == "BND"){
            }
            else{
            }
        }


    }

    for (auto x : name_to_seq){
        if (contig_to_variants.find(x.first) == contig_to_variants.end()){
            sequence_elem c_seq;
            c_seq.name = std::to_string(++global_max_id);
            c_seq.sequence.assign(x.second);
            gg.add_sequence(c_seq);

            // label path
            path_elem pp;
            pp.name = x.first;
            pp.rank = 1;
            pp.source_name = c_seq.name;
            pp.is_reverse = false;
            stringstream p_cig;
            p_cig << c_seq.sequence.length() << "M";
            pp.cigar = p_cig.str();
            gg.add_path(pp.source_name, pp);
        }
        
    }

    cout << gg;

    for (auto x : name_to_seq){
        delete [] x.second;
    }


    return 0;
}
