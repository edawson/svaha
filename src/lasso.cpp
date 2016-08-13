#include <iostream>
#include <getopt.h>
#include <vector>
#include <unordered_map>
#include <zlib.h>
#include <fstream>

#include "rodeo.hpp"
#include "wrangler.hpp"
#include "kseq.hpp"
#include "gfakluge.hpp"
#include "Variant.h"
#include "vcfparse.hpp"
#include <map>


using namespace std;
using namespace gfak;
using namespace rodeo;

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
    // parse fasta reference files using kseq
    parse_fastas(ref_files, name_to_seq, name_to_length); 

    // parse vcf files using VCFlib
    //vcflib::VariantCallFile v_file;
    //v_file.open(var_files[0]);
    //vcflib::Variant var;
    //while (v_file.getNextVariant(var)){
    //    cerr << var << endl;
    //}

    map<string, vector<Node*> > cont_nodes;
    map<string, vector<Edge*> > cont_edges;
    map<string, char*>::iterator it;
    for (it = name_to_seq.begin(); it != name_to_seq.end(); it++){
        string name = it->first;
        char* seq = it->second;
        int len = name_to_length[name];
        for (int i = 0; i < len; i++){
//            Node* n = new Node();
        }
    }

    vcfparse::Variant var;

    std::ifstream infile(var_files[0]);
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line[0] != '#'){
            var = vcfparse::parse_line(line);
            //cout << vcfparse::to_string(var);
        }
    }
    Wrangler wrg;




        return 0;
}
