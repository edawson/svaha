#include <iostream>
#include <functional>
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
    vector<Node*> next;
    vector<Node*> prev;
};

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

//    map<string, map<int, Node*> > pos_to_node;
    map<string, vector<Node*> > cont_nodes;
    map<string, char*>::iterator it;
    for (it = name_to_seq.begin(); it != name_to_seq.end(); it++){
        string name = it->first;
        char* seq = it->second;
        int len = name_to_length[name];
        for (int i = 0; i < len; i++){
            Node* n = new Node();
            n->sequence = seq[i];
            n->id = i;
            cont_nodes[name].push_back(n);
        }
    }

    for (auto x : cont_nodes){
        for (int i = 1; i < x.second.size() - 1; i++){
            x.second[i]->prev.push_back( x.second[i-1] );
            x.second[i]->next.push_back( x.second[i+1] );
        }
    }

    cerr << cont_nodes.size() << " Contigs parsed and turned into nodes." << endl;

    vector<vcfparse::Variant> variants;
    vcfparse::Variant var;

    std::ifstream infile(var_files[0]);
    std::string line;
    while (std::getline(infile, line))
    {   
        if (line[0] != '#'){
            variants.push_back(vcfparse::parse_line(line));
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
    #pragma omp parallel for
    for (int i = 0; i < variants.size(); i++){
       vcfparse::Variant var = variants[i];
       string contig = var.seq;
       int pos = var.pos;
       if (cont_nodes.find(contig) != cont_nodes.end()){
           vector<Node*> con = cont_nodes[contig];
           if (pos < con.size()){
                string svtype = var.info["SVTYPE"];
                std::string::size_type sz;
                if (svtype == "DEL"){
                    Node * start = con[ var.pos - 1 ];
                    Node * end = con[stoi( var.info["END"], &sz) - 1 ];
                    #pragma omp critical
                    cerr << "Creating edge between " << start->id << " and " << end->id << endl;


                                        //start->next.push_back(end);
                    //end->prev.push_back(start);

                }
                else if (svtype == "INV"){

                }
                else if (svtype == "DUP"){

                }
           }
           else{
                #pragma omp critical
                cerr << "Error: position greater than contig length." << endl
                 << "Contig length: " << con.size() << "\tPosition: " << pos << endl;
                exit(1);
           }
       }
       else{
#pragma omp critical
           cerr << "Error: vcf file contains variants that are on contigs not present in the reference." << endl
               << "Are you sure the VCF comes from the right reference? Could contigs be named differently (e.g. CHR1 instread of 1)?" << endl
               << "exiting" << endl;
           exit(1);
       }
    }

    /**
     * To transform this into GFA, we make an S entry for each node
     * and an L entry for each of its NEXT nodes. Each contig is a 
     * P entry ... we'll have to make sure the label these right.
     *
     * Also, we may want to collapse nodes. One way to do this is to pick a node at random, check
     * if it is of length == 1 and check its in/out degree. Then, loop down its neighbors in both directions
     * until one of these conditions fails (either the node is already collapsed [len > 1 ] or the node
     * has in/outdegree greater than 1/1.
     *
     */


    for (auto x : cont_nodes){
        for (int i = 0; i < x.second.size(); i++){
            delete x.second[i];
        }
    }

    for (auto x : name_to_seq){
        delete [] x.second;
    }



        return 0;
}
