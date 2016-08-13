#include "wrangler.hpp"

using namespace std;
using namespace rodeo;
using namespace gfak;
namespace rodeo{
    struct Genome* Wrangler::from_fasta_and_vcf(string fasta_file, string vcf_file){

    }

    /**
     * struct Genome* Wrangler::from_parsed(map<string, char*> seqs, map<string, int> lengths,
     *                                      vector<vcflib::Variant> variants
     *
     */

    Node* Wrangler::pos_to_node(int pos){
        return pos_to_node_map[pos];
    }

    Node* Wrangler::compact(Node* n, int num_nodes){

    }

    vector<Node*> Wrangler::compact(vector<Node*> nodes){

    }

    // These are complex, comprosing both links and sequences AND PATHS
    string Wrangler::to_gfa(struct Contig* c){

    }
    
    // These are sequence (S) lines
    string Wrangler::to_gfa(struct Node* n){

    }

    // These are link (L) lines
    string Wrangler::to_gfa(struct Edge* e){

    }

    // These have contigs (L/S/P) and header lines (H)
    string Wrangler::to_gfa(struct Genome* g){

    }
}
