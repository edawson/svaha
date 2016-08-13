#ifndef WRANGLER_HPP
#define WRANGLER_HPP

#include <string>
#include <iostream>
#include <unordered_map>
#include "rodeo.hpp"
#include "gfakluge.hpp"
#include "Variant.h"

using namespace std;
using namespace gfak;
namespace rodeo{
    class Wrangler{


        public:
            struct Genome* from_fasta_and_vcf(string fasta_file, string vcf_file);

            Node* pos_to_node(int pos);

            Node* compact(Node*, int num_nodes);

            vector<Node*> compact(vector<Node*> nodes);

            string to_gfa(struct Contig* c);

            string to_gfa(struct Node* n);

            string to_gfa(struct Edge* e);

            string to_gfa(struct Genome* g);

        private:
            unordered_map<int, struct Node*> pos_to_node_map;


    };
}
#endif
