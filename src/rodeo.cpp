#include "rodeo.hpp"

using namespace std;
namespace rodeo{
    struct Node{
        string seq;
        string variant_id;
        vector<Node*> next;
        vector<Node*> prev;
        vector<Edge*> edges_in;
        vector<Edge*> edges_out;
        string path;
    };   

    struct Edge{
        Edge* partner;
        Node* next;
        Node* prev;
        string variant_id;
    };

    struct Contig{
        string name;
        vector<Node*> nodes;
        vector<Edge*> edges;
        int length_helper = -1;
        inline int length(void){
            if (length_helper == -1){
                for (auto n : nodes){
                    length_helper += n->seq.length();
                }
            }
            return length_helper;
        };
        int begin;
        vector<string> paths;
    };

    struct Genome{
        string name;
        int num_contigs(void);
        vector<Contig> contigs;
    };
}
