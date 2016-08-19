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

    map<string, map<int, vector<int> > > contig_to_node_to_edges;
    map<string, map<int, vector<pair<int, string> > > > contig_to_node_to_endpoint_svtype;
    map<string, vector<int> > contig_to_breakpoints;
    //map<string, map<int, vector<tuple<bool, int, bool> > > >
    map<string, map<int, vector<pair<int, bool> > > > contig_to_node_to_endpoint_isreverse;

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
    for (int i = 0; i < variants.size(); i++){
        vcfparse::Variant var = variants[i];
        string contig = var.seq;
        int pos = var.pos;
        //vector<Node*> con = cont_nodes[contig];
        //if (pos < con.size()){
        string svtype = var.info["SVTYPE"];
        std::string::size_type sz;
        if (svtype == "DEL"){

            // These are easy:
            // simply look to make sure they have a defined length, then 
            // create an edge between the start node and node [start + len]
            // This even works if the lengths are negative.
            //
            // NB: we must subtract one from the position since genomic coordinates are 1-based and our
            // internal representations are 0-based
            if (var.info.find("SVLEN") != var.info.end()){
                contig_to_node_to_edges[var.seq][var.pos - 1].push_back( var.pos - 1 + stoi( var.info["SVLEN"], &sz) );
                //contig_to_node_to_endpoint_svtype[ var.seq ][ var.pos - 1 ].push_back( var.pos - 1 + stoi( var.info["SVLEN"], &sz) );

            }
        }
        else if (svtype == "INV"){
                // A tiny bit more complicated than deletions 
                // check if we have a defined length
                // then, we must create TWO edges
                // 1 from [ start - 1 ] to [ length ]
                // and one from [ length + 1 ] to [ start ]
                //
                // NB: as for deletions, we must subtract 1 from the posiiton due to 1/0-based difference
                // TODO weird stuff might happen if we are at the tip of a graph
            if (var.info.find("SVLEN") != var.info.end()){
                // Yes, these can be simplified.
                // But I've left the -1 there as a reminder that it's an offset
                contig_to_node_to_edges[var.seq][var.pos - 1 - 1 ].push_back( var.pos - 1 + stoi( var.info["SVLEN"], &sz) );
                //contig_to_node_to_edges[var.seq][var.pos - 1 ].push_back( var.pos - 1 + stoi( var.info["SVLEN"], &sz) );
                contig_to_node_to_edges[var.seq][var.pos - 1].push_back( var.pos - 1 + 1 + stoi( var.info["SVLEN"], &sz) );
                //contig_to_node_to_edges[var.seq][var.pos].push_back( var.pos + stoi( var.info["SVLEN"], &sz) );
            }

        }
        else if (svtype == "DUP"){
            // Two choices here
            // 1. Cyclic graph: simply create an edge from the end to the start of the repeated sequence
            // 2. Acyclic graph: unroll the sequence the number of times it is repeated,
            //      and add nodes for EACH number of possible times it is repeated. Then make edges
            //      from the node before the repeat to the node after the repeat passing through each of these.
        }
        else if (svtype == "INS"){
            // These aren't too bad, but we'll need to be able to grab the variant back from the functions below.
            // so that we can insert the sequence into the graph at the right spot.

        }
        else if (svtype == "BND"){
            // If these are paired, simply create edges
            // If they are not there isn't much we can do.
            // We should use the remapping abilities of VG to verify which breakends are valid.

        }
        else if (svtype == "TRANS"){
            // For these, we create FOUR edges (if they are reciprocal)
            // and TWO edges if not.
        }
        /* TODO: we should come up with a way to 
         * handle SNPs. They are relatively easy, but our VCF parsing and handling may have to change.
         *
         */
    }

    map<string, vector<Node*> > cont_to_nodes;

    map<string, vector<int> > cont_to_breakpoints;

    map<string, map<int, vector<int> > >::iterator it;
    for (it = contig_to_node_to_edges.begin(); it != contig_to_node_to_edges.end(); it++){
        for (auto node_brpts : it->second){
            cont_to_breakpoints[it->first].push_back(node_brpts.first);
            for (int i = 0; i < node_brpts.second.size(); i++){
                cont_to_breakpoints[it->first].push_back(node_brpts.second[i]);
            }
        }
        cont_to_breakpoints[it->first].push_back(name_to_length[it->first]);

        std::sort(cont_to_breakpoints[it->first].begin(), cont_to_breakpoints[it->first].end());
    }



    map<int, Node*> id_to_node;


    //map<string, vector<Node*> > cont_nodes;
    for (it = contig_to_node_to_edges.begin(); it != contig_to_node_to_edges.end(); it++){
        char* seq = name_to_seq[it->first];
        vector<int> breakpoints = cont_to_breakpoints[it->first];

        vector<Node*> con; 
        int n_start = 0;
        for (int i = 0; i < breakpoints.size(); i++){
           Node * nn = new Node();
           nn->sequence.assign(seq + n_start, breakpoints[i] - n_start);
           nn->id = n_start;
           nn->path = it->first;

           id_to_node[nn->id] = nn;

           if (con.size() > 0){
                con[con.size() - 1]->next.push_back(nn);
                nn->prev.push_back(con[con.size() - 1]);
           }
           con.push_back(nn);
           n_start = breakpoints[i];
        }
        cont_to_nodes[it->first] = con;
    }

    map<string, vector<Node*> >::iterator jt;
    for (jt = cont_to_nodes.begin(); jt != cont_to_nodes.end(); jt++){
        for (int i = 0; i < jt->second.size(); i++){
            if ( (contig_to_node_to_edges[jt->first][ (jt->second[i])->id ]).size() > 0){
                for (auto bp : (contig_to_node_to_edges[jt->first][ (jt->second[i])->id ])){
                    cerr << "Edge from " << jt->second[i]->id << " to " << bp << endl;
                    if (jt->second[i]->prev.size() == 1){
                        ((jt->second[i])->prev[0])->next.push_back( id_to_node[bp] );
                    }

                }
            }
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

    std::function<vector<Node*>(vector<Node*>)> compact = [](vector<Node*> contig){
        // Rules:
        // A node M may be merged with its next node N IFF:
        //      1. The outdegree of M is == 1
        //      2. The indegree of N is == 1
        //      3. id[M] < id[N] and the two IDs are consecutive

        // Algorithm: 
        // vector<Node *> ret;
        // Node* merged_node = contig[0];
        // for node in contig(1, contig.size]:
        //   if (outdegree(node - 1) == 1 and indegree(node) == 1 and
        //       id(node) - 1 == id(node)):
        //           merged_node->sequence += node.sequence
        //           merged_node->next = node.next
        //           delete (node)
        //           ret.push_back(merged_node)
        //   else:
        //       merged_node = node
        //       ret.push_back(merged_node);
        vector<Node*> ret;
        //Node* merged = memcpy(contig[0];
        Node* merged = new Node();
        merged->sequence = contig[0]->sequence;
        merged->path = contig[0]->path;
        merged->id = contig[0]->id;
        merged->next = contig[0]->next;
        ret.push_back(merged);
        //cerr << contig[0]->next.size();
        for (int i = 1; i < contig.size(); i++){
            Node* curr = contig[ i ];
            merged = ret[ ret.size () - 1 ];

            if ( merged->next.size() == 1 && curr->prev.size() == 1 &&
                    merged->path == curr->path){
                merged->sequence += curr->sequence;
                merged->next = curr->next;
                //merged->id = curr->id;
                delete curr;
            }
            else{
                //merged = curr;
                //ret.push_back(merged);
                ret.push_back(curr);
            }
        }
        return ret;
    };

    /*
       for (auto conti : cont_nodes){
//vector<Node*> orig = conti.second;
cont_nodes [ conti.first ] = compact(conti.second);

}
*/
GFAKluge og;
og.set_version();
for (auto conti : cont_to_nodes){
    string c_name = conti.first;
    vector<Node*> c_nodes = conti.second;
    for (int i = 0; i < c_nodes.size(); i++){
        Node * n = c_nodes[i];
        sequence_elem seq_el;
        seq_el.sequence = n->sequence;
        seq_el.name = std::to_string(n->id + 1);
        og.add_sequence(seq_el);

        //cerr << n->next.size() << endl;
        for (int next_ind = 0; next_ind < n->next.size(); next_ind++){
            link_elem link_el;
            link_el.source_name = std::to_string(n->id + 1);
            link_el.sink_name = std::to_string((n->next[next_ind])->id + 1);
            link_el.source_orientation_forward = true;
            link_el.sink_orientation_forward = true;
            //link_el.pos = 0;
            link_el.cigar = "0M";

            og.add_link(link_el.source_name, link_el);
        }

        if (!n->path.empty()){
            path_elem p_elem;
            p_elem.name = n->path;
            p_elem.source_name = std::to_string(n->id + 1);
            p_elem.rank = i + 1;
            p_elem.is_reverse = false;
            stringstream p_cig;
            p_cig << n->sequence.length() << "M";
            p_elem.cigar = p_cig.str();

            og.add_path(p_elem.source_name, p_elem);
        }
    }


}

cout << og;




for (auto x : cont_to_nodes){
    for (int i = 0; i < x.second.size(); i++){
        delete x.second[i];
    }
}

for (auto x : name_to_seq){
    delete [] x.second;
}



return 0;
}
