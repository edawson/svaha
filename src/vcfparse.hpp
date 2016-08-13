#ifndef VCFPARSE_HPP
#define VCFPARSE_HPP

#include <string>
#include <vector>
#include <map>
#include <stdlib.h>

//CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  REBC_UA0193_A90G-06A-11D_metastasis
//1   726500  1   N   <DEL>   .   .   SVTYPE=DEL;STRANDS=+-:4;SVLEN=-46;END=726546;CIPOS=-9,8;CIEND=-9,8;CIPOS95=0,0;CIEND95=0,0;SU=4;PE=0;SR=4   GT:SU:PE:SR ./.:4:0:4
//
//std::vector<std::string> split(string s, char delim);
//            std::string join(std::vector<std::string> splits, std::string glue);
//

using namespace std;
namespace vcfparse{

    // Borrow from
    //     //http://stackoverflow.com/questions/236129/split-a-string-in-c
    //         // Thanks StackOverflow!
    inline vector<string> split(string s, char delim){
        vector<string> ret;
        stringstream sstream(s);
        string temp;
        while(getline(sstream, temp, delim)){
            ret.push_back(temp);
        }
        return ret;
    };

    string join(vector<string> splits, string glue){
        string ret = "";
        for (int i = 0; i < splits.size(); i++){
            if (i != 0){
                ret += glue;
            }
            ret += splits[i];
        }

        return ret;
    }

    struct Variant{
        string seq;
        int pos;
        int sv_be_alpha;
        int sv_be_beta;
        string id;
        string ref;
        vector<string> alts;
        int qual;
        vector<string> filter;
        map<string, string> info;
    };

    inline Variant parse_line(string line){
        Variant ret;
        vector<string> t_splits = split(line, '\t');
        // First six fields are fixed
        //try{
        ret.seq = t_splits[0];
        ret.pos = atoi(t_splits[1].c_str());
        ret.id = t_splits[2];
        ret.ref = t_splits[3];
        
        for (auto a : split(t_splits[4], ',')){
            ret.alts.push_back(a);
        }
        
        // Filter field
        string filter = t_splits[8];
        // INFO Field - gets crazy here. Lots of ';' delimited values.
        string info_str = t_splits[9];
        // info_splits = split(info, ';');

        // Format field.

        // Sample specific ludicrous field(s)
        // These are again ';' delimited I think.

        //}
        //catch(){
        //    cerr << "VCF malformed" << endl;
        //}
        return ret;
    };

    inline vector<Variant> parse_file(string infile){

    };

    inline string to_string(struct Variant v){
        stringstream ret;
        ret << v.seq << "\t" << v.pos << "\t" << v.id << "\t" << v.ref;

        ret << endl;

        return ret.str();

    };


}
#endif
