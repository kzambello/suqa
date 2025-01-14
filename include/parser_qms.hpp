#include <string>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;

struct arg_list{
    double beta = 0.0;
    double t = 0.0;
    double U = 0.0;
    double mu = 0.0;
    int metro_steps = 0;
    int reset_each = 0;
    int ene_qbits = 0;
    string outfile = "";
    int max_reverse_attempts = 100;
    unsigned long long int seed = 0;
    double ene_min = 0.0;
    double ene_max = 1.0;
    int pe_steps = 10;
    int thermalization = 100;
    bool record_reverse = false;
    int nbatches = 20;
    int ene_threshold = 4;

    friend ostream& operator<<(ostream& o, const arg_list& al);
};

ostream& operator<<(ostream& o, const arg_list& al){
    o<<"beta: "<<al.beta<<endl;
    o<<"t: "<<al.t<<endl;
    o<<"U: "<<al.U<<endl;
    o<<"mu: "<<al.mu<<endl;
    o<<"metro steps: "<<al.metro_steps<<endl;
    o<<"reset each: "<<al.reset_each<<endl;
    o<<"num E qbits "<<al.ene_qbits<<endl;
    o<<"max reverse attempts: "<<al.max_reverse_attempts<<endl;
    o<<"seed: "<<al.seed<<endl;
    o<<"out datafile: "<<al.outfile<<endl;
    o<<"min energy: "<<al.ene_min<<endl;
    o<<"max energy: "<<al.ene_max<<endl;
    o<<"steps of PE evolution: "<<al.pe_steps<<endl;
    o<<"thermalization: "<<al.thermalization<<endl;
    o<<"record reverse: "<<al.record_reverse<<endl;
    o<<"number of batches: "<<al.nbatches<<endl;
    o<<"ene threshold: "<<al.ene_threshold<<endl;
    return o;
}

void parse_arguments(arg_list& args, int argc, char** argv){
    int fixed_args = 8;
    map<string,int> argmap;
    map<int,string> argmap_inv;
    char *end;
    int base_strtoull = 10;

    // fixed arguments
    args.beta = stod(argv[1],NULL);
    args.t = stod(argv[2],NULL);
    args.U = stod(argv[3],NULL);
    args.mu = stod(argv[4],NULL);
    args.metro_steps = atoi(argv[5]);
    args.reset_each = atoi(argv[6]);
    args.ene_qbits = atoi(argv[7]);
    args.outfile = argv[8];

    // floating arguments
    for(int i = fixed_args+1; i < argc; ++i){
        argmap[argv[i]]=i;
        argmap_inv[i]=argv[i];
    }
    int tmp_idx;

    // flagged arguments without value

    // (bool) record_reverse
    tmp_idx = argmap["--record-reverse"];
    if(tmp_idx>=fixed_args){
        args.record_reverse = true;
    }


    // flagged arguments without value

    // (int) max_reverse_attempts
    tmp_idx = argmap["--max-reverse"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--max-reverse' flag"; 
       
       args.max_reverse_attempts = atoi(argmap_inv[tmp_idx+1].c_str()); 
    }

    // (unsigned long long) seed
    tmp_idx = argmap["--seed"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--seed' flag"; 
       
       args.seed = strtoull(argmap_inv[tmp_idx+1].c_str(), &end, base_strtoull); 
    }

    // (double) ene_min
    tmp_idx = argmap["--ene-min"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--ene-min' flag"; 
       
       args.ene_min = stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }

    // (double) ene_max
    tmp_idx = argmap["--ene-max"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--ene-max' flag"; 
       
       args.ene_max = stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }

    // (int) pe_steps
    tmp_idx = argmap["--PE-steps"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--PE-steps' flag"; 
       
       args.pe_steps = stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }


    // (int) thermalization
    tmp_idx = argmap["--thermalization"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--thermalization' flag"; 
       
       args.thermalization = stod(argmap_inv[tmp_idx+1].c_str(), NULL); 
    }

    // (int) nbatches
    tmp_idx = argmap["--nbatches"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--nbatches' flag"; 
       
       args.nbatches = atoi(argmap_inv[tmp_idx+1].c_str()); 
    }

    // (int) ene-threshold
    tmp_idx = argmap["--ene-threshold"];
    if(tmp_idx>=fixed_args){
       if(tmp_idx+1>= argc)
           throw "ERROR: set value after '--ene-threshold' flag"; 
       
       args.ene_threshold = atoi(argmap_inv[tmp_idx+1].c_str()); 
    }


    // argument checking
    if(args.beta <= 0.0){
        throw "ERROR: argument <beta> invalid";
    }

    if(args.t <= 0.0){
        throw "ERROR: argument <t> invalid";
    }

    if(args.U <= 0.0){
        throw "ERROR: argument <U> invalid";
    }

    //if(args.mu <= 0.0){
    //    throw "ERROR: argument <mu> invalid";
    //}

    if(args.metro_steps <= 0){
        throw "ERROR: argument <metro steps> invalid";
    }

    if(args.reset_each <=0){
        throw "ERROR: argument <reset each> non positive";
    }
    
    if(args.ene_qbits <=0){
        throw "ERROR: argument <num ene qbits> non positive";
    }

    if(args.outfile == ""){
        throw "ERROR: argument <output file path> empty";
    }

    if(args.max_reverse_attempts <=0){
        throw "ERROR: argument <max reverse attempts> non positive";
    }

    if(args.pe_steps <=0){
        throw "ERROR: argument <steps of PE evolution> non positive";
    }

    if(args.nbatches <=0){
        throw "ERROR: argument <num batches> non positive";
    }

    if(args.ene_threshold <=0){
        throw "ERROR: argument <ene threshold> non positive";
    }
}
