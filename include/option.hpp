
#include <boost/program_options.hpp>
#include "global.hpp"

void options(unsigned argc, char **argv) {
    /////////
    // initialize the parameters here
    /////////
    
    using namespace boost::program_options;
    options_description desc{"Options"};
    desc.add_options() 
        ("data,d", value<std::string>()->default_value(""), "input data path")
        ("size,s", value<unsigned>()->default_value(256), "subgraph size")
        ("round,r", value<unsigned>()->default_value(3), "rounds")
        ("thread,t", value<unsigned>()->default_value(20), "threads")
        ("vertex,v", value<unsigned>()->default_value(MASK::MAX_UINT), "root vertex")
        ("iter,i", value<unsigned>()->default_value(20), "iterations")
        ("dynamic,y", value<bool>()->default_value(true), "dynamic split or not (1|0)")
        ("output,o", value<bool>()->default_value(false), "output the graph as *.mix or not (1|0)")
        ("verify,f", value<bool>()->default_value(false), "verify or not (1|0)")
        ;

    
    boost::program_options::variables_map vm;
    try
    {
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, desc), vm);
        boost::program_options::notify(vm);
    } catch (error& e) {
        std::cerr << "ERROR: " << e.what() << '\n' << '\n' << desc << '\n';
        exit(1);
    }

    params::input_file = vm["data"].as<std::string>();
    if (params::input_file.empty()) {
        std::cout << desc << '\n';
        exit(1);
    } 
    params::subgraph_size = vm["size"].as<unsigned>();
    params::iters = vm["iter"].as<unsigned>();
    params::rounds = vm["round"].as<unsigned>();
    params::threads  = vm["thread"].as<unsigned>();
    params::root_vertex = vm["vertex"].as<unsigned>();
    params::dynamic = vm["dynamic"].as<bool>();
    params::output = vm["output"].as<bool>();
    params::verify = vm["verify"].as<bool>();
    std::cout << "--------experimental setting--------"<<'\n';
    std::cout << "threads: " << params::threads 
              << " rounds: " << params::rounds
              << " iterations: " << params::iters 
              << " sub_size: " << params::subgraph_size << "KB (" << params::subgraph_size*1024/sizeof(VerxId) << ")"
              << '\n'
              << ", dynamic split: " << Bool_Str[params::dynamic]
              << ", verify: " << Bool_Str[params::verify]
              <<'\n';
            
}

