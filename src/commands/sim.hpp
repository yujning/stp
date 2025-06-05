/* also: Advanced Logic Synthesis and Optimization tool
 * Copyright (C) 2019- Ningbo University, Ningbo, China */

/**
 * @file imgff.hpp
 *
 * @brief construct a fanout-free implication logic network
 *
 * @author Zhufei Chu
 * @since  0.1
 */

#ifndef SIM_HPP
#define SIM_HPP
#include <chrono>
#include <alice/alice.hpp>
#include "../include/io/bench_parser.hpp"
#include "../include/algorithms/stp_circuit.hpp"
#include "../include/algorithms/normalize.hpp"
#include "../include/algorithms/excute.hpp"

using namespace stp;



std::string remove_quotes(const std::string& str) 
{
    std::string result = str;

    if (!result.empty() && result.front() == '"' && result.back() == '"') {
        result.erase(result.begin());
        result.pop_back();
    }
    return result;
}


namespace alice
{
    //./your_tool sim --aig your_file.bench --verbose
    class sim_command : public command
    {
    public:
        explicit sim_command(const environment::ptr &env) : command(env, "sim")
        {
            add_flag( "--verbose", "verbose output" );
            add_flag( "--lut, -l",  "using lut network" );
            add_flag( "--aig, -a",  "using aig network" );
            add_flag("--cuda, -c", "using cuda");
            add_option("filename", filename ,"input file name", true);
        }
        
        protected:
        void execute()
        {
            auto start = std::chrono::high_resolution_clock::now();
            if (filename.size()==0)
            {
                std::cout << "please specify the file " << std::endl;
                return;
            }
            // get path of input file
            std::string file_path = remove_quotes(filename);

            std::ifstream ifs(file_path);

            if(!ifs.good())
            {
                std::cout << "can't open file " << file_path << std::endl;
                return;
            }

            if ( is_set("lut") || is_set("-l") )
            {
                stp_circuit c;
                bench_reader parser;

                if ( !parser.parse( ifs, c ) )
                {
                    std::cout << "can't parse file" << file_path << std::endl;
                    return;
                }
                c.update_levels();

                std::cout << "***************************************" << std::endl;
                circuit_normalize_impl cn( c, is_set( "verbose" ) );

                if (is_set("cuda") || is_set("-c"))
                {
                    #ifdef ENABLE_CUDA  
                    Get_Total_Thread_Num();
                    std::string m3 = cn.run_str(true);
                    std::cout << "cuda method\n";
                    std::cout << m3 << "\n";
                    #else
                    std::cout << "can't find cuda" << std::endl;
                    #endif
                }
                else
                {
                    std::string m3 = cn.run_str(false);
                    std::cout << "rsult: \n";
                    std::cout << m3 << "\n";
                }
            }
            else
            {
              std::cout << "please specify the network type" << std::endl;
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            std::cout << "time: " << time << " us\n"<< std::endl;

        }
        
    private:
        std::string filename{};
    };
    ALICE_ADD_COMMAND(sim, "New Command");
}

#endif
