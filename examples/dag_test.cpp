#include <stp/bench_reader.hpp>
#include <stp/stp_dag2nmx.hpp>
#include <stp/stp_simulator.hpp>
using namespace stp;

//exapples test：  ./examples/dag_test ../benchmarks/c17_lut.bench
int main(int argc, char **argv)
{
	if (argc < 2)
	{
		std::cout << "no input file specified" << std::endl;
		return 1;
	}
	std::ifstream ifs(argv[1]);
	if (!ifs.good())
	{
		std::cout << "can't open file" << argv[1] << std::endl;
		return 1;
	}
	stp_circuit c;
	bench_reader parser;
	if (!parser.parse(ifs, c))
	{
		std::cout << "can't parse file" << argv[1] << std::endl;
		return 1;
	}

	// c.print_circuit();
	
	std::cout << "***************************************" << std::endl;
	circuit_normalize_impl cn(c, false);
	std::string m1 = cn.run_str(false);
	std::string m2 = cn.run_str(true);
	std::cout << "old methon\n";
	std::cout << m1 << "\n";
	std::cout << "new methon\n";
	std::cout << m2 << "\n";

	stp_simulator sim(c);
	std::cout << "simulation result\n";
	std::cout << sim.run() << "\n";

	return 0;
}