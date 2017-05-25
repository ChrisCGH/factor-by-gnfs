#include "SparseMatrix.h"
#include "MemoryMappedFile.h"
#include "convert.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace
{
std::string DependenciesFile;
std::string RelationSetsFile;
std::string RelationsFile;
long int FirstDependency = -1;
long int LastDependency = -1;

void usage()
{
    std::cerr << "Usage: processDependencies -r relations -d dependencies -rsi relationsets [ -f firstdependency ] [ -l lastdependency ] -help" << std::endl;
    std::exit(-1);
}

void init(int argc, char* argv[])
{
    try
    {
        int arg = 1;
        while (arg < argc)
        {
            if (strcmp(argv[arg], "-dependencies") == 0 ||
                    strcmp(argv[arg], "-d") == 0)
            {
                ++arg;
                DependenciesFile = argv[arg];
            }
            else if (strcmp(argv[arg], "-f") == 0)
            {
                ++arg;
                FirstDependency = std::atoi(argv[arg]);
            }
            else if (strcmp(argv[arg], "-l") == 0)
            {
                ++arg;
                LastDependency = std::atoi(argv[arg]);
            }
            else if (strcmp(argv[arg], "-r") == 0)
            {
                ++arg;
                RelationsFile = argv[arg];
            }
            else if (strcmp(argv[arg], "-rsi") == 0 ||
                     strcmp(argv[arg], "-relation_sets_input_file") == 0)
            {
                ++arg;
                RelationSetsFile = argv[arg];
            }
            else if (strcmp(argv[arg], "-help") == 0 ||
                     strcmp(argv[arg], "-h") == 0)
            {
                usage();
            }
            ++arg;
        }
        if (RelationSetsFile.empty() || DependenciesFile.empty())
        {
            usage();
        }
        if (FirstDependency < -1 || LastDependency < -1)
        {
            usage();
        }
        if (LastDependency == -1 || LastDependency < FirstDependency)
        {
            LastDependency = FirstDependency;
        }
    }
    catch (std::exception& e)
    {
        std::cerr << "init() failed : " << e.what() << std::endl;
        throw;
    }
    catch (...)
    {
        std::cerr << "init() failed : unknown exception" << std::endl;
        throw;
    }
}

std::vector<std::fstream*> rat_dependency_files;
std::vector<std::fstream*> alg_dependency_files;

void initialise_dependency_files(size_t dependency_count)
{
    if (FirstDependency < 0)
    {
        rat_dependency_files.resize(dependency_count);
        alg_dependency_files.resize(dependency_count);
        for (size_t i = 0; i < dependency_count; ++i)
        {
            std::ostringstream oss;
            oss << RelationsFile << ".rel" << i;
            alg_dependency_files[i] = new std::fstream(oss.str().c_str(), std::ios::out);
            oss << ".rat";
            rat_dependency_files[i] = new std::fstream(oss.str().c_str(), std::ios::out);
        }
    }
    else
    {
        size_t dependency_count = LastDependency - FirstDependency + 1;
        rat_dependency_files.resize(dependency_count);
        alg_dependency_files.resize(dependency_count);
        for (size_t dep = FirstDependency; (long int)dep <= LastDependency; ++dep)
        {
            std::ostringstream oss;
            oss << RelationsFile << ".rel" << dep;
            alg_dependency_files[dep - FirstDependency] = new std::fstream(oss.str().c_str(), std::ios::out);
            oss << ".rat";
            rat_dependency_files[dep - FirstDependency] = new std::fstream(oss.str().c_str(), std::ios::out);
        }
    }
}

void close_dependency_files()
{
    for (size_t i = 0; i < rat_dependency_files.size(); ++i)
    {
        *alg_dependency_files[i] << "!" << std::endl;
        *rat_dependency_files[i] << "!" << std::endl;
        delete alg_dependency_files[i];
        delete rat_dependency_files[i];
    }
}

#if 0
void parse(const std::string& str, std::string& ab, std::string& alg_primes, std::string& rat_primes)
{
    // str is of form:
    // a b : algebraic primes : rational primes :
    // 6924553 49395 : 7/2 43/39 3209/2549 6079/1642 9173/7520 48479/4437 101839/70977 201809/77963 207433/79821 534913/467630 1041253/275531 6000047/4249308 : 2 21701 100673 826799 16975523 30659723 :
    // split out by : and remove /n from the algebraic primes
    std::string::size_type colon_pos = str.find(':', 0);
    ab = str.substr(0, colon_pos - 1);

    std::string::size_type ap_pos = colon_pos + 2;
    colon_pos = str.find(':', ap_pos);
    alg_primes = str.substr(ap_pos, colon_pos - ap_pos - 1);

    std::string::size_type slash_pos = alg_primes.find('/', 0);
    std::string::size_type space_pos = alg_primes.find(' ', slash_pos);
    while (slash_pos != std::string::npos)
    {
        if (space_pos != std::string::npos)
        {
            alg_primes.erase(slash_pos, space_pos - slash_pos);
        }
        else
        {
            alg_primes.erase(slash_pos, space_pos);
        }
        slash_pos = alg_primes.find('/', 0);
        space_pos = alg_primes.find(' ', slash_pos);
    }

    std::string::size_type rp_pos = colon_pos + 2;
    colon_pos = str.find(':', rp_pos);
    rat_primes = str.substr(rp_pos, colon_pos - rp_pos - 1);
}
#endif

void add_relation_to_dependency_files(size_t dep, const std::string& str)
{
    //std::string ab;
    //std::string alg_primes;
    //std::string rat_primes;
    long long int a;
    long long int b;
    char* as;
    char* rs;
    //parse(str, ab, alg_primes, rat_primes);
    //*(alg_dependency_files[dep]) << ab << " " << alg_primes << std::endl;
    //*(rat_dependency_files[dep]) << ab << " " << rat_primes << std::endl;
    Convert::parse_FBGNFS(str, a, b, as, rs);
    *(alg_dependency_files[dep]) << a << " " << b << " " << as << std::endl;
    *(rat_dependency_files[dep]) << a << " " << b << " " << rs << std::endl;
}
};

int main(int argc, char* argv[])
{
    try
    {
        init(argc, argv);

        // relation_sets maps relation sets to relations
        SparseMatrix2 relation_sets(RelationSetsFile);
        std::cerr << "relation_sets is (" << relation_sets.rows() << "x" << relation_sets.cols() << ")" << std::endl;

        // dependencies maps relations sets to dependencies
        BitMatrix dependencies;
        //std::fstream dep_file(DependenciesFile.c_str(), std::ios::in);
        MemoryMappedFile dep_file(DependenciesFile.c_str());
        dependencies.readTranspose(dep_file);
        std::cerr << "dependencies is (" << dependencies.rows() << "x" << dependencies.cols() << ")" << std::endl;

        // relations_dependencies maps relations to dependencies
        BitMatrix relations_dependencies;
        multiplyt(relation_sets, dependencies, relations_dependencies);
        std::cerr << "relations_dependencies is (" << relations_dependencies.rows() << "x" << relations_dependencies.cols() << ")" << std::endl;

        initialise_dependency_files(relations_dependencies.cols());
        MemoryMappedFile relations(RelationsFile.c_str());

        std::string str;
        size_t relation = 0L;
        while (getline(relations, str))
        {
            unsigned long int drow = relations_dependencies.row_[relation];
            if (FirstDependency <0)
            {
                for (size_t dep = 0; dep < relations_dependencies.cols(); ++dep)
                {
                    if (BitOperations::bitSet(dep, drow))
                    {
                        add_relation_to_dependency_files(dep, str);
                    }
                }
            }
            else
            {
                for (size_t dep = FirstDependency; (long int)dep <= LastDependency && dep < relations_dependencies.cols(); ++dep)
                {
                    if (BitOperations::bitSet(dep, drow))
                    {
                        add_relation_to_dependency_files(dep - FirstDependency, str);
                    }
                }
            }
            ++relation;
        }

        close_dependency_files();

    }
    catch (...)
    {
        std::cerr << "Exception caught" << std::endl;
    }
}
