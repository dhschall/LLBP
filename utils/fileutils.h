
#ifndef __TraceReader__
#define __TraceReader__
/* code */



#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <zlib.h>

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>


#include "common.h"



struct branchTrace
{
    // "tick","pc","type","target","taken","misp","BTBmiss"
    /* data */
    uint64_t tick;
    uint64_t pc;
    BrType type;
    uint64_t target;
    bool taken;
    bool misp;
    bool btbMiss;
    uint32_t skipped;
};

bool parseLine(std::string line, branchTrace& bt);

class TraceReader
{
    private:
        std::string fileName;
        // std::string delimeter;
        bool binary;
        std::ifstream file;

    public:
        TraceReader(std::string filename, bool binary = false);
        TraceReader() {};

        bool open(std::string filename, bool binary = false);

        ~TraceReader();

        bool getNextBranch(branchTrace& bt);

};



class TraceReaderGz
{
    private:
        std::string fileName;
        // std::string delimeter;
        bool binary;
        gzFile file;

    public:
        TraceReaderGz(std::string filename, bool binary = false);
        TraceReaderGz() {};

        bool open(std::string filename, bool binary = false);

        ~TraceReaderGz();

        bool getNextBranch(branchTrace& bt);

};



class TraceReaderCSVGz
{
    protected:
        std::string fileName;
        // std::string delimeter;
        bool binary = false;
        gzFile file = nullptr;
        bool restart = false;

        static const int bufferSize = 1024;
        char buffer[bufferSize];

        bool refillBuffer();

        std::string reminingData;

    public:
        TraceReaderCSVGz(std::string filename);
        TraceReaderCSVGz() {};

        bool open(std::string filename);

        virtual ~TraceReaderCSVGz();

        bool getNextLine(std::string& line);

};

class GoogleTrace : public TraceReaderCSVGz
{
    bool parseGoogleLine(std::string line, branchTrace& bt);

  public:
    // GoogleTrace() : TraceReaderCSVGz(filename) {};
    bool getNextBranch(branchTrace& bt);
};



class ChampSimTrace : public TraceReaderCSVGz
{

    // special registers that help us identify branches
    const static char REG_STACK_POINTER = 6;
    const static char REG_FLAGS = 25;
    const static char REG_INSTRUCTION_POINTER = 26;


    // instruction format
    const static std::size_t NUM_INSTR_DESTINATIONS_SPARC = 4;
    const static std::size_t NUM_INSTR_DESTINATIONS = 2;
    const static std::size_t NUM_INSTR_SOURCES = 4;


    struct input_instr {
        // instruction pointer or PC (Program Counter)
        unsigned long long ip;

        // branch info
        unsigned char is_branch;
        unsigned char branch_taken;

        unsigned char destination_registers[NUM_INSTR_DESTINATIONS]; // output registers
        unsigned char source_registers[NUM_INSTR_SOURCES];           // input registers

        unsigned long long destination_memory[NUM_INSTR_DESTINATIONS]; // output memory
        unsigned long long source_memory[NUM_INSTR_SOURCES];           // input memory
    };


    bool covertBranch(branchTrace& bt);


    FILE *trace_file;
    std::string gunzip_command;
    input_instr cur_instr, next_instr;

    const bool restart = true;
    bool reopen(std::string gzfilename);


  public:
    // GoogleTrace() : TraceReaderCSVGz(filename) {};
    bool getNextBranch(branchTrace& bt);
    ~ChampSimTrace() override;

    bool open(const std::string filename);
};






#endif //__TraceReader__
