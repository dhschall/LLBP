

#include "fileutils.h"

#include <zlib.h>

#include <boost/format.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// ------------------------ TraceReader ------------------------

TraceReader::TraceReader(std::string filename, bool binary)
{
    open(filename, binary);
}

bool
TraceReader::open(std::string _filename, bool _binary)
{
    fileName = _filename;
    binary = _binary;

    file.open(_filename, std::ios::in);
    if (!file.is_open()) {
        std::cerr << "Error opening the file for reading." << std::endl;
        return false;
    }
    return true;
}

TraceReader::~TraceReader() { file.close(); }

bool TraceReader::getNextBranch(branchTrace& bt) {
    return (file.read((char*)&bt, sizeof(branchTrace))) ? true : false;
}



// ------------------------ TraceReaderGz ------------------------

TraceReaderGz::TraceReaderGz(std::string _filename, bool _binary)
{
    open(_filename, _binary);
}

bool
TraceReaderGz::open(std::string _filename, bool _binary)
{
    fileName = _filename;
    binary = _binary;

    file = gzopen(_filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Error opening the file for reading." << std::endl;
        return false;
    }
    return true;
}

TraceReaderGz::~TraceReaderGz() { gzclose(file); }


bool TraceReaderGz::getNextBranch(branchTrace& bt) {
    // return (file.read((char*)&bt, sizeof(branchTrace)))
    //             ? true : false;
    return (gzread(file, reinterpret_cast<void*>(&bt), sizeof(branchTrace)) > 0)
               ? true
               : false;
}





// ------------------------ TraceReaderCSVGz ------------------------

TraceReaderCSVGz::TraceReaderCSVGz(std::string _filename)
{
    open(_filename);
}

bool
TraceReaderCSVGz::open(std::string _filename)
{
    fileName = _filename;

    file = gzopen(_filename.c_str(), "rb");
    if (!file) {
        std::cerr << "Error opening the file for reading." << std::endl;
        return false;
    }
    return true;
}

TraceReaderCSVGz::~TraceReaderCSVGz()
{
    if (file) gzclose(file);
}

bool TraceReaderCSVGz::refillBuffer()
{
    int bytesRead = gzread(file, buffer, bufferSize - 1);
    if (bytesRead == 0) {
        return false; // End of file
    }

    if (bytesRead < 0) {
        std::cerr << "Error reading the gzipped file." << std::endl;
        return false;
    }

    buffer[bytesRead] = '\0'; // Null-terminate the buffer

    reminingData += std::string(buffer);

    return true;
}

bool TraceReaderCSVGz::getNextLine(std::string& line)
{
    size_t found=reminingData.find('\n');
    if (found == std::string::npos) {
        // No next line found
        if (!refillBuffer()) {
            return false;
        }
        found=reminingData.find('\n');
    }

    // Read the next line
    line = reminingData.substr(0, found);
    reminingData = reminingData.substr(found+1);
    return true;
}

// ------------------------ Google Tracer ------------------------

bool GoogleTrace::parseGoogleLine(std::string line, branchTrace& bt) {
    if (line.find("0x") == std::string::npos) {
        // Skip the header line
        return false;
    }

    std::istringstream iss(line);
    std::string token;

    // pc
    std::getline(iss, token, ',');
    bt.pc = std::stoull(token, nullptr, 16);

    // type
    std::getline(iss, token, ',');
    if (token == "conditional_jump")
        bt.type = BrType::DirectCond;
    else if (token == "direct_jump")
        bt.type = BrType::DirectUncond;
    else if (token == "indirect_jump")
        bt.type = BrType::IndirectUncond;
    else if (token == "indirect_call")
        bt.type = BrType::CallIndirect;
    else if (token == "direct_call")
        bt.type = BrType::CallDirect;
    else if (token == "return")
        bt.type = BrType::Return;
    else
        bt.type = BrType::NoBranch;

    // taken
    std::getline(iss, token, ',');
    bt.taken = std::stoi(token);

    // target
    std::getline(iss, token, ',');
    // convert hex to decimal
    bt.target = std::stoull(token, nullptr, 16);

    bt.skipped = 5;

    return true;
}

bool GoogleTrace::getNextBranch(branchTrace& bt) {
    std::string line;
    if (!getNextLine(line)) {
        return false;
    }

    return parseGoogleLine(line, bt);
}




// ------------------------ ChamSim Tracer ------------------------

bool ChampSimTrace::open(std::string filename)
{
    std::string last_dot = filename.substr(filename.find_last_of("."));
    std::string decomp_program;

    std::ifstream testfile(filename);
    if (!testfile.good())
    {
        std::cerr << "TRACE FILE NOT FOUND" << std::endl;
        return false;
    }

    if (last_dot[1] == 'g') // gzip format
        decomp_program = "gzip";
    else if (last_dot[1] == 'x') // xz
        decomp_program = "xz";
    else {
        std::cout << "ChampSim does not support traces other than gz or xz compression!" << std::endl;
        return false;
    }

    gunzip_command = decomp_program + " -dc " + filename;

    trace_file = popen(gunzip_command.c_str(), "r");
    if (trace_file == NULL) {
        std::cout << "\n*** Trace file not found: " << filename << " ***\n\n";
        return false;
    }
    if (!fread(&next_instr, sizeof(input_instr), 1, trace_file)) {
        std::cout << "\n*** Cannot read from file: " << filename << " ***\n\n";
        return false;
    }
    return true;
}

bool ChampSimTrace::reopen(std::string gzfilename)
{
    // reached end of file for this trace
    std::cout << "*** Reached end of trace. Repeating trace: " << gzfilename << std::endl;

    // close the trace file and re-open it
    pclose(trace_file);

    trace_file = popen(gunzip_command.c_str(), "r");
    if (trace_file == NULL) {
        std::cout << "\n*** CANNOT REOPEN TRACE FILE ***\n\n";
        return false;
    }
    if (!fread(&next_instr, sizeof(input_instr), 1, trace_file)) {
        std::cout << "\n*** Cannot read after reopen trace ***\n\n";
        return false;
    }
    return true;
}


bool ChampSimTrace::getNextBranch(branchTrace& bt) {

    int n_skip = 0;
    while (true) {
        /* code */
        cur_instr = next_instr;
        if (!fread(&next_instr, sizeof(input_instr), 1, trace_file)) {

            if (!restart) {
                return false;
            }
            if (!reopen(gunzip_command)) {
                return false;
            }
        }
        if (cur_instr.is_branch) {
            break;
        }
        n_skip++;
    }

    bt.skipped = n_skip;
    return covertBranch(bt);
}


ChampSimTrace::~ChampSimTrace() { pclose(trace_file); }


bool ChampSimTrace::covertBranch(branchTrace& bt) {


    int num_reg_ops = 0, num_mem_ops = 0;

    bt.pc = cur_instr.ip;
    bt.taken = cur_instr.branch_taken;

    bool reads_sp = false;
    bool writes_sp = false;
    bool reads_flags = false;
    bool reads_ip = false;
    bool writes_ip = false;
    bool reads_other = false;

    for (uint32_t i = 0; i < NUM_INSTR_DESTINATIONS; i++) {
        switch (cur_instr.destination_registers[i])
        {
        case 0:
            break;
        case REG_STACK_POINTER:
            writes_sp = true;
            break;
        case REG_INSTRUCTION_POINTER:
            writes_ip = true;
            break;
        default:
            break;
        }
    }

    for (int i = 0; i < NUM_INSTR_SOURCES; i++) {
        switch (cur_instr.source_registers[i])
        {
        case 0:
            break;
        case REG_STACK_POINTER:
            reads_sp = true;
            break;
        case REG_FLAGS:
            reads_flags = true;
            break;
        case REG_INSTRUCTION_POINTER:
            reads_ip = true;
            break;
        default:
            reads_other = true;
            break;
        }

    }
    // determine what kind of branch this is, if any
    if (!reads_sp && !reads_flags && writes_ip && !reads_other)
    {
        // direct jump
        bt.taken = 1;
        bt.type = BrType::DirectUncond;
    }
    else if (!reads_sp && !reads_flags && writes_ip && reads_other)
    {
        // indirect branch
        bt.taken = 1;
        bt.type = BrType::IndirectUncond;
    }
    else if (!reads_sp && reads_ip && !writes_sp && writes_ip && reads_flags && !reads_other)
    {
        // conditional branch
        bt.taken = cur_instr.branch_taken;
        bt.type = BrType::DirectCond;
    }
    else if (reads_sp && reads_ip && writes_sp && writes_ip && !reads_flags && !reads_other)
    {
        // direct call
        bt.taken = 1;
        bt.type = BrType::CallDirect;
    }
    else if (reads_sp && reads_ip && writes_sp && writes_ip && !reads_flags && reads_other)
    {
        // indirect call
        bt.taken = 1;
        bt.type = BrType::CallIndirect;
    }
    else if (reads_sp && !reads_ip && writes_sp && writes_ip)
    {
        // return
        bt.taken = 1;
        bt.type = BrType::Return;
    }
    else if (writes_ip)
    {
        bt.taken = cur_instr.branch_taken;
        bt.type = BrType::NoBranch;
    }

    if (bt.taken) {
        bt.target = next_instr.ip;
    } else {
        bt.target = 0;
    }
    return true;
}


