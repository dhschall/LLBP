#pragma once

#include <inttypes.h>
#include <stdlib.h>
#include <vector>
#include <string>



//==========================================================
// History Code:
// The code is based on the "Dynamically Sizing the TAGE Branch Predictor"
// paper by Stephen Pruett submitted to the CBP-5 workshop.



// [PPM, page 4] Discusses how to build a low latency FoldedHistory register
struct HistoryRegister {
   public:
    uint32_t size;
    uint32_t head;
    std::vector<bool> history;
    long long history_l;

    void init(uint32_t s) {
        size = s;
        history.resize(size);

        for (uint32_t i = 0; i < size; ++i) {
            history[i] = false;
        }
        history_l = 0;

        head = 0;
    }

    HistoryRegister() {}

    HistoryRegister(uint32_t s) { init(s); }

    void push(bool p) {
        head = (head + 1) % size;
        history[head] = p;

        history_l <<= 1;
        history_l += (p & 0x1);
    }

    bool operator[](const uint32_t i) {
        uint32_t index = (head + size - i) % size;
        assert(index < size);
        return history[index];
    }

    void print() {
        printf("History");
        for (uint32_t i = 0; i < size; ++i) {
            printf("%d, ", (bool)history[(head - i) % size]);
        }
        printf("\n");
    }

    long long getHistory() { return history_l; }

    uint32_t getSize() { return size; }
};

struct PCHistoryRegister {
   public:
    uint32_t size;
    uint32_t head;
    typedef std::pair<uint64_t,bool> entry_t;
    std::vector<entry_t> history;
    long long history_l;

    void init(uint32_t s) {
        size = s;
        history.resize(size);

        for (uint32_t i = 0; i < size; ++i) {
            history[i] = entry_t(0,false);
        }
        history_l = 0;

        head = 0;
    }

    PCHistoryRegister() {}

    PCHistoryRegister(uint32_t s) { init(s); }

    void push(uint64_t pc, bool t) {
        head = (head + 1) % size;
        history[head] = entry_t(pc,t);
    }

    entry_t operator[](const uint32_t i) {
        uint32_t index = (head + size - i) % size;
        assert(index < size);
        return history[index];
    }

    void print() {
        printf("History");
        for (uint32_t i = 0; i < size; ++i) {
            auto e = history[(head - i) % size];
            printf("%lu:%d ",e.first, (bool)e.second);
        }
        printf("\n");
    }

    std::string toStr() {
        std::string s = "";
        for (uint32_t i = 0; i < size; ++i) {
            auto e = history[(head - i) % size];
            // s += std::to_string(e.first) + "," + std::to_string(e.second) + ",";
            s += std::to_string(e.first) + ",";
        }
        return s;
    }


    uint32_t getSize() { return size; }
};


class FoldedHistory {
   private:
    uint32_t inputWidth;   // size of emulated history register
    uint32_t outputWidth;  // size of folded register
    uint32_t
        maxOutputWidth;  // first width register is set to. Used to calc size.
    int32_t remainder;
    int32_t value;
    HistoryRegister* ghr;

    FoldedHistory() {}

   public:
    FoldedHistory(HistoryRegister* g, uint32_t iw, uint32_t ow) {
        inputWidth = iw;
        outputWidth = ow;
        maxOutputWidth = outputWidth;
        ghr = g;

        // using a 32-bit integer as register
        //  -need an extra bit, so max is 31 bits...
        assert(outputWidth < 32);
        assert(outputWidth != 0);
        remainder = inputWidth % outputWidth;
        value = 0;
    }

    // Expectation is that FoldedHistory push is called
    //   after HistoryRegister push
    void update() {
        // input bit most recent shifted into ghr
        bool inBit = (*ghr)[0];

        // Shift in bit
        value = (value << 1) | (inBit ? 0x01 : 0x00);

        // Fold shifted-out bit in
        value = value ^ (value >> outputWidth);
        value = value & ((1 << outputWidth) - 1);

        // Get bit to shift out from Global History
        bool outputBit = (*ghr)[inputWidth];
        int32_t outputValue = (outputBit) ? (0x01 << (remainder)) : 0x0;

        // Shift out bit
        value = value ^ outputValue;
    }

    inline int32_t get() { return value; }

    void reset() { value = 0; }

    uint32_t getSize() { return maxOutputWidth; }
};




// [PPM, page 4] Discusses how to build a low latency FoldedHistory register
struct HistoryRegisterFast {
   public:
    const uint32_t size;
    uint32_t head;

    uint8_t* history;
    long long history_l;

    HistoryRegisterFast(uint32_t s)
        : size(s)
    {
        history = new uint8_t[size]();

        for (uint32_t i = 0; i < size; ++i) {
            history[i] = false;
        }
        history_l = 0;

        head = 0;
    }

    void push(bool p) {
        head--;
        history[head & (size - 1)] = p;
        history_l = (history_l << 1) + (p & 0x1);
    }

    inline bool operator[](const uint32_t i) {
        return history[(head + i) & (size - 1)];
    }

    void print() {
        printf("History");
        for (uint32_t i = 0; i < size; ++i) {
            printf("%d, ", (bool)history[(head - i) % size]);
        }
        printf("\n");
    }

    long long getHistory() { return history_l; }

    uint32_t getSize() { return size; }
};


class FoldedHistoryFast {
   private:
    const uint32_t inputWidth;   // size of emulated history register
    const uint32_t outputWidth;  // size of folded register
    const uint32_t
        maxOutputWidth;  // first width register is set to. Used to calc size.
    const int32_t remainder;
    HistoryRegisterFast& ghr;  // Reference to global history register


   public:
    int32_t value;

    FoldedHistoryFast(HistoryRegisterFast& g, uint32_t iw, uint32_t ow)
        : inputWidth(iw),
          outputWidth(ow),
          maxOutputWidth(ow),
          remainder(iw % ow),
          ghr(g)
    {
        // using a 32-bit integer as register
        //  -need an extra bit, so max is 31 bits...
        assert(outputWidth < 32);
        assert(outputWidth != 0);
        value = 0;
    }

    // Expectation is that FoldedHistory push is called
    //   after HistoryRegister push
    void update() {



        // // Shift in new bit
        // value = (value << 1) | ghr[0];

        // // Get bit to shift out from Global History
        // value = value ^ ghr[inputWidth] << outputWidth;

        // // Fold and mask
        // value = value ^ (value >> outputWidth);
        // value = value & ((1 << outputWidth) - 1);




        // input bit most recent shifted into ghr
        bool inBit = ghr[0];

        // Shift in bit
        value = (value << 1) | (inBit ? 0x01 : 0x00);







        // Fold shifted-out bit in
        value = value ^ (value >> outputWidth);
        value = value & ((1 << outputWidth) - 1);

        // Get bit to shift out from Global History
        bool outputBit = ghr[inputWidth];
        int32_t outputValue = (outputBit) ? (0x01 << (remainder)) : 0x0;

        // Shift out bit
        value = value ^ outputValue;
    }

    inline int32_t get() { return value; }

    void reset() { value = 0; }

    uint32_t getSize() { return maxOutputWidth; }
};








// utility class for index computation
// this is the cyclic shift register for folding
// a long global history into a smaller number of bits; see P. Michaud's
// PPM-like predictor at CBP-1
class folded_history {
    public:
    unsigned comp;
    int CLENGTH;
    int OLENGTH;
    int OUTPOINT;
    int histbufferlen;

    folded_history() {}

    void init(int original_length, int compressed_length, int _histbufferlen) {
        comp = 0;
        OLENGTH = original_length;
        CLENGTH = compressed_length;
        OUTPOINT = OLENGTH % CLENGTH;
        histbufferlen = _histbufferlen;
    }

    void update(uint8_t *h, int PT) {
        comp = (comp << 1) ^ h[PT & (histbufferlen - 1)];
        comp ^= h[(PT + OLENGTH) & (histbufferlen - 1)] << OUTPOINT;
        comp ^= (comp >> CLENGTH);
        comp = (comp) & ((1 << CLENGTH) - 1);
    }
};

class Bentry  // TAGE bimodal table entry
{
    public:
    int8_t ctr;
    int8_t hyst;
    int8_t pred;
    uint64_t pc;
    uint64_t id;

    Bentry() {
        ctr = -1;
        pred = 0;
        hyst = 1;
        pc = 0;
        id = 0;
    }
};

class Gentry  // TAGE global table entry
{
    public:
    int8_t ctr;
    uint tag;
    int8_t u;
    int correct;
    int incorrect;
    int useful;
    uint64_t pc;
    int hlen;
    int idx;
    uint64_t key;
    uint64_t ctx;

    Gentry() {
        ctr = 0;
        u = 0;
        tag = 0;
        correct = -1;
        incorrect = 0;
        useful = 0;
        pc = 0;
        key = 0;
        ctx = 0;
    }
};




// FTL++


//history register
class HISTORY {
private:
	int MAXHIST; // constant value

	bool *bhr; // branch history

public:
	void init(int HLENGTH) {
		MAXHIST = HLENGTH+1; // set constant value
		bhr = new bool[MAXHIST];
		for (int i = 0; i < MAXHIST; i++){
			bhr[i] = false;
		}
	}

	bool read(int n) {
		assert(n < MAXHIST);
		return bhr[n];
	}

	uint32_t read(int from, int n) {
		assert(n < MAXHIST);
		int r = 0;
		for (int i = from; i < n; ++i) {
			r ^= bhr[i] << ((i - from) % 32);
		}
		return r;
	}
	//push
	void push(bool t) {
		for(int i=MAXHIST-1; i>=0; i--) {
			bhr[i] = bhr[i-1];
		}
 		bhr[0] = t;
	}
};

//folded history
template<int WIDTH>
class FOLDHIST {
private:
	//path history length
	static const int PathHistoryLength = 32;

private:
	uint32_t start[2]; // constant value
	uint32_t end[2]; // constant value
	uint32_t pos[3]; // constant value

	uint32_t comp; // folded history

public:

	void init(int s, int e,int s2,int e2) {
		comp = 0;
		start[0] = s;
		start[1] = s2;
		end[0] = e;
		end[1] = e2;
		pos[0] = 0;
		pos[1] = (pos[0] + end[0] - start[0]) % WIDTH;
		pos[2] = (pos[1] + end[1] - start[1]) % WIDTH;
		assert(pos[0] < WIDTH);
		assert(pos[1] < WIDTH);
		assert(pos[2] < WIDTH);
	}
	void init(int s, int e) {
		comp = 0;
		start[0] = s;
		start[1] = s;
		end[0] = e;
		end[1] = e;
		pos[0] = 0;
		pos[1] = (pos[0] + end[0] - start[0]) % WIDTH;
		pos[2] = (pos[1] + end[1] - start[1]) % WIDTH;
		assert(pos[0] < WIDTH);
		assert(pos[1] < WIDTH);
		assert(pos[2] < WIDTH);
	}

	void init(int l) {
		init(0, l,0,(l>=PathHistoryLength)?PathHistoryLength:l);
	}

	uint32_t read(uint32_t pc) {
		assert(comp >= 0);
		assert(comp < (1 << WIDTH));

		pc &= (1 << WIDTH) - 1;
		return pc ^ comp;
	}
	uint32_t read(uint32_t pc, int rot) {
		assert(comp >= 0);
		assert(comp < (1 << WIDTH));
		uint32_t r = rot % WIDTH;
		uint32_t p = pc & ((1 << WIDTH) - 1);
		p = (p << r);
		p = p ^ (p >> WIDTH);
		p = (p & ((1 << WIDTH) - 1));
		return p ^ comp;
	}
	void update(class HISTORY &ghr, class HISTORY &phr) {
		comp = (comp << 1);
		comp |= (comp >> WIDTH);
		comp &= ((1 << WIDTH) - 1);
		comp ^= (ghr.read(start[0]) ? 1 : 0) << pos[0];
		comp ^= (ghr.read(end[0]) ? 1 : 0) << pos[1];
		comp ^= (phr.read(start[1]) ? 1 : 0) << pos[1];
		comp ^= (phr.read(end[1]) ? 1 : 0) << pos[2];
	}
};
