#pragma once

#include <inttypes.h>


// up-down saturating counter
inline void ctrupdate(int8_t& ctr, bool taken, int nbits) {
    if (taken) {
        if (ctr < ((1 << (nbits - 1)) - 1)) ctr++;
    } else {
        if (ctr > -(1 << (nbits - 1))) ctr--;
    }
}

// up-down saturating counter
template<typename T>
inline void ctrupdate(T& ctr, bool up, int nbits) {
    if (up) {
        if (ctr < ((1 << (T)nbits) - 1)) ctr++;
    } else {
        if (ctr > 0) ctr--;
    }
}

enum {LowConf = 0, MedConf = 1, HighConf = 2};
inline unsigned compConf(int8_t ctr, const int cwidth) {
    if (cwidth < 2)
        return HighConf;
    // Two bit counters saturate at +1 and -2
    if (cwidth < 3)
        return ((ctr == -2) || (ctr == 1)) ? HighConf : LowConf;

    if (abs (2 * ctr + 1) >= (1 << cwidth) - 1)
        return HighConf;
    if (abs (2 * ctr + 1) >= (1 << (cwidth - 1)) - 1)
        return MedConf;
    return LowConf;
}

inline int8_t saturate(int8_t ctr, int nbits) {
	if (ctr > ((1 << (nbits - 1)) - 1)) return ((1 << (nbits - 1)) - 1);
	if (ctr < -(1 << (nbits - 1))) return -(1 << (nbits - 1));
	return ctr;
}

inline int center(int8_t ctr) {
    return 2 * ctr + 1;
}

#define CUMAX(x) ((1 << (x)) - 1)

//counter
//MAX : max value
//MIN : min value
template<typename T, int MAX, int MIN>
class COUNTER {
private:
	T ctr;
public:
	T read() {
		return ctr;
	}

	bool pred() {
		return ctr >= 0;
	}
	bool satmax(){
		return ctr == MAX;
	}
	bool satmin(){
		return ctr == MIN;
	}
	void write(T v) {
		assert(v <= MAX);
		assert(v >= MIN);
		ctr = v;
	}
	void add(T d) {
		ctr = ctr + d;
		if (ctr > MAX){
			ctr = MAX;
		}else if (ctr < MIN){
			ctr = MIN;
		}
	}
	void update(bool incr) {
		if (incr) {
			if (ctr < MAX)
				ctr = ctr + 1;
		} else {
			if (ctr > MIN)
				ctr = ctr - 1;
		}
	}
};
//signed integer counter
template<int WIDTH>
class SCOUNTER : public COUNTER<int32_t,((1<<(WIDTH-1))-1),(-(1<<(WIDTH-1)))>{
};
//unsigned integer counter
template<int WIDTH>
class UCOUNTER : public COUNTER<int32_t,((1<<(WIDTH))-1),0>{
};

