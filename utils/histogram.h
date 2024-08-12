

// S
#pragma once

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <string>
#include <boost/format.hpp>

#include "error.h"


template <typename T>
class Histogram {

    // Define upper and lower bounds
    const T lower;
    const T upper;
    const int bins;



    // Bucket struct
    struct Bucket {
        std::string mark;
        T count;
        double frequency;
    };

    // Containers to hold the buckets and counts
    // The number of buckets is bins + 2 because there
    // is one extra bucket on each side to hold values
    // that are outside the range
    std::vector<Bucket> buckets;

    // Bucket size
    const T bs;

    // total number of elements inserted
    int samples;
    T max;
    T min;
    T sum;


  public:

    Histogram(T lower, T upper, int bins)
        : lower(lower), upper(upper), bins(bins),
          bs((upper - lower) / static_cast<T>(bins)),
          samples(0), max(0), min(0), sum(0)
    {

        panic_if(bs <= 0, "Bucket size must be positive");

        buckets.resize(bins + 2);
        {
            std::ostringstream oss;
            oss << "< " << std::fixed << std::setprecision(1) << lower;
            buckets[0] = Bucket{oss.str(), 0, 0};
        }
        {
            std::ostringstream oss;
            oss << "> " << std::fixed << std::setprecision(1) << upper;
            buckets[bins + 1] = Bucket{oss.str(), 0, 0};
        }

        for (int i = 1; i <= bins; ++i) {
            std::ostringstream oss;
            if (bs == 1) {
                oss << (lower + bs * static_cast<T>(i-1));
                buckets[i] = Bucket{oss.str(), 0, 0};
            } else {
                oss << std::fixed << std::setprecision(1)
                    << (lower + bs * static_cast<T>(i-1)) << "-"
                    << (lower + bs * static_cast<T>(i));
                buckets[i] = Bucket{oss.str(), 0, 0};
            }
        }
    }

    void insert(T value, int count=1) {
        samples += count;
        sum += value * count;
        if (value > max) {
            max = value;
        }
        if (value < min) {
            min = value;
        }

        if (value < lower) {
            buckets[0].count += count;
            return;
        } else if (value > upper) {
            buckets[bins + 1].count += count;
            return;
        }

        int idx = static_cast<int>((value - lower) / bs) + 1;
        buckets[idx].count += count;
    }

    void insert(const std::vector<T>& values) {
        for (auto v : values) {
            insert(v);
        }
    }

    void reset() {
        samples = 0;
        max = 0;
        min = 0;
        sum = 0;
        for (auto& b : buckets) {
            b.count = 0;
        }
    }


    std::string print(bool perc=false, bool cdf=false, int width = 40) {
        T _max = 0;
        for (const auto& b : buckets) {
            if (b.count > _max) {
                _max = b.count;
            }
        }

        std::ostringstream res;

        res << "N:" << samples
            << " Min:" << min
            << " Max:" << max
            << " Sum:" << sum
            << " Avg:" << std::fixed << std::setprecision(2) << getAvg()
            << "\n";


        int cum = 0;
        for (const auto& b : buckets) {
            int barLen = _max > 0 ? static_cast<int>((b.count * width + _max / 2) / _max) : 0;

            res << std::left << std::setw(10) << b.mark
                << " [" << std::setw(4) << b.count << "]\t";

            if (perc) {
                res << std::fixed << std::setprecision(1) << (100.0 * b.count / samples) << "%\t";
            }

            if (cdf) {
                cum += b.count;
                res << std::fixed << std::setprecision(1) << (100.0 * cum / samples) << "%\t";
            }

            res << "|" << std::string(barLen, '*') << "\n";

        }
        return res.str();
    }

    std::string printCDF() {
        return print(false, true);
    }

    double getAvg() {
        return static_cast<double>(sum) / static_cast<double>(samples);
    }

    T getMax() {
        return max;
    }

    T getMin() {
        return min;
    }

};

