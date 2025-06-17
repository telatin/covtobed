#include <benchmark/benchmark.h>
#include <random>
#include <queue>

// Include only the classes we need for benchmarking
#include <queue>
#include <vector>
#include <iostream>
#include <map>
#include <cstdlib>
#include <cassert>

using namespace std;

typedef uint32_t DepthType;

// Copy the classes we need for benchmarking
struct CovEnd {
    int end;
    bool rev;
    bool operator<(const CovEnd &o) const {
        return end > o.end;
    }
};

struct Coverage {
    DepthType f = 0, r = 0;
    operator DepthType() const { return f + r; }
    void inc(bool rev=false) {
        if (rev) ++r;
        else ++f;
    }
    void dec(bool rev=false) {
        if (rev) --r;
        else --f;
    }
    bool equal(const Coverage &o, bool stranded) const {
        if (stranded)
            return f == o.f && r == o.r;
        else
            return f + r == o.f + o.r;
    }
};

struct CoordinateInterval {
    int start, end;
    bool empty() const { return end <= start; }
    size_t length() const { return empty() ? 0 : end - start; }
    operator bool() const { return !empty(); }
    bool operator<(const CoordinateInterval &i) const { 
        return start < i.start || (start == i.start && end < i.end); 
    }
    bool operator<<(const CoordinateInterval &i) const { return end <= i.start; }
    CoordinateInterval operator*(const CoordinateInterval &o) const { 
        return {max(start, o.start), min(end, o.end)}; 
    }
    pair<CoordinateInterval, CoordinateInterval> operator-(const CoordinateInterval &o) const {
        return {{start, min(end, o.start)}, {max(start, o.end), end}};
    }
};

// Benchmark Coverage class operations
static void BM_CoverageIncrement(benchmark::State& state) {
    Coverage cov;
    for (auto _ : state) {
        cov.inc(false);
        benchmark::DoNotOptimize(cov);
    }
}
BENCHMARK(BM_CoverageIncrement);

static void BM_CoverageDecrement(benchmark::State& state) {
    Coverage cov;
    // Pre-populate with some coverage
    for (int i = 0; i < 1000; ++i) {
        cov.inc(i % 2 == 0);
    }
    
    for (auto _ : state) {
        cov.dec(false);
        benchmark::DoNotOptimize(cov);
        // Re-increment to avoid going negative
        cov.inc(false);
    }
}
BENCHMARK(BM_CoverageDecrement);

static void BM_CoverageEqual(benchmark::State& state) {
    Coverage cov1, cov2;
    cov1.f = 100; cov1.r = 50;
    cov2.f = 75;  cov2.r = 75;
    
    for (auto _ : state) {
        bool result = cov1.equal(cov2, false);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(BM_CoverageEqual);

// Benchmark priority queue operations (core algorithm)
static void BM_PriorityQueueOperations(benchmark::State& state) {
    const int num_operations = state.range(0);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 10000);
    
    for (auto _ : state) {
        std::priority_queue<CovEnd> pq;
        
        // Add random end positions
        for (int i = 0; i < num_operations; ++i) {
            pq.push({dis(gen), i % 2 == 0});
        }
        
        // Process all elements
        while (!pq.empty()) {
            auto top = pq.top();
            pq.pop();
            benchmark::DoNotOptimize(top);
        }
    }
}
BENCHMARK(BM_PriorityQueueOperations)->Range(100, 10000);

// Benchmark interval operations
static void BM_IntervalIntersection(benchmark::State& state) {
    CoordinateInterval a = {100, 500};
    CoordinateInterval b = {300, 700};
    
    for (auto _ : state) {
        auto intersection = a * b;
        benchmark::DoNotOptimize(intersection);
    }
}
BENCHMARK(BM_IntervalIntersection);

static void BM_IntervalDifference(benchmark::State& state) {
    CoordinateInterval a = {100, 500};
    CoordinateInterval b = {200, 300};
    
    for (auto _ : state) {
        auto diff = a - b;
        benchmark::DoNotOptimize(diff);
    }
}
BENCHMARK(BM_IntervalDifference);

// Benchmark coverage calculation simulation
static void BM_CoverageCalculationSimulation(benchmark::State& state) {
    const int num_alignments = state.range(0);
    const int ref_length = 10000;
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> pos_dis(0, ref_length - 100);
    std::uniform_int_distribution<> len_dis(50, 200);
    
    for (auto _ : state) {
        std::priority_queue<CovEnd> coverage_ends;
        Coverage coverage;
        std::vector<std::pair<int, int>> alignments;
        
        // Generate random alignments
        for (int i = 0; i < num_alignments; ++i) {
            int start = pos_dis(gen);
            int length = len_dis(gen);
            alignments.push_back({start, start + length});
        }
        
        // Sort by start position (simulating sorted BAM)
        std::sort(alignments.begin(), alignments.end());
        
        // Simulate coverage calculation
        int last_pos = 0;
        size_t align_idx = 0;
        
        while (align_idx < alignments.size() || !coverage_ends.empty()) {
            int next_change = ref_length;
            
            if (align_idx < alignments.size()) {
                next_change = std::min(next_change, alignments[align_idx].first);
            }
            if (!coverage_ends.empty()) {
                next_change = std::min(next_change, coverage_ends.top().end);
            }
            
            // Process alignments starting at this position
            while (align_idx < alignments.size() && 
                   alignments[align_idx].first == next_change) {
                coverage_ends.push({alignments[align_idx].second, false});
                coverage.inc(false);
                align_idx++;
            }
            
            // Process alignments ending at this position
            while (!coverage_ends.empty() && 
                   coverage_ends.top().end == next_change) {
                coverage.dec(coverage_ends.top().rev);
                coverage_ends.pop();
            }
            
            last_pos = next_change;
        }
        
        benchmark::DoNotOptimize(coverage);
    }
}
BENCHMARK(BM_CoverageCalculationSimulation)->Range(100, 5000);

// Memory allocation benchmark
static void BM_MemoryAllocation(benchmark::State& state) {
    const int num_allocs = state.range(0);
    
    for (auto _ : state) {
        std::vector<Coverage> coverages;
        coverages.reserve(num_allocs);
        
        for (int i = 0; i < num_allocs; ++i) {
            coverages.emplace_back();
            coverages.back().f = i;
            coverages.back().r = i * 2;
        }
        
        benchmark::DoNotOptimize(coverages);
    }
}
BENCHMARK(BM_MemoryAllocation)->Range(1000, 100000);

BENCHMARK_MAIN();