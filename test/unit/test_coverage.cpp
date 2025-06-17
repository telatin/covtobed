#include <catch2/catch_test_macros.hpp>

// Only include the classes we need for testing
#include <queue>
#include <vector>
#include <iostream>
#include <map>
#include <cstdlib>

using namespace std;

typedef uint32_t DepthType;

// Copy only the Coverage and CovEnd classes for unit testing
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

TEST_CASE("Coverage class basic operations", "[coverage]") {
    Coverage cov;
    
    SECTION("Initial state should be zero") {
        REQUIRE(cov.f == 0);
        REQUIRE(cov.r == 0);
        REQUIRE(static_cast<DepthType>(cov) == 0);
    }
    
    SECTION("Forward strand increment") {
        cov.inc(false);
        REQUIRE(cov.f == 1);
        REQUIRE(cov.r == 0);
        REQUIRE(static_cast<DepthType>(cov) == 1);
        
        // Multiple increments
        cov.inc(false);
        cov.inc(false);
        REQUIRE(cov.f == 3);
        REQUIRE(static_cast<DepthType>(cov) == 3);
    }
    
    SECTION("Reverse strand increment") {
        cov.inc(true);
        REQUIRE(cov.f == 0);
        REQUIRE(cov.r == 1);
        REQUIRE(static_cast<DepthType>(cov) == 1);
        
        // Multiple increments
        cov.inc(true);
        cov.inc(true);
        REQUIRE(cov.r == 3);
        REQUIRE(static_cast<DepthType>(cov) == 3);
    }
    
    SECTION("Mixed strand operations") {
        cov.inc(false);  // f=1, r=0
        cov.inc(true);   // f=1, r=1
        cov.inc(false);  // f=2, r=1
        
        REQUIRE(cov.f == 2);
        REQUIRE(cov.r == 1);
        REQUIRE(static_cast<DepthType>(cov) == 3);
    }
    
    SECTION("Decrement operations") {
        // Setup initial coverage
        cov.inc(false);
        cov.inc(false);
        cov.inc(true);
        REQUIRE(static_cast<DepthType>(cov) == 3);
        
        // Test decrements
        cov.dec(false);
        REQUIRE(cov.f == 1);
        REQUIRE(cov.r == 1);
        REQUIRE(static_cast<DepthType>(cov) == 2);
        
        cov.dec(true);
        REQUIRE(cov.f == 1);
        REQUIRE(cov.r == 0);
        REQUIRE(static_cast<DepthType>(cov) == 1);
    }
}

TEST_CASE("Coverage equality comparisons", "[coverage]") {
    Coverage cov1, cov2;
    
    SECTION("Equal empty coverage") {
        REQUIRE(cov1.equal(cov2, false) == true);
        REQUIRE(cov1.equal(cov2, true) == true);
    }
    
    SECTION("Unstranded equality") {
        cov1.inc(false);  // f=1, r=0
        cov2.inc(true);   // f=0, r=1
        
        // Total coverage is same (1), should be equal in unstranded mode
        REQUIRE(cov1.equal(cov2, false) == true);
        // But different in stranded mode
        REQUIRE(cov1.equal(cov2, true) == false);
    }
    
    SECTION("Stranded equality") {
        cov1.inc(false);  // f=1, r=0
        cov1.inc(true);   // f=1, r=1
        
        cov2.inc(false);  // f=1, r=0
        cov2.inc(true);   // f=1, r=1
        
        REQUIRE(cov1.equal(cov2, true) == true);
        REQUIRE(cov1.equal(cov2, false) == true);
    }
    
    SECTION("Different strand distributions") {
        cov1.f = 3; cov1.r = 1;  // Total: 4
        cov2.f = 2; cov2.r = 2;  // Total: 4
        
        REQUIRE(cov1.equal(cov2, false) == true);   // Same total
        REQUIRE(cov1.equal(cov2, true) == false);   // Different distribution
    }
}

TEST_CASE("CovEnd priority queue ordering", "[covend]") {
    SECTION("Basic ordering") {
        CovEnd end1 = {100, false};
        CovEnd end2 = {200, false};
        CovEnd end3 = {150, true};
        
        // Priority queue uses operator< which returns end > o.end
        // So smaller positions have higher priority (come out first)
        // Test the actual operator< implementation: return end > o.end
        REQUIRE(!(end1 < end2));  // end1.end > end2.end: 100 > 200 = false, so end1 < end2 is false
        REQUIRE(end2 < end1);     // end2.end > end1.end: 200 > 100 = true, so end2 < end1 is true
        REQUIRE(!(end1 < end3));  // end1.end > end3.end: 100 > 150 = false, so end1 < end3 is false
        REQUIRE(end3 < end1);     // end3.end > end1.end: 150 > 100 = true, so end3 < end1 is true
    }
    
    SECTION("Priority queue behavior simulation") {
        std::priority_queue<CovEnd> pq;
        
        pq.push({300, false});
        pq.push({100, true});
        pq.push({200, false});
        
        // Should come out in ascending order of end position
        REQUIRE(pq.top().end == 100);
        pq.pop();
        REQUIRE(pq.top().end == 200);
        pq.pop();
        REQUIRE(pq.top().end == 300);
    }
}