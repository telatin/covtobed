#include <catch2/catch_test_macros.hpp>
#include "../../interval.h"

TEST_CASE("CoordinateInterval basic operations", "[interval]") {
    SECTION("Empty interval detection") {
        CoordinateInterval empty1 = {10, 10};
        CoordinateInterval empty2 = {10, 5};  // end < start
        CoordinateInterval valid = {10, 20};
        
        REQUIRE(empty1.empty() == true);
        REQUIRE(empty2.empty() == true);
        REQUIRE(valid.empty() == false);
        
        // Boolean conversion
        REQUIRE(!empty1);
        REQUIRE(!empty2);
        REQUIRE(valid);
    }
    
    SECTION("Length calculation") {
        CoordinateInterval interval1 = {10, 20};
        CoordinateInterval interval2 = {100, 150};
        CoordinateInterval empty = {10, 10};
        
        REQUIRE(interval1.length() == 10);
        REQUIRE(interval2.length() == 50);
        REQUIRE(empty.length() == 0);
    }
    
    SECTION("Ordering comparison") {
        CoordinateInterval a = {10, 20};
        CoordinateInterval b = {30, 40};
        CoordinateInterval c = {10, 30};  // Same start, different end
        
        REQUIRE(a < b);      // 10 < 30
        REQUIRE(a < c);      // Same start (10), but 20 < 30
        REQUIRE(!(b < a));   // 30 > 10
    }
    
    SECTION("Disjointness check") {
        CoordinateInterval a = {10, 20};
        CoordinateInterval b = {20, 30};  // Adjacent (touching)
        CoordinateInterval c = {15, 25};  // Actually overlapping with a
        CoordinateInterval d = {25, 35};  // Disjoint from a
        
        REQUIRE(a << b);     // a.end <= b.start: 20 <= 20, true (disjoint/adjacent)
        REQUIRE(!(a << c));  // a.end <= c.start: 20 <= 15, false (not disjoint, overlapping)
        REQUIRE(a << d);     // a.end <= d.start: 20 <= 25, true (disjoint)
    }
    
    SECTION("Intersection operation") {
        CoordinateInterval a = {10, 30};
        CoordinateInterval b = {20, 40};
        CoordinateInterval c = {35, 45};  // No overlap
        
        auto intersection1 = a * b;
        REQUIRE(intersection1.start == 20);
        REQUIRE(intersection1.end == 30);
        REQUIRE(intersection1.length() == 10);
        
        auto intersection2 = a * c;
        REQUIRE(intersection2.empty());  // No overlap
    }
    
    SECTION("Difference operation") {
        CoordinateInterval a = {10, 50};
        CoordinateInterval b = {20, 30};  // Contained within a
        
        auto diff = a - b;
        
        // Should return two intervals: [10,20) and [30,50)
        REQUIRE(diff.first.start == 10);
        REQUIRE(diff.first.end == 20);
        REQUIRE(diff.second.start == 30);
        REQUIRE(diff.second.end == 50);
    }
}

TEST_CASE("Interval with reference names", "[interval]") {
    SECTION("Basic interval creation") {
        Interval interval("chr1", 100, 200);
        
        REQUIRE(interval.ref == "chr1");
        REQUIRE(interval.start == 100);
        REQUIRE(interval.end == 200);
        REQUIRE(interval.length() == 100);
    }
    
    SECTION("Interval ordering with different references") {
        Interval a("chr1", 100, 200);
        Interval b("chr2", 50, 150);   // Different chromosome
        Interval c("chr1", 150, 250);  // Same chromosome, later position
        
        // Should order by reference name first
        REQUIRE(a < b);  // "chr1" < "chr2"
        REQUIRE(a < c);  // Same ref, 100 < 150
        REQUIRE(!(b < a));
    }
    
    SECTION("Interval with same reference") {
        Interval a("chr1", 100, 200);
        Interval b("chr1", 150, 250);
        Interval c("chr1", 100, 150);  // Same start, different end
        
        REQUIRE(a < b);  // 100 < 150
        REQUIRE(c < a);  // Same start (100), end comparison: 150 < 200, so c < a
        REQUIRE(!(a < c));  // Same start (100), end comparison: 200 > 150, so a is NOT < c
    }
}

TEST_CASE("NamedInterval operations", "[named_interval]") {
    SECTION("Named interval creation") {
        NamedInterval interval;
        interval.ref = "chr1";
        interval.start = 100;
        interval.end = 200;
        interval.name = "region1";
        
        REQUIRE(interval.ref == "chr1");
        REQUIRE(interval.name == "region1");
        REQUIRE(interval.length() == 100);
    }
    
    SECTION("Empty name handling") {
        NamedInterval interval;
        interval.ref = "chr1";
        interval.start = 100;
        interval.end = 200;
        interval.name = "";
        
        REQUIRE(interval.name.empty());
    }
}