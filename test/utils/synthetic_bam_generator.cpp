#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <api/BamWriter.h>
#include <api/BamAlignment.h>

using namespace BamTools;
using namespace std;

class SyntheticBAMGenerator {
private:
    mt19937 gen;
    
public:
    SyntheticBAMGenerator(unsigned int seed = 42) : gen(seed) {}
    
    // Generate BAM with specific coverage pattern
    bool generateCoveragePattern(const string& filename, 
                               const string& ref_name,
                               int ref_length,
                               const vector<pair<pair<int, int>, int>>& coverage_regions) {
        
        BamWriter writer;
        if (!writer.Open(filename)) {
            cerr << "Failed to open " << filename << " for writing" << endl;
            return false;
        }
        
        // Create reference data
        RefVector references;
        RefData ref_data;
        ref_data.RefName = ref_name;
        ref_data.RefLength = ref_length;
        references.push_back(ref_data);
        
        // Set header
        SamHeader header;
        header.Version = "1.4";
        header.SortOrder = "coordinate";
        
        writer.SetReferences(references);
        if (!writer.WriteHeader(header)) {
            cerr << "Failed to write header" << endl;
            return false;
        }
        
        // Generate alignments for each coverage region
        int read_id = 1;
        for (const auto& region : coverage_regions) {
            int start = region.first.first;
            int end = region.first.second;
            int coverage = region.second;
            int read_length = 100;
            
            // Generate overlapping reads to achieve desired coverage
            for (int i = 0; i < coverage; ++i) {
                BamAlignment alignment;
                alignment.Name = "read_" + to_string(read_id++);
                alignment.RefID = 0;  // First (and only) reference
                alignment.Position = start + (i * (end - start - read_length) / max(1, coverage - 1));
                alignment.MapQuality = 30;
                alignment.AlignmentFlag = 0;  // Not paired, mapped, forward strand
                alignment.InsertSize = 0;
                alignment.MateRefID = -1;
                alignment.MatePosition = -1;
                
                // Create CIGAR string (simple match)
                CigarOp cigar_op;
                cigar_op.Type = 'M';
                cigar_op.Length = read_length;
                alignment.CigarData.push_back(cigar_op);
                
                // Create dummy sequence and quality
                alignment.QueryBases = string(read_length, 'A');
                alignment.Qualities = string(read_length, '!');  // Quality 0
                
                if (!writer.WriteAlignment(alignment)) {
                    cerr << "Failed to write alignment" << endl;
                    return false;
                }
            }
        }
        
        writer.Close();
        return true;
    }
    
    // Generate stranded coverage data
    bool generateStrandedData(const string& filename,
                            const string& ref_name,
                            int ref_length) {
        BamWriter writer;
        if (!writer.Open(filename)) {
            return false;
        }
        
        RefVector references;
        RefData ref_data;
        ref_data.RefName = ref_name;
        ref_data.RefLength = ref_length;
        references.push_back(ref_data);
        
        SamHeader header;
        header.Version = "1.4";
        header.SortOrder = "coordinate";
        
        writer.SetReferences(references);
        writer.WriteHeader(header);
        
        // Generate reads on both strands
        int read_id = 1;
        
        // Forward strand reads (position 100-200, coverage 3)
        for (int i = 0; i < 3; ++i) {
            BamAlignment alignment;
            alignment.Name = "forward_" + to_string(read_id++);
            alignment.RefID = 0;
            alignment.Position = 100 + i * 10;
            alignment.MapQuality = 30;
            alignment.AlignmentFlag = 0;  // Forward strand
            
            CigarOp cigar_op;
            cigar_op.Type = 'M';
            cigar_op.Length = 50;
            alignment.CigarData.push_back(cigar_op);
            
            alignment.QueryBases = string(50, 'A');
            alignment.Qualities = string(50, '!');
            
            writer.WriteAlignment(alignment);
        }
        
        // Reverse strand reads (position 150-250, coverage 2)
        for (int i = 0; i < 2; ++i) {
            BamAlignment alignment;
            alignment.Name = "reverse_" + to_string(read_id++);
            alignment.RefID = 0;
            alignment.Position = 150 + i * 10;
            alignment.MapQuality = 30;
            alignment.AlignmentFlag = 16;  // Reverse strand
            
            CigarOp cigar_op;
            cigar_op.Type = 'M';
            cigar_op.Length = 50;
            alignment.CigarData.push_back(cigar_op);
            
            alignment.QueryBases = string(50, 'T');
            alignment.Qualities = string(50, '!');
            
            writer.WriteAlignment(alignment);
        }
        
        writer.Close();
        return true;
    }
    
    // Generate physical coverage data (paired-end)
    bool generatePhysicalCoverageData(const string& filename,
                                    const string& ref_name,
                                    int ref_length) {
        BamWriter writer;
        if (!writer.Open(filename)) {
            return false;
        }
        
        RefVector references;
        RefData ref_data;
        ref_data.RefName = ref_name;
        ref_data.RefLength = ref_length;
        references.push_back(ref_data);
        
        SamHeader header;
        header.Version = "1.4";
        header.SortOrder = "coordinate";
        
        writer.SetReferences(references);
        writer.WriteHeader(header);
        
        // Generate paired-end reads
        int read_id = 1;
        uniform_int_distribution<> pos_dist(100, ref_length - 500);
        uniform_int_distribution<> insert_dist(200, 400);
        
        for (int i = 0; i < 50; ++i) {
            int start_pos = pos_dist(gen);
            int insert_size = insert_dist(gen);
            int mate_pos = start_pos + insert_size - 100;
            
            // First read in pair
            BamAlignment read1;
            read1.Name = "pair_" + to_string(read_id);
            read1.RefID = 0;
            read1.Position = start_pos;
            read1.MapQuality = 30;
            read1.AlignmentFlag = 99;  // Paired, proper pair, first in pair
            read1.InsertSize = insert_size;
            read1.MateRefID = 0;
            read1.MatePosition = mate_pos;
            
            CigarOp cigar_op;
            cigar_op.Type = 'M';
            cigar_op.Length = 100;
            read1.CigarData.push_back(cigar_op);
            
            read1.QueryBases = string(100, 'A');
            read1.Qualities = string(100, '!');
            
            // Second read in pair
            BamAlignment read2;
            read2.Name = "pair_" + to_string(read_id);
            read2.RefID = 0;
            read2.Position = mate_pos;
            read2.MapQuality = 30;
            read2.AlignmentFlag = 147;  // Paired, proper pair, second in pair, reverse
            read2.InsertSize = -insert_size;
            read2.MateRefID = 0;
            read2.MatePosition = start_pos;
            
            read2.CigarData.push_back(cigar_op);
            read2.QueryBases = string(100, 'T');
            read2.Qualities = string(100, '!');
            
            writer.WriteAlignment(read1);
            writer.WriteAlignment(read2);
            
            read_id++;
        }
        
        writer.Close();
        return true;
    }
    
    // Generate edge case data
    bool generateEdgeCaseData(const string& filename,
                            const string& ref_name,
                            int ref_length,
                            const string& case_type) {
        BamWriter writer;
        if (!writer.Open(filename)) {
            return false;
        }
        
        RefVector references;
        RefData ref_data;
        ref_data.RefName = ref_name;
        ref_data.RefLength = ref_length;
        references.push_back(ref_data);
        
        SamHeader header;
        header.Version = "1.4";
        header.SortOrder = "coordinate";
        
        writer.SetReferences(references);
        writer.WriteHeader(header);
        
        if (case_type == "high_coverage") {
            // Generate very high coverage region
            for (int i = 0; i < 1000; ++i) {
                BamAlignment alignment;
                alignment.Name = "high_cov_" + to_string(i);
                alignment.RefID = 0;
                alignment.Position = 1000 + (i % 100);
                alignment.MapQuality = 30;
                alignment.AlignmentFlag = 0;
                
                CigarOp cigar_op;
                cigar_op.Type = 'M';
                cigar_op.Length = 100;
                alignment.CigarData.push_back(cigar_op);
                
                alignment.QueryBases = string(100, 'A');
                alignment.Qualities = string(100, '!');
                
                writer.WriteAlignment(alignment);
            }
        } else if (case_type == "single_read") {
            // Single read only
            BamAlignment alignment;
            alignment.Name = "single_read";
            alignment.RefID = 0;
            alignment.Position = 1000;
            alignment.MapQuality = 30;
            alignment.AlignmentFlag = 0;
            
            CigarOp cigar_op;
            cigar_op.Type = 'M';
            cigar_op.Length = 100;
            alignment.CigarData.push_back(cigar_op);
            
            alignment.QueryBases = string(100, 'A');
            alignment.Qualities = string(100, '!');
            
            writer.WriteAlignment(alignment);
        }
        
        writer.Close();
        return true;
    }
};

// Command line interface
int main(int argc, char* argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <output_file> <pattern_type> [options]" << endl;
        cout << "Pattern types:" << endl;
        cout << "  coverage <start> <end> <depth> - Generate specific coverage pattern" << endl;
        cout << "  stranded - Generate stranded coverage data" << endl;
        cout << "  physical - Generate physical coverage data" << endl;
        cout << "  edge <case_type> - Generate edge case data" << endl;
        return 1;
    }
    
    string output_file = argv[1];
    string pattern_type = argv[2];
    
    SyntheticBAMGenerator generator;
    
    if (pattern_type == "coverage" && argc >= 6) {
        int start = stoi(argv[3]);
        int end = stoi(argv[4]);
        int depth = stoi(argv[5]);
        
        vector<pair<pair<int, int>, int>> regions;
        regions.push_back({{start, end}, depth});
        
        if (generator.generateCoveragePattern(output_file, "chr1", 10000, regions)) {
            cout << "Generated coverage pattern BAM: " << output_file << endl;
        } else {
            cerr << "Failed to generate BAM file" << endl;
            return 1;
        }
    } else if (pattern_type == "stranded") {
        if (generator.generateStrandedData(output_file, "chr1", 10000)) {
            cout << "Generated stranded BAM: " << output_file << endl;
        } else {
            cerr << "Failed to generate BAM file" << endl;
            return 1;
        }
    } else if (pattern_type == "physical") {
        if (generator.generatePhysicalCoverageData(output_file, "chr1", 10000)) {
            cout << "Generated physical coverage BAM: " << output_file << endl;
        } else {
            cerr << "Failed to generate BAM file" << endl;
            return 1;
        }
    } else if (pattern_type == "edge" && argc >= 4) {
        string case_type = argv[3];
        if (generator.generateEdgeCaseData(output_file, "chr1", 10000, case_type)) {
            cout << "Generated edge case BAM: " << output_file << endl;
        } else {
            cerr << "Failed to generate BAM file" << endl;
            return 1;
        }
    } else {
        cerr << "Invalid arguments" << endl;
        return 1;
    }
    
    return 0;
}