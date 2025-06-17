#include <queue>
#include <vector>
#include <iostream>
#include <map>
#include <cstdlib>
#include <cassert>
#include <api/BamMultiReader.h>
#include <api/BamAlignment.h>
#include "OptionParser.h"
#include "interval.h"


using namespace BamTools;
using namespace std;

typedef uint32_t DepthType; // type for depth of coverage, kept it small
const char ref_char = '>';  // reference prefix for "counts" output

const string VERSION = "%prog 1.4.0"
	"\nCopyright (C) 2014-2019 Giovanni Birolo and Andrea Telatin\n"
	"https://github.com/telatin/covtobed - License MIT"
	".\n"
	"This is free software: you are free to change and redistribute it.\n"
	"There is NO WARRANTY, to the extent permitted by law.";

#define debug if(false)

// Class that stores info about the end position of alignments, used in the alignment queue
struct CovEnd {
	PositionType end;
	bool rev;
	// order for queuing
	bool operator<(const CovEnd &o) const {
		return end > o.end;
	}
};

// Class for input handling: reads and filters alignments from bams
class Input {
	public:
		BamMultiReader input_bams;
		const int min_mapq;
		const int discard_invalid_alignments;
		// open all files
		Input(const vector<string> &paths, const int q, const int v) : min_mapq(q), discard_invalid_alignments(v) {
			if (paths.empty()) {
				// no input files, use standard input
				// 1.3.0 - provide feedback unless $COVTOBED_QUIET is set to 1
				// Check environment variable COVTOBED_QUIET
				if (getenv("COVTOBED_QUIET") == NULL) {
					cerr << "Reading from STDIN... [Ctrl+C to exit; 'covtobed -h' for help]" << endl;
				}
				
				if (!input_bams.OpenFile("-"))
					throw string("cannot read BAM from standard input, are you piping a BAM file?");
					//throw input_bams.GetErrorString();
			} else {
				for (const auto &path : paths)
					if (!input_bams.OpenFile(path))
						throw "cannot open '" + path + "'. File not found, or not a BAM file.";
						//throw input_bams.GetErrorString();
			}
		}

		// get next good alignment (if any)
		bool get_next_alignment(BamAlignment & alignment) {
			bool more_alignments, good_alignment;
			do {
				debug cerr << "[M] Read  on ref#" << alignment.RefID << " pos:" << alignment.Position << 
				    "\n\t| Is mapped? " << alignment.IsMapped() << " | AlignmentFlag:" << alignment.AlignmentFlag << endl;
				more_alignments = input_bams.GetNextAlignmentCore(alignment);
				if (discard_invalid_alignments) {
					good_alignment = alignment.IsMapped() && alignment.MapQuality >= min_mapq
									&& !alignment.IsDuplicate() && !alignment.IsFailedQC() && alignment.IsPrimaryAlignment();

					// 1.2.0 ProperPair
					if (alignment.IsPaired() && ! alignment.IsProperPair() ) {
					    good_alignment = false;
					}
				} else {
					good_alignment = alignment.IsMapped() && alignment.MapQuality >= min_mapq;
				}
			} while (more_alignments && !good_alignment);
			return more_alignments;
		}
		vector<RefData> get_ref_data() { return input_bams.GetReferenceData(); }
		int get_ref_id(const string &ref) { return input_bams.GetReferenceID(ref); }
};

// Class that stores info about the current coverage
struct Coverage {
	DepthType f = 0, r = 0;

	operator DepthType() const { return f + r; }

	void inc(bool rev=false) {
		if (rev)
			++r;
		else
			++f;
	}
	void dec(bool rev=false) {
		if (rev)
			--r;
		else
			--f;
	}
	bool equal(const Coverage &o, bool stranded) const {
	    	if (stranded)
		    return f == o.f && r == o.r;
		else
		    return f + r == o.f + o.r;
	}
};

// Class for output handling: writes coverage in the specified format
class Output {
	public:
		enum Format {BED, COUNTS};

		// class constructor
		Output(ostream *o, const char *f, bool s=false, int m=0, int x=100000, int l=1) : out(o), format(parse_format(f)), strands(s), mincov(m), maxcov(x), minlen(l) {
		}

		// write interval to bed
		void operator() (const Interval &i, const Coverage &c) {
			// can the last interval be extended with the same coverage?
			if (i.ref == last_interval.ref && i.start == last_interval.end && last_coverage.equal(c, strands))
				// extend previous interval
				last_interval.end = i.end;
			else {
				// output previous interval
				write(last_interval, last_coverage);
				if (i.ref != last_interval.ref)
					// new reference
					switch(format) {
						case BED:
							break;
						case COUNTS:
							*out << ref_char << i.ref << endl;
							break;
					}
				last_interval = i;
				last_coverage = c;
			}
		}
		~Output() {
			write(last_interval, last_coverage);
		}
	private:
		void write(const Interval &i, const Coverage &c) {
			if (i and c >= mincov and c < maxcov and (i.end - i.start >= minlen) )  { // interval not empty plus user constraints
				switch(format) {
					case Format::BED:
						*out << i.ref << '\t' << i.start << '\t' << i.end << '\t';
						if (c >= 0) {
							write_coverage(c);
						}
						break;
					case Format::COUNTS:
						write_coverage(c);
						*out << '\t' << i.length();
						break;
				}
				*out << endl;
			}
		}
		void write_coverage(const Coverage &c) {
			if (strands)
				*out << c.f << '\t' << c.r;
			else
				*out  << static_cast<DepthType>(c);
		}
		static Format parse_format(const char *format_str) {
			const string s = format_str;
			if (s == "bed")
				return Output::BED;
			if (s == "counts")
				return Output::COUNTS;
			throw string("unkown format specification: \"") + format_str + "\"";
		}


		ostream *out;
		const Format format;
		const bool strands;
		const int mincov;
		const int maxcov;
		const int minlen;
		Interval last_interval;
		Coverage last_coverage;
};

int main(int argc, char *argv[]) {
	// general options

	optparse::OptionParser parser = optparse::OptionParser().description("Computes coverage from alignments").usage("%prog [options] [BAM]...").version(VERSION);
	//parser.add_option("-v") .action("version") .help("prints program version");

	// input options
	parser.add_option("--physical-coverage").action("store_true").set_default("0").help("compute physical coverage (needs paired alignments in input)");
	parser.add_option("-q", "--min-mapq").metavar("MINQ").type("int").set_default("0").help("skip alignments whose mapping quality is less than MINQ (default: %default)");
	parser.add_option("-m", "--min-cov").metavar("MINCOV").type("int").set_default("0").help("print BED feature only if the coverage is bigger than (or equal to) MINCOV (default: %default)");
	parser.add_option("-x", "--max-cov").metavar("MAXCOV").type("int").set_default("100000").help("print BED feature only if the coverage is lower than MAXCOV (default: %default)");
	parser.add_option("-l", "--min-len").metavar("MINLEN").type("int").set_default("1").help("print BED feature only if its length is bigger (or equal to) than MINLELN (default: %default)");
	parser.add_option("-z", "--min-ctg-len").metavar("MINCTGLEN").type("int").help("Skip reference sequences (contigs) shorter than this value");
	parser.add_option("-d", "--discard-invalid-alignments").action("store_true").set_default("0").help("skip duplicates, failed QC, and non primary alignment, minq>0 (or user-defined if higher) (default: %default)");
	parser.add_option("--keep-invalid-alignments").action("store_true").set_default("0").help("Keep duplicates, failed QC, and non primary alignment, min=0 (or user-defined if higher) (default: %default)");

	// output options
	parser.add_option("--output-strands").action("store_true").set_default("0").help("output coverage and stats separately for each strand");
	vector<string> choices = {"bed", "counts"};
	parser.add_option("--format").choices(choices.begin(), choices.end()).set_default("bed").help("output format");

	// parse arguments
	optparse::Values options = parser.parse_args(argc, argv);

	const bool physical_coverage = options.get("physical_coverage");
	bool only_valid              = options.get("discard_invalid_alignments"); 
	const bool keep_invalid      = options.get("keep_invalid_alignments"); 
	const int  minimum_coverage  = options.get("min_cov");
	const int  maximum_coverage  = options.get("max_cov");
	const int  minimum_length    = options.get("min_len");
	const int  minimum_contig_len= options.get("min_ctg_len");
	int min_mapq                 = options.get("min_mapq");

	// since 1.3.2 deprecate --discard-invalid-alignments

	if (only_valid and keep_invalid) {
		cerr << "ERROR: --discard-invalid-alignments and --keep-invalid-alignments are incompatible." << endl << "In the future --discard-invalid-alignments will be the default." << endl;
		exit(1);
	} else if (!only_valid and !keep_invalid) {
		if (getenv("COVTOBED_QUIET") == NULL) {
			cerr << "WARNING: --discard-invalid-alignments in the future will be activated by default." << endl;
		}
		// only_valid = false by default as legacy behaviour
	} else if (only_valid) {
		if (getenv("COVTOBED_QUIET") == NULL) {
			cerr << "WARNING: --discard-invalid-alignments will be automatically enabled in the future." << endl;
		}
		only_valid = true;
	} else if (keep_invalid) {
		// this if statement will remain: in the future only_valid will be the default
		only_valid = false;
	}

	if (only_valid and !min_mapq) {
		min_mapq = 1;
	} else {
		min_mapq = options.get("min_mapq");
	}


	try {
		// open input and output
		Input input(parser.args(), min_mapq, only_valid);
		Output output(&cout,
			static_cast<const char *>(options.get("format")), 
			static_cast<bool>(options.get("output_strands")), 
			minimum_coverage, maximum_coverage, minimum_length);


		// main alignment parsing loop
		BamAlignment alignment;
		bool more_alignments = input.get_next_alignment(alignment);
		for (const auto &ref : input.get_ref_data()) { // loop on reference
			// init new reference data
			const int ref_id = input.get_ref_id(ref.RefName);
			if (ref.RefLength <= minimum_contig_len) {
				continue;
			}
			debug cerr << "[R] Reference: " << ref_id << endl;
			PositionType last_pos = 0;
			priority_queue<CovEnd> coverage_ends;
			Coverage coverage;
			bool more_alignments_for_ref = more_alignments && alignment.RefID == ref_id;

			// loop current reference
			do {
				// find next possible coverage change
				const PositionType next_change = more_alignments_for_ref ? 
					(coverage_ends.empty() ? alignment.Position : min(alignment.Position, coverage_ends.top().end)) :
					(coverage_ends.empty() ? ref.RefLength : coverage_ends.top().end);
				debug cerr << "[-] Coverage is " << coverage << " up to " << next_change << endl;

				// output coverage
				assert(coverage_ends.size() == coverage);

				// check unsorted 1.3.4
				if  (  last_pos > next_change ) {
					throw string("Position going backward, is the BAM sorted? last_pos=" + 
						to_string(last_pos) + " next_change=" + to_string(next_change));
				}
				output({ref.RefName, last_pos, next_change}, coverage);
				

				// increment coverage with alignments that start here
				while (more_alignments_for_ref && next_change == alignment.Position) {
					if (physical_coverage) {
						if (alignment.InsertSize > 0) {
						        debug cerr << "   [phy] pos:" << alignment.Position << " size:" << alignment.InsertSize << endl;
							coverage_ends.push({alignment.Position + alignment.InsertSize, alignment.IsReverseStrand()});
							coverage.inc(alignment.IsReverseStrand());
						}
					} else {
						coverage_ends.push({alignment.GetEndPosition(), alignment.IsReverseStrand()});
						coverage.inc(alignment.IsReverseStrand());
					}
					more_alignments = input.get_next_alignment(alignment);
					more_alignments_for_ref = more_alignments && alignment.RefID == ref_id;
				}
				// decrement coverage with alignments that end here
				while (!coverage_ends.empty() && next_change == coverage_ends.top().end) {
					coverage.dec(coverage_ends.top().rev);
					coverage_ends.pop();
				}

				debug cerr << "[<] Coverage is " << coverage << " from " << next_change << endl;
				last_pos = next_change;
				

			} while (last_pos != ref.RefLength);
			debug cerr << "[_] Completed at " << ref.RefName << ":" << last_pos << " [coverage:"  << coverage_ends.size() << "]" << endl;
			// reference ended
			if (!coverage_ends.empty()) {
			    cerr << "Coverage is not zero at the end of " << ref.RefName << endl;
			    cerr << "Try samtools fixmate on the input file" << endl;
			}
			// 1.2.0 -- removed: assert(coverage_ends.empty());
		
		}
		//assert(!more_alignments);
		if (more_alignments) {
			throw string("Unexpected alignment found, is the BAM sorted?");			
		}
	} catch (const string &msg) {
		parser.error(msg);
	}
	return 0;
}
