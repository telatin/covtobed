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

typedef uint16_t DepthType; // type for depth of coverage, kept it small
const char ref_char = '>';

const string VERSION = "%prog 0.7"
	"\nCopyright (C) 2014-2019 Giovanni Birolo and Andrea Telatin\n"
	"License MIT"
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

		// open all files
		Input(const vector<string> &paths, const int q) : min_mapq(q) {
			if (paths.empty()) {
				// no input files, use standard input
				cerr << "Reading from STDIN... [try 'covtobed -h' for options]" << endl;
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
				debug cerr << "[M] " << alignment.Name << ":" << alignment.Position << " | Is mapped? " << alignment.IsMapped() << endl;
				more_alignments = input_bams.GetNextAlignmentCore(alignment);
				good_alignment = alignment.IsMapped() && alignment.MapQuality >= min_mapq;
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
};

// Class for input handling: reads and filters alignments from bams
class Output {
	public:
		enum Format {BED, COUNTS};

		// class constructor
		Output(ostream *o, Format f, bool s=false, int m=0, int x=100000, int l=1) : 
			out(o), format(f), strands(s), mincov(m), maxcov(x), minlen(l) {}

		// write interval to bed
		void operator() (const Interval &i, const Coverage &c) {
			// can the last interval be extended with the same coverage?
			if (i.ref == last_interval.ref && i.start == last_interval.end && last_coverage == c)
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
						if (c < 1000 and c >= 0) {
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
	parser.add_option("-v") .action("version") .help("prints program version");

	// input options
	parser.add_option("--physical-coverage").action("store_true").set_default("0").help("compute physical coverage (needs paired alignments in input)");
	parser.add_option("-q", "--min-mapq").metavar("MINQ").type("int").set_default("0").help("skip alignments whose mapping quality is less than MINQ (default: %default)");
	parser.add_option("-m", "--min-cov").metavar("MINCOV").type("int").set_default("0").help("print BED feature only if the coverage is bigger than (or equal to) MINCOV (default: %default)");
	parser.add_option("-x", "--max-cov").metavar("MAXCOV").type("int").set_default("100000").help("print BED feature only if the coverage is lower than MAXCOV (default: %default)");
	parser.add_option("-l", "--min-len").metavar("MINLEN").type("int").set_default("1").help("print BED feature only if its length is bigger (or equal to) than MINLELN (default: %default)");
	//parser.add_option("-f", "--flatten").action("store_true").set_default("0").help("Flatten adjacent BED features (usually when specifying -m/-x)");

	// output options
	parser.add_option("--output-strands").action("store_true").set_default("0").help("outputs coverage and stats separately for each strand");
	//parser.add_option("--target").action("store_true").set_default("0").help("outputs coverage and stats separately for each strand");
	vector<string> choices = {"bed", "counts"};
	parser.add_option("--format").choices(choices.begin(), choices.end()).set_default("bed").help("output format");

	// parse arguments
	optparse::Values options = parser.parse_args(argc, argv);

	const bool physical_coverage = options.get("physical_coverage");
	const int  minimum_coverage  = options.get("min_cov");
	const int  maximum_coverage  = options.get("max_cov");
	const int  minimum_length    = options.get("min_len");

	// get format
	Output::Format f;
	const string format_str = static_cast<const char *>(options.get("format"));

	if (format_str == "bed")
		f = Output::BED;
	else if (format_str == "counts")
		f = Output::COUNTS;
	else {
		parser.error("bad output format specification");
		return 1;
	}

	// open input and output
	try {

		Output output(&cout, f, options.get("output_strands"), minimum_coverage, maximum_coverage, minimum_length);
		Input input(parser.args(), options.get("min_mapq"));


		// main alignment parsing loop
		BamAlignment alignment;
		bool more_alignments = input.get_next_alignment(alignment);
		for (const auto &ref : input.get_ref_data()) { // loop on reference
			// init new reference data
			const int ref_id = input.get_ref_id(ref.RefName);
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
				output({ref.RefName, last_pos, next_change}, coverage);

				// increment coverage with alignments that start here
				while (more_alignments_for_ref && next_change == alignment.Position) {
					if (physical_coverage) {
						if (alignment.InsertSize > 0) {
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
			debug cerr << "[_] Completed at " << ref.RefName << ":" << last_pos << endl;
			// reference ended
			assert(coverage_ends.empty());
		}
		assert(!more_alignments);
	} catch (string msg) {
		parser.error(msg);
	}
	return 0;
}
