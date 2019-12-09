#ifndef __INTERVAL_H__
#define __INTERVAL_H__
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <algorithm>
//#include <cstdlib>
#include <cassert>

using namespace std;

// types
typedef int PositionType;
//typedef uint16_t DepthType;
//typedef size_t CountType;

struct CoordinateInterval {
	PositionType start, end;
	bool empty() const { return end <= start; }
	size_t length() const { return empty() ? 0 : end - start; }
	operator bool() const { return !empty(); }
	// ordering
	bool operator<(const CoordinateInterval &i) const { return start < i.start || (start == i.start && end < i.end); }
	// disjointness
	bool operator<<(const CoordinateInterval &i) const { return end <= i.start; }
	// intersection
	CoordinateInterval operator*(const CoordinateInterval &o) const { return {max(start, o.start), min(end, o.end)}; }
	/*
	tuple<CoordinateInterval, CoordinateInterval, CoordinateInterval> split(const CoordinateInterval &o) const {
		return {{start, min(end, o.start)}, 
		}
		*/

	// difference, returns two (possibly empty) intervals
	pair<CoordinateInterval, CoordinateInterval> operator-(const CoordinateInterval &o) const {
		assert(!o.empty());

		//debug_inter cerr << *this << " - " << o << " = " << CoordinateInterval{start, min(end, o.start)} << ", " << CoordinateInterval{max(start, o.end), end} << endl;
		return {{start, min(end, o.start)}, {max(start, o.end), end}};
	}
};

ostream &operator<<(ostream &o, const CoordinateInterval &i) { return o << i.start << '-' << i.end; }

bool compare_ref_name(const string &r1, const string &r2) {
	//int prefix_len;
	//for (prefix_len = 0; prefix_len < min(r1.size(), r2.size()); ++prefix_len) {
	//FIXME implement
	return r1 < r2;
}

struct Interval : public CoordinateInterval {
	string ref;

	Interval(const string &r="", PositionType s=0, PositionType e=0) : CoordinateInterval{s, e}, ref(r) {}
	//Interval(istream &s) { read_bed(s); }
	// ordering
	bool operator<(const Interval &o) const { 
		return compare_ref_name(ref, o.ref) || 
			(ref == o.ref && CoordinateInterval::operator<(o));
	}
	
	ostream &print_bed(ostream &s) const { return s << ref << '\t' << start << '\t' << end; }
	istream &read_bed(istream &s) { return s >> ref >> start >> end; }
	ostream &print_igv(ostream &s) const { return s << ref << ':' << start << '-' << end; }
};
ostream &operator<< (ostream &s, const Interval &i) { return i.print_bed(s); }
istream &operator>> (istream &s, Interval &i) { return i.read_bed(s); }


struct NamedInterval : public Interval {
	string name;

	//Interval(const string &r, const CoordinateInterval &i, const string &n="") : CoordinateInterval(i), ref(r), name(n) {}

	/*
	ostream &write_bed(ostream &o, const Coverage &c) const {
		o << ref << '\t' << start << '\t' << end << '\t' << (name.empty() ? "." : name) << '\t';
		if (output_strands)
			o << c.forward_coverage() << '\t' << c.reverse_coverage();
		else
			o << c.coverage();
		return o << endl;
	}
	ostream &write_bed(ostream &o, const int score) const {
		return o << ref << '\t' << start << '\t' << end << '\t' << (name.empty() ? "." : name) << '\t' << score << endl;
	}
	ostream &write_bed(ostream &o) const {
		o << ref << '\t' << start << '\t' << end;
		if (!name.empty())
			o << '\t' << name;
		return o << endl;
	}

	ostream &write_bed(ostream &o, const string &cols) const {
	o << ref << '\t' << start << '\t' << end << '\t' << (name.empty() ? "." : name);
	if (!cols.empty())
	o << '\t' << cols;
	return o << endl;
	}
	 */
	//friend ostream &operator<<(ostream &o, const Interval &i) {
	//	return o << i.ref << ':' << i.start << '-' << i.end;
	//}
	ostream &print_bed(ostream &s) const { return s << ref << '\t' << start << '\t' << end << '\t' << (name.empty() ? "." : name); }
	istream &read_bed(istream &s) { return s >> ref >> start >> end >> name; }
};

/*
template <class I>
ostream &operator<< (ostream &s, const I &i) { return i.print_bed(s); }
template <class I>
istream &operator>> (istream &s, I &i) { return i.read_bed(s); }
*/
ostream &operator<< (ostream &s, const NamedInterval &i) { return i.print_bed(s); }
istream &operator>> (istream &s, NamedInterval &i) { return i.read_bed(s); }

class Intervals {
	public:
		Intervals() {}

		Intervals(istream &i) {
			read_bed(i);
			sort();
		}

		// load intervals from bed file
		bool read_bed(istream &input) {
			string line;
			for (int line_count = 0; getline(input, line); ++line_count) {
				//debug_bed cerr << "read line " << line << endl;
				istringstream line_input(line);
				string ref;
				NamedInterval i;
				//if (line_input >> ref >> i.start >> i.end) {
				if (line_input >> i) {
					if (line_input >> i.name) {
						if (i.name == ".")
							i.name.clear();
						if (!i.name.empty())
							_has_names = true;
					}
					intervals_by_ref[ref].push_back(i);
					//debug_bed cerr << "read interval '" << i.name << "' " << ref << ":" << i.start << "-" << i.end << endl;
				} else {
					ostringstream err;
					cerr << "bad format at line " << line_count << " of target file";
				}
				input.clear();
				//input.ignore(10000, '\n');

			}
			if (!input.eof()) {
				cerr << "error reading target file" << endl;
				exit(1);
			}
			//debug cerr << "loaded " << count() << " intervals in " << intervals_by_ref.size() << " references" << endl;

			return true;
		}

		void sort() {
			for (auto &intervals : intervals_by_ref)
				std::sort(intervals.second.begin(), intervals.second.end());
		}


		bool has_names() const { return _has_names; }

		size_t count() const {
			size_t acc = 0;
			for (const auto &ref_intervals : intervals_by_ref)
				acc += ref_intervals.second.size();
			return acc;
		}



		// we are assuming that coverage_interval is more to the right than in previous calls
		/*vector<Interval> intersections(const Interval &o) {
			vector<Interval> intersections;


			return intersections;
		}*/

		/*
		vector<Interval> update(Coverage<DepthType> coverage, const Interval &coverage_interval) {
			vector<Interval> intersections;
			if (intervals_by_ref.empty()) {
				// empty target behaves as whole genome target
				intersections.push_back(coverage_interval);
			} else {
				// reset starting interval if the reference changed
				if (coverage_interval.ref != last_ref) {
					last_first_interval = 0;
					last_ref = coverage_interval.ref;
				}
				auto &intervals = intervals_by_ref[last_ref];

				// go backward if we are too advanced, not needed if coverage_interval is increasing
				*/ /* FIXME this does not work anyway, because intervals are not sorted by end coordinate
				   should remember max interval length, then we can find a bound
				   while (last_first_interval > 0 && !(intervals[last_first_interval] << coverage_interval)) {
				   debug_inter cerr << "go back interval " << ref << ":" << intervals[last_first_interval] << endl;
				   --last_first_interval;
				   }
				 */ /*

				debug_inter cerr << "looking for interval " << coverage_interval << endl;

				// advance target interval in the given ref until we hit coverage_interval
				while (last_first_interval < intervals.size() && intervals[last_first_interval] << coverage_interval) {
					debug_inter cerr << "advance interval " << last_ref << ":" << intervals[last_first_interval] << endl;
					++last_first_interval;
				}

				CoordinateInterval out_of_target_right = coverage_interval;
				// advance target interval in the given ref until we do not hit coverage_interval anymore
				for (size_t i = last_first_interval; i < intervals.size() && !(coverage_interval << intervals[i]); ++i) {
					CoordinateInterval out_of_target_left;
					// the left interval is sure to be outside the target since we are proceeding rightwards
					debug_inter cerr << "intersecting interval " << last_ref << ":" << intervals[i] << endl;
					tie(out_of_target_left, out_of_target_right) = out_of_target_right - intervals[i];
					if (out_of_target_left) {
						outside_stats.add_coverage(coverage, out_of_target_left);
						debug_inter cerr << "outside interval (left) " << last_ref << ":" << out_of_target_left << endl;
					}

					const auto intersection = coverage_interval*intervals[i];
					debug_inter cerr << "intersection " << last_ref << ":" << intersection << endl;
					intervals[i].stats.add_coverage(coverage, intersection);
					intersections.emplace_back(last_ref, intersection, intervals[i].name);
				}

				if (out_of_target_right) {
					outside_stats.add_coverage(coverage, out_of_target_right);
					debug_inter cerr << "outside interval (right) " << last_ref << ":" << out_of_target_right << endl;
				}
			}

			return intersections;
		}
		*/
	private:
		map<string, vector<CoordinateInterval> > intervals_by_ref;
		string last_ref;
		bool _has_names = false;
		size_t last_first_interval = 0;
};


class BEDOutput {
	public:
		BEDOutput(const char *p) : path(p), use_stdout(path == "-") {
			if (!use_stdout)
				file_out.open(p, ofstream::out);
			last_interval.start = last_interval.end = 0;
		}
		// write interval to bed
		// checks if we can extend the last interval
		// cols is appended after the interval name if not empty
		void operator() (const NamedInterval &i, const string &cols="") {
			if (writable()) {
				// check for extension
				if (i.ref == last_interval.ref && i.start == last_interval.end && 
						i.name == last_interval.name && cols == last_cols)
					// extend previous interval
					last_interval.end = i.end;
				else {
					// output previous interval
					if (!last_interval.ref.empty())
						write_bed(last_interval, last_cols);
					last_interval = i;
					last_cols = cols;
				}
			}
		}
		~BEDOutput() {
			if (writable()) {
				if (!last_interval.ref.empty())
					write_bed(last_interval, last_cols);
			}
		}
	private:
		bool writable() const {
			return use_stdout || file_out.is_open();
		}
		void write_bed(const NamedInterval &i, const string &cols="") {
			ostream *out = use_stdout ? &cout : &file_out;
			*out << i.ref << '\t' << i.start << '\t' << i.end;
			if (!cols.empty())
				*out << '\t' << (i.name.empty() ? "." : i.name) << '\t' << cols;
			else if (!i.name.empty())
				*out << '\t' << i.name;
			*out << endl;
		}
		
		const string path;
		const bool use_stdout;
		ofstream file_out;


		NamedInterval last_interval;
		string last_cols;
};

class SimpleOutput {
	public:
		SimpleOutput(const char *p) : path(p), use_stdout(path == "-") {
			if (!use_stdout)
				file_out.open(p, ofstream::out);
		}
		// write interval to bed
		// checks if we can extend the last interval
		// cols is appended after the interval name if not empty
		void operator() (const NamedInterval &i, const string &cols="") {
			if (writable()) {
				ostream *out = use_stdout ? &cout : &file_out;
				for (int b = i.start; b < i.end; ++b)
					if (cols.empty())
						*out << i.ref << '\t' << b << endl;
					else
						*out << i.ref << '\t' << b << '\t' << cols << endl;
			}
		}
	private:
		bool writable() const {
			return use_stdout || file_out.is_open();
		}
		
		const string path;
		const bool use_stdout;
		ofstream file_out;
};

#endif /*__INTERVAL_H__*/
