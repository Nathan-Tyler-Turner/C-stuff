#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <stdexcept>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <list>
#include <utility>
#include <stdexcept>
#include <deque>
#include <regex>
#include <iostream>
#include <mutex>
#include <thread>
template<typename K, typename V>
void print_map(std::unordered_map<K, V> const &m)
{
    for (auto const &pair: m) {
        std::cout << "{" << pair.first << ": " << pair.second << "}\n";
    }
}

std::mutex mtx;

#include "util.hpp"

void usage(char *progname);
std::vector<std::string> extractHitAccessions(const std::string& filename);
void parseXML(const std::string& filename, const std::string& outputFilename);
std::string get_common_name(const std::string& taxonomy_id);
void createDictionary(std::unordered_map<std::string, std::string>& taxid2name_dictionary, std::string names_filename);
std::string get_common_name(const std::unordered_map<std::string, std::string>& taxid2name_dictionary, const std::string& key);
std::vector<std::string> extractHitAccessionsMultiBLAST(const std::string& filename);
void parse_dmnd_output(const std::string& filename, std::vector<std::vector<std::string>>& result);
void parse_prot_acc2taxid_multithread(int threadNum, const std::string& filename, std::unordered_map<std::string, std::string>& dictionary, std::mutex& mergeMutex);
std::string processFileforBLASTdbcmd(const std::vector<std::string>& accession_nos, const std::string& sequence, const std::string& db_path, const std::string& outputFilename);
std::string executeCommand(const std::string& command);
std::vector<std::string> splitString(const std::string& str, char delimiter);

int main(int argc, char** argv) {

    std::unordered_map<Accession, TaxonId> acc2taxid;
	//std::cout << "Allocating memory..." << std::endl;
	//acc2taxid.reserve(4932678007);
    TaxTree nodes;
    std::unordered_map<TaxonId, TaxonName> node2name;
	//std::unordered_map<Accession, TaxonId> prot_acc2taxid_dict;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string acc2taxid_filename = "";
	std::string in1_filename = "";
	std::string out_filename;

	//These are only required if using BLAST command line tool
	std::string filename = "";
    std::string outputFilename = "parsed_blast_output.tsv";

	std::unordered_map<std::string, std::string> taxid2name_dictionary;
	//std::unordered_map<std::string, std::string> prot_acc2taxid_dict;
	
	
	bool full_path = false;
	bool specified_ranks = false;
    bool verbose = false;
	std::string multiblast = "false";
	std::string ranks_arg;
	std::string sequence;
	std::string db_path;

	// Variables related to multithreaded parse_prot_acc2tax_multithread function, not required if not generating dictionary
	std::vector<std::thread> threads;
    const int numThreads = 98;
    const std::string directoryPath = "processed/";

    std::unordered_map<std::string, std::string> mergedDictionary;
    std::mutex mergeMutex; // Mutex for thread-safe access to mergedDictionary

	//use the parseXML function to write a tsv file with the blast data
	//parseXML(filename, outputFilename);

	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hva:f:s:i:o:r:pd:")) != -1) {
		switch (c) {
			case 'h':
				usage(argv[0]);
				break;
			case 'v':
				verbose = true;
				break;
			//case 'a':
			//	acc2taxid_filename = optarg;
			//	break;
			case 'd':
				db_path = optarg;
				break;
			case 'p':
				full_path = true;
				break;
			case 'o':
				out_filename = optarg;
				break;
			case 'f':
				multiblast = optarg;
				break;
			case 's':
				sequence = optarg;
				break;
			case 'i':
				in1_filename = optarg;
				break;
			case 'r': {
				specified_ranks = true;
				ranks_arg = optarg;
				break;
			}
			default:
				usage(argv[0]);
		}
	}

	nodes_filename = db_path + "/FASTA/nodes.dmp";
	names_filename = db_path + "/FASTA/names.dmp";
	//Extract accession numbers from BLAST xml or multiBLAST table
	std::cout << "Generating taxid2commonname dictionary from..." << names_filename << std::endl;

	createDictionary(taxid2name_dictionary, names_filename);
	std::vector<std::string> accession_nos;
	if (multiblast == "true") {
		std::cout << "MultiBLAST table selected as input: " << in1_filename << std::endl;
		accession_nos = extractHitAccessionsMultiBLAST(in1_filename);
	} else if (multiblast == "false"){
		std::cout << "BLAST .xml file selected as input: " << in1_filename << std::endl;
		accession_nos = extractHitAccessions(in1_filename);
	} else if (multiblast == "diamond") {
		std::cout << "Diamond output selected as input: " << in1_filename << std::endl;
		
	}
	// Print the extracted accession numbers
    //for (const auto& accession : accession_nos) {
   //     std::cout << accession << std::endl;
    //}
	if(out_filename.length() == 0) { error("Error: Please specify the name of the output file, using the -o option."); usage(argv[0]); }
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
//	if(acc2taxid_filename.length() == 0) { error("Please specify the location of the nucl_gb.accession2taxid file, using the -a option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }

    std::ifstream nodes_file;
    nodes_file.open(nodes_filename);
    if(!nodes_file.is_open()) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
    if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
    parseNodesDmp(nodes, nodes_file);
    nodes_file.close();

    std::ifstream names_file;
    names_file.open(names_filename);
    //if(!names_file.is_open()) { error("Could not open file " + names_filename); usage(argv[0]); }
    if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
    parseNamesDmp(node2name, names_file);
    names_file.close();


    std::ifstream in1_file;
    in1_file.open(in1_filename);
    if(!in1_file.is_open()) {  error("Could not open file " + in1_filename); exit(EXIT_FAILURE); }

	std::ifstream acc2taxid_file;
	if (multiblast != "diamond") {
//		acc2taxid_file.open(acc2taxid_filename);
//		if(!acc2taxid_file.is_open()) { error("Could not open file " + acc2taxid_filename); exit(EXIT_FAILURE); }
//		if(verbose) std::cerr << "Reading accession to taxon id map from file " << acc2taxid_filename << std::endl;
		if (sequence == "nucleotide_threaded") {
			parse_accession2taxid(acc2taxid,acc2taxid_file);
//			acc2taxid_file.close();
		} else if (sequence == "protein_threaded") {
			for (int i = 0; i < numThreads; ++i) {
				std::string filename = directoryPath + "subfile_" + std::to_string(i) + ".txt";
				threads.emplace_back(parse_prot_acc2taxid_multithread, i, filename, std::ref(mergedDictionary), std::ref(mergeMutex));
			}
			// Join all the threads
			for (std::thread& t : threads) {
				t.join();
			}
			// Print the merged dictionary contents
			//std::cout << "Merged Dictionary Contents:" << std::endl;
			//for (const auto& pair : mergedDictionary) {
				//std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
				//continue;
			//}
		} else {
			//parse_accession2taxid(acc2taxid,acc2taxid_file);
			std::cout << "accession2taxid dictionary unnecesarry, skipping..." << std::endl;
		}
	}

    std::vector<std::vector<std::string>> result;

    //processAccessions(in1_filename, acc2taxid, node2name, result);
    //if(!out_file.is_open()) {  error("Could not open file " + out_filename + " for writing"); exit(EXIT_FAILURE); }
    if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;

    if(verbose) std::cerr << "Processing " << in1_filename <<"..." << "\n";
	
	//look at the acc2taxid unordered map
	//int count = 0;
	//for (const auto& pair : node2name) {
	//std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
	//count++;
	//if (count == 10)
	//	break;
    //}
	std::string taxid;
	if (multiblast != "diamond") {
		for (const auto& accession : accession_nos) {
			Accession acc;
			if (multiblast == "false") {
				//append a ".1" to the accession number to match format found in nucl_gb.accession2taxid
				acc = accession + ".1";			
			} else {
				acc = accession ;
				std::cout << "Using normal acc no: " << acc << std::endl;
			}
			if (sequence == "protein" || sequence == "nucleotide") {
				//taxid = mergedDictionary[acc];
				//std::string filename = in1_filename;
				//std::string sequence_1 = sequence;
				std::cout << "Using " << sequence << " pipeline" << std::endl;
    			std::string firstTaxid = processFileforBLASTdbcmd(accession_nos, sequence, db_path, out_filename);
    			std::cout << "First taxid: " << firstTaxid << std::endl;
				TaxonId id = std::stoull(firstTaxid);
				//auto it = mergedDictionary.find(acc);
				//taxid = it->second;
    			//std::cout << "Value of key '" << acc << "' is: " << id << std::endl;
				auto it_name = node2name.find(id);
				if (it_name == node2name.end()) {
					std::cerr << "Warning: Taxon ID " << id << " for accession " << acc << " is not contained in names.dmp file " << names_filename << ". Skipping.\n";
					continue;
				}
				std::string name = it_name->second;

				std::vector<std::string> entry;
				entry.push_back(acc);
				entry.push_back(std::to_string(id));
				entry.push_back(name);
				result.push_back(entry);
			} else {
				auto it_id = acc2taxid.find(acc);
				if(it_id == acc2taxid.end()) {
					std::cerr << "Warning: Accession " << acc << " is not found in "<< acc2taxid_filename << ". Skipping.\n";
					continue;
				}
				TaxonId id = it_id->second;
				//std::cout << "ID:" << id << std::endl;

				auto it_name = node2name.find(id);
				if (it_name == node2name.end()) {
					std::cerr << "Warning: Taxon ID " << id << " for accession " << acc << " is not contained in names.dmp file " << names_filename << ". Skipping.\n";
					continue;
				}
				std::string name = it_name->second;

				std::vector<std::string> entry;
				entry.push_back(acc);
				entry.push_back(std::to_string(id));
				entry.push_back(name);
				result.push_back(entry);
			}
		}
	} else if (multiblast == "diamond") {
		parse_dmnd_output(in1_filename, result);
	}

    in1_file.close();
    //out_file.close();
    //for (const auto& entry : result) {
    //    for (const auto& value : entry) {
    //        std::cout << value << "\t";
    //    }
    //    std::cout << "\n";          
    //}
	
//################################################taxonid2namestart###################################
	std::unordered_map<uint64_t,uint64_t> nodes_1;
	std::unordered_map<uint64_t, std::string> node2name_1;
	std::unordered_map<uint64_t, std::string> node2rank_1;
	std::unordered_map<uint64_t, std::string> node2path_1;

	std::list<std::string> ranks_list;
	std::set<std::string> ranks_set;

	bool filter_unclassified = false;

	// --------------------- START ------------------------------------------------------------------


	/* parse user-supplied rank list into list and set */
	if(ranks_arg.length() > 0) {
		size_t begin = 0;
		size_t pos = -1;
		std::string rankname;
		while((pos = ranks_arg.find(",",pos+1)) != std::string::npos) {
			rankname = ranks_arg.substr(begin,(pos - begin));
			if(rankname.length()==0 || rankname==",") { begin=pos+1; continue; }
			ranks_list.emplace_back(rankname);
			ranks_set.emplace(rankname);
			begin = pos+1;
		}
		rankname = ranks_arg.substr(begin);
		if(!(rankname.length()==0 || rankname==",")) {
			ranks_set.emplace(rankname);
			ranks_list.emplace_back(rankname);
		}
	}

	/* read nodes.dmp */
	//std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl;
	 }
	if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmpWithRank(nodes,node2rank_1,nodes_file);
	nodes_file.close();

	/* read names.dmp */
	//std::ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { error("Could not open file " + names_filename);}
	if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
	parseNamesDmp(node2name,names_file);
	names_file.close();

	std::ostream * out_stream;
	if(out_filename.length()>0) {
		if(verbose) std::cerr << "Output file: " << out_filename << std::endl;
		std::ofstream * output_file = new std::ofstream();
		output_file->open(out_filename);
		if(!output_file->is_open()) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
		out_stream = output_file;
	}
	else {
		out_stream = &std::cout;
	}

	std::string line;
	//write column headers
	*out_stream << "Hit_accession" << ',' << "taxid" << ',' << "name" << ',' << "taxonomy" << ',' << "common_name" << '\n';
	std::cout << "Looking up taxonomy info..." << std::endl;
	for (const auto& entry : result) {
		std::string acc = entry[0];
		TaxonId taxonid = std::stoul(entry[1]);
		//std::cout << "START" << std::to_string(taxonid) << "END" << std::endl;
		std::string common_name = get_common_name(taxid2name_dictionary, std::to_string(taxonid));
		//if (common_name != "") {
		//	std::cout << std::to_string(taxonid) << " common name: " << common_name << std::endl;
		//}
		//std::cout << acc << std::endl;
		line = entry[0] +"," + entry[1] + "," + entry[2];
		std::size_t dotPos = entry[0].find('.');
		std::string raw_acc_no = entry[0].substr(0, dotPos);
		//std::cout << "UPDATED" << node2path.at(taxonid) << std::endl;
		//TaxonId taxonid;
		try {
			taxonid = std::stoul(entry[1]);
			//std::cout << "taxid: " << taxonid << std::endl;
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad taxon id (out of range error) in line: " << line << std::endl;
			continue;
		}
		
		if(nodes.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
			continue;
		}
		if(node2name.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not found in file "<< names_filename << ".\n";
			continue;
		}

		if(full_path || specified_ranks) {
			if(node2path_1.count(taxonid)>0) { // look if path is already saved
				*out_stream << raw_acc_no +"," + entry[1] + "," + entry[2] + ',' << node2path_1.at(taxonid) << ',' << common_name << "\n";
				//std::cout << raw_acc_no +"," + entry[1] + "," + entry[2] + ',' << node2path_1.at(taxonid) << common_name << std::endl;
				continue;
			} //else {
				//std::cout << "Path not already saved" << std::endl;
				//std::cout << taxonid << std::endl;
				//std::cout << line << node2path.at(taxonid) << std::endl;
			//}
			std::deque<std::string> lineage; // for full_path
			std::map<std::string,std::string> curr_rank_values;
			if(specified_ranks) { //set the values for all specified ranks to NA, which will be overwritten by the actual values if they are found
				for(auto it : ranks_list) {
					curr_rank_values.emplace(it,"NA");
				}
			}
			//  go from leaf to root starting at taxonid and gather values for ranks
			uint64_t id = taxonid;
			while(nodes.count(id)>0 && id != nodes.at(id)) {
				std::string taxon_name;
				if(specified_ranks) {
					if(node2rank_1.count(id)==0 || node2rank_1.at(id)=="no rank") {  // no rank name
						id = nodes.at(id);
						continue;
					}
					std::string rank_name = node2rank_1.at(id);
					if(ranks_set.count(rank_name)==0) { // rank name is not in specified list of ranks
						id = nodes.at(id);
						continue;
					}
					taxon_name = getTaxonNameFromId(node2name, id, names_filename);
					curr_rank_values[rank_name] = taxon_name;
				}
				else { //full path
					taxon_name = getTaxonNameFromId(node2name, id, names_filename);
					lineage.emplace_front(taxon_name);
				}
				id = nodes.at(id);
			} // end while

			// assemble lineage into one string
			std::string lineage_text;
			if(specified_ranks) {
				for(auto it : ranks_list) {
					lineage_text += curr_rank_values[it];
					lineage_text += "; ";
				}
			}
			else { // full path
				for(auto  itl : lineage) {
					lineage_text += itl;
					lineage_text += "; ";
				}
			}
			// now lineage_text contains the final lineage for the taxon
			node2path_1.emplace(taxonid,lineage_text);
			*out_stream << raw_acc_no +"," + entry[1] + "," + entry[2] + ',' << lineage_text << ',' << common_name << "\n";
			//std::cout <<  raw_acc_no +"," + entry[1] + "," + entry[2] + "," << lineage_text << common_name << std::endl;
		}
		else {
			*out_stream << raw_acc_no +"," + entry[1] + "," + entry[2] + ',' << getTaxonNameFromId(node2name, taxonid, names_filename) << ',' << common_name << "\n";
			//std::cout << raw_acc_no +"," + entry[1] + "," + entry[2] + ',' << getTaxonNameFromId(node2name, taxonid, names_filename) << common_name << "\n";
		}
	}  // end while getline

	out_stream->flush();
	if(out_filename.length()>0) {
		((std::ofstream*)out_stream)->close();
		delete ((std::ofstream*)out_stream);
	}

	return 0;

}

void usage(char *progname) {
	fprintf(stderr, "Copyright 2018 Peter Menzel\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -a nucl_gb.accession2taxid -i blast.out -o blast2krona.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -d PATH       Path to blast databases (nt/nr/FASTA subfolders must be present) \n");
//	fprintf(stderr, "   -a FILENAME   Name of accession2taxid file\n");
	fprintf(stderr, "   -s nucleotide|protein  identify which database was searched (nt or nr)\n");
	fprintf(stderr, "   -f true or false, refers to whether or not the input is a MultiBLAST table or not \n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	exit(EXIT_FAILURE);
}

std::vector<std::string> extractHitAccessions(const std::string& filename) {
    std::vector<std::string> accession_nos;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Failed to open file: " << filename << std::endl;
        return accession_nos;
    }

    std::string line;
    std::string xmlContent;
    while (std::getline(file, line)) {
        xmlContent += line;
    }

    std::regex hitAccessionRegex("<Hit_accession>(.*?)</Hit_accession>");
    std::smatch match;

    auto pos = xmlContent.cbegin();
    auto end = xmlContent.cend();
    while (std::regex_search(pos, end, match, hitAccessionRegex)) {
        accession_nos.push_back(match[1].str());
        pos = match.suffix().first;
    }

    file.close();
    return accession_nos;
}

void parseXML(const std::string& filename, const std::string& outputFilename) {
    std::ifstream xmlFile(filename);
    if (!xmlFile) {
        std::cerr << "Error opening XML file: " << filename << std::endl;
        return;
    }

    std::ofstream outputFile(outputFilename);
    if (!outputFile) {
        std::cerr << "Error opening output file: " << outputFilename << std::endl;
        return;
    }

    std::string line;
    std::vector<std::pair<std::string, std::string>> hitData; // Use a pair to store tag name and value
    bool insideHitTag = false;
    bool headersWritten = false; // Flag to track if headers are already written

    while (std::getline(xmlFile, line)) {
        if (line.find("<Hit>") != std::string::npos) {
            insideHitTag = true;
            hitData.clear();
        }
        else if (line.find("</Hit>") != std::string::npos) {
            insideHitTag = false;

            // Write column headers as the first line in the TSV file if not already written
            if (!headersWritten) {
                for (const auto& data : hitData) {
                    outputFile << data.first << ",";
                }
                outputFile << std::endl;
                headersWritten = true;
            }

            // Write data values in subsequent lines
            for (const auto& data : hitData) {
                outputFile << data.second << ",";
            }
            outputFile << std::endl;
        }
        else if (insideHitTag && line.find("<") != std::string::npos && line.find(">") != std::string::npos) {
            std::string elementName = line.substr(line.find("<") + 1, line.find(">") - line.find("<") - 1);
            std::string elementValue = line.substr(line.find(">") + 1, line.rfind("<") - line.find(">") - 1);

            hitData.push_back(std::make_pair(elementName, elementValue));
        }
    }

    xmlFile.close();
    outputFile.close();
    std::cout << "Output file generated: " << outputFilename << std::endl;
}

std::string get_common_name(const std::string& taxonomy_id, std::string names_filename) {
    std::ifstream names_file(names_filename);
    std::string line;

    while (std::getline(names_file, line)) {
        std::istringstream iss(line);
        std::string field;
        std::vector<std::string> fields;

        while (std::getline(iss, field, '\t'))
            fields.push_back(field);

		//std::cout << fields[0] << "END" << std::endl;
		
        if (fields[0] == taxonomy_id) {
			std::cout << "Executed" << std::endl;
			std::cout << fields[1] << std::endl;
		}            
    }

    return "";
}

void createDictionary(std::unordered_map<std::string, std::string>& taxid2name_dictionary, std::string names_filename) {
	std::cout << "file" << names_filename;
    std::ifstream file(names_filename);  // Replace "data.txt" with your file name
    std::string line;
    bool is_common;

    if (!file) {
        std::cerr << "Failed to open the file." << std::endl;
        return;
    }
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string value;
        std::string value2;

        // Extract the value before the first pipe
        if (std::getline(iss, value, '|')) {
            value.erase(value.find_last_not_of(" \t") + 1);  // Trim trailing whitespace
            //std::cout << "Value 1: " << value << std::endl;
            // Extract the value between the first two pipes
            if (std::getline(iss, value2, '|')) {
                value2.erase(0, value2.find_first_not_of(" \t"));  // Trim leading whitespace
                value2.erase(value2.find_last_not_of(" \t") + 1);  // Trim trailing whitespace
                //std::cout << "Value 2: " << value2 << std::endl;
                //value2 = value;
            }            
        }

        // Extract the value between the last two pipes
        std::string value3;
        while (std::getline(iss, value3, '|')) {
            value3.erase(0, value3.find_first_not_of(" \t"));  // Trim leading whitespace
            value3.erase(value3.find_last_not_of(" \t") + 1);  // Trim trailing whitespace
            //value3 = value;
            if (value3 == "common name") {
                is_common = true;
                //std::cout << "Common name found for " << value << "\t" << value3 << std::endl;
                //std::cout << std::boolalpha << is_common << std::endl;
            } else {
                is_common = false;
                //std::cout << "Common name NOT found for " << value << std::endl;
            } // Store the last extracted value
        }
     

        //create kv pair using the array as a value
        //std::cout << std::noboolalpha << is_common << std::endl;
        if (is_common) {
            taxid2name_dictionary[value] = value2;
            //std::cout << "Common name found for " << value << " : " << value2 << std::endl;
        }
    }
    
    file.close();
}

std::string get_common_name(const std::unordered_map<std::string, std::string>& taxid2name_dictionary, const std::string& key) {
    auto it = taxid2name_dictionary.find(key);
    if (it != taxid2name_dictionary.end()) {
		//std::cout << "common name found!" << std::endl;
        return it->second;
    } else {
        return "No common name";  // Return empty string if key not found
    }
}

std::vector<std::string> extractHitAccessionsMultiBLAST(const std::string& filename) {
    std::vector<std::string> accession_nos;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cout << "Failed to open file: " << filename << std::endl;
        return accession_nos;
    }

    std::string line;
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> headers;
        std::string header;
        
        // Extract headers from the first line
        while (std::getline(iss, header, '\t')) {
            headers.push_back(header);
        }

        // Find the index of the "Accession (E-value)" header
        auto it = std::find(headers.begin(), headers.end(), "Accession (bit score)");
        if (it == headers.end()) {
            std::cout << "Column header not found: Accession (bit score)" << std::endl;
            return accession_nos;
        }
        int columnIndex = std::distance(headers.begin(), it);

        // Read the rest of the lines and extract values under the "Accession (E-value)" column
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::vector<std::string> values;
            std::string value;

            int colIndex = 0;
            while (std::getline(iss, value, '\t')) {
                if (colIndex == columnIndex) {
                    accession_nos.push_back(value);
                    break;  // No need to continue reading this line
                }
                colIndex++;
            }
        }
    }

    file.close();
    return accession_nos;
}

void parse_dmnd_output(const std::string& filename, std::vector<std::vector<std::string>>& result) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    // Skip the first three lines
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    // Read the remaining lines and extract the entries
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<std::string> entry;
        std::string value;
        while (std::getline(iss, value, '\t')) {
            entry.push_back(value);
        }

        // Add the entry to the result
        result.push_back(entry);
    }

    file.close();
}

void parse_prot_acc2taxid_multithread(int threadNum, const std::string& filename, std::unordered_map<std::string, std::string>& dictionary, std::mutex& mergeMutex) {
    std::ifstream file(filename);
    std::string line;
    int lineNumber = 1;
    std::unordered_map<std::string, std::string> threadDictionary; // Separate dictionary for each thread

    while (std::getline(file, line) && lineNumber <= 50000000) {
        std::istringstream iss(line);
        std::string key, value;
        if (std::getline(iss, key, '\t') || std::getline(iss, key, ' ')) {
            std::getline(iss, value);
            std::lock_guard<std::mutex> lock(mtx); // Lock to ensure thread-safe output
            lineNumber++;
            //std::cout << "Thread " << threadNum << ": Line " << lineNumber++ << ": " << key << " -> " << value << std::endl;
            //std::cout << "";
            threadDictionary[key] = value;
        }
    }
    //std::cout << threadNum << " has finished" << std::endl;

    file.close();

    // Merge the thread's dictionary into the main mergedDictionary
    std::lock_guard<std::mutex> lock(mergeMutex); // Lock to ensure thread-safe access to mergedDictionary
    dictionary.insert(threadDictionary.begin(), threadDictionary.end());
}

std::vector<std::string> splitString(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::istringstream iss(str);
    std::string token;
    while (std::getline(iss, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string executeCommand(const std::string& command) {
    std::string output;
    FILE* pipe = popen(command.c_str(), "r");
    if (pipe) {
        char buffer[128];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }
        pclose(pipe);
    } else {
        std::cerr << "Error executing command: " << command << std::endl;
    }
    return output;
}

std::string processFileforBLASTdbcmd(const std::vector<std::string>& accession_nos, const std::string& sequence, const std::string& db_path, const std::string& outputFilename) {
    std::ofstream outputFile(outputFilename);
    if (!outputFile.is_open()) {
        std::cerr << "Error opening output file: " << outputFilename << std::endl;
        return "";  // Return an empty string to indicate an error
    }
	
	for (const auto& value : accession_nos) {
        std::string command;
        if (sequence == "protein") {
            command = "blastdbcmd -db " + db_path + "/nr/nr/nr -entry \"" + value + "\" -outfmt %T";
        } else if (sequence == "nucleotide") {
            command = "blastdbcmd -db " + db_path + "/nt/nt/nt -entry \"" + value + "\" -outfmt %T";
        }

        std::string output = executeCommand(command);
        std::vector<std::string> taxids = splitString(output, '\n');
        if (!taxids.empty()) {
            // Perform any additional processing if needed
            outputFile << "Accession: " << value << " - TaxID: " << taxids[0] << std::endl;
        } else {
            outputFile << "No taxid found for accession: " << value << std::endl;
        }
        //std::cout << "-----------" << std::endl;
    }
	
	outputFile.close();
		// Return an empty string to indicate success
        return "";
    }