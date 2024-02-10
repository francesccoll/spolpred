/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   SpolPred.cpp

   Copyright (C) 2012 Francesc Coll (Francesc.Coll@lshtm.ac.uk)
   <http://www.gnu.org/licenses/lgpl-3.0.html>

   This file contains the whole software code.

   SpolPred is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License, version 3 (GPL-3.0)
   as published by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   SpolPred is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with SpolPred; if not, see http://www.opensource.org/licenses/gpl-3.0.html

   CODE STRUCTURE

    0. Function declarations
    1. Main function
	1.1 Default options and thresholds
	1.2 Saving user-defined options and thresholds
	1.3 Printing user-defined options and thresholds
	1.4 Defining spacers
	1.5 Scanning the fastq file
		1.5.1 Read processing and read-spacer matching
		1.5.2 Checking stop and detail options
		1.5.3 Printing detailed information on the screen
	1.6 Obtaining octal code from read-spacer matches
	1.7 Printing the resulting octal code and writing the output into a file	
    2. Function definitions
	2.1 seq_to_bits function
	2.2 read_to_bits function
	2.3 spacers_to_bits function
	2.4 comparison function
	2.5 octalcode function
	2.6 writeresults function

 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <cstring>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#define S_NUM 43

using namespace std;

/*
    0. FUNCTION DECLARATIONS
*/
unsigned int seq_to_bits(string seq, int from, int to);
void read_to_bits(string read, unsigned int **bit_subreads, int length);
void spacers_to_bits(string spacers_str[], unsigned int bit_spacers[][2]);
void comparison(unsigned int **bit_subreads, unsigned int bit_spacers[][2], int spacer_occur[], int spacer_occur1[], int length);
void octalcode(int spacer_01[], int spacer_octal[]);
void writeresults(int spacer_octal[], int spacer_occur[], int spacer_occur1[], const char *sample_name, const char *output_name);


/*
    1. MAIN FUNCTION
*/
int main(int argc, char *argv[]) {
    
    bool file_opened = true;
    /*
        1.1 Default options and thresholds
    */
    bool rev = false;
    bool stop_option = false;
    bool stop = false;
    int screen_thr = 50;
    bool details = true;
    int match_thr = 4;
    int length = 75;
    const char *output_name = "output.txt";
    
    int rm = argc%2;
    if ((argc < 2) || (argc > 14) || (rm != 0)) {
        cerr << "Incorrect number of arguments." << "\nUsage: SpolPred <reads.fastq> [-l <ReadLength>] [-b <d/r>] [-o <output.txt>] [-d <on/off>] [-s <on/off>] [-a <ScreeningThreshold>] [-m <MatchThreshold>]" << endl;
        exit(1);
    }

    /* Checking the file can be opened */
    ifstream file;
    file.open (argv[1]);
    if ( (file.rdstate() & ifstream::failbit ) != 0 ) {
    cerr << "Error opening fastq file\n";
    file_opened = false;
    }
    file.close();

    if (file_opened == true ) {
    /*
        1.2 Saving user-defined options and thresholds
    */
    for (int i = 0; i < argc; i++) {
        if (strncmp(argv[i], "-o", 2) == 0) {
            output_name = argv[i + 1];
        }
        if (strncmp(argv[i], "-l", 2) == 0) {
		length = atoi(argv[i + 1]);
		if (length == 0){
			cout << "\nERROR: Chosen length (-l) is not a valid value.";
			cout << "\nSpolpred stops\n";
		        exit(1);
		    } else if ((length < 35) || (length > 1000)) {
		        cout << "\nWARNING: chosen length (" << length << ") out of range. Using default length (75 bp)\n";
		        length = 75;
		    }
        }
        if (strncmp(argv[i], "-d", 2) == 0) {
            if ((strncmp(argv[i + 1], "off", 3) == 0) || (strncmp(argv[i + 1], "OFF", 3) == 0)) {
                details = false;
            } else if ((strncmp(argv[i + 1], "on", 2) == 0) || (strncmp(argv[i + 1], "ON", 2) == 0)) {
                details = true;
            } else {
		cout << "\nWARNING: Details Option (-d) should be either 'on' or 'off'. Using default value 'on'\n";
		details = true;
	    }
        }
        if (strncmp(argv[i], "-s", 2) == 0) {
            if ((strncmp(argv[i + 1], "off", 3) == 0) || (strncmp(argv[i + 1], "OFF", 3) == 0)){
                stop_option = false;
            } else if ((strncmp(argv[i + 1], "on", 2) == 0) || (strncmp(argv[i + 1], "ON", 2) == 0)) {
                stop_option = true;
            } else {
	        cout << "\nWARNING: Stop Option (-s) should be either 'on' or 'off'. Using default value 'off'\n";
	        stop_option = false;
 	      }
        }
        if (strncmp(argv[i], "-b", 2) == 0) {
            if ((strncmp(argv[i + 1], "d", 1) == 0) || (strncmp(argv[i + 1], "D", 1) == 0)){
                rev = false;
            } else if ((strncmp(argv[i + 1], "r", 1) == 0) || (strncmp(argv[i + 1], "R", 1) == 0)) {
                rev = true;
            } else {
	        cout << "\nWARNING: Mode option (-b) should be either 'd' (direct) or 'r' (reverse). Using default value 'd'\n";
	        rev = false;
 	      }
        }
        if (strncmp(argv[i], "-a", 2) == 0) {
            if (stop_option == false) {
                cout << "\nWARNING: Screening Stop option must be switch on first ([-s on])." << "\nUsing -a default value: " << screen_thr << "\n";
            } else {
                screen_thr = atoi(argv[i + 1]);
		if (screen_thr == 0) {
			cout << "\nERROR: Screening Threshold (-a) is not a valid value.";
			cout << "\nSpolpred stops\n";
		        exit(1);
		    } else if (screen_thr < 0) {
		        cout << "\nWARNING: Screening Threshold (-a): " << screen_thr << " should be a positive value.";
			cout << "Using default value (50)\n";
		        screen_thr = 50;
		    }
            }
        }
        if (strncmp(argv[i], "-m", 2) == 0) {
            match_thr = atoi(argv[i + 1]);
	    if (match_thr == 0) {
			cout << "\nERROR: Matching Threshold (-m) is not a valid value.";
			cout << "\nSpolpred stops\n";
		        exit(1);
	    } else if (match_thr < 0) {
		        cout << "\nWARNING: Matching Threshold (-m): " << match_thr << " should be a positive value. ";
			cout << "Using default value (4)\n";
		        match_thr = 4;
	    }
        }
    }

    /*
        1.3 Printing user-defined options and thresholds
    */
    cout << "\nSample name:" << "\t" << argv[1];
    cout << "\nOutput name:" << "\t" << output_name;
    cout << "\nRead Length:" << "\t" << length;
    if(rev == false) {
    cout << "\nMode:" << "\t" << "direct reads";
    }
    if(rev == true) {
    cout << "\nMode:" << "\t" << "reverse reads";
    }
    cout << "\nNumber of read-spacer matches below which absence is assigned:" << "\t" << match_thr;
    if (stop_option == true) {
        cout << "\nStop Screening option: ON\nThreshold = " << screen_thr << " spacer-read matches on average.";
    } else {
        cout << "\nStop Reads Screening: OFF";
    }
    cout << "\n\n";

    /*
        1.4 Defining spacers
    */
    std::string spacers_str[S_NUM] = { "" };

    if(rev == false){
	    spacers_str[0] = "ATAGAGGGTCGCCGGCTCTGGATCA";
	    spacers_str[1] = "CCTCATGCTTGGGCGACAGCTTTTG";
	    spacers_str[2] = "CCGTGCTTCCAGTGATCGCCTTCTA";
	    spacers_str[3] = "ACGTCATACGCCGACCAATCATCAG";
	    spacers_str[4] = "TTTTCTGACCACTTGTGCGGGATTA";
	    spacers_str[5] = "CGTCGTCATTTCCGGCTTCAATTTC";
	    spacers_str[6] = "GAGGAGAGCGAGTACTCGGGGCTGC";
	    spacers_str[7] = "CGTGAAACCGCCCCCAGCCTCGCCG";
	    spacers_str[8] = "ACTCGGAATCCCATGTGCTGACAGC";
	    spacers_str[9] = "TCGACACCCGCTCTAGTTGACTTCC";
	    spacers_str[10] = "GTGAGCAACGGCGGCGGCAACCTGG";
	    spacers_str[11] = "ATATCTGCTGCCCGCCCGGGGAGAT";
	    spacers_str[12] = "GACCATCATTGCCATTCCCTCTCCC";
	    spacers_str[13] = "GGTGTGATGCGGATGGTCGGCTCGG";
	    spacers_str[14] = "CTTGAATAACGCGCAGTGAATTTCG";
	    spacers_str[15] = "CGAGTTCCCGTCAGCGTCGTAAATC";
	    spacers_str[16] = "GCGCCGGCCCGCGCGGATGACTCCG";
	    spacers_str[17] = "CATGGACCCGGGCGAGCTGCAGATG";
	    spacers_str[18] = "TAACTGGCTTGGCGCTGATCCTGGT";
	    spacers_str[19] = "TTGACCTCGCCAGGAGAGAAGATCA";
	    spacers_str[20] = "TCGATGTCGATGTCCCAATCGTCGA";
	    spacers_str[21] = "ACCGCAGACGGCACGATTGAGACAA";
	    spacers_str[22] = "AGCATCGCTGATGCGGTCCAGCTCG";
	    spacers_str[23] = "CCGCCTGCTGGGTGAGACGTGCTCG";
	    spacers_str[24] = "GATCAGCGACCACCGCACCCTGTCA";
	    spacers_str[25] = "CTTCAGCACCACCATCATCCGGCGC";
	    spacers_str[26] = "GGATTCGTGATCTCTTCCCGCGGAT";
	    spacers_str[27] = "TGCCCCGGCGTTTAGCGATCACAAC";
	    spacers_str[28] = "AAATACAGGCTCCACGACACGACCA";
	    spacers_str[29] = "GGTTGCCCCGCGCCCTTTTCCAGCC";
	    spacers_str[30] = "TCAGACAGGTTCGCGTCGATCAAGT";
	    spacers_str[31] = "GACCAAATAGGTATCGGCGTGTTCA";
	    spacers_str[32] = "GACATGACGGCGGTGCCGCACTTGA";
	    spacers_str[33] = "AAGTCACCTCGCCCACACCGTCGAA";
	    spacers_str[34] = "TCCGTACGCTCGAAACGCTTCCAAC";
	    spacers_str[35] = "CGAAATCCAGCACCACATCCGCAGC";
	    spacers_str[36] = "CGCGAACTCGTCCACAGTCCCCCTT";
	    spacers_str[37] = "CGTGGATGGCGGATGCGTTGTGCGC";
	    spacers_str[38] = "GACGATGGCCAGTAAATCGGCGTGG";
	    spacers_str[39] = "CGCCATCTGTGCCTCATACAGGTCC";
	    spacers_str[40] = "GGAGCTTTCCGGCTTCTATCAGGTA";
	    spacers_str[41] = "ATGGTGGGACATGGACGAGCGCGAC";
	    spacers_str[42] = "CGCAGAATCGCACCGGGTGCGGGAG";
    }
    if(rev == true) {
	spacers_str[0] = "TGATCCAGAGCCGGCGACCCTCTAT";
	spacers_str[1] = "CAAAAGCTGTCGCCCAAGCATGAGG";
	spacers_str[2] = "TAGAAGGCGATCACTGGAAGCACGG";
	spacers_str[3] = "CTGATGATTGGTCGGCGTATGACGT";
	spacers_str[4] = "TAATCCCGCACAAGTGGTCAGAAAA";
	spacers_str[5] = "GAAATTGAAGCCGGAAATGACGACG";
	spacers_str[6] = "GCAGCCCCGAGTACTCGCTCTCCTC";
	spacers_str[7] = "CGGCGAGGCTGGGGGCGGTTTCACG";
	spacers_str[8] = "GCTGTCAGCACATGGGATTCCGAGT";
	spacers_str[9] = "GGAAGTCAACTAGAGCGGGTGTCGA";
	spacers_str[10] = "CCAGGTTGCCGCCGCCGTTGCTCAC";
	spacers_str[11] = "ATCTCCCCGGGCGGGCAGCAGATAT";
	spacers_str[12] = "GGGAGAGGGAATGGCAATGATGGTC";
	spacers_str[13] = "CCGAGCCGACCATCCGCATCACACC";
	spacers_str[14] = "CGAAATTCACTGCGCGTTATTCAAG";
	spacers_str[15] = "GATTTACGACGCTGACGGGAACTCG";
	spacers_str[16] = "CGGAGTCATCCGCGCGGGCCGGCGC";
	spacers_str[17] = "CATCTGCAGCTCGCCCGGGTCCATG";
	spacers_str[18] = "ACCAGGATCAGCGCCAAGCCAGTTA";
	spacers_str[19] = "TGATCTTCTCTCCTGGCGAGGTCAA";
	spacers_str[20] = "TCGACGATTGGGACATCGACATCGA";
	spacers_str[21] = "TTGTCTCAATCGTGCCGTCTGCGGT";
	spacers_str[22] = "CGAGCTGGACCGCATCAGCGATGCT";
	spacers_str[23] = "CGAGCACGTCTCACCCAGCAGGCGG";
	spacers_str[24] = "TGACAGGGTGCGGTGGTCGCTGATC";
	spacers_str[25] = "GCGCCGGATGATGGTGGTGCTGAAG";
	spacers_str[26] = "ATCCGCGGGAAGAGATCACGAATCC";
	spacers_str[27] = "GTTGTGATCGCTAAACGCCGGGGCA";
	spacers_str[28] = "TGGTCGTGTCGTGGAGCCTGTATTT";
	spacers_str[29] = "GGCTGGAAAAGGGCGCGGGGCAACC";
	spacers_str[30] = "ACTTGATCGACGCGAACCTGTCTGA";
	spacers_str[31] = "TGAACACGCCGATACCTATTTGGTC";
	spacers_str[32] = "TCAAGTGCGGCACCGCCGTCATGTC";
	spacers_str[33] = "TTCGACGGTGTGGGCGAGGTGACTT";
	spacers_str[34] = "GTTGGAAGCGTTTCGAGCGTACGGA";
	spacers_str[35] = "GCTGCGGATGTGGTGCTGGATTTCG";
	spacers_str[36] = "AAGGGGGACTGTGGACGAGTTCGCG";
	spacers_str[37] = "GCGCACAACGCATCCGCCATCCACG";
	spacers_str[38] = "CCACGCCGATTTACTGGCCATCGTC";
	spacers_str[39] = "GGACCTGTATGAGGCACAGATGGCG";
	spacers_str[40] = "TACCTGATAGAAGCCGGAAAGCTCC";
	spacers_str[41] = "GTCGCGCTCGTCCATGTCCCACCAT";
	spacers_str[42] = "CTCCCGCACCCGGTGCGATTCTGCG";
    }

    unsigned int bit_spacers[S_NUM][2];
    int spacer_occur[S_NUM] = { 0 }; // Array to track the number of read-spacer matches
    int spacer_occur1[S_NUM] = { 0 }; //Allowing one mismatch

    spacers_to_bits(spacers_str, bit_spacers); /*String spacers to bit spacers, i.e. filling unsigned int bit_spacer[S_NUM][2]*/

    /*
        1.5 Scanning the fastq file
    */
    ifstream finReads;
    finReads.open(argv[1]);
    string line;
    double proc = 0;
    double step = 100000;
    int mean = 0;

    while (!finReads.eof() && (stop == false)) {
        getline(finReads, line);
        getline(finReads, line);
        string read = line;
        proc++;
        /*
            1.5.1 Read processing and read-spacer matching
        */
        if (read.length() >= length) {
            unsigned int **bit_subreads = 0;
            bit_subreads = new unsigned int *[length - 25];
            for (int i = 0; i < (length - 25); i++) {
                bit_subreads[i] = new unsigned int[2];
            }

            read_to_bits(read, bit_subreads, length); // unsigned int bit_subreads[length-25][2];
            comparison(bit_subreads, bit_spacers, spacer_occur, spacer_occur1, length);

            for (int i = 0; i < (length - 25); i++) {
                delete [] bit_subreads[i];
            }
            delete [] bit_subreads;
        } /*else {
           //stop = true;
	   cout << "Error: Chosen length longer than the real one: ";
	   cout << read.length();
	   cout << "\n";
	}*/
        /*
            1.5.2 Checking stop and detail options
        */
        if (fmod(proc, step) == 0) {
            if (stop_option == true) {
                mean = 0;   // Average number of matches per 'present' spacer. Value used to stop read screening.
                int total = 0;
                int count = 0;
                for (int i = 0; i < 43; i++) {
                    if (spacer_occur1[i] > match_thr) { // match_thr: matching threshold
                        total = total + spacer_occur1[i];
                        count++;
                    }
                    mean = (double)total / count;
                    if (mean > (double)screen_thr) { // Breaking the while loop and therefore read screening
                        stop = true;
                    }
                }
            }

            /*
                1.5.3 Printing detailed information on the screen
            */
            if (details == true) {
                cout << "\n" << proc << " sequences processed";
                cout << "\nExact match";
                for (int j = 0; j < S_NUM; j++) {
                    cout << " " << spacer_occur[j];
                }
                cout << "\nOne mismatch";
                for (int j = 0; j < S_NUM; j++) {
                    cout << " " << spacer_occur1[j];
                }
                if (mean > 0) {
                    cout << "\nAverage matches: " << mean;
                }
            }
        }

        /* Discarding third and forth lines */
        getline(finReads, line);
        getline(finReads, line);
    }
    finReads.close();

    /*
        1.6 Obtaining octal code from read-spacer matches
    */
    int spacer_01[S_NUM] = { 0 };
    for (int i = 0; i < S_NUM; i++) {
        if (spacer_occur1[i] > match_thr) {
            spacer_01[i] = 1;
        }
    }
    int spacer_octal[15];
    octalcode(spacer_01, spacer_octal);

    /*
        1.7 Printing the resulting octal code and writing the output into a file
    */
    cout << "\nOctal code: ";
    for (int j = 0; j < 15; j++) {
        cout << " " << spacer_octal[j];
    }
    cout << "\n\n";
    writeresults(spacer_octal, spacer_occur, spacer_occur1, argv[1], output_name);

    } // end of if (file_opened == true )
    return 0;
}


/*
    2. FUNCTION DEFINTIONS
*/

/*
    2.1 seq_to_bits function
    It returns a 32-bit unsigned int representing an up to 16-bp subsequence from a given DNA sequence.
    From and to are positions at the DNA sequence, not array indices, therefore from=1 points at the first nucleotide.
*/
unsigned int seq_to_bits(string seq, int from, int to) {
    string sequence = seq;
    unsigned int liValue = 0;
    for (int i = from - 1; i < to; i++) {
        switch (sequence.at(i)) {
        case 'C':
            liValue = liValue << 2;
            liValue = liValue + 1;
            break;
        case 'G':
            liValue = liValue << 2;
            liValue = liValue + 2;
            break;
        case 'A':
            liValue = liValue << 2;
            liValue = liValue + 0;
            break;
        case 'T':
            liValue = liValue << 2;
            liValue = liValue + 3;
            break;
        }
    }
    return liValue;
}

/*
    2.2 read_to_bits function
    The function splits an input read (string) into 16-bp overlapping subreads, which are 1-bp shifted
    with respect to the previous one. Those sequences, in the form of 32-bit unsigned integers,
    are stored into the array passed in as the second argument.
*/
void read_to_bits(string read, unsigned int **bit_subreads, int length) {
    string r = read;
    for (int i = 1; i <= (length - 25); i++) {
        unsigned int libsr16 = 0;  // libsr = local integer bit subread
        unsigned int libsr9 = 0;
        libsr16 = seq_to_bits(r, i, i + 15);
        libsr9 = seq_to_bits(r, i + 16, i + 24);
        bit_subreads[i - 1][0] = libsr16;
        bit_subreads[i - 1][1] = libsr9;
    }
    return;
}

/*
    2.3 spacers_to_bits function
    This function splits each spacer (string) into two unsigned integers.
*/
void spacers_to_bits(string spacers_str[], unsigned int bit_spacers[][2]) {
    for (int i = 0; i < S_NUM; i++) {
        string s = spacers_str[i];
        unsigned int libsr16 = 0;
        unsigned int libsr9 = 0;
        libsr16 = seq_to_bits(s, 1, 16);
        libsr9 = seq_to_bits(s, 17, 25);
        bit_spacers[i][0] = libsr16;
        bit_spacers[i][1] = libsr9;
    }
    return;
}


/*
    2.4 comparison function
    This function tracks the number of read-spacer matches (allowing up to 1 SNP)
*/
void comparison(unsigned int **bit_subreads, unsigned int bit_spacers[][2], int spacer_occur[], int spacer_occur1[], int length) {
    unsigned int a, b, k, error, maxError;
    maxError = 2;   // Maximum allowed error 2 bits --> one SNP

    for (int i = 0; i < S_NUM; i++) {
        for (int j = 0; j < (length - 25); j++) {
            error = 0;
            a = bit_spacers[i][0] ^ bit_subreads[j][0]; // Comparing the first 16 positions
            for (k = 0; k < 16; k++) {
                if ((a & 3) > 0) {
                    error++;
                    if (error > maxError) {
                        break;
                    }
                }
                a = a >> 2;
            }
            if (error < maxError) { // If no more than 1 SNP occurs
                b = bit_spacers[i][1] ^ bit_subreads[j][1]; // Comparing the remaining 9 positions
                for (k = 0; k < 16; k++) {
                    if ((b & 3) > 0) {
                        error++;
                        if (error > maxError) {
                            break;
                        }
                    }
                    b = b >> 2;
                }
                if (error < maxError) {
                    spacer_occur1[i]++; // Only one SNP
                    if (error == 0) { // If no SNPs are found
                        spacer_occur[i]++;
                    }
                }
            }
        }
    }
    return;
}

/*
    2.5 octalcode function
    This function converts absence-presence code into the octal code
*/
void octalcode(int spacer_01[], int spacer_octal[]) {
    for (int i = 0; i < 40 ; i = i + 3) {
        int j = i;
        int k;
        if (spacer_01[j] == 0) {
            if (spacer_01[j + 1] == 0) {
                if (spacer_01[j + 2] == 0) {
                    k = j / 3;
                    spacer_octal[k] = 0;
                    // 000 = 0
                } else {
                    k = j / 3;
                    spacer_octal[k] = 1;
                    // 001 = 1
                }
            } else if (spacer_01[j + 1] == 1) {
                if (spacer_01[j + 2] == 0) {
                    k = j / 3;
                    spacer_octal[k] = 2;
                    // 010 = 2
                } else {
                    k = j / 3;
                    spacer_octal[k] = 3;
                    // 011 = 3
                }
            }
        } else if (spacer_01[j] == 1) {
            if (spacer_01[j + 1] == 0) {
                if (spacer_01[j + 2] == 0) {
                    k = j / 3;
                    spacer_octal[k] = 4;
                    // 100 = 4
                } else {
                    k = j / 3;
                    spacer_octal[k] = 5;
                    // 101 = 5
                }
            } else if (spacer_01[j + 1] == 1) {
                if (spacer_01[j + 2] == 0) {
                    k = j / 3;
                    spacer_octal[k] = 6;
                    // 110 = 6
                } else {
                    k = j / 3;
                    spacer_octal[k] = 7;
                    // 111 = 7
                }
            }
        }
    }

    if (spacer_01[42] == 0) {
        spacer_octal[14] = 0;
    } else {
        spacer_octal[14] = 1;
    }
    return;
}

/*
    2.6 writeresults function
    This function writes the results (octal code and number of read-spacer matches)
*/
void writeresults(int spacer_octal[], int spacer_occur[], int spacer_occur1[], const char *sample_name, const char *output_name) {
    string o_code = "";
    for (int j = 0; j < 15; j++) {
        std::ostringstream sin;
        sin << spacer_octal[j];
        std::string val = sin.str();
        o_code.append(val);
    }

    string spc_occur = "";
    for (int i = 0; i < S_NUM; i++) {
        std::ostringstream sin; // int to string
        sin << spacer_occur[i];
        std::string val = sin.str();
        spc_occur.append(val);
        spc_occur.append(" ");
    }
    spc_occur = spc_occur.substr(0, spc_occur.size() - 1); // Removing last space

    string spc_occur1 = "";
    for (int i = 0; i < S_NUM; i++) {
        std::ostringstream sin1; // int to string
        sin1 << spacer_occur1[i];
        std::string val1 = sin1.str();
        spc_occur1.append(val1);
        spc_occur1.append(" ");
    }
    spc_occur1 = spc_occur1.substr(0, spc_occur1.size() - 1);

    ofstream rfile; // Writing the results into a text file
    rfile.open(output_name, ios::out | ios::app);

    string rline = "";
    rline.append(sample_name);
    rline.append("\t");
    rline.append(o_code);
    rline.append("\t");
    rline.append(spc_occur);
    rline.append("\t");
    rline.append(spc_occur1);
    rline.append("\n");
    rfile << rline; // Writing the result line into the results file
    rfile.close();

    return;
}




