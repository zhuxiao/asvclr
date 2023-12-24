#include "consistency.h"

int max(int a, int b, int c) {
    return max(a, max(b, c));
}
void needleman_wunschOp(const string& seq1, const string& seq2, int32_t match_score, int32_t mismatch_score, int32_t gap_penalty, string& seq1_new, string& seq2_new) {
    int64_t rowsNum = seq1.size() + 1, colsNum = seq2.size() + 1, arrSize;
    int32_t i, j, * score_matrix;
    int32_t score, diag_score, up_score;

    // Initialize the score matrix
    arrSize = rowsNum * colsNum;
    score_matrix = (int32_t*)calloc(arrSize, sizeof(int32_t));
    if (score_matrix == NULL) {
        cerr << "line=" << __LINE__ << ", rowsNum=" << rowsNum << ", colsNum=" << colsNum << ", cannot allocate memory, error!" << endl;
        exit(1);
    }

    // Initialize the first row and first column
    for (i = 1; i < rowsNum; i++) score_matrix[i * colsNum] = score_matrix[(i - 1) * colsNum] + gap_penalty;
    for (j = 1; j < colsNum; j++) score_matrix[j] = score_matrix[j - 1] + gap_penalty;

    // Fill in the rest of the score matrix
    int32_t match, delete_gap, insert_gap;
    for (i = 1; i < rowsNum; i++) {
        for (j = 1; j < colsNum; j++) {
            match = score_matrix[(i - 1) * colsNum + j - 1] + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score);
            delete_gap = score_matrix[(i - 1) * colsNum + j] + gap_penalty;
            insert_gap = score_matrix[i * colsNum + j - 1] + gap_penalty;
            score_matrix[i * colsNum + j] = max(match, delete_gap, insert_gap);
        }
    }

    // Traceback to find the alignment
    i = rowsNum - 1, j = colsNum - 1;
    while (i > 0 && j > 0) {
        score = score_matrix[i * colsNum + j];
        diag_score = score_matrix[(i - 1) * colsNum + j - 1];
        //left_score = score_matrix[i*colsNum+j-1];
        up_score = score_matrix[(i - 1) * colsNum + j];

        if (score == diag_score + (seq1[i - 1] == seq2[j - 1] ? match_score : mismatch_score)) {
            seq1_new = seq1[i - 1] + seq1_new;
            seq2_new = seq2[j - 1] + seq2_new;
            i--;
            j--;
        }
        else if (score == up_score + gap_penalty) {
            seq1_new = seq1[i - 1] + seq1_new;
            seq2_new = "-" + seq2_new;
            i--;
        }
        else {
            seq1_new = "-" + seq1_new;
            seq2_new = seq2[j - 1] + seq2_new;
            j--;
        }
    }

    while (i > 0) {
        seq1_new = seq1[i - 1] + seq1_new;
        seq2_new = "-" + seq2_new;
        i--;
    }

    while (j > 0) {
        seq1_new = "-" + seq1_new;
        seq2_new = seq2[j - 1] + seq2_new;
        j--;
    }

    free(score_matrix);
}

double calculate_consistency(const string& seq1, const string& seq2) {

    int matching_chars = 0;
    int seq_length = max(seq1.length(), seq2.length());

    for (int32_t i = 0; i < seq_length; ++i) {
        if (i < (int32_t)seq1.length() && i < (int32_t)seq2.length() && seq1[i] == seq2[i]) {
            matching_chars++;
        }
    }

    double relief = 0;

    int consecutive_gap_count = 0;
    for (int i = 0; i <= seq_length; ++i) {
        if (i < (int)seq1.length() && i < (int)seq2.length() && (seq1[i] == '-' || seq2[i] == '-')) {
            consecutive_gap_count++;
        }
        else {
            if (consecutive_gap_count >= 3) {   //>
                relief += (consecutive_gap_count * RELIEF_FACTOR);
            }
            consecutive_gap_count = 0;
        }
    }

    double consistency = static_cast<double>(matching_chars + relief) / seq_length;
    return consistency;
}

size_t customHashFunction(const string& kmer, size_t kmerSize) {
    size_t hashValue = 0;
    size_t base = 4; 

    for (size_t i = 0; i < kmer.length(); ++i) {
        char baseChar = kmer[i];

        size_t baseValue;
        switch (baseChar) {
        case 'A':
            baseValue = 0;
            break;
        case 'C':
            baseValue = 1;
            break;
        case 'G':
            baseValue = 2;
            break;
        case 'T':
            baseValue = 3;
            break;
        default:
            //cerr << "Error: Invalid character in kmer." << "   kmer " << kmer << endl;
            hashValue = INT_MAX;
            return hashValue;
        }
        hashValue += baseValue * pow(base, kmerSize - i - 1);
    }
    return hashValue;
}

vector<Minimizer> findMinimizers(const string& sequence, size_t windowSize, size_t kmerSize) {
    vector<Minimizer> minimizers;

    for (size_t i = 0; i <= sequence.length() - windowSize - kmerSize + 1; ++i) {
        string window = sequence.substr(i, windowSize + kmerSize - 1);
        size_t minHash = INT_MAX;
        string minKmer = "";
        size_t pos;

        for (size_t j = 0; j <= window.size() - kmerSize; ++j) {
            string kmer = window.substr(j, kmerSize);
            size_t hashValue = customHashFunction(kmer, kmerSize);

            if (hashValue < minHash) {
                minHash = hashValue;
                minKmer = kmer;
                pos = i + j;

            }
        }
        Minimizer minimizer(pos, minHash, minKmer, false);
        if (!minimizers.empty()) {
            const Minimizer& lastMinimizer = minimizers.back();
            if (minHash != lastMinimizer.hashValue or pos != lastMinimizer.position)
                minimizers.push_back(minimizer);
        }
        else {
            minimizers.push_back(minimizer);
        }
    }

    return minimizers;
}

void retainCommonHashValuesMinimizers(vector<Minimizer>& containerA, vector<Minimizer>& containerB) {

    unordered_multiset<size_t> setB;
    for (const auto& element : containerB) {
        setB.insert(element.hashValue);
    }

    containerA.erase(remove_if(containerA.begin(), containerA.end(),
        [&setB](const Minimizer& element) {
            return setB.find(element.hashValue) == setB.end();
        }),
        containerA.end());

    unordered_multiset<size_t> setA;
    for (const auto& element : containerA) {
        setA.insert(element.hashValue);
    }

    containerB.erase(remove_if(containerB.begin(), containerB.end(),
        [&setA](const Minimizer& element) {
            return setA.find(element.hashValue) == setA.end();
        }),
        containerB.end());
}

void findSimilarityPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1) {

    int record_pos = -1;
    if (minimizers.empty() or minimizers1.empty())  return;
    if (minimizers.size() <= minimizers1.size()) {
        for (int i = 0; i < (int)minimizers.size(); i++) {
            int start = max(0, i - 5);
            int end = min(i + 5, (int)minimizers1.size() - 1);

            if (record_pos != -1 and i != 0) {
                if (minimizers[record_pos].matchflag and abs((int)minimizers[i].position - (int)minimizers[record_pos].position) <= 15) continue;
            }
            //another
            int MinPos = INT_MAX;
            int j_subscript = -1;

            for (int j = start; j <= end; j++) {


                if (minimizers1[j].matchflag == true)   continue;
                if ((minimizers[i].hashValue == minimizers1[j].hashValue) and (!minimizers[i].matchflag and !minimizers1[j].matchflag)) {
                    if (abs((int)minimizers1[i].position - (int)minimizers[j].position) < MinPos)
                    {
                        MinPos = abs((int)minimizers1[i].position - (int)minimizers[j].position);
                        j_subscript = j;
                    }
                }
            }
            //another
            if (j_subscript != -1 and MinPos < KMERSIZE) {
                minimizers[i].matchflag = minimizers1[j_subscript].matchflag = true;
                record_pos = i;
            }
        }
    }
    else {
        for (int i = 0; i < (int)minimizers1.size(); i++) {

            int start = max(0, i - 5);
            int end = min(i + 5, (int)minimizers.size() - 1);

            if (record_pos != -1 and i != 0) {
                if (minimizers1[record_pos].matchflag and abs((int)minimizers1[i].position - (int)minimizers1[record_pos].position) <= 15) continue;
            }

            //another
            int MinPos = INT_MAX;
            int j_subscript = -1;
            for (int j = start; j <= end; j++) {

                if (minimizers[j].matchflag == true)   continue;

                    //another
                if ((minimizers1[i].hashValue == minimizers[j].hashValue) and (!minimizers1[i].matchflag and !minimizers[j].matchflag)) {
                    if (abs((int)minimizers1[i].position - (int)minimizers[j].position) < MinPos)
                    {
                        MinPos = abs((int)minimizers1[i].position - (int)minimizers[j].position);
                        j_subscript = j;
                    }
                    
                }
            }
            //another
            if (j_subscript != -1 and MinPos < KMERSIZE) {
                minimizers1[i].matchflag = minimizers[j_subscript].matchflag = true;
                record_pos = i;
            }
        }
    }
}

void FilterDissimilaPosMinimizers(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1) {

    minimizers.erase(remove_if(minimizers.begin(), minimizers.end(),
        [](const Minimizer& m) { return !m.matchflag; }),
        minimizers.end());

    minimizers1.erase(remove_if(minimizers1.begin(), minimizers1.end(),
        [](const Minimizer& m) { return !m.matchflag; }),
        minimizers1.end());
}

void LongSequenceSplitAlignment(vector<Minimizer>& minimizers, vector<Minimizer>& minimizers1, int k, string& AlignSeq, string& AlignSeq1, string sequence, string sequence1) {

    int Exact_len, Exact_len1;
    string seq = "", seq1 = "", temporarily_seq = "", temporarily_seq1 = "";
    if (minimizers.size() != minimizers1.size() or (minimizers.empty() or minimizers1.empty())) return;
    int Length = minimizers.size(), sequece_len = sequence.size(), seqence_len1 = sequence1.size();
    for (int i = 0; i <= Length; i++) {
        if (i == 0) {
            seq = sequence.substr(0, minimizers[i].position);
            seq1 = sequence.substr(0, minimizers1[i].position);
        }
        else if (i == Length) {
            Exact_len = sequece_len - (minimizers[i - 1].position + k);
            Exact_len1 = seqence_len1 - (minimizers1[i - 1].position + k);
            seq = sequence.substr(minimizers[i - 1].position + k, Exact_len); 

            seq1 = sequence1.substr(minimizers1[i - 1].position + k, Exact_len1);
        }
        else {
            Exact_len = minimizers[i].position - (minimizers[i - 1].position + k);
            Exact_len1 = minimizers1[i].position - (minimizers1[i - 1].position + k);
            seq = sequence.substr(minimizers[i - 1].position + k, Exact_len);
            seq1 = sequence1.substr(minimizers1[i - 1].position + k, Exact_len1);

        }

        needleman_wunschOp(seq, seq1, MATCHSCORE, MISMATCHSCORE, GAPPENALTY, temporarily_seq, temporarily_seq1);


        if (i != Length) {
            temporarily_seq = temporarily_seq + minimizers[i].sequence;
            temporarily_seq1 = temporarily_seq1 + minimizers1[i].sequence;
            AlignSeq += temporarily_seq;
            AlignSeq1 += temporarily_seq1;

        }
        else if (i == Length) {
            AlignSeq += temporarily_seq;
            AlignSeq1 += temporarily_seq1;

        }

        seq = "", seq = "", temporarily_seq = "", temporarily_seq1 = "";
    }

}

double computeVarseqConsistency(const string& seq, const string& seq1) {
    string  seq1_new, seq2_new, AlignSeq, AlignSeq1;
    double result;
    if (seq.size() < MAX_SEQ_LEN or seq1.size() < MAX_SEQ_LEN) {
        needleman_wunschOp(seq, seq1, MATCHSCORE, MISMATCHSCORE, GAPPENALTY, seq1_new, seq2_new);
        result = calculate_consistency(seq1_new, seq2_new);
    }
    else {
        vector<Minimizer> minimizers = findMinimizers(seq, WINDOWSIZE, KMERSIZE);
        vector<Minimizer> minimizers1 = findMinimizers(seq1, WINDOWSIZE, KMERSIZE);
        retainCommonHashValuesMinimizers(minimizers, minimizers1);
        findSimilarityPosMinimizers(minimizers, minimizers1);
        FilterDissimilaPosMinimizers(minimizers, minimizers1);
        LongSequenceSplitAlignment(minimizers, minimizers1, KMERSIZE, AlignSeq, AlignSeq1, seq, seq1);
        if (minimizers.empty() or minimizers1.empty()) {
            needleman_wunschOp(seq, seq1, MATCHSCORE, MISMATCHSCORE, GAPPENALTY, seq1_new, seq2_new);
            result = calculate_consistency(seq1_new, seq2_new);
        }
        else
            result = calculate_consistency(AlignSeq, AlignSeq1);
    }
    return result;
}

//void Readfile(const string file_path, string& sequence) {
//    ifstream file(file_path);
//
//    if (file.is_open()) {
//        string line;
//        while (getline(file, line)) {
//            sequence += line;
//        }
//        cout << "File Content:\n" << sequence << endl;
//        file.close();
//    }
//    else {
//        cerr << "Unable to open file: " << file_path << endl;
//    }
//}

//int main() {
//
//	string seq1, seq2;
//    //const string file_path = "C:\\Users\\Administrator\\Desktop\\seq\\seq21690.txt"; //seq32199.txt //seq40739.txt \\seq.txt
//    //const string file_path1 = "C:\\Users\\Administrator\\Desktop\\seq\\seq121690.txt";   //seq132199.txt //seq140739.txt \\seq1.txt
//    //Readfile(file_path, seq1);
//    //Readfile(file_path1, seq2);
//	double consistency;
//	consistency = computeVarseqConsistency(seq1,seq2);
//    cout << "result:" << consistency << endl;
//	return 0;
//}
