/*
 * msa2db
 * Converts multiple MSA files (FASTA/A3M format) to an msaDB database
 */
#include "FileUtil.h"
#include "DBWriter.h"
#include "Debug.h"
#include "Util.h"
#include "KSeqWrapper.h"
#include "FastSort.h"

#ifdef OPENMP
#include <omp.h>
#endif

int msa2db(int argc, const char **argv, const Command& command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_VARIADIC, 0);
    par.printParameters(command.cmd, argc, argv, *command.params);

    std::vector<std::string> filenames(par.filenames);
    std::string dataFile = filenames.back();
    filenames.pop_back();

    if (filenames.empty()) {
        Debug(Debug::ERROR) << "No input files provided\n";
        EXIT(EXIT_FAILURE);
    }

    // Handle TSV input file
    if (Util::endsWith(".tsv", filenames[0])) {
        if (filenames.size() > 1) {
            Debug(Debug::ERROR) << "Only one tsv file can be given\n";
            EXIT(EXIT_FAILURE);
        }
        std::string tsv = filenames.back();
        filenames.pop_back();

        FILE* file = FileUtil::openFileOrDie(tsv.c_str(), "r", true);
        char* line = NULL;
        size_t len = 0;
        ssize_t read;
        while ((read = getline(&line, &len, file)) != -1) {
            if (line[read - 1] == '\n') {
                line[read - 1] = '\0';
                read--;
            }
            filenames.push_back(line);
        }
        free(line);
        fclose(file);
    }

    // Consistent order
    SORT_SERIAL(filenames.begin(), filenames.end(), [](const std::string &a, const std::string &b) {
        return FileUtil::baseName(a) < FileUtil::baseName(b);
    });
    
    for (size_t i = 0; i < filenames.size(); i++) {
        if (FileUtil::directoryExists(filenames[i].c_str()) == true) {
            Debug(Debug::ERROR) << "File " << filenames[i] << " is a directory\n";
            EXIT(EXIT_FAILURE);
        }
    }

    std::string indexFile = dataFile + ".index";

    DBWriter writer(dataFile.c_str(), indexFile.c_str(), 1, par.compressed, Parameters::DBTYPE_MSA_DB);
    writer.open();

    // Create lookup file
    std::string lookupFile = dataFile + ".lookup";
    FILE* lookupFileHandle = FileUtil::openAndDelete(lookupFile.c_str(), "w");
    if (lookupFileHandle == NULL) {
        Debug(Debug::ERROR) << "Cannot open lookup file " << lookupFile << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    // Create source file
    std::string sourceFile = dataFile + ".source";
    FILE *source = fopen(sourceFile.c_str(), "w");
    if (source == NULL) {
        Debug(Debug::ERROR) << "Cannot open " << sourceFile << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Converting MSA files\n";
    Debug::Progress progress;

    unsigned int entries_num = 0;
    std::string result;
    result.reserve(10 * 1024 * 1024);
    std::string lookupBuffer;
    char buffer[4096];
    DBReader<unsigned int>::LookupEntry lookupEntry;

    for (size_t fileIdx = 0; fileIdx < filenames.size(); fileIdx++) {
        progress.updateProgress();
        
        std::string sourceName = FileUtil::baseName(filenames[fileIdx]);
        
        // Write to source file
        size_t len = snprintf(buffer, sizeof(buffer), "%zu\t%s\n", fileIdx, sourceName.c_str());
        size_t written = fwrite(buffer, sizeof(char), len, source);
        if (written != len) {
            Debug(Debug::ERROR) << "Cannot write to source file " << sourceFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        // Read the MSA file using KSeqWrapper
        KSeqWrapper* kseq = KSeqFactory(filenames[fileIdx].c_str());
        if (kseq == NULL) {
            Debug(Debug::ERROR) << "Cannot open file " << filenames[fileIdx] << "\n";
            EXIT(EXIT_FAILURE);
        }

        result.clear();
        std::string firstSeqName;
        bool hasEntries = false;

        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry &e = kseq->entry;
            if (e.name.l == 0) {
                Debug(Debug::WARNING) << "Invalid entry in file " << filenames[fileIdx] << "\n";
                continue;
            }

            hasEntries = true;

            // Store the first sequence name as identifier
            if (firstSeqName.empty()) {
                firstSeqName = std::string(e.name.s, e.name.l);
            }

            // Write FASTA format
            result.append(">");
            result.append(e.name.s, e.name.l);
            if (e.comment.l > 0) {
                result.append(" ");
                result.append(e.comment.s, e.comment.l);
            }
            result.append("\n");
            result.append(e.sequence.s, e.sequence.l);
            result.append("\n");
        }
        delete kseq;

        if (!hasEntries) {
            Debug(Debug::WARNING) << "File " << filenames[fileIdx] << " is empty or invalid, skipping\n";
            continue;
        }

        // Write the MSA entry to the database
        unsigned int id = par.identifierOffset + entries_num;
        writer.writeData(result.c_str(), result.length(), id, 0);

        // Write lookup entry
        lookupEntry.id = id;
        lookupEntry.entryName = firstSeqName.empty() ? sourceName : firstSeqName;
        lookupEntry.fileNumber = fileIdx;
        lookupBuffer.clear();
        DBReader<unsigned int>::lookupEntryToBuffer(lookupBuffer, lookupEntry);
        written = fwrite(lookupBuffer.c_str(), sizeof(char), lookupBuffer.size(), lookupFileHandle);
        if (written != lookupBuffer.size()) {
            Debug(Debug::ERROR) << "Cannot write to lookup file " << lookupFile << "\n";
            EXIT(EXIT_FAILURE);
        }

        entries_num++;
    }

    Debug(Debug::INFO) << "\n";

    if (fclose(source) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << sourceFile << "\n";
        EXIT(EXIT_FAILURE);
    }

    if (fclose(lookupFileHandle) != 0) {
        Debug(Debug::ERROR) << "Cannot close file " << lookupFile << "\n";
        EXIT(EXIT_FAILURE);
    }

    writer.close(true);

    if (entries_num == 0) {
        Debug(Debug::ERROR) << "No valid MSA files found in input\n";
        EXIT(EXIT_FAILURE);
    }

    Debug(Debug::INFO) << "Created MSA database with " << entries_num << " entries\n";

    return EXIT_SUCCESS;
}
