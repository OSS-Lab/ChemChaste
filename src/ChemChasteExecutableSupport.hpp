/*
#ifndef EXECUTABLESUPPORT_HPP_
#define EXECUTABLESUPPORT_HPP_

#include <string>
#include "FileFinder.hpp"
#include "OutputFileHandler.hpp" // For out_stream

/**
 * Various helpful static methods for people writing their own executables
 * within the Chaste framework.
 *
 * Most executables will just need to call StandardStartup as the first thing in their
 * main() function, and FinalizePetsc before quitting.  The other methods allow you to
 * fine-tune what output is presented to users.
 */
class ExecutableSupport
{
public:
    /**
     * Initialise PETSc from the command line arguments.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void InitializePetsc(int* pArgc, char*** pArgv);

    /**
     * Display Chaste's copyright information.
     */
    static void ShowCopyright();

    /**
     * Output extra diagnostics when Chaste is launched in parallel.
     */
    static void ShowParallelLaunching();

    /**
     * Set the directory to which the files created by WriteMachineInfoFile,
     * and WriteProvenanceInfoFile will be written.
     * By default they will write to the CHASTE_TEST_OUTPUT folder.
     *
     * @param rOutputDirectory  the directory to write to
     */
    static void SetOutputDirectory(const std::string& rOutputDirectory);

    /**
     * Write to log files information about the machine that ran the code.
     * Each process will output its own file.
     *
     * @param fileBaseName base name of the file to write to
     */
    static void WriteMachineInfoFile(std::string fileBaseName);

    /**
     * Write to log files provenance information about this executable.
     * Each process will output its own file, named provenance_info.n.txt.
     */
    static void WriteProvenanceInfoFile();

    /**
     * Get information about library and compiler versions in XML-like format.
     *
     * @param rInfo a string to populate with library info.
     */
    static void GetBuildInfo(std::string& rInfo);

    /**
     * Call InitializePetsc, ShowCopyright, then ShowParallelLaunching.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StandardStartup(int* pArgc, char*** pArgv);

    /**
     * Call InitializePetsc, then ShowParallelLaunching. Omit ShowCopyright.
     *
     * @param pArgc  pointer to the number of arguments
     * @param pArgv  pointer to the argument list
     */
    static void StartupWithoutShowingCopyright(int* pArgc, char*** pArgv);

    /**
     * Display an error message to the user, on stderr.
     *
     * @param rMessage  the message to display
     * @param masterOnly  whether only the master process should display the error
     */
    static void PrintError(const std::string& rMessage, bool masterOnly=false);

    /**
     * Display an informative message to the user, on stdout.
     * Message is only displayed by master process
     * @param rMessage  the message to display
     */
    static void Print(const std::string& rMessage);

    /**
     * Shut down PETSc so we exit cleanly.
     */
    static void FinalizePetsc();

    /**
     * Standard exit codes for executables to return from main():
     * successful termination.
     */
    static const int EXIT_OK = 0;

    /**
     * Standard exit codes for executables to return from main():
     * exception thrown during execution.
     */
    static const int EXIT_ERROR = 1;

    /**
     * Standard exit codes for executables to return from main():
     * bad arguments passed on command line.
     */
    static const int EXIT_BAD_ARGUMENTS = 2;

private:
    /** The output directory to put machine provenance information into. */
    static FileFinder mOutputDirectory;
};

#endif /* EXECUTABLESUPPORT_HPP_ */
*/