package org.nextprot.api.blast.service;

import com.fasterxml.jackson.annotation.JsonInclude;
import org.nextprot.api.blast.controller.SystemCommandExecutor;
import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.commons.utils.ExceptionWithReason;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.nio.file.Files;
import java.util.List;
import java.util.Objects;
import java.util.logging.Logger;

/**
 * A base object for blast suite programs
 */
public abstract class BlastProgram<I, O, C extends BlastProgram.Config> {

    private static final Logger LOGGER = Logger.getLogger(BlastProgram.class.getName());

    private final String name;
    private final C config;

    protected BlastProgram(String name, C config) {

        Objects.requireNonNull(config, "missing blast configuration");

        this.name = name;
        this.config = config;
    }

    public O run(I input) throws ExceptionWithReason {

        try {
            // pre process
            File fastaFile = constructFastaFile(input);
            List<String> commandLine = buildCommandLine(config, fastaFile);
            preConfig(input, config);

            O out = process(commandLine);

            // post process
            tearDownFastaFile(fastaFile);
            postConfig(config);

            return out;
        } catch (Exception e) {

            throw new NextProtException("could not run "+name, e);
        }
    }

    private File constructFastaFile(I input) throws Exception {

        File tmpQueryFastaFile = File.createTempFile(name, ".fasta");

        PrintWriter pw = new PrintWriter(tmpQueryFastaFile);
        writeFastaInput(pw, input);
        pw.close();

        LOGGER.info("create temporary file "+tmpQueryFastaFile.getName());

        return tmpQueryFastaFile;
    }

    private O process(List<String> commandLine) throws IOException, InterruptedException, ExceptionWithReason {

        SystemCommandExecutor commandExecutor = new SystemCommandExecutor(commandLine);
        commandExecutor.executeCommand();

        String stderr = commandExecutor.getStandardErrorFromCommand().toString();

        if (!stderr.isEmpty()) {

            ExceptionWithReason ewr = new ExceptionWithReason();
            ewr.getReason().addCause(name+ " exception", stderr.replace("\n", " "));
            ewr.getReason().setMessage("Error while executing "+name);

            throw ewr;
        }

        return buildOutputFromStdout(commandExecutor.getStandardOutputFromCommand().toString());
    }

    private void tearDownFastaFile(File fastaFile) throws Exception {

        if (!config.isDebugMode()) {

            LOGGER.info("delete temporary file "+fastaFile.getName());
            Files.deleteIfExists(fastaFile.toPath());
        }
    }

    protected void preConfig(I input, C config) { }

    protected void postConfig(C config) {

        config.setDebugMode(null);
        config.unsetPathes();
    }

    /**
     * Build command line to execute by blast program
     * @param config input configuration
     * @param fastaFile the input fasta file
     * @return command line list
     */
    protected abstract List<String> buildCommandLine(C config, File fastaFile);

    /**
     * Build output object from command stdout
     * @param stdout the command output
     * @return an object containing output
     */
    protected abstract O buildOutputFromStdout(String stdout) throws IOException;

    /**
     * Write fasta entries from input object
     *
     * @param pw the writer object
     * @param input the input to extract entries from
     */
    protected abstract void writeFastaInput(PrintWriter pw, I input);

    /**
     * Configuration object for program execution
     */
    @JsonInclude(JsonInclude.Include.NON_NULL)
    public static class Config implements Serializable {

        private String binPath;
        private String nextprotBlastDbPath;
        private Boolean isDebugMode = false;

        public Config(String binPath, String nextprotBlastDbPath) {

            if (binPath == null)
                throw new NextProtException("Internal error: bin path is missing");

            if (nextprotBlastDbPath == null)
                throw new NextProtException("Internal error: nextprot blast db path is missing");

            this.binPath = binPath;
            this.nextprotBlastDbPath = nextprotBlastDbPath;
        }

        public String getNextprotBlastDbPath() {
            return nextprotBlastDbPath;
        }

        public String getBinPath() {
            return binPath;
        }

        public Boolean isDebugMode() {
            return isDebugMode;
        }

        public void setDebugMode(Boolean debugMode) {
            isDebugMode = debugMode;
        }

        public void unsetPathes() {

            binPath = null;
            nextprotBlastDbPath = null;
        }
    }

    /**
     * Write a fasta entry
     * @param pw
     * @param header
     * @param sequence
     */
    static void writeFastaEntry(PrintWriter pw, String header, String sequence) {

        pw.write(">");
        pw.write(header);
        pw.write("\n");
        pw.write(sequence);
        pw.write("\n");
    }
}
