package org.nextprot.api.core.service.export.io;

import org.nextprot.api.core.dao.EntityName;
import org.nextprot.api.core.domain.ChromosomeReport;
import org.nextprot.api.core.domain.EntryReport;
import org.nextprot.api.core.domain.Overview;
import org.nextprot.api.core.domain.ProteinExistenceLevel;
import org.nextprot.api.core.service.OverviewService;
import org.nextprot.api.core.service.export.HPPChromosomeReportWriter;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;


/**
 * Writes a {@code ChromosomeReport} in TXT format
 *
 * Created by fnikitin on 19.04.17.
 */
public class HPPChromosomeReportTXTWriter implements HPPChromosomeReportWriter {

    private PrintWriter writer;
    private final OverviewService overviewService;

    public HPPChromosomeReportTXTWriter(OutputStream os, OverviewService overviewService) {

        this.overviewService = overviewService;
        this.writer = new PrintWriter(os);
    }

    @Override
    public void write(ChromosomeReport chromosomeReport) throws IOException {

        writer.write(String.format(buildHeaderFormat(),
                Arrays.asList("neXtProt AC", "Gene name(s)", "Protein existence", "Proteomics", "Antibody").toArray()));

        Map<String, EntryReport> groupedByAccession = chromosomeReport.getEntryReports().stream()
                .collect(Collectors.toMap(
                        EntryReport::getAccession,
                        er -> er,
                        (er1, er2) -> er1) // 1. keep one entry report for each entry accession
                );

        List<String> sortedAccessions = new ArrayList<>(groupedByAccession.keySet()).stream()
                .sorted()                  // 2. order by entry accession
                .collect(Collectors.toList());

        for (String accession : sortedAccessions) {

            EntryReport er = groupedByAccession.get(accession);
            write(er, overviewService.findOverviewByEntry(er.getAccession()));
        }

        writer.close();
    }

    private void write(EntryReport report, Overview overview) throws IOException {

        writer.write(String.format(buildRowFormat(), extractValues(report, overview).toArray()));
    }

    private List<String> extractValues(EntryReport entryReport, Overview overview) {

        return Arrays.asList(
                entryReport.getAccession(),
                getMainEntityNames(overview.getGeneNames()),
                ProteinExistenceLevel.valueOfString(entryReport.getProteinExistence()).getDescription(),
                (entryReport.isProteomics()) ? "yes" : "no",
                (entryReport.isAntibody()) ? "yes" : "no"
        );
    }

    private static String buildHeaderFormat() {

        return "%-12s %-12s %-25s %-10s %-8s%n";
    }

    private static String buildRowFormat() {

        return "%-12s %-12s %-25s %-10s %-8s%n";
    }

    private static String getMainEntityNames(List<EntityName> entityNameList) {

        return entityNameList.stream()
                .filter(EntityName::isMain)
                .map(EntityName::getName)
                .collect(Collectors.joining(";"));
    }
}
