package org.nextprot.api.core.service.export;

import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.core.domain.ChromosomeReport;
import org.nextprot.api.core.service.ChromosomeReportExportService;
import org.nextprot.api.core.service.ChromosomeReportService;
import org.nextprot.api.core.service.ChromosomeReportSummaryService;
import org.nextprot.api.core.service.OverviewService;
import org.nextprot.api.core.service.export.format.NextprotMediaType;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.stereotype.Service;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@Service
public class ChromosomeReportExportServiceImpl implements ChromosomeReportExportService {

	@Autowired
	private ChromosomeReportService chromosomeReportService;
	@Autowired
	private ChromosomeReportSummaryService chromosomeReportSummaryService;
	@Autowired
	private OverviewService overviewService;

	@Override
	public void exportChromosomeEntryReport(String chromosome, NextprotMediaType nextprotMediaType, OutputStream os) throws IOException {

		Optional<ChromosomeReportWriter> writer = ChromosomeReportWriter.valueOf(nextprotMediaType, os);

		if (writer.isPresent()) {
			writer.get().write(chromosomeReportService.reportChromosome(chromosome));
		}
		else {
			throw new NextProtException("cannot export chromosome "+chromosome+": " + "unsupported "+nextprotMediaType+" format");
		}
	}

	@Override
	public void exportHPPChromosomeEntryReport(String chromosome, NextprotMediaType nextprotMediaType, OutputStream os) throws IOException {

		Optional<HPPChromosomeReportWriter> writer = HPPChromosomeReportWriter.valueOf(nextprotMediaType, os, overviewService);

		if (writer.isPresent()) {
			writer.get().write(chromosomeReportService.reportChromosome(chromosome));
		}
		else {
			throw new NextProtException("cannot export hpp chromosome "+chromosome+": " + "unsupported "+nextprotMediaType+" format");
		}
	}

	@Override
	public void exportHPPChromosomeEntryReportCountByProteinExistence(OutputStream os) throws IOException {

		PrintWriter writer = new PrintWriter(os);

		Map<String, ChromosomeReport.Summary> summaries = chromosomeReportSummaryService.getChromosomeSummaries();

		writer.write(Stream.of(
				"chromosome	",
				"entry count",
				"protein level (PE1)",
				"transcript level (PE2)",
				"homology (PE3)",
				"predicted (PE4)",
				"uncertain (PE5)",
				"awaiting protein validation (P2+P3+P4)"
		).collect(Collectors.joining("\t")));

		writer.write("\n");

		for (String chromosome : ChromosomeReportService.getChromosomeNames()) {

			ChromosomeReport.Summary summary = summaries.get(chromosome);

			writer.write(Stream.of(
					chromosome,
					String.valueOf(summary.getEntryCount()),
					String.valueOf(summary.getProteinLevelEntryCount()),
					String.valueOf(summary.getTranscriptLevelEntryCount()),
					String.valueOf(summary.getHomologyLevelEntryCount()),
					String.valueOf(summary.getPredictedLevelEntryCount()),
					String.valueOf(summary.getUncertainLevelEntryCount()),
					String.valueOf(summary.getTranscriptLevelEntryCount()+summary.getHomologyLevelEntryCount()+summary.getPredictedLevelEntryCount())
			).collect(Collectors.joining("\t")));

			writer.write("\n");
		}

		writer.close();
	}
}
