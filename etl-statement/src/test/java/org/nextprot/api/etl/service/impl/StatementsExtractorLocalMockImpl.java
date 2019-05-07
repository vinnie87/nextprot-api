package org.nextprot.api.etl.service.impl;

import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.etl.StatementSource;
import org.nextprot.api.etl.service.StatementDictionary;
import org.nextprot.api.etl.service.StatementExtractorService;
import org.nextprot.commons.statements.Statement;
import org.nextprot.commons.statements.reader.JsonReader;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.Set;

public class StatementsExtractorLocalMockImpl implements StatementExtractorService {

	@Override
	public Collection<Statement> getStatementsFromJsonFile(StatementSource source, String release, String jsonFilename) throws IOException {

		StatementDictionary sd = new StatementDictionary();
		String content = sd.getStatements(jsonFilename);
		String removedComments = content.replaceAll("((['\"])(?:(?!\\2|\\\\).|\\\\.)*\\2)|\\/\\/[^\\n]*|\\/\\*(?:[^*]|\\*(?!\\/))*\\*\\/", "$1");

		return new JsonReader(source.getSpecifications()).readStatements(new ByteArrayInputStream(removedComments.getBytes(StandardCharsets.UTF_8)));
	}

	@Override
	public Set<Statement> getStatementsForSource(StatementSource source, String release) {
		throw new NextProtException("Method not supported");
	}
}