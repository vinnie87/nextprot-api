package org.nextprot.api.web;

import org.junit.Before;
import org.junit.Test;
import org.nextprot.api.commons.utils.XMLPrettyPrinter;
import org.nextprot.api.web.dbunit.base.mvc.WebIntegrationBaseTest;
import org.nextprot.api.web.service.ExportService;
import org.nextprot.api.web.service.impl.writer.EntryStreamWriter;
import org.nextprot.api.web.service.impl.writer.EntryXMLStreamWriter;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ActiveProfiles;

import javax.xml.XMLConstants;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;
import javax.xml.validation.Validator;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;

import static org.junit.Assert.fail;

@ActiveProfiles()
public class XSDValidationTest extends WebIntegrationBaseTest {


	@Autowired
	private ExportService exportService;

	@Before
	public void clearRepository() {
		exportService.clearRepository();
	}

	@Test
	public void shouldValidateXMLFilewithXSD() {

		
		Schema schema;
		try {

			SchemaFactory factory = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
			schema = factory.newSchema(new StreamSource(new File("src/main/webapp/nextprot-export-v2.xsd")));

			File f = new File("tmp.xml");
			StreamSource xmlFile = new StreamSource(f);
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
		
            EntryStreamWriter<?> writer = new EntryXMLStreamWriter(baos, "entry");
			exportService.streamResults(writer, "entry", Arrays.asList(new String[] { "NX_Q15858" }));
			//exportService.streamResults(writer, "entry", Arrays.asList(new String[] { "NX_Q6PIU2" })); 

			XMLPrettyPrinter prettyPrinter = new XMLPrettyPrinter();
			//System.err.println(baos.toString());
			
			String prettyXml = prettyPrinter.prettify(baos.toString());
			PrintWriter out = new PrintWriter(f);
			out.print(prettyXml);
			out.close();
			
			// instance document
			Validator validator = schema.newValidator();
			// validate the DOM tree
			validator.validate(xmlFile);

			f.delete();

		} catch (Exception e) {
			e.printStackTrace();
			fail();
			
		}

	}
}
