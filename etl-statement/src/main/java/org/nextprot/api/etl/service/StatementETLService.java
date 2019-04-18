package org.nextprot.api.etl.service;


import org.nextprot.api.etl.NextProtSource;

import java.io.IOException;

public interface StatementETLService {
    
	String etlStatements(NextProtSource source, String release, boolean load) throws IOException;
}
