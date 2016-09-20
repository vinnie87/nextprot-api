package org.nextprot.api.web.controller;

import org.jsondoc.core.annotation.Api;
import org.jsondoc.core.annotation.ApiMethod;
import org.jsondoc.core.annotation.ApiPathParam;
import org.jsondoc.core.pojo.ApiVerb;
import org.nextprot.api.commons.service.MasterIdentifierService;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.http.MediaType;
import org.springframework.stereotype.Controller;
import org.springframework.web.bind.annotation.PathVariable;
import org.springframework.web.bind.annotation.RequestMapping;
import org.springframework.web.bind.annotation.RequestMethod;
import org.springframework.web.bind.annotation.ResponseBody;

import java.util.ArrayList;
import java.util.List;

@Controller
@Api(name = "Master Identifiers", description = "Retrieves nextProt idenfitiers")
public class MasterIdentifierController {

	@Autowired	private MasterIdentifierService masterIdentifierService;

	@ApiMethod(path = "/master-identifiers", verb = ApiVerb.GET, description = "Retrieve all neXtProt entry accession numbers", produces = MediaType.APPLICATION_JSON_VALUE)
	@RequestMapping(value = "/master-identifiers", method = { RequestMethod.GET })
	public List<String> masterIdentifiers() {
		return new ArrayList<>(masterIdentifierService.findUniqueNames());
	}
	
	
	@ApiMethod(path = "/master-identifiers/chromosome/{chromosome}", verb = ApiVerb.GET, description = "Retrieve all neXtProt entry accession numbers of the given chromosome", produces = MediaType.APPLICATION_JSON_VALUE)
	@RequestMapping(value = "/master-identifiers/chromosome/{chromosome}", method = { RequestMethod.GET })
	public List<String> masterIdentifiersPerChromosome(
			@ApiPathParam(name = "chromosome", description = "The chromosome number or name (X,Y..)",  allowedvalues = { "master-isoform-mapping"}) @PathVariable("chromosome")  String chromosome) {
		return new ArrayList<>(masterIdentifierService.findUniqueNamesOfChromosome(chromosome));
	}

	
	@ApiMethod(path = "/master-identifiers/gene/{geneName}", verb = ApiVerb.GET, description = "Retrieves the entry accession number(s) corresponding to the given gene name", produces = MediaType.APPLICATION_JSON_VALUE)
	@RequestMapping(value = "/master-identifiers/gene/{geneName}", method = { RequestMethod.GET })
	@ResponseBody
	public List<String> masterIdentifierByGeneName(
			@ApiPathParam(name = "geneName", description = "The gene name",  allowedvalues = { "INSR"}) @PathVariable("geneName")  String geneName) {
		return new ArrayList<>(masterIdentifierService.findEntryAccessionByGeneName(geneName));
	}

}
