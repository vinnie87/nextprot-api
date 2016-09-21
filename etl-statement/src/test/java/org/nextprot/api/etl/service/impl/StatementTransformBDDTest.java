package org.nextprot.api.etl.service.impl;

import static org.junit.Assert.fail;
import static org.nextprot.api.commons.constants.AnnotationCategory.VARIANT;

import java.util.Set;
import java.util.function.Predicate;

import org.junit.Assert;
import org.junit.Test;
import org.nextprot.api.commons.constants.AnnotationCategory;
import org.nextprot.api.commons.exception.NextProtException;
import org.nextprot.api.commons.utils.StringUtils;
import org.nextprot.api.etl.statement.StatementETLBaseUnitTest;
import org.nextprot.commons.statements.Statement;
import org.nextprot.commons.statements.StatementField;

import com.hp.hpl.jena.sparql.modify.request.Target;
import com.nextprot.api.annotation.builder.statement.TargetIsoformSerializer;



public class StatementTransformBDDTest extends StatementETLBaseUnitTest {


	/**
	 * It is not allowed to have a subject composed by variants in different genes 
	 */
	@Test
	public void shouldThrowAnExceptionWhenMultipleMutantsAreLocatedOnDifferentGenes() {
	
		try {
			StatementsExtractorLocalMockImpl sle = new StatementsExtractorLocalMockImpl();
			Set<Statement> rawStatements = sle.getStatementsForSourceForGeneName(null, "msh2-msh6-multiple-mutants-on-different-genes");

			statementETLServiceMocked.transformStatements(rawStatements);
			
			fail();
			
		}catch(NextProtException e){
			
			Assert.assertEquals("Mixing iso numbers for subjects is not allowed", e.getMessage());
			Assert.assertEquals(NextProtException.class, e.getClass());
			
		}
	
	}
	
	static class AnnotationCategoryPredicate implements Predicate<Statement>{
		
		private AnnotationCategory category = null;
		public AnnotationCategoryPredicate(AnnotationCategory category){
			this.category = category;
		}

		@Override
		public boolean test(Statement s) {
			String sCat = s.getValue(StatementField.ANNOTATION_CATEGORY);
			AnnotationCategory sCategory = AnnotationCategory.getDecamelizedAnnotationTypeName(StringUtils.camelToKebabCase(sCat));
			return sCategory.equals(category);
		}
	
	}
	
	/**
	 * Specification: A variant should always be propagated to all possible isoforms.
	 * If a variant can not be propagated to an isoform, then the 'phenotypic variation' should not exist for that isoform neiher.
	 * However the object annotation should be propagate to all isoforms (if not positional)
	 * 
	 * In this test we check if the propagation of the variant MSH6-p.Thr1219Asp is not propagated to the isoform 2, because the Isoform 2 can not contain this exons.
	 * Additionally we check that the phenotypic variation annotation is well propagated according to the variant.
	 * Finally we check that the object is propagate to all isoforms.
	 * 
	 */
	@Test
	public void shouldPropagateVariantsOnlyToMappableIsoforms() {

		StatementsExtractorLocalMockImpl sle = new StatementsExtractorLocalMockImpl();
		Set<Statement> rawStatements = sle.getStatementsForSourceForGeneName(null, "msh6-variant-on-iso1-but-not-on-iso2");

		//Variant 
		Set<Statement> mappedStatements =statementETLServiceMocked.transformStatements(rawStatements);
		
		Statement variantMappedStatement = mappedStatements.stream().filter(new AnnotationCategoryPredicate(VARIANT)).findFirst().orElseThrow(RuntimeException::new);
		
		String variantMappedStatementIsoformsJson = variantMappedStatement.getValue(StatementField.TARGET_ISOFORMS);
		String isoPropagationsWithoutIsoform2 = "[{\"isoformAccession\":\"NX_P52701-1\",\"specificity\":\"BY_DEFAULT\",\"begin\":1219,\"end\":1219,\"name\":\"MSH6-isoGTBP-N-p.Thr1219Asp\"},{\"isoformAccession\":\"NX_P52701-3\",\"specificity\":\"BY_DEFAULT\",\"begin\":1089,\"end\":1089,\"name\":\"MSH6-iso3-p.Thr1089Asp\"},{\"isoformAccession\":\"NX_P52701-4\",\"specificity\":\"BY_DEFAULT\",\"begin\":917,\"end\":917,\"name\":\"MSH6-iso4-p.Thr917Asp\"}]";
		
		Assert.assertEquals(variantMappedStatementIsoformsJson, isoPropagationsWithoutIsoform2);
	
		//Phenotypic variation
		Statement phentypicMappedStatement = mappedStatements.stream().filter(new AnnotationCategoryPredicate(AnnotationCategory.PHENOTYPIC_VARIATION)).findFirst().orElseThrow(RuntimeException::new);
		String phenotypicMappedStatementIsoformJson = phentypicMappedStatement.getValue(StatementField.TARGET_ISOFORMS);
		Assert.assertEquals(TargetIsoformSerializer.deSerializeFromJsonString(phenotypicMappedStatementIsoformJson).size(), 3);
		
		String phenotypicWithoutIsoform2 = "[{\"isoformAccession\":\"NX_P52701-1\",\"specificity\":\"BY_DEFAULT\",\"begin\":null,\"end\":null,\"name\":\"MSH6-isoGTBP-N-p.Thr1219Asp\"},{\"isoformAccession\":\"NX_P52701-3\",\"specificity\":\"BY_DEFAULT\",\"begin\":null,\"end\":null,\"name\":\"MSH6-iso3-p.Thr1089Asp\"},{\"isoformAccession\":\"NX_P52701-4\",\"specificity\":\"BY_DEFAULT\",\"begin\":null,\"end\":null,\"name\":\"MSH6-iso4-p.Thr917Asp\"}]";
		
		Assert.assertEquals(TargetIsoformSerializer.deSerializeFromJsonString(phenotypicWithoutIsoform2).size(), 3);
		Assert.assertEquals(phenotypicMappedStatementIsoformJson, phenotypicWithoutIsoform2);
		

		//Object
		Statement objectStatement = mappedStatements.stream().filter(new AnnotationCategoryPredicate(AnnotationCategory.GO_MOLECULAR_FUNCTION)).findFirst().orElseThrow(RuntimeException::new);
		String objectMappedStatementIsoformJson = phentypicMappedStatement.getValue(StatementField.TARGET_ISOFORMS);
		
		//TODO what to do??? System.err.println(objectMappedStatementIsoformJson);
	}
	
	
	
	@Test
	public void shouldThrowAnExceptionIfFeatureNameDoesNotCorrespondToNextprotAccession() {
		
		
	}

}
