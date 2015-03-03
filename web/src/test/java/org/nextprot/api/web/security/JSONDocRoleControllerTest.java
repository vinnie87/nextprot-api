package org.nextprot.api.web.security;

import static org.junit.Assert.assertTrue;
import static org.junit.Assert.assertFalse;
import static org.springframework.test.web.servlet.request.MockMvcRequestBuilders.get;
import static org.springframework.test.web.servlet.result.MockMvcResultMatchers.status;

import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import org.junit.Test;
import org.nextprot.api.web.dbunit.base.mvc.MVCBaseSecurityTest;
import org.springframework.http.MediaType;

/**
 * Tests GET, PUT, POST, DELETE for 3 different scenarios (anonymous, owner and other logged user) 
 * @author dteixeira
 *
 */
public class JSONDocRoleControllerTest extends MVCBaseSecurityTest {

	@Test
	public void sheldonShouldBeAbleToSeeHisSuperGeniousQuery() throws Exception {

		String sheldonToken = generateTokenWithExpirationDate("Sheldon", 1, TimeUnit.DAYS, Arrays.asList("ROLE_USER"));

		String responseString = this.getJSONDocByUser(sheldonToken);

		// Admin group does not exist
		assertFalse(this.isMatchRegExpGroup(responseString, "Admin"));
		// User and "" groups exist
		assertTrue(this.isMatchRegExpGroup(responseString, "User"));
		assertTrue(this.isMatchRegExpGroup(responseString, ""));
		
		// Check presence of User subgroups
		assertTrue(this.containsWithKeyValue(responseString, "name", "User"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Application"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Protein Lists"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Queries"));
	}

	@Test
	public void adminShouldBeAbleToSeeAllData() throws Exception {

		String adminToken = generateTokenWithExpirationDate("Admin", 1, TimeUnit.DAYS, Arrays.asList("ROLE_ADMIN"));

		String responseString = this.getJSONDocByUser(adminToken);

		// All groups exist
		assertTrue(this.isMatchRegExpGroup(responseString, "Admin"));
		assertTrue(this.isMatchRegExpGroup(responseString, "User"));
		assertTrue(this.isMatchRegExpGroup(responseString, ""));

		// Check presence of User subgroups
		assertTrue(this.containsWithKeyValue(responseString, "name", "User"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Application"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Protein Lists"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Queries"));
	}

	@Test
	public void anonymousShouldBeAbleToSeeSimpleData() throws Exception {

		String adminToken = generateTokenWithExpirationDate("Anonymous", 1, TimeUnit.DAYS, Arrays.asList("ROLE_ANONYMOUS"));

		String responseString = this.getJSONDocByUser(adminToken);

		// Admin group does not exist
		assertFalse(this.isMatchRegExpGroup(responseString, "Admin"));
		// User and "" groups exist 
		assertTrue(this.isMatchRegExpGroup(responseString, "User"));
		assertTrue(this.isMatchRegExpGroup(responseString, ""));

		// Check presence/absence of User subgroups
		assertFalse(this.containsWithKeyValue(responseString, "name", "User"));
		assertFalse(this.containsWithKeyValue(responseString, "name", "User Application"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Protein Lists"));
		assertTrue(this.containsWithKeyValue(responseString, "name", "User Queries"));

		// Check that does not contain any "modification" verbs
		assertFalse(this.containsWithKeyValue(responseString, "verb", "POST"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "PATCH"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "PUT"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "DELETE"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "HEAD"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "OPTIONS"));
		assertFalse(this.containsWithKeyValue(responseString, "verb", "TRACE"));
	}

	/**
	 * Get MVC mock for jsondoc with the providen user.
	 */
	private String getJSONDocByUser(String user) throws Exception {
		return this.mockMvc.perform(get("/jsondoc").contentType(MediaType.APPLICATION_JSON).
				header("Authorization", "Bearer " + user).accept(MediaType.APPLICATION_JSON)).
				andExpect(status().isOk()).andReturn().getResponse().getContentAsString();
	}

	/**
	 * Returns true if and only if the provided string contains the specified string formed by key and value 
	 * (for instance, '"name":"User Application"').
	 */
	private boolean containsWithKeyValue(String string, String key, String value) {
		return string.contains("\"" + key + "\":\"" + value +"\"");
	}
	
	/**
	 * Returns true if and only if the provided string contains the specified string of a JSONDoc group 
	 * (for instance, '"Admin":[' not succeeded by ']').
	 */
	private boolean isMatchRegExpGroup(String string, String groupName) {
		return  string.matches(".*\""+groupName+"\":\\[[^\\]].*");
	}
}
