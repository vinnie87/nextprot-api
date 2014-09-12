package org.nextprot.api.dbunit;

import static org.springframework.test.web.servlet.setup.MockMvcBuilders.webAppContextSetup;

import org.junit.Before;
import org.nextprot.api.user.domain.UserApplication;
import org.nextprot.api.user.security.UserApplicationKeyGenerator;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.security.web.FilterChainProxy;
import org.springframework.test.annotation.DirtiesContext;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.web.WebAppConfiguration;
import org.springframework.test.web.servlet.MockMvc;
import org.springframework.web.context.WebApplicationContext;

/**
 * Base class for dbunit tests using the spring-test-dbunit framework http://springtestdbunit.github.io/
 * Transactions are rollback and dev profile is activated by default
 * Dev profile includes database connection to the dev database
 * 
 * @RunWith(SpringJUnit4ClassRunner.class)
 * @ContextConfiguration("classpath:api-servlet-test.xml")
 * @ActiveProfiles("test")
 * @author dteixeira
 */

@WebAppConfiguration
@ContextConfiguration("classpath:META-INF/spring/web-context.xml")
@DirtiesContext
public abstract class MVCBaseIntegrationTest extends AbstractIntegrationBaseTest {

	@Autowired
	protected UserApplicationKeyGenerator keyGenerator;
	
	@Autowired
	protected WebApplicationContext wac;

	@Autowired
    private FilterChainProxy springSecurityFilterChain;
	
	protected MockMvc mockMvc;

	@Before
	public void setup() {
		this.mockMvc = webAppContextSetup(this.wac).addFilters(this.springSecurityFilterChain).build();
	}
	
	
	protected String generateTestToken(){
		
		UserApplication ua = new UserApplication();
		ua.setName("unit-test-application");
		
		return keyGenerator.generateToken(ua);
	}
	

}
