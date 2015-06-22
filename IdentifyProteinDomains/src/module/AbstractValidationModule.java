/**
 * 
 */
package module;

import java.util.Map;
import java.util.Set;

import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public abstract class AbstractValidationModule extends Thread implements IValidationModule {
	
	protected Map<String,Set<PutativeDomain>> putativeDomainsByProt;
	protected Map<String,PutativeDomain> mapIdPutativeDomain;
	
	public AbstractValidationModule(Map<String,Set<PutativeDomain>> putativeDomainsByProt, Map<String,PutativeDomain> mapIdPutativeDomain) {
		super();
		this.putativeDomainsByProt = putativeDomainsByProt;
		this.mapIdPutativeDomain = mapIdPutativeDomain;
	}
	
	public abstract void run();
	
}
