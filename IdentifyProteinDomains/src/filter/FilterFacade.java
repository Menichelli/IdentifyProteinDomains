/**
 * 
 */
package filter;

import java.lang.reflect.Method;
import java.util.HashSet;
import java.util.Set;

import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class FilterFacade {

	private FilterFacade() {}
	
	public static Set<PutativeDomain> applyFilter(Class<? extends Ifilter> filter, Set<PutativeDomain> domains) throws Exception {
		Set<PutativeDomain> ret = new HashSet<PutativeDomain>();
		PutativeDomain tmp;
		for(PutativeDomain domain : domains) {
			Method m = filter.getMethod("filter", PutativeDomain.class);
			Ifilter obj = (Ifilter)(filter.getMethod("getInstance",new Class[0]).invoke(null));
			tmp = (PutativeDomain) m.invoke(obj, domain);
			if(tmp!=null) ret.add(tmp);
		}
		return ret;
	}

}
