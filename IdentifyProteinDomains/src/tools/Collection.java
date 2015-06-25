/**
 * 
 */
package tools;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class Collection {
	
	private Collection() {}
	
	public static <E> int intersectionSize(Set<E> set1, Set<E> set2) {
		Set<E> tmp = new HashSet<E>();
		tmp.addAll(set1);
		tmp.addAll(set2);
		return set1.size() + set2.size() - tmp.size();
	}
	
	public static Map<String, Set<PutativeDomain>> clonePutativeDomains(Map<String, Set<PutativeDomain>> map) {
		Map<String, Set<PutativeDomain>> ret = new HashMap<String, Set<PutativeDomain>>();
		
		Set<PutativeDomain> tmp = new HashSet<PutativeDomain>();
		for(String key : map.keySet()) {
			tmp.clear();
			for(PutativeDomain val : map.get(key)) {
				tmp.add((PutativeDomain)val.clone());
			}
			ret.put(key,tmp);
		}
		
		return ret;
	}
}
