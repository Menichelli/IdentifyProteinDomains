/**
 * 
 */
package tools;

import java.util.HashSet;
import java.util.Set;

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
	
}
