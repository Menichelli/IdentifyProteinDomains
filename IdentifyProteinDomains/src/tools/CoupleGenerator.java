/**
 * 
 */
package tools;

import java.util.HashSet;
import java.util.Set;

import model.CouplePfamPutative;
import model.CouplePutativePutative;
import model.PfamFamily;
import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class CoupleGenerator {

	private CoupleGenerator() {}

	/**
	 * Pour les couples Pfam-Putative
	 * @param pDomains
	 * @param pfamFamilies
	 * @return
	 */
	public static Set<CouplePfamPutative> getCouplePfamPutative(Set<PutativeDomain> pDomains, Set<PfamFamily> pfamFamilies) {
		Set<CouplePfamPutative> ret = new HashSet<CouplePfamPutative>();

		for(PutativeDomain putativeDomain : pDomains) {
			for(PfamFamily pf : pfamFamilies) {
				if(pf.isFoundOn(putativeDomain.getQueryName()+"_"+putativeDomain.getQuerySpecies())) {
					ret.add(new CouplePfamPutative(pf, putativeDomain));
				}
			}
		}
		return ret;
	}

	/**
	 * Pour les couples Putative-Putative
	 * @param pDomains
	 * @return
	 */
	public static Set<CouplePutativePutative> getCouplePutativePutative(Set<PutativeDomain> pDomains) {
		Set<CouplePutativePutative> ret = new HashSet<CouplePutativePutative>();

		for(PutativeDomain p1 : pDomains) {
			for(PutativeDomain p2 : pDomains) {
				if(p1.getQueryName().equals(p2.getQueryName()) && p1.getDomainStart()!=p2.getDomainStart()) { //check if they are found on the same protein and if they are differents
					ret.add(new CouplePutativePutative(p1, p2));
				}
			}
		}

		return ret;
	}

}
