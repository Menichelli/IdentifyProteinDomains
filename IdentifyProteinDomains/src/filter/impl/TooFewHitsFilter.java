/**
 * 
 */
package filter.impl;

import filter.Ifilter;
import global.Global;
import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class TooFewHitsFilter implements Ifilter {

	private static Ifilter instance;
	
	private TooFewHitsFilter() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new TooFewHitsFilter():instance;
		return instance;
	}
	
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		if(domain.getBlastHits().size()>=Global.NB_SEQ_INTERSECT) {
			ret = domain;
		}
		return ret;
	}
}
