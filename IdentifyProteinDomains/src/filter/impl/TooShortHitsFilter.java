/**
 * 
 */
package filter.impl;

import java.util.HashSet;
import java.util.Set;

import model.BlastHit;
import model.PutativeDomain;
import filter.Ifilter;

/**
 * @author christophe
 *
 */
public class TooShortHitsFilter implements Ifilter {
	
	private static Ifilter instance;
	
	private TooShortHitsFilter() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new TooShortHitsFilter():instance;
		return instance;
	}
	
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		
		Set<BlastHit> hits = domain.getBlastHits();
		domain.setBlastHits(new HashSet<BlastHit>());
		
		for(BlastHit hit : hits) {
			if((hit.getqEnd()-hit.getqStart()) >= (domain.getDomainEnd()-domain.getDomainStart()*.5)) domain.addBlastHit(hit);
		}
		
		ret = domain;
		
		return ret;
	}
	
}
