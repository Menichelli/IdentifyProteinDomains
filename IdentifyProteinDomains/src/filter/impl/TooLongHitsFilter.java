/**
 * 
 */
package filter.impl;

import java.util.HashSet;
import java.util.Set;

import model.BlastHit;
import model.PutativeDomain;
import filter.Ifilter;
import global.Global;

/**
 * @author christophe
 *
 */
public class TooLongHitsFilter implements Ifilter {
	
	private static Ifilter instance;
	
	private TooLongHitsFilter() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new TooLongHitsFilter():instance;
		return instance;
	}
	
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		
		Set<BlastHit> hits = domain.getBlastHits();
		domain.setBlastHits(new HashSet<BlastHit>());
		
		for(BlastHit hit : hits) {
			if(isIn(hit, domain)) domain.addBlastHit(hit);
		}
		
		ret = domain;
		
		return ret;
	}

	
	
	
	//--------------------------------------------------------------------------------------------//
	//--------------------------------------------------------------------------------------------//
	//--------------------------------------------------------------------------------------------//
	
	/**
	 * Determine if the blasthit is in a given putative domain
	 * @param bh
	 * @param pd
	 * @return
	 */
	private static boolean isIn(BlastHit bh, PutativeDomain pd) {
		boolean ret = false;

		if(getOverlapSize(bh,pd) >= (bh.getqEnd()-bh.getqStart())*Global.OVERLAP_RATE_MIN) { //if overlap by at least X%
			ret = true;
		}

		return ret;
	}

	private static int getOverlapSize(BlastHit bh, PutativeDomain pd) {
		int ret = 0;
		if(bh.getqStart() < pd.getDomainStart() && bh.getqEnd() > pd.getDomainEnd()) {
			ret = pd.getDomainEnd() - pd.getDomainStart() + 1;
		} else if(bh.getqStart() < pd.getDomainEnd() && bh.getqEnd() > pd.getDomainEnd()) {
			ret = pd.getDomainEnd() - Math.max(bh.getqStart(),pd.getDomainStart()) + 1;//
		} else if(bh.getqStart() < pd.getDomainStart() && bh.getqEnd() > pd.getDomainStart()) {
			ret = Math.min(bh.getqEnd(),pd.getDomainEnd()) - pd.getDomainStart() + 1;//
		} else if(bh.getqStart() > pd.getDomainStart() && bh.getqEnd() < pd.getDomainEnd()) {
			ret = bh.getqEnd() - bh.getqStart() + 1;
		} else if(bh.getqStart()==pd.getDomainStart()) {
			ret = Math.min(bh.getqEnd(), pd.getDomainEnd()) - bh.getqStart() + 1;
		} else if(bh.getqEnd()==pd.getDomainEnd()) {
			ret = bh.getqEnd() - Math.max(bh.getqStart(), pd.getDomainStart()) +1;
		} else if(bh.getqStart()==pd.getDomainEnd() || bh.getqEnd()==pd.getDomainStart()) {
			ret = 1;
		}
		return ret;
	}
	
	
}
