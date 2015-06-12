/**
 * 
 */
package filter.impl;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import filter.Ifilter;
import model.BlastHit;
import model.PutativeDomain;

/**
 * @author christophe
 *
 */
public class BestEvalueFilter implements Ifilter {

	private static Ifilter instance;

	private BestEvalueFilter() {}

	public static Ifilter getInstance() {
		instance=(instance==null)?new BestEvalueFilter():instance;
		return instance;
	}

	@Override
	public PutativeDomain filter(PutativeDomain domain) {
		Set<BlastHit> hits = domain.getBlastHits();
		List<BlastHit> listHits = new ArrayList<BlastHit>(hits);
		//sort by ascending order
		Collections.sort(listHits, new Comparator<BlastHit>() {
			public int compare(BlastHit o1, BlastHit o2) {
				if(o1.geteValue() > o2.geteValue()) return -1;
				if(o1.geteValue() < o2.geteValue()) return 1;
				return 0;
			}
		});

		//clean the putative domain
		domain.setBlastHits(new HashSet<BlastHit>());

		//replace hits in the putative domain
		for(BlastHit bh : listHits) {
			if(!isThereOverlap(bh, domain.getBlastHits())) {
				domain.addBlastHit(bh);
			}
		}
		
		//garde que les hits qui chevauchent la meilleur position
//		//reload hits
//		hits = domain.getBlastHits();
//		int bestPos = domain.getBestPosition();
//
//		domain.setBlastHits(new HashSet<BlastHit>());
//
//		//keep hits overlapping the best position
//		for(BlastHit bh : hits) {
//			if(bh.getqStart()<=bestPos && bh.getqEnd()>=bestPos) {
//				domain.addBlastHit(bh);
//			}
//		}

		return domain;
	}

	/**
	 * Is there an overlap of this hit with already added hits
	 * @param bh1
	 * @param set
	 * @return
	 */
	private boolean isThereOverlap(BlastHit h1, Set<BlastHit> set) {
		boolean ret = false;
		for(BlastHit h2 : set) {
			if(h1.getSubjectName().equals(h2.getSubjectName())) {
				if((h1.getsStart() < h2.getsStart() && h1.getsEnd() > h2.getsEnd()) 
						|| (h1.getsStart() < h2.getsEnd() && h1.getsEnd() > h2.getsEnd())
						|| (h1.getsStart() < h2.getsEnd() && h1.getsEnd() > h2.getsEnd())
						|| (h1.getsStart() < h2.getsStart() && h1.getsEnd() > h2.getsStart())
						|| (h1.getsStart() > h2.getsStart() && h1.getsEnd() < h2.getsEnd())
						|| (h1.getsStart()==h2.getsStart())
						|| (h1.getsEnd()==h2.getsEnd())
						|| (h1.getsStart()==h2.getsEnd() || h1.getsEnd()==h2.getsStart())) {
					ret = true;
				}
			}
		}
		return ret;
	}
}
