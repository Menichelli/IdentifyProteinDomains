package model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class PutativeDomain implements Cloneable {

	private String queryName;
	private String querySpecies;
	private int domainStart;
	private int domainEnd;
	private Set<BlastHit> blastHits;
	private boolean hasBeenCertified;
	
	private final static int BEST_POSITION_BOUND_MARGIN = 10;

	public PutativeDomain(String queryName, String querySpecies, int domStart, int domEnd) {
		super();
		this.queryName = queryName;
		this.querySpecies = querySpecies;
		this.domainStart = domStart;
		this.domainEnd = domEnd;
		this.blastHits = new HashSet<BlastHit>();
		this.hasBeenCertified = false;
	}

	/**
	 * @return the queryName
	 */
	public String getQueryName() {
		return queryName;
	}
	/**
	 * @param queryName the queryName to set
	 */
	public void setQueryName(String queryName) {
		this.queryName = queryName;
	}
	/**
	 * @return the querySpecies
	 */
	public String getQuerySpecies() {
		return querySpecies;
	}
	/**
	 * @param querySpecies the querySpecies to set
	 */
	public void setQuerySpecies(String querySpecies) {
		this.querySpecies = querySpecies;
	}
	/**
	 * @return the hitStart
	 */
	public int getDomainStart() {
		return domainStart;
	}
	/**
	 * @param domainStart the hitStart to set
	 */
	public void setDomainStart(int domainStart) {
		this.domainStart = domainStart;
	}
	/**
	 * @return the hitEnd
	 */
	public int getDomainEnd() {
		return domainEnd;
	}
	/**
	 * @param domainEnd the hitEnd to set
	 */
	public void setDomainEnd(int domainEnd) {
		this.domainEnd = domainEnd;
	}

	public Set<BlastHit> getBlastHits() {
		return this.blastHits;
	}

	public void setBlastHits(Set<BlastHit> set) {
		this.blastHits = (set==null)?new HashSet<BlastHit>():set;
	}

	public void addBlastHit(BlastHit bh) {
		this.blastHits.add(bh);
	}

	public String getIdentifier() {
		return queryName+"_"+querySpecies+"/"+domainStart+"-"+domainEnd;
	}

	public String getFileIdentifier() {
		return queryName+"_"+querySpecies+"_"+domainStart+"-"+domainEnd;
	}

	/**
	 * Nom de toutes les especes sur lesquelles notre putative a hit
	 * @return
	 */
	public Set<String> getHitSpecies() {
		Set<String> ret = new HashSet<String>();

		for(BlastHit bh : blastHits) {
			ret.add(bh.getSubjectSpecies());
		}

		return ret;
	}

	/**
	 * Permet de recuperer tous les noms de prot ayant un hit couvrant un residu particulier
	 * @param position
	 * @return
	 */
	public Set<String> getProteinsCoveringResidue(int position) {
		Set<String> ret = new HashSet<String>();

		for(BlastHit bh : this.blastHits) {
			if(position <= bh.getqEnd() && position >= bh.getqStart()) {
				ret.add(bh.getSubjectName()+"_"+bh.getSubjectSpecies());
			}
		}

		return ret;
	}

	public Set<BlastHit> getHitsCoveringResidue(int position) {
		Set<BlastHit> ret = new HashSet<BlastHit>();

		for(BlastHit bh : this.blastHits) {
			if(position <= bh.getqEnd() && position >= bh.getqStart()) {
				ret.add(bh);
			}
		}

		return ret;
	}

	public int getBestPosition() {
		int ret = 0;

		//compute coverage
		int[] coverage = new int[domainEnd-domainStart+1-2*BEST_POSITION_BOUND_MARGIN];
		int bestCoverage=0,tmp;
		for(int index = 0; index < coverage.length; index++) {
			tmp = getHitsCoveringResidue(index+domainStart+BEST_POSITION_BOUND_MARGIN).size();
			bestCoverage=(tmp>bestCoverage)?tmp:bestCoverage;
			coverage[index] = tmp;
		}

		List<Integer> bestPositions = new ArrayList<Integer>();
		for(int index = 0; index < coverage.length; index++) {
			if(coverage[index]==bestCoverage) bestPositions.add(new Integer(index+domainStart+BEST_POSITION_BOUND_MARGIN));
		}

		Collections.sort(bestPositions, new Comparator<Integer>() {
			public int compare(Integer o1, Integer o2) {
				if(o1>o2) return 1;
				if(o1<o2) return -1;
				return 0;
			}
		});

		int quantile50 = (int)(bestPositions.size()/2);
		ret = bestPositions.get(quantile50);
		
		return ret;
	}
	
	public void certify() throws Exception {
		if(this.hasBeenCertified) throw new Exception("Already certified!");
		this.hasBeenCertified = true;
	}
	
	public boolean hasBeenCertified() {
		return this.hasBeenCertified;
	}
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return String.format(
				"PutativeDomain [queryName=%s, querySpecies=%s, hitStart=%s, hitEnd=%s, nbHits=%d]",
				queryName, querySpecies, domainStart, domainEnd, getHitsCoveringResidue(getBestPosition()).size());
	}
	
	public Object clone() {
		PutativeDomain ret = new PutativeDomain(this.queryName, this.querySpecies, this.domainStart, this.domainEnd);
		
		for(BlastHit bh : this.blastHits) {
			ret.addBlastHit((BlastHit)bh.clone());
		}
		
		return ret;
	}
	
}
