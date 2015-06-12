package tools;

import filter.FilterFacade;
import filter.impl.BestEvalueFilter;
import filter.impl.BoundShrinker;
import filter.impl.LowComplexityFilter;
import filter.impl.TooFewHitsFilter;
import filter.impl.TooLongHitsFilter;
import filter.impl.TooSmallDomainFilter;
import global.Global;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import model.BlastHit;
import model.PutativeDomain;
import model.Protein;

public class HitsGathering {

	private HitsGathering() {}

	public static Set<PutativeDomain> gatherHits(List<BlastHit> listHits, Protein protein, Set<PutativeDomain> forbiddenRegions) throws Exception {
		Set<PutativeDomain> ret = new HashSet<PutativeDomain>();
		
		if(forbiddenRegions==null) forbiddenRegions = new HashSet<PutativeDomain>(); //on first iteration
		
		//compute coverage
		int[] coverage = new int[protein.getSequence().length()];
		for(int i = 0; i < coverage.length; i++) coverage[i] = 0;

		for(BlastHit bh : listHits) {
			if(!Global.EXCLUDE_SPECIES_SET.contains(bh.getSubjectSpecies()) && canIuse(bh, forbiddenRegions)) {
				for(int i=bh.getqStart()-1; i<=bh.getqEnd()-1; i++) {
					coverage[i]++;
				}
			}
		}
		
//		//debug==
//		System.out.print("\n"+protein.getProtName()+"\ndata <- c(");
//		for(int i=0; i<coverage.length; i++) {
//			System.out.print(coverage[i]+",");
//		}
//		System.out.println(")");
//		//==debug
		
		//cut in putative domains
		Set<PutativeDomain> pDomains = new HashSet<PutativeDomain>();
		PutativeDomain currentPutativeDomain = null;
		for(int index = 0; index < coverage.length; index++) {
			if(coverage[index]>=Global.NUMBER_OF_HITS_MIN) {
				if(currentPutativeDomain==null) currentPutativeDomain = new PutativeDomain(protein.getProtName(),protein.getSpecies(),index+1, index+1);
				else currentPutativeDomain.setDomainEnd(currentPutativeDomain.getDomainEnd()+1);
			} else {
				if(currentPutativeDomain==null) continue;
				else {
					pDomains.add(currentPutativeDomain);
					currentPutativeDomain = null;
				}
			}
		}
		if(currentPutativeDomain!=null) pDomains.add(currentPutativeDomain);
		
		//supprime les domaines trop petits
		pDomains = FilterFacade.applyFilter(TooSmallDomainFilter.class, pDomains);
		
		//ajoute tous les hits dans les putative domains correspondant
		for(BlastHit bh : listHits) {
			for(PutativeDomain pd : pDomains) {
				if(isIn(bh, pd)) pd.addBlastHit(bh);
			}
		}
		
		//supprime les hits redondants
		pDomains = FilterFacade.applyFilter(BestEvalueFilter.class, pDomains);
		pDomains = FilterFacade.applyFilter(TooFewHitsFilter.class, pDomains);
		
		//raccourci les bornes, verifie si les hits sont encore bien dedans
		pDomains = FilterFacade.applyFilter(TooLongHitsFilter.class, pDomains);
		pDomains = FilterFacade.applyFilter(TooFewHitsFilter.class, pDomains);
		pDomains = FilterFacade.applyFilter(BoundShrinker.class, pDomains);
		pDomains = FilterFacade.applyFilter(TooSmallDomainFilter.class, pDomains);
		
		//low complexity filter
		pDomains = FilterFacade.applyFilter(LowComplexityFilter.class, pDomains);
		
		ret.addAll(pDomains);
		return ret;
	}

	
	
	/**
	 * Can I use this hit when computing the coverage
	 * @param bh
	 * @param pd
	 * @return
	 */
	private static boolean canIuse(BlastHit bh, Set<PutativeDomain> forbiddenRegions) {
		boolean ret = true;
		
		for(PutativeDomain domain : forbiddenRegions) {
			if(getOverlapSize(bh, domain)>0) {
				ret = false;
				break;
			}
		}
		
		return ret;
	}
	

	/**
	 * Determine if the blasthit is in a given putative domain
	 * @param bh
	 * @param pd
	 * @return
	 */
	private static boolean isIn(BlastHit bh, PutativeDomain pd) {
		boolean ret = false;

		if(getOverlapSize(bh, pd) >= (bh.getqEnd()-bh.getqStart())*Global.OVERLAP_RATE_MIN) { //if overlap by at least X%
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
