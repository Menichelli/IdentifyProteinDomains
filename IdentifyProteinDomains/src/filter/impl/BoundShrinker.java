package filter.impl;

import filter.Ifilter;
import global.Global;
import model.PutativeDomain;

public class BoundShrinker implements Ifilter {

	private static Ifilter instance;
	
	private BoundShrinker() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new BoundShrinker():instance;
		return instance;
	}

	@Override
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		boolean modif = true;
		PutativeDomain tmp = domain;
		while(modif) {
			modif = false;
			if(tmp.getDomainEnd()<=tmp.getDomainStart()) break;
			if(tmp.getHitsCoveringResidue(tmp.getDomainStart()).size() <= Global.NB_SEQ_INTERSECT) {
				tmp.setDomainStart(tmp.getDomainStart()+1);
				modif = true;
			}			
			if(tmp.getHitsCoveringResidue(tmp.getDomainEnd()).size() <= Global.NB_SEQ_INTERSECT) {
				tmp.setDomainEnd(tmp.getDomainEnd()-1);
				modif = true;
			}
		}
		ret = tmp;
		return ret;
	}

}
