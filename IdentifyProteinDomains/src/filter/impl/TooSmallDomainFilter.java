package filter.impl;

import filter.Ifilter;
import global.Global;
import model.PutativeDomain;

public class TooSmallDomainFilter implements Ifilter {

	private static Ifilter instance;
	
	private TooSmallDomainFilter() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new TooSmallDomainFilter():instance;
		return instance;
	}
	
	@Override
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		
		if((domain.getDomainEnd()-domain.getDomainStart())>=Global.SIZE_HITS_MIN) {
			ret = domain;
		}
		
		return ret;
	}

}
