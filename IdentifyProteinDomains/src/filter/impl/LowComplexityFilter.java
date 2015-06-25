/**
 * 
 */
package filter.impl;

import java.util.HashSet;
import java.util.Set;

import model.Protein;
import model.PutativeDomain;
import filter.Ifilter;
import global.Global;

/**
 * @author christophe
 *
 */
public class LowComplexityFilter implements Ifilter {
	
	private static Ifilter instance;
	
	private final static int WINDOWSIZE = 20;
	private final static int STEPSIZE = 1;
	
	private LowComplexityFilter() {}
	
	public static Ifilter getInstance() {
		instance = (instance==null)?new LowComplexityFilter():instance;
		return instance;
	}
	
	public PutativeDomain filter(PutativeDomain domain) {
		PutativeDomain ret = null;
		
		Protein p = Global.PROTEOME_AIMED.getProteinByName(domain.getQueryName()+"_"+domain.getQuerySpecies());

		String sequence = "";
		for(int index=domain.getDomainStart(); index<domain.getDomainEnd()-1; index++) {
			sequence += p.getSequence().charAt(index);
		}
		
		int step = 0;
		String subseq;
		boolean isValid = true;
		int beginIndex = 0;
		while(beginIndex<sequence.length()-WINDOWSIZE-2) {
			beginIndex=step*STEPSIZE;
			subseq = sequence.substring(beginIndex, beginIndex+WINDOWSIZE);
			if(isLowComplexity(subseq)) {
				isValid = false;
				break;
			}
			step++;
		}

		if(isValid) ret = domain;
		//else System.err.println("LCF REJECTED >> "+sequence);

		return ret;
	}

	private boolean isLowComplexity(String sequence) {
		boolean ret = false;
		Set<Character> chars = new HashSet<Character>();
		for(char c : sequence.toCharArray()) {chars.add(c);}
		ret = (Math.log(chars.size())/Math.log(2) < 2.2); //Algo SEG
		return ret;
	}

}
